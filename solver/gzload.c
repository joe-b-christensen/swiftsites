/*
 *  SwiftSites: gzip input loader with two-thread split decoding.
 *  GPL-3.0 (part of SwiftSites, derived from snp-sites v2.5.1).
 *
 *  Sequential path: libdeflate, one gzip member after another while the next two bytes are the gzip magic; anything
 *  else after a member is ignored (zlib gzread's trailing-garbage rule).
 *
 *  Split path (one large member, >= 2 threads):
 *   1. find_boundary(): starting at the middle of the deflate data, try bit offsets until one starts a dynamic
 *      Huffman block whose header is fully valid (complete codes, end-of-block symbol present), which decodes to its
 *      end-of-block, and which is followed by another valid block header. Call it X.
 *   2. Thread 0 decodes from the start with libdeflate, stopping when a block starts exactly at bit X. Reaching X at
 *      a block start proves X is a real block boundary of this stream, so the second half below is the true
 *      continuation; if thread 0 does not stop at X the split is abandoned and the member is decoded sequentially.
 *   3. Thread 1 decodes from X with a small marker-aware inflate (16-bit output: 0..255 bytes, 256+w for byte w of
 *      the unknown 32 KB window before X) until, at a block boundary, its last 32 KB of output contain no markers;
 *      it then writes those 32 KB as real bytes and continues with libdeflate using them as a preset window.
 *   4. Markers are resolved from the last 32 KB of thread 0's output, thread 1's bytes are moved to follow thread 0's,
 *      and CRC32 and ISIZE from the gzip trailer are checked (any mismatch -> sequential decoding).
 */
#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <stdio.h>
#include <zlib.h>
#include "libdeflate/libdeflate.h"
#include "gzload.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/* libdeflate split extension (fields at the end of struct libdeflate_decompressor; see deflate_decompress.c).
   Accessed through a small accessor compiled inside libdeflate to avoid duplicating the struct layout. */
void libdeflate_sw_set(struct libdeflate_decompressor *d, size_t start_bit, size_t dict_len, size_t stop_bit);
int libdeflate_sw_stopped(struct libdeflate_decompressor *d);

static unsigned char *big_alloc2(size_t n)
{
  void *m = mmap(NULL, n, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
  if (m == MAP_FAILED) return NULL;
  madvise(m, n, MADV_HUGEPAGE);
  return m;
}

/* ------------------------------------------------------------------ gzip header */

/* Returns the offset of the deflate data, or 0 if the header is not a plain-parsable gzip header. */
static size_t gzip_header_len(const unsigned char *c, size_t n)
{
  if (n < 18 || c[0] != 31 || c[1] != 139 || c[2] != 8) return 0;
  unsigned flg = c[3];
  size_t p = 10;
  if (flg & 0xE0) return 0;
  if (flg & 4) { if (p + 2 > n) return 0; size_t xlen = c[p] | (c[p + 1] << 8); p += 2 + xlen; }
  if (flg & 8) { while (p < n && c[p]) p++; p++; }
  if (flg & 16) { while (p < n && c[p]) p++; p++; }
  if (flg & 2) p += 2;
  return p < n ? p : 0;
}

/* ------------------------------------------------------------------ marker-aware inflate (slow, simple) */

typedef struct { const unsigned char *c; size_t nbits; size_t pos; } bits_t;

static inline uint32_t peek(const bits_t *b, unsigned k) /* k <= 25 */
{
  size_t byte = b->pos >> 3; unsigned sh = b->pos & 7;
  uint64_t v = 0;
  size_t nbytes = (b->nbits + 7) >> 3;
  for (unsigned i = 0; i < 5; i++) if (byte + i < nbytes) v |= (uint64_t) b->c[byte + i] << (8 * i);
  return (uint32_t) ((v >> sh) & ((1u << k) - 1));
}
static inline uint32_t getb(bits_t *b, unsigned k) { uint32_t v = peek(b, k); b->pos += k; return v; }

/* Huffman table indexed by 15 reversed bits: entry = symbol | len << 9 (len 0 = invalid). */
typedef struct { uint16_t t[1 << 15]; int empty; } huff_t;

static uint32_t rev(uint32_t code, unsigned len) { uint32_t r = 0; for (unsigned i = 0; i < len; i++) { r = (r << 1) | (code & 1); code >>= 1; } return r; }

/* Builds the table; returns 0 if the code is over-subscribed, or incomplete unless allow_single (one code of
   length 1, as zlib allows for length/distance codes), or (for distance codes, allow_empty) all lengths zero. */
static int huff_build(huff_t *h, const uint8_t *lens, int n, int allow_single, int allow_empty)
{
  int count[16] = {0};
  for (int i = 0; i < n; i++) count[lens[i]]++;
  count[0] = 0;
  int total = 0, maxlen = 0;
  for (int l = 1; l < 16; l++) { total += count[l]; if (count[l]) maxlen = l; }
  memset(h->t, 0, sizeof h->t);
  h->empty = total == 0;
  if (total == 0) return allow_empty;
  long left = 1;
  for (int l = 1; l < 16; l++) { left <<= 1; left -= count[l]; if (left < 0) return 0; }
  if (left > 0 && !(allow_single && maxlen == 1)) return 0;
  uint32_t next[16]; uint32_t code = 0;
  for (int l = 1; l < 16; l++) { code = (code + (uint32_t) count[l - 1]) << 1; next[l] = code; }
  for (int s = 0; s < n; s++) {
    unsigned l = lens[s];
    if (!l) continue;
    uint32_t r = rev(next[l]++, l);
    for (uint32_t k = r; k < (1u << 15); k += 1u << l) h->t[k] = (uint16_t) (s | (l << 9));
  }
  return 1;
}

/* Precode (code-length code): max 7-bit codes, so a 128-entry table. Must be complete (as zlib requires). */
typedef struct { uint16_t t[128]; } pre_t;
static int pre_build(pre_t *h, const uint8_t *lens)
{
  int count[8] = {0};
  for (int i = 0; i < 19; i++) count[lens[i]]++;
  count[0] = 0;
  long left = 1;
  for (int l = 1; l < 8; l++) { left <<= 1; left -= count[l]; if (left < 0) return 0; }
  if (left != 0) return 0;
  uint32_t next[8]; uint32_t code = 0;
  for (int l = 1; l < 8; l++) { code = (code + (uint32_t) count[l - 1]) << 1; next[l] = code; }
  for (int sy = 0; sy < 19; sy++) {
    unsigned l = lens[sy];
    if (!l) continue;
    uint32_t r = rev(next[l]++, l);
    for (uint32_t k = r; k < 128; k += 1u << l) h->t[k] = (uint16_t) (sy | (l << 9));
  }
  return 1;
}

static inline int hdecode(bits_t *b, const huff_t *h)
{
  uint16_t e = h->t[peek(b, 15)];
  unsigned l = e >> 9;
  if (!l || b->pos + l > b->nbits) return -1;
  b->pos += l;
  return e & 511;
}

static const uint16_t LBASE[29] = {3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258};
static const uint8_t LEXT[29] = {0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0};
static const uint16_t DBASE[30] = {1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577};
static const uint8_t DEXT[30] = {0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13};

typedef struct {
  huff_t lit, dist;
  uint16_t *out; size_t n, cap; size_t last_marker_end; /* index after the last marker written */
} m16_t;

/* Reads a block header at b->pos and builds tables. type: 0 stored, 1 fixed, 2 dynamic. Returns 0 if invalid. */
static int read_header(bits_t *b, m16_t *m, int *final, int *type, uint32_t *stored_len)
{
  if (b->pos + 3 > b->nbits) return 0;
  *final = (int) getb(b, 1);
  *type = (int) getb(b, 2);
  if (*type == 3) return 0;
  if (*type == 0) {
    b->pos = (b->pos + 7) & ~(size_t) 7;
    if (b->pos + 32 > b->nbits) return 0;
    uint32_t len = getb(b, 16), nlen = getb(b, 16);
    if ((len ^ 0xFFFF) != nlen) return 0;
    *stored_len = len;
    return 1;
  }
  uint8_t lens[320];
  if (*type == 1) {
    for (int i = 0; i < 144; i++) lens[i] = 8;
    for (int i = 144; i < 256; i++) lens[i] = 9;
    for (int i = 256; i < 280; i++) lens[i] = 7;
    for (int i = 280; i < 288; i++) lens[i] = 8;
    huff_build(&m->lit, lens, 288, 0, 0);
    for (int i = 0; i < 30; i++) lens[i] = 5;
    huff_build(&m->dist, lens, 30, 1, 0);
    return 1;
  }
  if (b->pos + 14 > b->nbits) return 0;
  int nlen = (int) getb(b, 5) + 257, ndist = (int) getb(b, 5) + 1, ncode = (int) getb(b, 4) + 4;
  if (nlen > 286 || ndist > 30) return 0;
  static const uint8_t ORD[19] = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};
  uint8_t pl[19] = {0};
  if (b->pos + 3u * (unsigned) ncode > b->nbits) return 0;
  for (int i = 0; i < ncode; i++) pl[ORD[i]] = (uint8_t) getb(b, 3);
  pre_t pre;
  if (!pre_build(&pre, pl)) return 0;
  int k = 0;
  while (k < nlen + ndist) {
    uint16_t e = pre.t[peek(b, 7)];
    unsigned el = e >> 9;
    if (!el) return 0;
    b->pos += el;
    int sym = e & 511;
    if (sym < 16) lens[k++] = (uint8_t) sym;
    else {
      int rep, val = 0;
      if (sym == 16) { if (k == 0) return 0; val = lens[k - 1]; rep = 3 + (int) getb(b, 2); }
      else if (sym == 17) rep = 3 + (int) getb(b, 3);
      else rep = 11 + (int) getb(b, 7);
      if (k + rep > nlen + ndist) return 0;
      while (rep--) lens[k++] = (uint8_t) val;
    }
    if (b->pos > b->nbits) return 0;
  }
  if (lens[256] == 0) return 0;
  if (!huff_build(&m->dist, lens + nlen, ndist, 1, 1)) return 0;
  if (!huff_build(&m->lit, lens, nlen, 1, 0)) return 0;
  return 1;
}

static int m_put(m16_t *m, uint16_t v)
{
  if (m->n == m->cap) {
    size_t nc = m->cap ? m->cap * 2 : (1u << 20);
    uint16_t *o = realloc(m->out, nc * sizeof(uint16_t));
    if (!o) return 0;
    m->out = o; m->cap = nc;
  }
  m->out[m->n++] = v;
  if (v >= 256) m->last_marker_end = m->n;
  return 1;
}

/* Decodes the body of a block whose header was just read. Returns 1 at end of block, 0 on invalid data. */
static int decode_body(bits_t *b, m16_t *m, int type, uint32_t stored_len)
{
  if (type == 0) {
    if (b->pos + 8u * stored_len > b->nbits) return 0;
    for (uint32_t i = 0; i < stored_len; i++) if (!m_put(m, (uint16_t) getb(b, 8))) return 0;
    return 1;
  }
  for (;;) {
    int sym = hdecode(b, &m->lit);
    if (sym < 0) return 0;
    if (sym < 256) { if (!m_put(m, (uint16_t) sym)) return 0; continue; }
    if (sym == 256) return 1;
    sym -= 257;
    if (sym >= 29) return 0;
    unsigned len = LBASE[sym] + getb(b, LEXT[sym]);
    if (m->dist.empty) return 0;
    int ds = hdecode(b, &m->dist);
    if (ds < 0 || ds >= 30) return 0;
    size_t dist = DBASE[ds] + getb(b, DEXT[ds]);
    if (b->pos > b->nbits) return 0;
    if (dist > m->n + 32768) return 0;
    for (unsigned i = 0; i < len; i++) {
      long src = (long) m->n - (long) dist;
      uint16_t v = src >= 0 ? m->out[src] : (uint16_t) (256 + (32768 + src));
      if (!m_put(m, v)) return 0;
    }
  }
}

/* Is there a plausible block boundary at bit x? (valid dynamic header, block decodes, next header valid) */
static int try_boundary(const unsigned char *c, size_t nbits, size_t x, m16_t *m)
{
  bits_t b = { c, nbits, x };
  int final, type; uint32_t sl = 0;
  if (((peek(&b, 3) >> 1) & 3) != 2 || (peek(&b, 1) & 1)) return 0; /* non-final dynamic block only */
  if (!read_header(&b, m, &final, &type, &sl)) return 0;
  m->n = 0; m->last_marker_end = 0;
  if (!decode_body(&b, m, type, sl)) return 0;
  if (m->n < 1024) return 0;
  return read_header(&b, m, &final, &type, &sl);
}

/* ------------------------------------------------------------------ split decoding */

typedef struct {
  const unsigned char *c; size_t dstart, dbytes; /* deflate data */
  unsigned char *buf; size_t cap, G;              /* output buffer; thread 1 writes from G */
  size_t X;                                       /* split bit (relative to deflate start) */
  /* results */
  int ok0, ok1; size_t A_len, B_len, slow_len, end_byte; uint16_t *slow;
  uint32_t crcA, crcBfast; /* CRC32 of thread 0's bytes and of thread 1's libdeflate part */
} split_t;

static void thread0(split_t *s)
{
  struct libdeflate_decompressor *d = libdeflate_alloc_decompressor();
  libdeflate_sw_set(d, 0, 0, s->X);
  size_t ain = 0, aout = 0;
  enum libdeflate_result rc = libdeflate_deflate_decompress_ex(d, s->c + s->dstart, s->dbytes, s->buf, s->G, &ain, &aout);
  s->ok0 = rc == LIBDEFLATE_SUCCESS && libdeflate_sw_stopped(d) && aout >= 32768;
  s->A_len = aout;
  libdeflate_free_decompressor(d);
  if (s->ok0) s->crcA = libdeflate_crc32(0, s->buf, aout);
}

static void thread1(split_t *s)
{
  s->ok1 = 0;
  m16_t *m = calloc(1, sizeof *m);
  bits_t b = { s->c + s->dstart, s->dbytes * 8, s->X };
  int final = 0, type; uint32_t sl = 0;
  for (;;) {
    if (!read_header(&b, m, &final, &type, &sl)) goto out;
    if (!decode_body(&b, m, type, sl)) goto out;
    if (final) { /* the whole second part was decoded slowly */
      s->slow = m->out; s->slow_len = m->n; s->B_len = m->n; s->crcBfast = 0;
      s->end_byte = (b.pos + 7) >> 3;
      if (s->G + m->n > s->cap) goto out;
      s->ok1 = 1; free(m); return;
    }
    if (m->n >= 32768 && m->n - m->last_marker_end >= 32768) break;
    if (m->n > (64u << 20)) goto out; /* window dependence does not die out: give up */
  }
  /* switch to libdeflate with the last 32 KB as a preset window */
  size_t sl_n = m->n;
  if (s->G + sl_n > s->cap) goto out;
  for (size_t i = sl_n - 32768; i < sl_n; i++) s->buf[s->G + i] = (unsigned char) m->out[i];
  {
    struct libdeflate_decompressor *d = libdeflate_alloc_decompressor();
    libdeflate_sw_set(d, b.pos, 32768, (size_t) -1);
    size_t ain = 0, aout = 0;
    enum libdeflate_result rc = libdeflate_deflate_decompress_ex(d, s->c + s->dstart, s->dbytes, s->buf + s->G + sl_n,
                                                                 s->cap - s->G - sl_n, &ain, &aout);
    libdeflate_free_decompressor(d);
    if (rc != LIBDEFLATE_SUCCESS) goto out;
    s->crcBfast = libdeflate_crc32(0, s->buf + s->G + sl_n, aout);
    s->slow = m->out; s->slow_len = sl_n; s->B_len = sl_n + aout; s->end_byte = ain;
    s->ok1 = 1; free(m); return;
  }
out:
  free(m->out); free(m);
}

static int pgz_min_bytes(void)
{
  const char *e = getenv("SWIFTSITES_PGZ_MIN"); /* test hook: force the split path on small inputs */
  return e ? atoi(e) : (32 << 20);
}

/* Split-decode the first member. Returns 1 and sets out, outlen and consumed on success, 0 to fall back. */
static int split_member(const unsigned char *c, size_t n, unsigned char **out, size_t *outlen, size_t *consumed,
                        size_t *outcap)
{
  size_t h = gzip_header_len(c, n);
  if (!h || n < h + 8 + 4096) return 0;
  uint32_t isize; memcpy(&isize, c + n - 4, 4);
  if (isize < (uint32_t) pgz_min_bytes()) return 0;
  size_t dbytes = n - h - 8; /* deflate data upper bound (single member assumed; verified below) */
  double ts0 = omp_get_wtime();
  m16_t *m = calloc(1, sizeof *m);
  size_t start = (dbytes / 2) * 8, X = 0;
  size_t limit = start + (size_t) 8 * (8u << 20);
  for (size_t x = start; x < limit && x < dbytes * 8 - 64; x++)
    if (try_boundary(c + h, dbytes * 8, x, m)) { X = x; break; }
  free(m->out); free(m);
  if (getenv("SWIFTSITES_PGZ_DEBUG")) fprintf(stderr, "PGZ search %.4f s (X-start=%zu bits)\n", omp_get_wtime() - ts0, X - start);
  if (!X) return 0;
  split_t s; memset(&s, 0, sizeof s);
  s.c = c; s.dstart = h; s.dbytes = dbytes; s.X = X;
  double frac = (double) X / (double) (dbytes * 8);
  size_t est = (size_t) (frac * (double) isize);
  s.G = est + est / 20 + (4u << 20);
  s.cap = s.G + (size_t) isize + (64u << 20);
  s.buf = big_alloc2(s.cap);
  if (!s.buf) return 0;
  double tt0 = omp_get_wtime(), t0e = 0, t1e = 0;
#pragma omp parallel sections num_threads(2)
  {
#pragma omp section
    { thread0(&s); t0e = omp_get_wtime(); }
#pragma omp section
    { thread1(&s); t1e = omp_get_wtime(); }
  }
  if (getenv("SWIFTSITES_PGZ_DEBUG"))
    fprintf(stderr, "PGZ X=%zu/%zu thread0 %.3f s (A=%zu ok=%d) thread1 %.3f s (slow=%zu B=%zu ok=%d)\n", X, dbytes * 8,
            t0e - tt0, s.A_len, s.ok0, t1e - tt0, s.slow_len, s.B_len, s.ok1);
  int ok = s.ok0 && s.ok1 && s.A_len <= s.G && s.A_len >= 32768;
  double tr0 = omp_get_wtime();
  if (ok) {
    /* resolve markers from thread 0's last 32 KB, then move thread 1's bytes behind thread 0's */
    const unsigned char *win = s.buf + s.A_len - 32768;
    unsigned char *dst = s.buf + s.G;
    for (size_t i = 0; i < s.slow_len; i++) {
      uint16_t v = s.slow[i];
      dst[i] = v < 256 ? (unsigned char) v : win[v - 256];
    }
    uint32_t crcB = (uint32_t) crc32_combine(libdeflate_crc32(0, dst, s.slow_len), s.crcBfast, (z_off_t) (s.B_len - s.slow_len));
    /* Move thread 1's bytes down by gap = G - A_len, in chunks of cs = floor(gap/2). Chunk i's destination overlaps
       only the sources of chunks i-2 and i-3 (i-3 by at most one byte, when gap is odd), so chunks are processed in
       rounds {2k, 2k+1} with a barrier between rounds: both earlier chunks finished in previous rounds, and no later
       round writes over a source that has not been read yet. */
    size_t gap = s.G - s.A_len;
    if (gap >= (1u << 16)) {
      size_t cs = gap / 2;
      long nch = (long) ((s.B_len + cs - 1) / cs);
#pragma omp parallel num_threads(2)
      {
        int nt = omp_get_num_threads(), tid = omp_get_thread_num();
        for (long k = 0; 2 * k < nch; k++) {
          for (long i = 2 * k + (nt == 1 ? 0 : tid); i < 2 * k + 2 && i < nch; i += (nt == 1 ? 1 : 2)) {
            size_t a = (size_t) i * cs, len = s.B_len - a < cs ? s.B_len - a : cs;
            memmove(s.buf + s.A_len + a, s.buf + s.G + a, len);
          }
#pragma omp barrier
        }
      }
    } else if (gap) memmove(s.buf + s.A_len, s.buf + s.G, s.B_len);
    size_t total = s.A_len + s.B_len;
    size_t tpos = h + s.end_byte; /* gzip trailer */
    ok = tpos + 8 <= n;
    if (ok) {
      uint32_t crc, isz;
      memcpy(&crc, c + tpos, 4); memcpy(&isz, c + tpos + 4, 4);
      uint32_t cc = (uint32_t) crc32_combine(s.crcA, crcB, (z_off_t) s.B_len);
      ok = cc == crc && isz == (uint32_t) total;
      if (getenv("SWIFTSITES_PGZ_DEBUG")) fprintf(stderr, "PGZ resolve+move+crc %.4f s (gap %zu)\n", omp_get_wtime() - tr0, s.G - s.A_len);
      if (ok) { *out = s.buf; *outlen = total; *consumed = tpos + 8; *outcap = s.cap; }
    }
  }
  free(s.slow);
  if (!ok) munmap(s.buf, s.cap);
  return ok;
}

/* ------------------------------------------------------------------ public entry */

int gz_split_used = 0;

int gz_load(const unsigned char *c, size_t csize, unsigned char **outp, size_t *outlenp)
{
  int threads = 1;
#ifdef _OPENMP
  threads = omp_get_max_threads();
#endif
  uint32_t isize = 0;
  if (csize >= 4) memcpy(&isize, c + csize - 4, 4);
  unsigned char *out = NULL; size_t len = 0, off = 0, cap = 0;
  gz_split_used = 0;
  if (threads >= 2 && split_member(c, csize, &out, &len, &off, &cap)) {
    gz_split_used = 1;
    if (!(csize - off >= 2 && c[off] == 31 && c[off + 1] == 139)) { *outp = out; *outlenp = len; return 1; }
    /* more members follow: continue sequentially below */
  } else {
    cap = (size_t) isize + (1u << 20);
    out = big_alloc2(cap);
    if (!out) return 0;
  }
  struct libdeflate_decompressor *d = libdeflate_alloc_decompressor();
  int ok = 1;
  for (;;) {
    size_t ain = 0, aout = 0;
    enum libdeflate_result rc = libdeflate_gzip_decompress_ex(d, c + off, csize - off, out + len, cap - len, &ain, &aout);
    if (rc == LIBDEFLATE_INSUFFICIENT_SPACE) {
      size_t ncap = cap * 2;
      unsigned char *n2 = big_alloc2(ncap);
      if (!n2) { ok = 0; break; }
      memcpy(n2, out, len);
      munmap(out, cap);
      out = n2; cap = ncap;
      continue;
    }
    if (rc != LIBDEFLATE_SUCCESS) { ok = 0; break; }
    off += ain; len += aout;
    if (csize - off >= 2 && c[off] == 31 && c[off + 1] == 139) continue;
    break;
  }
  libdeflate_free_decompressor(d);
  if (!ok) { munmap(out, cap); return 0; }
  *outp = out; *outlenp = len;
  return 1;
}
