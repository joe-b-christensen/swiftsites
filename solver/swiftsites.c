/*
 *  SwiftSites: a byte-identical, faster drop-in for snp-sites v2.5.1.
 *
 *  Derived from snp-sites (Wellcome Trust Sanger Institute, 2013-2016; Page et al., Microbial Genomics 2016),
 *  GPL-3.0. This program is free software; you can redistribute it and/or modify it under the terms of the GNU
 *  General Public License as published by the Free Software Foundation; either version 3 of the License, or (at
 *  your option) any later version.
 *
 *  Design: the input is read once into memory (mmap for plain files, one decompression for gzip). An in-memory
 *  re-implementation of the 2012 kseq parser used by snp-sites produces the record list with identical semantics
 *  (names end at the first isspace byte, '>'/'@' start records, '+' switches to FASTQ, empty lines skipped, the
 *  "strip a trailing \r only when the accumulated length is > 1" rule, -2 on bad FASTQ quality ends the input).
 *  Single-line sequences are used in place; multi-line sequences are copied into an arena.
 *  Column states are computed in column blocks over all records; outputs are generated straight from the record
 *  slices (FASTA/PHYLIP row gathers, VCF from 64-column transposed tiles).
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <getopt.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/uio.h>
#include <sched.h>
#include <zlib.h>
#include <immintrin.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "libdeflate/libdeflate.h"
#include "gzload.h"

#define MAXIMUM_NUMBER_OF_ALT_BASES 30

/* Per-phase wall-clock timers, printed to stderr when SWIFTSITES_PROFILE is set (used for the paper's breakdown). */
#include <time.h>
static int g_prof = 0; static double g_t_last = 0;
static double now_s(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }
static void prof_mark(const char *phase)
{
  if (!g_prof) return;
  double t = now_s();
  fprintf(stderr, "PROFILE %-10s %.4f\n", phase, t - g_t_last);
  g_t_last = t;
}

typedef struct {
  const unsigned char *name; size_t name_len;
  const unsigned char *seq; size_t len;
} rec_t;

static inline unsigned char up(unsigned char c) { return (unsigned char) (c - (((unsigned char) (c - 'a') < 26) ? 32 : 0)); }
static inline int unk(unsigned char c) { return c == 'N' || c == 'n' || c == '-' || c == '?'; }
static inline int ks_isspace(unsigned char c) { return c == ' ' || (c >= '\t' && c <= '\r'); }

static void die(const char *msg) { fprintf(stderr, "%s", msg); fflush(stderr); exit(EXIT_FAILURE); }

/* ------------------------------------------------------------------ input */

typedef struct { const unsigned char *buf; size_t n; } input_t;

static int is_gzip_file(int fd)
{
  unsigned char m[2];
  ssize_t k = pread(fd, m, 2, 0);
  return k == 2 && m[0] == 31 && m[1] == 139;
}

static void load_gzip_zlib(const char *fn, input_t *in)
{
  gzFile z = gzopen(fn, "r");
  if (!z) die("ERROR: cannot open input\n");
  gzbuffer(z, 1 << 20);
  size_t cap = 1 << 26, n = 0;
  unsigned char *b = malloc(cap);
  for (;;) {
    if (cap - n < (1u << 24)) { cap *= 2; b = realloc(b, cap); }
    int k = gzread(z, b + n, (unsigned) (cap - n > (1u << 30) ? (1u << 30) : cap - n));
    if (k < 0) die("ERROR: gzip decompression failed\n");
    if (k == 0) break;
    n += (size_t) k;
  }
  gzclose(z);
  in->buf = b; in->n = n;
}


/* gzip input: see gzload.c (libdeflate, two-thread split decoding for one large member). Returns 0 if libdeflate
   rejects the data, in which case the caller falls back to zlib. */
static int load_gzip_libdeflate(int fd, size_t csize, input_t *in)
{
  unsigned char *c = mmap(NULL, csize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0);
  if (c == MAP_FAILED) return 0;
  unsigned char *out = NULL; size_t len = 0;
  int ok = gz_load(c, csize, &out, &len);
  munmap(c, csize);
  if (!ok) return 0;
  in->buf = out; in->n = len;
  return 1;
}

static void load_input(const char *fn, input_t *in, int force_private)
{
  int fd = open(fn, O_RDONLY);
  if (fd < 0) die("ERROR: cannot open input alignment file\n");
  struct stat st;
  if (fstat(fd, &st) != 0) die("ERROR: cannot stat input\n");
  size_t n = (size_t) st.st_size;
  if (is_gzip_file(fd)) {
    if (S_ISREG(st.st_mode) && load_gzip_libdeflate(fd, n, in)) { close(fd); return; }
    close(fd); load_gzip_zlib(fn, in); return;
  }
  if (n == 0 || !S_ISREG(st.st_mode) || force_private) {
    /* read(2) into private memory (non-regular files, or an output path that aliases the input) */
    size_t cap = n ? n + 1 : (1 << 20), len = 0;
    unsigned char *b = malloc(cap);
    for (;;) {
      if (len == cap) { cap *= 2; b = realloc(b, cap); }
      ssize_t k = read(fd, b + len, cap - len);
      if (k < 0) die("ERROR: read failed\n");
      if (k == 0) break;
      len += (size_t) k;
    }
    close(fd);
    in->buf = b; in->n = len;
    return;
  }
  /* map, then populate the page tables from all threads at once (MADV_POPULATE_READ, Linux 5.14+; takes the mmap
     lock shared); if that is unavailable, populate from this thread as before */
  void *m = mmap(NULL, n, PROT_READ, MAP_PRIVATE, fd, 0);
  if (m == MAP_FAILED) die("ERROR: mmap failed\n");
  close(fd);
#ifndef MADV_POPULATE_READ
#define MADV_POPULATE_READ 22
#endif
  int popfail = 0;
  if (n >= (64u << 20)) {
#pragma omp parallel reduction(|:popfail)
    {
      size_t T = (size_t) omp_get_num_threads(), t = (size_t) omp_get_thread_num();
      size_t a = n / T * t & ~(size_t) 4095, e = t + 1 == T ? n : n / T * (t + 1) & ~(size_t) 4095;
      if (e > a && madvise((char *) m + a, e - a, MADV_POPULATE_READ) != 0) popfail = 1;
    }
  } else popfail = madvise(m, n, MADV_POPULATE_READ) != 0;
  if (popfail) { /* older kernel: touch every page (the effect of MAP_POPULATE) */
    const volatile unsigned char *pv = m; unsigned char acc = 0;
    for (size_t i = 0; i < n; i += 4096) acc ^= pv[i];
    (void) acc;
  }
  in->buf = m; in->n = n;
}

/* ------------------------------------------------------------------ kseq-equivalent parser */

typedef struct {
  const unsigned char *b; size_t n, p; int last_char;
  unsigned char *arena; size_t arena_used, arena_cap;
  unsigned char *qual; size_t qual_cap;
} kq_t;

static unsigned char *arena_alloc_tail(kq_t *k, size_t need)
{
  if (!k->arena) {
    k->arena_cap = k->n + 64;
    k->arena = mmap(NULL, k->arena_cap, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    if (k->arena == MAP_FAILED) die("ERROR: out of memory\n");
    madvise(k->arena, k->arena_cap, MADV_HUGEPAGE);
  }
  if (k->arena_cap - k->arena_used < need) die("ERROR: arena overflow\n");
  return k->arena + k->arena_used;
}

/* Returns sequence length (>= 0), -1 at end of input, -2 on a FASTQ quality error (ends the input, as in kseq). */
static long kq_read(kq_t *k, rec_t *r)
{
  const unsigned char *b = k->b;
  size_t n = k->n, p = k->p;
  int c;
  if (k->last_char == 0) {
    while (p < n && b[p] != '>' && b[p] != '@') p++;
    if (p >= n) { k->p = n; return -1; }
    k->last_char = b[p]; p++;
  }
  /* name: ks_getuntil(KS_SEP_SPACE) */
  if (p >= n) { k->p = p; return -1; }
  size_t ns = p;
  while (p < n && !ks_isspace(b[p])) p++;
  r->name = b + ns; r->name_len = p - ns;
  if (p < n) { c = b[p]; p++; } else c = 0;
  if (c != '\n' && p < n) { /* comment: rest of the line */
    const unsigned char *q = memchr(b + p, '\n', n - p);
    p = q ? (size_t) (q - b) + 1 : n;
  }
  /* sequence lines */
  const unsigned char *single = NULL; size_t len = 0;
  unsigned char *dst = NULL; /* non-NULL once copying into the arena */
  c = -1;
  while (p < n) {
    c = b[p++];
    if (c == '>' || c == '+' || c == '@') break;
    if (c == '\n') { c = -1; continue; }
    size_t ls = p - 1, le; int strip_ok;
    if (p >= n) { le = n; strip_ok = 0; } /* ks_getuntil2 returns -1 at EOF: no \r strip */
    else {
      const unsigned char *q = memchr(b + p, '\n', n - p);
      le = q ? (size_t) (q - b) : n;
      p = q ? le + 1 : n;
      strip_ok = 1;
    }
    if (!single && !dst) { /* first content line: keep in place */
      single = b + ls; len = le - ls;
      if (strip_ok && len > 1 && single[len - 1] == '\r') len--;
    } else {
      if (!dst) { /* second content line: switch to the arena */
        dst = arena_alloc_tail(k, (n - ls) + len + 1);
        memcpy(dst, single, len);
        single = NULL;
      }
      memcpy(dst + len, b + ls, le - ls);
      len += le - ls;
      if (strip_ok && len > 1 && dst[len - 1] == '\r') len--;
    }
    c = -1;
  }
  if (c == '>' || c == '@') k->last_char = c;
  if (dst) { r->seq = dst; k->arena_used += len; } else r->seq = single ? single : (const unsigned char *) "";
  r->len = len;
  if (c != '+') { k->p = p; return (long) len; }
  /* FASTQ: skip the rest of the '+' line, then read quality lines until qual.l >= seq.l */
  while (p < n && b[p] != '\n') p++;
  if (p >= n) { k->p = p; return -2; }
  p++;
  size_t ql = 0;
  for (;;) {
    if (p >= n) break; /* ks_getuntil2 returns -1 */
    const unsigned char *q = memchr(b + p, '\n', n - p);
    size_t le = q ? (size_t) (q - b) : n;
    size_t add = le - p;
    if (k->qual_cap < ql + add + 1) { k->qual_cap = (ql + add + 1) * 2; k->qual = realloc(k->qual, k->qual_cap); }
    memcpy(k->qual + ql, b + p, add);
    ql += add;
    p = q ? le + 1 : n;
    if (ql > 1 && k->qual[ql - 1] == '\r') ql--;
    if (!(ql < len)) break;
  }
  k->last_char = 0;
  k->p = p;
  if (ql != len) return -2;
  return (long) len;
}

/* ------------------------------------------------------------------ fast layout path */

/* Most alignments are either one line per sequence or wrapped at one fixed width. parse_fast() walks record 1 line
   by line to get its length L and line width W, then predicts every later record's extent from (L, W) and checks
   only its header and the byte after its predicted end ('>' or '@' or end of input). validate_layout() then checks
   the predicted interior: every expected line end is '\n', the record holds no other '\n', no line starts with
   '>', '+' or '@' and no line ends with '\r'. Under those conditions kseq would return exactly these L content
   bytes for each record, so the records are used in place without copying. Anything else (CRLF, blank lines,
   FASTQ, ragged wrapping, unequal lengths, a missing final newline) returns 0 and the exact parser runs. */
static int parse_fast(const unsigned char *b, size_t n, rec_t **Rp, size_t *nrp, size_t *Lp, size_t *Wp)
{
  size_t p = 0;
  while (p < n && b[p] != '>' && b[p] != '@') p++;
  if (p >= n) return 0;
  size_t cap = 1024, nr = 0, L = 0, W = 0, nlines = 0;
  rec_t *R = malloc(cap * sizeof(rec_t));
  while (p < n) {
    /* b[p] is '>' or '@' */
    p++;
    size_t ns = p;
    while (p < n && !ks_isspace(b[p])) p++;
    if (p >= n) goto fail;
    size_t ne = p;
    if (b[p++] != '\n') {
      const unsigned char *q = memchr(b + p, '\n', n - p);
      if (!q) goto fail;
      p = (size_t) (q - b) + 1;
    }
    size_t ss = p;
    if (nr == 0) {
      /* measure record 1 exactly */
      size_t lp = p, total = 0, k = 0;
      for (;;) {
        if (lp >= n) break;
        unsigned char fc = b[lp];
        if (fc == '>' || fc == '@') break;
        if (fc == '+' || fc == '\n') goto fail;
        const unsigned char *q = memchr(b + lp, '\n', n - lp);
        if (!q) goto fail;
        size_t ll = (size_t) (q - b) - lp;
        if (b[lp + ll - 1] == '\r') goto fail;
        if (k == 0) W = ll;
        else if (total % W != 0 || ll > W) goto fail; /* every line but the last must be exactly W */
        total += ll; lp += ll + 1; k++;
        if (ll < W) { /* short line: must be the last one */
          if (lp < n && b[lp] != '>' && b[lp] != '@') goto fail;
          break;
        }
      }
      if (k == 0 || total == 0) goto fail;
      L = total; nlines = k;
      if (nlines == 1) W = L;
      p = lp;
    } else {
      size_t end = ss + L + nlines; /* one '\n' per line */
      if (end > n || b[end - 1] != '\n') goto fail;
      if (b[ss + (W < L ? W : L)] != '\n') goto fail;
      if (end < n && b[end] != '>' && b[end] != '@') goto fail;
      p = end;
    }
    if (nr == cap) { cap *= 2; R = realloc(R, cap * sizeof(rec_t)); }
    R[nr].name = b + ns; R[nr].name_len = ne - ns; R[nr].seq = b + ss; R[nr].len = L;
    nr++;
  }
  *Rp = R; *nrp = nr; *Lp = L; *Wp = W;
  return 1;
fail:
  free(R);
  return 0;
}

/* ------------------------------------------------------------------ column states */

/* f[j]: first non-unknown base (upper-cased) of column j, 'N' if none; var[j] != 0: some later non-unknown base
   differs from f[j]; np[j] != 0 (pure mode): some byte is not ACGTacgt. Processed in column blocks so the state
   arrays stay in L1; AVX2 kernel with a scalar tail (identical logic). */


static inline void scan_scalar(const unsigned char *s, size_t len, unsigned char *f, unsigned char *var,
                               unsigned char *np)
{
  for (size_t j = 0; j < len; j++) {
    unsigned char c = s[j];
    unsigned char u = up(c);
    unsigned char known = (u != 'N') & (u != '-') & (u != '?');
    unsigned char fj = f[j];
    fj = ((fj == 'N') & known) ? u : fj;
    f[j] = fj;
    var[j] |= known & (u != fj);
    if (np) np[j] |= (u != 'A') & (u != 'C') & (u != 'G') & (u != 'T');
  }
}

#ifdef __AVX2__
static inline void scan_vec(const unsigned char *s, unsigned char *f, unsigned char *var, unsigned char *np)
{
  const __m256i ka = _mm256_set1_epi8('a'), k25 = _mm256_set1_epi8(25), k32 = _mm256_set1_epi8(32);
  const __m256i kN = _mm256_set1_epi8('N'), kD = _mm256_set1_epi8('-'), kQ = _mm256_set1_epi8('?');
  const __m256i ones = _mm256_set1_epi8(-1);
  __m256i c = _mm256_loadu_si256((const __m256i *) s);
  __m256i t = _mm256_sub_epi8(c, ka);
  __m256i lo = _mm256_cmpeq_epi8(_mm256_min_epu8(t, k25), t);
  __m256i u = _mm256_sub_epi8(c, _mm256_and_si256(lo, k32));
  __m256i unkm = _mm256_or_si256(_mm256_or_si256(_mm256_cmpeq_epi8(u, kN), _mm256_cmpeq_epi8(u, kD)),
                                 _mm256_cmpeq_epi8(u, kQ));
  __m256i fj = _mm256_loadu_si256((const __m256i *) f);
  __m256i set = _mm256_andnot_si256(unkm, _mm256_cmpeq_epi8(fj, kN));
  fj = _mm256_blendv_epi8(fj, u, set);
  _mm256_storeu_si256((__m256i *) f, fj);
  /* var |= known & (u != f)  ==  ~(unk | (u == f)) */
  __m256i same = _mm256_or_si256(unkm, _mm256_cmpeq_epi8(u, fj));
  __m256i v = _mm256_loadu_si256((const __m256i *) var);
  _mm256_storeu_si256((__m256i *) var, _mm256_or_si256(v, _mm256_andnot_si256(same, ones)));
  if (np) {
    __m256i pure = _mm256_or_si256(_mm256_or_si256(_mm256_cmpeq_epi8(u, _mm256_set1_epi8('A')), _mm256_cmpeq_epi8(u, _mm256_set1_epi8('C'))),
                                   _mm256_or_si256(_mm256_cmpeq_epi8(u, _mm256_set1_epi8('G')), _mm256_cmpeq_epi8(u, _mm256_set1_epi8('T'))));
    __m256i q = _mm256_loadu_si256((const __m256i *) np);
    _mm256_storeu_si256((__m256i *) np, _mm256_or_si256(q, _mm256_andnot_si256(pure, ones)));
  }
}
#endif

/* -c -b without -r, -v or -C: the output columns are exactly those without a non-ACGT byte, so only np is needed
   (f stays 'N' and var 0, which nothing reads in that mode). One OR-accumulate per 32 bytes: (c | 0x20) is one of
   "acgt" exactly for the bytes A C G T a c g t. */
static int g_np_only;
#ifdef __AVX2__
static inline void scan_vec_np(const unsigned char *s, unsigned char *np)
{
  const __m256i k20 = _mm256_set1_epi8(0x20);
  __m256i c = _mm256_or_si256(_mm256_loadu_si256((const __m256i *) s), k20);
  __m256i pure = _mm256_or_si256(_mm256_or_si256(_mm256_cmpeq_epi8(c, _mm256_set1_epi8('a')), _mm256_cmpeq_epi8(c, _mm256_set1_epi8('c'))),
                                 _mm256_or_si256(_mm256_cmpeq_epi8(c, _mm256_set1_epi8('g')), _mm256_cmpeq_epi8(c, _mm256_set1_epi8('t'))));
  __m256i q = _mm256_loadu_si256((const __m256i *) np);
  _mm256_storeu_si256((__m256i *) np, _mm256_or_si256(q, _mm256_andnot_si256(pure, _mm256_set1_epi8(-1))));
}
#endif

/* Same update, and also report whether any byte in s[0..len) is '\n' (fused layout validation). */
static int scan_span_nl(const unsigned char *s, size_t len, unsigned char *f, unsigned char *var, unsigned char *np)
{
  int bad = 0;
#ifdef __AVX2__
  if (len >= 32 && g_np_only && np) {
    const __m256i knl = _mm256_set1_epi8('\n');
    __m256i acc = _mm256_setzero_si256();
    size_t j = 0;
    for (; j + 32 <= len; j += 32) {
      acc = _mm256_or_si256(acc, _mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i *) (s + j)), knl));
      scan_vec_np(s + j, np + j);
    }
    if (j < len) {
      j = len - 32;
      acc = _mm256_or_si256(acc, _mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i *) (s + j)), knl));
      scan_vec_np(s + j, np + j);
    }
    return !_mm256_testz_si256(acc, acc);
  }
  if (len >= 32) {
    const __m256i knl = _mm256_set1_epi8('\n');
    __m256i acc = _mm256_setzero_si256();
    size_t j = 0;
    for (; j + 32 <= len; j += 32) {
      acc = _mm256_or_si256(acc, _mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i *) (s + j)), knl));
      scan_vec(s + j, f + j, var + j, np ? np + j : NULL);
    }
    if (j < len) {
      j = len - 32;
      acc = _mm256_or_si256(acc, _mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i *) (s + j)), knl));
      scan_vec(s + j, f + j, var + j, np ? np + j : NULL);
    }
    return !_mm256_testz_si256(acc, acc);
  }
#endif
  for (size_t j = 0; j < len; j++) bad |= s[j] == '\n';
  scan_scalar(s, len, f, var, np);
  return bad;
}

/* Update columns [0,len) of f/var/np from bytes s[0..len). The update is idempotent for a repeated (column, byte)
   pair, so the last vector may overlap the previous one instead of running a scalar tail. */
static void scan_span(const unsigned char *s, size_t len, unsigned char *f, unsigned char *var, unsigned char *np)
{
#ifdef __AVX2__
  if (len >= 32 && g_np_only && np) {
    size_t j = 0;
    for (; j + 32 <= len; j += 32) scan_vec_np(s + j, np + j);
    if (j < len) { j = len - 32; scan_vec_np(s + j, np + j); }
    return;
  }
  if (len >= 32) {
    size_t j = 0;
    for (; j + 32 <= len; j += 32) scan_vec(s + j, f + j, var + j, np ? np + j : NULL);
    if (j < len) { j = len - 32; scan_vec(s + j, f + j, var + j, np ? np + j : NULL); }
    return;
  }
#endif
  scan_scalar(s, len, f, var, np);
}

static void scan_block(const unsigned char *s, size_t j0, size_t j1, unsigned char *f, unsigned char *var,
                       unsigned char *np)
{
  scan_span(s + j0, j1 - j0, f + j0, var + j0, np ? np + j0 : NULL);
}

/* Copy len bytes (a line of a wrapped record) with 32-byte vector moves; the last move may overlap. */
static inline void copy_small(unsigned char *d, const unsigned char *s, size_t len)
{
#ifdef __AVX2__
  if (len >= 32) {
    size_t j = 0;
    for (; j + 32 <= len; j += 32) _mm256_storeu_si256((__m256i *) (d + j), _mm256_loadu_si256((const __m256i *) (s + j)));
    if (j < len) _mm256_storeu_si256((__m256i *) (d + len - 32), _mm256_loadu_si256((const __m256i *) (s + len - 32)));
    return;
  }
#endif
  memcpy(d, s, len);
}

/* In-place ASCII upper-casing (toupper in the C locale), AVX2. */
static void upper_inplace(unsigned char *d, size_t n)
{
  size_t j = 0;
#ifdef __AVX2__
  const __m256i ka = _mm256_set1_epi8('a'), k25 = _mm256_set1_epi8(25), k32 = _mm256_set1_epi8(32);
  for (; j + 32 <= n; j += 32) {
    __m256i c = _mm256_loadu_si256((const __m256i *) (d + j));
    __m256i t = _mm256_sub_epi8(c, ka);
    __m256i lo = _mm256_cmpeq_epi8(_mm256_min_epu8(t, k25), t);
    _mm256_storeu_si256((__m256i *) (d + j), _mm256_sub_epi8(c, _mm256_and_si256(lo, k32)));
  }
#endif
  for (; j < n; j++) d[j] = up(d[j]);
}

/* Exact-parser path: contiguous records. Column blocks are independent, so they are split across threads. */
static void scan_columns(const rec_t *R, size_t nr, size_t L, int pure, unsigned char *f, unsigned char *var,
                         unsigned char *np)
{
  memset(f, 'N', L); memset(var, 0, L); if (pure) memset(np, 0, L);
  const size_t B = 8192;
  long nb = (long) ((L + B - 1) / B);
#pragma omp parallel for schedule(static)
  for (long bi = 0; bi < nb; bi++) {
    size_t j0 = (size_t) bi * B, j1 = j0 + B < L ? j0 + B : L;
    for (size_t r = 0; r < nr; r++) scan_block(R[r].seq, j0, j1, f, var, pure ? np : NULL);
  }
}

/* Merge column states of a later record group (f1, var1, np1) into an earlier one (f, var, np). Exact: the first
   known base of the union is the earlier group's if it has one; the union is variable if either group is, or if both
   have a first known base and they differ; np is a plain OR. */
static void merge_states(unsigned char *f, unsigned char *var, unsigned char *np, const unsigned char *f1,
                         const unsigned char *var1, const unsigned char *np1, size_t L)
{
  size_t j = 0;
#ifdef __AVX2__
  const __m256i kN = _mm256_set1_epi8('N'), ones = _mm256_set1_epi8(-1);
  for (; j + 32 <= L; j += 32) {
    __m256i a = _mm256_loadu_si256((const __m256i *) (f + j)), b = _mm256_loadu_si256((const __m256i *) (f1 + j));
    __m256i an = _mm256_cmpeq_epi8(a, kN), bn = _mm256_cmpeq_epi8(b, kN);
    __m256i both = _mm256_andnot_si256(_mm256_or_si256(an, bn), ones);
    __m256i differ = _mm256_andnot_si256(_mm256_cmpeq_epi8(a, b), both);
    __m256i v = _mm256_or_si256(_mm256_or_si256(_mm256_loadu_si256((const __m256i *) (var + j)),
                                                _mm256_loadu_si256((const __m256i *) (var1 + j))), differ);
    _mm256_storeu_si256((__m256i *) (var + j), v);
    _mm256_storeu_si256((__m256i *) (f + j), _mm256_blendv_epi8(a, b, an));
    if (np) _mm256_storeu_si256((__m256i *) (np + j), _mm256_or_si256(_mm256_loadu_si256((const __m256i *) (np + j)),
                                                                       _mm256_loadu_si256((const __m256i *) (np1 + j))));
  }
#endif
  for (; j < L; j++) {
    unsigned char a = f[j], b = f1[j];
    var[j] = (unsigned char) (var[j] | var1[j] | ((a != 'N' && b != 'N' && a != b) ? 0xFF : 0));
    if (a == 'N') f[j] = b;
    if (np) np[j] |= np1[j];
  }
}

/* Short alignments (state for all columns fits in L2): each thread scans a contiguous range of records front to
   back (long sequential reads instead of thousands of short per-block streams) into private state, and the group
   states are merged in record order. Same validation as scan_columns_layout. Returns 0 if the layout was wrong. */
static int scan_records_layout(const rec_t *R, size_t nr, size_t L, size_t W, int pure, unsigned char *f,
                               unsigned char *var, unsigned char *np)
{
  int T = 1;
#ifdef _OPENMP
  T = omp_get_max_threads();
#endif
  if ((size_t) T > nr) T = (int) (nr ? nr : 1);
  unsigned char **F = calloc((size_t) T, sizeof(unsigned char *)), **V = calloc((size_t) T, sizeof(unsigned char *)),
                **P = calloc((size_t) T, sizeof(unsigned char *));
  int *bad = calloc((size_t) T, sizeof(int));
  F[0] = f; V[0] = var; P[0] = pure ? np : NULL;
#pragma omp parallel num_threads(T)
  {
    int tid = 0, nt = 1;
#ifdef _OPENMP
    tid = omp_get_thread_num(); nt = omp_get_num_threads();
#endif
    for (int g = tid; g < T; g += nt) {
      if (g > 0) { F[g] = malloc(L + 64); V[g] = malloc(L + 64); P[g] = pure ? malloc(L + 64) : NULL; }
      memset(F[g], 'N', L); memset(V[g], 0, L); if (pure) memset(P[g], 0, L);
      size_t r0 = nr * (size_t) g / (size_t) T, r1 = nr * (size_t) (g + 1) / (size_t) T;
      int b = 0;
      if (W >= L) {
        for (size_t r = r0; r < r1; r++) {
          const unsigned char *q = R[r].seq;
          b |= scan_span_nl(q, L, F[g], V[g], P[g]);
          b |= q[0] == '>' || q[0] == '+' || q[0] == '@' || q[L - 1] == '\r' || q[L] != '\n';
        }
      } else {
        size_t lpb = 8192 / W ? 8192 / W : 1, B = lpb * W;
        unsigned char *tmp = aligned_alloc(64, (B + 127) / 64 * 64);
        for (size_t r = r0; r < r1; r++) {
          const unsigned char *src = R[r].seq;
          for (size_t j0 = 0; j0 < L; j0 += B) {
            size_t j1 = j0 + B < L ? j0 + B : L;
            unsigned char *d = tmp;
            for (size_t j = j0; j < j1; j += W, src += W + 1) {
              size_t len = j1 - j < W ? j1 - j : W;
              copy_small(d, src, len);
              d += len;
              b |= src[0] == '>' || src[0] == '+' || src[0] == '@' || src[len - 1] == '\r' || src[len] != '\n';
            }
            b |= scan_span_nl(tmp, j1 - j0, F[g] + j0, V[g] + j0, P[g] ? P[g] + j0 : NULL);
          }
        }
        free(tmp);
      }
      bad[g] = b;
    }
  }
  int anybad = 0;
  for (int g = 0; g < T; g++) anybad |= bad[g];
  for (int g = 1; g < T; g++) {
    merge_states(f, var, pure ? np : NULL, F[g], V[g], P[g], L);
    free(F[g]); free(V[g]); free(P[g]);
  }
  free(F); free(V); free(P); free(bad);
  return !anybad;
}

/* Fast-layout path: records contiguous (W >= L) or wrapped at width W (column j at seq + j + j/W). Also validates
   the predicted layout (see parse_fast): no '\n' inside any line, '\n' at every line end, no line starting with
   '>', '+' or '@', no line ending in '\r'. Returns 0 if the layout was wrong (the caller then re-parses exactly). */
static int scan_columns_layout(const rec_t *R, size_t nr, size_t L, size_t W, int pure, unsigned char *f,
                               unsigned char *var, unsigned char *np)
{
  /* column states are initialised block by block inside the parallel loops */
  int bad = 0;
  if (W >= L) {
    const size_t B = 8192;
    long nb = (long) ((L + B - 1) / B);
#pragma omp parallel for schedule(static) reduction(|:bad)
    for (long bi = 0; bi < nb; bi++) {
      size_t j0 = (size_t) bi * B, j1 = j0 + B < L ? j0 + B : L;
      memset(f + j0, 'N', j1 - j0); memset(var + j0, 0, j1 - j0); if (pure) memset(np + j0, 0, j1 - j0);
      for (size_t r = 0; r < nr; r++)
        bad |= scan_span_nl(R[r].seq + j0, j1 - j0, f + j0, var + j0, pure ? np + j0 : NULL);
    }
    for (size_t r = 0; r < nr; r++) {
      const unsigned char *q = R[r].seq;
      bad |= q[0] == '>' || q[0] == '+' || q[0] == '@' || q[L - 1] == '\r' || q[L] != '\n';
    }
    return !bad;
  }
  /* wrapped: copy each block's lines into a contiguous L1 buffer, then run the contiguous kernel on it */
  size_t lpb = 8192 / W ? 8192 / W : 1, B = lpb * W;
  long nb = (long) ((L + B - 1) / B);
#pragma omp parallel reduction(|:bad)
  {
    unsigned char *tmp = aligned_alloc(64, (B + 127) / 64 * 64);
#pragma omp for schedule(static)
    for (long bi = 0; bi < nb; bi++) {
      size_t j0 = (size_t) bi * B, j1 = j0 + B < L ? j0 + B : L;
      size_t line0 = j0 / W;
      memset(f + j0, 'N', j1 - j0); memset(var + j0, 0, j1 - j0); if (pure) memset(np + j0, 0, j1 - j0);
      for (size_t r = 0; r < nr; r++) {
        const unsigned char *src = R[r].seq + line0 * (W + 1);
        unsigned char *d = tmp;
        for (size_t j = j0; j < j1; j += W, src += W + 1) {
          size_t len = j1 - j < W ? j1 - j : W;
          copy_small(d, src, len);
          d += len;
          bad |= src[0] == '>' || src[0] == '+' || src[0] == '@' || src[len - 1] == '\r' || src[len] != '\n';
        }
        bad |= scan_span_nl(tmp, j1 - j0, f + j0, var + j0, pure ? np + j0 : NULL);
      }
    }
    free(tmp);
  }
  return !bad;
}

/* -C with unequal sequence lengths: snp-sites' first pass reads the first L bytes of kseq's reused sequence buffer
   for every record. For a shorter record those bytes are its sequence, its NUL terminator, then whatever earlier
   records left in the buffer. All positions < L were written by record 1, so this is deterministic; emulate it. */
static void scan_columns_shadow(const rec_t *R, size_t nr, size_t L, unsigned char *f, unsigned char *var)
{
  unsigned char *shadow = malloc(L + 1), *eff = malloc(L + 1);
  memcpy(shadow, R[0].seq, L); shadow[L] = 0;
  memset(f, 'N', L); memset(var, 0, L);
  for (size_t r = 0; r < nr; r++) {
    size_t l = R[r].len;
    const unsigned char *s;
    if (l >= L) s = R[r].seq;
    else {
      memcpy(eff, R[r].seq, l); eff[l] = 0;
      if (l + 1 < L) memcpy(eff + l + 1, shadow + l + 1, L - l - 1);
      s = eff;
    }
    scan_block(s, 0, L, f, var, NULL);
    memcpy(shadow, R[r].seq, l < L ? l : L);
    if (l < L) shadow[l] = 0;
  }
  free(shadow); free(eff);
}

/* ------------------------------------------------------------------ outputs */

typedef struct {
  const rec_t *R; size_t nr;
  const int *pos; size_t S;          /* output columns */
  const uint32_t *off;               /* byte offset of each output column inside a record's sequence bytes */
  const uint32_t *seg_off, *seg_len; size_t nseg; int use_segs; /* contiguous copy segments for dense gathers */
  /* 16-byte windows for medium-density gathers: load 16 raw bytes at win_src[w], pshufb with win_shuf[w] packs the
     window's output columns to the front, store, advance by win_cnt[w] */
  const uint32_t *win_src; const uint8_t *win_cnt; const __m128i *win_shuf; size_t nwin; int use_wins;
  const unsigned char *f; size_t L;
} snps_t;

static size_t name_len_c(const rec_t *r) { return strnlen((const char *) r->name, r->name_len); }

static void gather_row(const snps_t *z, const unsigned char *s, unsigned char *dst)
{
  if (z->use_wins) { /* may write up to 15 bytes past dst + S; callers leave slack and overwrite it */
    const __m128i ka = _mm_set1_epi8('a'), k25 = _mm_set1_epi8(25), k32 = _mm_set1_epi8(32);
    unsigned char *d = dst;
    for (size_t w = 0; w < z->nwin; w++) {
      __m128i v = _mm_shuffle_epi8(_mm_loadu_si128((const __m128i *) (s + z->win_src[w])), z->win_shuf[w]);
      __m128i t = _mm_sub_epi8(v, ka);
      __m128i lo = _mm_cmpeq_epi8(_mm_min_epu8(t, k25), t);
      _mm_storeu_si128((__m128i *) d, _mm_sub_epi8(v, _mm_and_si128(lo, k32)));
      d += z->win_cnt[w];
    }
    return;
  }
  if (z->use_segs) {
    unsigned char *d = dst;
    for (size_t k = 0; k < z->nseg; k++) { memcpy(d, s + z->seg_off[k], z->seg_len[k]); d += z->seg_len[k]; }
  } else {
    const uint32_t *off = z->off;
    for (size_t j = 0; j < z->S; j++) dst[j] = s[off[j]];
  }
  upper_inplace(dst, z->S);
}

static unsigned char *thp_alloc(size_t n)
{
  void *m = mmap(NULL, n, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
  if (m == MAP_FAILED) die("ERROR: out of memory\n");
  madvise(m, n, MADV_HUGEPAGE);
  return m;
}

/* Large never-freed working arrays (column states, output column tables): transparent huge pages, so first-touch
   faults are 2 MB at a time. */
static unsigned char *big_alloc(size_t n)
{
  if (n < (2u << 20)) { unsigned char *p = malloc(n); if (!p) die("ERROR: out of memory\n"); return p; }
  return thp_alloc(n);
}

static void write_all(int fd, const char *p, size_t n)
{
  while (n) {
    ssize_t k = write(fd, p, n);
    if (k <= 0) { perror("write"); exit(EXIT_FAILURE); }
    p += k; n -= (size_t) k;
  }
}

/* writev with partial-write handling, at most 1024 segments per call */
static void writev_all(int fd, struct iovec *iov, size_t cnt)
{
  while (cnt > 0) {
    int c = cnt < 1024 ? (int) cnt : 1024;
    ssize_t k = writev(fd, iov, c);
    if (k <= 0) { perror("write"); exit(EXIT_FAILURE); }
    while (cnt > 0 && (size_t) k >= iov->iov_len) { k -= (ssize_t) iov->iov_len; iov++; cnt--; }
    if (cnt > 0 && k > 0) { iov->iov_base = (char *) iov->iov_base + k; iov->iov_len -= (size_t) k; }
  }
}

static int open_out(const char *name)
{
  int fd = open(name, O_WRONLY | O_CREAT | O_TRUNC, 0666);
  if (fd < 0) { fprintf(stderr, "ERROR: cannot open %s for writing: %s\n", name, strerror(errno)); exit(EXIT_FAILURE); }
  return fd;
}

/* Output files of known approximate size are preallocated (fallocate, FALLOC_FL_KEEP_SIZE: the file size still
   grows only with the bytes written, so contents are unaffected); ext4 then skips per-page delayed-allocation work
   in write(), measured 25% faster for a 361 MB file. Only regular files; any failure is ignored. At close the file
   is truncated to its written length, which releases preallocated blocks beyond the end. */
static int prealloc_fds[8], n_prealloc;
static void prealloc(int fd, size_t bytes)
{
  struct stat st;
  if (bytes < (1u << 20) || n_prealloc >= 8 || fstat(fd, &st) != 0 || !S_ISREG(st.st_mode)) return;
  if (fallocate(fd, FALLOC_FL_KEEP_SIZE, 0, (off_t) bytes) == 0) prealloc_fds[n_prealloc++] = fd;
}
static void close_out(int fd)
{
  for (int i = 0; i < n_prealloc; i++)
    if (prealloc_fds[i] == fd) {
      off_t e = lseek(fd, 0, SEEK_CUR);
      if (e >= 0 && ftruncate(fd, e) != 0) { /* blocks beyond the end stay allocated; contents unaffected */ }
      prealloc_fds[i] = -1;
    }
  close(fd);
}

/* Output pipeline: the file is cut into chunks of a few MB; threads format chunks into their own buffers and write
   them in chunk order (OpenMP ordered), so one thread's write overlaps the others' formatting. Chunking does not
   depend on the thread count, and the bytes are the same for any thread count. */
typedef size_t (*fmt_fn)(const void *ctx, size_t chunk, char *buf, unsigned char *scratch);

static void run_chunked(int fd, size_t nchunks, size_t maxbytes, size_t scratch_bytes, fmt_fn fmt, const void *ctx)
{
#pragma omp parallel
  {
    char *buf = malloc(maxbytes + 64);
    unsigned char *scratch = scratch_bytes ? malloc(scratch_bytes + 64) : NULL;
#pragma omp for ordered schedule(dynamic, 1)
    for (long c = 0; c < (long) nchunks; c++) {
      size_t len = fmt(ctx, (size_t) c, buf, scratch);
#pragma omp ordered
      write_all(fd, buf, len);
    }
    free(buf); free(scratch);
  }
}

/* Several output files at once: the chunks of all files form one task list, interleaved by relative progress, so
   the two threads usually fill different files (writes to different files run in parallel in the kernel; writes
   to one file are serialised by its inode lock). A thread formats its chunk, waits for its file's turn (the
   previous chunk of that file written), writes, and passes the turn on. Tasks are handed out in list order and each
   file's chunks appear in order, so the lowest unfinished task can always proceed. */
/* A job fills one file, or two files (fd2 >= 0) whose chunks come from one formatting pass (fmt2 writes the second
   file's chunk into buf2): FASTA and PHYLIP rows share one gather per record. */
typedef size_t (*fmt2_fn)(const void *ctx, size_t chunk, char *buf, char *buf2, size_t *len2, unsigned char *scratch);
/* iov2: fmt2 leaves the second file's chunk as an iovec list in scratch (count in *len2), pointing at its own
   prefixes in buf2 and at row bytes already formatted in buf, so long rows are not copied twice */
typedef struct { int fd; size_t nchunks, maxbytes, scratch_bytes; fmt_fn fmt; const void *ctx; int fd2; fmt2_fn fmt2; int iov2; } job_t;

static inline void wait_turn(volatile size_t *turn, size_t c)
{
  int spins = 0;
  while (__atomic_load_n(turn, __ATOMIC_ACQUIRE) != c) { if (++spins < 64) _mm_pause(); else sched_yield(); }
}
static inline void pass_turn(volatile size_t *turn, size_t c) { __atomic_store_n(turn, c + 1, __ATOMIC_RELEASE); }

static void run_jobs(job_t *jobs, int njobs)
{
  if (njobs == 1 && jobs[0].fd2 < 0) { run_chunked(jobs[0].fd, jobs[0].nchunks, jobs[0].maxbytes, jobs[0].scratch_bytes, jobs[0].fmt, jobs[0].ctx); return; }
  size_t ntasks = 0, maxb = 0, maxs = 0;
  for (int i = 0; i < njobs; i++) {
    ntasks += jobs[i].nchunks;
    if (jobs[i].maxbytes > maxb) maxb = jobs[i].maxbytes;
    if (jobs[i].scratch_bytes > maxs) maxs = jobs[i].scratch_bytes;
  }
  int *tj = malloc((ntasks + 1) * sizeof(int));
  size_t *tc = malloc((ntasks + 1) * sizeof(size_t)), nextc[8] = {0};
  for (size_t t = 0; t < ntasks; t++) {
    int best = -1; double bp = 2.0;
    for (int i = 0; i < njobs; i++)
      if (nextc[i] < jobs[i].nchunks) {
        double pr = (nextc[i] + 0.5) / (double) jobs[i].nchunks;
        if (pr < bp) { bp = pr; best = i; }
      }
    tj[t] = best; tc[t] = nextc[best]++;
  }
  volatile size_t turns[8] = {0}, turns2[8] = {0};
  int any2 = 0;
  for (int i = 0; i < njobs; i++) any2 |= jobs[i].fd2 >= 0;
#pragma omp parallel
  {
    char *buf = malloc(maxb + 64), *buf2 = any2 ? malloc(maxb + 64) : NULL;
    unsigned char *scratch = maxs ? malloc(maxs + 64) : NULL;
#pragma omp for schedule(dynamic, 1)
    for (long t = 0; t < (long) ntasks; t++) {
      job_t *jb = &jobs[tj[t]];
      size_t len2 = 0;
      size_t len = jb->fd2 >= 0 ? jb->fmt2(jb->ctx, tc[t], buf, buf2, &len2, scratch) : jb->fmt(jb->ctx, tc[t], buf, scratch);
      wait_turn(&turns[tj[t]], tc[t]);
      write_all(jb->fd, buf, len);
      pass_turn(&turns[tj[t]], tc[t]);
      if (jb->fd2 >= 0) {
        wait_turn(&turns2[tj[t]], tc[t]);
        if (jb->iov2) writev_all(jb->fd2, (struct iovec *) scratch, len2);
        else write_all(jb->fd2, buf2, len2);
        pass_turn(&turns2[tj[t]], tc[t]);
      }
    }
    free(buf); free(buf2); free(scratch);
  }
  free(tj); free(tc);
}

#define CHUNK_TARGET (2u << 20)

typedef struct {
  const snps_t *z; int phylip, output_reference;
  const size_t *rstart; /* chunk c covers records [rstart[c], rstart[c+1]) */
  size_t maxrows; int iov;  /* most records in a chunk; PHYLIP chunk as iovecs into the FASTA chunk (fmt_rows2) */
} rows_ctx_t;

static size_t fmt_rows(const void *vctx, size_t c, char *buf, unsigned char *scratch)
{
  (void) scratch;
  const rows_ctx_t *x = vctx; const snps_t *z = x->z;
  char *d = buf;
  if (c == 0) {
    if (x->phylip) d += sprintf(d, "%d %d\n", (int) z->nr + (x->output_reference ? 1 : 0), (int) z->S);
    if (x->output_reference) {
      const char *h = x->phylip ? "pseudo_reference_sequence\t" : ">pseudo_reference_sequence\n";
      size_t hl = strlen(h); memcpy(d, h, hl); d += hl;
      for (size_t j = 0; j < z->S; j++) d[j] = (char) z->f[z->pos[j]];
      d += z->S; *d++ = '\n';
    }
  }
  for (size_t r = x->rstart[c]; r < x->rstart[c + 1]; r++) {
    size_t nl = name_len_c(&z->R[r]);
    if (!x->phylip) *d++ = '>';
    memcpy(d, z->R[r].name, nl); d += nl;
    *d++ = x->phylip ? '\t' : '\n';
    gather_row(z, z->R[r].seq, (unsigned char *) d); d += z->S; *d++ = '\n';
  }
  return (size_t) (d - buf);
}

/* FASTA rows into buf and PHYLIP rows into buf2 from one gather per record */
static size_t fmt_rows2(const void *vctx, size_t c, char *buf, char *buf2, size_t *len2, unsigned char *scratch)
{
  const rows_ctx_t *x = vctx; const snps_t *z = x->z;
  char *d = buf, *p = buf2;
  if (x->iov) {
    struct iovec *v = (struct iovec *) scratch; size_t nv = 0;
    if (c == 0) {
      char *p0 = p;
      p += sprintf(p, "%d %d\n", (int) z->nr + (x->output_reference ? 1 : 0), (int) z->S);
      if (x->output_reference) {
        memcpy(p, "pseudo_reference_sequence\t", 26); p += 26;
        memcpy(d, ">pseudo_reference_sequence\n", 27); d += 27;
        v[nv].iov_base = p0; v[nv].iov_len = (size_t) (p - p0); nv++;
        v[nv].iov_base = d; v[nv].iov_len = z->S + 1; nv++;
        for (size_t j = 0; j < z->S; j++) d[j] = (char) z->f[z->pos[j]];
        d += z->S; *d++ = '\n';
      } else { v[nv].iov_base = p0; v[nv].iov_len = (size_t) (p - p0); nv++; }
    }
    for (size_t r = x->rstart[c]; r < x->rstart[c + 1]; r++) {
      size_t nl = name_len_c(&z->R[r]);
      *d++ = '>'; memcpy(d, z->R[r].name, nl); d += nl; *d++ = '\n';
      v[nv].iov_base = p; memcpy(p, z->R[r].name, nl); p += nl; *p++ = '\t'; v[nv].iov_len = nl + 1; nv++;
      gather_row(z, z->R[r].seq, (unsigned char *) d);
      v[nv].iov_base = d; v[nv].iov_len = z->S + 1; nv++;
      d += z->S; *d++ = '\n';
    }
    *len2 = nv;
    return (size_t) (d - buf);
  }
  if (c == 0) {
    p += sprintf(p, "%d %d\n", (int) z->nr + (x->output_reference ? 1 : 0), (int) z->S);
    if (x->output_reference) {
      memcpy(d, ">pseudo_reference_sequence\n", 27); d += 27;
      memcpy(p, "pseudo_reference_sequence\t", 26); p += 26;
      for (size_t j = 0; j < z->S; j++) d[j] = p[j] = (char) z->f[z->pos[j]];
      d += z->S; *d++ = '\n'; p += z->S; *p++ = '\n';
    }
  }
  for (size_t r = x->rstart[c]; r < x->rstart[c + 1]; r++) {
    size_t nl = name_len_c(&z->R[r]);
    *d++ = '>'; memcpy(d, z->R[r].name, nl); d += nl; *d++ = '\n';
    memcpy(p, z->R[r].name, nl); p += nl; *p++ = '\t';
    gather_row(z, z->R[r].seq, (unsigned char *) d);
    memcpy(p, d, z->S);
    d += z->S; *d++ = '\n'; p += z->S; *p++ = '\n';
  }
  *len2 = (size_t) (p - buf2);
  return (size_t) (d - buf);
}

static size_t rows_bytes(const snps_t *z)
{
  size_t b = z->S + 64;
  for (size_t r = 0; r < z->nr; r++) b += name_len_c(&z->R[r]) + z->S + 3;
  return b;
}
static size_t vcf_bytes(const snps_t *z)
{
  size_t b = 256 + z->S * (2 * z->nr + 48);
  for (size_t r = 0; r < z->nr; r++) b += name_len_c(&z->R[r]) + 1;
  return b;
}

static job_t prepare_rows(int fd, const snps_t *z, int output_reference, int phylip)
{
  size_t *rstart = malloc((z->nr + 2) * sizeof(size_t));
  size_t nchunks = 0, acc = 0, maxb = 0, pre = z->S + 128;
  rstart[0] = 0;
  for (size_t r = 0; r < z->nr; r++) {
    size_t rb = name_len_c(&z->R[r]) + z->S + 3; /* ">" name "\n" row "\n" */
    if (acc && acc + rb > CHUNK_TARGET) { if (acc > maxb) maxb = acc; rstart[++nchunks] = r; acc = 0; }
    acc += rb + (nchunks == 0 && acc == 0 ? pre : 0);
  }
  if (acc > maxb) maxb = acc;
  rstart[++nchunks] = z->nr;
  if (z->nr == 0) { nchunks = 1; maxb = pre; }
  rows_ctx_t *x = malloc(sizeof *x);
  x->z = z; x->phylip = phylip; x->output_reference = output_reference; x->rstart = rstart;
  x->maxrows = 0; x->iov = 0;
  for (size_t q = 0; q < nchunks; q++) if (rstart[q + 1] - rstart[q] > x->maxrows) x->maxrows = rstart[q + 1] - rstart[q];
  job_t j = { fd, nchunks, maxb + pre, 0, fmt_rows, x, -1, NULL };
  return j;
}

/* One VCF data row. b[] holds the column's upper-cased bases in sample order. ALT symbols are the distinct values
   != REF in order of first appearance, unknowns (N n - ?) written as '*'. A vector pass finds first appearances:
   each 32-byte block is compared against the set of raw byte values already seen (REF, the ALT bytes and, once '*'
   is in ALT, all four unknown bytes); only blocks with an unseen byte drop to scalar code, lowest lane first. */
static char *vcf_row(char *o, const unsigned char *b, int snp_location, size_t n, unsigned char ref)
{
  if (ref == '\0') { fprintf(stderr, "Couldnt get the reference base for coordinate %d\n", snp_location); exit(EXIT_FAILURE); }
  /* per-thread lookup tables, all-zero between rows: only the entries a row sets are cleared at its end */
  static __thread unsigned char altidx[256], known[256];
  unsigned char kv[64]; int nk = 0;
  unsigned char alts[32]; int na = 0;
  kv[nk++] = ref; known[ref] = 1;
  /* 4-bit hash h(v) = (v >> 1) & 15 is collision-free on A C G T N - ? (and R Y W K D H): while every raw value
     present in the row has its own slot, one pshufb replaces the per-value compare chains. b[] is upper-cased, so
     'n' never occurs and is not hashed. Empty slots hold a sentinel whose own hash differs from the slot. */
  unsigned char hslot[16]; int hash_ok = 1; unsigned hused = 0;
  for (int k = 0; k < 16; k++) hslot[k] = (unsigned char) (((k + 1) & 15) << 1);
#define HASH_ADD(cc) do { unsigned char h_ = (unsigned char) (((cc) >> 1) & 15); \
    if (hused & (1u << h_)) { if (hslot[h_] != (cc)) hash_ok = 0; } else { hused |= 1u << h_; hslot[h_] = (cc); } } while (0)
  HASH_ADD(ref);
#define ADD_RAW(cc) do { unsigned char c_ = (cc); \
    unsigned char m_ = unk(c_) ? '*' : c_; \
    if (!altidx[m_]) { \
      if (na >= MAXIMUM_NUMBER_OF_ALT_BASES) \
        die("Unexpectedly large number of alternative bases found between sequences.  Please check input file is not corrupted\n\n"); \
      alts[na++] = m_; altidx[m_] = (unsigned char) na; } \
    if (unk(c_)) { static const unsigned char U4[4] = {'N', 'n', '-', '?'}; \
      for (int q_ = 0; q_ < 4; q_++) if (!known[U4[q_]]) { known[U4[q_]] = 1; kv[nk++] = U4[q_]; if (U4[q_] != 'n') HASH_ADD(U4[q_]); } } \
    else if (!known[c_]) { known[c_] = 1; kv[nk++] = c_; HASH_ADD(c_); } } while (0)
  size_t i = 0;
#ifdef __AVX2__
  const __m256i k0f = _mm256_set1_epi8(0x0F);
  while (hash_ok && i + 32 <= n) {
    __m256i v = _mm256_loadu_si256((const __m256i *) (b + i));
    __m256i h = _mm256_and_si256(_mm256_srli_epi16(v, 1), k0f);
    __m256i sv = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *) hslot));
    unsigned mask = ~(unsigned) _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_shuffle_epi8(sv, h), v));
    if (!mask) { i += 32; continue; }
    ADD_RAW(b[i + __builtin_ctz(mask)]); /* lowest unseen lane first; then re-check this vector */
  }
  for (; i + 32 <= n; i += 32) {
    __m256i v = _mm256_loadu_si256((const __m256i *) (b + i));
    __m256i m = _mm256_cmpeq_epi8(v, _mm256_set1_epi8((char) kv[0]));
    for (int k = 1; k < nk; k++) m = _mm256_or_si256(m, _mm256_cmpeq_epi8(v, _mm256_set1_epi8((char) kv[k])));
    unsigned mask = ~(unsigned) _mm256_movemask_epi8(m);
    while (mask) {
      int lane = __builtin_ctz(mask);
      int before = nk;
      ADD_RAW(b[i + lane]);
      for (int k = before; k < nk; k++) m = _mm256_or_si256(m, _mm256_cmpeq_epi8(v, _mm256_set1_epi8((char) kv[k])));
      mask = ~(unsigned) _mm256_movemask_epi8(m);
    }
  }
#endif
  for (; i < n; i++) if (!known[b[i]]) ADD_RAW(b[i]);
#undef ADD_RAW
#undef HASH_ADD
  if (na >= MAXIMUM_NUMBER_OF_ALT_BASES) abort(); /* snp-sites: assert(number_of_alt_bases < MAX) */
  { /* "1\t<POS>\t.\t<REF>\t" without sprintf */
    char dg[16]; int nd = 0; unsigned v = (unsigned) snp_location + 1;
    do { dg[nd++] = (char) ('0' + v % 10); v /= 10; } while (v);
    *o++ = '1'; *o++ = '\t';
    while (nd) *o++ = dg[--nd];
    *o++ = '\t'; *o++ = '.'; *o++ = '\t'; *o++ = (char) ref; *o++ = '\t';
  }
  if (na == 0) *o++ = '.';
  else { *o++ = (char) alts[0]; for (int k = 1; k < na; k++) { *o++ = ','; *o++ = (char) alts[k]; } }
  memcpy(o, "\t.\t.\t.\tGT\t", 10); o += 10;
  if (na >= 10) { /* 2-digit allele indices: scalar */
    for (size_t i2 = 0; i2 < n; i2++) {
      unsigned char m = unk(b[i2]) ? '*' : b[i2];
      if (m == ref) *o++ = '0';
      else if (altidx[m]) { int v = altidx[m]; if (v >= 10) *o++ = (char) ('0' + v / 10); *o++ = (char) ('0' + v % 10); }
      else *o++ = '.';
      *o++ = i2 + 1 != n ? '\t' : '\n';
    }
    for (int k = 0; k < na; k++) altidx[alts[k]] = 0;
    for (int k = 0; k < nk; k++) known[kv[k]] = 0;
    return o;
  }
  /* single-character genotype code for each raw value present */
  unsigned char code[64] = {0};
  for (int k = 0; k < nk; k++) {
    unsigned char m = unk(kv[k]) ? '*' : kv[k];
    code[k] = (m == ref) ? '0' : (altidx[m] ? (unsigned char) ('0' + altidx[m]) : '.');
  }
  size_t need = n * 2;
  unsigned char *p = (unsigned char *) o;
  size_t i2 = 0;
#ifdef __AVX2__
  if (hash_ok) { /* one pshufb per 32 samples: code table indexed by h(v) */
    unsigned char ctab[16] = {0};
    for (int k = 0; k < nk; k++) if (kv[k] != 'n') ctab[(kv[k] >> 1) & 15] = code[k];
    const __m256i tab = _mm256_set1_epi8('\t'), k0f2 = _mm256_set1_epi8(0x0F);
    const __m256i cv = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *) ctab));
    for (; i2 + 32 <= n; i2 += 32) {
      __m256i v = _mm256_loadu_si256((const __m256i *) (b + i2));
      __m256i out = _mm256_shuffle_epi8(cv, _mm256_and_si256(_mm256_srli_epi16(v, 1), k0f2));
      out = _mm256_permute4x64_epi64(out, 0xD8);
      _mm256_storeu_si256((__m256i *) (p + 2 * i2), _mm256_unpacklo_epi8(out, tab));
      _mm256_storeu_si256((__m256i *) (p + 2 * i2 + 32), _mm256_unpackhi_epi8(out, tab));
    }
  } else if (nk <= 12) {
    const __m256i tab = _mm256_set1_epi8('\t');
    __m256i c0 = _mm256_set1_epi8((char) code[0]);
    for (; i2 + 32 <= n; i2 += 32) {
      __m256i v = _mm256_loadu_si256((const __m256i *) (b + i2));
      __m256i out = c0;
      for (int k = 1; k < nk; k++)
        out = _mm256_blendv_epi8(out, _mm256_set1_epi8((char) code[k]), _mm256_cmpeq_epi8(v, _mm256_set1_epi8((char) kv[k])));
      out = _mm256_permute4x64_epi64(out, 0xD8);
      _mm256_storeu_si256((__m256i *) (p + 2 * i2), _mm256_unpacklo_epi8(out, tab));
      _mm256_storeu_si256((__m256i *) (p + 2 * i2 + 32), _mm256_unpackhi_epi8(out, tab));
    }
  }
#endif
  if (i2 < n) {
    unsigned char lut[256];
    for (int k = 0; k < nk; k++) lut[kv[k]] = code[k];
    for (; i2 < n; i2++) { p[2 * i2] = lut[b[i2]]; p[2 * i2 + 1] = '\t'; }
  }
  p[need - 1] = '\n';
  for (int k = 0; k < na; k++) altidx[alts[k]] = 0;
  for (int k = 0; k < nk; k++) known[kv[k]] = 0;
  return o + need;
}

typedef struct { const snps_t *z; size_t snps_per_chunk; } vcf_ctx_t;

/* 16x16 byte transpose: dst row i, byte k = src row k, byte i. SSE2 unpack network (8/16/32/64-bit stages). */
static inline void transpose16(const unsigned char *src, size_t ss, unsigned char *dst, size_t ds)
{
  __m128i a[16], b[16], c4[4][4], c8[2][8];
  for (int i = 0; i < 16; i++) a[i] = _mm_loadu_si128((const __m128i *) (src + i * ss));
  for (int p = 0; p < 8; p++) { b[2 * p] = _mm_unpacklo_epi8(a[2 * p], a[2 * p + 1]); b[2 * p + 1] = _mm_unpackhi_epi8(a[2 * p], a[2 * p + 1]); }
  for (int q = 0; q < 4; q++) {
    c4[q][0] = _mm_unpacklo_epi16(b[4 * q], b[4 * q + 2]);
    c4[q][1] = _mm_unpackhi_epi16(b[4 * q], b[4 * q + 2]);
    c4[q][2] = _mm_unpacklo_epi16(b[4 * q + 1], b[4 * q + 3]);
    c4[q][3] = _mm_unpackhi_epi16(b[4 * q + 1], b[4 * q + 3]);
  }
  for (int h = 0; h < 2; h++)
    for (int g = 0; g < 4; g++) {
      c8[h][2 * g] = _mm_unpacklo_epi32(c4[2 * h][g], c4[2 * h + 1][g]);
      c8[h][2 * g + 1] = _mm_unpackhi_epi32(c4[2 * h][g], c4[2 * h + 1][g]);
    }
  for (int m = 0; m < 8; m++) {
    _mm_storeu_si128((__m128i *) (dst + (size_t) (2 * m) * ds), _mm_unpacklo_epi64(c8[0][m], c8[1][m]));
    _mm_storeu_si128((__m128i *) (dst + (size_t) (2 * m + 1) * ds), _mm_unpackhi_epi64(c8[0][m], c8[1][m]));
  }
}

static size_t fmt_vcf(const void *vctx, size_t c, char *buf, unsigned char *T)
{
  const vcf_ctx_t *x = vctx; const snps_t *z = x->z;
  char *o = buf;
  const size_t TB = 64, n = z->nr;
  if (c == 0) {
    o += sprintf(o, "##fileformat=VCFv4.1\n##contig=<ID=1,length=%i>\n", (int) z->L);
    o += sprintf(o, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    o += sprintf(o, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (size_t r = 0; r < n; r++) { *o++ = '\t'; size_t nl = name_len_c(&z->R[r]); memcpy(o, z->R[r].name, nl); o += nl; }
    *o++ = '\n';
  }
  size_t a = c * x->snps_per_chunk, e = a + x->snps_per_chunk < z->S ? a + x->snps_per_chunk : z->S;
  for (size_t j0 = a; j0 < e; j0 += TB) {
    size_t k = e - j0 < TB ? e - j0 : TB;
    const int *pos = z->pos + j0;
    const uint32_t *off = z->off + j0;
    if (k == TB && off[TB - 1] - off[0] == TB - 1) {
      /* 64 contiguous columns: copy 64 bytes per record into a sample-major tile G, then 16x16 SSE transposes */
      unsigned char *G = T + TB * n + 64;
      size_t o0 = off[0], r0 = 0;
      for (size_t r = 0; r < n; r++) {
        const unsigned char *s = z->R[r].seq + o0;
        _mm256_storeu_si256((__m256i *) (G + r * TB), _mm256_loadu_si256((const __m256i *) s));
        _mm256_storeu_si256((__m256i *) (G + r * TB + 32), _mm256_loadu_si256((const __m256i *) (s + 32)));
      }
      for (; r0 + 16 <= n; r0 += 16)
        for (size_t t0 = 0; t0 < TB; t0 += 16) transpose16(G + r0 * TB + t0, TB, T + t0 * n + r0, n);
      for (; r0 < n; r0++)
        for (size_t t = 0; t < TB; t++) T[t * n + r0] = G[r0 * TB + t];
    } else {
      for (size_t r = 0; r < n; r++) {
        const unsigned char *s = z->R[r].seq;
        for (size_t t = 0; t < k; t++) T[t * n + r] = s[off[t]];
      }
    }
    upper_inplace(T, k * n);
    for (size_t t = 0; t < k; t++) o = vcf_row(o, T + t * n, pos[t], n, z->f[pos[t]]);
  }
  return (size_t) (o - buf);
}

static job_t prepare_vcf(int fd, const snps_t *z)
{
  size_t n = z->nr, row_typ = 2 * n + 40;
  size_t spc = CHUNK_TARGET / row_typ / 64 * 64;
  if (spc < 64) spc = 64;
  size_t nchunks = (z->S + spc - 1) / spc;
  size_t header = 256;
  for (size_t r = 0; r < n; r++) header += name_len_c(&z->R[r]) + 1;
  vcf_ctx_t *x = malloc(sizeof *x);
  x->z = z; x->snps_per_chunk = spc;
  job_t j = { fd, nchunks, header + spc * (3 * n + 160), 128 * n + 256, fmt_vcf, x, -1, NULL };
  return j;
}

/* Large-n path (VCF requested, n >= 512). Strided VCF tile gathers from the 4 KB-page input cost about one TLB
   miss per record per tile, which dominates when n is large. Instead:
   Phase A (parallel over units of 16-record blocks, ordered writes): gather each record's output columns once
     (sequential reads), write its FASTA/PHYLIP rows from that gathered tile, and SSE-transpose 16x16 blocks of the
     tile into an SNP-major matrix TF (S x n bytes, transparent huge pages, so no TLB pressure).
   Phase B (parallel over SNP chunks, ordered writes): VCF rows straight from TF rows. */

typedef struct { const snps_t *z; const unsigned char *TF; size_t snps_per_chunk, stride; } tfvcf_ctx_t;

static size_t fmt_vcf_tf(const void *vctx, size_t c, char *buf, unsigned char *scratch)
{
  (void) scratch;
  const tfvcf_ctx_t *x = vctx; const snps_t *z = x->z;
  char *o = buf;
  const size_t n = z->nr;
  if (c == 0) {
    o += sprintf(o, "##fileformat=VCFv4.1\n##contig=<ID=1,length=%i>\n", (int) z->L);
    o += sprintf(o, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    o += sprintf(o, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (size_t r = 0; r < n; r++) { *o++ = '\t'; size_t nl = name_len_c(&z->R[r]); memcpy(o, z->R[r].name, nl); o += nl; }
    *o++ = '\n';
  }
  size_t a = c * x->snps_per_chunk, e = a + x->snps_per_chunk < z->S ? a + x->snps_per_chunk : z->S;
  for (size_t j = a; j < e; j++) o = vcf_row(o, x->TF + j * x->stride, z->pos[j], n, z->f[z->pos[j]]);
  return (size_t) (o - buf);
}

static void write_big_n(const char *vcf_fn, const char *phy_fn, const char *fa_fn, int want_phy, int want_fa,
                        int output_reference, const snps_t *z)
{
  const size_t n = z->nr, S = z->S;
  int fdv = open_out(vcf_fn);
  int fdp = want_phy ? open_out(phy_fn) : -1;
  int fdf = want_fa ? open_out(fa_fn) : -1;
  prealloc(fdv, vcf_bytes(z));
  if (want_phy || want_fa) { size_t rb = rows_bytes(z); if (fdp >= 0) prealloc(fdp, rb); if (fdf >= 0) prealloc(fdf, rb); }
  /* matrix rows padded to whole cache lines: a full 64-record block then writes exactly one aligned line per row,
     assembled in L1 and stored with non-temporal stores (no read-for-ownership of the matrix) */
  const size_t ns = (n + 63) & ~(size_t) 63;
  unsigned char *TF = thp_alloc(S * ns + 64);
  size_t rbu = (CHUNK_TARGET / (S + 64)) / 64 * 64;
  if (rbu < 64) rbu = 64;
  long nunits = (long) ((n + rbu - 1) / rbu);
  /* Phase A rows: a unit's gathered rows stay in Gt (stride S + 1, each followed by its '\n'), and the FASTA and
     PHYLIP chunks are written from there with writev, behind per-record name prefixes. */
  size_t maxname = 0;
  for (size_t r = 0; r < n; r++) { size_t nl = name_len_c(&z->R[r]); if (nl > maxname) maxname = nl; }
  const size_t GS = S + 1;
  volatile size_t fa_turn = 0, phy_turn = 0;
#pragma omp parallel
  {
    unsigned char *Gt = malloc(rbu * GS + 64);
    unsigned char *tile = aligned_alloc(64, 16 * 64);
    char *hf = malloc(S + 128), *hp = malloc(S + 128);              /* unit 0: header and pseudo-reference rows */
    char *fpre = malloc(rbu * (maxname + 2) + 64), *ppre = malloc(rbu * (maxname + 2) + 64);
    struct iovec *fv = malloc((2 * rbu + 4) * sizeof(struct iovec)), *pv = malloc((2 * rbu + 4) * sizeof(struct iovec));
#pragma omp for schedule(dynamic, 1)
    for (long u = 0; u < nunits; u++) {
      size_t r0 = (size_t) u * rbu, r1 = r0 + rbu < n ? r0 + rbu : n;
      size_t nf = 0, np_ = 0;
      if (u == 0) {
        char *fo = hf, *po = hp;
        if (want_phy) po += sprintf(po, "%d %d\n", (int) n + (output_reference ? 1 : 0), (int) S);
        if (output_reference) {
          if (want_fa) { memcpy(fo, ">pseudo_reference_sequence\n", 27); fo += 27; for (size_t j = 0; j < S; j++) fo[j] = (char) z->f[z->pos[j]]; fo += S; *fo++ = '\n'; }
          if (want_phy) { memcpy(po, "pseudo_reference_sequence\t", 26); po += 26; for (size_t j = 0; j < S; j++) po[j] = (char) z->f[z->pos[j]]; po += S; *po++ = '\n'; }
        }
        if (fo > hf) { fv[nf].iov_base = hf; fv[nf].iov_len = (size_t) (fo - hf); nf++; }
        if (po > hp) { pv[np_].iov_base = hp; pv[np_].iov_len = (size_t) (po - hp); np_++; }
      }
      for (size_t rb = r0; rb < r1; rb += 64) {
        size_t k = r1 - rb < 64 ? r1 - rb : 64;
        unsigned char *Gb = Gt + (rb - r0) * GS;
        for (size_t i = 0; i < k; i++) { gather_row(z, z->R[rb + i].seq, Gb + i * GS); Gb[i * GS + S] = '\n'; }
        size_t kf = k / 16 * 16; /* records in full 16-blocks */
        size_t j0 = 0;
        if (k == 64) {
          for (; j0 + 16 <= S; j0 += 16) {
            for (size_t i0 = 0; i0 < 64; i0 += 16) transpose16(Gb + i0 * GS + j0, GS, tile + i0, 64);
            for (size_t m = 0; m < 16; m++) {
              __m256i *d = (__m256i *) (TF + (j0 + m) * ns + rb);
              _mm256_stream_si256(d, _mm256_load_si256((const __m256i *) (tile + m * 64)));
              _mm256_stream_si256(d + 1, _mm256_load_si256((const __m256i *) (tile + m * 64 + 32)));
            }
          }
        } else {
          for (; j0 + 16 <= S; j0 += 16)
            for (size_t i0 = 0; i0 < kf; i0 += 16) transpose16(Gb + i0 * GS + j0, GS, TF + j0 * ns + rb + i0, ns);
        }
        for (size_t j = j0; j < S; j++) for (size_t i = 0; i < kf; i++) TF[j * ns + rb + i] = Gb[i * GS + j];
        for (size_t j = 0; j < S; j++) for (size_t i = kf; i < k; i++) TF[j * ns + rb + i] = Gb[i * GS + j];
      }
      char *fo = fpre, *po = ppre;
      for (size_t r = r0; r < r1; r++) {
        const rec_t *R = &z->R[r];
        size_t nl = name_len_c(R);
        unsigned char *row = Gt + (r - r0) * GS;
        if (want_fa) {
          fv[nf].iov_base = fo; *fo++ = '>'; memcpy(fo, R->name, nl); fo += nl; *fo++ = '\n'; fv[nf].iov_len = nl + 2; nf++;
          fv[nf].iov_base = row; fv[nf].iov_len = GS; nf++;
        }
        if (want_phy) {
          pv[np_].iov_base = po; memcpy(po, R->name, nl); po += nl; *po++ = '\t'; pv[np_].iov_len = nl + 1; np_++;
          pv[np_].iov_base = row; pv[np_].iov_len = GS; np_++;
        }
      }
      /* per-file turns: one thread can write FASTA unit u+1 while another writes PHYLIP unit u */
      wait_turn(&fa_turn, (size_t) u);
      if (want_fa) writev_all(fdf, fv, nf);
      pass_turn(&fa_turn, (size_t) u);
      wait_turn(&phy_turn, (size_t) u);
      if (want_phy) writev_all(fdp, pv, np_);
      pass_turn(&phy_turn, (size_t) u);
    }
    _mm_sfence(); /* order the non-temporal matrix stores before Phase B reads them */
    free(Gt); free(tile); free(hf); free(hp); free(fpre); free(ppre); free(fv); free(pv);
  }
  prof_mark("rows+TF");
  size_t row_typ = 2 * n + 40;
  size_t spc = CHUNK_TARGET / row_typ;
  if (spc < 1) spc = 1;
  size_t nchunks = (S + spc - 1) / spc;
  size_t header = 256;
  for (size_t r = 0; r < n; r++) header += name_len_c(&z->R[r]) + 1;
  tfvcf_ctx_t x = { z, TF, spc, ns };
  run_chunked(fdv, nchunks, header + spc * (3 * n + 160), 0, fmt_vcf_tf, &x);
  close_out(fdv); if (fdp >= 0) close_out(fdp); if (fdf >= 0) close_out(fdf);
}

/* ------------------------------------------------------------------ driver */

static void strip_directory_from_filename(const char *input_filename, char *output_filename)
{
  int i, end_index = 0, last_forward_slash_index = -1;
  for (i = 0; i < FILENAME_MAX; i++) {
    if (input_filename[i] == '/') last_forward_slash_index = i;
    if (input_filename[i] == '\0' || input_filename[i] == '\n') { end_index = i; break; }
  }
  int current_index = 0;
  for (i = last_forward_slash_index + 1; i < end_index; i++) output_filename[current_index++] = input_filename[i];
  output_filename[current_index] = '\0';
}

static int same_file(const char *a, const struct stat *in_st)
{
  struct stat st;
  if (stat(a, &st) != 0) return 0;
  return st.st_dev == in_st->st_dev && st.st_ino == in_st->st_ino;
}

typedef struct {
  int fasta, vcf, phylip, reference, constant, pure, mono;
  const char *input; const char *output;
} opts_t;

static int run(const opts_t *op)
{
  char vcf_fn[FILENAME_MAX + 16], phy_fn[FILENAME_MAX + 16], fa_fn[FILENAME_MAX + 32];
  int want_fasta = op->fasta || (!op->vcf && !op->phylip && !op->fasta);
  int cso = -1; /* -C output, opened before reading as in snp-sites */
  if (op->constant) cso = open_out(op->output);
  else {
    char base[FILENAME_MAX], nodir[FILENAME_MAX];
    strip_directory_from_filename(op->input, nodir);
    strncpy(base, nodir, FILENAME_MAX - 1); base[FILENAME_MAX - 1] = 0;
    if (op->output[0] != '\0') { strncpy(base, op->output, FILENAME_MAX - 1); base[FILENAME_MAX - 1] = 0; }
    int suffix = (op->vcf + op->phylip + op->fasta) > 1 || op->output[0] == '\0';
    snprintf(vcf_fn, sizeof vcf_fn, "%s%s", base, suffix ? ".vcf" : "");
    snprintf(phy_fn, sizeof phy_fn, "%s%s", base, suffix ? ".phylip" : "");
    snprintf(fa_fn, sizeof fa_fn, "%s%s", base, suffix ? ".snp_sites.aln" : "");
  }
  /* If an output path is the input file itself, read the input into private memory first. */
  int force_private = 0;
  struct stat in_st;
  if (stat(op->input, &in_st) == 0) {
    if (op->constant) force_private = same_file(op->output, &in_st);
    else force_private = (op->vcf && same_file(vcf_fn, &in_st)) || (op->phylip && same_file(phy_fn, &in_st)) ||
                         (want_fasta && same_file(fa_fn, &in_st));
  }
  input_t in;
  load_input(op->input, &in, force_private);
  prof_mark("read");

  rec_t *R = NULL;
  size_t nr = 0, L = 0, W = 0;
  int fast = parse_fast(in.buf, in.n, &R, &nr, &L, &W);
  if (!fast) {
    kq_t k; memset(&k, 0, sizeof k);
    k.b = in.buf; k.n = in.n;
    size_t cap = 1024;
    nr = 0;
    R = malloc(cap * sizeof(rec_t));
    for (;;) {
      if (nr == cap) { cap *= 2; R = realloc(R, cap * sizeof(rec_t)); }
      long l = kq_read(&k, &R[nr]);
      if (l < 0) break;
      nr++;
    }
    L = nr ? R[0].len : 0;
    W = L; /* contiguous */
  }
  prof_mark("parse");
  int unequal = 0;
  for (size_t r = 1; r < nr; r++) if (R[r].len != L) { unequal = r; break; }
  if (unequal && !op->constant) {
    fprintf(stderr, "Alignment %s contains sequences of unequal length. Expected length is %i but got %i in sequence %.*s\n\n",
            op->input, (int) L, (int) R[unequal].len, (int) name_len_c(&R[unequal]), (const char *) R[unequal].name);
    exit(EXIT_FAILURE);
  }
  int pure = op->constant ? 0 : op->pure, mono = op->constant ? 0 : op->mono;
  g_np_only = pure && mono && !op->reference && !op->vcf;
  unsigned char *f = big_alloc(L + 1), *var = big_alloc(L + 1), *np = pure ? big_alloc(L + 1) : NULL;
  if (unequal) scan_columns_shadow(R, nr, L, f, var);
  else if (fast) {
    int lay_ok = L <= 262144 ? scan_records_layout(R, nr, L, W, pure, f, var, np)
                             : scan_columns_layout(R, nr, L, W, pure, f, var, np);
    if (!lay_ok) {
      /* predicted layout was wrong: start again with the exact parser */
      free(R);
      kq_t k; memset(&k, 0, sizeof k);
      k.b = in.buf; k.n = in.n;
      size_t cap = 1024;
      nr = 0;
      R = malloc(cap * sizeof(rec_t));
      for (;;) {
        if (nr == cap) { cap *= 2; R = realloc(R, cap * sizeof(rec_t)); }
        long l = kq_read(&k, &R[nr]);
        if (l < 0) break;
        nr++;
      }
      L = nr ? R[0].len : 0;
      W = L;
      fast = 0;
      for (size_t r = 1; r < nr; r++) if (R[r].len != L) { unequal = (int) r; break; }
      if (unequal && !op->constant) {
        fprintf(stderr, "Alignment %s contains sequences of unequal length\n\n", op->input);
        exit(EXIT_FAILURE);
      }
      f = realloc(f, L + 1); var = realloc(var, L + 1); if (pure) np = realloc(np, L + 1);
      if (unequal) scan_columns_shadow(R, nr, L, f, var);
      else scan_columns(R, nr, L, pure, f, var, np);
    }
  } else scan_columns(R, nr, L, pure, f, var, np);

  prof_mark("scan");
  /* classify (parallel over fixed 64 Ki-column chunks, AVX2 masks): output columns, -C constant-site counts */
  size_t S = 0;
  long counts[4] = {0, 0, 0, 0};
  const size_t CK = 65536;
  long nck = (long) ((L + CK - 1) / CK);
  size_t *ccount = calloc((size_t) nck + 1, sizeof(size_t));
  long *ccnt4 = calloc((size_t) nck * 4 + 4, sizeof(long));
  int *pos = (int *) big_alloc((L + 32) * sizeof(int));
  for (int pass = 0; pass < 2; pass++) {
#pragma omp parallel for schedule(static)
    for (long ci = 0; ci < nck; ci++) {
      size_t j0 = (size_t) ci * CK, j1 = j0 + CK < L ? j0 + CK : L;
      size_t k = pass ? ccount[ci] : 0, cnt = 0;
      long c4[4] = {0, 0, 0, 0};
      for (size_t j = j0; j < j1; j += 32) {
        size_t m = j1 - j < 32 ? j1 - j : 32;
        uint32_t outm = 0, am = 0, cm = 0, gm = 0, tm = 0;
        if (m == 32) {
          __m256i fv = _mm256_loadu_si256((const __m256i *) (f + j));
          __m256i vv = _mm256_loadu_si256((const __m256i *) (var + j));
          __m256i zero = _mm256_setzero_si256();
          __m256i excl = pure ? _mm256_xor_si256(_mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i *) (np + j)), zero), _mm256_set1_epi8(-1))
                              : _mm256_cmpeq_epi8(fv, _mm256_set1_epi8('#'));
          __m256i isv = _mm256_or_si256(_mm256_xor_si256(_mm256_cmpeq_epi8(vv, zero), _mm256_set1_epi8(-1)),
                                        _mm256_cmpeq_epi8(fv, _mm256_set1_epi8('>')));
          __m256i outv = _mm256_andnot_si256(excl, mono ? _mm256_set1_epi8(-1) : isv);
          outm = (uint32_t) _mm256_movemask_epi8(outv);
          if (op->constant) {
            am = (uint32_t) _mm256_movemask_epi8(_mm256_cmpeq_epi8(fv, _mm256_set1_epi8('A')));
            cm = (uint32_t) _mm256_movemask_epi8(_mm256_cmpeq_epi8(fv, _mm256_set1_epi8('C')));
            gm = (uint32_t) _mm256_movemask_epi8(_mm256_cmpeq_epi8(fv, _mm256_set1_epi8('G')));
            tm = (uint32_t) _mm256_movemask_epi8(_mm256_cmpeq_epi8(fv, _mm256_set1_epi8('T')));
          }
        } else {
          for (size_t t = 0; t < m; t++) {
            unsigned char fj = f[j + t];
            int ex = pure ? np[j + t] != 0 : (fj == '#');
            int iv = !ex && (var[j + t] || fj == '>');
            if (iv || (mono && !ex)) outm |= 1u << t;
            am |= (uint32_t) (fj == 'A') << t; cm |= (uint32_t) (fj == 'C') << t;
            gm |= (uint32_t) (fj == 'G') << t; tm |= (uint32_t) (fj == 'T') << t;
          }
        }
        if (!pass) {
          cnt += (size_t) __builtin_popcount(outm);
          if (op->constant) {
            c4[0] += __builtin_popcount(am & ~outm); c4[1] += __builtin_popcount(cm & ~outm);
            c4[2] += __builtin_popcount(gm & ~outm); c4[3] += __builtin_popcount(tm & ~outm);
          }
        } else if (outm == 0xFFFFFFFFu) {
          __m256i base = _mm256_set1_epi32((int) j);
          _mm256_storeu_si256((__m256i *) (pos + k), _mm256_add_epi32(base, _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7)));
          _mm256_storeu_si256((__m256i *) (pos + k + 8), _mm256_add_epi32(base, _mm256_setr_epi32(8, 9, 10, 11, 12, 13, 14, 15)));
          _mm256_storeu_si256((__m256i *) (pos + k + 16), _mm256_add_epi32(base, _mm256_setr_epi32(16, 17, 18, 19, 20, 21, 22, 23)));
          _mm256_storeu_si256((__m256i *) (pos + k + 24), _mm256_add_epi32(base, _mm256_setr_epi32(24, 25, 26, 27, 28, 29, 30, 31)));
          k += 32;
        } else {
          while (outm) { pos[k++] = (int) (j + (size_t) __builtin_ctz(outm)); outm &= outm - 1; }
        }
      }
      if (!pass) { ccount[ci] = cnt; for (int q = 0; q < 4; q++) ccnt4[ci * 4 + q] = c4[q]; }
    }
    if (!pass) { /* exclusive prefix sums over chunks */
      for (long ci = 0; ci < nck; ci++) {
        size_t c = ccount[ci]; ccount[ci] = S; S += c;
        for (int q = 0; q < 4; q++) counts[q] += ccnt4[ci * 4 + q];
      }
      ccount[nck] = S;
      if (S == 0 || op->constant) break;
    }
  }
  if (S == 0) die("Warning: No SNPs were detected so there is nothing to output.\n");
  if (op->constant) {
    char t[128];
    int m = snprintf(t, sizeof t, "%d,%d,%d,%d\n", (int) counts[0], (int) counts[1], (int) counts[2], (int) counts[3]);
    write_all(cso, t, (size_t) m);
    close(cso);
    return 0;
  }
  /* byte offsets of the output columns (pos + line index), and contiguous copy segments (consecutive columns on the
     same line), computed per chunk in parallel and merged at chunk boundaries */
  uint32_t *off = (uint32_t *) big_alloc((S + 32) * sizeof(uint32_t));
  uint32_t *seg_off = (uint32_t *) big_alloc((S + 1) * sizeof(uint32_t)), *seg_len = (uint32_t *) big_alloc((S + 1) * sizeof(uint32_t));
  size_t *cseg = calloc((size_t) nck + 1, sizeof(size_t));
#pragma omp parallel for schedule(static)
  for (long ci = 0; ci < nck; ci++) {
    size_t a = ccount[ci], e = ccount[ci + 1];
    if (a == e) continue;
    size_t line = W >= L ? 0 : (size_t) pos[a] / W, line_end = W >= L ? (size_t) -1 : (line + 1) * W;
    size_t ns = 0; /* segments are written at index a.. (worst case one per column) */
    /* a segment continues while the byte offset grows by exactly 1 (pos is strictly increasing, so this also
       means consecutive columns on the same line) */
    uint32_t prev = 0, run = 0;
    for (size_t j = a; j < e; j++) {
      size_t pj = (size_t) pos[j];
      while (pj >= line_end) { line++; line_end += W; }
      uint32_t oj = (uint32_t) (pj + line);
      off[j] = oj;
      if (j == a || oj != prev + 1) { if (ns) seg_len[a + ns - 1] = run; seg_off[a + ns] = oj; ns++; run = 0; }
      run++; prev = oj;
    }
    if (ns) seg_len[a + ns - 1] = run;
    cseg[ci] = ns;
  }
  size_t nseg = 0;
  for (long ci = 0; ci < nck; ci++) {
    size_t a = ccount[ci];
    for (size_t q = 0; q < cseg[ci]; q++) {
      uint32_t so = seg_off[a + q], sl = seg_len[a + q];
      size_t jfirst = a; /* first column of this chunk: merge only the chunk's first segment */
      if (q == 0 && nseg && so == seg_off[nseg - 1] + seg_len[nseg - 1] && pos[jfirst] == pos[jfirst - 1] + 1) seg_len[nseg - 1] += sl;
      else { seg_off[nseg] = so; seg_len[nseg] = sl; nseg++; }
    }
  }
  /* 16-byte windows (see snps_t): greedy from the first unassigned column; a window must stay inside the record's
     raw bytes (raw length = L plus its line breaks), so a window near the end is shifted back. */
  size_t rawlen = W >= L ? L : L + (L - 1) / W;
  uint32_t *win_src = malloc((S + 1) * sizeof(uint32_t));
  uint8_t *win_cnt = malloc(S + 1);
  __m128i *win_shuf = aligned_alloc(16, (S + 1) * sizeof(__m128i));
  size_t nwin = 0;
  int use_segs = nseg * 16 <= S, use_wins = 0;
  if (!use_segs && rawlen >= 16) {
    for (size_t j = 0; j < S;) {
      size_t st = off[j];
      if (st + 16 > rawlen) st = rawlen - 16;
      uint8_t sh[16]; int cnt = 0;
      while (j < S && off[j] < st + 16) { sh[cnt++] = (uint8_t) (off[j] - st); j++; }
      for (int q = cnt; q < 16; q++) sh[q] = 0x80;
      win_src[nwin] = (uint32_t) st; win_cnt[nwin] = (uint8_t) cnt;
      win_shuf[nwin] = _mm_loadu_si128((const __m128i *) sh);
      nwin++;
    }
    use_wins = nwin * 2 <= S; /* worth it if windows hold >= 2 columns on average */
  }
  snps_t z = { R, nr, pos, S, off, seg_off, seg_len, nseg, use_segs, win_src, win_cnt, win_shuf, nwin, use_wins, f, L };
  if (g_prof) fprintf(stderr, "PROFILE gather: S=%zu nseg=%zu nwin=%zu mode=%s\n", S, nseg, nwin, use_wins ? "windows" : use_segs ? "segments" : "index");
  prof_mark("classify");
  /* large n with scattered output columns: rows + transposed matrix path (see write_big_n) */
  size_t contig_tiles = 0, ntiles = S / 64;
  for (size_t t = 0; t < ntiles; t++) contig_tiles += off[t * 64 + 63] - off[t * 64] == 63;
  int big_n = op->vcf && nr >= 512 && contig_tiles * 2 < ntiles;
  if (big_n) { write_big_n(vcf_fn, phy_fn, fa_fn, op->phylip, want_fasta, op->reference, &z); prof_mark("vcf"); return 0; }
  /* files are created in snp-sites' order (VCF, PHYLIP, FASTA), then filled (concurrently when several) */
  job_t jobs[3]; int fds[4], nj = 0, close_extra = -1;
  if (op->vcf) { fds[nj] = open_out(vcf_fn); jobs[nj] = prepare_vcf(fds[nj], &z); nj++; }
  if (op->phylip && want_fasta) { /* one job, one gather per record, two files */
    int fp = open_out(phy_fn), ff = open_out(fa_fn);
    fds[nj] = ff; jobs[nj] = prepare_rows(ff, &z, op->reference, 0);
    jobs[nj].fd2 = fp; jobs[nj].fmt2 = fmt_rows2;
    if (z.S >= 1024) { /* long rows: PHYLIP written from the FASTA chunk's row bytes (writev) */
      rows_ctx_t *rx = (rows_ctx_t *) jobs[nj].ctx;
      rx->iov = 1; jobs[nj].iov2 = 1; jobs[nj].scratch_bytes = (2 * rx->maxrows + 8) * sizeof(struct iovec);
    }
    nj++;
    close_extra = fp;
  } else {
    if (op->phylip) { fds[nj] = open_out(phy_fn); jobs[nj] = prepare_rows(fds[nj], &z, op->reference, 1); nj++; }
    if (want_fasta) { fds[nj] = open_out(fa_fn); jobs[nj] = prepare_rows(fds[nj], &z, op->reference, 0); nj++; }
  }
  if (op->vcf) prealloc(fds[0], vcf_bytes(&z));
  if (op->phylip || want_fasta) {
    size_t rb = rows_bytes(&z);
    for (int i = op->vcf ? 1 : 0; i < nj; i++) prealloc(fds[i], rb);
    if (close_extra >= 0) prealloc(close_extra, rb);
  }
  run_jobs(jobs, nj);
  for (int i = 0; i < nj; i++) close_out(fds[i]);
  if (close_extra >= 0) close_out(close_extra);
  prof_mark("outputs");
  return 0;
}

static void print_usage(void)
{
  printf("Usage: snp-sites [-rmvpcbhV] [-o output_filename] <file>\n");
  printf("This program finds snp sites from a multi FASTA alignment file.\n");
  printf(" -r     output internal pseudo reference sequence\n");
  printf(" -m     output a multi fasta alignment file (default)\n");
  printf(" -v     output a VCF file\n");
  printf(" -p     output a phylip file\n");
  printf(" -o STR specify an output filename [STDOUT]\n");
  printf(" -c     only output columns containing exclusively ACGT\n");
  printf(" -C     only output count of constant sites (suitable for IQ-TREE -fconst) and nothing else\n");
  printf(" -b     output monomorphic sites, used for BEAST\n");
  printf(" -h     this help message\n");
  printf(" -V     print version and exit\n");
  printf(" <file> input alignment file which can optionally be gzipped\n\n");
  printf("Example: creating files for BEAST\n");
  printf(" snp-sites -cb -o outputfile.aln inputfile.aln\n\n");
  printf("If you use this program, please cite:\n");
  printf("\"SNP-sites: rapid efficient extraction of SNPs from multi-FASTA alignments\",\n");
  printf("Andrew J. Page, Ben Taylor, Aidan J. Delaney, Jorge Soares, Torsten Seemann, Jacqueline A. Keane, Simon R. Harris,\n");
  printf("Microbial Genomics 2(4), (2016). http://dx.doi.org/10.1099/mgen.0.000056\n");
}

int main(int argc, char **argv)
{
  static char output_filename[FILENAME_MAX] = {"/dev/stdout"};
  static char multi_fasta_filename[FILENAME_MAX] = {""};
  opts_t op; memset(&op, 0, sizeof op);
  int c;
  if (getenv("SWIFTSITES_PROFILE")) { g_prof = 1; g_t_last = now_s(); }
  while ((c = getopt(argc, argv, "mvrbpcCo:V")) != -1)
    switch (c) {
      case 'm': op.fasta = 1; break;
      case 'v': op.vcf = 1; break;
      case 'V': printf("%s %s\n", "snp-sites", "2.5.1"); return 0;
      case 'p': op.phylip = 1; break;
      case 'r': op.reference = 1; break;
      case 'c': op.pure = 1; break;
      case 'C': op.constant = 1; break;
      case 'b': op.mono = 1; break;
      case 'o': strncpy(output_filename, optarg, FILENAME_MAX); output_filename[FILENAME_MAX - 1] = 0; break;
      case 'h': print_usage(); exit(EXIT_SUCCESS);
      default: op.fasta = 1;
    }
  if (optind < argc) {
    if (access(argv[optind], F_OK) == -1) {
      fprintf(stderr, "ERROR: cannot access input alignment file '%s'\n", argv[optind]);
      fflush(stderr);
      exit(EXIT_FAILURE);
    }
    strncpy(multi_fasta_filename, argv[optind], FILENAME_MAX); multi_fasta_filename[FILENAME_MAX - 1] = 0;
    op.input = multi_fasta_filename; op.output = output_filename;
    run(&op);
  } else {
    print_usage();
  }
  exit(EXIT_SUCCESS);
}
