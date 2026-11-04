/* SwiftSites: gzip input loader (GPL-3.0). */
#ifndef SWIFT_GZLOAD_H
#define SWIFT_GZLOAD_H
#include <stddef.h>
/* Decodes a gzip file image c[0..csize) into fresh memory (*out, *outlen). Returns 0 if libdeflate rejects the data
   (the caller then falls back to zlib). Uses two-thread split decoding for one large member when possible. */
int gz_load(const unsigned char *c, size_t csize, unsigned char **out, size_t *outlen);
extern int gz_split_used; /* diagnostics: set when the last gz_load used the split path */
#endif
