#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
LD=libdeflate/lib
# libdeflate's decompressor is measurably faster at -O2 than at -O3 with gcc 13
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
for f in deflate_decompress gzip_decompress utils crc32 x86/cpu_features; do
  gcc -O2 -march=native -Wall -Ilibdeflate -c "$LD/$f.c" -o "$T/$(basename "$f").o"
done
gcc -O3 -march=native -fopenmp -Wall -Ilibdeflate -o swiftsites swiftsites.c gzload.c "$T"/*.o -lz
