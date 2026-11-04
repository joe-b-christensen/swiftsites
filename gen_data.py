#!/usr/bin/env python3
"""SwiftSites: dataset specs and fuzz-case generator.

Benchmark sets are mostly REAL public alignments (data/, see data/README.txt), plus two simulated
sets that exercise what the real ones lack (lower-case soft-masking, '?' and every IUPAC code, wrapped lines, gzip):
  * kind "real":    an xz-compressed alignment in data/, decompressed as is;
  * kind "real_gz": the same, re-compressed with gzip (level 6, mtime 0) so the gzip input path is benchmarked;
  * kind "file":    a file shipped in the snp-sites repo (ext/snp-sites/...);
  * kind "sim":     simulated with a fixed seed (core-genome-like: a reference, clonal lineages, N runs, gaps).
Fuzz cases are tiny adversarial FASTA/FASTQ files (optionally gzipped) plus random snp-sites flags.
"""
import gzip
import random

import numpy as np

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
IUPAC = b"RYSWKMBDHV"

# name -> spec. `args` are the snp-sites flags of the benchmark (besides `-o OUT` and the input file).
DEV_SPECS = {
    "real_vch_beast": dict(kind="real", src="data/vcholerae_n100.fa.xz", args=["-c", "-b"]),
    "real_mtb_all4": dict(kind="real", src="data/mtb_n1788_snps.fa.xz", args=["-m", "-v", "-p", "-r"]),
    "real_sars2_8596_all4": dict(kind="real", src="data/sarscov2_n8596.fa.xz", args=["-m", "-v", "-p", "-r"]),
    "real_vch_gz_vcf_phylip": dict(kind="real_gz", src="data/vcholerae_n100.fa.xz", args=["-v", "-p"]),
    "sim_softmasked_default": dict(kind="sim", n=1200, L=250_000, seed=31, wrap=70, lower=0.25, args=[]),
}
HOLDOUT_SPECS = {
    "h_real_vch_mono_vcf": dict(kind="real", src="data/vcholerae_n100.fa.xz", args=["-b", "-v"]),
    "h_real_sars2_global_vcf_ref": dict(kind="real", src="data/sarscov2_global_n4637.fa.xz", args=["-v", "-r"]),
    "h_real_mtb_phylip": dict(kind="real", src="data/mtb_n1788_snps.fa.xz", args=["-p"]),
    "h_real_sars2_europe_mono": dict(kind="real", src="data/sarscov2_europe_n2012.fa.xz",
                                     args=["-b", "-p", "-r"]),
    "h_real_salmonella_gz_vcf_phylip": dict(kind="file", src="ext/snp-sites/example_data/salmonella_serovars_core_genes.aln.gz",
                                            args=["-v", "-p"]),
    "h_sim_gz_default_vcf": dict(kind="sim", n=600, L=400_000, seed=47, wrap=60, lower=0.1, gz=True,
                                        args=["-m", "-v"]),
}


def simulate(spec):
    rng = np.random.default_rng(spec["seed"] * 11 + 3)
    n, L = spec["n"], spec["L"]
    ref = BASES[rng.choice(4, size=L, p=[0.26, 0.24, 0.24, 0.26])]
    lineages = []
    for _ in range(12):
        s = ref.copy()
        k = int(rng.integers(200, 1500))
        pos = rng.integers(0, L, k)
        s[pos] = BASES[rng.integers(0, 4, k)]
        lineages.append(s)
    seqs = []
    for i in range(n):
        s = lineages[int(rng.integers(len(lineages)))].copy()
        k = int(rng.poisson(40))
        pos = rng.integers(0, L, k)
        s[pos] = BASES[rng.integers(0, 4, k)]
        for _ in range(int(rng.integers(0, 6))):  # N runs (missing coverage) and gap runs (deletions)
            a = int(rng.integers(0, L - 3000))
            s[a:a + int(rng.integers(50, 3000))] = ord("N") if rng.random() < 0.7 else ord("-")
        k = int(rng.poisson(3))
        if k:
            s[rng.integers(0, L, k)] = np.frombuffer(IUPAC + b"?", dtype=np.uint8)[rng.integers(0, 11, k)]
        if rng.random() < spec.get("lower", 0):  # soft-masked stretch
            a = int(rng.integers(0, L - 20000))
            seg = s[a:a + 20000]
            s[a:a + 20000] = np.where((seg >= 65) & (seg <= 90), seg + 32, seg)
        seqs.append(s)
    names = [f"isolate_{spec['seed']}_{i:05d}" for i in range(n)]
    return names, seqs


def write_dataset(spec, path):
    names, seqs = simulate(spec)
    wrap = spec.get("wrap", 0)
    opener = (lambda p: gzip.GzipFile(p, "wb", compresslevel=6, mtime=0)) if spec.get("gz") else \
        (lambda p: open(p, "wb"))
    with opener(path) as f:
        for nm, s in zip(names, seqs):
            f.write(b">" + nm.encode() + b"\n")
            b = s.tobytes()
            if wrap:
                for i in range(0, len(b), wrap):
                    f.write(b[i:i + wrap] + b"\n")
            else:
                f.write(b + b"\n")
    return (len(seqs), len(seqs[0]))


# ---------------------------------------------------------------- fuzz

CORE = b"ACGT"
EXTRA = b"acgtNn-?RYSWKMBDHVryswkmbdhv"
WEIRD = b"#*.XxUu0~"  # rare characters with special or odd treatment ('#' excludes a column in -c/-b mode logic)
FLAG_POOL = ["-m", "-v", "-p", "-r", "-c", "-b"]


def _records(rng, n, L, extra_p, names_mode, equal):
    weird_p = 0.01 if rng.random() < 0.1 else 0.0
    base = bytes(rng.choice(CORE) for _ in range(L))
    recs, names = [], []
    for i in range(n):
        if names_mode == "dup" and i and rng.random() < 0.4:
            nm = names[rng.randrange(len(names))]
        elif names_mode == "odd" and rng.random() < 0.4:
            nm = rng.choice(["", "a b c", "x\ty", "q,r", "long_" * rng.randint(10, 300), "ü-name", "dup"])
        else:
            nm = f"s{i}"
        names.append(nm)
        s = bytearray(base)
        for _ in range(rng.choice([0, 0, 1, 2, 3, max(1, L // 10), L])):
            if not s:
                break
            j = rng.randrange(len(s))
            r = rng.random()
            s[j] = rng.choice(WEIRD) if r < weird_p else rng.choice(EXTRA) if r < extra_p else rng.choice(CORE)
        if rng.random() < 0.1 and s:  # N or gap runs
            a = rng.randrange(len(s))
            b = min(len(s), a + rng.randint(1, max(1, L // 3)))
            s[a:b] = (b"N" if rng.random() < 0.5 else b"-") * (b - a)
        if rng.random() < 0.1:  # soft-masking
            s = bytearray(bytes(s).lower())
        if not equal and rng.random() < 0.3:
            s = s[:max(0, len(s) + rng.randint(-4, 4))] if rng.random() < 0.5 else s + b"A" * rng.randint(1, 4)
        recs.append((nm, bytes(s)))
    if recs and rng.random() < 0.08:  # a literal '#' early in a column (it collides with snp-sites' internal marker)
        j = 0 if rng.random() < 0.7 else rng.randrange(len(recs))
        s = bytearray(recs[j][1])
        for _ in range(rng.randint(1, 3)):
            if s:
                s[rng.randrange(len(s))] = ord("#")
        recs[j] = (recs[j][0], bytes(s))
    return recs


def _render(recs, rng):
    fastq = rng.random() < 0.03
    wrap = 0 if fastq else rng.choice([0, 0, 0, 1, 7, 60])
    nl = b"\r\n" if rng.random() < 0.04 else b"\n"
    buf = bytearray()
    for nm, s in recs:
        if fastq:
            buf += b"@" + nm.encode() + nl + s + nl + b"+" + nl + b"I" * len(s) + nl
            continue
        buf += b">" + nm.encode() + (b" a comment" if rng.random() < 0.05 else b"") + nl
        if rng.random() < 0.04:
            buf += nl
        if wrap and s:
            for i in range(0, len(s), wrap):
                buf += s[i:i + wrap] + nl
        else:
            buf += s + nl
    if rng.random() < 0.04:
        buf = buf.rstrip(b"\r\n")
    return bytes(buf)


def fuzz_case(rng: random.Random):
    """Return (data, suffix, flags). Flags contain the literal token {OUT} for the output name, or no -o at all."""
    r = rng.random()
    if r < 0.02:
        data = b""
    elif r < 0.04:
        data = rng.choice([b">", b">a\n", b">a\nACGT\n", b"ACGT\n", b">a\nAC\n>b\nAG\n", b">a\nA\n>b\nC\n>c\nN\n"])
    else:
        n = rng.choice([1, 2, 2, 3, 4, 5, 7, 10, 16, 33])
        L = rng.choice([1, 2, 3, 8, 31, 32, 33, 63, 64, 65, 100, 127, 128, 129, 200, 256, 300])
        extra_p = rng.choice([0.0, 0.0, 0.05, 0.2, 0.6])
        recs = _records(rng, n, L, extra_p, rng.choice(["plain", "plain", "dup", "odd"]), rng.random() > 0.03)
        data = _render(recs, rng)
    suffix = ".fa"
    if rng.random() < 0.1:
        data = gzip.compress(data, compresslevel=rng.choice([1, 6, 9]), mtime=0)
        suffix = ".fa.gz"
    if rng.random() < 0.06:
        flags = ["-C"]
    else:
        k = rng.choice([0, 0, 1, 1, 1, 2, 2, 3, 4])
        chosen = rng.sample(FLAG_POOL, k)
        if chosen and rng.random() < 0.2:  # combined short options, e.g. -cb
            flags = ["-" + "".join(f[1] for f in chosen)]
        else:
            flags = list(chosen)
    o = rng.random()
    if o < 0.85:
        flags += ["-o", "{OUT}"]
    elif o < 0.9:
        flags += ["-o", ""]  # empty: outputs are named after the input file, in the working directory
    # else: no -o, output goes to /dev/stdout
    return data, suffix, flags
