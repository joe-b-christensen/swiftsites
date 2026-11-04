# swiftsites

swiftsites is a drop-in replacement for [snp-sites](https://github.com/sanger-pathogens/snp-sites) 2.5.1, which reduces
a whole-genome multiple sequence alignment to its variable columns and writes them as multi-FASTA, VCF or relaxed
PHYLIP. It takes the same command line and writes the same files, byte for byte, much faster.

On four real benchmark sets (*Vibrio cholerae*, *Mycobacterium tuberculosis*, 8,596 SARS-CoV-2 genomes) it was 43×
faster than snp-sites with two threads and 26× with one; on five held-out real sets it was 39× faster with two threads.
On seven real inputs it was 5.0× faster than SNPick 1.0.2, the fastest alternative, whose output differs from snp-sites'
on three of them (writing variable-site FASTA, two threads). These are geometric means on a cloud VM with 2 vCPUs;
snp-sites has no thread option. One *V. cholerae* set took snp-sites 11.9 s and swiftsites 0.16 s. On a 12-thread AVX2
laptop swiftsites was 16× faster on the same four sets. Details are in the preprint (see Citation).

## Build

Requirements: Linux on x86-64 with AVX2, gcc with OpenMP, zlib headers (`zlib1g-dev` on Debian/Ubuntu). Tested with
gcc 13.3.0 on Ubuntu 24.04.

```
bash solver/build.sh        # produces solver/swiftsites
```

The build uses `-march=native`, so the binary is tuned for the build machine's CPU. CPUs without AVX2 and other
architectures are not supported.

## Usage

swiftsites takes the same options as snp-sites 2.5.1 and writes the same files:

```
solver/swiftsites -c -b -o core.aln core_full.aln        # pure-ACGT columns, including monomorphic ones (BEAST)
solver/swiftsites -m -v -p -r -o snps core_full.aln.gz   # snps.snp_sites.aln, snps.vcf, snps.phylip
solver/swiftsites -C -o fconst.txt core_full.aln         # constant-site counts for IQ-TREE -fconst
```

| Option | Meaning (as in snp-sites) |
|---|---|
| `-m` | multi-FASTA alignment of the variable columns (default) |
| `-v` | VCF |
| `-p` | relaxed PHYLIP |
| `-r` | add the internal pseudo-reference sequence |
| `-c` | only columns that contain exclusively A, C, G, T |
| `-b` | also output monomorphic columns (for BEAST) |
| `-C` | only the constant-site counts |
| `-o` | output name; with several formats the suffixes `.vcf`, `.phylip`, `.snp_sites.aln` are appended; `-o ""` names outputs after the input file |
| `-V` | prints `snp-sites 2.5.1`, so pipelines that check the version keep working |

- Every output file and standard output are byte-identical to snp-sites 2.5.1 with the same options, at every thread
  count. Where snp-sites fails on an input, swiftsites fails as well. To use it as a drop-in, put it on the `PATH` under
  the name `snp-sites`.
- `OMP_NUM_THREADS` sets the number of threads (default: all CPUs).
- Input: FASTA or FASTQ, plain or gzip-compressed. swiftsites holds the whole input in memory: its peak was 64 to
  527 MiB on the benchmark sets, against 16 to 488 MiB for snp-sites.
- swiftsites copies snp-sites' unusual behaviours, such as its handling of a literal `#` or `>` as the first known base
  of a column. Sequence names of 2,048 bytes or more are out of scope. Alignments longer than 2^31 columns and corrupt
  gzip input beyond the random-test cases were not tested.

## Check it yourself

```
python3 prepare.py              # builds snp-sites 2.5.1 (commit 657e74d), checks it against upstream's expected
                                # outputs, unpacks data/, generates the simulated sets and records snp-sites' outputs
                                # (needs gcc, git, zlib, numpy, network; about 10 min, 2 GB of disk)
python3 evaluate.py             # builds solver/; 300 random inputs and five benchmark sets compared byte for byte,
                                # plus timing; prints "gate: PASS" when every byte matches
python3 evaluate.py --holdout   # the same with the six holdout sets
```

Timings depend on the machine; the comparison does not. `data/` holds the five real alignments used in the preprint
(see `data/README.txt`); `prepare.py` fetches the *Salmonella* alignment from the snp-sites repository.

## Licence

GPL-3.0-or-later (`LICENSE`), as a derivative of snp-sites. Vendored code: libdeflate 1.26 (MIT, decompression sources
only, with small marked changes); see `LICENSE` and `solver/LICENSE_NOTE.txt`. The alignments in `data/` keep their own
licences (see `data/README.txt`).

## Citation

Christensen J. SwiftSites extracts variable sites 43-fold faster than snp-sites and 5-fold faster than SNPick, with
byte-identical output. Preprint, 2026. Please also cite snp-sites (Page AJ et al.,
Microbial Genomics 2:e000056, 2016). See `CITATION.cff`.
