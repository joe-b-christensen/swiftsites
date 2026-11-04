The five real alignments used in the SwiftSites preprint (about 5.5 MB in all): FASTA multiple sequence alignments,
xz-compressed. Sequences are unmodified: no re-alignment, masking, trimming or renaming. The only processing is
decompression and recompression, plus for sarscov2_n8596 the concatenation of the seven Nextstrain regional builds,
keeping the first occurrence of each sequence name. prepare.py decompresses them into bench/data/. Every set is used
whole.

  file                          sha256
  vcholerae_n100.fa.xz          37e63fba4f3b58f2d9f46ab6c17007415436f2ae63a2ea597b7c4a4d544c3880
  mtb_n1788_snps.fa.xz          8c05f6e0c4fa67bc7030504673a02f8eae4988f46958e2da73aede8023664d72
  sarscov2_europe_n2012.fa.xz   087efa164ccbc2df6f1372cbc596ab26fb9b9cca98852a27a3ac0418f10f09e8
  sarscov2_global_n4637.fa.xz   b3cca372df11e21d3ddac3a6b507d98c9057c76d3f8e31094ddf66b7b753b91a
  sarscov2_n8596.fa.xz          b205e5c601ed5cc17431e4e4ad27505e62075719332271f9a9d97b253fcba7b1

How they are used (gen_data.py DEV_SPECS and HOLDOUT_SPECS):
  development  real_vch_beast                   vcholerae_n100         -c -b
               real_mtb_all4                    mtb_n1788_snps         -m -v -p -r
               real_sars2_8596_all4             sarscov2_n8596         -m -v -p -r
               real_vch_gz_vcf_phylip           vcholerae_n100 (gzip)  -v -p
  holdout      h_real_vch_mono_vcf              vcholerae_n100         -b -v
               h_real_sars2_global_vcf_ref      sarscov2_global_n4637  -v -r
               h_real_mtb_phylip                mtb_n1788_snps         -p
               h_real_sars2_europe_mono         sarscov2_europe_n2012  -b -p -r
               h_real_salmonella_gz_vcf_phylip  see below              -v -p
The V. cholerae and M. tuberculosis alignments appear in both groups with different options, so different code paths
and outputs. gen_data.py also defines one simulated set per group, generated from fixed seeds when prepare.py runs.

Salmonella: prepare.py takes example_data/salmonella_serovars_core_genes.aln.gz from its clone of the snp-sites
repository at the pinned tag (Salmonella enterica serovar core-gene alignment from the snp-sites authors; GPL-3.0
repository content). It is not stored here.

1. vcholerae_n100.fa.xz (100 sequences x 4,047,835 columns; 3.4 MB xz, 412 MB plain)
   Vibrio cholerae O1, global, Snippy core-genome (whole-genome, reference-length) alignment, unmasked. Contains N
   (about 3.7%) and gaps. prepare.py also derives a gzip level-6 (mtime 0) copy for real_vch_gz_vcf_phylip.
   Source: Wilkinson E. "Global Vibrio Cholerae O1 dataset". Zenodo, 2025-09-01. doi:10.5281/zenodo.16909942, file
   cholera_global_n100_2025.09.01.fasta. Licence CC BY 4.0.

2. mtb_n1788_snps.fa.xz (1,788 sequences x 88,562 columns; 0.9 MB xz, 158 MB plain)
   Mycobacterium tuberculosis, a diverse multi-lineage set, Snippy v4.5.1 against H37Rv, filtered to variable sites.
   Source: Mixao V, Pinto M, Sobral D, Di Pasquale A, Gomes JP, Borges V. "Multiple Sequence Alignment of a diverse
   dataset with 1788 Mycobacterium tuberculosis isolates" (ReporTree benchmarking). Zenodo, 2023-03-27.
   doi:10.5281/zenodo.7772652, file MTb_original_align_profile.fasta. Licence CC BY 4.0.
   Underlying data: Walker TM et al. Lancet Infect Dis 2015;15:1193-1202, doi:10.1016/S1473-3099(15)00062-6; ENA
   PRJNA282721. ReporTree: Mixao V et al. Genome Med 2023;15:43.

3. sarscov2_europe_n2012.fa.xz (2,012 x 29,903; 0.17 MB xz, 60 MB plain)
4. sarscov2_global_n4637.fa.xz (4,637 x 29,903; 0.35 MB xz, 139 MB plain)
5. sarscov2_n8596.fa.xz (8,596 x 29,903; union of the seven regional builds; 0.66 MB xz, 257 MB plain)
   SARS-CoV-2 genomes aligned to Wuhan-Hu-1 by the Nextstrain ncov workflow (Nextclade alignment), as published in
   Nextstrain's open (GenBank-derived) builds, downloaded on 11 October 2026:
   https://data.nextstrain.org/files/ncov/open/{global,europe,africa,asia,north-america,oceania,south-america}/aligned.fasta.xz
   They contain N runs (about 3-4%), gaps and IUPAC ambiguity codes. Nextstrain regenerates these builds regularly, so
   the exact sequence sets are pinned by the files here, not by the URL.
   Licence and credit: Nextstrain open data are derived from GenBank/INSDC (no use restrictions). Cite Hadfield J et
   al. "Nextstrain: real-time tracking of pathogen evolution", Bioinformatics 2018;34:4121-4123, and Aksamentov I et
   al. "Nextclade", JOSS 2021;6(67):3773.

To download them again from the sources:
  bash data/tools/fetch_real_data.sh OUT_DIR
The V. cholerae and M. tuberculosis records are versioned, so those sequences are unchanged. The SARS-CoV-2 files will
differ from the sets above.
