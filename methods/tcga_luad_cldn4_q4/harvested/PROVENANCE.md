# Provenance — TCGA-LUAD CLDN4 Q4 extra

| File | n | Source |
|---|---:|---|
| `hiseqv2_cldn4_cd8_cd274.tsv` | 515 primaries | Extracted from UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2.gz` (CLDN4, CD8A, CD8B, CD274 only) |
| `../../../../data/tcga_luad_cldn4_q4/MDACC_estimate_LUAD_RNAseqV2.txt` | 515 primaries + 2 `−02` | https://bioinformatics.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt |

HiSeqV2 URL: `https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap/HiSeqV2.gz`

The full 29 MB HiSeqV2 matrix is not committed. The slim four-gene table is.

TJ7 vs CD8 is taken as given from PR #69 / #90 and is not re-scored here.
