# Harvested TCGA-LUAD CLDN4 / CXCL slice

- Matrix: UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` log2(RSEM norm_count+1).
- URL: `https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap/HiSeqV2.gz`
- Genes kept: CLDN4, CXCL9, CXCL10, CXCL13, CXCR3.
- Samples: primary tumors (`-01`) only; one row per 15-character barcode (mean if duplicated).
- ESTIMATE: official MDACC `lung_adenocarcinoma_RNAseqV2.txt` committed under `data/tcga_luad_cldn4_cxcl/`.
- Full HiSeqV2.gz is not committed (~29 MB). Re-download with `methods/tcga_luad_cldn4_cxcl/download.py` if the slim table is missing.
