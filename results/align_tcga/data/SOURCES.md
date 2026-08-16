# Data sources (all public, open access)

Raw matrices are NOT committed to git (see `.gitignore`); regenerate them with
`scripts/align_tcga/download_data.sh`.

| File | Source | Description | Unit |
|------|--------|-------------|------|
| `LUAD.HiSeqV2.gz` | UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` | TCGA LUAD gene expression (Illumina HiSeq, RSEM) | log2(norm_count+1) |
| `LUSC.HiSeqV2.gz` | UCSC Xena `TCGA.LUSC.sampleMap/HiSeqV2` | TCGA LUSC gene expression (Illumina HiSeq, RSEM) | log2(norm_count+1) |
| `Aran_CPE_purity.xlsx` | Aran, Sirota & Butte, *Nat Commun* 2015 (`ncomms9971`), Supplementary Data 1 | Per-sample tumour purity: ESTIMATE, ABSOLUTE, LUMP, IHC, CPE | fraction 0–1 |

Barcode matching: expression columns are 15-char barcodes (e.g. `TCGA-69-7978-01`);
the purity table uses 16-char barcodes with a vial letter (e.g. `TCGA-69-7978-01A`),
trimmed to 15 chars for the join. Only primary tumours (barcode suffix `-01`) are used.
