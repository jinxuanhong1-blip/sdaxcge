# Data sources (all public, open access)

Raw expression / purity matrices are NOT committed (see `.gitignore`).
Regenerate them with `scripts/align_tcga/download_data.sh`.

| File | Source | Description | Unit |
|------|--------|-------------|------|
| `LUAD.HiSeqV2.gz` | UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` | TCGA LUAD gene expression (Illumina HiSeq, RSEM) | log2(norm_count+1) |
| `LUSC.HiSeqV2.gz` | UCSC Xena `TCGA.LUSC.sampleMap/HiSeqV2` | TCGA LUSC gene expression (Illumina HiSeq, RSEM) | log2(norm_count+1) |
| `Aran_CPE_purity.xlsx` | Aran, Sirota & Butte, *Nat Commun* 2015 (`ncomms9971`), Supplementary Data 1 | Per-sample tumour purity: ESTIMATE, ABSOLUTE, LUMP, IHC, CPE | fraction 0–1 |
| `GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz` | GEO GSE253564 (Altorki et al. *Nat Commun* 2024) | Neoadjuvant durvalumab ± SBRT, **pre-treatment** FPKM | FPKM |
| `GSE248378_Durva_Post_FPKMs.txt.gz` | GEO GSE248378 (same trial) | **Post-treatment / resected** FPKM | FPKM |
| `SI_geneset.gmt` | Yoshihara et al. ESTIMATE R package `inst/extdata` | 141 stromal + 141 immune genes | — |
| `common_genes.txt` | same ESTIMATE package | 10,412 genes used by `filterCommonGenes` | — |
| `sample_estimate.gct` | same ESTIMATE package | Official sample scores for the Python-port self-test | — |

OncoSG LUAD 2020 expression + `PURITY` are pulled live from the public
cBioPortal API (`studyId=luad_oncosg_2020`); nothing is cached in `data/`.

Barcode matching (TCGA): expression columns are 15-char barcodes
(e.g. `TCGA-69-7978-01`); the purity table uses 16-char barcodes with a vial
letter, trimmed to 15 chars. Only primary tumours (suffix `-01`) are used.
