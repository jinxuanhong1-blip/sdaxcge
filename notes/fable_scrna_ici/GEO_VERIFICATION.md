# GEO verification — human ICI scRNA-seq slice

Verified 2026-08-16 against NCBI GEO (accession viewer + family SOFT files
from `ftp.ncbi.nlm.nih.gov`). Processed matrices only; nothing >2 GB.

## GSE207422 — neoadjuvant anti-PD-1 + chemotherapy NSCLC (MPR)

- Title: "Tumor microenvironment remodeling after neoadjuvant immunotherapy in
  non-small cell lung cancer revealed by single-cell RNA sequencing"
- Paper: Hu et al., *Genome Medicine* 2023, PMID 36869384, PMC10024472.
- Status: Public 2023-03-21. Platform GPL24676 (NovaSeq 6000, Homo sapiens).
- 15 scRNA samples (BD_immune01–15) from 15 patients: 3 pre-treatment biopsies
  + 12 post-treatment resections. Authors group pCR with MPR (post-tx MPR n=4,
  NMPR n=8). P01 pre-tx is NE.
- Contact: Peng Zhang, Shanghai Pulmonary Hospital. Raw data: GSA HRA001033
  (not on GEO).
- Files used (both << 2 GB):
  - `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz` (175.5 MB): dense genes ×
    cells UMI matrix, **92,330 cells**, 24,292 genes, barcodes
    `BD_immuneXX_<id>`.
  - `GSE207422_NSCLC_scRNAseq_metadata.xlsx` (11.1 KB): Patient, Resource,
    PD-1 antibody, chemo, Pathologic Response (MPR / NMPR / pCR / NE),
    Residual Tumor, RECIST.
- **No per-cell annotation is on GEO.** Hu et al. called malignant epithelium
  with **CopyKAT** (fibroblast + endothelial reference). Those IDs are not in
  the paper supplements or `Junjie-Hu/NSCLC-immunotherapy`.
- **Third-party labels used as the primary malignant set:** DRMref (Liu et al.
  *NAR* 2024, https://ccsm.uth.edu/DRMref/) Seurat objects
  `GSE207422_Tor` (8,690 cells) + `GSE207422_Sin` (22,187 cells) = 30,877
  post-treatment cells, 2,051 labeled `Malignant cells`. Marker-based 16-type
  annotation, **not** CopyKAT. All 30,877 barcodes match the GEO UMI matrix.
  Pre-treatment biopsies P01/P05/P08 are absent from DRMref.
- A marker-based epithelial-minus-normal-lung filter is kept as sensitivity
  (see WRITEUP). CopyKAT was not re-run.

## GSE205335 — lung cancer on ICI (RECIST)

- Title: "Single-cell transcriptome profiles of tumor tissues from lung cancer
  patients receiving immune checkpoint inhibitors."
- Status: Public 2024-11-05. Platforms GPL16791 / GPL24676.
- 33 GSM (GSM6210624–GSM6210656) from 26 patients: LN / liver / effusion /
  lung / bronchus tumors plus Normal Lung / Normal LN / Normal Brain.
- SOFT characteristics: patient, tissue, stage, cancer subtype (ADC / SQ /
  SCLC / NUT), RECIST (PR / SD / PD / NE), 10x 3' or 5'.
- Raw data: EGA EGAD00001008703 (not on GEO).
- Files used (both << 2 GB):
  - `GSE205335_Lung_IO_UMI_matrix.rds.gz` (499.5 MB download): **double-gzipped**;
    one `zcat` yields a readable RDS `dgCMatrix` of **33,714 genes × 96,505 cells**.
  - `GSE205335_Lung_IO_CellIdentity.txt.gz` (718.7 KB): author labels including
    `lineage.sub = Malignant cells` (28,512 cells) and `lineage.total = T/NK cells`.
- GSM → `orig.ident` reconstructed from title + platform (e.g. "P1006 EBUS_06"
  + 3' → `EBUS-06-3P`) and matched **33/33**
  (`results/fable_scrna_ici/gse205335_sample_table.tsv`).

## Size policy

Largest download: 523.7 MB (GSE205335 RDS). No FASTQ/BAM. No file >2 GB.
