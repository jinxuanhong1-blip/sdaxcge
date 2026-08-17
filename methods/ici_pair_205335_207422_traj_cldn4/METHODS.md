# Methods — ICI-pair CLDN4-only trajectory (GSE205335 + GSE207422)

## Scope

Public processed UMI only. Two ICI scRNA sets with outcome labels. Patient is
the inferential unit. Cell-level statistics are not reported as claims.

Excluded by design: GSE148071, GSE179994, any 7-pool / joint Harmony object,
TACSTD2∩CLDN4 dual-high, EGA FASTQ, GSA-Human HRA001033.

## Matrices

- **GSE205335** (Ahn/Lee, eLife 2024): `GSE205335_Lung_IO_UMI_matrix.rds.gz`
  + author `GSE205335_Lung_IO_CellIdentity.txt.gz` + GEO SOFT characteristics.
  Epithelium = `lineage.total == Epithelial cells` on non-normal tissues.
  Malignant = `lineage.sub == Malignant cells`. Normal LN / lung / brain dropped.
- **GSE207422** (Hu, *Genome Med* 2023): `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz`
  + sample xlsx. Author CopyKAT / DRMref barcodes are not public. Epithelium =
  Hu canonical marker argmax (EPCAM/KRT vs T/NK/B/myeloid/…). A3-malignant-like
  = epithelial AND zero UMI of SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3. Leftover
  epithelium (those normal-lung genes > 0) is kept for the DPT root.
  Pre-treatment biopsies (P01/P05/P08) are out of the graph.

Per-patient cap: 400 epithelial cells (seed 1 / 3), leftover cells protected.

## Trajectory (per dataset)

QC: ≥200 genes, ≥500 UMI, genes in ≥10 cells. log1p(CP10k). 3000 HVG
(seurat_v3 if available), 30 PCs, 30 neighbors, Leiden 0.6, PAGA, diffusion
map, DPT (10 DCs).

**No Harmony.** Each dataset is its own graph. Platforms differ (GSE205335
mixed 3′/5′ 10x biopsies/effusions; GSE207422 10x post-resection).

Root: leftover epithelium **in the giant PAGA component**, CLDN4 at or below the
leftover second tertile, AT2 near the 80th percentile (or max AT2 if that
percentile is 0). Never CLDN4-high. Leftover cells are protected at extract
time but still capped (one leftover-heavy patient cannot dominate).

Barrier/keratin score: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR. **CLDN4 out.**

## Outcomes

- GSE205335 primary: RECIST **PR vs PD**. SD and NE excluded. R vs NR
  (PR vs SD+PD) and ADC+SQ-only are sensitivity.
- GSE207422 primary: **MPR vs NMPR** on the 12 post-treatment patients
  (pCR P06 = MPR).
- Stacked binary: PR+MPR vs PD+NMPR. Raw DPT is not pooled. Within-dataset
  rank and z are used. Fisher-z DerSimonian–Laird of the two per-dataset
  CLDN4–DPT Spearman ρ is descriptive (k=2).

Eligible Spearman n: ≥10 epithelial cells on the analysis object. Malignant-only
sensitivity: ≥10 malignant / A3-malignant-like cells.

## Software

scanpy PAGA/DPT (Wolf 2019; Haghverdi 2016). Slingshot is not required.
