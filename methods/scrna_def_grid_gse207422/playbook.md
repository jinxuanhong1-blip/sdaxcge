# Methods — GSE207422 public-UMI definition grid

ADDITIVE slice. Public GEO UMI only. Author CopyKAT / epithelium RDS barcodes are not on GEO, GitHub, TISCH2, or CELLxGENE; that object is not reconstructed. User A3 slide numbers (malignant TACSTD2 higher in NMPR than MPR; per-patient Spearman vs T/NK fraction ρ ≈ −0.40 to −0.50) are taken as given and are not re-argued. This playbook defines an exhaustive public grid and reports every honest ρ / p / n.

## Data

- GEO `GSE207422` processed UMI: `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz` (24,292 genes × 92,330 barcodes). User PPT ~90,652; the deposited matrix is 92,330. No public 90,652 subset is used.
- Sample sheet: `GSE207422_NSCLC_scRNAseq_metadata.xlsx` (15 rows). Barcode prefix `BD_immuneXX`.
- Primary unit: **12 post-treatment surgery samples**. MPR n=4 including pCR P06 (paper Fig. 1); NMPR n=8. Three pre-treatment biopsies are excluded from NMPR vs MPR and from the Spearman grid.
- DRMref third-party labels (Liu et al., *NAR* 2024): 30,877 / 92,330 post-treatment barcodes, 16 types, 2,051 `Malignant cells`. Marker-based, not Hu CopyKAT. Reused from the public DRMref Seurat `meta.data` already extracted in this repo.
- Raw GSA-Human HRA001033 FASTQ is not used.

## Malignant definitions

Author CopyKAT IDs: **not present**. Row is omitted from numeric cells and noted on the figure.

| Name | Rule |
|---|---|
| `hu_markers_mean` | Lineage = argmax mean `log1p(CP10k)` of Hu-style canonical markers (T, NK, B, Plasma, Myeloid, Neutrophil, Mast, Epithelial, Fibroblast, Endothelial). T vs NK broken by CD3E when scores are close. Malignant = Epithelial AND normal-lung score below the epithelial 60th percentile. Normal-lung panel: SFTPA1/2, SFTPB, SFTPC, AGER, SCGB1A1, SCGB3A1/2, TPPP3, FOXJ1, CAPS. |
| `hu_markers_pctpos` | Same lineages, but scores are the fraction of panel genes with UMI > 0. Malignant = Epithelial AND zero UMI on the normal-lung panel. |
| `epcam_krt_pos` | EPCAM UMI ≥ 1 AND (KRT8 ≥ 1 OR KRT18 ≥ 1 OR KRT19 ≥ 1). No normal-lung subtract. |
| `epcam_krt_mean` | Mean `log1p(CP10k)` of {EPCAM, KRT8, KRT18, KRT19} ≥ the median of Hu-mean epithelial cells, with a floor of 0.40, and PTPRC UMI < 1. (An all-cell percentile is near zero because the module is zero-inflated.) |
| `copykat_stromal_p95` | Window-smoothed expression CNV (not the R `copykat` package). Diploid reference = marker fibroblasts + endothelia. `log1p(CP10k)`, subtract **median** of reference cells, order by chromosome/start, sliding window = 25 genes, no cross-chromosome smoothing. Score = sum \|smoothed\|. Malignant = Hu-mean epithelial AND score > stromal 95th percentile. |
| `infercnv_stromal_p95` | Same genes/window/reference. Subtract **mean** of reference cells (inferCNV-style). Score = mean \|smoothed\|. Malignant = Hu-mean epithelial AND score > stromal 95th percentile. Not the R `infercnv` HMM. |
| `drmref_public` | DRMref `Malignant cells` on the 30,877 annotated barcodes. Unannotated cells are not malignant. |
| `drmref_like` | Rebuild: argmax of 12 DRMref-like type modules (Malignant / CD8+ T / CD4+ T / NK / B / Plasma / Mono-Macro / Neutrophils / pDCs / Mast / Fibroblasts / Endothelial). Malignant = type `Malignant cells` AND normal-lung score < 0.25. This is a marker rebuild, not the DRMref RDS. |

Hu epithelial panel: EPCAM, KRT8, KRT18, KRT19, KRT7, CDH1. T: CD3D, CD3E, CD3G, CD2, TRAC. NK: NKG7, GNLY, KLRD1, KLRF1, NCR1. Remaining lineages as in `scripts/analyze.py`.

## Immune definitions

Fractions use **all cells in the sample** as the denominator (patient-level composition).

| Name | Rule |
|---|---|
| `T_NK` | Hu-mean lineage T or NK |
| `CD8_only` | lineage T AND (CD8A ≥ 1 OR CD8B ≥ 1) |
| `CD8_NK` | `CD8_only` OR lineage NK |
| `CXCL13_pos` | (T or NK) AND CXCL13 UMI ≥ 1 |
| `cyto_high_T` | lineage T AND cytotoxicity module (`GZMB`, `GZMA`, `PRF1`, `IFNG`, `NKG7`) ≥ the 75th percentile of all T cells (pooled) |

## TACSTD2 scores (inside the malignant set of that definition)

| Name | Rule |
|---|---|
| `mean_log1p` | mean `log1p(CP10k)` of TACSTD2 |
| `pct_pos` | fraction of malignant cells with TACSTD2 UMI > 0 |
| `ucell_module` | mean UCell of {TACSTD2, CLDN4, EPCAM}. Python implementation of Andreatta & Carmona 2021: competition rank = 1 + n genes with strictly higher UMI; ranks > `maxRank=1500` set to 1501; score = 1 − U / (n_sig × maxRank), clipped to [0, 1]. |

Samples with zero malignant cells under a definition contribute NaN and drop from that test. n is the number of remaining post-treatment patients.

## Tests

- **vs immune:** Spearman of the sample-level TACSTD2 score vs each immune fraction. `scipy.stats.spearmanr`. Report ρ, p, n.
- **NMPR vs MPR:** exact two-sided Wilcoxon by enumerating all C(n, n_NMPR) label assignments. One-sided p for NMPR mean ≥ observed is also stored. pCR P06 = MPR.
- Highlight: ρ ≤ −0.35, or NMPR−MPR Δ > 0. These are flags on the public grid, not a retuning of the slide.

## Outputs

All under `methods/scrna_def_grid_gse207422/`:

- `grid_spearman.tsv`, `grid_nmpr_mpr.tsv` — full grid
- `highlight_rho_le_neg035.tsv`, `highlight_nmpr_gt_mpr.tsv`
- `per_sample.tsv`, `sample_metadata.tsv`, `cell_calls.tsv.gz`
- `fig_def_grid.png` — extra figure of the grid
- `REPORT.md` — every ρ / p / n, with highlighted rows
- `summary.json`, `cnv_info.json`
