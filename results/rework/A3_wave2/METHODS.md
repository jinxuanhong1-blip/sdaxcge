# Methods — A3 wave-2 (GSE207422)

## Data (public only)

- GEO `GSE207422` processed UMI matrix: 24,292 genes × **92,330** barcodes (`GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz`, 176 MB). User PPT said ~90,652; the deposited post-QC matrix is 92,330.
- Sample sheet: 15 rows (3 pre-treatment biopsies + 12 post-treatment resections). MPR n=4 including pCR P06; NMPR n=8. Residual tumor % is the pathology column, not a cell label.
- Author CopyKAT / epithelium RDS objects are **not public** (confirmed in wave-1 `results/rework/A3_annot/`).
- DRMref third-party Seurat labels (30,877 / 92,330 cells; 12 post-treatment samples) are reused as a T/NK and malignant sensitivity, not as author CopyKAT.

Raw GSA-Human HRA001033 FASTQ was not used.

## Lineage

Canonical marker scores on `log1p(CP10k)` (Hu Fig. 1B / Methods plus standard lung genes). Assigned lineage = argmax. T vs NK broken by CD3E when scores are close. This is a marker reconstruction, not the unpublished Seurat object.

## CopyKAT-like / inferCNV-like aneuploidy (primary malignant call)

Actual `copykat` / `infercnv` R packages were **not** run (no R/Bioconductor in this environment; full HMM inferCNV on ~11k epithelial cells is out of scope). Instead a window-smoothed expression-CNV score that follows the authors' stated rule:

1. Restrict scoring to epithelial cells; diploid reference = marker fibroblasts + endothelia (author: “stromal cells as normal reference”). If stromal n < 80, T + myeloid cells are added and that is recorded.
2. Keep genes with genomic coordinates (`gene_chr.tsv`, 15,320 symbols), mean UMI > 0.05 and ≥20 expressing cells among the keep-set.
3. `log1p(CP10k)`, subtract the **median of reference cells** per gene.
4. Order by chromosome, start. Sliding-window mean, window = 25 genes (CopyKAT default `win.size`), no smoothing across chromosome boundaries. ≥5 genes/chromosome.
5. Per-cell score = **sum of |smoothed values|** (author: “sum of calculated CNV for each gene per cell”).
6. **Primary malignant** = epithelial AND score > 95th percentile of the reference-cell scores.
7. Sensitivities (all reported): reference p90; reference mean+2 SD; 1-D 2-means cut on epithelial scores; p95 AND below the epithelial 75th percentile of a normal-lung marker score; DRMref `Malignant cells`; all epithelial; CNV-low (diploid-like) epithelial.

This is a CopyKAT-**like** score on the public UMI. It is not the unpublished `copykat_res.rds`.

## TACSTD2 metrics (malignant-only)

Per post-treatment sample, among cells called malignant by a given definition:

- mean `log1p(CP10k)`
- mean `log1p(UMI)`
- **%positive at UMI ≥1, ≥2, ≥3**

Samples with zero malignant cells contribute NaN and drop out of that test. A sensitivity requires ≥5 malignant cells.

## T/NK definitions

| Name | Rule | Denominator |
|---|---|---|
| `lineage_TNK` | assigned T or NK | all cells in the sample |
| `umi_CD3E_CD8A_NKG7` | CD3E≥1 OR CD8A≥1 OR NKG7≥1 | all cells |
| `umi_CD3E_or_NKG7notCD3E` | CD3E≥1 (T) or (NKG7≥1 and CD3E=0) (NK) | all cells |
| `drmref_TNK` | DRMref CD8+ T + CD4+ T + NK | DRMref-annotated cells in that sample |

Author cell-type IDs do not exist on GEO.

## Tests

- **NMPR vs MPR:** exact two-sided Wilcoxon by enumerating all ways to assign the observed number of MPR labels. pCR = MPR. Unit = patient/sample, never cell.
- **vs T/NK:** Spearman of the sample-level TACSTD2 metric vs each T/NK fraction.
- **Residual on epithelial fraction:** OLS residual of TACSTD2 ~ epithelial fraction, then Spearman vs T/NK. Also partial Spearman via rank residuals of both variables on epithelial fraction.

## Match / mismatch (pre-specified)

- NMPR>MPR: **MATCH** if Δ(NMPR−MPR)>0 and exact p<0.05; **PARTIAL_direction_NS** if Δ>0 and p≥0.05; **MISMATCH** if Δ≤0.
- vs T/NK: **MATCH** if ρ ∈ [−0.50, −0.40]; **NEAR** if ρ ∈ [−0.55, −0.35]; else **MISMATCH**.

User claim recovered only if a pre-specified primary test is MATCH on both halves.
