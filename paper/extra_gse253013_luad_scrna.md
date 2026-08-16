# Extra scRNA figure · GSE253013 treatment-naïve LUAD

**Placement.** Extra LUAD TME scRNA panel. Do **not** replace or re-analyze user A3 (GSE207422 malignant TACSTD2 vs T/NK), which is taken as given.

## Honest cohort identity

GSE253013 (Sze, Xiang et al., *Cancer Res* 2024, PMID 38335304) is public LUAD scRNA-seq of **9 treatment-naïve** patients (256,379 cells after author QC; **89** GSM 10x lanes: tumor and adjacent non-tumor lung). The paper discusses immunotherapy resistance as motivation; the sequenced tissues were **not** on-treatment ICI. GEO and the Seurat `treatment` field are **None**. **No public MPR / pCR / R / NR labels.** This panel is an extra LUAD TME scRNA figure, not a response-stratified neoadjuvant ICI replication.

The series is lung tissue, not blood. After a real gene-name check of the GEO processed object (`GSE253013_all_luad_garnett_temp.rds.gz`), **TACSTD2 and CLDN4 are both present**.

## Results (patient unit)

Tumor lanes only. Malignant-like = marker epithelial minus normal-lung markers. Eligible: ≥10 malignant-like and ≥20 T/NK cells (all 9 tumor patients).

| Contrast | n | ρ | p |
|---|---:|---:|---:|
| Malignant-like TACSTD2 mean log1p UMI vs T/NK fraction | 9 | −0.72 | 0.030 |
| Malignant-like TACSTD2 mean log1p(CP10k) vs T/NK | 9 | −0.78 | 0.013 |
| Malignant-like TACSTD2 %pos vs T/NK | 9 | −0.83 | 0.0053 |
| All-epithelial TACSTD2 mean log1p vs T/NK | 9 | −0.12 | 0.77 |
| Author-Epithelial TACSTD2 vs author T-cell fraction | 9 | −0.58 | 0.099 |
| Malignant-like TACSTD2 vs T/NK, patients with ≥50 malignant-like cells | 5 | −0.60 | 0.28 |
| Malignant-like CLDN4 mean log1p UMI vs T/NK | 9 | −0.33 | 0.38 |
| Malignant-like CLDN4 mean log1p(CP10k) vs T/NK | 9 | −0.55 | 0.12 |
| All-epithelial CLDN4 mean log1p vs T/NK | 9 | −0.20 | 0.61 |
| Response / MPR / R vs TACSTD2 or CLDN4 | 0 | — | no public labels |

Both genes are epithelial-restricted (paired Wilcoxon, n=9): median %pos malignant-like TACSTD2 61.3 vs T/NK 0.59 (p=0.0020); CLDN4 67.4 vs 0.58 (p=0.0020).

## Suggested results text

In an independent public LUAD scRNA-seq series (GSE253013; 9 treatment-naïve patients, 89 10x lanes; not an ICI-response cohort), malignant-like *TACSTD2* was inversely associated with the per-patient tumor T/NK fraction (Spearman ρ=−0.72, p=0.030, n=9). The same direction held for log1p(CP10k) (ρ=−0.78, p=0.013) and percent-positive cells (ρ=−0.83, p=0.0053). All-epithelial *TACSTD2* was not associated with T/NK fraction (ρ=−0.12, p=0.77). *CLDN4* showed a weaker, non-significant inverse trend (ρ=−0.33, p=0.38). No major-pathologic-response or radiographic-response labels are released with this series. Three of nine tumors have fewer than 50 malignant-like cells; restricting to the five tumors with ≥50 such cells left the TACSTD2 trend but removed significance (ρ=−0.60, p=0.28, n=5). T/NK fractions can be inflated by CD45 sorting (author protocol: CD45+/CD45− from 6/9 tumors).

## Methods (paper methods)

Processed counts were taken from the GEO supplementary Seurat RDS (series matrix has no expression). A streaming XDR parser extracted a lineage/target gene panel plus author metadata (`methods/gse253013_*.py`). Lineages were assigned by argmax of mean log1p marker scores. Malignant-like epithelium was epithelial cells with near-zero normal-lung markers (*SFTPA2*, *AGER*, *SCGB1A1*, *SCGB3A1*, *TPPP3*). Per-patient tumor T/NK fraction used T (*CD3D*/*CD3E*/*CD2*) plus NK (*NKG7*/*GNLY*/*FGFBP2*). Library size for CP10k was Seurat `nCount_RNA`. Associations used Spearman rank correlation on patients. GSE207422 was not re-analyzed.

## Figure

`results/gse253013/fig_gse253013_tacstd2_cldn4_vs_tnk.png`

**Extra scRNA.** Malignant-like *TACSTD2* and *CLDN4* versus per-patient T/NK fraction in GSE253013 treatment-naïve LUAD (tumor lanes; 256,379 cells, 9 patients). No public ICI response / MPR labels. Spearman, patient unit: *TACSTD2* n=9, ρ=−0.72, p=0.030; *CLDN4* n=9, ρ=−0.33, p=0.38.
