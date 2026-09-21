# Methods

Concordant-4 only: GSE123902 (13 donors), GSE131907 (21 samples), GSE205335 (22 patients), GSE189357 (9 patients). Not GSE148071, GSE127465, GSE154826, GSE200563, or E-MTAB-13526. The unit list, cell counts, and `frac_tnk = n_tnk / n_cells` are the locked table in `data/locked_units.tsv` (n = 65). Cell counts are not n.

Compartments are the same as the locked Seurat run.

| Cohort | Malignant | T/NK |
| --- | --- | --- |
| GSE123902, GSE189357 | (EPCAM or KRT8 or KRT18 or KRT19) > 0 and PTPRC = 0 | (CD3D or CD3E or CD8A or NKG7 or GNLY or KLRD1) > 0, and not malignant |
| GSE131907 | author `Cell_subtype == Malignant cells` | author `Cell_type` in {T lymphocytes, NK cells} |
| GSE205335 | author `lineage.sub == Malignant cells` | author `lineage.total == T/NK cells` |

GSE123902 uses the 13 tumor captures only. GSE205335 drops captures whose tissue starts with "Normal". GSE131907 keeps the 21 tumor-bearing samples already in the locked table. Expression is the full malignant-cell count vector, not a UMAP cap.

## Endpoint

Within each cohort, Spearman of the malignant score versus `frac_tnk`. The four cohort coefficients are pooled with DerSimonian–Laird on Fisher z, the same estimator as the locked CLDN4 result. A specification is kept only when all four cohorts have a finite Spearman and N = 65.

The run stops unless malignant CLDN4 % of cells with UMI > 0 matches the locked patient table (max absolute difference < 10⁻⁴) and the pooled ρ matches −0.5311678045689989.

## What was searched

ELF3 alone (28 scores): % of malignant cells with UMI ≥ k for k in {1, 2, 3, 4, 5, 6, 8, 10, 12, 15}; % with CP10k ≥ t for t in {0.25, 0.5, 1, 2, 3, 5, 8, 10, 15}; % at or above the cohort quantile q of CP10k for q in {0.50, 0.60, 0.70, 0.75, 0.80, 0.90}; mean log1p(UMI); mean log1p(CP10k); mean log1p among detected cells; mean log1p in the within-unit upper quartile.

Four-gene module. The core is ELF3, TACSTD2, and CLDN4. The pre-specified fourth gene is CLDN7. The same score families were repeated after replacing CLDN7 with each gene that is present in every locked unit: CDH1, CDH3, CLDN1, CLDN3, EPCAM, F11R, GRHL1, KLF5, KRT7, KRT8, KRT18, KRT19, MUC1, OCLN, TJP1. GRHL2, SPDEF, and OVOL2 are missing from at least one unit and were not used as the fourth gene. Module summaries, all computed inside cohort, are: mean of z-scores, mean of ranks, first principal component sign-aligned to ELF3, the % of malignant cells in which all four genes pass the cut, and the % in which at least three pass the cut.

A module row is dropped when any member is constant inside any cohort. A constant member is removed by z-scoring, so the row would no longer be a four-gene score in that cohort. This removes cohort-median calls whose threshold is zero because more than half of the malignant cells are exact zeros (TACSTD2 and CLDN7 in parts of this matrix).

## What is reported

Two concordant maxima, both required to be negative in all four cohorts:

1. ELF3. The marginal score with the largest |ρ|.
2. The four-gene module. The pre-specified CLDN7 membership, and separately the largest |ρ| across fourth-gene swaps.

The unconstrained maximum (a cohort is allowed to be positive) is stored in `results/tables/winners.tsv`. It is not the headline when it differs from the concordant maximum. Stacked within-cohort Q4 versus Q1 rank-biserial correlations are companions. They are not the quantity that was maximized.

p-values on a searched maximum are descriptive. The grid has 1,030 specifications (28 ELF3, 1,002 module rows after the constant-member filter).

## Software

Python (numpy, scipy, pandas, matplotlib). GSE205335 is a double-gzipped dgCMatrix and is subset with R (`Matrix`). R does not run the correlations.
