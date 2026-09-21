# TROP2–CLDN4 protein on Gygi CCLE: the histology maximum is 0.85 (n=17)

**Maximum in the locked histology / NSCLC / partial grid: Spearman partial ρ = 0.846, n = 17.**

That estimate is TACSTD2 (TROP2) versus CLDN4 protein in DepMap 24Q4 NSCLC lines whose model was derived from a primary tumor, after a rank partial correlation for Oncotree subtype. It is not n = 118. Nothing in this grid has n = 118.

The number that rounds to the prior **0.69** is still the unadjusted S1 Lung complete-case correlation: **ρ = 0.693, n = 45**.

| Analysis | n | Spearman ρ | 95% CI | p | Same-size null |
|---|---:|---:|---|---:|---:|
| S1 Lung, unadjusted (the 0.69) | 45 | 0.693 | 0.48–0.82 | 1.3×10⁻⁷ | — |
| NSCLC | 35 | 0.727 | 0.51–0.84 | 7.7×10⁻⁷ | 0.21 |
| LUAD | 22 | 0.782 | 0.52–0.90 | 1.7×10⁻⁵ | 0.15 |
| NSCLC, metastatic | 18 | 0.600 | 0.16–0.85 | 0.0085 | 0.75 |
| NSCLC, primary tumor | 17 | 0.821 | 0.52–0.94 | 5.3×10⁻⁵ | 0.081 |
| **NSCLC primary \| Oncotree subtype** | **17** | **0.846** | **0.55–0.96** | **7.1×10⁻⁵** | — |
| LUAD, primary tumor | 11 | 0.873 | 0.52–0.99 | 4.5×10⁻⁴ | 0.056 |
| S1 Lung, CLDN4 peptides ≥ 2 | 25 | 0.829 | 0.60–0.93 | 3.0×10⁻⁷ | 0.023 |

LUAD primary (0.873) is higher and is reported. It is not the maximum used here because the rule requires n ≥ 15. LUSC (0.750, n = 7) and SCLC (0.643, n = 7) are below that line and are not significant at 0.05.

Benjamini–Hochberg q for the 0.846 partial, among the 20 grid rows that meet the n / df / p rule, is 9.5×10⁻⁵. The bootstrap interval still covers 0.69. A random 17-line subset of the same 45 lung complete cases reaches ρ ≥ 0.821 about 8% of the time, so the primary-NSCLC point estimate is the top of this grid and is not a rare draw under subset noise.

Leave-one-out on the 17 primary NSCLC lines keeps the unadjusted ρ between 0.785 and 0.856. Dropping any one line does not create the correlation.

## What was held fixed

Protein matrix: Nusinow et al., *Cell* 2020, `protein_quant_current_normalized.csv.gz` (SHA256 `b72a9ff3…1be2a7ed`, same file as the earlier n=45 / ρ=0.69 recompute). Sample info is Gygi Table S1. Histology and primary versus metastatic labels are DepMap Public 24Q4 `Model.csv`.

S1 Tissue = Lung is 77 lines. CLDN4 is quantified in 45 of them. TACSTD2 is quantified in all 77. Pairwise complete cases are the analysis n. No imputation.

The maximum is the largest Spearman in a list written in `scripts/depmap_trop2_cldn4_maxrho/analyze.py` before the sort:

- Named Oncotree slices of those 45 lines (lung, NSCLC, LUAD, LUSC, SCLC, NSCLC without LUSC, NSCLC without LUAD).
- The DepMap `PrimaryOrMetastasis` split inside lung, NSCLC, and LUAD. Both sides are in the table. Metastatic NSCLC is the weaker side (0.600).
- Partials for Oncotree subtype, Oncotree primary disease, and primary versus metastatic. Method: rank the variables, residualize by least squares, then Pearson. Rare subtype levels with fewer than 3 lines are pooled as Other. Residual df = n − k − 2 must be at least 10.

Epithelial partials are computed and cannot win. On all 45 lung lines, adjusting for EPCAM lowers ρ from 0.693 to 0.553, and adjusting for CDH1 lowers it to 0.389. The co-expression is partly the shared epithelial program. VIM moves it only to 0.640.

Holding histology fixed inside all 45 lines does not raise ρ (subtype partial 0.703). The increase is the primary-NSCLC slice (0.821), plus 0.025 from the subtype partial (0.846).

## Peptide depth, separate from the histology maximum

CLDN4 in this matrix is a low-peptide protein. In the 45 lung complete cases the plex-level peptide count is 1, 2, 3, or 5. Restricting to plexes with at least 2 CLDN4 peptides gives **ρ = 0.829, n = 25**. Only 2.3% of random 25-line subsets of the same 45 reach a ρ that high. A subtype partial inside that set stays at 0.830.

That filter is not selecting high CLDN4. Peptide count versus TACSTD2 abundance is ρ = 0.019 (p = 0.90), and versus CLDN4 abundance is ρ = 0.017 (p = 0.91). Lines quantified on a single peptide are noisier, and removing them raises the correlation. n = 11 at the ≥3-peptide cut (ρ = 0.836) is below the n ≥ 15 rule, and its same-size null p is 0.12.

## 中文

公开 Gygi/CCLE 蛋白表上，TROP2（TACSTD2）对 CLDN4 的组织学网格最大值是 **原发瘤来源 NSCLC、按 Oncotree 亚型做偏相关，ρ = 0.846，n = 17**。原先能四舍五入到 **0.69** 的仍是 S1 肺系成对完整病例 **ρ = 0.693，n = 45**。S1 肺系一共 77 株，没有哪一条预设定切片的 n 是 118。LUAD 原发瘤 ρ = 0.873 但 n = 11，不够 n ≥ 15 的入选线。用 EPCAM 或 CDH1 做偏相关会把 ρ 降到 0.55 或 0.39。CLDN4 肽段数 ≥ 2 的灵敏度分析是 ρ = 0.829，n = 25。没有填补，也没有为了抬高 ρ 删细胞系。

## Reproduce

```bash
python3 scripts/depmap_trop2_cldn4_maxrho/download.py
python3 scripts/depmap_trop2_cldn4_maxrho/analyze.py
```
