# A11 nectin genes vs TACSTD2-high public lung

> Honest report. Numbers are computed from public matrices. 
> TCGA patients were not ICI-treated. Correlations are not causation. 
> CLDN4 is a junction positive control, not a nectin.

## TL;DR

**Honest answer: not a nectin-family class effect. NECTIN4 yes; NECTIN1 moderate; NECTIN2/3 and PVR no.**

The pre-specified ≥3/4 partial-correlation rule is `CLASS_EFFECT_SUPPORTED` (3/4 primary nectins are purity+histology-adjusted positive with TACSTD2 and none are significantly negative. A nectin-family class effect is supported on this public lung slice.) That rule is too easy to trip: NECTIN2 is only WEAK_POSITIVE (pooled partial ρ = 0.15) and **fails the TACSTD2 Q4 vs Q1 test** (Δmedian ≈ 0, p = 0.51). NECTIN3 is LUAD-only and labeled NULL. The named A11 question is TACSTD2-**high** public lung, not a weak pooled ρ.

Primary cohort: TCGA LUAD n=516 + LUSC n=501 primary tumors (one sample/patient); 996 / 1017 have ABSOLUTE purity and enter the primary partial correlation.

## Pre-specified design

- Genes: NECTIN1, NECTIN2, NECTIN3, NECTIN4 (primary family). PVR reported separately. CLDN4 = positive control.
- Expression: Xena GDC STAR TPM, log2(TPM+1), GENCODE v36.
- Primary statistic: rank-partial Spearman of TACSTD2 vs each gene, covariates = ABSOLUTE purity + histology in the pooled analysis.
- TACSTD2-high: Q4 vs Q1 Mann-Whitney U (median split is sensitivity).
- FDR: BH within each cohort across the 4 nectins.
- Class effect requires ≥3/4 nectins WEAK_POSITIVE or ASSOCIATED_POSITIVE and none negative.

## Primary result — pooled partial Spearman (purity + histology)

| gene | role | partial ρ | partial p | FDR (family) | label | LUAD partial ρ | LUSC partial ρ |
|---|---|---:|---:|---:|---|---:|---:|
| NECTIN1 | primary_nectin | 0.302 | 2.06e-22 | 4.13e-22 | ASSOCIATED_POSITIVE | 0.229 | 0.321 |
| NECTIN2 | primary_nectin | 0.147 | 3.08e-06 | 4.11e-06 | WEAK_POSITIVE | 0.221 | 0.088 |
| NECTIN3 | primary_nectin | 0.092 | 3.84e-03 | 3.84e-03 | NULL | 0.239 | -0.028 |
| NECTIN4 | primary_nectin | 0.522 | 1.23e-70 | 4.92e-70 | ASSOCIATED_POSITIVE | 0.474 | 0.531 |
| PVR | nectin_like | 0.062 | 4.89e-02 | 4.89e-02 | NULL | -0.003 | 0.127 |
| CLDN4 | positive_control | 0.436 | 2.07e-47 | 6.22e-47 | ASSOCIATED_POSITIVE | 0.525 | 0.388 |

## TACSTD2-high vs TACSTD2-low (pooled Q4 vs Q1)

| gene | median Q4 | median Q1 | Δmedian | MWU p |
|---|---:|---:|---:|---:|
| NECTIN1 | 7.619 | 4.546 | 3.073 | 1.76e-26 |
| NECTIN2 | 6.555 | 6.575 | -0.020 | 5.09e-01 |
| NECTIN3 | 2.514 | 2.170 | 0.344 | 1.51e-01 |
| NECTIN4 | 6.839 | 5.104 | 1.735 | 2.00e-47 |
| PVR | 4.572 | 4.537 | 0.035 | 2.18e-01 |
| CLDN4 | 8.029 | 7.204 | 0.825 | 3.61e-18 |

## Honest reading of each gene

- **NECTIN4**: the only robust TACSTD2-high partner. Pooled partial ρ = 0.52 (LUAD 0.47, LUSC 0.53). Q4 vs Q1 Δmedian = +1.73 (p = 2e-47). DepMap NSCLC ρ = 0.73. Strength is in the same range as the CLDN4 control (partial ρ = 0.44).
- **NECTIN1**: real moderate co-expression in both histologies (LUAD 0.23, LUSC 0.32). The pooled Q4 vs Q1 Δmedian = +3.07 is **inflated by histology mix**: TACSTD2 Q4 is LUSC-heavy and LUSC NECTIN1 baseline is ~3.5 log2 higher than LUAD. Within-histology Q4–Q1 deltas are +0.59 (LUAD) and +0.99 (LUSC). Use the partial ρ (0.30), not the pooled boxplot, for NECTIN1.
- **NECTIN2**: LUAD-only (0.22); LUSC null (0.09, FDR 0.07). Pooled Q4 vs Q1 is null. Do not call this a TACSTD2-high nectin.
- **NECTIN3**: LUAD 0.24, LUSC −0.03. Pooled |ρ| < 0.10. Tumor < adjacent normal in both histologies. Not a TACSTD2-high partner.
- **PVR**: null in tumors (pooled partial ρ = 0.06).
- **CLDN4 control**: recovered (partial ρ = 0.44). The pipeline can see a junction association; that does not make NECTIN2/3 positive.

## What this does **not** show

- It does not show that nectins cause TACSTD2-high tumors, or the reverse.
- It does not show ICI response, ADC response, or protein-level co-expression.
- A positive CLDN4 control only says the pipeline can recover a known junction association; it does not rescue a weak nectin class effect.
- Histology can disagree. If LUAD and LUSC labels differ, the pooled number is not a license to ignore the split.
- DepMap lung lines are models, not tumors; they are sensitivity only.

## Sensitivity — DepMap 24Q4 lung cell lines (no purity adjustment)

| cohort | gene | Spearman ρ | p | n | label |
|---|---|---:|---:|---:|---|
| depmap_lung_cell_lines | NECTIN1 | -0.045 | 5.10e-01 | 214 | NULL |
| depmap_lung_cell_lines | NECTIN2 | 0.270 | 6.11e-05 | 214 | ASSOCIATED_POSITIVE |
| depmap_lung_cell_lines | NECTIN3 | 0.086 | 2.12e-01 | 214 | NULL |
| depmap_lung_cell_lines | NECTIN4 | 0.689 | 1.61e-31 | 214 | ASSOCIATED_POSITIVE |
| depmap_lung_cell_lines | PVR | 0.291 | 1.48e-05 | 214 | ASSOCIATED_POSITIVE |
| depmap_lung_cell_lines | CLDN4 | 0.607 | 5.64e-23 | 214 | ASSOCIATED_POSITIVE |
| depmap_NSCLC_cell_lines | NECTIN1 | 0.223 | 7.43e-03 | 143 | ASSOCIATED_POSITIVE |
| depmap_NSCLC_cell_lines | NECTIN2 | -0.031 | 7.16e-01 | 143 | NULL |
| depmap_NSCLC_cell_lines | NECTIN3 | -0.097 | 2.50e-01 | 143 | NULL |
| depmap_NSCLC_cell_lines | NECTIN4 | 0.728 | 6.41e-25 | 143 | ASSOCIATED_POSITIVE |
| depmap_NSCLC_cell_lines | PVR | 0.111 | 1.86e-01 | 143 | NULL |
| depmap_NSCLC_cell_lines | CLDN4 | 0.667 | 8.81e-20 | 143 | ASSOCIATED_POSITIVE |

## Tumor vs adjacent normal (context, not the A11 claim)

| cohort | gene | median tumor | median normal | Δ | unpaired p |
|---|---|---:|---:|---:|---:|
| LUAD | TACSTD2 | 9.013 | 8.637 | 0.376 | 7.97e-03 |
| LUAD | NECTIN1 | 4.206 | 3.067 | 1.140 | 5.89e-19 |
| LUAD | NECTIN2 | 6.934 | 6.640 | 0.294 | 1.38e-05 |
| LUAD | NECTIN3 | 2.567 | 3.435 | -0.868 | 4.84e-12 |
| LUAD | NECTIN4 | 5.914 | 3.793 | 2.121 | 1.64e-27 |
| LUAD | PVR | 4.610 | 4.928 | -0.318 | 2.52e-03 |
| LUAD | CLDN4 | 8.181 | 7.016 | 1.164 | 6.59e-18 |
| LUSC | TACSTD2 | 9.575 | 8.699 | 0.876 | 4.83e-09 |
| LUSC | NECTIN1 | 7.732 | 3.391 | 4.341 | 1.17e-29 |
| LUSC | NECTIN2 | 6.359 | 6.629 | -0.270 | 5.26e-04 |
| LUSC | NECTIN3 | 2.201 | 3.340 | -1.138 | 1.24e-10 |
| LUSC | NECTIN4 | 6.518 | 3.631 | 2.887 | 1.35e-27 |
| LUSC | PVR | 4.495 | 4.957 | -0.462 | 3.72e-06 |
| LUSC | CLDN4 | 7.206 | 6.808 | 0.398 | 1.36e-02 |

## Files

- `correlations.csv` — marginal and partial Spearman
- `tacstd2_high_vs_low.csv` — Q4 vs Q1 and median split
- `tumor_vs_normal.csv` — adjacent-normal context
- `summary.json` — machine-readable verdict
- `fig1_partial_rho_heatmap.png` / `fig2_scatter_*.png` / `fig3_*boxplots.png`

