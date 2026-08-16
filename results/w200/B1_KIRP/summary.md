# CLDN4 vs TACSTD2 co-expression in TCGA-KIRP (honest ranking)

**Question:** Is `CLDN4` the top co-expression partner of `TACSTD2` among surface genes?

**Answer:** NO -- CLDN4 is NOT the single top surface-gene co-expression partner of TACSTD2.

- Cohort: TCGA-KIRP, 290 primary-tumour samples (sample-type ['01']).
- Surface-gene universe ranked: 2608 (10 zero-variance genes dropped; correlation undefined).
- Primary metric: Spearman correlation.
- Expression: UCSC Xena TCGA.KIRP.sampleMap/HiSeqV2 (log2 norm_count+1).
- Surfaceome: Bausch-Fluck et al. 2018 in-silico surfaceome (table S3).

| metric | value |
| --- | --- |
| CLDN4 Spearman rank | **#57 of 2608** (97.853th percentile) |
| CLDN4 Spearman rho | 0.435 (FDR q=2.68e-13) |
| CLDN4 Pearson rank | #80 of 2608 |
| CLDN4 Pearson rho | 0.386 |
| Actual #1 (Spearman) | MUC1 (rho=0.651) |

## Top 15 surface-gene partners of TACSTD2 (Spearman)

| rank | gene | spearman_rho | pearson_rho | FDR q |
| --- | --- | --- | --- | --- |
| 1 | MUC1 | 0.651 | 0.614 | 7.12e-33 |
| 2 | CLDN7 | 0.621 | 0.458 | 3.52e-29 |
| 3 | KCNS1 | 0.610 | 0.610 | 4.57e-28 |
| 4 | PROM1 | 0.609 | 0.474 | 4.57e-28 |
| 5 | MST1R | 0.608 | 0.599 | 4.57e-28 |
| 6 | GPR110 | 0.597 | 0.579 | 8.49e-27 |
| 7 | EFNB2 | 0.587 | 0.599 | 8.86e-26 |
| 8 | EPHA1 | 0.584 | 0.574 | 1.80e-25 |
| 9 | ITGB6 | 0.557 | 0.580 | 1.08e-22 |
| 10 | TPBG | 0.541 | 0.546 | 3.47e-21 |
| 11 | GJB3 | 0.541 | 0.482 | 3.47e-21 |
| 12 | ERBB2 | 0.538 | 0.524 | 6.06e-21 |
| 13 | CD9 | 0.536 | 0.483 | 1.02e-20 |
| 14 | EREG | 0.527 | 0.530 | 6.44e-20 |
| 15 | EMP1 | 0.523 | 0.491 | 1.34e-19 |

## Focus gene (`CLDN4`) — not in the top 15

| rank | gene | spearman_rho | pearson_rho | FDR q |
| --- | --- | --- | --- | --- |
| 57 | CLDN4 | 0.435 | 0.386 | 2.68e-13 |
