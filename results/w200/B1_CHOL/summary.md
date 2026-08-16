# CLDN4 vs TACSTD2 co-expression in TCGA-CHOL (honest ranking)

**Question:** Is `CLDN4` the top co-expression partner of `TACSTD2` among surface genes?

**Answer:** NO -- CLDN4 is NOT the single top surface-gene co-expression partner of TACSTD2, and it is not in the FDR-significant set.

- Cohort: TCGA-CHOL, 36 primary-tumour samples (sample-type ['01']). Small n.
- Surface-gene universe: 2394 rankable surfaceome genes present in the matrix (anchor removed; 224 zero-variance genes dropped).
- Primary metric: Spearman correlation.

| metric | value |
| --- | --- |
| CLDN4 Spearman rank | **#153 of 2394** (93.651th percentile) |
| CLDN4 Spearman rho | 0.412 (FDR q=1.61e-01) |
| CLDN4 Pearson rank | #146 of 2394 |
| CLDN4 Pearson rho | 0.418 |
| Actual #1 (Spearman) | FZD6 (rho=0.744) |
| Surface-abundance rank, CLDN4 | **#49 of 2619** (median log2=12.50) |
| Surface-abundance rank, TACSTD2 | **#299 of 2619** (median log2=10.15) |

Pairwise TACSTD2–CLDN4 Spearman rho = 0.412 (bootstrap 95% CI 0.09–0.68, p=0.012). The interval is wide because n=36.

## Top 15 surface-gene partners of TACSTD2 (Spearman)

| rank | gene | spearman_rho | pearson_rho | FDR q |
| --- | --- | --- | --- | --- |
| 1 | FZD6 | 0.744 | 0.696 | 4.72e-04 |
| 2 | MYADM | 0.709 | 0.718 | 1.55e-03 |
| 3 | QSOX1 | 0.662 | 0.657 | 8.35e-03 |
| 4 | EPHB3 | 0.656 | 0.639 | 8.35e-03 |
| 5 | ADAM8 | 0.633 | 0.571 | 1.64e-02 |
| 6 | ENTPD3 | 0.626 | 0.681 | 1.64e-02 |
| 7 | SORT1 | 0.615 | 0.581 | 1.64e-02 |
| 8 | ITGA2 | 0.613 | 0.608 | 1.64e-02 |
| 9 | ITGA11 | 0.610 | 0.609 | 1.64e-02 |
| 10 | SLC41A1 | 0.608 | 0.530 | 1.64e-02 |
| 11 | MANSC1 | 0.607 | 0.600 | 1.64e-02 |
| 12 | SLC9A6 | 0.606 | 0.606 | 1.64e-02 |
| 13 | PVRL4 | 0.600 | 0.494 | 1.79e-02 |
| 14 | CD47 | 0.599 | 0.619 | 1.79e-02 |
| 15 | IGSF9 | 0.596 | 0.538 | 1.79e-02 |
