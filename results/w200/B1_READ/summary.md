# CLDN4 vs TACSTD2 co-expression in TCGA-READ (honest ranking)

**Question:** Is `CLDN4` the top co-expression partner of `TACSTD2` among surface genes?

**Answer:** NO -- CLDN4 is not a TACSTD2 surface-gene partner in this cohort. Spearman rho is -0.007 (FDR q=9.92e-01); rank #1579 of 2618 is mid-pack / null, not a near-miss.

- Cohort: TCGA-READ, 94 primary-tumour samples (sample-type ['01']; matrix is 20530 genes × 105 samples).
- Surface-gene universe: 2618 surfaceome genes present in the matrix (anchor removed). 169 genes had undefined Spearman (zero variance across the cohort) and sort to the bottom.
- Primary metric: Spearman correlation.
- `TACSTD2` log2(norm_count+1): mean 8.15, range 3.58–13.72.
- `CLDN4` log2(norm_count+1): mean 14.04, range 12.36–15.20 (constitutively high; little dynamic range).

| metric | value |
| --- | --- |
| CLDN4 Spearman rank | **#1579 of 2618** (39.725th percentile) |
| CLDN4 Spearman rho | -0.007 (FDR q=9.92e-01) |
| CLDN4 Pearson rank | #951 of 2618 |
| CLDN4 Pearson rho | 0.060 |
| Actual #1 (Spearman) | PVRL4 (rho=0.604) |

## Top 15 surface-gene partners of TACSTD2 (Spearman)

| rank | gene | spearman_rho | pearson_rho | FDR q |
| --- | --- | --- | --- | --- |
| 1 | PVRL4 | 0.604 | 0.614 | 2.72e-07 |
| 2 | LYPD3 | 0.465 | 0.463 | 2.00e-03 |
| 3 | AMIGO2 | 0.464 | 0.462 | 2.00e-03 |
| 4 | TSPAN14 | 0.459 | 0.466 | 2.00e-03 |
| 5 | GPR87 | 0.441 | 0.412 | 3.70e-03 |
| 6 | ALPP | 0.440 | 0.515 | 3.70e-03 |
| 7 | MSLN | 0.436 | 0.408 | 4.02e-03 |
| 8 | SLC26A9 | 0.421 | 0.549 | 7.17e-03 |
| 9 | DRD1 | 0.399 | 0.403 | 1.41e-02 |
| 10 | GRIK2 | 0.397 | 0.323 | 1.41e-02 |
| 11 | LY6G6C | 0.396 | 0.365 | 1.41e-02 |
| 12 | MFI2 | 0.395 | 0.430 | 1.41e-02 |
| 13 | CD164L2 | 0.395 | 0.347 | 1.41e-02 |
| 14 | DPCR1 | 0.394 | 0.336 | 1.41e-02 |
| 15 | ALPPL2 | 0.391 | 0.484 | 1.50e-02 |
