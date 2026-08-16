# CLDN4 vs TACSTD2 co-expression in TCGA-STAD (honest ranking)

**Question:** Is `CLDN4` the top co-expression partner of `TACSTD2` among surface genes?

**Answer:** NO -- CLDN4 is NOT the single top surface-gene co-expression partner of TACSTD2.

- Cohort: TCGA-STAD, 415 primary-tumour samples (sample-type ['01']).
- Surface-gene universe: 2618 surfaceome genes present in the matrix (anchor removed). 1 surfaceome gene(s) had zero variance and were left unranked.
- Ranked genes (finite Spearman ρ): 2617.
- Primary metric: Spearman correlation.

| metric | value |
| --- | --- |
| CLDN4 Spearman rank | **#33 of 2617** (98.777th percentile) |
| CLDN4 Spearman rho | 0.310 (FDR q=8.68e-09) |
| CLDN4 Pearson rank | #26 of 2617 |
| CLDN4 Pearson rho | 0.319 |
| Actual #1 (Spearman) | PVRL4 (rho=0.578) |

## Top 15 surface-gene partners of TACSTD2 (Spearman)

| rank | gene | spearman_rho | pearson_rho | FDR q |
| --- | --- | --- | --- | --- |
| 1 | PVRL4 | 0.578 | 0.507 | 6.48e-35 |
| 2 | TM4SF1 | 0.438 | 0.382 | 9.53e-18 |
| 3 | TGFA | 0.416 | 0.383 | 6.79e-16 |
| 4 | GPR110 | 0.405 | 0.361 | 5.79e-15 |
| 5 | LYPD3 | 0.396 | 0.375 | 2.44e-14 |
| 6 | GPR87 | 0.393 | 0.389 | 3.86e-14 |
| 7 | ITGB6 | 0.388 | 0.341 | 8.97e-14 |
| 8 | GJB3 | 0.384 | 0.324 | 1.62e-13 |
| 9 | TNFRSF21 | 0.376 | 0.370 | 6.61e-13 |
| 10 | SLC26A9 | 0.373 | 0.381 | 1.01e-12 |
| 11 | LY6G6C | 0.368 | 0.356 | 2.11e-12 |
| 12 | ITGA3 | 0.367 | 0.315 | 2.56e-12 |
| 13 | SERINC2 | 0.366 | 0.315 | 2.56e-12 |
| 14 | PCDH1 | 0.362 | 0.372 | 5.27e-12 |
| 15 | RHBDL2 | 0.359 | 0.361 | 8.38e-12 |
