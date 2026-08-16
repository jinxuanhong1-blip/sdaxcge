# CLDN4 vs TACSTD2 co-expression in TCGA-CESC (honest ranking)

**Question:** Is `CLDN4` the top co-expression partner of `TACSTD2` among surface genes?

**Answer:** NO -- CLDN4 is NOT the single top surface-gene co-expression partner of TACSTD2.

- Cohort: TCGA-CESC, 303 primary-tumour samples (sample-type ['01']).
- Surface-gene universe: 2575 surfaceome genes with defined Spearman ρ (anchor removed; 43 constant/undefined excluded).
- Primary metric: Spearman correlation.

| metric | value |
| --- | --- |
| CLDN4 Spearman rank | **#193 of 2575** (92.544th percentile) |
| CLDN4 Spearman rho | 0.279 (FDR q=5.68e-06) |
| CLDN4 Pearson rank | #204 of 2575 |
| CLDN4 Pearson rho | 0.258 |
| Actual #1 (Spearman) | PVRL4 (rho=0.690) |

## Top 15 surface-gene partners of TACSTD2 (Spearman)

| rank | gene | spearman_rho | pearson_rho | FDR q |
| --- | --- | --- | --- | --- |
| 1 | PVRL4 | 0.690 | 0.764 | 1.01e-40 |
| 2 | DUOXA1 | 0.640 | 0.593 | 2.83e-33 |
| 3 | CLCA4 | 0.625 | 0.581 | 2.70e-31 |
| 4 | TMPRSS11D | 0.618 | 0.576 | 1.92e-30 |
| 5 | SDC1 | 0.604 | 0.623 | 7.48e-29 |
| 6 | LYPD3 | 0.603 | 0.636 | 8.79e-29 |
| 7 | RHBDL2 | 0.602 | 0.615 | 1.14e-28 |
| 8 | GJB5 | 0.592 | 0.626 | 1.43e-27 |
| 9 | MUC21 | 0.588 | 0.550 | 4.40e-27 |
| 10 | DSG3 | 0.585 | 0.566 | 9.30e-27 |
| 11 | GPR115 | 0.582 | 0.603 | 1.62e-26 |
| 12 | GPR87 | 0.569 | 0.613 | 5.17e-25 |
| 13 | DSC3 | 0.568 | 0.564 | 5.21e-25 |
| 14 | THBD | 0.558 | 0.527 | 6.52e-24 |
| 15 | BDKRB1 | 0.557 | 0.543 | 7.75e-24 |

## Histology split (robustness, not the primary ranking)

CESC mixes cervical squamous carcinoma and endocervical adenocarcinoma. The primary ranking above pools all primary tumours (BRCA-analog protocol). The pooled CLDN4 rank is **diluted by adenocarcinoma**, where the TACSTD2–CLDN4 link is weak and not significant.

### squamous (n=250)

- `CLDN4` Spearman rank **#12 of 2560** (99.57th percentile), ρ=0.491 (FDR q=2.92e-14).
- Actual #1: MUC21 (ρ=0.601).

### adeno (n=52)

- `CLDN4` Spearman rank **#317 of 2464** (87.175th percentile), ρ=0.203 (FDR q=3.99e-01).
- Actual #1: PVRL4 (ρ=0.725).

