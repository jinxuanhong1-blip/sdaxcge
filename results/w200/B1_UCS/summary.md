# CLDN4 vs TACSTD2 co-expression in TCGA-UCS (honest ranking)

**Question:** Is `CLDN4` the top co-expression partner of `TACSTD2` among surface genes in TCGA-UCS?

**Answer:** NO -- CLDN4 is NOT the single top surface-gene co-expression partner of TACSTD2.

- Cohort: TCGA-UCS, 57 primary-tumour samples (sample-type 01; one per patient). All cases are mixed Müllerian / carcinosarcoma.
- Surface-gene universe: 2671 surfaceome genes present in the STAR TPM matrix (anchor removed).
- Primary metric: Spearman correlation on log2(TPM+1).
- Small-n note: n=57. CLDN4 surface-rank bootstrap median 8 (95% 1–28); fraction of resamples with rank 1 = 0.09.

| metric | value |
| --- | --- |
| CLDN4 Spearman rank | **#3 of 2671** (99.925th percentile) |
| CLDN4 Spearman rho | 0.835 (p=7.15e-16, FDR q=6.37e-13) |
| CLDN4 Spearman rho 95% CI | 0.705 – 0.913 |
| CLDN4 Pearson rank | #1 of 2671 |
| CLDN4 Pearson rho | 0.863 |
| Actual #1 (Spearman) | CRB3 (rho=0.851) |
| Partial rho given ABSOLUTE purity | 0.819 (n=56) |

## Top 15 surface-gene partners of TACSTD2 (Spearman)

| rank | gene | spearman_rho | pearson_rho | FDR q |
| --- | --- | --- | --- | --- |
| 1 | CRB3 | 0.851 | 0.860 | 1.32e-13 |
| 2 | EPHA1 | 0.835 | 0.843 | 6.37e-13 |
| 3 | CLDN4  <-- focus | 0.835 | 0.863 | 6.37e-13 |
| 4 | SCNN1A | 0.833 | 0.832 | 6.48e-13 |
| 5 | SLC44A4 | 0.831 | 0.832 | 7.00e-13 |
| 6 | MFSD6L | 0.829 | 0.805 | 7.25e-13 |
| 7 | CLDN7 | 0.827 | 0.825 | 8.36e-13 |
| 8 | F11R | 0.826 | 0.810 | 8.66e-13 |
| 9 | VTCN1 | 0.816 | 0.826 | 3.05e-12 |
| 10 | PRSS8 | 0.815 | 0.856 | 3.39e-12 |
| 11 | CLDN3 | 0.814 | 0.814 | 3.54e-12 |
| 12 | PROM2 | 0.808 | 0.812 | 6.86e-12 |
| 13 | MPZL3 | 0.803 | 0.795 | 1.17e-11 |
| 14 | CDH1 | 0.799 | 0.802 | 1.87e-11 |
| 15 | GPRC5A | 0.791 | 0.778 | 4.13e-11 |
