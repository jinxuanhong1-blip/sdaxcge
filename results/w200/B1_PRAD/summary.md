# CLDN4 vs TACSTD2 co-expression in TCGA-PRAD (honest ranking)

**Question:** Is `CLDN4` the top co-expression partner of `TACSTD2` among surface genes?

**Answer:** POINT ESTIMATE YES, STABILITY NO -- CLDN4 is Spearman #1 in this sample, but the lead is not stable (rank-1 in only 29% of bootstraps; Pearson rank #6).

- Cohort: TCGA-PRAD, 497 primary-tumour samples (sample-type ['01']).
- Surface-gene universe: 2614 surfaceome genes with defined correlation (anchor removed; 4 zero-variance genes dropped).
- Primary metric: Spearman correlation.

| metric | value |
| --- | --- |
| CLDN4 Spearman rank | **#1 of 2614** (100.0th percentile) |
| CLDN4 Spearman rho | 0.383 (FDR q=1.25e-15) |
| CLDN4 Pearson rank | #6 of 2614 |
| CLDN4 Pearson rho | 0.377 |
| Actual #1 (Spearman) | CLDN4 (rho=0.383) |

## Top 15 surface-gene partners of TACSTD2 (Spearman)

| rank | gene | spearman_rho | pearson_rho | FDR q |
| --- | --- | --- | --- | --- |
| 1 | CLDN4  <-- focus | 0.383 | 0.377 | 1.25e-15 |
| 2 | SLC44A4 | 0.382 | 0.378 | 1.25e-15 |
| 3 | CLDN3 | 0.373 | 0.345 | 6.79e-15 |
| 4 | PVRL2 | 0.366 | 0.339 | 2.38e-14 |
| 5 | SEMA4B | 0.343 | 0.376 | 1.11e-12 |
| 6 | GPR108 | 0.339 | 0.357 | 2.09e-12 |
| 7 | EFNA1 | 0.338 | 0.339 | 2.37e-12 |
| 8 | RPN1 | 0.333 | 0.303 | 5.00e-12 |
| 9 | SLC37A1 | 0.323 | 0.363 | 3.05e-11 |
| 10 | SERINC2 | 0.318 | 0.331 | 7.12e-11 |
| 11 | LYPD3 | 0.314 | 0.367 | 1.19e-10 |
| 12 | SPINT2 | 0.314 | 0.308 | 1.19e-10 |
| 13 | PCDH1 | 0.312 | 0.298 | 1.50e-10 |
| 14 | FURIN | 0.311 | 0.283 | 1.72e-10 |
| 15 | CLDN7 | 0.310 | 0.213 | 1.77e-10 |

## Honest caveats

- Bulk-tumour mRNA co-expression (mixes tumour/stroma/immune), not protein or single-cell. Xena HiSeqV2 is log2(norm_count+1), not TPM.
- Primary metric is Spearman, matching the BRCA B1 analog. Pearson can disagree: here Pearson rank is #6 (#1 is PVRL4).
- The Spearman lead over #2 SLC44A4 is only 0.0007. Treat '#1' as a near-tie, not a unique winner.
- Bootstrap (n=1000, seed 20260816): CLDN4 Spearman rho 95% CI [0.307, 0.456]; it is rank #1 in **29.0%** of resamples. Most frequent usurpers: SLC44A4 (256/1000), CLDN3 (163/1000), PVRL2 (95/1000).
- Surface universe is the 2018 in-silico surfaceome (legacy HGNC symbols, e.g. PVRL4 = NECTIN4, PVRL2 = NECTIN2). Four zero-variance surface genes were dropped rather than ranked.
- One primary-tumour aliquot per patient in this matrix (n=497 code-01; 52 adjacent-normal and 1 metastatic aliquot exist and were not used).
