# CLDN4 vs TACSTD2 co-expression in TCGA-UCEC (honest ranking)

**Question:** Is `CLDN4` the top co-expression partner of `TACSTD2` among surface genes?

**Answer:** NO -- CLDN4 is NOT the single top surface-gene co-expression partner of TACSTD2.

- Cohort: TCGA-UCEC, 549 primary-tumour samples (sample-type ['01']; 545 unique patients).
- Surface-gene universe: 2751 surfaceome genes with defined Spearman (anchor removed; 7 zero-variance genes dropped).
- Primary metric: Spearman correlation.
- Expression: Xena GDC STAR FPKM-UQ, log2(fpkm-uq+1).

| metric | value |
| --- | --- |
| CLDN4 Spearman rank | **#22 of 2751** (99.237th percentile) |
| CLDN4 Spearman rho | 0.314 (p=5.14e-14, FDR q=6.33e-12) |
| CLDN4 Pearson rank | #4 of 2751 |
| CLDN4 Pearson rho | 0.442 |
| Actual #1 (Spearman) | NECTIN4 (rho=0.465) |

## Top 25 surface-gene partners of TACSTD2 (Spearman)

| rank | gene | spearman_rho | pearson_rho | FDR q |
| --- | --- | --- | --- | --- |
| 1 | NECTIN4 | 0.465 | 0.574 | 2.05e-27 |
| 2 | SLC52A3 | 0.412 | 0.448 | 9.98e-21 |
| 3 | TGFA | 0.410 | 0.377 | 1.19e-20 |
| 4 | PLB1 | 0.361 | 0.302 | 1.64e-15 |
| 5 | CLCA4 | 0.359 | 0.261 | 2.25e-15 |
| 6 | SLC34A2 | 0.358 | 0.469 | 2.38e-15 |
| 7 | BACE2 | 0.356 | 0.412 | 2.83e-15 |
| 8 | CDCP1 | 0.347 | 0.372 | 1.75e-14 |
| 9 | EMP2 | 0.343 | 0.355 | 3.76e-14 |
| 10 | PODXL | 0.340 | 0.360 | 6.58e-14 |
| 11 | ITGB4 | 0.337 | 0.384 | 1.18e-13 |
| 12 | PLAUR | 0.327 | 0.316 | 8.86e-13 |
| 13 | MPZL2 | 0.326 | 0.411 | 9.20e-13 |
| 14 | ECE1 | 0.325 | 0.302 | 1.06e-12 |
| 15 | CD59 | 0.323 | 0.326 | 1.41e-12 |
| 16 | PSCA | 0.323 | 0.290 | 1.41e-12 |
| 17 | TPBG | 0.323 | 0.280 | 1.41e-12 |
| 18 | RAET1L | 0.322 | 0.240 | 1.62e-12 |
| 19 | MUC1 | 0.320 | 0.408 | 2.17e-12 |
| 20 | UNC93A | 0.317 | 0.222 | 4.21e-12 |
| 21 | SLC28A3 | 0.315 | 0.301 | 5.13e-12 |
| 22 | CLDN4  <-- focus | 0.314 | 0.442 | 6.33e-12 |
| 23 | GJB3 | 0.314 | 0.312 | 6.33e-12 |
| 24 | SCNN1G | 0.313 | 0.322 | 6.96e-12 |
| 25 | RHBDL2 | 0.312 | 0.338 | 7.78e-12 |

## Sensitivities (same Spearman, not cherry-picked)

- **symbol_only_no_ensembl_fallback:** CLDN4 rank #21 of 2678 (rho=0.314); #1 is SLC52A3 (rho=0.412).
- **one_sample_per_patient:** CLDN4 rank #25 of 2751 (rho=0.312); #1 is NECTIN4 (rho=0.462).

## Surfaceome mapping

- Table S3 rows: 2882; unique UniProt symbols: 2799.
- Matched by 2018 UniProt symbol: 2769.
- Recovered by Ensembl id → GENCODE v36 symbol: 73 (e.g. TMEM30C→TMEM30CP, KIAA0922→TMEM131L, PCNXL2→PCNX2, HIDE1→C19orf38, KIAA1324L→ELAPOR2, UPK3BL→UPK3BL1, C1orf233→FNDC10, PPAP2A→PLPP1).

## Caveats

- Bulk-tumour mRNA co-expression (mixes tumour/stroma/immune), not protein or single-cell.
- UCEC mixes endometrioid and serous histologies; this ranking is not histology-stratified.
- Surfaceome symbols are 2018 UniProt names; renamed genes (PVRL4→NECTIN4) enter via Ensembl ids.
- Rank is among surfaceome genes only, not the full transcriptome.
- Pearson and Spearman disagree on exact rank (CLDN4 is higher by Pearson); Spearman is primary.
