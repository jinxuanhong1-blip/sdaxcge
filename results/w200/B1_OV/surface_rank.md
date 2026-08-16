# CLDN4 vs TACSTD2 co-expression in TCGA-OV (honest surface-gene ranking)

**Question:** Is `CLDN4` the top co-expression partner of `TACSTD2` among surface genes?

**Answer:** NO — CLDN4 is NOT the single top surface-gene co-expression partner of TACSTD2.

- Cohort: TCGA-OV, **n = 421** primary-tumour samples (sample-type `01`, one per patient).
- Surface-gene universe: **2,677** Bausch-Fluck 2018 surfaceome genes present in the STAR-TPM matrix with non-zero variance (anchor removed). Table S3 had 2,799 unique UniProt gene names; 113 were absent from GENCODE v36 symbols in this matrix; 8 had zero variance (rho undefined) and were dropped: CLDN22, GP1BB, ICAM4, LRRC24, OR14I1, OR4N4, OR5M9, SERINC4.
- Primary metric: Spearman correlation of log2(TPM+1).
- Expression: Xena GDC hub `TCGA-OV.star_tpm` (GENCODE v36).

| metric | value |
| --- | --- |
| n | **421** |
| CLDN4 Spearman ρ | **0.320** (p = 1.64e-11; BH-FDR q = 1.51e-09) |
| CLDN4 Spearman rank | **#27 of 2,677** (99.029th percentile) |
| CLDN4 Pearson r | 0.407 (rank #10 of 2,677) |
| Actual #1 (Spearman) | TGFA (ρ = 0.470) |

CLDN4 is a **high** partner (top ~1%) but not #1. The same pair was #4 of 2,618 in the B1_BRCA analog (ρ = 0.348). The OV ρ matches the pairwise TACSTD2–CLDN4 Spearman already reported in `cohort_summary.json` (0.32).

## Top 15 surface-gene partners of TACSTD2 (Spearman)

| rank | gene | spearman_rho | pearson_r | FDR q |
| --- | --- | --- | --- | --- |
| 1 | TGFA | 0.470 | 0.445 | 3.95e-21 |
| 2 | UNC93B1 | 0.450 | 0.438 | 3.04e-19 |
| 3 | ITGB8 | 0.431 | 0.431 | 1.44e-17 |
| 4 | EVA1C | 0.408 | 0.393 | 1.58e-15 |
| 5 | CXCL16 | 0.403 | 0.465 | 3.95e-15 |
| 6 | SLC37A1 | 0.397 | 0.415 | 9.78e-15 |
| 7 | PROM2 | 0.375 | 0.408 | 6.22e-13 |
| 8 | CDCP1 | 0.357 | 0.407 | 1.34e-11 |
| 9 | GPR132 | 0.356 | 0.342 | 1.48e-11 |
| 10 | NIPAL2 | 0.355 | 0.337 | 1.55e-11 |
| 11 | GJB3 | 0.349 | 0.388 | 4.08e-11 |
| 12 | SERINC2 | 0.348 | 0.374 | 4.60e-11 |
| 13 | PIEZO1 | 0.342 | 0.321 | 1.11e-10 |
| 14 | HRH1 | 0.340 | 0.306 | 1.48e-10 |
| 15 | GPR157 | 0.339 | 0.348 | 1.70e-10 |
| 27 | CLDN4  ← focus | 0.320 | 0.407 | 1.51e-09 |

Full ranking: `coexpression_TACSTD2_surfaceome.csv`. Top-200 window: `top200.csv`. Headline row: `surface_rank_headline.csv`.
