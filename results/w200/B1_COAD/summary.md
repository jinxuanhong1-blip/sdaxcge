# CLDN4 vs TACSTD2 co-expression in TCGA-COAD (honest ranking)

**Question:** Is `CLDN4` the top co-expression partner of `TACSTD2` among surface genes?

**Answer:** NO -- CLDN4 is NOT the single top surface-gene co-expression partner of TACSTD2.

- Cohort: TCGA-COAD, 458 primary-tumour samples (GDC `Primary Tumor`, one aliquot per case).
- Expression: GDC STAR - Counts `tpm_unstranded`, log2(TPM+1).
- Surface-gene universe: 2672 surfaceome genes present in the matrix (anchor removed).
- Primary metric: Spearman correlation.

| metric | value |
| --- | --- |
| CLDN4 Spearman rank | **#1876 of 2672** (29.828th percentile) |
| CLDN4 Spearman rho | 0.007 (FDR q=9.26e-01) |
| CLDN4 Pearson rank | #1716 of 2672 |
| CLDN4 Pearson rho | 0.029 |
| Actual #1 (Spearman) | LYPD3 (rho=0.389) |

## Top 15 surface-gene partners of TACSTD2 (Spearman)

| rank | gene | spearman_rho | pearson_rho | FDR q |
| --- | --- | --- | --- | --- |
| 1 | LYPD3 | 0.389 | 0.412 | 1.49e-14 |
| 2 | MSLN | 0.379 | 0.374 | 5.20e-14 |
| 3 | PTPRU | 0.358 | 0.356 | 1.77e-12 |
| 4 | PSCA | 0.357 | 0.319 | 1.77e-12 |
| 5 | ALPP | 0.357 | 0.311 | 1.77e-12 |
| 6 | GJB4 | 0.355 | 0.351 | 2.04e-12 |
| 7 | GJB3 | 0.347 | 0.318 | 7.51e-12 |
| 8 | ADAM8 | 0.327 | 0.338 | 2.50e-10 |
| 9 | GJB5 | 0.321 | 0.318 | 6.13e-10 |
| 10 | CHRNB1 | 0.313 | 0.307 | 1.87e-09 |
| 11 | TLR5 | 0.305 | 0.318 | 6.48e-09 |
| 12 | SEMA4A | 0.304 | 0.304 | 7.14e-09 |
| 13 | CALHM2 | 0.293 | 0.302 | 3.48e-08 |
| 14 | VTCN1 | 0.289 | 0.277 | 5.70e-08 |
| 15 | VSIG2 | 0.288 | 0.283 | 6.26e-08 |

## Expression context (why the rank is this low)

`CLDN4` is constitutively high in COAD primary tumours (median log2(TPM+1) = 9.15; detected in 100% of samples). `TACSTD2` is moderate and more variable (median 4.53, range 0.29–10.54). The two vectors do not co-vary: Spearman ρ = 0.007, FDR q = 9.26e-01. This is a null co-expression result, not a near-miss of #1.

## Independent matrix (UCSC Xena HiSeqV2)

Same question on the legacy Xena `TCGA.COAD.sampleMap/HiSeqV2` matrix (log2(norm_count+1), sample-type `01`, n=286, 2510 ranked surface genes):

- CLDN4 Spearman rank **#2136 of 2510**, ρ = -0.058 (FDR q=4.73e-01).
- Actual #1: PVRL4 (ρ = 0.448).

The GDC STAR-Counts and Xena HiSeqV2 matrices agree: CLDN4 is not a TACSTD2 surface-gene partner in COAD.

## Caveats

- Bulk-tumour mRNA co-expression (mixes tumour / stroma / immune); not protein and not single-cell.
- Surfaceome gene symbols are the 2018 UniProt names (e.g. `PVRL4` = `NECTIN4`).
- 13 surfaceome genes with zero variance in this cohort were dropped (undefined correlation), not ranked.
- One Primary Tumor aliquot per case; technical-replicate files were not averaged.
- This COAD null does not contradict the BRCA B1 analog (CLDN4 #4, ρ≈0.35) or the PAAD analog (ρ≈0.71); those are different tissues.
- TCGA-COAD is untreated archival RNA-seq; it is not an ICI or ADC outcome cohort.
