# CLDN4 vs TACSTD2 co-expression in TCGA-BRCA (honest ranking)

**Question:** Is `CLDN4` the top co-expression partner of `TACSTD2` among surface genes?

**Answer:** NO — CLDN4 ranks **#4 of 2618** surface genes (Spearman ρ = 0.348). #1 is EFNA1. Genes ahead of CLDN4: EFNA1, PVRL4, EFNA4.

- Main cohort: 1097 primary tumours, one sample per patient (all-vial n = 1097 gives the same rank #4).
- Pair-wise TACSTD2–CLDN4 Spearman ρ = 0.348 [0.29, 0.40].
- Transcriptome-wide rank: #30 of 17659.
- Claudin-family rank: #1 (top claudin = CLDN4).
- Purity-partial ρ = 0.343 (unadjusted on the same n = 1047: 0.346).

See `README.md` in this folder for the full bilingual write-up.

## Top 15 surface-gene partners (one-per-patient, Spearman)

| rank | gene | spearman_rho | pearson_rho | FDR q |
| --- | --- | --- | --- | --- |
| 1 | EFNA1 | 0.429 | 0.394 | 6.96e-47 |
| 2 | PVRL4 | 0.351 | 0.405 | 4.86e-30 |
| 3 | EFNA4 | 0.350 | 0.348 | 4.86e-30 |
| 4 | CLDN4  <-- focus | 0.348 | 0.386 | 1.11e-29 |
| 5 | TM4SF1 | 0.341 | 0.325 | 1.59e-28 |
| 6 | MPZL2 | 0.340 | 0.347 | 1.69e-28 |
| 7 | EPHA2 | 0.340 | 0.274 | 1.69e-28 |
| 8 | ITGB4 | 0.328 | 0.341 | 1.85e-26 |
| 9 | CD151 | 0.325 | 0.317 | 5.15e-26 |
| 10 | SLC5A1 | 0.322 | 0.303 | 2.16e-25 |
| 11 | MUC16 | 0.321 | 0.305 | 2.93e-25 |
| 12 | GPR87 | 0.320 | 0.300 | 3.65e-25 |
| 13 | FOLR1 | 0.316 | 0.292 | 1.33e-24 |
| 14 | RHBDL2 | 0.313 | 0.307 | 4.08e-24 |
| 15 | ABCA4 | 0.308 | 0.298 | 2.35e-23 |
