# CLDN4-high malignant signature: concordant-4 + TCGA, then inverse immune

Additive public analysis. Bulk correlations are not a spatial exclusion test, and this folder does not rewrite the locked concordant-4 T/NK result (ρ = −0.53, n = 65, PR 503) or the locked single-gene OncoSG / GSE10072 / GSE11969 / GSE248378 CLDN4 tables. GSE10072, GSE11969, and GSE248378 are not reopened.

## Signature

Source is the locked malignant pseudobulk OLS on log2(TMM-CPM+1), contrast `q4q1_combined`, n_Q1 = 18, n_Q4 = 16 (PR 503 `de_all.tsv`). Positive logFC is higher in CLDN4-high malignant cells. GSE189357 has no separate per-cohort contrast in that table (Q4 n = 2); stability uses GSE123902, GSE131907, and GSE205335.

Pre-specified keep rules:

- concordant-4: logFC > 0.5, FDR < 0.05, logFC > 0 in at least 2 tested per-cohort contrasts
- held out: CLDN4, TACSTD2, KRT8, KRT18, KRT19
- held out: IFN, MHC-I/APM, and chemokine family genes, plus a fixed T/NK/B/myeloid and HLA list, so the score is not an immune program
- TCGA-LUAD primary tumors with a called ABSOLUTE purity: partial Spearman of the gene vs CLDN4 given KRT8 + KRT18 + KRT19 + purity, BH FDR < 0.05 on this candidate list, partial ρ > 0

| step | n genes |
|---|---:|
| combined DE genes | 15779 |
| FDR < 0.05 and logFC > 0.5 | 272 |
| after hold-outs and 2-cohort direction | 270 |
| present in TCGA tumor+purity matrix | 264 |
| **signature** (partial ρ > 0, FDR < 0.05) | **221** |

TCGA filter n = **490** primary tumors (one `-01` sample per patient, called ABSOLUTE purity). This filter does not use CD8A or ImmuneScore. TCGA immune correlations below are in-sample and are not validation.

Normal-lung / AT2-club markers inside the signature: **1** / 221 (SFTA3). They stay in the primary score. A sensitivity drops them.

Top genes by TCGA partial ρ (full list: `tables/signature_genes.tsv`):

| gene | concordant-4 logFC | concordant-4 FDR | cohorts up | TCGA partial ρ | TCGA FDR | lung marker |
|---|---:|---:|---:|---:|---:|---|
| CLDN3 | +1.88 | 0.0171 | 3/3 | +0.600 | 1.81e-46 |  |
| FAXC | +1.27 | 0.0403 | 3/3 | +0.425 | 1.14e-20 |  |
| CASD1 | +1.32 | 0.0316 | 3/3 | +0.410 | 3.92e-19 |  |
| RSBN1L | +0.65 | 0.0414 | 3/3 | +0.393 | 1.37e-17 |  |
| ARL8A | +0.62 | 0.0475 | 3/3 | +0.382 | 1.21e-16 |  |
| DHX40 | +0.85 | 0.0474 | 3/3 | +0.382 | 1.21e-16 |  |
| PPP1R9A | +1.73 | 0.00833 | 3/3 | +0.380 | 1.48e-16 |  |
| NUDT16L1 | +0.77 | 0.0391 | 3/3 | +0.379 | 1.66e-16 |  |
| DUSP8 | +1.40 | 0.0259 | 3/3 | +0.378 | 1.87e-16 |  |
| PAK6 | +1.33 | 0.0343 | 3/3 | +0.375 | 3.20e-16 |  |
| RPS6KA3 | +0.89 | 0.0358 | 2/3 | +0.363 | 3.67e-15 |  |
| NUBP2 | +0.60 | 0.0129 | 3/3 | +0.355 | 1.43e-14 |  |
| FNIP2 | +1.22 | 0.0306 | 3/3 | +0.354 | 1.77e-14 |  |
| GDE1 | +0.83 | 0.0431 | 3/3 | +0.348 | 5.51e-14 |  |
| ZMYM3 | +0.84 | 0.0204 | 3/3 | +0.346 | 7.63e-14 |  |

## OncoSG QC (locked single gene, not a new result)

Same public z-score matrix and published PURITY as PR 333 (n = 169). ImmuneScore for this QC row is the mean of the deposited z-scores of the A1 8 genes, not ESTIMATE.

| pair | this run ρ | PR 333 ρ | |Δ| |
|---|---:|---:|---:|
| CLDN4 vs CD8A (unadj) | -0.416 | -0.416 | 0.000 |
| CLDN4 vs CD8A (partial) | -0.285 | -0.285 | 0.000 |
| CLDN4 vs ImmuneScore (unadj) | -0.432 | -0.432 | 0.000 |
| CLDN4 vs ImmuneScore (partial) | -0.308 | -0.308 | 0.000 |

QC matches PR 333 within 0.015 on all four locked rhos.

## Inverse immune

Signature score = mean of within-cohort z-scores of signature genes. ImmuneScore = mean of within-cohort z-scores of the same 8 genes. Spearman, two-sided. Partial correlation is Pearson of rank residuals. GEO cohorts have no published purity; the composition control there is a 9-gene stromal mean (FAP, COL1A1/1A2/3A1, DCN, LUM, PDGFRA, TAGLN, ACTA2), which is not ESTIMATE and not ABSOLUTE.

GSE273377 is one FFPE exome-capture stage I LUAD study with a discovery stratum (GPL30173) and a validation stratum (GPL16791). Samples with `passed qc: FALSE` are out. The two strata are not two studies. GSE233774 uses tumor columns only (paracancerous `N*` out). GSE282774 is pN2 LUAD FPKM.

| cohort | role | n | genes used | vs CD8A ρ (p) | vs ImmuneScore ρ (p) | vs CD8A partial | vs ImmuneScore partial |
|---|---|---:|---:|---|---|---|---|
| OncoSG | validation | 169 | 215/221 | -0.578 (1.90e-16) | -0.615 (6.23e-19) | -0.416 (2.00e-08) purity | -0.472 (1.06e-10) purity |
| GSE273377 discovery | validation_stratum | 103 | 221/221 | -0.381 (7.29e-05) | -0.391 (4.34e-05) | -0.350 (0.000308) stromal | -0.362 (0.000187) stromal |
| GSE273377 validation | validation_stratum | 60 | 221/221 | -0.418 (0.000882) | -0.430 (0.000604) | -0.461 (0.000242) stromal | -0.457 (0.000272) stromal |
| GSE233774 tumor | validation | 30 | 217/221 | -0.332 (0.0733) | -0.440 (0.0149) | -0.100 (0.608) stromal | -0.249 (0.192) stromal |
| GSE282774 | validation | 58 | 221/221 | -0.212 (0.111) | -0.296 (0.0239) | -0.202 (0.131) stromal | -0.279 (0.0358) stromal |
| GSE288479 solid | below_n_floor | 8 | 217/221 | +0.167 (0.693) | +0.024 (0.955) | NA | NA |
| TCGA-LUAD (not validation) | discovery_echo | 490 | 221/221 | -0.190 (2.35e-05) | -0.224 (5.18e-07) | -0.082 (0.0702) purity | -0.111 (0.0139) purity |

Partials: OncoSG and the TCGA echo use published purity or ABSOLUTE. GEO rows use the stromal score. A second partial, given CLDN4, asks whether the program still tracks immune genes after the locked single gene.

### OncoSG partials (primary external cohort)

| endpoint | unadj ρ (p) | partial \| PURITY | partial \| CLDN4 | partial \| PURITY+CLDN4 | partial \| KRT8/18/19 |
|---|---|---|---|---|---|
| CD8A | -0.578 (1.90e-16) | -0.416 (2.00e-08) | -0.472 (1.06e-10) | -0.347 (4.25e-06) | -0.618 (7.17e-19) |
| ImmuneScore | -0.615 (6.23e-19) | -0.472 (1.06e-10) | -0.513 (1.22e-12) | -0.403 (6.54e-08) | -0.649 (3.49e-21) |
| CD274 | -0.486 (2.05e-11) | -0.332 (1.09e-05) | -0.380 (3.79e-07) | -0.267 (0.000476) | -0.498 (8.48e-12) |
| IMSIG_T_cells | -0.597 (1.05e-17) | -0.408 (3.90e-08) | -0.458 (4.46e-10) | -0.298 (9.19e-05) | -0.630 (1.02e-19) |

On OncoSG, the signature vs CD8A is inverse at -0.578 (p = 1.90e-16). After published PURITY it is inverse at -0.416 (p = 2.00e-08). After CLDN4 it is inverse at -0.472 (p = 1.06e-10). After PURITY and CLDN4 together it is inverse at -0.347 (p = 4.25e-06).

TCGA is not a validation cohort. In-sample, the signature vs CD8A is inverse at -0.190 (p = 2.35e-05). After ABSOLUTE purity it is -0.082 (p = 0.0702), which does not clear 0.05.

Partial correlation given CLDN4 in each validation stratum (same question outside OncoSG):

| cohort | endpoint | n | ρ \| CLDN4 | p |
|---|---|---:|---:|---:|
| OncoSG | CD8A | 169 | -0.472 | 1.06e-10 |
| OncoSG | ImmuneScore | 169 | -0.513 | 1.22e-12 |
| GSE273377 discovery | CD8A | 103 | -0.435 | 5.01e-06 |
| GSE273377 discovery | ImmuneScore | 103 | -0.464 | 8.92e-07 |
| GSE273377 validation | CD8A | 60 | -0.305 | 0.019 |
| GSE273377 validation | ImmuneScore | 60 | -0.387 | 0.00244 |
| GSE233774 tumor | CD8A | 30 | -0.241 | 0.207 |
| GSE233774 tumor | ImmuneScore | 30 | -0.352 | 0.0614 |
| GSE282774 | CD8A | 58 | -0.197 | 0.142 |
| GSE282774 | ImmuneScore | 58 | -0.283 | 0.0327 |

### Cross-study combination

Independent studies in the combination: OncoSG, GSE273377 (inverse-variance of its two strata), GSE282774, GSE233774. GSE288479 (8 solid regions) is below the n = 20 floor and is not in the combination. TCGA is not in the combination.

| endpoint | model | k | n sum | ρ | p | I² |
|---|---|---:|---:|---:|---:|---:|
| CD8A | DL random | 4 | 420 | -0.409 | 1.82e-05 | 71% |
| CD8A, GSE273377 strata only | fixed | 2 | 163 | -0.394 | 1.73e-07 | 0% |
| ImmuneScore | DL random | 4 | 420 | -0.460 | 7.99e-07 | 71% |
| ImmuneScore, GSE273377 strata only | fixed | 2 | 163 | -0.406 | 6.93e-08 | 0% |

Cross-study unadjusted combination: CD8A DL random ρ = -0.409 (p = 1.82e-05, I² = 71%); ImmuneScore DL random ρ = -0.460 (p = 7.99e-07, I² = 71%). GSE273377 enters as one study (its two strata combined first). I² is about 70% because OncoSG is stronger than the smaller GEO sets. GSE282774 CD8A and GSE233774 CD8A do not clear 0.05 on the 221-gene score; their ImmuneScore rows do. The combination is inverse, and it is not the same magnitude in every cohort.

## Sensitivities

| score | cohort | endpoint | n | ρ | p |
|---|---|---|---:|---:|---:|
| primary | OncoSG | CD8A | 169 | -0.578 | 1.90e-16 |
| primary | OncoSG | ImmuneScore | 169 | -0.615 | 6.23e-19 |
| primary | GSE282774 | CD8A | 58 | -0.212 | 0.111 |
| primary | GSE282774 | ImmuneScore | 58 | -0.296 | 0.0239 |
| primary | GSE233774 tumor | CD8A | 30 | -0.332 | 0.0733 |
| primary | GSE233774 tumor | ImmuneScore | 30 | -0.440 | 0.0149 |
| primary | GSE273377 discovery | CD8A | 103 | -0.381 | 7.29e-05 |
| primary | GSE273377 discovery | ImmuneScore | 103 | -0.391 | 4.34e-05 |
| primary | GSE273377 validation | CD8A | 60 | -0.418 | 0.000882 |
| primary | GSE273377 validation | ImmuneScore | 60 | -0.430 | 0.000604 |
| drop lung markers | OncoSG | CD8A | 169 | -0.579 | 1.54e-16 |
| drop lung markers | OncoSG | ImmuneScore | 169 | -0.616 | 5.09e-19 |
| drop lung markers | GSE282774 | CD8A | 58 | -0.207 | 0.119 |
| drop lung markers | GSE282774 | ImmuneScore | 58 | -0.297 | 0.0237 |
| drop lung markers | GSE233774 tumor | CD8A | 30 | -0.299 | 0.108 |
| drop lung markers | GSE233774 tumor | ImmuneScore | 30 | -0.421 | 0.0206 |
| drop lung markers | GSE273377 discovery | CD8A | 103 | -0.378 | 8.13e-05 |
| drop lung markers | GSE273377 discovery | ImmuneScore | 103 | -0.388 | 5.03e-05 |
| drop lung markers | GSE273377 validation | CD8A | 60 | -0.417 | 0.000905 |
| drop lung markers | GSE273377 validation | ImmuneScore | 60 | -0.429 | 0.000635 |
| top 30 by TCGA partial ρ | OncoSG | CD8A | 169 | -0.478 | 5.01e-11 |
| top 30 by TCGA partial ρ | OncoSG | ImmuneScore | 169 | -0.540 | 3.70e-14 |
| top 30 by TCGA partial ρ | GSE282774 | CD8A | 58 | -0.312 | 0.0169 |
| top 30 by TCGA partial ρ | GSE282774 | ImmuneScore | 58 | -0.356 | 0.00616 |
| top 30 by TCGA partial ρ | GSE233774 tumor | CD8A | 30 | -0.519 | 0.00329 |
| top 30 by TCGA partial ρ | GSE233774 tumor | ImmuneScore | 30 | -0.565 | 0.00113 |
| top 30 by TCGA partial ρ | GSE273377 discovery | CD8A | 103 | -0.253 | 0.00986 |
| top 30 by TCGA partial ρ | GSE273377 discovery | ImmuneScore | 103 | -0.269 | 0.00602 |
| top 30 by TCGA partial ρ | GSE273377 validation | CD8A | 60 | -0.365 | 0.00415 |
| top 30 by TCGA partial ρ | GSE273377 validation | ImmuneScore | 60 | -0.307 | 0.0169 |

OncoSG vs CD8A after dropping the fixed normal-lung markers is inverse at -0.579 (p = 1.54e-16). The top-30 score is inverse at -0.478 (p = 5.01e-11). Primary signature size is 221; no-lung size is 220.

## Null on OncoSG

OncoSG, 400 random gene sets of size 221, seed 20260921. Draws exclude the signature, the keratin covariates, CLDN4, TACSTD2, CD274, the ImmuneScore genes, and the fixed lineage list, so the null is not filled with T/NK markers. Observed signature vs CD8A ρ = -0.578. Fraction of random sets as low or lower: 0.000. Null median ρ = +0.168. Concordant-4 genes that failed the TCGA filter (n = 42 scored): vs CD8A ρ = -0.399 (p = 7.56e-08, n = 169).

## Sets that were screened and not used as validation

| set | why it is not in the combination |
|---|---|
| GSE10072, GSE11969, GSE248378 | locked single-gene bulks; not reopened |
| TCGA-LUAD | used to choose genes |
| GSE288479 solid component | n = 8 patients, below the floor of 20; appendix row only |
| GSE271259 | series matrix: 5 lung tumors and 35 brain metastases. Brain immune context is not this claim. Not scored |
| GSE319666 and cell-line or adjacent-lung series from the 2023–2026 GEO pass | not invasive LUAD tumor bulk |

## What this does not say

A negative bulk correlation is not the CosMx exclusion result and is not an ICI-resistance result. These cohorts are surgical LUAD. The signature was not trained on immune labels, but it was trained on CLDN4, so agreement with CLDN4 is expected. Agreement with CD8A beyond CLDN4 is the extra question, and it is answered in the partial-|CLDN4 row.

## Reproduce

```bash
python3 methods/cldn4_high_malignant_signature/analyze.py
```

