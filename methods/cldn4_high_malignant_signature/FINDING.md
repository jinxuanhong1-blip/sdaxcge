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


## Weaker GEO cohorts: size, ssGSEA, ESTIMATE, histology

The primary 221-gene z-mean result is unchanged. GSE282774 CD8A ρ = −0.212 (p = 0.11) and GSE233774 CD8A ρ = −0.332 (p = 0.073) on that score. This section is the requested grid. A star marks an inverse association with p < 0.05. Cells that do not clear 0.05 are left in the tables.

ssGSEA of the full 221-gene set is inverse with p < 0.05 for CD8A and ImmuneScore in both cohorts: GSE282774 ssGSEA n = 221 vs CD8A -0.486 (0.000109)*; GSE282774 ssGSEA n = 221 vs ImmuneScore -0.509 (4.41e-05)*; GSE233774 ssGSEA n = 221 vs CD8A -0.420 (0.0209)*; GSE233774 ssGSEA n = 221 vs ImmuneScore -0.578 (0.000818)*. Shortening the z-mean score also clears unadjusted CD8A at sizes 20–100 in GSE282774 and at sizes 10–100 in GSE233774. Size 221 z-mean CD8A still does not. After ESTIMATE StromalScore, GSE282774 ssGSEA vs CD8A stays inverse (-0.447 (0.000494)* at size 221). GSE233774 ssGSEA vs CD8A does not clear after StromalScore (best -0.343 (0.0685) at size 50). Restricting GSE233774 to IAC (n = 22) does not clear CD8A at any size or score (best unadjusted ssGSEA size 50 -0.397 (0.0674)). ImmuneScore on that IAC subset still clears (ssGSEA size 221 -0.554 (0.00748)*). GSE282774 has no histologic subtype to filter.

GSE233774 pathologic diagnosis, from `GSE233774_Pathological_and_radiological_information.xlsx`: IAC 22, MIA 3, AIS 3, AAH 2. All 30 tumor columns have a diagnosis. This is AAH / AIS / MIA / IAC, not a WHO lepidic/acinar/solid label. GSE282774's series matrix and expression file contain only `disease: pN2 LUAD`. There is no subtype column to filter.

ESTIMATE scores use package 1.0.13 gene sets (Yoshihara et al., Nat Commun 2013): common-gene filter, then ssGSEA of the 141-gene stromal set and the 141-gene immune set. ESTIMATEScore is their sum. Tumor purity from `cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)` is Affymetrix-only in that package and is not used here. Partial correlation on ESTIMATEScore removes immune signal as well as stroma, because the score contains ImmuneScore. StromalScore is the composition control that is not the endpoint.

ESTIMATE coverage: GSE282774 all pN2 LUAD: n = 58, common genes 9885, stromal genes 136/141, immune genes 140/141. GSE233774 all tumors: n = 30, common genes 8327, stromal genes 135/141, immune genes 139/141. GSE233774 IAC+MIA: n = 25, common genes 8327, stromal genes 135/141, immune genes 139/141. GSE233774 IAC only: n = 22, common genes 8327, stromal genes 135/141, immune genes 139/141.

### Unadjusted, all samples

| cohort | subset | method | endpoint | n = 10 | n = 20 | n = 30 | n = 50 | n = 100 | n = 221 |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| GSE282774 | all pN2 LUAD | z-mean | CD8A | -0.249 (0.0592) | -0.276 (0.0357)* | -0.312 (0.0169)* | -0.291 (0.0267)* | -0.260 (0.0488)* | -0.212 (0.111) |
| GSE282774 | all pN2 LUAD | z-mean | ImmuneScore | -0.326 (0.0126)* | -0.348 (0.00746)* | -0.356 (0.00616)* | -0.346 (0.00783)* | -0.327 (0.0122)* | -0.296 (0.0239)* |
| GSE282774 | all pN2 LUAD | ssGSEA | CD8A | -0.577 (2.11e-06)* | -0.592 (1.00e-06)* | -0.567 (3.41e-06)* | -0.602 (5.76e-07)* | -0.567 (3.49e-06)* | -0.486 (0.000109)* |
| GSE282774 | all pN2 LUAD | ssGSEA | ImmuneScore | -0.629 (1.24e-07)* | -0.657 (2.14e-08)* | -0.586 (1.34e-06)* | -0.622 (1.87e-07)* | -0.602 (5.65e-07)* | -0.509 (4.41e-05)* |
| GSE233774 | all tumors | z-mean | CD8A | -0.464 (0.00983)* | -0.403 (0.0273)* | -0.519 (0.00329)* | -0.475 (0.00799)* | -0.412 (0.0236)* | -0.332 (0.0733) |
| GSE233774 | all tumors | z-mean | ImmuneScore | -0.591 (0.00059)* | -0.524 (0.00293)* | -0.565 (0.00113)* | -0.519 (0.00329)* | -0.528 (0.00269)* | -0.440 (0.0149)* |
| GSE233774 | all tumors | ssGSEA | CD8A | -0.438 (0.0156)* | -0.423 (0.0197)* | -0.526 (0.00285)* | -0.523 (0.00299)* | -0.496 (0.00528)* | -0.420 (0.0209)* |
| GSE233774 | all tumors | ssGSEA | ImmuneScore | -0.576 (0.000876)* | -0.591 (0.00059)* | -0.617 (0.000279)* | -0.586 (0.000672)* | -0.637 (0.000154)* | -0.578 (0.000818)* |

### Partial correlation given ESTIMATE StromalScore, all samples

| cohort | subset | method | endpoint | n = 10 | n = 20 | n = 30 | n = 50 | n = 100 | n = 221 |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| GSE282774 | all pN2 LUAD | z-mean | CD8A | -0.152 (0.259) | -0.185 (0.169) | -0.230 (0.0846) | -0.203 (0.13) | -0.172 (0.2) | -0.127 (0.345) |
| GSE282774 | all pN2 LUAD | z-mean | ImmuneScore | -0.201 (0.133) | -0.230 (0.0853) | -0.245 (0.0661) | -0.229 (0.0861) | -0.215 (0.109) | -0.192 (0.153) |
| GSE282774 | all pN2 LUAD | ssGSEA | CD8A | -0.534 (1.89e-05)* | -0.550 (9.37e-06)* | -0.524 (2.88e-05)* | -0.563 (5.08e-06)* | -0.527 (2.49e-05)* | -0.447 (0.000494)* |
| GSE282774 | all pN2 LUAD | ssGSEA | ImmuneScore | -0.560 (5.86e-06)* | -0.596 (9.81e-07)* | -0.523 (2.98e-05)* | -0.566 (4.45e-06)* | -0.551 (8.94e-06)* | -0.460 (0.00032)* |
| GSE233774 | all tumors | z-mean | CD8A | -0.249 (0.194) | -0.128 (0.509) | -0.332 (0.0787) | -0.263 (0.168) | -0.173 (0.37) | -0.135 (0.485) |
| GSE233774 | all tumors | z-mean | ImmuneScore | -0.430 (0.0198)* | -0.314 (0.0976) | -0.390 (0.0363)* | -0.317 (0.0934) | -0.336 (0.0746) | -0.268 (0.159) |
| GSE233774 | all tumors | ssGSEA | CD8A | -0.204 (0.288) | -0.162 (0.402) | -0.343 (0.0686) | -0.343 (0.0685) | -0.297 (0.118) | -0.229 (0.233) |
| GSE233774 | all tumors | ssGSEA | ImmuneScore | -0.405 (0.0291)* | -0.425 (0.0216)* | -0.472 (0.00981)* | -0.425 (0.0215)* | -0.502 (0.00555)* | -0.436 (0.0182)* |

### Partial correlation given ESTIMATEScore, all samples

This adjustment includes the immune ssGSEA, so a lost CD8 association here is not evidence against the unadjusted result.

| cohort | subset | method | endpoint | n = 10 | n = 20 | n = 30 | n = 50 | n = 100 | n = 221 |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| GSE282774 | all pN2 LUAD | z-mean | CD8A | +0.096 (0.477) | +0.045 (0.738) | -0.051 (0.707) | -0.026 (0.85) | -0.004 (0.974) | +0.015 (0.91) |
| GSE282774 | all pN2 LUAD | z-mean | ImmuneScore | +0.073 (0.588) | +0.023 (0.865) | -0.047 (0.727) | -0.039 (0.775) | -0.038 (0.782) | -0.049 (0.717) |
| GSE282774 | all pN2 LUAD | ssGSEA | CD8A | -0.318 (0.0161)* | -0.362 (0.00563)* | -0.405 (0.0018)* | -0.458 (0.000336)* | -0.435 (0.000711)* | -0.399 (0.00213)* |
| GSE282774 | all pN2 LUAD | ssGSEA | ImmuneScore | -0.321 (0.015)* | -0.402 (0.00195)* | -0.401 (0.00201)* | -0.464 (0.000276)* | -0.470 (0.000226)* | -0.430 (0.000845)* |
| GSE233774 | all tumors | z-mean | CD8A | -0.179 (0.352) | -0.071 (0.713) | -0.266 (0.163) | -0.219 (0.253) | -0.125 (0.52) | -0.084 (0.663) |
| GSE233774 | all tumors | z-mean | ImmuneScore | -0.403 (0.0303)* | -0.294 (0.122) | -0.363 (0.0529) | -0.307 (0.105) | -0.322 (0.089) | -0.249 (0.193) |
| GSE233774 | all tumors | ssGSEA | CD8A | -0.118 (0.542) | -0.093 (0.632) | -0.274 (0.15) | -0.282 (0.138) | -0.228 (0.234) | -0.170 (0.378) |
| GSE233774 | all tumors | ssGSEA | ImmuneScore | -0.371 (0.0477)* | -0.395 (0.0342)* | -0.443 (0.0161)* | -0.400 (0.0315)* | -0.474 (0.00944)* | -0.413 (0.0259)* |

### GSE233774 histology subsets, unadjusted

| cohort | subset | method | endpoint | n = 10 | n = 20 | n = 30 | n = 50 | n = 100 | n = 221 |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| GSE233774 | IAC+MIA | z-mean | CD8A | -0.400 (0.0476)* | -0.359 (0.0778) | -0.457 (0.0217)* | -0.442 (0.0268)* | -0.361 (0.0764) | -0.248 (0.231) |
| GSE233774 | IAC+MIA | z-mean | ImmuneScore | -0.578 (0.00249)* | -0.533 (0.00607)* | -0.568 (0.00303)* | -0.512 (0.00884)* | -0.493 (0.0123)* | -0.382 (0.0593) |
| GSE233774 | IAC+MIA | ssGSEA | CD8A | -0.392 (0.0524) | -0.375 (0.0644) | -0.484 (0.0143)* | -0.486 (0.0137)* | -0.461 (0.0204)* | -0.377 (0.0633) |
| GSE233774 | IAC+MIA | ssGSEA | ImmuneScore | -0.610 (0.00121)* | -0.597 (0.00163)* | -0.644 (0.000515)* | -0.609 (0.00123)* | -0.635 (0.000643)* | -0.577 (0.00253)* |
| GSE233774 | IAC only | z-mean | CD8A | -0.311 (0.159) | -0.225 (0.313) | -0.375 (0.0851) | -0.335 (0.128) | -0.278 (0.21) | -0.202 (0.368) |
| GSE233774 | IAC only | z-mean | ImmuneScore | -0.542 (0.00925)* | -0.446 (0.0377)* | -0.496 (0.0188)* | -0.424 (0.0492)* | -0.442 (0.0394)* | -0.407 (0.06) |
| GSE233774 | IAC only | ssGSEA | CD8A | -0.302 (0.172) | -0.286 (0.196) | -0.396 (0.0682) | -0.397 (0.0674) | -0.351 (0.11) | -0.310 (0.16) |
| GSE233774 | IAC only | ssGSEA | ImmuneScore | -0.571 (0.00553)* | -0.549 (0.00809)* | -0.592 (0.00368)* | -0.554 (0.00748)* | -0.560 (0.00677)* | -0.554 (0.00748)* |

### Where the inverse association clears 0.05

The grid has 288 cells. 132 are inverse with p < 0.05. That count is not a family-wise error rate. The primary test is still the 221-gene z-mean.

**GSE282774, all pN2 LUAD.**
- CD8A: z-mean clears at size 20, 30, 50, 100; z-mean | ESTIMATE StromalScore does not clear (best -0.230 (0.0846) at size 30); z-mean | ESTIMATEScore does not clear (best -0.051 (0.707) at size 30); ssGSEA clears at size 10, 20, 30, 50, 100, 221; ssGSEA | ESTIMATE StromalScore clears at size 10, 20, 30, 50, 100, 221; ssGSEA | ESTIMATEScore clears at size 10, 20, 30, 50, 100, 221.
- ImmuneScore: z-mean clears at size 10, 20, 30, 50, 100, 221; z-mean | ESTIMATE StromalScore does not clear (best -0.245 (0.0661) at size 30); z-mean | ESTIMATEScore does not clear (best -0.049 (0.717) at size 221); ssGSEA clears at size 10, 20, 30, 50, 100, 221; ssGSEA | ESTIMATE StromalScore clears at size 10, 20, 30, 50, 100, 221; ssGSEA | ESTIMATEScore clears at size 10, 20, 30, 50, 100, 221.

**GSE233774, all tumors.**
- CD8A: z-mean clears at size 10, 20, 30, 50, 100; z-mean | ESTIMATE StromalScore does not clear (best -0.332 (0.0787) at size 30); z-mean | ESTIMATEScore does not clear (best -0.266 (0.163) at size 30); ssGSEA clears at size 10, 20, 30, 50, 100, 221; ssGSEA | ESTIMATE StromalScore does not clear (best -0.343 (0.0685) at size 50); ssGSEA | ESTIMATEScore does not clear (best -0.282 (0.138) at size 50).
- ImmuneScore: z-mean clears at size 10, 20, 30, 50, 100, 221; z-mean | ESTIMATE StromalScore clears at size 10, 30; z-mean | ESTIMATEScore clears at size 10; ssGSEA clears at size 10, 20, 30, 50, 100, 221; ssGSEA | ESTIMATE StromalScore clears at size 10, 20, 30, 50, 100, 221; ssGSEA | ESTIMATEScore clears at size 10, 20, 30, 50, 100, 221.

**GSE233774, IAC only.**
- CD8A: z-mean does not clear (best -0.375 (0.0851) at size 30); z-mean | ESTIMATE StromalScore does not clear (best -0.250 (0.275) at size 30); z-mean | ESTIMATEScore does not clear (best -0.182 (0.431) at size 30); ssGSEA does not clear (best -0.397 (0.0674) at size 50); ssGSEA | ESTIMATE StromalScore does not clear (best -0.284 (0.213) at size 50); ssGSEA | ESTIMATEScore does not clear (best -0.224 (0.328) at size 50).
- ImmuneScore: z-mean clears at size 10, 20, 30, 50, 100; z-mean | ESTIMATE StromalScore clears at size 10; z-mean | ESTIMATEScore does not clear (best -0.398 (0.0737) at size 10); ssGSEA clears at size 10, 20, 30, 50, 100, 221; ssGSEA | ESTIMATE StromalScore clears at size 10, 20, 30, 50, 100, 221; ssGSEA | ESTIMATEScore clears at size 30.

**GSE233774, IAC+MIA.**
- CD8A: z-mean clears at size 10, 30, 50; z-mean | ESTIMATE StromalScore does not clear (best -0.267 (0.207) at size 30); z-mean | ESTIMATEScore does not clear (best -0.203 (0.343) at size 50); ssGSEA clears at size 30, 50, 100; ssGSEA | ESTIMATE StromalScore does not clear (best -0.313 (0.137) at size 50); ssGSEA | ESTIMATEScore does not clear (best -0.249 (0.241) at size 50).
- ImmuneScore: z-mean clears at size 10, 20, 30, 50, 100; z-mean | ESTIMATE StromalScore clears at size 10, 30; z-mean | ESTIMATEScore does not clear (best -0.384 (0.0636) at size 10); ssGSEA clears at size 10, 20, 30, 50, 100, 221; ssGSEA | ESTIMATE StromalScore clears at size 10, 20, 30, 50, 100, 221; ssGSEA | ESTIMATEScore clears at size 10, 30, 50, 100.

Full grid: `tables/weak_geo_grid.tsv`. Diagnosis map: `tables/gse233774_diagnosis.tsv`.

