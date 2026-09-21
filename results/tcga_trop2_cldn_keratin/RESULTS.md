# TCGA: TACSTD2 tracks CLDN4 and CLDN7 after keratin adjustment

Numbers in this file are written by `scripts/tcga_trop2_cldn_keratin/run_analysis.py`. They are not transcribed by hand.

The immune-specification sweep (gene sets, KRT5/6 versus KRT8/18/19, histology, ESTIMATE and ABSOLUTE purity, CLDN4 quantiles, Spearman versus Pearson residualization, and NHEJ/STING/IFN modules) is in `SWEEP.md`. The locked CLDN4-versus-KRT8 surface-gene screen is not rerun there.

## 中文摘要

在 TCGA 上皮癌原发灶（Xena GDC STAR log2(TPM+1)，样本类型 01）里，校正 KRT8+KRT18+KRT19 之后，TACSTD2–CLDN4 偏相关正向 19/21，随机效应 meta ρ=0.315 (95% CI 0.251 to 0.376), p=3.82e-20, I²=89.1%, cohorts=21。TACSTD2–CLDN7 偏相关正向 19/21，meta ρ=0.242 (95% CI 0.178 to 0.304), p=4.35e-13, I²=88.1%, cohorts=21。肺加角蛋白漏斗（LUAD、LUSC、BRCA、CESC、KIRC、STAD、BLCA、PAAD）里这两对都是 8/8 与 8/8。不正向的上皮癌队列在英文 Coexpression 一节按表列出（结直肠 COAD、READ）。

免疫排斥（同一角蛋白校正，预指定队列为 LUAD、LUSC 加上已锁定的角蛋白漏斗：LUAD、BRCA、CESC、KIRC、STAD、BLCA、PAAD）看的是与 CD8 分数的偏相关，负值表示角蛋白校正后仍更少 CD8。TACSTD2–CD8 负向 6/8，meta ρ=-0.069 (95% CI -0.135 to -0.002), p=0.0427, I²=76.1%, cohorts=8。CLDN4–CD8 负向 7/8，meta ρ=-0.081 (95% CI -0.138 to -0.023), p=0.0062, I²=68.1%, cohorts=8。CLDN7–CD8 负向 7/8，meta ρ=-0.109 (95% CI -0.173 to -0.045), p=0.0009, I²=74.9%, cohorts=8。

CLDN4–KRT8（两边只校正 KRT18/KRT19）在这套原发灶 STAR 偏相关里，漏斗七队列 7/7 为正，BRCA 是 0.202。这不是已锁定表面分子筛选里的 BRCA −0.071，也不重估该筛选的 q 值；锁定数字保持原样。

## Question

Does Trop2 (`TACSTD2`) track the claudin barrier genes `CLDN4` and `CLDN7` in TCGA epithelial carcinomas, and is that program still associated with lower CD8 / cytotoxic signal after a simple-epithelial keratin adjustment?

## Data

| item | choice |
| --- | --- |
| expression | UCSC Xena GDC hub `TCGA-<COHORT>.star_tpm.tsv.gz`, log2(TPM+1), GENCODE v36 |
| samples | primary solid tumor only (barcode field 4 starts with `01`); replicate aliquots averaged to the patient |
| gene map | `gencode.v36.annotation.gtf.gene.probemap` |
| keratin covariates | KRT8, KRT18, KRT19 |
| CD8 score | mean of CD8A and CD8B |
| cytotoxic score | mean of GZMA, GZMB, PRF1, NKG7 |
| purity sensitivity | PanCanAtlas ABSOLUTE purity, added as a fourth covariate where the patient has a call |
| not used | adjacent normal, metastatic samples, ICI-treated cohorts, spatial data |

Pan-cancer coexpression uses these epithelial projects: BLCA, BRCA, CESC, CHOL, COAD, ESCA, HNSC, KICH, KIRC, KIRP, LIHC, LUAD, LUSC, OV, PAAD, PRAD, READ, STAD, THCA, UCEC, UCS.

Keratin-adjusted immune exclusion is restricted to LUAD, LUSC, and the locked keratin funnel (LUAD, BRCA, CESC, KIRC, STAD, BLCA, PAAD). LUAD is in both the lung pair and the funnel; it is counted once.

This is not an ICI cohort. A negative partial correlation with a CD8 score is a bulk RNA association after a keratin control. It is not a spatial exclusion measurement.

## Methods

Marginal association is Spearman correlation. Keratin-adjusted association is partial Spearman: Pearson correlation of rank residuals after linear regression on an intercept plus the ranks of the covariates. The t test uses df = n − 2 − k. Confidence intervals use the Fisher z variance 1/(n − 3 − k).

Benjamini–Hochberg q-values are computed inside each pre-specified family:

- Coexpression primary family: TACSTD2–CLDN4 and TACSTD2–CLDN7 partial tests across all epithelial cohorts.
- Immune primary family: TACSTD2, CLDN4, and CLDN7 versus the CD8 score, partial tests, lung + funnel cohorts only.
- Cytotoxic family: the same three predictors versus the cytotoxic score, same cohorts.
- Epithelial comparators (EPCAM, MUC1) are reported and FDR-controlled in their own families, not mixed into the primary q-values.

Cross-cohort summary is a DerSimonian–Laird random-effects model on Fisher z. I² is reported next to the pooled ρ. Per-cohort signs are the primary read when I² is large.

The locked concordance model is different on purpose: CLDN4 versus KRT8, both residualized only on KRT18 and KRT19. KRT8 is the outcome, so it is not a covariate in that one test.

## Coexpression

TACSTD2–CLDN4 keratin-partial ρ is positive in 19/21 epithelial cohorts. Random-effects meta, all epithelial: ρ=0.315 (95% CI 0.251 to 0.376), p=3.82e-20, I²=89.1%, cohorts=21.

TACSTD2–CLDN7 keratin-partial ρ is positive in 19/21 epithelial cohorts. Random-effects meta, all epithelial: ρ=0.242 (95% CI 0.178 to 0.304), p=4.35e-13, I²=88.1%, cohorts=21.

Inside LUAD/LUSC plus the keratin funnel, TACSTD2–CLDN4 partial is positive in 8/8 and TACSTD2–CLDN7 partial is positive in 8/8.

The epithelial cohorts whose TACSTD2–CLDN4 keratin-partial ρ is not positive are: COAD, READ. For TACSTD2–CLDN7 they are: COAD, READ.

Keratin adjustment removes the shared KRT8/18/19 axis. The Trop2–claudin partial ρ stays positive in the lung pair and in every keratin-funnel cohort. In some cohorts the partial ρ is larger than the marginal ρ. COAD and READ are the epithelial projects where the partial ρ is null or slightly negative.

Same keratin model, epithelial comparators, all-epithelial random-effects meta: TACSTD2–EPCAM ρ=0.072 (95% CI -0.007 to 0.150), p=0.0751, I²=91.7%, cohorts=21; TACSTD2–MUC1 ρ=0.217 (95% CI 0.161 to 0.272), p=1.61e-13, I²=84.4%, cohorts=21. After KRT8/18/19, the CLDN4 association stays larger than the EPCAM association. That comparison does not reopen the locked surface-gene ranking.

| cohort | tier | n | TACSTD2_CLDN4_marginal | TACSTD2_CLDN4_partial | q_CLDN4 | TACSTD2_CLDN7_marginal | TACSTD2_CLDN7_partial | q_CLDN7 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LUAD | lung | 516 | 0.539 | 0.405 | 7.92e-21 | 0.140 | 0.110 | 0.0166 |
| LUSC | lung | 501 | 0.391 | 0.396 | 2.37e-19 | 0.372 | 0.357 | 9.86e-16 |
| BLCA | keratin_funnel | 406 | 0.566 | 0.507 | 2.35e-26 | 0.427 | 0.344 | 3.70e-12 |
| BRCA | keratin_funnel | 1095 | 0.330 | 0.313 | 4.48e-25 | 0.223 | 0.202 | 5.43e-11 |
| CESC | keratin_funnel | 304 | 0.284 | 0.277 | 2.16e-06 | 0.150 | 0.244 | 2.97e-05 |
| KIRC | keratin_funnel | 533 | 0.149 | 0.189 | 2.04e-05 | 0.241 | 0.169 | 0.0001 |
| PAAD | keratin_funnel | 178 | 0.709 | 0.515 | 1.07e-12 | 0.483 | 0.167 | 0.0325 |
| STAD | keratin_funnel | 412 | 0.326 | 0.215 | 1.92e-05 | 0.166 | 0.041 | 0.4372 |
| CHOL | extended_epithelial | 35 | 0.522 | 0.553 | 0.0015 | 0.287 | 0.370 | 0.0436 |
| COAD | extended_epithelial | 458 | 0.009 | -0.048 | 0.3377 | 0.060 | -0.012 | 0.8230 |
| ESCA | extended_epithelial | 184 | -0.159 | 0.165 | 0.0325 | -0.220 | 0.091 | 0.2535 |
| HNSC | extended_epithelial | 520 | 0.489 | 0.425 | 4.88e-23 | 0.492 | 0.481 | 1.13e-29 |
| KICH | extended_epithelial | 66 | 0.340 | 0.372 | 0.0035 | 0.354 | 0.416 | 0.0011 |
| KIRP | extended_epithelial | 290 | 0.537 | 0.289 | 1.38e-06 | 0.661 | 0.479 | 3.68e-17 |
| LIHC | extended_epithelial | 371 | 0.401 | 0.162 | 0.0026 | 0.275 | 0.161 | 0.0027 |
| OV | extended_epithelial | 421 | 0.324 | 0.279 | 1.60e-08 | 0.280 | 0.230 | 3.99e-06 |
| PRAD | extended_epithelial | 497 | 0.442 | 0.336 | 5.87e-14 | 0.332 | 0.270 | 2.86e-09 |
| READ | extended_epithelial | 166 | 0.031 | -0.040 | 0.6409 | 0.072 | -6.20e-04 | 0.9937 |
| THCA | extended_epithelial | 505 | 0.676 | 0.412 | 4.51e-21 | 0.558 | 0.283 | 2.85e-10 |
| UCEC | extended_epithelial | 545 | 0.405 | 0.329 | 1.48e-14 | 0.279 | 0.195 | 9.06e-06 |
| UCS | extended_epithelial | 57 | 0.835 | 0.637 | 5.17e-07 | 0.827 | 0.601 | 3.13e-06 |

EPCAM and MUC1 use the same KRT8/18/19 adjustment. Their per-cohort rows are in `tables/coexpression.tsv`; the pooled estimates are in the sentence above. A positive TACSTD2–CLDN partial ρ does not make CLDN4 the top surface partner. The locked surface-gene ranking is unchanged.

## Keratin-adjusted immune association

Negative ρ means higher predictor, lower immune score, after KRT8+KRT18+KRT19.

TACSTD2 versus CD8 score: negative in 6/8 lung+funnel cohorts. Meta: ρ=-0.069 (95% CI -0.135 to -0.002), p=0.0427, I²=76.1%, cohorts=8.

CLDN4 versus CD8 score: negative in 7/8. Meta: ρ=-0.081 (95% CI -0.138 to -0.023), p=0.0062, I²=68.1%, cohorts=8.

CLDN7 versus CD8 score: negative in 7/8. Meta: ρ=-0.109 (95% CI -0.173 to -0.045), p=0.0009, I²=74.9%, cohorts=8.

TACSTD2 versus cytotoxic score (GZMA, GZMB, PRF1, NKG7): negative in 5/8. Meta: ρ=-0.061 (95% CI -0.129 to 0.006), p=0.0762, I²=77.1%, cohorts=8.

| cohort | n | TACSTD2_marginal | TACSTD2_partial | q_TACSTD2 | CLDN4_marginal | CLDN4_partial | q_CLDN4 | CLDN7_partial | q_CLDN7 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LUAD | 516 | -0.113 | -0.104 | 0.0508 | -0.099 | -0.085 | 0.1100 | -0.267 | 2.00e-08 |
| LUSC | 501 | -0.302 | -0.221 | 7.31e-06 | -0.077 | -0.018 | 0.7320 | -0.046 | 0.3719 |
| BRCA | 1095 | -0.056 | -0.039 | 0.2813 | 0.010 | 0.034 | 0.3442 | -0.120 | 0.0005 |
| CESC | 304 | 0.069 | 0.118 | 0.0891 | -0.062 | -0.082 | 0.2381 | 4.93e-04 | 0.9932 |
| KIRC | 533 | -0.024 | -0.047 | 0.3513 | -0.108 | -0.066 | 0.2233 | -0.020 | 0.7320 |
| STAD | 412 | -0.132 | -0.079 | 0.2007 | -0.159 | -0.179 | 0.0017 | -0.121 | 0.0432 |
| BLCA | 406 | -0.200 | -0.158 | 0.0059 | -0.191 | -0.145 | 0.0123 | -0.070 | 0.2381 |
| PAAD | 178 | -0.139 | 0.029 | 0.7320 | -0.276 | -0.168 | 0.0620 | -0.249 | 0.0043 |

q-values in this immune table are from the CD8 primary family (TACSTD2, CLDN4, CLDN7 across the eight cohorts).

## Concordance with the locked CLDN4–KRT8 test

Both CLDN4 and KRT8 are residualized on KRT18 and KRT19 only. The locked public statement is that 6 of the 7 funnel cohorts are positive and BRCA is negative (partial ρ −0.071), with screen q-values 0.035 in LUAD and 0.0035 in STAD. This script does not rerun that surface-gene screen, so it does not emit those q-values.

In this run the funnel sign count is 7/7 positive, and the BRCA partial ρ is 0.202. That BRCA estimate is positive on primary-tumor STAR counts. It does not reproduce the locked screen's BRCA partial ρ of −0.071, and it does not replace it. LUAD and STAD single-test p-values in the table below are not the locked screen q-values (0.035 and 0.0035).

| cohort | n | rho_marginal | rho_partial | p_partial |
| --- | --- | --- | --- | --- |
| LUAD | 516 | 0.332 | 0.081 | 0.0655 |
| LUSC | 501 | 0.371 | 0.012 | 0.7887 |
| BRCA | 1095 | 0.222 | 0.202 | 1.48e-11 |
| CESC | 304 | 0.247 | 0.100 | 0.0826 |
| KIRC | 533 | 0.387 | 0.410 | 6.73e-23 |
| STAD | 412 | 0.510 | 0.126 | 0.0108 |
| BLCA | 406 | 0.493 | 0.106 | 0.0327 |
| PAAD | 178 | 0.604 | 0.042 | 0.5794 |
| CHOL | 35 | 0.195 | 0.006 | 0.9753 |
| COAD | 458 | 0.307 | 0.202 | 1.37e-05 |
| ESCA | 184 | 0.638 | 0.236 | 0.0013 |
| HNSC | 520 | 0.278 | -0.041 | 0.3561 |
| KICH | 66 | 0.434 | 0.169 | 0.1817 |
| KIRP | 290 | 0.629 | 0.410 | 4.21e-13 |
| LIHC | 371 | 0.391 | 0.155 | 0.0028 |
| OV | 421 | 0.205 | 0.052 | 0.2907 |
| PRAD | 497 | 0.250 | 0.176 | 7.93e-05 |
| READ | 166 | 0.289 | 0.227 | 0.0035 |
| THCA | 505 | 0.515 | 0.143 | 0.0013 |
| UCEC | 545 | 0.364 | 0.313 | 8.53e-14 |
| UCS | 57 | 0.745 | 0.262 | 0.0536 |

## Sensitivities

ABSOLUTE purity added on top of KRT8/18/19, lung + funnel, CD8 score. Patients without a purity call are dropped in this block only.

| cohort | x | n | rho_partial | p_partial |
| --- | --- | --- | --- | --- |
| LUAD | TACSTD2 | 503 | -0.056 | 0.2091 |
| LUAD | CLDN4 | 503 | -0.016 | 0.7162 |
| LUAD | CLDN7 | 503 | -0.150 | 0.0007 |
| LUSC | TACSTD2 | 493 | -0.233 | 1.79e-07 |
| LUSC | CLDN4 | 493 | -0.032 | 0.4823 |
| LUSC | CLDN7 | 493 | -0.040 | 0.3828 |
| BRCA | TACSTD2 | 1046 | -0.089 | 0.0040 |
| BRCA | CLDN4 | 1046 | 0.043 | 0.1615 |
| BRCA | CLDN7 | 1046 | -0.061 | 0.0485 |
| CESC | TACSTD2 | 291 | 0.112 | 0.0584 |
| CESC | CLDN4 | 291 | -0.054 | 0.3630 |
| CESC | CLDN7 | 291 | -0.002 | 0.9694 |
| KIRC | TACSTD2 | 498 | -0.058 | 0.1949 |
| KIRC | CLDN4 | 498 | 0.025 | 0.5737 |
| KIRC | CLDN7 | 498 | 0.012 | 0.7819 |
| STAD | TACSTD2 | 400 | -0.097 | 0.0543 |
| STAD | CLDN4 | 400 | -0.133 | 0.0082 |
| STAD | CLDN7 | 400 | -0.048 | 0.3440 |
| BLCA | TACSTD2 | 395 | -0.159 | 0.0016 |
| BLCA | CLDN4 | 395 | -0.176 | 0.0005 |
| BLCA | CLDN7 | 395 | -0.021 | 0.6826 |
| PAAD | TACSTD2 | 158 | 0.047 | 0.5602 |
| PAAD | CLDN4 | 158 | -0.114 | 0.1595 |
| PAAD | CLDN7 | 158 | -0.195 | 0.0154 |

LUSC squamous-keratin sensitivity: the same partial tests with covariates KRT5+KRT6A+KRT14 instead of KRT8/18/19. LUSC is a squamous carcinoma, so the simple-epithelial keratins are not its dominant keratin program. This row is a sensitivity, not a replacement of the pre-specified model.

| x | y | n | rho_marginal | rho_partial | p_partial |
| --- | --- | --- | --- | --- | --- |
| TACSTD2 | CLDN4 | 501 | 0.391 | 0.460 | 1.82e-27 |
| TACSTD2 | CLDN7 | 501 | 0.372 | 0.412 | 8.60e-22 |
| TACSTD2 | CD8_score | 501 | -0.302 | -0.138 | 0.0020 |
| CLDN4 | CD8_score | 501 | -0.077 | -0.116 | 0.0097 |
| CLDN7 | CD8_score | 501 | -0.099 | -0.111 | 0.0134 |

## Reading

TACSTD2 stays positively associated with CLDN4 and with CLDN7 after KRT8/18/19 adjustment in LUAD, LUSC, and every keratin-funnel cohort. The epithelial exceptions are the colorectal projects named above. The CD8 result is smaller and uneven: the pooled partial ρ is negative, I² is high, and cohort q-values show the signal is concentrated (LUSC for TACSTD2, LUAD for CLDN7, BLCA and STAD for parts of the trio) rather than shared by all eight cohorts. CESC is the cohort where TACSTD2’s keratin-adjusted CD8 association is positive.

What this does not say: it does not say CLDN4 is the top TACSTD2 surface partner, it does not say the association is spatial immune exclusion, and it does not say the genes rise after checkpoint blockade. Adjacent normal tissue was not mixed into the correlations. Private KL single-cell matrices are not in this repository and were not used.

## Reproduce

```bash
python3 scripts/tcga_trop2_cldn_keratin/test_stats.py
python3 scripts/tcga_trop2_cldn_keratin/run_analysis.py
```

The script streams each cohort matrix from the Xena GDC hub, caches patient-level gene tables under `/tmp/tcga_trop2` (or `TCGA_TROP2_CACHE`), and rewrites this file plus the TSV tables and figures. Provenance for the extracted tables, not the full matrices, is in `tables/provenance.json`.

## Gene IDs used

| symbol | ensembl |
| --- | --- |
| TACSTD2 | ENSG00000184292.7 |
| CLDN4 | ENSG00000189143.9 |
| CLDN7 | ENSG00000181885.18 |
| CLDN1 | ENSG00000163347.6 |
| CLDN3 | ENSG00000165215.6 |
| F11R | ENSG00000158769.18 |
| EPCAM | ENSG00000119888.11 |
| MUC1 | ENSG00000185499.16 |
| KRT8 | ENSG00000170421.12 |
| KRT18 | ENSG00000111057.11 |
| KRT19 | ENSG00000171345.13 |
| KRT5 | ENSG00000186081.12 |
| KRT6A | ENSG00000205420.11 |
| KRT14 | ENSG00000186847.6 |
| CD8A | ENSG00000153563.16 |
| CD8B | ENSG00000172116.23 |
| CD3D | ENSG00000167286.10 |
| CD3E | ENSG00000198851.10 |
| GZMA | ENSG00000145649.8 |
| GZMB | ENSG00000100453.13 |
| PRF1 | ENSG00000180644.8 |
| NKG7 | ENSG00000105374.10 |

## Sample counts

| cohort | tier | n_01_columns | n_patients |
| --- | --- | --- | --- |
| BLCA | keratin_funnel | 409 | 406 |
| BRCA | keratin_funnel | 1106 | 1095 |
| CESC | keratin_funnel | 304 | 304 |
| CHOL | extended_epithelial | 35 | 35 |
| COAD | extended_epithelial | 471 | 458 |
| ESCA | extended_epithelial | 184 | 184 |
| HNSC | extended_epithelial | 520 | 520 |
| KICH | extended_epithelial | 66 | 66 |
| KIRC | keratin_funnel | 537 | 533 |
| KIRP | extended_epithelial | 290 | 290 |
| LIHC | extended_epithelial | 371 | 371 |
| LUAD | lung | 528 | 516 |
| LUSC | lung | 501 | 501 |
| OV | extended_epithelial | 422 | 421 |
| PAAD | keratin_funnel | 178 | 178 |
| PRAD | extended_epithelial | 501 | 497 |
| READ | extended_epithelial | 166 | 166 |
| STAD | keratin_funnel | 412 | 412 |
| THCA | extended_epithelial | 505 | 505 |
| UCEC | extended_epithelial | 549 | 545 |
| UCS | extended_epithelial | 57 | 57 |

