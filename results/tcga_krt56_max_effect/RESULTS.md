# TCGA TACSTD2 / CLDN4 / CLDN7 versus CD3 / CD8 after KRT5/6

Primary tumors only (sample type 01), one row per patient, Xena GDC STAR log2(TPM+1). Eight cohorts: LUAD, LUSC, BRCA, CESC, KIRC, STAD, BLCA, PAAD. This file is written by `scripts/tcga_krt56_max_effect/sweep.py` from the tables next to it.

## Objective

The grid below was declared in the script before it was evaluated. A specification is eligible only when TACSTD2, CLDN4, and CLDN7 each have a negative random-effects pooled effect and each is negative in at least 7 of the 8 cohorts. The reported maximum maximizes the minimum of those three absolute pooled effects. Correlations use a Fisher-z DerSimonian-Laird pool. Quartile and tail contrasts use Cliff's delta (negative means the upper predictor group has the lower immune score), pooled on the delta scale with the Cliff 1993 variance.

The maximum is the extreme of this search. Selecting it pushes |effect| upward relative to a single pre-specified test. `q_worst` is Benjamini-Hochberg on the worst of the three meta p-values across every specification in that family, eligible or not. It does not undo that selection.

Every covariate set contains KRT5 and a KRT6 gene, except `squamous_PC1` (first principal component of KRT5, KRT6A, KRT6B, and KRT14, oriented with KRT5) and `stratum_only` (allowed only inside a KRT5/6 median or tertile split, where the split itself is the KRT5/6 control). Cytotoxic-only gene sets are not in the grid.

## Data

| cohort | n_patients | n_with_ESTIMATE |
| --- | --- | --- |
| LUAD | 516 | 515 |
| LUSC | 501 | 501 |
| BRCA | 1095 | 1093 |
| CESC | 304 | 304 |
| KIRC | 533 | 533 |
| STAD | 412 | 411 |
| BLCA | 406 | 406 |
| PAAD | 178 | 178 |

Continuous specifications: 9483. Quantile specifications: 56898. Eligible continuous: 2590. Eligible quantile: 9171.

Immune scores: CD3D, CD3E, CD3G, CD8A, CD8B; arithmetic mean, within-sample z-mean, rank-mean, and first principal component of CD3 (CD3D/E/G), CD8 (CD8A/B), CD3+CD8, the Danaher T-cell set (CD6, CD3D, CD3E, SH2D1A, TRAT1, CD3G), signature-H (15 T-cell genes), and Tsig (CD3D/E/G, CD2, CD247, LCK). CD8_mean is the two-gene Danaher CD8 score.

Strata, defined inside each cohort: all patients; KRT5/6 score (mean of KRT5, KRT6A, KRT6B) at or above the median, below the median, top tertile, bottom tertile; ESTIMATE purity at or above the median and below the median. ESTIMATE is the Yoshihara cosine transform of the MD Anderson RNAseqV2 ESTIMATE score.

## Pipeline check against the previous KRT5/6 panel

Partial Spearman of CD3_mean and CD8_mean on KRT5+KRT6A+KRT6B+KRT14, all patients, was recomputed and compared with the previous grid at 6 decimal places. Maximum absolute difference across 48 cohort-level correlations: 4.897e-07.

That previous CD8_mean specification, which is inside this grid, has minimum |pooled ρ| = 0.112 (TACSTD2, CLDN4, CLDN7 pooled ρ -0.117, -0.112, -0.142). TACSTD2 is the gene that misses a negative sign, in CESC. CD3_mean on the same covariates has a smaller minimum |ρ| because TACSTD2 is weaker.

## Continuous maximum

Among eligible continuous specifications restricted to all patients, the minimum |pooled ρ| has median 0.099, 95th percentile 0.112, and maximum 0.119 (595 eligible).

All-patient maximum: `CD3CD8_z | KRT5_6_14_TP63 | all | spearman`.

| predictor | pooled | ci_low | ci_high | p | I2 | negative_cohorts | cohorts_in_pool | fixed_effect |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TACSTD2 | -0.119 | -0.176 | -0.062 | 4.98e-05 | 0.681 | 7 | 8 | -0.121 |
| CLDN4 | -0.134 | -0.194 | -0.073 | 1.68e-05 | 0.716 | 8 | 8 | -0.110 |
| CLDN7 | -0.158 | -0.241 | -0.073 | 0.0003 | 0.858 | 8 | 8 | -0.154 |

Worst meta p = 0.0003. Sweep q = 0.0027.

| cohort | n | TACSTD2 | CLDN4 | CLDN7 |
| --- | --- | --- | --- | --- |
| LUAD | 516 | -0.169 | -0.160 | -0.318 |
| LUSC | 501 | -0.131 | -0.137 | -0.146 |
| BRCA | 1095 | -0.144 | -0.020 | -0.174 |
| CESC | 304 | 0.087 | -0.114 | -0.006 |
| KIRC | 533 | -0.083 | -0.131 | -0.003 |
| STAD | 412 | -0.107 | -0.090 | -0.062 |
| BLCA | 406 | -0.124 | -0.132 | -0.181 |
| PAAD | 178 | -0.305 | -0.373 | -0.385 |

Across the full continuous grid, including KRT5/6 and purity strata, 2590 specifications are eligible. Median minimum |ρ| = 0.100, 95th percentile = 0.128, maximum = 0.152.

Grid maximum: `CD8A | KRT5_6_14_stroma | krt_bot_tertile | pearson_winsor`.

| predictor | pooled | ci_low | ci_high | p | I2 | negative_cohorts | cohorts_in_pool | fixed_effect |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TACSTD2 | -0.152 | -0.222 | -0.079 | 4.27e-05 | 0.379 | 7 | 8 | -0.151 |
| CLDN4 | -0.153 | -0.226 | -0.078 | 6.60e-05 | 0.416 | 7 | 8 | -0.133 |
| CLDN7 | -0.169 | -0.223 | -0.115 | 1.54e-09 | 0.000 | 8 | 8 | -0.169 |

Worst meta p = 6.60e-05. Sweep q = 0.0017.

| cohort | n | TACSTD2 | CLDN4 | CLDN7 |
| --- | --- | --- | --- | --- |
| LUAD | 172 | -0.104 | -0.173 | -0.202 |
| LUSC | 167 | -0.220 | -0.138 | -0.128 |
| BRCA | 365 | -0.133 | 0.007 | -0.173 |
| CESC | 102 | -0.215 | -0.265 | -0.228 |
| KIRC | 178 | -0.141 | -0.195 | -0.158 |
| STAD | 138 | 0.004 | -0.130 | -0.062 |
| BLCA | 136 | -0.342 | -0.220 | -0.192 |
| PAAD | 60 | -0.004 | -0.277 | -0.302 |

Largest eligible continuous specifications, ordered by the minimum |pooled ρ|:

| outcome | covariates | stratum | method | min_abs | TACSTD2 | CLDN4 | CLDN7 | neg_T2 | neg_C4 | neg_C7 | worst_p | q_worst |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CD8A | KRT5_6_14_stroma | krt_bot_tertile | pearson_winsor | 0.152 | -0.152 | -0.153 | -0.169 | 7 | 7 | 8 | 6.60e-05 | 0.0017 |
| CD8A | KRT5_6_14_TP63 | krt_bot_tertile | spearman | 0.150 | -0.150 | -0.176 | -0.171 | 7 | 7 | 8 | 0.0006 | 0.0040 |
| CD3CD8_z | KRT5_6_14_TP63 | krt_bot_tertile | spearman | 0.150 | -0.150 | -0.180 | -0.184 | 8 | 7 | 8 | 0.0002 | 0.0021 |
| CD3CD8_pc1 | KRT5_6_14_TP63 | krt_bot_tertile | spearman | 0.150 | -0.150 | -0.180 | -0.183 | 8 | 7 | 8 | 0.0001 | 0.0021 |
| CD3CD8_mean | KRT5_6_14_TP63 | krt_bot_tertile | spearman | 0.148 | -0.148 | -0.179 | -0.183 | 8 | 7 | 8 | 0.0001 | 0.0021 |
| CD3G | KRT5_6_14_TP63 | krt_bot_tertile | spearman | 0.146 | -0.146 | -0.171 | -0.169 | 8 | 8 | 7 | 0.0002 | 0.0025 |
| CD8A | KRT5_6_14_stroma | krt_bot_tertile | pearson | 0.146 | -0.149 | -0.146 | -0.160 | 7 | 7 | 8 | 0.0003 | 0.0026 |
| CD8A | KRT5_6_14_stroma | krt_bot_tertile | spearman | 0.145 | -0.145 | -0.157 | -0.157 | 7 | 7 | 8 | 0.0005 | 0.0037 |
| CD8_mean | KRT5_6_14_TP63 | krt_bot_tertile | spearman | 0.145 | -0.145 | -0.155 | -0.170 | 7 | 7 | 8 | 0.0012 | 0.0061 |
| CD3CD8_rank | KRT5_6_14_TP63 | krt_bot_tertile | spearman | 0.144 | -0.144 | -0.176 | -0.183 | 8 | 7 | 8 | 0.0001 | 0.0021 |
| DanaherT_mean | KRT5_6_14_TP63 | krt_bot_tertile | spearman | 0.144 | -0.144 | -0.181 | -0.190 | 8 | 8 | 8 | 4.36e-05 | 0.0016 |
| CD3CD8_z | KRT5_6_14_stroma | krt_bot_tertile | spearman | 0.144 | -0.144 | -0.157 | -0.168 | 7 | 7 | 8 | 0.0001 | 0.0020 |

## Cliff's delta maximum

Quartiles are the outer 25% (cut 0.25). The grid also contains outer 20% and outer 10% cuts, because those tail contrasts are where |Cliff's delta| can grow when the association is monotone. Groups are formed either on raw log2(TPM+1) (`expression`) or on the same residual used for the immune score (`residual`). The immune score is always the residual after the named covariates. Arms smaller than 12 patients are not tested.

Eligible quartile contrasts in all patients: 979. Median minimum |pooled δ| = 0.135, 95th percentile = 0.166, maximum = 0.185.

All-patient quartile maximum: `CD8_mean | KRT5_6_14_TP63 | all | spearman | cut 0.25 | residual`.

| predictor | pooled | ci_low | ci_high | p | I2 | negative_cohorts | cohorts_in_pool | fixed_effect |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TACSTD2 | -0.185 | -0.285 | -0.086 | 0.0003 | 0.735 | 7 | 8 | -0.194 |
| CLDN4 | -0.193 | -0.295 | -0.092 | 0.0002 | 0.742 | 8 | 8 | -0.166 |
| CLDN7 | -0.194 | -0.302 | -0.086 | 0.0004 | 0.778 | 7 | 8 | -0.202 |

Worst meta p = 0.0004. Sweep q = 0.0065.

| cohort | n | n_low | n_high | TACSTD2 | CLDN4 | CLDN7 |
| --- | --- | --- | --- | --- | --- | --- |
| LUAD | 516 | 129 | 129 | -0.249 | -0.297 | -0.367 |
| LUSC | 501 | 126 | 126 | -0.246 | -0.185 | -0.174 |
| BRCA | 1095 | 274 | 274 | -0.237 | -0.012 | -0.252 |
| CESC | 304 | 76 | 76 | 0.176 | -0.108 | 0.012 |
| KIRC | 533 | 134 | 134 | -0.115 | -0.202 | -0.029 |
| STAD | 412 | 103 | 103 | -0.117 | -0.057 | -0.038 |
| BLCA | 406 | 102 | 102 | -0.245 | -0.271 | -0.205 |
| PAAD | 178 | 45 | 45 | -0.451 | -0.491 | -0.510 |

Across every cut and stratum, 9171 quantile specifications are eligible. Median minimum |pooled δ| = 0.160, 95th percentile = 0.212, maximum = 0.266.

Grid maximum: `CD3D | KRT5_6_14_S100A8_A9 | krt_low | spearman | cut 0.10 | residual`.

| predictor | pooled | ci_low | ci_high | p | I2 | negative_cohorts | cohorts_in_pool | fixed_effect |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TACSTD2 | -0.266 | -0.402 | -0.130 | 0.0001 | 0.302 | 7 | 7 | -0.266 |
| CLDN4 | -0.271 | -0.393 | -0.150 | 1.28e-05 | 0.149 | 7 | 7 | -0.268 |
| CLDN7 | -0.271 | -0.382 | -0.159 | 2.00e-06 | 0.000 | 7 | 7 | -0.271 |

Worst meta p = 0.0001. Sweep q = 0.0033.

| cohort | n | n_low | n_high | TACSTD2 | CLDN4 | CLDN7 |
| --- | --- | --- | --- | --- | --- | --- |
| LUAD | 258 | 26 | 26 | -0.272 | -0.104 | -0.399 |
| LUSC | 250 | 25 | 25 | -0.274 | -0.443 | -0.280 |
| BRCA | 547 | 55 | 55 | -0.227 | -0.106 | -0.237 |
| CESC | 152 | 16 | 16 | -0.609 | -0.234 | -0.188 |
| KIRC | 266 | 27 | 27 | -0.067 | -0.471 | -0.147 |
| STAD | 206 | 21 | 21 | -0.016 | -0.234 | -0.102 |
| BLCA | 203 | 21 | 21 | -0.370 | -0.351 | -0.478 |
| PAAD | 89 | 9 | 9 | NA | NA | NA |

Largest eligible Cliff's delta specifications:

| outcome | covariates | stratum | method | cut | predictor_cut | min_abs | TACSTD2 | CLDN4 | CLDN7 | neg_T2 | neg_C4 | neg_C7 | worst_p | q_worst |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CD3D | KRT5_6_14_S100A8_A9 | krt_low | spearman | 0.100 | residual | 0.266 | -0.266 | -0.271 | -0.271 | 7 | 7 | 7 | 0.0001 | 0.0033 |
| CD8A | KRT5_6_14_KRT13 | all | spearman | 0.100 | residual | 0.265 | -0.265 | -0.272 | -0.325 | 8 | 8 | 8 | 0.0004 | 0.0061 |
| CD8A | KRT5_6_14_CD68 | krt_low | spearman | 0.100 | residual | 0.257 | -0.257 | -0.268 | -0.279 | 7 | 7 | 7 | 0.0002 | 0.0040 |
| CD3D | KRT5_6_14_KRT13 | all | spearman | 0.100 | residual | 0.255 | -0.255 | -0.286 | -0.334 | 8 | 8 | 8 | 9.21e-05 | 0.0029 |
| SigH_rank | KRT5_6_14_CD68 | pur_low | pearson | 0.200 | residual | 0.255 | -0.291 | -0.258 | -0.255 | 7 | 7 | 7 | 0.0013 | 0.0116 |
| SigH_rank | KRT5_6_14_CD68 | pur_low | pearson_winsor | 0.200 | residual | 0.254 | -0.273 | -0.261 | -0.254 | 7 | 7 | 7 | 0.0014 | 0.0124 |
| CD8_mean | KRT5_6_14_CD68 | all | spearman | 0.100 | residual | 0.253 | -0.253 | -0.294 | -0.315 | 8 | 8 | 7 | 0.0006 | 0.0075 |
| SigH_z | KRT5_6_14_CD68 | pur_low | pearson | 0.200 | residual | 0.252 | -0.286 | -0.256 | -0.252 | 7 | 7 | 7 | 0.0011 | 0.0105 |
| SigH_pc1 | KRT5_6_14_CD68 | pur_low | pearson | 0.200 | residual | 0.252 | -0.286 | -0.254 | -0.252 | 7 | 7 | 7 | 0.0016 | 0.0131 |
| CD8_mean | KRT5_6_14_KRT13 | all | spearman | 0.100 | residual | 0.249 | -0.262 | -0.249 | -0.307 | 8 | 8 | 8 | 0.0007 | 0.0084 |
| CD3CD8_mean | KRT5_6_14_KRT13 | all | spearman | 0.100 | residual | 0.249 | -0.249 | -0.282 | -0.352 | 8 | 8 | 8 | 5.51e-05 | 0.0022 |
| SigH_z | KRT5_6_14_CD68 | pur_low | pearson_winsor | 0.200 | residual | 0.248 | -0.269 | -0.257 | -0.248 | 7 | 7 | 7 | 0.0011 | 0.0102 |

## Reading

The strongest all-patient partial correlation that keeps all three predictors negative in at least 7 cohorts is CD3CD8_z | KRT5_6_14_TP63 | all | spearman, with pooled ρ TACSTD2 -0.119, CLDN4 -0.134, CLDN7 -0.158 (minimum |ρ| 0.119). Allowing the pre-declared strata raises the minimum |pooled ρ| to 0.152 at CD8A | KRT5_6_14_stroma | krt_bot_tertile | pearson_winsor. The strongest all-patient quartile contrast on the same rule has pooled Cliff's delta TACSTD2 -0.185, CLDN4 -0.193, CLDN7 -0.194 (minimum |δ| 0.185) at CD8_mean | KRT5_6_14_TP63 | all | spearman | cut 0.25 | residual. The largest Cliff's delta in the full cut-and-stratum grid has minimum |δ| 0.266 at CD3D | KRT5_6_14_S100A8_A9 | krt_low | spearman | cut 0.10 | residual. That Cliff maximum does not include every cohort: an arm below the pre-set floor of 12 patients is left untested and is not counted as negative. The next Cliff specification keeps a negative sign in all 8 cohorts: CD8A | KRT5_6_14_KRT13 | all | spearman | cut 0.10 | residual, minimum |δ| 0.265 (TACSTD2 -0.265, CLDN4 -0.272, CLDN7 -0.325). On all patients, the balanced continuous effect moves only from the previous CD8_mean minimum |ρ| of 0.112 to 0.119. The three pooled correlations remain small. The continuous grid maximum minimum |ρ| is 0.152, and it is a KRT5/6-stratum result rather than an all-patient result.

These are bulk RNA associations after a keratin adjustment. They are not a spatial exclusion measurement, not a knockdown result, and not an estimate of what a single pre-specified test would have returned.
