# TCGA TACSTD2 versus CD8/CD3 after keratin

Primary tumors only (sample type 01), one row per patient, Xena GDC STAR log2(TPM+1). Eight cohorts: LUAD, LUSC, BRCA, CESC, KIRC, STAD, BLCA, PAAD. This file is written by `scripts/tcga_tacstd2_keratin_cd8/sweep.py` from the tables next to it.

## Objective

The grid is declared in `sweep.py` before the matrices are scored. Immune scores are CD3D, CD3E, CD3G, CD8A, CD8B, and the mean, z-mean, rank-mean, and first principal component of CD3, of CD8, and of CD3+CD8. Covariates are basal keratins (KRT5/6/14, the broader basal set, and that set plus TP63), simple keratins (KRT8/18/19, plus KRT7 and KRT20), z-means and principal components of those sets, the basal-plus-simple set, and the same keratin scores with ESTIMATE purity. Strata are all patients, basal and simple keratin median and tertile splits, and an ESTIMATE median split.

A keratin-adjusted specification has at least one covariate. It is eligible when the random-effects pooled association is negative, at least 6 of the 8 cohorts are negative, and at least 6 cohorts enter the pool. Untested cohorts are not counted as negative. The continuous objective maximizes |pooled ρ|. The Cliff objective maximizes |pooled δ|. Negative δ means the upper TACSTD2 group has the lower immune score. Correlations use a Fisher-z DerSimonian-Laird pool. Deltas use a DerSimonian-Laird pool on the delta scale with the Cliff 1993 variance.

The maximum is the extreme of this search. Selecting it pushes |effect| upward relative to a single pre-specified test. `q` is Benjamini-Hochberg on the meta p across every specification in that family, eligible or not.

`stratum_only` drops the covariate inside a keratin split. Those rows are stored and are excluded from the keratin-adjusted maximum, because that maximum is a partial association.

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

Continuous specifications: 10506. Quantile specifications: 63036. Keratin-adjusted eligible continuous: 6995. Keratin-adjusted eligible quantile: 32009.

Partial Spearman of CD8_mean and CD3_mean on KRT5+KRT6A+KRT6B+KRT14, all patients, was recomputed and compared with the previous KRT5/6 grid at 6 decimal places. Maximum absolute difference across 16 cohort-level TACSTD2 correlations: 4.153e-07.

## Pre-specified all-patient Spearman panel

These rows are fixed in the script. They are the keratin-adjusted CD8/CD3 associations on all patients, before any search for a larger effect.

| outcome | covariates | pooled | ci_low | ci_high | p | I2 | neg | cohorts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CD8A | KRT5_6_14 | -0.112 | -0.176 | -0.046 | 0.0008 | 0.753 | 7 | 8 |
| CD8A | KRT8_18_19 | -0.061 | -0.124 | 0.002 | 0.0581 | 0.734 | 6 | 8 |
| CD8A | both_programs | -0.073 | -0.124 | -0.021 | 0.0055 | 0.593 | 7 | 8 |
| CD8A | keratin_z | -0.077 | -0.145 | -0.009 | 0.0255 | 0.772 | 7 | 8 |
| CD8A | simple_z | -0.077 | -0.147 | -0.007 | 0.0307 | 0.787 | 7 | 8 |
| CD8A | basal_z | -0.111 | -0.180 | -0.041 | 0.0018 | 0.787 | 7 | 8 |
| CD8_mean | KRT5_6_14 | -0.117 | -0.177 | -0.056 | 0.0002 | 0.719 | 7 | 8 |
| CD8_mean | KRT8_18_19 | -0.069 | -0.135 | -0.002 | 0.0427 | 0.761 | 6 | 8 |
| CD8_mean | both_programs | -0.080 | -0.131 | -0.028 | 0.0024 | 0.596 | 7 | 8 |
| CD8_mean | keratin_z | -0.084 | -0.149 | -0.019 | 0.0118 | 0.757 | 7 | 8 |
| CD8_mean | simple_z | -0.087 | -0.159 | -0.015 | 0.0185 | 0.802 | 7 | 8 |
| CD8_mean | basal_z | -0.116 | -0.183 | -0.049 | 0.0008 | 0.772 | 7 | 8 |
| CD3D | KRT5_6_14 | -0.109 | -0.156 | -0.061 | 9.13e-06 | 0.540 | 7 | 8 |
| CD3D | KRT8_18_19 | -0.067 | -0.119 | -0.014 | 0.0136 | 0.619 | 6 | 8 |
| CD3D | both_programs | -0.086 | -0.122 | -0.049 | 5.61e-06 | 0.236 | 7 | 8 |
| CD3D | keratin_z | -0.077 | -0.130 | -0.024 | 0.0047 | 0.632 | 7 | 8 |
| CD3D | simple_z | -0.067 | -0.131 | -0.001 | 0.0456 | 0.753 | 6 | 8 |
| CD3D | basal_z | -0.108 | -0.163 | -0.052 | 0.0002 | 0.665 | 7 | 8 |
| CD3_mean | KRT5_6_14 | -0.102 | -0.155 | -0.047 | 0.0003 | 0.642 | 7 | 8 |
| CD3_mean | KRT8_18_19 | -0.053 | -0.109 | 0.004 | 0.0700 | 0.672 | 6 | 8 |
| CD3_mean | both_programs | -0.068 | -0.108 | -0.027 | 0.0010 | 0.350 | 7 | 8 |
| CD3_mean | keratin_z | -0.063 | -0.124 | -8.90e-04 | 0.0468 | 0.722 | 6 | 8 |
| CD3_mean | simple_z | -0.058 | -0.126 | 0.011 | 0.1003 | 0.777 | 6 | 8 |
| CD3_mean | basal_z | -0.100 | -0.163 | -0.037 | 0.0020 | 0.741 | 7 | 8 |
| CD3CD8_z | KRT5_6_14 | -0.112 | -0.166 | -0.057 | 7.38e-05 | 0.654 | 7 | 8 |
| CD3CD8_z | KRT8_18_19 | -0.062 | -0.123 | -3.38e-04 | 0.0488 | 0.719 | 6 | 8 |
| CD3CD8_z | both_programs | -0.076 | -0.118 | -0.033 | 0.0006 | 0.425 | 7 | 8 |
| CD3CD8_z | keratin_z | -0.074 | -0.136 | -0.011 | 0.0205 | 0.730 | 6 | 8 |
| CD3CD8_z | simple_z | -0.073 | -0.143 | -0.002 | 0.0446 | 0.792 | 6 | 8 |
| CD3CD8_z | basal_z | -0.111 | -0.174 | -0.047 | 0.0007 | 0.748 | 7 | 8 |

## Continuous maximum

Eligible all-patient Spearman specifications: 284. Median |ρ| 0.080, 95th percentile 0.116, maximum 0.120.

All-patient Spearman maximum:

`CD8_mean | KRT5_6_14_TP63 | all | spearman`.

| pooled | ci_low | ci_high | p | q | I2 | negative_cohorts | cohorts_in_pool | fixed_effect | k |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -0.120 | -0.180 | -0.058 | 0.0002 | 0.0019 | 0.724 | 7 | 8 | -0.123 | 5 |

Non-negative cohorts: CESC.

| cohort | n | effect | ci_low | ci_high | p |
| --- | --- | --- | --- | --- | --- |
| LUAD | 516 | -0.162 | -0.245 | -0.076 | 0.0002 |
| LUSC | 501 | -0.152 | -0.237 | -0.065 | 0.0007 |
| BRCA | 1095 | -0.146 | -0.204 | -0.087 | 1.28e-06 |
| CESC | 304 | 0.099 | -0.015 | 0.210 | 0.0886 |
| KIRC | 533 | -0.102 | -0.186 | -0.017 | 0.0190 |
| STAD | 412 | -0.053 | -0.150 | 0.044 | 0.2837 |
| BLCA | 406 | -0.150 | -0.245 | -0.053 | 0.0026 |
| PAAD | 178 | -0.305 | -0.435 | -0.164 | 4.38e-05 |

Eligible keratin-adjusted specifications in the full continuous grid: 6995. Median |ρ| 0.082, 95th percentile 0.140, maximum 0.170.

Full-grid maximum:

`CD3D | KRT5_6_14_TP63 | simple_top_tertile | pearson_winsor`.

| pooled | ci_low | ci_high | p | q | I2 | negative_cohorts | cohorts_in_pool | fixed_effect | k |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -0.170 | -0.238 | -0.100 | 2.69e-06 | 0.0002 | 0.356 | 8 | 8 | -0.172 | 5 |

Every tested cohort is negative.

| cohort | n | effect | ci_low | ci_high | p |
| --- | --- | --- | --- | --- | --- |
| LUAD | 172 | -0.122 | -0.269 | 0.031 | 0.1171 |
| LUSC | 167 | -0.089 | -0.240 | 0.066 | 0.2619 |
| BRCA | 365 | -0.166 | -0.265 | -0.064 | 0.0015 |
| CESC | 102 | -0.187 | -0.373 | 0.013 | 0.0662 |
| KIRC | 178 | -0.373 | -0.494 | -0.237 | 4.41e-07 |
| STAD | 138 | -0.149 | -0.312 | 0.021 | 0.0861 |
| BLCA | 136 | -0.139 | -0.304 | 0.033 | 0.1121 |
| PAAD | 60 | -0.036 | -0.299 | 0.231 | 0.7927 |

Largest eligible continuous specifications:

| outcome | covariates | stratum | method | abs_effect | effect | neg | cohorts | p | q |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CD3D | KRT5_6_14_TP63 | simple_top_tertile | pearson_winsor | 0.170 | -0.170 | 8 | 8 | 2.69e-06 | 0.0002 |
| CD3G | KRT5_6_14_TP63 | pur_low | pearson_winsor | 0.169 | -0.169 | 7 | 8 | 6.78e-05 | 0.0012 |
| CD3G | KRT5_6_14_TP63 | pur_low | spearman | 0.168 | -0.168 | 7 | 8 | 7.49e-05 | 0.0013 |
| CD3D | KRT5_6_14_TP63 | simple_top_tertile | pearson | 0.166 | -0.166 | 8 | 8 | 2.72e-06 | 0.0002 |
| CD3G | KRT5_6AB | pur_low | spearman | 0.164 | -0.164 | 7 | 8 | 0.0002 | 0.0020 |
| CD3CD8_pc1 | KRT5_6AB | pur_low | spearman | 0.164 | -0.164 | 7 | 8 | 5.64e-05 | 0.0011 |
| CD3CD8_z | KRT5_6AB | pur_low | spearman | 0.164 | -0.164 | 7 | 8 | 5.50e-05 | 0.0011 |
| CD3G | KRT5_6_14 | pur_low | pearson_winsor | 0.164 | -0.164 | 7 | 8 | 0.0002 | 0.0021 |
| CD3CD8_mean | KRT5_6AB | pur_low | spearman | 0.164 | -0.164 | 7 | 8 | 6.21e-05 | 0.0012 |
| CD3CD8_pc1 | KRT5_6_14 | pur_low | spearman | 0.163 | -0.163 | 7 | 8 | 8.59e-05 | 0.0014 |
| CD3CD8_z | KRT5_6_14 | pur_low | spearman | 0.163 | -0.163 | 7 | 8 | 8.40e-05 | 0.0013 |
| CD3CD8_pc1 | KRT5_KRT6A | pur_low | spearman | 0.162 | -0.162 | 7 | 8 | 2.83e-05 | 0.0008 |

Largest positive pooled ρ with at least 6 cohorts in the pool: 0.056 at `CD3E | KRT5_6_14_ESTIMATE | basal_top_tertile | pearson_winsor` (p=0.2007, 5 cohorts positive, 8 cohorts in the pool).

## Cliff's delta

Quartiles are the outer 25%. The immune score is the keratin residual. `residual` cuts TACSTD2 on that residual; `expression` cuts TACSTD2 on log2(TPM+1). Arms below 12 patients are untested.

Pre-specified all-patient quartile contrasts, Spearman residuals, cut 0.25:

| outcome | covariates | pooled | ci_low | ci_high | p | I2 | neg | cohorts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CD8A | KRT5_6_14 | -0.168 | -0.266 | -0.070 | 0.0008 | 0.728 | 7 | 8 |
| CD8A | KRT8_18_19 | -0.093 | -0.181 | -0.005 | 0.0377 | 0.647 | 6 | 8 |
| CD8A | both_programs | -0.106 | -0.179 | -0.034 | 0.0042 | 0.482 | 7 | 8 |
| CD8A | keratin_z | -0.117 | -0.220 | -0.013 | 0.0270 | 0.748 | 7 | 8 |
| CD8_mean | KRT5_6_14 | -0.173 | -0.265 | -0.082 | 0.0002 | 0.684 | 7 | 8 |
| CD8_mean | KRT8_18_19 | -0.105 | -0.204 | -0.007 | 0.0355 | 0.719 | 6 | 8 |
| CD8_mean | both_programs | -0.127 | -0.190 | -0.064 | 7.33e-05 | 0.313 | 7 | 8 |
| CD8_mean | keratin_z | -0.119 | -0.222 | -0.017 | 0.0226 | 0.745 | 7 | 8 |
| CD3D | KRT5_6_14 | -0.160 | -0.231 | -0.089 | 9.23e-06 | 0.462 | 7 | 8 |
| CD3D | KRT8_18_19 | -0.106 | -0.180 | -0.031 | 0.0053 | 0.500 | 6 | 8 |
| CD3D | both_programs | -0.126 | -0.183 | -0.069 | 1.54e-05 | 0.179 | 7 | 8 |
| CD3D | keratin_z | -0.118 | -0.200 | -0.036 | 0.0047 | 0.589 | 6 | 8 |
| CD3_mean | KRT5_6_14 | -0.152 | -0.235 | -0.070 | 0.0003 | 0.605 | 7 | 8 |
| CD3_mean | KRT8_18_19 | -0.081 | -0.161 | -0.001 | 0.0463 | 0.565 | 6 | 8 |
| CD3_mean | both_programs | -0.094 | -0.158 | -0.031 | 0.0038 | 0.325 | 7 | 8 |
| CD3_mean | keratin_z | -0.092 | -0.186 | 0.002 | 0.0547 | 0.687 | 6 | 8 |
| CD3CD8_z | KRT5_6_14 | -0.164 | -0.247 | -0.081 | 0.0001 | 0.608 | 7 | 8 |
| CD3CD8_z | KRT8_18_19 | -0.094 | -0.179 | -0.009 | 0.0302 | 0.618 | 6 | 8 |
| CD3CD8_z | both_programs | -0.114 | -0.169 | -0.058 | 6.25e-05 | 0.142 | 7 | 8 |
| CD3CD8_z | keratin_z | -0.107 | -0.203 | -0.011 | 0.0284 | 0.704 | 6 | 8 |

Eligible all-patient Spearman residual quartile contrasts: 273. Median |δ| 0.118, 95th percentile 0.171, maximum 0.185.

All-patient Spearman residual quartile maximum:

`CD8_mean | KRT5_6_14_TP63 | all | spearman | cut 0.25 | residual`.

| pooled | ci_low | ci_high | p | q | I2 | negative_cohorts | cohorts_in_pool | fixed_effect | k |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -0.185 | -0.285 | -0.086 | 0.0003 | 0.0042 | 0.735 | 7 | 8 | -0.194 | 5 |

Non-negative cohorts: CESC.

| cohort | n | n_low | n_high | effect | ci_low | ci_high | p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| LUAD | 516 | 129 | 129 | -0.249 | -0.386 | -0.111 | 0.0004 |
| LUSC | 501 | 126 | 126 | -0.246 | -0.383 | -0.108 | 0.0005 |
| BRCA | 1095 | 274 | 274 | -0.237 | -0.330 | -0.143 | 6.74e-07 |
| CESC | 304 | 76 | 76 | 0.176 | -0.005 | 0.357 | 0.0568 |
| KIRC | 533 | 134 | 134 | -0.115 | -0.253 | 0.023 | 0.1027 |
| STAD | 412 | 103 | 103 | -0.117 | -0.274 | 0.041 | 0.1465 |
| BLCA | 406 | 102 | 102 | -0.245 | -0.399 | -0.092 | 0.0017 |
| PAAD | 178 | 45 | 45 | -0.451 | -0.660 | -0.241 | 2.44e-05 |

Eligible keratin-adjusted quantile specifications: 32009. Median |δ| 0.138, 95th percentile 0.261, maximum 0.441.

Full-grid Cliff maximum:

`CD8_mean | simple_z | simple_top_tertile | pearson | cut 0.10 | residual`.

| pooled | ci_low | ci_high | p | q | I2 | negative_cohorts | cohorts_in_pool | fixed_effect | k |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -0.441 | -0.667 | -0.215 | 0.0001 | 0.0031 | 0.666 | 6 | 6 | -0.461 | 1 |

Every tested cohort is negative. Untested cohorts, not counted as negative: CESC, PAAD.

| cohort | n | n_low | n_high | effect | ci_low | ci_high | p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| LUAD | 172 | 18 | 18 | -0.488 | -0.842 | -0.134 | 0.0069 |
| LUSC | 167 | 17 | 17 | -0.121 | -0.520 | 0.277 | 0.5514 |
| BRCA | 365 | 37 | 37 | -0.160 | -0.423 | 0.103 | 0.2331 |
| CESC | 102 | 11 | 11 | NA | NA | NA | NA |
| KIRC | 178 | 18 | 18 | -0.716 | -0.980 | -0.452 | 1.01e-07 |
| STAD | 138 | 14 | 14 | -0.357 | -0.762 | 0.048 | 0.0839 |
| BLCA | 136 | 14 | 14 | -0.724 | -1.011 | -0.438 | 7.00e-07 |
| PAAD | 60 | 6 | 6 | NA | NA | NA | NA |

Largest |δ| among eligible specifications that keep all 8 cohorts in the pool:

`CD3G | KRT5_6_14_TP63 | pur_low | pearson | cut 0.20 | residual`.

| pooled | ci_low | ci_high | p | q | I2 | negative_cohorts | cohorts_in_pool | fixed_effect | k |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -0.293 | -0.455 | -0.130 | 0.0004 | 0.0054 | 0.775 | 7 | 8 | -0.298 | 5 |

Non-negative cohorts: STAD.

| cohort | n | n_low | n_high | effect | ci_low | ci_high | p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| LUAD | 257 | 52 | 52 | -0.217 | -0.433 | -1.81e-06 | 0.0500 |
| LUSC | 250 | 50 | 50 | -0.406 | -0.612 | -0.200 | 0.0001 |
| BRCA | 546 | 110 | 110 | -0.322 | -0.468 | -0.177 | 1.37e-05 |
| CESC | 152 | 31 | 31 | -0.330 | -0.614 | -0.046 | 0.0229 |
| KIRC | 266 | 54 | 54 | -0.199 | -0.415 | 0.017 | 0.0707 |
| STAD | 205 | 41 | 41 | 0.193 | -0.054 | 0.440 | 0.1252 |
| BLCA | 203 | 41 | 41 | -0.305 | -0.544 | -0.066 | 0.0123 |
| PAAD | 89 | 18 | 18 | -0.747 | -0.988 | -0.506 | 1.29e-09 |

Largest eligible Cliff specifications:

| outcome | covariates | stratum | method | cut | predictor_cut | abs_effect | effect | neg | cohorts | p | q |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CD8_mean | simple_z | simple_top_tertile | pearson | 0.10 | residual | 0.441 | -0.441 | 6 | 6 | 0.0001 | 0.0031 |
| CD8_z | simple_z | simple_top_tertile | pearson | 0.10 | residual | 0.437 | -0.437 | 6 | 6 | 0.0002 | 0.0034 |
| CD8_pc1 | simple_z | simple_top_tertile | pearson | 0.10 | residual | 0.437 | -0.437 | 6 | 6 | 0.0002 | 0.0034 |
| CD8B | simple_z | simple_top_tertile | pearson | 0.10 | residual | 0.431 | -0.431 | 6 | 6 | 0.0008 | 0.0082 |
| CD8_z | KRT5_6_14_TP63 | simple_top_tertile | pearson_winsor | 0.10 | residual | 0.430 | -0.430 | 6 | 6 | 5.41e-08 | 9.74e-05 |
| CD8_pc1 | KRT5_6_14_TP63 | simple_top_tertile | pearson_winsor | 0.10 | residual | 0.430 | -0.430 | 6 | 6 | 5.41e-08 | 9.74e-05 |
| CD8_rank | simple_z | simple_top_tertile | pearson | 0.10 | residual | 0.428 | -0.428 | 6 | 6 | 0.0001 | 0.0031 |
| CD8_mean | KRT5_6_14_TP63 | simple_top_tertile | pearson_winsor | 0.10 | residual | 0.425 | -0.425 | 6 | 6 | 5.87e-08 | 1.00e-04 |
| CD3D | KRT5_6_14_TP63 | simple_top_tertile | pearson_winsor | 0.10 | residual | 0.422 | -0.422 | 6 | 6 | 1.24e-09 | 3.30e-05 |
| CD8B | KRT5_6_14_TP63 | simple_top_tertile | pearson_winsor | 0.10 | residual | 0.420 | -0.420 | 6 | 6 | 1.89e-06 | 0.0005 |
| CD3CD8_z | KRT5_6_14_TP63 | simple_top_tertile | pearson_winsor | 0.10 | residual | 0.419 | -0.419 | 6 | 6 | 1.38e-09 | 3.30e-05 |
| CD3CD8_z | simple_z | simple_top_tertile | pearson | 0.10 | residual | 0.419 | -0.419 | 6 | 6 | 1.07e-05 | 0.0011 |

## Reading

On all patients, partial Spearman of CD8_mean on KRT5+KRT6A+KRT6B+KRT14 is pooled ρ -0.117 (7/8 negative, p=0.0002). The same CD8_mean score on KRT8+KRT18+KRT19 is pooled ρ -0.069 (6/8 negative, p=0.0427). The largest eligible all-patient Spearman |ρ| is 0.120 at `CD8_mean | KRT5_6_14_TP63 | all | spearman` (pooled ρ -0.120, 95% CI -0.180 to -0.058, p=0.0002, q=0.0019, 7/8 negative). The keratin-adjusted continuous maximum is |ρ| 0.170 at `CD3D | KRT5_6_14_TP63 | simple_top_tertile | pearson_winsor` (pooled ρ -0.170, 95% CI -0.238 to -0.100, p=2.69e-06, q=0.0002, 8/8 negative, I²=0.356). The same outcome, covariates, and stratum without winsorization are Spearman pooled ρ -0.159 and Pearson pooled ρ -0.166, both with 8/8 and 8/8 cohorts negative. The all-patient Spearman residual quartile maximum is |δ| 0.185 at `CD8_mean | KRT5_6_14_TP63 | all | spearman | cut 0.25 | residual` (pooled δ -0.185, 7/8 negative, p=0.0003). The keratin-adjusted Cliff maximum is |δ| 0.441 at `CD8_mean | simple_z | simple_top_tertile | pearson | cut 0.10 | residual` (pooled δ -0.441, 95% CI -0.667 to -0.215, p=0.0001, q=0.0031, 6/8 negative).

These are bulk RNA associations in TCGA primary tumors after the named keratin covariates. The maximum is the extreme of the pre-declared grid.

```bash
python3 scripts/tcga_tacstd2_keratin_cd8/test_stats.py
python3 scripts/tcga_tacstd2_keratin_cd8/sweep.py
```
