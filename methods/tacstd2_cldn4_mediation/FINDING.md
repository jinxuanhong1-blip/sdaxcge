# TACSTD2–immune partial on CLDN4, and mediation by the locked 221-gene signature

Public bulk RNA only. OncoSG is East-Asian surgical LUAD (Chen et al., Nat Genet 2020). The GEO cohorts are the PR 590 LUAD bulks. These correlations are not a spatial exclusion result and not an ICI-response result. GSE10072, GSE11969, and GSE248378 stay closed. The 221-gene list is not refit.

## Question

TACSTD2 and CLDN4 are positively correlated in these tumors, and each is inversely correlated with CD8A. The locked signature is the CLDN4-high malignant program from PR 590 (concordant-4 Q4 vs Q1, then TCGA partial vs CLDN4 given KRT8/18/19 and ABSOLUTE purity). CLDN4 and TACSTD2 are not members of the 221 genes. Two quantities are estimated:

1. Partial Spearman of TACSTD2 vs CD8A and vs ImmuneScore given CLDN4.
2. Single-mediator model: TACSTD2 → mediator → immune, with mediator = the 221-gene score. The CLDN4 gene is the comparison mediator.

## Estimators

Spearman is two-sided. Partial Spearman is the Pearson correlation of rank residuals, with t degrees of freedom n−2−k and Fisher interval using SE = 1/√(n−3−k). That is the OncoSG estimator from PR 139 / PR 333 / PR 590.

Mediation uses average ranks, then z-scores those ranks, then OLS. With one mediator and the same covariates in every equation, and no TACSTD2×mediator interaction, the indirect effect a×b equals the total coefficient c minus the direct coefficient c′. Uncertainty for a×b and for the ratio a×b/c is a percentile interval from 5000 row bootstraps. The bootstrap p for a×b is twice the fraction of bootstrap draws on the other side of zero. A Sobel p that ignores covariance between a and b is stored in the table and is not the interval used in the text. Proportion mediated is a×b/c. It is unstable when c is near zero; the indirect effect is the quantity that is meta-analyzed.

ImmuneScore is the mean of within-cohort z-scores of CD8A, GZMA, GZMB, IFNG, EOMES, CXCL9, CXCL10, and TBX21, the same 8 genes as PR 590. OncoSG also reports the mean of the deposited z-scores, which is the PR 139 immune score, as a QC row. Composition adjustment is published PURITY on OncoSG and the PR 590 9-gene stromal mean on GEO (FAP, COL1A1, COL1A2, COL3A1, DCN, LUM, PDGFRA, TAGLN, ACTA2). GEO has no ABSOLUTE purity.

A cohort is called partial mediation when the total effect, the a path, and the indirect effect all have bootstrap intervals excluding zero, in the observed directions (TACSTD2 positively associated with the mediator; total and indirect effects inverse), the proportion is between 0 and 1, and the direct effect stays inverse. If the direct-effect interval includes zero, the call is indirect with a direct interval that includes zero. `total_not_inverse` means the total-effect bootstrap interval is not entirely below zero; the point estimate can still be negative. That is a statistical description of these tumors. It is not an intervention.

## Locked inputs and QC

Signature file SHA-256 `9d398e71f8a2569d3c5c7bb940577a24eaca1f19db0a1d589c84e13a0bac6a3a`, 221 genes, CLDN4/TACSTD2/KRT8/KRT18/KRT19/CD8A/CD274 absent. OncoSG matrix columns: 169. Signature genes present: OncoSG 215/221.

| check | this run | locked | |Δ| |
|---|---:|---:|---:|
| TACSTD2 CD8A unadj | -0.380 | -0.380 | 0.000 |
| TACSTD2 CD8A purity | -0.309 | -0.309 | 0.000 |
| TACSTD2 ImmuneScore_deposited unadj | -0.387 | -0.387 | 0.000 |
| TACSTD2 ImmuneScore_deposited purity | -0.318 | -0.318 | 0.000 |
| CLDN4 CD8A unadj | -0.416 | -0.416 | 0.000 |
| CLDN4 CD8A purity | -0.285 | -0.285 | 0.000 |
| CLDN4 TACSTD2 unadj | +0.505 | +0.505 | 0.000 |
| signature CD8A unadj | -0.578 | -0.578 | 0.000 |

Every locked rho is within 0.015 of the published value. The partials and mediation rows below use this same matrix and this same 221-gene score.

## Sample sizes

| cohort | n | signature genes used | ImmuneScore genes | composition |
|---|---:|---:|---:|---|
| OncoSG | 169 | 215/221 | 8/8 | PURITY |
| GSE273377 discovery | 103 | 221/221 | 8/8 | stromal |
| GSE273377 validation | 60 | 221/221 | 8/8 | stromal |
| GSE282774 | 58 | 221/221 | 8/8 | stromal |
| GSE233774 tumor | 30 | 217/221 | 7/8 | stromal |

GSE233774 ImmuneScore uses 7 of 8 genes: IFNG is absent from the FPKM file. The other cohorts have 8 of 8.

GSE273377 is one stage I LUAD study with a discovery stratum (GPL30173) and a validation stratum (GPL16791). Samples with `passed qc: FALSE` are out. GSE233774 uses tumor columns only. GSE282774 is pN2 LUAD. GSE273377 enters cross-study combinations as one study after its two strata are combined.

## Partial Spearman

Primary covariate is CLDN4. The signature column is the same partial with the 221-gene score as the covariate. Unadjusted TACSTD2–CD8A on OncoSG is the PR 139 anchor.

| cohort | endpoint | n | unadjusted ρ (p) | partial \| CLDN4 | partial \| signature | partial \| CLDN4+signature | partial \| composition | partial \| composition+CLDN4 |
|---|---|---:|---|---|---|---|---|---|
| OncoSG | CD8A | 169 | -0.380 (3.38e-07) | -0.217 (0.00471) | -0.284 (0.000193) | -0.228 (0.003) | -0.309 (4.69e-05) | -0.208 (0.00696) |
| OncoSG | ImmuneScore | 169 | -0.387 (1.98e-07) | -0.217 (0.00463) | -0.289 (0.000141) | -0.233 (0.00241) | -0.318 (2.70e-05) | -0.208 (0.00702) |
| GSE273377 discovery | CD8A | 103 | +0.039 (0.696) | -0.023 (0.82) | +0.154 (0.122) | +0.009 (0.929) | -0.028 (0.782) | -0.047 (0.642) |
| GSE273377 discovery | ImmuneScore | 103 | +0.096 (0.332) | +0.014 (0.887) | +0.223 (0.0245) | +0.054 (0.595) | +0.030 (0.763) | -0.008 (0.935) |
| GSE273377 validation | CD8A | 60 | -0.423 (0.000753) | -0.265 (0.0425) | -0.389 (0.00231) | -0.327 (0.0123) | -0.442 (0.00046) | -0.279 (0.0342) |
| GSE273377 validation | ImmuneScore | 60 | -0.330 (0.0101) | -0.265 (0.0424) | -0.284 (0.0293) | -0.350 (0.00701) | -0.355 (0.00578) | -0.304 (0.0203) |
| GSE282774 | CD8A | 58 | -0.206 (0.121) | -0.205 (0.126) | -0.131 (0.331) | -0.163 (0.229) | -0.205 (0.126) | -0.196 (0.148) |
| GSE282774 | ImmuneScore | 58 | -0.179 (0.179) | -0.158 (0.242) | -0.062 (0.648) | -0.092 (0.499) | -0.178 (0.184) | -0.139 (0.307) |
| GSE233774 tumor | CD8A | 30 | -0.437 (0.0158) | -0.341 (0.0699) | -0.364 (0.0523) | -0.327 (0.0889) | -0.266 (0.163) | -0.397 (0.0363) |
| GSE233774 tumor | ImmuneScore | 30 | -0.439 (0.0152) | -0.313 (0.0984) | -0.340 (0.0708) | -0.296 (0.126) | -0.268 (0.16) | -0.358 (0.0614) |

OncoSG TACSTD2 vs CD8A is -0.380 (p = 3.38e-07, n = 169). Partial on CLDN4 it is -0.217 (p = 0.00471). Partial on the 221-gene score it is -0.284 (p = 0.000193). Partial on CLDN4 and the score together it is -0.228 (p = 0.003). Partial on published PURITY it is -0.309 (p = 4.69e-05). Partial on PURITY and CLDN4 it is -0.208 (p = 0.00696).

OncoSG TACSTD2 vs ImmuneScore is -0.387 (p = 1.98e-07). Partial on CLDN4 it is -0.217 (p = 0.00463). Partial on the 221-gene score it is -0.289 (p = 0.000141). Partial on PURITY and CLDN4 it is -0.208 (p = 0.00702).

Cross-study combination (OncoSG, GSE273377 as one study, GSE282774, GSE233774):

| endpoint | covariate | model | k | n sum | ρ | 95% CI | p | I² |
|---|---|---|---:|---:|---:|---|---:|---:|
| CD8A | unadjusted | DL random | 4 | 420 | -0.292 | -0.405 to -0.171 | 4.24e-06 | 34% |
| CD8A | CLDN4 | fixed | 4 | 420 | -0.190 | -0.282 to -0.095 | 0.000108 | 0% |
| CD8A | signature | DL random | 4 | 420 | -0.206 | -0.307 to -0.101 | 0.000151 | 13% |
| ImmuneScore | unadjusted | DL random | 4 | 420 | -0.268 | -0.427 to -0.094 | 0.0029 | 65% |
| ImmuneScore | CLDN4 | fixed | 4 | 420 | -0.176 | -0.268 to -0.080 | 0.000356 | 0% |
| ImmuneScore | signature | DL random | 4 | 420 | -0.166 | -0.325 to +0.003 | 0.0548 | 60% |

Cross-study TACSTD2 vs CD8A partial on CLDN4 is -0.190 (fixed, p = 0.000108, I² = 0%, n sum = 420). The ImmuneScore partial on CLDN4 is -0.176 (p = 0.000356, I² = 0%). Partial on the 221-gene score, CD8A is -0.206 (p = 0.000151, I² = 13%) and ImmuneScore is -0.166 (p = 0.0548, I² = 60%).

## Mediation

X is TACSTD2. Y is CD8A or ImmuneScore. The primary mediator is the 221-gene score. Coefficients are standardized rank regressions, so the no-covariate total effect c equals the Spearman ρ.

| cohort | endpoint | mediator | covariates | n | a | b | c | c′ | a×b (95% CI) | boot p | proportion (95% CI) | call |
|---|---|---|---|---:|---:|---:|---:|---:|---|---:|---|---|
| OncoSG | CD8A | signature | none | 169 | +0.273 | -0.512 | -0.380 | -0.241 | -0.140 (-0.228 to -0.056) | 0.0004 | +0.367 (+0.160 to +0.637) | partial |
| OncoSG | CD8A | signature | composition | 169 | +0.157 | -0.370 | -0.269 | -0.211 | -0.058 (-0.131 to -0.005) | 0.0288 | +0.216 (+0.021 to +0.535) | partial |
| OncoSG | CD8A | signature | CLDN4 | 169 | +0.034 | -0.485 | -0.229 | -0.212 | -0.016 (-0.102 to +0.069) | 0.724 | +0.072 (-0.602 to +0.545) | indirect_not_supported |
| OncoSG | CD8A | CLDN4 | none | 169 | +0.505 | -0.301 | -0.380 | -0.229 | -0.152 (-0.254 to -0.067) | 0.0008 | +0.399 (+0.173 to +0.791) | partial |
| OncoSG | CD8A | CLDN4 | composition | 169 | +0.444 | -0.164 | -0.269 | -0.196 | -0.073 (-0.160 to -0.004) | 0.0396 | +0.271 (+0.014 to +0.761) | partial |
| OncoSG | ImmuneScore | signature | none | 169 | +0.273 | -0.550 | -0.387 | -0.237 | -0.150 (-0.243 to -0.064) | 0.0004 | +0.387 (+0.182 to +0.653) | partial |
| OncoSG | ImmuneScore | signature | composition | 169 | +0.157 | -0.430 | -0.280 | -0.212 | -0.068 (-0.145 to -0.004) | 0.0332 | +0.242 (+0.017 to +0.560) | partial |
| OncoSG | ImmuneScore | signature | CLDN4 | 169 | +0.034 | -0.523 | -0.227 | -0.209 | -0.018 (-0.106 to +0.078) | 0.693 | +0.078 (-0.611 to +0.555) | indirect_not_supported |
| OncoSG | ImmuneScore | CLDN4 | none | 169 | +0.505 | -0.317 | -0.387 | -0.227 | -0.160 (-0.266 to -0.076) | <0.0004 | +0.413 (+0.187 to +0.799) | partial |
| OncoSG | ImmuneScore | CLDN4 | composition | 169 | +0.444 | -0.188 | -0.280 | -0.196 | -0.083 (-0.169 to -0.018) | 0.0132 | +0.298 (+0.061 to +0.756) | partial |
| GSE273377 discovery | CD8A | signature | none | 103 | +0.259 | -0.419 | +0.039 | +0.148 | -0.109 (-0.211 to -0.023) | 0.0092 | -2.788 (-16.779 to +18.415) | total_not_inverse |
| GSE273377 discovery | CD8A | signature | composition | 103 | +0.295 | -0.352 | -0.026 | +0.078 | -0.104 (-0.199 to -0.028) | 0.0016 | +3.978 (-15.281 to +15.201) | total_not_inverse |
| GSE273377 discovery | CD8A | signature | CLDN4 | 103 | +0.087 | -0.458 | -0.030 | +0.010 | -0.040 (-0.166 to +0.055) | 0.419 | +1.355 (-4.985 to +5.199) | total_not_inverse |
| GSE273377 discovery | CD8A | CLDN4 | none | 103 | +0.638 | +0.107 | +0.039 | -0.030 | +0.069 (-0.083 to +0.250) | 0.379 | +1.759 (-13.139 to +12.757) | total_not_inverse |
| GSE273377 discovery | CD8A | CLDN4 | composition | 103 | +0.622 | +0.049 | -0.026 | -0.056 | +0.030 (-0.111 to +0.193) | 0.681 | -1.162 (-10.391 to +9.537) | total_not_inverse |
| GSE273377 discovery | ImmuneScore | signature | none | 103 | +0.259 | -0.446 | +0.096 | +0.212 | -0.116 (-0.225 to -0.024) | 0.0088 | -1.199 (-15.083 to +12.763) | total_not_inverse |
| GSE273377 discovery | ImmuneScore | signature | composition | 103 | +0.295 | -0.376 | +0.028 | +0.139 | -0.111 (-0.208 to -0.030) | 0.0012 | -3.962 (-19.228 to +19.726) | total_not_inverse |
| GSE273377 discovery | ImmuneScore | signature | CLDN4 | 103 | +0.087 | -0.490 | +0.018 | +0.061 | -0.043 (-0.179 to +0.057) | 0.42 | -2.328 (-7.600 to +6.977) | total_not_inverse |
| GSE273377 discovery | ImmuneScore | CLDN4 | none | 103 | +0.638 | +0.122 | +0.096 | +0.018 | +0.078 (-0.076 to +0.256) | 0.316 | +0.809 (-9.275 to +10.776) | total_not_inverse |
| GSE273377 discovery | ImmuneScore | CLDN4 | composition | 103 | +0.622 | +0.061 | +0.028 | -0.010 | +0.038 (-0.109 to +0.202) | 0.603 | +1.349 (-11.988 to +10.186) | total_not_inverse |
| GSE273377 validation | CD8A | signature | none | 60 | +0.181 | -0.353 | -0.423 | -0.359 | -0.064 (-0.188 to +0.025) | 0.189 | +0.151 (-0.093 to +0.502) | indirect_not_supported |
| GSE273377 validation | CD8A | signature | composition | 60 | +0.325 | -0.386 | -0.460 | -0.335 | -0.126 (-0.288 to -0.011) | 0.024 | +0.273 (+0.026 to +0.714) | partial |
| GSE273377 validation | CD8A | signature | CLDN4 | 60 | -0.162 | -0.362 | -0.315 | -0.374 | +0.059 (-0.052 to +0.189) | 0.299 | -0.186 (-2.450 to +1.334) | total_not_inverse |
| GSE273377 validation | CD8A | CLDN4 | none | 60 | +0.625 | -0.173 | -0.423 | -0.315 | -0.108 (-0.351 to +0.112) | 0.326 | +0.255 (-0.250 to +1.164) | indirect_not_supported |
| GSE273377 validation | CD8A | CLDN4 | composition | 60 | +0.705 | -0.132 | -0.460 | -0.367 | -0.093 (-0.360 to +0.174) | 0.472 | +0.202 (-0.365 to +1.004) | indirect_not_supported |
| GSE273377 validation | ImmuneScore | signature | none | 60 | +0.181 | -0.383 | -0.330 | -0.261 | -0.069 (-0.198 to +0.026) | 0.18 | +0.210 (-0.151 to +0.790) | indirect_not_supported |
| GSE273377 validation | ImmuneScore | signature | composition | 60 | +0.325 | -0.419 | -0.370 | -0.233 | -0.136 (-0.288 to -0.015) | 0.0172 | +0.369 (+0.043 to +1.063) | indirect_direct_crosses_zero |
| GSE273377 validation | ImmuneScore | signature | CLDN4 | 60 | -0.162 | -0.474 | -0.333 | -0.409 | +0.077 (-0.060 to +0.237) | 0.3 | -0.231 (-2.142 to +0.566) | total_not_inverse |
| GSE273377 validation | ImmuneScore | CLDN4 | none | 60 | +0.625 | +0.004 | -0.330 | -0.333 | +0.003 (-0.231 to +0.237) | 0.969 | -0.008 (-1.067 to +1.038) | indirect_not_supported |
| GSE273377 validation | ImmuneScore | CLDN4 | composition | 60 | +0.705 | +0.075 | -0.370 | -0.423 | +0.053 (-0.221 to +0.314) | 0.721 | -0.144 (-1.313 to +0.755) | indirect_not_supported |
| GSE282774 | CD8A | signature | none | 58 | +0.424 | -0.152 | -0.206 | -0.142 | -0.064 (-0.230 to +0.060) | 0.313 | +0.312 (-2.510 to +3.261) | total_not_inverse |
| GSE282774 | CD8A | signature | composition | 58 | +0.421 | -0.142 | -0.205 | -0.145 | -0.060 (-0.231 to +0.064) | 0.36 | +0.292 (-2.565 to +3.313) | total_not_inverse |
| GSE282774 | CD8A | signature | CLDN4 | 58 | +0.314 | -0.166 | -0.272 | -0.220 | -0.052 (-0.219 to +0.035) | 0.309 | +0.192 (-1.045 to +1.849) | total_not_inverse |
| GSE282774 | CD8A | CLDN4 | none | 58 | +0.659 | +0.100 | -0.206 | -0.272 | +0.066 (-0.203 to +0.292) | 0.634 | -0.321 (-5.181 to +3.948) | total_not_inverse |
| GSE282774 | CD8A | CLDN4 | composition | 58 | +0.662 | +0.085 | -0.205 | -0.261 | +0.057 (-0.218 to +0.277) | 0.66 | -0.276 (-4.603 to +3.734) | total_not_inverse |
| GSE282774 | ImmuneScore | signature | none | 58 | +0.424 | -0.269 | -0.179 | -0.065 | -0.114 (-0.292 to +0.011) | 0.0736 | +0.637 (-3.832 to +5.985) | total_not_inverse |
| GSE282774 | ImmuneScore | signature | composition | 58 | +0.421 | -0.249 | -0.177 | -0.072 | -0.105 (-0.281 to +0.017) | 0.0984 | +0.594 (-3.204 to +6.177) | total_not_inverse |
| GSE282774 | ImmuneScore | signature | CLDN4 | 58 | +0.314 | -0.279 | -0.209 | -0.121 | -0.088 (-0.263 to +0.018) | 0.142 | +0.420 (-2.632 to +4.249) | total_not_inverse |
| GSE282774 | ImmuneScore | CLDN4 | none | 58 | +0.659 | +0.045 | -0.179 | -0.209 | +0.030 (-0.236 to +0.245) | 0.812 | -0.165 (-4.698 to +6.387) | total_not_inverse |
| GSE282774 | ImmuneScore | CLDN4 | composition | 58 | +0.662 | +0.011 | -0.177 | -0.184 | +0.007 (-0.257 to +0.230) | 0.991 | -0.039 (-4.693 to +5.325) | total_not_inverse |
| GSE233774 tumor | CD8A | signature | none | 30 | +0.345 | -0.205 | -0.437 | -0.366 | -0.071 (-0.226 to +0.086) | 0.355 | +0.162 (-0.274 to +1.038) | indirect_not_supported |
| GSE233774 tumor | CD8A | signature | composition | 30 | +0.207 | -0.040 | -0.231 | -0.223 | -0.008 (-0.128 to +0.113) | 0.907 | +0.036 (-1.191 to +1.281) | total_not_inverse |
| GSE233774 tumor | CD8A | signature | CLDN4 | 30 | +0.125 | -0.219 | -0.428 | -0.401 | -0.027 (-0.159 to +0.164) | 0.831 | +0.064 (-1.086 to +1.227) | total_not_inverse |
| GSE233774 tumor | CD8A | CLDN4 | none | 30 | +0.647 | -0.013 | -0.437 | -0.428 | -0.008 (-0.556 to +0.237) | 0.908 | +0.019 (-0.876 to +2.809) | indirect_not_supported |
| GSE233774 tumor | CD8A | CLDN4 | composition | 30 | +0.487 | +0.371 | -0.231 | -0.411 | +0.181 (-0.194 to +0.352) | 0.257 | -0.782 (-6.687 to +7.055) | total_not_inverse |
| GSE233774 tumor | ImmuneScore | signature | none | 30 | +0.345 | -0.328 | -0.439 | -0.326 | -0.113 (-0.279 to +0.038) | 0.159 | +0.258 (-0.142 to +1.046) | indirect_not_supported |
| GSE233774 tumor | ImmuneScore | signature | composition | 30 | +0.207 | -0.177 | -0.231 | -0.195 | -0.037 (-0.199 to +0.055) | 0.508 | +0.158 (-0.960 to +1.513) | total_not_inverse |
| GSE233774 tumor | ImmuneScore | signature | CLDN4 | 30 | +0.125 | -0.335 | -0.387 | -0.345 | -0.042 (-0.238 to +0.156) | 0.656 | +0.108 (-1.325 to +1.973) | total_not_inverse |
| GSE233774 tumor | ImmuneScore | CLDN4 | none | 30 | +0.647 | -0.080 | -0.439 | -0.387 | -0.052 (-0.535 to +0.193) | 0.694 | +0.118 (-0.635 to +2.200) | indirect_not_supported |
| GSE233774 tumor | ImmuneScore | CLDN4 | composition | 30 | +0.487 | +0.287 | -0.231 | -0.371 | +0.140 (-0.236 to +0.370) | 0.381 | -0.603 (-7.283 to +7.954) | total_not_inverse |

Primary OncoSG paths, no composition covariate:

- CD8A via 221-gene score: a = +0.273 (CI +0.117 to +0.414), b = -0.512 (CI -0.634 to -0.379), c = -0.380, c′ = -0.241 (CI -0.382 to -0.101), a×b = -0.140 (CI -0.228 to -0.056, boot p = 0.0004), proportion = +0.367 (CI +0.160 to +0.637). Call: partial.
- CD8A via CLDN4: a = +0.505 (CI +0.370 to +0.623), b = -0.301 (CI -0.457 to -0.146), c = -0.380, c′ = -0.229 (CI -0.386 to -0.056), a×b = -0.152 (CI -0.254 to -0.067, boot p = 0.0008), proportion = +0.399 (CI +0.173 to +0.791). Call: partial.
- ImmuneScore via 221-gene score: a = +0.273 (CI +0.123 to +0.425), b = -0.550 (CI -0.668 to -0.415), c = -0.387, c′ = -0.237 (CI -0.370 to -0.100), a×b = -0.150 (CI -0.243 to -0.064, boot p = 0.0004), proportion = +0.387 (CI +0.182 to +0.653). Call: partial.
- ImmuneScore via CLDN4: a = +0.505 (CI +0.368 to +0.623), b = -0.317 (CI -0.475 to -0.166), c = -0.387, c′ = -0.227 (CI -0.387 to -0.057), a×b = -0.160 (CI -0.266 to -0.076, boot p < 0.0004), proportion = +0.413 (CI +0.187 to +0.799). Call: partial.

Cross-study indirect effect a×b for the 221-gene score, no composition covariate. GSE273377 strata are combined first. Variance is the bootstrap variance.

| endpoint | model | k | a×b | 95% CI | p | I² |
|---|---|---:|---:|---|---:|---:|
| CD8A | fixed | 4 | -0.100 | -0.149 to -0.052 | 5.00e-05 | 0% |
| ImmuneScore | fixed | 4 | -0.117 | -0.167 to -0.066 | 6.56e-06 | 0% |

The same indirect effect with the composition covariate in every equation:

| endpoint | model | k | a×b | 95% CI | p | I² |
|---|---|---:|---:|---|---:|---:|
| CD8A | fixed | 4 | -0.067 | -0.109 to -0.025 | 0.00165 | 0% |
| ImmuneScore | fixed | 4 | -0.086 | -0.131 to -0.040 | 0.000216 | 0% |

Signature as mediator with CLDN4 entered as a covariate (does the score carry an indirect path beyond the single gene):

| cohort | endpoint | a×b (95% CI) | boot p | call |
|---|---|---|---:|---|
| OncoSG | CD8A | -0.016 (-0.102 to +0.069) | 0.724 | indirect_not_supported |
| OncoSG | ImmuneScore | -0.018 (-0.106 to +0.078) | 0.693 | indirect_not_supported |
| GSE273377 discovery | CD8A | -0.040 (-0.166 to +0.055) | 0.419 | total_not_inverse |
| GSE273377 discovery | ImmuneScore | -0.043 (-0.179 to +0.057) | 0.42 | total_not_inverse |
| GSE273377 validation | CD8A | +0.059 (-0.052 to +0.189) | 0.299 | total_not_inverse |
| GSE273377 validation | ImmuneScore | +0.077 (-0.060 to +0.237) | 0.3 | total_not_inverse |
| GSE282774 | CD8A | -0.052 (-0.219 to +0.035) | 0.309 | total_not_inverse |
| GSE282774 | ImmuneScore | -0.088 (-0.263 to +0.018) | 0.142 | total_not_inverse |
| GSE233774 tumor | CD8A | -0.027 (-0.159 to +0.164) | 0.831 | total_not_inverse |
| GSE233774 tumor | ImmuneScore | -0.042 (-0.238 to +0.156) | 0.656 | total_not_inverse |

## Answer

Partialling CLDN4 leaves an inverse TACSTD2–CD8A association on OncoSG (-0.217, p = 0.00471, n = 169) and in the four-study combination (-0.190, fixed, p = 0.000108, I² = 0%). The ImmuneScore partial on CLDN4 is -0.217 on OncoSG (p = 0.00463) and -0.176 across studies (p = 0.000356, I² = 0%). Partialling the 221-gene score leaves OncoSG TACSTD2–CD8A at -0.284 (p = 0.000193) and the cross-study CD8A partial at -0.206 (p = 0.000151, I² = 13%).

Under the pre-specified call, the 221-gene score is a partial statistical mediator in 2 of 10 primary cohort×endpoint rows (both are OncoSG). A direct interval includes zero in 0 rows. A proportion outside 0–1 occurs in 0 rows. 8 primary rows fail the rule. On OncoSG the indirect path remains after published PURITY (CD8A a×b = -0.058, CI -0.131 to -0.005).

GSE273377 discovery is a different pattern. TACSTD2 vs CD8A Spearman is +0.039 (p = 0.696, n = 103), and vs ImmuneScore it is +0.096 (p = 0.332). The indirect path in that stratum is inverse (CD8A a×b = -0.109, CI -0.211 to -0.023) while the direct coefficient is positive (c′ = +0.148). That is an indirect path beside a total association that is not inverse. It is not evidence that the score explains an inverse TACSTD2–immune correlation in that stratum.

Averaging every indirect path, including that discovery stratum, gives a cross-study a×b of -0.100 for CD8A (fixed, 95% CI -0.149 to -0.052, p = 5.00e-05, I² = 0%) and -0.117 for ImmuneScore (95% CI -0.167 to -0.066, p = 6.56e-06, I² = 0%). With the composition covariate in every equation, the CD8A indirect effect is -0.067 (95% CI -0.109 to -0.025, p = 0.00165). Restricting to cohort-endpoints whose total-effect interval lies entirely below zero (OncoSG, GSE273377 validation, GSE233774 tumor for CD8A) gives a×b = -0.104 (k = 3, 95% CI -0.165 to -0.042, p = 0.000918, I² = 0%). The ImmuneScore restriction (OncoSG, GSE273377 validation, GSE233774 tumor) is -0.118 (k = 3, 95% CI -0.182 to -0.054, p = 0.000308, I² = 0%). Inside that restricted set, indirect intervals entirely below zero are: CD8A in OncoSG; ImmuneScore in OncoSG.

With CLDN4 as a covariate, the signature indirect path meets the rule in 0 of 10 cohort×endpoint rows. On OncoSG that indirect effect is -0.016 for CD8A (CI -0.102 to +0.069, call indirect_not_supported) and -0.018 for ImmuneScore (CI -0.106 to +0.078, call indirect_not_supported).

## What this does not say

The model has no exposure–mediator interaction and no unmeasured-confounding correction beyond the stated covariates. Bulk RNA from a surgical resection does not identify a cellular path from TACSTD2 to CD8A. The signature was built to track CLDN4 after keratin and purity adjustment, so a positive a path from TACSTD2 into the score is expected from the known TACSTD2–CLDN4 correlation. Agreement of the score with CD8A beyond CLDN4 is the extra quantity, and it is the CLDN4-covariate mediation block. TCGA was the filter that chose the 221 genes and is not in this mediation.

## Reproduce

```bash
python3 methods/tacstd2_cldn4_mediation/analyze.py
```

