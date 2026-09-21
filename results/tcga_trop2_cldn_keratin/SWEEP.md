# Immune-limb specification sweep

Numbers are written by `scripts/tcga_trop2_cldn_keratin/sweep_immune.py`.

The question for this sweep is which pre-set choice of immune gene set, keratin covariates, purity covariate, residualization, histology arm, or CLDN4/TACSTD2/CLDN7 quantile puts the barrier genes with a lower T, CD8, or cytotoxic score. Thesis-aligned means partial ρ < 0, or, for quantiles, Cliff's delta < 0 (Q4 lower than Q1).

The locked CLDN4-versus-KRT8 surface-gene ranking is not recomputed and is not overwritten.

## Grid

Continuous grid, eight cohorts (LUAD, LUSC, BRCA, CESC, KIRC, STAD, BLCA, PAAD), three predictors, 8 immune scores, 4 keratin sets, three purity settings, two residualizations. Completed rows: 4608. Thesis-aligned rows: 2849.

ESTIMATE purity is the Yoshihara cosine transform of the MD Anderson RNAseqV2 ESTIMATE score and covers all eight cohorts. ABSOLUTE purity is Aran et al. 2015. That table has no STAD or PAAD rows, and CESC has no ABSOLUTE call, so ABSOLUTE-adjusted tests drop those cohorts. A panel is eligible for the best-panel rank only when all eight cohorts still contribute to every predictor.

Immune scores: CD8A; CD8A+CD8B; CD3D/E/G; CD3D/E/G+CD2+CD247+LCK; Rooney CYT (GZMA+PRF1); cytotoxic four (GZMA, GZMB, PRF1, NKG7); cytotoxic eight (those four plus GNLY, GZMK, GZMH, CTSW); CD8 effector (CD8A, CD8B, GZMB, PRF1, NKG7, IFNG).

Keratin sets: none; KRT8+KRT18+KRT19; KRT5+KRT6A+KRT6B+KRT14; the union plus KRT17.

Spearman residualization ranks every variable, then correlates OLS residuals. Pearson residualization correlates OLS residuals of the log2(TPM+1) values.

LUAD and LUSC open clinical matrices are mostly histology NOS, so they are not split. Arms with n≥40: BRCA ductal vs lobular, CESC squamous vs adenocarcinoma, CESC keratinizing vs non-keratinizing squamous, STAD intestinal vs diffuse, BLCA papillary vs non-papillary.

## Best shared panel

A shared panel is one outcome, keratin set, purity covariate, and residualization, applied to all three predictors and all eight cohorts. Panels are ranked by the worst of the three random-effects meta p-values, and only panels with all three pooled ρ values negative are eligible. `q_worst` is Benjamini-Hochberg across every shared panel's worst p, including panels that are not thesis-aligned.

Shared panels: 192. Panels with all three pooled ρ < 0: 128. Eligible eight-cohort panels among those: 64.

Best shared panel: outcome `T_CD3`, keratin `squamous_KRT5_6`, purity `none`, residualization `pearson`.

TACSTD2 pooled ρ=-0.080 (p=9.23e-05, negative cohorts 6/8). CLDN4 pooled ρ=-0.127 (p=3.38e-09, 8/8). CLDN7 pooled ρ=-0.160 (p=3.28e-05, 7/8). Worst p=9.23e-05. Sweep q on that worst p=0.0122.

Per-cohort partial ρ for this panel:

| cohort | n | TACSTD2 | CLDN4 | CLDN7 |
| --- | --- | --- | --- | --- |
| LUAD | 516 | -0.044 | -0.118 | -0.302 |
| LUSC | 501 | -0.089 | -0.133 | -0.149 |
| BRCA | 1095 | -0.090 | -0.049 | -0.144 |
| CESC | 304 | 0.003 | -0.088 | 0.025 |
| KIRC | 533 | -0.097 | -0.124 | -0.049 |
| STAD | 412 | -0.094 | -0.169 | -0.116 |
| BLCA | 406 | -0.192 | -0.184 | -0.193 |
| PAAD | 178 | 0.040 | -0.249 | -0.365 |

Next shared panels with all three pooled ρ < 0, ordered by worst meta p:

| outcome | keratin | purity | method | rho_TACSTD2 | p_TACSTD2 | rho_CLDN4 | p_CLDN4 | rho_CLDN7 | p_CLDN7 | worst_p | q_worst |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| T_CD3 | squamous_KRT5_6 | none | pearson | -0.080 | 9.23e-05 | -0.127 | 3.38e-09 | -0.160 | 3.28e-05 | 9.23e-05 | 0.0122 |
| CD8A | squamous_KRT5_6 | none | pearson | -0.100 | 0.0001 | -0.129 | 2.17e-05 | -0.149 | 8.26e-05 | 0.0001 | 0.0122 |
| CD8_AB | squamous_KRT5_6 | none | spearman | -0.117 | 0.0002 | -0.112 | 8.51e-05 | -0.142 | 0.0002 | 0.0002 | 0.0122 |
| T_CD3 | squamous_KRT5_6 | none | spearman | -0.102 | 0.0003 | -0.131 | 1.52e-08 | -0.165 | 0.0001 | 0.0003 | 0.0122 |
| CD8_AB | squamous_KRT5_6 | none | pearson | -0.108 | 0.0003 | -0.128 | 6.06e-05 | -0.145 | 0.0001 | 0.0003 | 0.0133 |
| T_broad | squamous_KRT5_6 | none | spearman | -0.097 | 0.0004 | -0.123 | 6.91e-08 | -0.158 | 0.0002 | 0.0004 | 0.0135 |
| cytotoxic_8 | squamous_KRT5_6 | none | spearman | -0.107 | 0.0001 | -0.136 | 2.03e-05 | -0.142 | 0.0005 | 0.0005 | 0.0135 |
| cytotoxic_8 | squamous_KRT5_6 | none | pearson | -0.100 | 2.21e-06 | -0.135 | 4.21e-06 | -0.138 | 0.0006 | 0.0006 | 0.0137 |
| T_broad | squamous_KRT5_6 | none | pearson | -0.078 | 0.0008 | -0.118 | 1.64e-08 | -0.153 | 6.23e-05 | 0.0008 | 0.0146 |
| CD8A | squamous_KRT5_6 | none | spearman | -0.112 | 0.0008 | -0.117 | 2.53e-05 | -0.147 | 0.0002 | 0.0008 | 0.0146 |
| CD8_effector | squamous_KRT5_6 | none | pearson | -0.119 | 3.14e-06 | -0.139 | 3.69e-05 | -0.144 | 0.0008 | 0.0008 | 0.0146 |
| CD8_effector | squamous_KRT5_6 | none | spearman | -0.124 | 4.67e-05 | -0.130 | 0.0001 | -0.141 | 0.0013 | 0.0013 | 0.0203 |

The original specification (CD8A+CD8B, KRT8+KRT18+KRT19, no purity, Spearman) remains in the grid: TACSTD2 ρ=-0.069 (p=0.0427), CLDN4 ρ=-0.081 (p=0.0062), CLDN7 ρ=-0.109 (p=0.0009). Its worst p is 0.0427. The KRT5/6 panels above are the stronger shared specifications.

ESTIMATE-adjusted panels that keep all three pooled ρ values negative in all eight cohorts: 0.

## Strongest single cohort tests

These are the smallest two-sided p-values among continuous-grid tests with ρ < 0. The grid contains 4608 tests. `q_grid` is Benjamini-Hochberg across all 4608 continuous tests, not only the negative ones.

| cohort | predictor | outcome | keratin | purity | method | n | rho | p | q_grid |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LUAD | CLDN7 | T_CD3 | none | none | spearman | 516 | -0.355 | 9.96e-17 | 4.02e-13 |
| LUAD | CLDN7 | T_broad | none | none | spearman | 516 | -0.349 | 3.06e-16 | 6.17e-13 |
| LUAD | CLDN7 | T_CD3 | squamous_KRT5_6 | none | spearman | 516 | -0.341 | 1.92e-15 | 2.58e-12 |
| LUAD | CLDN7 | T_broad | squamous_KRT5_6 | none | spearman | 516 | -0.334 | 7.81e-15 | 7.87e-12 |
| LUAD | CLDN7 | CYT_Rooney | squamous_KRT5_6 | none | spearman | 516 | -0.329 | 2.25e-14 | 1.53e-11 |
| LUAD | CLDN7 | CYT_Rooney | none | none | spearman | 516 | -0.328 | 2.28e-14 | 1.53e-11 |
| LUAD | CLDN7 | T_CD3 | simple_KRT8_18_19 | none | spearman | 516 | -0.317 | 2.10e-13 | 9.40e-11 |
| BLCA | CLDN4 | CYT_Rooney | none | none | pearson | 406 | -0.353 | 2.50e-13 | 1.01e-10 |
| LUAD | CLDN7 | CYT_Rooney | squamous_KRT5_6 | none | pearson | 516 | -0.311 | 6.09e-13 | 1.72e-10 |
| LUAD | CLDN7 | CD8_effector | squamous_KRT5_6 | none | spearman | 516 | -0.311 | 6.40e-13 | 1.72e-10 |
| LUAD | CLDN7 | CYT_Rooney | simple_KRT8_18_19 | none | spearman | 516 | -0.309 | 8.41e-13 | 2.12e-10 |
| BRCA | CLDN7 | CD8A | none | none | spearman | 1095 | -0.213 | 9.65e-13 | 2.16e-10 |

## Histology arms

Histology-grid rows: 1800. The same BH is computed inside this family only.

| stratum | predictor | outcome | keratin | purity | method | n | rho | p | q_grid |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BLCA_nonpapillary | CLDN4 | CYT_Rooney | none | none | pearson | 269 | -0.362 | 9.86e-10 | 1.78e-06 |
| BRCA_ductal | CLDN7 | T_CD3 | none | none | spearman | 783 | -0.206 | 5.79e-09 | 5.21e-06 |
| BLCA_nonpapillary | CLDN4 | cytotoxic_4 | none | none | pearson | 269 | -0.337 | 1.42e-08 | 8.32e-06 |
| BRCA_ductal | CLDN7 | CYT_Rooney | none | none | spearman | 783 | -0.199 | 2.03e-08 | 8.32e-06 |
| BRCA_ductal | CLDN7 | cytotoxic_8 | none | none | spearman | 783 | -0.198 | 2.49e-08 | 8.32e-06 |
| BRCA_ductal | CLDN7 | T_CD3 | none | none | pearson | 783 | -0.197 | 2.77e-08 | 8.32e-06 |
| BRCA_ductal | CLDN7 | cytotoxic_8 | none | none | pearson | 783 | -0.191 | 7.51e-08 | 1.35e-05 |
| BLCA_nonpapillary | TACSTD2 | CYT_Rooney | squamous_KRT5_6 | none | pearson | 269 | -0.321 | 9.49e-08 | 1.52e-05 |
| BRCA_ductal | CLDN7 | cytotoxic_4 | none | none | spearman | 783 | -0.189 | 1.01e-07 | 1.52e-05 |
| BRCA_ductal | CLDN7 | T_CD3 | squamous_KRT5_6 | none | spearman | 783 | -0.187 | 1.56e-07 | 1.92e-05 |
| BRCA_ductal | CLDN7 | cytotoxic_4 | none | none | pearson | 783 | -0.185 | 1.82e-07 | 1.92e-05 |
| BLCA_nonpapillary | CLDN4 | CYT_Rooney | none | none | spearman | 269 | -0.311 | 1.84e-07 | 1.92e-05 |

## Quantiles

Q4 versus Q1 of the predictor. The outcome is either the raw score or the residual of that score after the named adjustment. Cliff's delta < 0 means the upper quartile has the lower immune score. Stouffer combines the eight cohorts. `q` is Benjamini-Hochberg across quantile specifications.

| predictor | outcome | adjustment | n_cohorts | n_negative | median_delta | stouffer_z | stouffer_p | q |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CLDN7 | T_CD3 | none | 8 | 8 | -0.250 | -9.099 | 9.14e-20 | 1.54e-17 |
| CLDN7 | T_CD3 | squamous_pearson | 8 | 8 | -0.221 | -8.788 | 1.52e-18 | 1.28e-16 |
| CLDN7 | T_CD3 | squamous_spearman | 8 | 8 | -0.212 | -8.368 | 5.87e-17 | 3.28e-15 |
| CLDN7 | CD8A | squamous_pearson | 8 | 7 | -0.180 | -8.084 | 6.29e-16 | 2.64e-14 |
| CLDN7 | CD8A | none | 8 | 8 | -0.184 | -8.037 | 9.24e-16 | 3.10e-14 |
| CLDN4 | CYT_Rooney | none | 8 | 8 | -0.221 | -7.994 | 1.30e-15 | 3.65e-14 |
| CLDN4 | CYT_Rooney | squamous_pearson | 8 | 8 | -0.228 | -7.752 | 9.05e-15 | 1.83e-13 |
| CLDN7 | CYT_Rooney | none | 8 | 7 | -0.167 | -7.747 | 9.39e-15 | 1.83e-13 |
| CLDN7 | CYT_Rooney | squamous_pearson | 8 | 7 | -0.170 | -7.742 | 9.79e-15 | 1.83e-13 |
| CLDN7 | CD8_AB | squamous_pearson | 8 | 7 | -0.184 | -7.702 | 1.34e-14 | 2.25e-13 |
| CLDN4 | cytotoxic_8 | none | 8 | 8 | -0.230 | -7.670 | 1.73e-14 | 2.64e-13 |
| CLDN7 | CD8A | squamous_spearman | 8 | 7 | -0.171 | -7.631 | 2.34e-14 | 3.27e-13 |

Per-cohort Cliff's delta for the strongest keratin-adjusted quantile specification of each predictor (negative = Q4 lower than Q1):

Columns are TACSTD2 vs CD8 effector after a KRT5/6 Spearman residual, CLDN4 vs Rooney CYT after a KRT5/6 Pearson residual, and CLDN7 vs CD3 after a KRT5/6 Pearson residual.

| cohort | TACSTD2 | CLDN4 | CLDN7 |
| --- | --- | --- | --- |
| LUAD | -0.240 | -0.232 | -0.498 |
| LUSC | -0.177 | -0.227 | -0.199 |
| BRCA | -0.199 | -0.055 | -0.243 |
| CESC | 0.088 | -0.291 | -0.016 |
| KIRC | -0.022 | -0.192 | -0.014 |
| STAD | -0.033 | -0.186 | -0.200 |
| BLCA | -0.323 | -0.228 | -0.273 |
| PAAD | -0.309 | -0.435 | -0.504 |

## NHEJ, STING, and IFN modules versus the barrier genes

IFN-negative is the immune-cold direction. NHEJ and STING are reported at the sign the data give. The block below is partial Spearman of CLDN4 versus each module after KRT8+KRT18+KRT19, which matches the original keratin model. The full module grid is `tables/module_partial.tsv`.

| cohort | outcome | n | rho | p |
| --- | --- | --- | --- | --- |
| LUAD | IFN | 516 | -0.091 | 0.0385 |
| LUSC | IFN | 501 | -0.036 | 0.4264 |
| BRCA | IFN | 1095 | 0.036 | 0.2364 |
| CESC | IFN | 304 | -0.165 | 0.0041 |
| KIRC | IFN | 533 | -0.021 | 0.6249 |
| STAD | IFN | 412 | -0.095 | 0.0537 |
| BLCA | IFN | 406 | -0.117 | 0.0189 |
| PAAD | IFN | 178 | -0.062 | 0.4149 |
| LUAD | NHEJ | 516 | 0.106 | 0.0162 |
| LUSC | NHEJ | 501 | -0.046 | 0.3075 |
| BRCA | NHEJ | 1095 | -0.024 | 0.4323 |
| CESC | NHEJ | 304 | 0.124 | 0.0318 |
| KIRC | NHEJ | 533 | 0.151 | 0.0005 |
| STAD | NHEJ | 412 | 0.177 | 0.0003 |
| BLCA | NHEJ | 406 | 0.024 | 0.6333 |
| PAAD | NHEJ | 178 | 0.047 | 0.5371 |
| LUAD | STING | 516 | 0.303 | 2.53e-12 |
| LUSC | STING | 501 | 0.034 | 0.4426 |
| BRCA | STING | 1095 | 0.181 | 1.79e-09 |
| CESC | STING | 304 | 0.159 | 0.0058 |
| KIRC | STING | 533 | -0.047 | 0.2776 |
| STAD | STING | 412 | 0.139 | 0.0047 |
| BLCA | STING | 406 | -0.045 | 0.3638 |
| PAAD | STING | 178 | 0.267 | 0.0004 |

Random-effects pool, CLDN4, same keratin model, Spearman:

- IFN: ρ=-0.062 (95% CI -0.111 to -0.012), p=0.0147, I²=56.6%, cohorts=8.
- STING: ρ=0.123 (95% CI 0.031 to 0.213), p=0.0088, I²=87.7%, cohorts=8.
- NHEJ: ρ=0.068 (95% CI 0.005 to 0.130), p=0.0342, I²=73.2%, cohorts=8.

## Reading

Use the best shared panel as the sweep's answer for a specification that moves TACSTD2, CLDN4, and CLDN7 together. Use `q_worst` and `q_grid` when citing a p-value from this search. Cohort-level minima and histology arms describe where the association is largest. They are part of the same search.

