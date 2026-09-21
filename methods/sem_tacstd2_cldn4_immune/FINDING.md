# Which graph fits: TACSTD2, CLDN4, and immune

Three linear Gaussian graphs were fit to the same three observed variables. Each restricted graph has 5 free covariance parameters, so AIC and BIC rank them by likelihood. The saturated graph adds the direct path (6 parameters) and is the reference for a likelihood-ratio test.

M1 is TACSTD2 → CLDN4 → immune, with no direct TACSTD2 → immune arrow. M2 is CLDN4 → TACSTD2 → immune, with no direct CLDN4 → immune arrow. M3 is independent: TACSTD2 and CLDN4 are uncorrelated, and each has its own arrow into the immune readout.

A lower BIC is a better description of this covariance. It is not a knockdown, not an instrument, and not by itself a causal order. Concordant-4 and TCGA are not entered into one likelihood.

## Concordant-4

Locked units only: GSE123902, GSE131907, GSE205335, and GSE189357 (n = 65). Primary variables are the malignant-cell detection percentages `pct_TACSTD2` and `pct_CLDN4`, and the unit T/NK fraction `frac_tnk`. Each variable is z-scored inside its dataset (maximum-likelihood standard deviation), then the four datasets are stacked. Percent positive is the locked CLDN4 measurement: the same table reproduces the locked DerSimonian–Laird Spearman.

On the concordant-4 primary specification (percent positive, n = 65), using one SEM on the within-dataset z-scores, the lowest BIC is **TACSTD2 → CLDN4 → immune** (BIC 539.4; ΔBIC 14.66 versus CLDN4 → TACSTD2 → immune (greater than 10); ΔBIC 17.95 versus independent (TACSTD2 ⊥ CLDN4, both → immune)). AIC selects the same graph. The saturated model does not improve BIC (ΔBIC 3.39; likelihood-ratio versus the winner χ² = 0.78 on 1 df, p = 0.3764). On the sum of the four dataset-specific BICs, with a separate 5-parameter model in each dataset, the lowest BIC is **TACSTD2 → CLDN4 → immune** (BIC 538.0; ΔBIC 8.75 versus CLDN4 → TACSTD2 → immune (a large gap); ΔBIC 46.80 versus independent (TACSTD2 ⊥ CLDN4, both → immune)). AIC selects the same graph. The saturated model does not improve BIC (ΔBIC 8.04; likelihood-ratio versus the winner χ² = 2.86 on 4 df, p = 0.5822).

## TCGA

UCSC Xena GDC STAR log2(TPM+1), primary solid tumor only, replicate aliquots averaged. Primary cohorts are the locked keratin funnel: LUAD, BRCA, CESC, KIRC, STAD, BLCA, PAAD (n = 3,444). LUSC is estimated and is not in the sum. The immune readout is the CD8 score, the mean of CD8A and CD8B. Inside each cohort, TACSTD2, CLDN4, and the CD8 score are residualized on KRT8, KRT18, and KRT19, then z-scored. The fit that decides the TCGA comparison is the sum of those cohort-specific likelihoods. A pooled z-score fit is reported beside it.

On the TCGA primary specification (CD8 score, KRT8/18/19 residuals, seven cohorts, n = 3,444), summing cohort-specific BICs, the lowest BIC is **TACSTD2 → CLDN4 → immune** (BIC 29129.8; ΔBIC 11.60 versus CLDN4 → TACSTD2 → immune (greater than 10); ΔBIC 321.96 versus independent (TACSTD2 ⊥ CLDN4, both → immune)). AIC selects the same graph. The saturated model does not improve BIC (ΔBIC 11.08; likelihood-ratio versus the winner χ² = 31.37 on 7 df, p = 5.31e-05). AIC still prefers the saturated model (AIC 28970.2 versus 28987.6 for the BIC winner). The extra direct-path parameters improve the likelihood enough for AIC and not enough for BIC. On the same residuals stacked after within-cohort z-scoring, as one pooled SEM, the lowest BIC is **TACSTD2 → CLDN4 → immune** (BIC 29012.0; ΔBIC 14.25 versus CLDN4 → TACSTD2 → immune (greater than 10); ΔBIC 322.82 versus independent (TACSTD2 ⊥ CLDN4, both → immune)). AIC selects the same graph. The saturated model does not improve BIC (ΔBIC 5.19; likelihood-ratio versus the winner χ² = 2.95 on 1 df, p = 0.0858).

## Mediation percentage

The percentage is the product-of-coefficients indirect path divided by the total effect, with the direct path retained. On z-scored variables this is a standardized path decomposition. The bootstrap resamples patients inside each dataset or cohort, repeats the z-scoring (and, for TCGA, the keratin residualization), and refits. There are 2000 replicates and the seed is 20260921. The interval is the 2.5 and 97.5 percentiles.

Concordant-4 primary: TACSTD2 → CLDN4 → immune is 173.8% (bootstrap 95% CI -1152.5 to 1532.0; standardized indirect -0.263, direct 0.112, total -0.151); direct and indirect paths have opposite signs (inconsistent mediation), so the percentage is not a share of a one-direction effect; 6.0% of bootstrap draws have a standardized total effect below 0.02 in absolute value, which stretches the ratio. CLDN4 → TACSTD2 → immune is -11.9% (bootstrap 95% CI -62.9 to 18.8; standardized indirect 0.056, direct -0.525, total -0.469); direct and indirect paths have opposite signs (inconsistent mediation), so the percentage is not a share of a one-direction effect.

TCGA primary, pooled within-cohort z-scores: TACSTD2 → CLDN4 → immune is 42.1% (bootstrap 95% CI 19.5 to 114.9; standardized indirect -0.022, direct -0.031, total -0.053); 3.0% of bootstrap draws have a standardized total effect below 0.02 in absolute value, which stretches the ratio. CLDN4 → TACSTD2 → immune is 11.1% (bootstrap 95% CI -1.1 to 28.2; standardized indirect -0.009, direct -0.074, total -0.083).

## Reading the percentage and the cohort split

On concordant-4, TACSTD2 → CLDN4 → T/NK splits into a standardized indirect path of -0.263 (bootstrap 95% CI -0.439 to -0.123) and a direct path of 0.112 (bootstrap 95% CI -0.135 to 0.339). The total effect is -0.151 (bootstrap 95% CI -0.394 to 0.136). The indirect path is negative and the direct path is positive, so they oppose each other. The ratio of the indirect path to that net total is 173.8%, and the bootstrap interval on the ratio is -1153 to 1532. The total-effect interval includes zero. The indirect-path interval does not. The ratio is not an estimate of a mediation share.

The reverse order, CLDN4 → TACSTD2 → T/NK, has indirect 0.056 (bootstrap 95% CI -0.086 to 0.192), direct -0.525 (bootstrap 95% CI -0.726 to -0.264), and total -0.469 (bootstrap 95% CI -0.665 to -0.236). The ratio is -11.9% (bootstrap 95% CI -62.9 to 18.8). The CLDN4 association with T/NK stays on the direct path. It is not carried by TACSTD2.

On the pooled TCGA CD8 residuals both paths in TACSTD2 → CLDN4 → CD8 are negative: indirect -0.022 (bootstrap 95% CI -0.034 to -0.011), direct -0.031 (bootstrap 95% CI -0.066 to 0.003), total -0.053 (bootstrap 95% CI -0.087 to -0.019). The ratio is 42.1% (bootstrap 95% CI 19.5 to 114.9). The interval extends above 100, so the split between the indirect path and the direct path is not pinned down. The reverse order is 11.1% (CI -1.1 to 28.2): most of the CLDN4 association with the CD8 score is direct, and that percentage's interval includes zero.

The TCGA sum is not a 7-cohort vote. Cohort-specific BIC gaps on the primary CD8 specification:

| cohort | n | BIC winner | ΔBIC to 2nd | T→C→immune % | C→T→immune % |
| --- | --- | --- | --- | --- | --- |
| LUAD | 516 | M1 | 2.02 | 33.3 | 20.4 |
| BRCA | 1095 | M2 | 0.70 | -28.1 (opposite signs) | -77.5 (unstable total; opposite signs) |
| CESC | 304 | M2 | 0.10 | -24.8 (opposite signs) | -25.9 (opposite signs) |
| KIRC | 533 | M2 | 6.32 | 2.8 | 74.6 |
| STAD | 412 | M1 | 9.50 | 54.4 | 3.1 |
| BLCA | 406 | M1 | 2.54 | 35.9 | 26.4 |
| PAAD | 178 | M3 | 6.96 | -28.2 (opposite signs) | -19.2 (opposite signs) |

LUSC alone, same CD8 score and KRT8/18/19 residual, prefers M2 by ΔBIC 33.62 (n = 501). Adding LUSC to the seven-cohort sum selects M2 (ΔBIC 22.03). LUSC was not in the pre-specified funnel.

Replacing KRT8/18/19 with KRT5+KRT6A+KRT6B+KRT14 on the same seven cohorts leaves the cohort-sum gap at ΔBIC 0.29 for M1, while the pooled SEM prefers M2 (ΔBIC 3.34). The pooled residual correlations are -0.071 for TACSTD2 after CLDN4 and -0.064 for CLDN4 after TACSTD2. Under that keratin specification the two directions are not separated.

Concordant-4 datasets, fit separately: GSE123902 M1 (ΔBIC to 2nd 2.65, n=13), GSE131907 M1 (ΔBIC to 2nd 1.59, n=21), GSE189357 M1 (ΔBIC to 2nd 1.66, n=9), GSE205335 M3 (ΔBIC to 2nd 2.24, n=22). GSE205335 is the one dataset whose BIC winner is not M1, and that gap is small. The pooled n=65 result is where the ΔBIC exceeds 10.

## Two-stage residual

Conditional-independence check: residualize the immune variable on the hypothesized mediator, then correlate that residual with the other gene. M1 says the TACSTD2 residual correlation should be the one near zero. M2 says the CLDN4 residual correlation should be the one near zero. M3 says the TACSTD2–CLDN4 correlation itself should be near zero.

Effect decomposition: residualize the mediator on the exposure, then regress the immune variable on the exposure and that residual. The exposure coefficient equals the total effect, and the residual coefficient equals the mediator coefficient. On these matrices the identity holds to numerical error (the script checks it before writing results).

Concordant-4 primary: Pearson residual correlation of TACSTD2 with immune after CLDN4 is 0.109 (p = 0.3858). Pearson residual correlation of CLDN4 with immune after TACSTD2 is -0.460 (p = 0.0001). Pearson correlation of TACSTD2 with CLDN4 is 0.500 (p = 2.19e-05). Within-stratum partial Spearman, DerSimonian–Laird: TACSTD2–immune | CLDN4 DL ρ = 0.216 (p = 0.1239, I² = 0%); CLDN4–immune | TACSTD2 DL ρ = -0.543 (p = 2.02e-05, I² = 0%); TACSTD2–CLDN4 DL ρ = 0.535 (p = 0.0058, I² = 56%).

TCGA primary pooled residuals: Pearson residual correlation of TACSTD2 with immune after CLDN4 is -0.029 (p = 0.0859). Pearson residual correlation of CLDN4 with immune after TACSTD2 is -0.071 (p = 3.38e-05). Pearson correlation of TACSTD2 with CLDN4 is 0.300 (p = 9.03e-73). Within-stratum partial Spearman, DerSimonian–Laird: TACSTD2–immune | CLDN4 DL ρ = -0.007 (p = 0.8480, I² = 78%); CLDN4–immune | TACSTD2 DL ρ = -0.093 (p = 0.0176, I² = 79%); TACSTD2–CLDN4 DL ρ = 0.318 (p = 3.54e-19, I² = 76%).

## Sensitivities

The primary specifications above were fixed first. The rows below use the same three graphs. A sensitivity that picks a different graph is listed as a disagreement, not averaged away.

| role | spec | estimator | n | BIC winner | ΔBIC 2nd | saturated ΔBIC | M1 mediation % | M2 mediation % |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| primary | c4_pct | pooled_within_z | 65 | M1 | 14.66 | 3.39 | 173.8 | -11.9 |
| sensitivity | c4_mean | pooled_within_z | 65 | M1 | 5.74 | 2.63 | 548.3 | -33.1 |
| sensitivity | c4_pct_logit | pooled_within_z | 65 | M1 | 17.04 | 2.49 | 215.7 | -15.5 |
| sensitivity | c4_pct_min30_malignant | pooled_within_z | 64 | M1 | 13.96 | 3.14 | 196.9 | -14.1 |
| primary | tcga_cd8_krt819 | cohort_bic_sum | 3444 | M1 | 11.60 | 11.08 | 42.1 | 11.1 |
| sensitivity | tcga_cd3_krt819 | cohort_bic_sum | 3444 | M1 | 20.42 | 16.81 | 107.8 | -0.7 |
| sensitivity | tcga_cytotoxic_krt819 | cohort_bic_sum | 3444 | M1 | 50.51 | 13.34 | 60.0 | 5.7 |
| sensitivity | tcga_cd8_none | cohort_bic_sum | 3444 | M1 | 16.85 | 7.43 | 41.2 | 19.9 |
| sensitivity | tcga_cd8_krt56 | cohort_bic_sum | 3444 | M1 | 0.29 | -1.77 | 28.0 | 32.5 |
| sensitivity | tcga_cd3_krt56 | cohort_bic_sum | 3444 | M1 | 23.41 | 15.73 | 51.7 | 14.9 |
| sensitivity | tcga_cd8_krt819_plus_lusc | cohort_bic_sum | 3945 | M2 | 22.03 | 2.33 | 23.1 | 24.8 |

TCGA primary, cohort-by-cohort BIC winners: M1 3, M2 3, M3 1 of 7 funnel cohorts. Per-cohort labels: LUAD M1, BRCA M2, CESC M2, KIRC M2, STAD M1, BLCA M1, PAAD M3.

Concordant-4 primary, dataset-by-dataset BIC winners: GSE123902 M1, GSE131907 M1, GSE189357 M1, GSE205335 M3.

## Calibration

Concordant-4 DerSimonian–Laird Spearman of malignant CLDN4 percent positive versus T/NK fraction is -0.531168 (p = 1.65e-05, I² = 0%, N = 65). That matches the locked value −0.5311678045689989 from the concordant-4 patient table.
The same meta-analysis on mean log1p is -0.394 for CLDN4 and 0.041 for TACSTD2; percent-positive TACSTD2 is -0.112.
LUAD partial Spearman of TACSTD2 versus the CD8 score after KRT8/18/19 is -0.103504 (published −0.103504), and CLDN4 versus the CD8 score is -0.084776 (published −0.084776). Patient counts match the earlier Xena extract (funnel n = 3,444).

## What this does not identify

Once both the TACSTD2–CLDN4 arrow and the direct arrow into immune are free, the saturated Gaussian model fits the 3×3 covariance exactly. The two directions of the gene–gene arrow are then the same likelihood, and AIC cannot choose between them. Only the restricted graphs, which drop one association, are separated by AIC and BIC.

Concordant-4 percent-positive and T/NK fraction are bounded. The linear SEM is a covariance approximation, which is why the mean log1p and logit-fraction sensitivities are in the table. TCGA CD8 is a bulk transcript score after a linear keratin residual, not a spatial immune count and not a purity-adjusted fraction.

