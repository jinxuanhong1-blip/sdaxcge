# Methods

Three observed variables, one Gaussian structural equation model per graph. Means are removed. Variances use the maximum-likelihood divisor n. The log-likelihood is the recursive factorization (each node given its parents), which matches the multivariate-normal likelihood of the implied covariance.

## Graphs

| code | edges | constraint |
|---|---|---|
| M1 | TACSTD2 → CLDN4 → immune | no direct TACSTD2 → immune |
| M2 | CLDN4 → TACSTD2 → immune | no direct CLDN4 → immune |
| M3 | TACSTD2 → immune ← CLDN4 | TACSTD2 independent of CLDN4 |
| saturated | gene–gene edge plus both arrows into immune | none (6 covariance parameters) |

M1, M2, and M3 each have 5 free parameters. AIC = −2 log L + 2k. BIC = −2 log L + k log n. AICc is reported for the single-sample fits. The likelihood-ratio statistic versus the saturated model is compared to χ² on the difference in parameter count.

ΔBIC between models is unchanged if each variable is rescaled, because every model is fit to the same transformed data. Z-scoring is used so the path coefficients are standardized. The script checks that the BIC winner on the z-scored matrix matches the winner on the unscaled residuals whenever the BIC gap is larger than 0.05.

## Mediation percentage

The reported percentage keeps the direct path. For exposure X, mediator M, and immune Y:

- M = a X + error
- Y = c′ X + b M + error
- indirect = a b
- total = c′ + a b
- percent = 100 × indirect / total

On centered ordinary least squares, c′ + a b equals the coefficient from Y on X alone. The two-stage residual inclusion regression (residualize M on X, then regress Y on X and that residual) returns the same total effect and the same b. The script refuses to write results if that identity fails.

A chain that omits the direct path is 100% mediated by construction. That figure is not used.

The percentage is not bounded to 0–100 when the direct and indirect paths have opposite signs. A standardized total effect below 0.02 in absolute value is marked unstable.

Uncertainty is a stratified bootstrap: resample patients inside each dataset or cohort, repeat z-scoring and (for TCGA) keratin residualization, and refit. 2,000 replicates, seed 20260921, percentile interval 2.5 to 97.5.

## Two-stage residual

Residualize the immune variable on the hypothesized mediator and correlate the residual with the other gene. M1 predicts that the TACSTD2 residual correlation is the one near zero. M2 predicts that for CLDN4. M3 predicts that the TACSTD2–CLDN4 correlation is near zero. These residual correlations are Pearson correlations on the analysis matrix. A DerSimonian–Laird partial Spearman across strata is reported beside them.

## Concordant-4 primary specification

Malignant percent positive for TACSTD2 and CLDN4, outcome `frac_tnk`, n = 65. Z-score inside each of the four datasets, then one SEM. The sum of four dataset-specific BICs is a second estimator, not the primary one (GSE189357 has 9 units).

Sensitivities, not used to choose the primary winner: mean log1p instead of percent positive; logit of the T/NK fraction; units with at least 30 malignant cells.

## TCGA primary specification

Funnel cohorts LUAD, BRCA, CESC, KIRC, STAD, BLCA, PAAD. CD8 score = mean(CD8A, CD8B). Within each cohort, residualize TACSTD2, CLDN4, and the CD8 score on KRT8, KRT18, and KRT19, then z-score. The primary comparison sums the cohort BICs. A pooled z-score SEM is reported next to it. LUSC is fit and is excluded from the sum.

Sensitivities: CD3 score = mean(CD3D, CD3E, CD3G); cytotoxic score = mean(GZMA, GZMB, PRF1, NKG7); no keratin residual; KRT5+KRT6A+KRT6B+KRT14; funnel plus LUSC.

## Software checks

`test_semcore.py` checks that the factorized likelihood matches the multivariate-normal likelihood, that mediation and the two-stage residual agree, that rescaling does not change ΔBIC, and that data simulated from each graph are assigned to that graph by BIC.
