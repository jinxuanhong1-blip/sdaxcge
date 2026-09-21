# Finding — CPTAC protein: TROP2, partial CLDN4, CD8A/MHC-I, EPCAM control

Public CPTAC TMT freeze v1.2 tumor protein. **LUAD** (Gillette *Cell* 2020, n=110) and **LSCC** (Satpathy *Cell* 2021, n=108) kept separate. TROP2 is the protein product of TACSTD2. No log2 abundance is imputed. Coefficients below are written by `analyze.py` from `tables/`; displayed values are rounded to 3 decimals.

All-tumor TROP2 versus CD8A and versus MHC-I: 0 of 4 bootstrap CIs lie entirely below 0, 0 entirely above 0, 4 cross 0. Parallel-model indirect effects on the shared complete cases: CLDN4 mediator, 0 of 4 entirely below 0; EPCAM mediator, 0 of 4 entirely below 0.

## What the paths show

TROP2 versus CD8A and MHC-I, using every tumor in which both proteins were quantified: LUAD TROP2 vs CD8A ρ=-0.038 (n=110, 95% CI -0.247 to +0.159, crosses 0); LUAD TROP2 vs MHC1 ρ=-0.083 (n=110, 95% CI -0.276 to +0.119, crosses 0); LSCC TROP2 vs CD8A ρ=-0.103 (n=108, 95% CI -0.288 to +0.085, crosses 0); LSCC TROP2 vs MHC1 ρ=+0.029 (n=108, 95% CI -0.172 to +0.224, crosses 0).

In LSCC the inverse CD8A and MHC-I associations are with CLDN4 and with EPCAM. CLDN4 vs CD8A ρ=-0.444 (n=78, 95% CI -0.626 to -0.235, entirely below 0); EPCAM vs CD8A ρ=-0.396 (n=108, 95% CI -0.564 to -0.205, entirely below 0); CLDN4 vs MHC-I ρ=-0.410 (n=78, 95% CI -0.586 to -0.197, entirely below 0); EPCAM vs MHC-I ρ=-0.385 (n=108, 95% CI -0.543 to -0.209, entirely below 0).

Each of those LSCC associations remains after the other epithelial protein, or TROP2, is held constant. Partial Spearman on the CLDN4-quantified tumors: CLDN4 vs CD8A given EPCAM partial ρ=-0.396 (n=78, 95% CI -0.597 to -0.163, q=0.0015, entirely below 0); CLDN4 vs CD8A given TROP2 partial ρ=-0.440 (n=78, 95% CI -0.615 to -0.222, q=6.7e-04, entirely below 0); EPCAM vs CD8A given CLDN4 partial ρ=-0.393 (n=78, 95% CI -0.595 to -0.152, q=0.0015, entirely below 0); CLDN4 vs MHC-I given EPCAM partial ρ=-0.356 (n=78, 95% CI -0.556 to -0.135, q=0.00413, entirely below 0); CLDN4 vs MHC-I given TROP2 partial ρ=-0.417 (n=78, 95% CI -0.591 to -0.207, q=0.00106, entirely below 0); EPCAM vs MHC-I given CLDN4 partial ρ=-0.408 (n=78, 95% CI -0.579 to -0.197, q=0.00126, entirely below 0).

The a path, TROP2 to the epithelial protein, is histology-specific. LUAD TROP2 vs CLDN4 ρ=+0.275 (n=79, 95% CI +0.069 to +0.473, entirely above 0); LUAD TROP2 vs EPCAM ρ=+0.372 (n=110, 95% CI +0.189 to +0.547, entirely above 0); LSCC TROP2 vs CLDN4 ρ=+0.081 (n=78, 95% CI -0.168 to +0.323, crosses 0); LSCC TROP2 vs EPCAM ρ=+0.085 (n=108, 95% CI -0.107 to +0.287, crosses 0); LUAD CLDN4 vs EPCAM ρ=+0.195 (n=79, 95% CI -0.030 to +0.401, crosses 0); LSCC CLDN4 vs EPCAM ρ=+0.221 (n=78, 95% CI -0.010 to +0.449, crosses 0). Max VIF of the three rank predictors is 1.12 in LUAD and 1.08 in LSCC.

Shared-sample indirect effects: 16 of 16 ab confidence intervals cross 0. LSCC b intervals entirely below 0: 8 of 8. LSCC TROP2→CLDN4 a intervals that cross 0: 4 of 4. LUAD CLDN4-mediator a intervals entirely above 0: 4 of 4. LUAD CLDN4-mediator b intervals that cross 0: 4 of 4. Primary q on the 16 indirect tests has minimum 0.559.

The one shared-sample total effect whose interval lies entirely below 0 is LUAD TROP2 vs MHC-I on the CLDN4-quantified tumors: c=-0.245 (n=79, 95% CI -0.448 to -0.013). The all-tumor TROP2–MHC-I interval crosses 0, so this is the CLDN4-complete subset, not the full LUAD series. TROP2 vs MHC-I given CLDN4 partial ρ=-0.269 (n=79, 95% CI -0.465 to -0.059, q=0.046, entirely below 0). The same contrast given CLDN4 and EPCAM partial ρ=-0.250 (n=79, 95% CI -0.448 to -0.038, q=0.0677, entirely below 0). The CLDN4 indirect effect on that subset is ab=+0.035 (95% CI -0.031 to +0.119), ab/c=-0.142. The product is positive because the CLDN4 b coefficient is positive. The inverse total association on this subset is the direct path. WES-purity rank adjustment of the parallel model on this subset (n=77) moves the total-effect interval to -0.412 to +0.030 (crosses 0) and the direct-effect interval to -0.440 to +4.7e-04 (crosses 0).

WES-purity rank adjustment: 16 of 16 indirect-effect intervals cross 0. LSCC parallel-model b coefficients under that adjustment: CD8A CLDN4 b=-0.348 (CI -0.558 to -0.120, entirely below 0); CD8A EPCAM b=-0.237 (CI -0.520 to +0.049, crosses 0); MHC-I CLDN4 b=-0.291 (CI -0.485 to -0.078, entirely below 0); MHC-I EPCAM b=-0.197 (CI -0.459 to +0.036, crosses 0). Unadjusted partial Spearman had already kept both LSCC epithelial proteins inverse; the purity-adjusted regression coefficient is a different estimand and is where EPCAM and CLDN4 separate.

The pre-specified mediation rule, which requires an inverse total TROP2 association and an inverse CLDN4 indirect effect that EPCAM does not copy, is met for 0 of 4 cohort–outcome pairs.

## Treatment-naive

LUAD phenotype has 183 columns. Therapy-like column names found: none.
LSCC phenotype has 193 columns. Therapy-like column names found: none.

Treatment-naive is the published cohort definition: prospectively collected, previously untreated surgical resections. There is no ICI arm, no on-treatment biopsy, and no response label in these files. A path coefficient here is not acquired resistance and not a drug effect.

## Estimands

Exposure X = TACSTD2 protein. Mediator under test M1 = CLDN4 protein. Control protein M2 = EPCAM, a second epithelial surface protein run through the same equations. Outcomes are CD8A protein and the MHC-I score used on the earlier page for this freeze: the mean of per-gene z-scores (ddof=0) of HLA-A, HLA-B, and HLA-C, requiring all three. B2M is not in that score. Ranks use average ties. Single-mediator z-rank OLS is M ~ X and Y ~ X + M. The product ab equals c − c′ on that scale, and a and c equal Spearman ρ. The parallel model is M1 ~ X, M2 ~ X, and Y ~ X + M1 + M2, so each indirect effect is adjusted for the other epithelial protein at the b step. Bootstrap: 2,000 resamples, seed `20260921`, ranks recomputed inside each resample, percentile 95% CI. Bootstrap p for ab is two-sided on the resampled products; a value of 0 is stored as 1/2000 and flagged. Primary q is Benjamini–Hochberg across the 16 shared-sample indirect tests (2 cohorts × 2 outcomes × 2 mediators × single and parallel). Partial Spearman is the Pearson correlation of rank residuals and has its own BH family (the covariate list in `PARTIAL_SPECS`). WES purity is a rank covariate in a sensitivity refit, not part of the primary rule.

The mediation rule, applied to the parallel model on tumors with TACSTD2, CLDN4, EPCAM, and the outcome all quantified, requires all of: CLDN4 a CI entirely above 0, CLDN4 b CI entirely below 0, CLDN4 ab CI entirely below 0, total c CI entirely below 0, and EPCAM ab CI not entirely below 0. Where a total-effect CI includes 0, ab/c is left in `tables/mediation_shared.tsv` and is not read as a fraction mediated.

## Reproduction of the eight earlier totals

Pairwise-complete Spearman ρ and analytic p for TACSTD2 and CLDN4 versus CD8A and MHC-I were compared with `summary.json` from the earlier protein page, which used these same freeze files. All eight match within 1e-8 on n, ρ, and p. That page's q and its WES-purity partial correlations are not re-adjudicated here.

| Cohort | X | Y | n | ρ | |Δρ| | Match |
| --- | --- | --- | --- | --- | --- | --- |
| LUAD | CLDN4 | CD8A | 79 | -0.098 | 0.0e+00 | yes |
| LUAD | CLDN4 | MHC1 | 79 | +0.049 | 0.0e+00 | yes |
| LUAD | TACSTD2 | CD8A | 110 | -0.038 | 0.0e+00 | yes |
| LUAD | TACSTD2 | MHC1 | 110 | -0.083 | 0.0e+00 | yes |
| LSCC | CLDN4 | CD8A | 78 | -0.444 | 0.0e+00 | yes |
| LSCC | CLDN4 | MHC1 | 78 | -0.410 | 0.0e+00 | yes |
| LSCC | TACSTD2 | CD8A | 108 | -0.103 | 0.0e+00 | yes |
| LSCC | TACSTD2 | MHC1 | 108 | +0.029 | 0.0e+00 | yes |

## Partial quantification

CLDN4 is the incompletely quantified protein. Mediation uses complete cases only. Tumors missing CLDN4 stay in the all-tumor TROP2 and EPCAM scatters and are dropped from every path that includes CLDN4. The B2M row is absent from both matrices, which is why the MHC-I score stays HLA-A/B/C.

| Cohort | Protein | Matrix row | Quantified | Missing | % missing |
| --- | --- | --- | --- | --- | --- |
| LUAD | TACSTD2 | ENSG00000184292.7 | 110 | 0 | 0.0 |
| LUAD | CLDN4 | ENSG00000189143.9 | 79 | 31 | 28.2 |
| LUAD | EPCAM | ENSG00000119888.11 | 110 | 0 | 0.0 |
| LUAD | CD8A | ENSG00000153563.15 | 110 | 0 | 0.0 |
| LUAD | HLA-A | ENSG00000206503.13 | 110 | 0 | 0.0 |
| LUAD | HLA-B | ENSG00000234745.11 | 110 | 0 | 0.0 |
| LUAD | HLA-C | ENSG00000204525.16 | 110 | 0 | 0.0 |
| LUAD | B2M | nan | 0 | 110 | 100.0 |
| LSCC | TACSTD2 | ENSG00000184292.7 | 108 | 0 | 0.0 |
| LSCC | CLDN4 | ENSG00000189143.9 | 78 | 30 | 27.8 |
| LSCC | EPCAM | ENSG00000119888.11 | 108 | 0 | 0.0 |
| LSCC | CD8A | ENSG00000153563.15 | 108 | 0 | 0.0 |
| LSCC | HLA-A | ENSG00000206503.13 | 108 | 0 | 0.0 |
| LSCC | HLA-B | ENSG00000234745.11 | 108 | 0 | 0.0 |
| LSCC | HLA-C | ENSG00000204525.16 | 108 | 0 | 0.0 |
| LSCC | B2M | nan | 0 | 108 | 100.0 |

Shared complete-case n (TACSTD2, CLDN4, EPCAM, and the outcome):

- LUAD CD8A: n=79
- LUAD MHC1: n=79
- LSCC CD8A: n=78
- LSCC MHC1: n=78

## All-tumor Spearman

These are pairwise-complete. They are not the mediation-sample total effects. Analytic Spearman p is the usual t approximation. The CI is this page's bootstrap. No FDR is applied to this descriptive table.

| Cohort | X | Y | n | ρ | CI low | CI high | p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| LUAD | TACSTD2 | CD8A | 110 | -0.038 | -0.247 | +0.159 | 0.691 |
| LUAD | TACSTD2 | MHC1 | 110 | -0.083 | -0.276 | +0.119 | 0.387 |
| LUAD | CLDN4 | CD8A | 79 | -0.098 | -0.328 | +0.145 | 0.391 |
| LUAD | CLDN4 | MHC1 | 79 | +0.049 | -0.190 | +0.279 | 0.665 |
| LUAD | EPCAM | CD8A | 110 | -0.135 | -0.327 | +0.069 | 0.16 |
| LUAD | EPCAM | MHC1 | 110 | -0.030 | -0.232 | +0.169 | 0.754 |
| LUAD | TACSTD2 | CLDN4 | 79 | +0.275 | +0.069 | +0.473 | 0.0143 |
| LUAD | TACSTD2 | EPCAM | 110 | +0.372 | +0.189 | +0.547 | 6.4e-05 |
| LUAD | CLDN4 | EPCAM | 79 | +0.195 | -0.030 | +0.401 | 0.0852 |
| LSCC | TACSTD2 | CD8A | 108 | -0.103 | -0.288 | +0.085 | 0.287 |
| LSCC | TACSTD2 | MHC1 | 108 | +0.029 | -0.172 | +0.224 | 0.764 |
| LSCC | CLDN4 | CD8A | 78 | -0.444 | -0.626 | -0.235 | 4.6e-05 |
| LSCC | CLDN4 | MHC1 | 78 | -0.410 | -0.586 | -0.197 | 2.0e-04 |
| LSCC | EPCAM | CD8A | 108 | -0.396 | -0.564 | -0.205 | 2.2e-05 |
| LSCC | EPCAM | MHC1 | 108 | -0.385 | -0.543 | -0.209 | 3.8e-05 |
| LSCC | TACSTD2 | CLDN4 | 78 | +0.081 | -0.168 | +0.323 | 0.483 |
| LSCC | TACSTD2 | EPCAM | 108 | +0.085 | -0.107 | +0.287 | 0.382 |
| LSCC | CLDN4 | EPCAM | 78 | +0.221 | -0.010 | +0.449 | 0.0519 |

## Mediation on the shared sample

a is the standardized rank coefficient, equal to Spearman ρ for the single-mediator and parallel a paths (those a paths are M ~ X, not adjusted for the other mediator). b and c′ are standardized rank regression coefficients. ab is their product. q is the 16-test primary BH on ab.

| Cohort | Y | Mediator | Model | n | a | b | c | c′ | ab | ab CI low | ab CI high | ab p | q |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LUAD | CD8A | CLDN4 | single | 79 | +0.275 | -0.104 | -0.005 | +0.023 | -0.029 | -0.106 | +0.048 | 0.44 | 0.559 |
| LUAD | CD8A | EPCAM | single | 79 | +0.223 | -0.134 | -0.005 | +0.024 | -0.030 | -0.097 | +0.033 | 0.368 | 0.559 |
| LUAD | CD8A | CLDN4 | parallel | 79 | +0.275 | -0.087 | -0.005 | +0.045 | -0.024 | -0.102 | +0.053 | 0.524 | 0.559 |
| LUAD | CD8A | EPCAM | parallel | 79 | +0.223 | -0.121 | -0.005 | +0.045 | -0.027 | -0.093 | +0.036 | 0.397 | 0.559 |
| LUAD | MHC1 | CLDN4 | single | 79 | +0.275 | +0.126 | -0.245 | -0.279 | +0.035 | -0.031 | +0.119 | 0.257 | 0.559 |
| LUAD | MHC1 | EPCAM | single | 79 | +0.223 | -0.078 | -0.245 | -0.227 | -0.017 | -0.083 | +0.043 | 0.581 | 0.581 |
| LUAD | MHC1 | CLDN4 | parallel | 79 | +0.275 | +0.140 | -0.245 | -0.261 | +0.039 | -0.030 | +0.121 | 0.22 | 0.559 |
| LUAD | MHC1 | EPCAM | parallel | 79 | +0.223 | -0.098 | -0.245 | -0.261 | -0.022 | -0.090 | +0.036 | 0.497 | 0.559 |
| LSCC | CD8A | CLDN4 | single | 78 | +0.081 | -0.439 | -0.097 | -0.062 | -0.035 | -0.146 | +0.065 | 0.503 | 0.559 |
| LSCC | CD8A | EPCAM | single | 78 | +0.173 | -0.438 | -0.097 | -0.021 | -0.076 | -0.202 | +0.026 | 0.185 | 0.559 |
| LSCC | CD8A | CLDN4 | parallel | 78 | +0.081 | -0.364 | -0.097 | -0.005 | -0.029 | -0.126 | +0.058 | 0.504 | 0.559 |
| LSCC | CD8A | EPCAM | parallel | 78 | +0.173 | -0.360 | -0.097 | -0.005 | -0.062 | -0.172 | +0.021 | 0.188 | 0.559 |
| LSCC | MHC1 | CLDN4 | single | 78 | +0.081 | -0.417 | +0.061 | +0.095 | -0.034 | -0.129 | +0.060 | 0.495 | 0.559 |
| LSCC | MHC1 | EPCAM | single | 78 | +0.173 | -0.478 | +0.061 | +0.144 | -0.083 | -0.220 | +0.027 | 0.161 | 0.559 |
| LSCC | MHC1 | CLDN4 | parallel | 78 | +0.081 | -0.332 | +0.061 | +0.159 | -0.027 | -0.104 | +0.051 | 0.495 | 0.559 |
| LSCC | MHC1 | EPCAM | parallel | 78 | +0.173 | -0.407 | +0.061 | +0.159 | -0.071 | -0.187 | +0.025 | 0.162 | 0.559 |

Parallel-model detail and the mediation rule:

**LSCC CD8A, n=78.** Total c -0.097 (CI -0.311 to +0.130, analytic p=0.398), which crosses 0. CLDN4 a +0.081 (CI -0.155 to +0.302). Single-mediator CLDN4 b -0.439 (CI -0.619 to -0.207), c′ -0.062 (CI -0.255 to +0.139), ab -0.035 (CI -0.146 to +0.065, p=0.503, q=0.559). Single-mediator EPCAM ab -0.076 (CI -0.202 to +0.026, q=0.559). Parallel CLDN4 b -0.364 (CI -0.554 to -0.147), ab -0.029 (CI -0.126 to +0.058, p=0.504, q=0.559). Parallel EPCAM b -0.360 (CI -0.565 to -0.121), ab -0.062 (CI -0.172 to +0.021, q=0.559). Parallel direct c′ -0.005 (CI -0.178 to +0.172). Max VIF of TACSTD2, CLDN4, and EPCAM ranks in the parallel outcome model: 1.08. Rule clauses met: CLDN4 b CI entirely below 0, EPCAM ab CI not entirely below 0. Rule clauses not met: CLDN4 a CI entirely above 0, CLDN4 ab CI entirely below 0, total c CI entirely below 0.

**LSCC MHC1, n=78.** Total c +0.061 (CI -0.184 to +0.298, analytic p=0.594), which crosses 0. CLDN4 a +0.081 (CI -0.150 to +0.298). Single-mediator CLDN4 b -0.417 (CI -0.593 to -0.206), c′ +0.095 (CI -0.121 to +0.301), ab -0.034 (CI -0.129 to +0.060, p=0.495, q=0.559). Single-mediator EPCAM ab -0.083 (CI -0.220 to +0.027, q=0.559). Parallel CLDN4 b -0.332 (CI -0.518 to -0.128), ab -0.027 (CI -0.104 to +0.051, p=0.495, q=0.559). Parallel EPCAM b -0.407 (CI -0.592 to -0.212), ab -0.071 (CI -0.187 to +0.025, q=0.559). Parallel direct c′ +0.159 (CI -0.033 to +0.341). Max VIF of TACSTD2, CLDN4, and EPCAM ranks in the parallel outcome model: 1.08. Rule clauses met: CLDN4 b CI entirely below 0, EPCAM ab CI not entirely below 0. Rule clauses not met: CLDN4 a CI entirely above 0, CLDN4 ab CI entirely below 0, total c CI entirely below 0.

**LUAD CD8A, n=79.** Total c -0.005 (CI -0.259 to +0.249, analytic p=0.962), which crosses 0. CLDN4 a +0.275 (CI +0.050 to +0.469). Single-mediator CLDN4 b -0.104 (CI -0.340 to +0.154), c′ +0.023 (CI -0.245 to +0.273), ab -0.029 (CI -0.106 to +0.048, p=0.44, q=0.559). Single-mediator EPCAM ab -0.030 (CI -0.097 to +0.033, q=0.559). Parallel CLDN4 b -0.087 (CI -0.322 to +0.166), ab -0.024 (CI -0.102 to +0.053, p=0.524, q=0.559). Parallel EPCAM b -0.121 (CI -0.352 to +0.130), ab -0.027 (CI -0.093 to +0.036, q=0.559). Parallel direct c′ +0.045 (CI -0.218 to +0.277). Max VIF of TACSTD2, CLDN4, and EPCAM ranks in the parallel outcome model: 1.12. Rule clauses met: CLDN4 a CI entirely above 0, EPCAM ab CI not entirely below 0. Rule clauses not met: CLDN4 b CI entirely below 0, CLDN4 ab CI entirely below 0, total c CI entirely below 0.

**LUAD MHC1, n=79.** Total c -0.245 (CI -0.448 to -0.013, analytic p=0.0299), which entirely below 0. CLDN4 a +0.275 (CI +0.051 to +0.463). Single-mediator CLDN4 b +0.126 (CI -0.111 to +0.350), c′ -0.279 (CI -0.480 to -0.058), ab +0.035 (CI -0.031 to +0.119, p=0.257, q=0.559). Single-mediator EPCAM ab -0.017 (CI -0.083 to +0.043, q=0.581). Parallel CLDN4 b +0.140 (CI -0.100 to +0.365), ab +0.039 (CI -0.030 to +0.121, p=0.22, q=0.559). Parallel EPCAM b -0.098 (CI -0.322 to +0.139), ab -0.022 (CI -0.090 to +0.036, q=0.559). Parallel direct c′ -0.261 (CI -0.471 to -0.047). Max VIF of TACSTD2, CLDN4, and EPCAM ranks in the parallel outcome model: 1.12. Rule clauses met: CLDN4 a CI entirely above 0, total c CI entirely below 0, EPCAM ab CI not entirely below 0. Rule clauses not met: CLDN4 b CI entirely below 0, CLDN4 ab CI entirely below 0.

## Partial Spearman

Partial ρ is the correlation of rank residuals, with analytic p from that Pearson correlation and a bootstrap CI. q is BH inside this covariate family only. A partial CI that excludes 0 is not promoted into the mediation rule.

| Cohort | Focal | Y | Given | n | partial ρ | CI low | CI high | p | q |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LUAD | TACSTD2 | CD8A | CLDN4 | 79 | +0.022 | -0.236 | +0.254 | 0.845 | 0.922 |
| LUAD | TACSTD2 | CD8A | EPCAM | 110 | +0.013 | -0.170 | +0.200 | 0.895 | 0.948 |
| LUAD | TACSTD2 | CD8A | CLDN4+EPCAM | 79 | +0.043 | -0.188 | +0.271 | 0.704 | 0.792 |
| LUAD | TACSTD2 | MHC1 | CLDN4 | 79 | -0.269 | -0.465 | -0.059 | 0.0166 | 0.046 |
| LUAD | TACSTD2 | MHC1 | EPCAM | 110 | -0.078 | -0.265 | +0.112 | 0.42 | 0.6 |
| LUAD | TACSTD2 | MHC1 | CLDN4+EPCAM | 79 | -0.250 | -0.448 | -0.038 | 0.0263 | 0.0677 |
| LUAD | CLDN4 | CD8A | TACSTD2 | 79 | -0.100 | -0.345 | +0.151 | 0.379 | 0.583 |
| LUAD | CLDN4 | CD8A | EPCAM | 79 | -0.075 | -0.324 | +0.168 | 0.511 | 0.614 |
| LUAD | CLDN4 | CD8A | TACSTD2+EPCAM | 79 | -0.083 | -0.330 | +0.161 | 0.466 | 0.614 |
| LUAD | CLDN4 | MHC1 | TACSTD2 | 79 | +0.125 | -0.106 | +0.336 | 0.272 | 0.516 |
| LUAD | CLDN4 | MHC1 | EPCAM | 79 | +0.077 | -0.155 | +0.297 | 0.502 | 0.614 |
| LUAD | CLDN4 | MHC1 | TACSTD2+EPCAM | 79 | +0.138 | -0.096 | +0.343 | 0.225 | 0.449 |
| LUAD | EPCAM | CD8A | CLDN4 | 79 | -0.112 | -0.339 | +0.118 | 0.327 | 0.56 |
| LUAD | EPCAM | CD8A | TACSTD2 | 110 | -0.130 | -0.320 | +0.078 | 0.176 | 0.396 |
| LUAD | EPCAM | CD8A | CLDN4+TACSTD2 | 79 | -0.118 | -0.336 | +0.127 | 0.301 | 0.542 |
| LUAD | EPCAM | MHC1 | CLDN4 | 79 | -0.141 | -0.365 | +0.078 | 0.214 | 0.449 |
| LUAD | EPCAM | MHC1 | TACSTD2 | 110 | +0.001 | -0.183 | +0.178 | 0.993 | 0.993 |
| LUAD | EPCAM | MHC1 | CLDN4+TACSTD2 | 79 | -0.098 | -0.330 | +0.133 | 0.388 | 0.583 |
| LSCC | TACSTD2 | CD8A | CLDN4 | 78 | -0.069 | -0.287 | +0.146 | 0.55 | 0.639 |
| LSCC | TACSTD2 | CD8A | EPCAM | 108 | -0.076 | -0.256 | +0.104 | 0.433 | 0.6 |
| LSCC | TACSTD2 | CD8A | CLDN4+EPCAM | 78 | -0.006 | -0.221 | +0.196 | 0.955 | 0.983 |
| LSCC | TACSTD2 | MHC1 | CLDN4 | 78 | +0.104 | -0.125 | +0.330 | 0.366 | 0.583 |
| LSCC | TACSTD2 | MHC1 | EPCAM | 108 | +0.067 | -0.128 | +0.249 | 0.488 | 0.614 |
| LSCC | TACSTD2 | MHC1 | CLDN4+EPCAM | 78 | +0.187 | -0.048 | +0.402 | 0.1 | 0.241 |
| LSCC | CLDN4 | CD8A | TACSTD2 | 78 | -0.440 | -0.615 | -0.222 | 5.6e-05 | 6.7e-04 |
| LSCC | CLDN4 | CD8A | EPCAM | 78 | -0.396 | -0.597 | -0.163 | 3.3e-04 | 0.0015 |
| LSCC | CLDN4 | CD8A | TACSTD2+EPCAM | 78 | -0.395 | -0.582 | -0.154 | 3.4e-04 | 0.0015 |
| LSCC | CLDN4 | MHC1 | TACSTD2 | 78 | -0.417 | -0.591 | -0.207 | 1.5e-04 | 0.00106 |
| LSCC | CLDN4 | MHC1 | EPCAM | 78 | -0.356 | -0.556 | -0.135 | 0.00138 | 0.00413 |
| LSCC | CLDN4 | MHC1 | TACSTD2+EPCAM | 78 | -0.368 | -0.554 | -0.161 | 9.1e-04 | 0.00299 |
| LSCC | EPCAM | CD8A | CLDN4 | 78 | -0.393 | -0.595 | -0.152 | 3.8e-04 | 0.0015 |
| LSCC | EPCAM | CD8A | TACSTD2 | 108 | -0.391 | -0.561 | -0.209 | 2.8e-05 | 5.6e-04 |
| LSCC | EPCAM | CD8A | CLDN4+TACSTD2 | 78 | -0.388 | -0.590 | -0.140 | 4.5e-04 | 0.00164 |
| LSCC | EPCAM | MHC1 | CLDN4 | 78 | -0.408 | -0.579 | -0.197 | 2.1e-04 | 0.00126 |
| LSCC | EPCAM | MHC1 | TACSTD2 | 108 | -0.389 | -0.546 | -0.210 | 3.1e-05 | 5.6e-04 |
| LSCC | EPCAM | MHC1 | CLDN4+TACSTD2 | 78 | -0.432 | -0.604 | -0.229 | 7.8e-05 | 7.0e-04 |

Contrasts that carry the partial-CLDN4 question:

**LSCC CD8A.** All-tumor TROP2 ρ=-0.103 (n=108, CI -0.288 to +0.085). TROP2 given CLDN4: partial ρ=-0.069 (n=78, CI -0.287 to +0.146, q=0.639). TROP2 given EPCAM: partial ρ=-0.076 (n=108, CI -0.256 to +0.104, q=0.6). TROP2 given CLDN4 and EPCAM: partial ρ=-0.006 (n=78, CI -0.221 to +0.196, q=0.983). CLDN4 given EPCAM: partial ρ=-0.396 (n=78, CI -0.597 to -0.163, q=0.0015). CLDN4 given TROP2: partial ρ=-0.440 (n=78, CI -0.615 to -0.222, q=6.7e-04).

**LSCC MHC1.** All-tumor TROP2 ρ=+0.029 (n=108, CI -0.172 to +0.224). TROP2 given CLDN4: partial ρ=+0.104 (n=78, CI -0.125 to +0.330, q=0.583). TROP2 given EPCAM: partial ρ=+0.067 (n=108, CI -0.128 to +0.249, q=0.614). TROP2 given CLDN4 and EPCAM: partial ρ=+0.187 (n=78, CI -0.048 to +0.402, q=0.241). CLDN4 given EPCAM: partial ρ=-0.356 (n=78, CI -0.556 to -0.135, q=0.00413). CLDN4 given TROP2: partial ρ=-0.417 (n=78, CI -0.591 to -0.207, q=0.00106).

**LUAD CD8A.** All-tumor TROP2 ρ=-0.038 (n=110, CI -0.247 to +0.159). TROP2 given CLDN4: partial ρ=+0.022 (n=79, CI -0.236 to +0.254, q=0.922). TROP2 given EPCAM: partial ρ=+0.013 (n=110, CI -0.170 to +0.200, q=0.948). TROP2 given CLDN4 and EPCAM: partial ρ=+0.043 (n=79, CI -0.188 to +0.271, q=0.792). CLDN4 given EPCAM: partial ρ=-0.075 (n=79, CI -0.324 to +0.168, q=0.614). CLDN4 given TROP2: partial ρ=-0.100 (n=79, CI -0.345 to +0.151, q=0.583).

**LUAD MHC1.** All-tumor TROP2 ρ=-0.083 (n=110, CI -0.276 to +0.119). TROP2 given CLDN4: partial ρ=-0.269 (n=79, CI -0.465 to -0.059, q=0.046). TROP2 given EPCAM: partial ρ=-0.078 (n=110, CI -0.265 to +0.112, q=0.6). TROP2 given CLDN4 and EPCAM: partial ρ=-0.250 (n=79, CI -0.448 to -0.038, q=0.0677). CLDN4 given EPCAM: partial ρ=+0.077 (n=79, CI -0.155 to +0.297, q=0.614). CLDN4 given TROP2: partial ρ=+0.125 (n=79, CI -0.106 to +0.336, q=0.516).

## Sensitivity

Mediator-complete fits drop the requirement that the other epithelial protein was quantified. WES-purity fits add purity rank to every equation on tumors that also have a purity value. q in the purity file is BH across those purity indirect tests only. Neither family replaces the shared-sample primary rule.

Mediator-complete single-mediator models:

| Cohort | Y | Mediator | n | a | b | c | ab | ab CI low | ab CI high | ab p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LUAD | CD8A | CLDN4 | 79 | +0.275 | -0.104 | -0.005 | -0.029 | -0.112 | +0.044 | 0.428 |
| LUAD | CD8A | EPCAM | 110 | +0.372 | -0.140 | -0.038 | -0.052 | -0.127 | +0.034 | 0.208 |
| LUAD | MHC1 | CLDN4 | 79 | +0.275 | +0.126 | -0.245 | +0.035 | -0.030 | +0.113 | 0.293 |
| LUAD | MHC1 | EPCAM | 110 | +0.372 | +0.001 | -0.083 | +3.4e-04 | -0.077 | +0.073 | 0.974 |
| LSCC | CD8A | CLDN4 | 78 | +0.081 | -0.439 | -0.097 | -0.035 | -0.144 | +0.066 | 0.502 |
| LSCC | CD8A | EPCAM | 108 | +0.085 | -0.390 | -0.103 | -0.033 | -0.121 | +0.046 | 0.414 |
| LSCC | MHC1 | CLDN4 | 78 | +0.081 | -0.417 | +0.061 | -0.034 | -0.137 | +0.062 | 0.506 |
| LSCC | MHC1 | EPCAM | 108 | +0.085 | -0.391 | +0.029 | -0.033 | -0.121 | +0.042 | 0.397 |

WES-purity rank-adjusted paths:

| Cohort | Y | Mediator | Model | n | a | b | c | ab | ab CI low | ab CI high | ab p | q |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LUAD | CD8A | CLDN4 | single | 77 | +0.256 | -0.066 | +0.079 | -0.017 | -0.103 | +0.050 | 0.592 | 0.803 |
| LUAD | CD8A | EPCAM | single | 77 | +0.128 | +0.002 | +0.079 | +2.6e-04 | -0.033 | +0.047 | 0.956 | 0.956 |
| LUAD | CD8A | CLDN4 | parallel | 77 | +0.256 | -0.067 | +0.079 | -0.017 | -0.104 | +0.049 | 0.585 | 0.803 |
| LUAD | CD8A | EPCAM | parallel | 77 | +0.128 | +0.010 | +0.079 | +0.001 | -0.032 | +0.050 | 0.918 | 0.956 |
| LUAD | MHC1 | CLDN4 | single | 77 | +0.256 | +0.131 | -0.199 | +0.034 | -0.028 | +0.119 | 0.308 | 0.803 |
| LUAD | MHC1 | EPCAM | single | 77 | +0.128 | -0.057 | -0.199 | -0.007 | -0.055 | +0.043 | 0.779 | 0.89 |
| LUAD | MHC1 | CLDN4 | parallel | 77 | +0.256 | +0.140 | -0.199 | +0.036 | -0.028 | +0.121 | 0.296 | 0.803 |
| LUAD | MHC1 | EPCAM | parallel | 77 | +0.128 | -0.074 | -0.199 | -0.010 | -0.057 | +0.040 | 0.727 | 0.89 |
| LSCC | CD8A | CLDN4 | single | 77 | +0.063 | -0.354 | -0.075 | -0.022 | -0.112 | +0.056 | 0.601 | 0.803 |
| LSCC | CD8A | EPCAM | single | 77 | +0.139 | -0.251 | -0.075 | -0.035 | -0.112 | +0.017 | 0.221 | 0.803 |
| LSCC | CD8A | CLDN4 | parallel | 77 | +0.063 | -0.348 | -0.075 | -0.022 | -0.109 | +0.054 | 0.602 | 0.803 |
| LSCC | CD8A | EPCAM | parallel | 77 | +0.139 | -0.237 | -0.075 | -0.033 | -0.105 | +0.013 | 0.197 | 0.803 |
| LSCC | MHC1 | CLDN4 | single | 77 | +0.063 | -0.296 | +0.093 | -0.019 | -0.086 | +0.051 | 0.558 | 0.803 |
| LSCC | MHC1 | EPCAM | single | 77 | +0.139 | -0.209 | +0.093 | -0.029 | -0.104 | +0.015 | 0.222 | 0.803 |
| LSCC | MHC1 | CLDN4 | parallel | 77 | +0.063 | -0.291 | +0.093 | -0.018 | -0.084 | +0.050 | 0.555 | 0.803 |
| LSCC | MHC1 | EPCAM | parallel | 77 | +0.139 | -0.197 | +0.093 | -0.027 | -0.095 | +0.011 | 0.198 | 0.803 |

## What this page does not do

It does not impute CLDN4. It does not pool LUAD with LSCC. It does not substitute RNA for a missing protein. It does not use ImmuneScore, GEP, or CIBERSORT as the outcome; those are already on other pages. It does not restate private PDX numbers. It does not turn a cross-sectional product of rank coefficients into a causal demonstration that TROP2 changes CD8 or MHC through CLDN4.

Figures: `figures/fig1_missingness.png`, `figures/fig2_indirect_forest.png`, `figures/fig3_partial_forest.png`, `figures/fig4_cd8a_scatter.png`.

