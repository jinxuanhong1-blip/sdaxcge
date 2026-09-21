# TACSTD2 → T/NK on concordant-4, and how much of it runs through CLDN4

ADDITIVE mechanism test. Locked concordant-4 only: GSE123902 (13 donors) + GSE131907 (21 samples) + GSE205335 (22 patients) + GSE189357 (9 patients). N = 65. Not GSE148071, GSE127465, GSE207422, GSE154826, GSE200563, or E-MTAB-13526.

The unit is the patient, donor, or sample already locked for malignant CLDN4 % positive versus T/NK. Cell counts are not n. p-values are descriptive. This is not a causal mediation, not a spatial exclusion result, and not a claim about private KL tumors.

## Pipeline check

Malignant CLDN4 fraction with UMI > 0 versus T/NK fraction, DerSimonian–Laird ρ = -0.531168 (p = 1.65e-05, I² = 0.0%). The locked published value is −0.5311678045689989. Cohort ρ: GSE123902 -0.659 (n=13), GSE131907 -0.522 (n=21), GSE205335 -0.435 (n=22), GSE189357 -0.600 (n=9).

## What was maximized

Pre-specified grid. Same malignant-cell summary for TACSTD2 and CLDN4, plus a cross of each gene against the locked UMI>0 fraction of the other. Summaries: UMI>0, ≥2, ≥3, ≥5, ≥10; CP10k ≥1, ≥5, ≥10; fraction above the cohort median UMI and above the cohort 75th percentile; mean log1p(UMI); mean log1p(CP10k); mean log1p among detected cells; 90th percentile on both scales; pseudobulk log1p(CP10k).

Models, each with a cohort fixed effect: raw scores; within-cohort z-scores of TACSTD2, CLDN4, and T/NK; within-cohort ranks; logit of proportions (clip 10⁻⁴) when both scores are fractions. Mediation is Baron–Kenny on that linear model. The product of coefficients equals the drop in the TACSTD2 coefficient. Proportion = a×b / c.

A row is tier A when all four of these hold: every cohort has TACSTD2–T/NK ρ < 0, CLDN4–T/NK ρ < 0, and TACSTD2–CLDN4 ρ > 0 with |ρ| < 0.95; the pooled Spearman partial of TACSTD2 also moves toward zero and stays negative; the pooled OLS path is a > 0, b < 0, c < 0, c′ < 0, and 0 < a×b/c < 1; and the same partial-attenuation pattern holds in the separate OLS of each cohort. Tier B drops only the per-cohort OLS requirement. Tier C drops the Spearman-attenuation requirement and keeps the sign pattern plus pooled OLS attenuation. The reported maximum is the largest mediation proportion inside the best non-empty tier, among rows whose pooled TACSTD2–T/NK ρ is ≤ −0.20. Ties break on lower I², then on |ρ| × proportion.

Grid size 89. Tier counts: none = 89.

No row cleared tier C with pooled TACSTD2–T/NK ρ ≤ −0.20. Every TACSTD2 summary in the grid is positively correlated with T/NK in GSE205335, so none is 4-cohort negative. The most negative pooled TACSTD2–T/NK ρ in the grid is the locked UMI>0 fraction.

Largest same-sign attenuation in the grid, without the 4-cohort sign rule: `pct_gt0|pb_log1p_cp10k|identity`. Proportion 0.328, total β -0.0871, direct β -0.0586, pooled TACSTD2–T/NK ρ -0.112, partial ρ -0.046. Cohorts with the per-cohort OLS attenuation pattern: 1/4. This row stays in the table as the grid maximum under a weaker rule. It is not the 4-cohort result.

The six tests below are the locked UMI>0 specification, which is the pre-specified score.

## Locked UMI > 0

Pre-specified row.

Specification `pct_gt0|pct_gt0|identity` (paired, tier none). N = 65.

### (1) Spearman

TACSTD2 vs T/NK: ρ = -0.112 (p = 0.4603, I² = 15.4%, -0.388 to 0.183).
CLDN4 vs T/NK: ρ = -0.531 (p = 1.65e-05, I² = 0.0%).
TACSTD2 vs CLDN4: ρ = 0.535 (p = 0.0058, I² = 56.3%).

### (2) Partial Spearman

TACSTD2–T/NK | CLDN4: ρ = 0.216 (p = 0.1239, I² = 0.0%, -0.060 to 0.462).
CLDN4–T/NK | TACSTD2: ρ = -0.543 (p = 2.02e-05, I² = 0.0%).
The partial Spearman changes sign relative to the marginal Spearman, so 1 − ρ_partial/ρ_total is not an attenuation fraction in (0, 1).

### (3) Nested OLS with cohort fixed effects

Total TACSTD2 β = -0.0871 (SE 0.1160, p = 0.4555).
Direct TACSTD2 β after CLDN4 = 0.1659 (SE 0.1225, p = 0.1808).
The TACSTD2 coefficient changes sign, from -0.0871 to 0.1659. The ratio a×b/c equals 2.904, which is outside (0, 1).

### (4) Mediation TACSTD2 → CLDN4 → T/NK

Baron–Kenny / product of coefficients, same cohort fixed effects. a (TACSTD2 → CLDN4) = 0.4945 (p = 1.12e-05). b (CLDN4 → T/NK | TACSTD2) = -0.5116 (p = 2.29e-04). Indirect a×b = -0.2530. Proportion mediated a×b/c = 2.904. Sobel p = 0.0024.
The OLS identity (c − c') − a×b has absolute gap 1.11e-16.
Cohort-stratified bootstrap, 5000 draws, 5000 finite. Proportion median 1.603, 95% percentile interval -20.572 to 24.473. Indirect 95% interval -0.4362 to -0.0954. Share of draws with proportion in (0, 1): 0.084. Share with the same partial-attenuation sign pattern: 0.083.
Leave-one-cohort-out proportion: drop GSE123902 → 10.384 (n=52), drop GSE131907 → 6.074 (n=44), drop GSE205335 → 1.270 (n=43), drop GSE189357 → 3.998 (n=56).

### (5) TACSTD2 × CLDN4 interaction

Product-term β = 0.0734 (SE 0.4763, p = 0.8781). Simple slope of TACSTD2 at the pooled 25th percentile of CLDN4 = 0.1465; at the 75th percentile = 0.1740.

### (6) CLDN4 Q1 vs Q4

| split | level | n | cohorts in Spearman meta | ρ | p | I² | OLS β | p |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| quartile | Q1 | 19 | 2 | 0.063 | 0.8764 | 0.0% | -0.1263 | 0.6186 |
| quartile | Q4 | 16 | 2 | 0.702 | 0.1432 | 42.0% | 0.2399 | 0.3024 |
| median | low | 34 | 4 | -0.292 | 0.1589 | 0.0% | -0.2069 | 0.1739 |
| median | high | 31 | 3 | 0.402 | 0.2806 | 61.3% | 0.3148 | 0.0923 |

Quartile Spearman meta uses cohorts with n ≥ 5 inside the bin. GSE189357 has 9 patients, so a quartile bin often falls under that floor.
CLDN4 low half, cohort Spearman of the TACSTD2 score versus T/NK: GSE123902 -0.429 (n=7), GSE131907 -0.136 (n=11), GSE205335 -0.409 (n=11), GSE189357 -0.100 (n=5).
CLDN4 high half, cohort Spearman of the TACSTD2 score versus T/NK: GSE123902 -0.543 (n=6), GSE131907 +0.552 (n=10), GSE205335 +0.727 (n=11), GSE189357 -0.800 (n=4).

### GLMM random intercept for cohort

Total β = -0.0871 (p = 0.4343). Direct β = 0.1659 (p = 0.1552). a = 0.4945, b = -0.5116. Proportion = 2.904. Random-intercept variance: base 0.0000, full 0.0000. Converged: base True, full True. A random-intercept variance of 0 leaves the coefficients equal to the cohort-fixed-effect OLS coefficients on this specification.

Cohort detail:

| cohort | n | TACSTD2–T/NK ρ | CLDN4–T/NK ρ | TACSTD2–CLDN4 ρ | partial TACSTD2 | partial CLDN4 | cohort OLS proportion |
|---|---:|---:|---:|---:|---:|---:|---:|
| GSE123902 | 13 | -0.324 | -0.659 | 0.302 | -0.174 | -0.623 | 0.663 |
| GSE131907 | 21 | -0.151 | -0.522 | 0.671 | 0.316 | -0.575 | 1.825 |
| GSE205335 | 22 | 0.208 | -0.435 | 0.217 | 0.345 | -0.503 | -0.058 |
| GSE189357 | 9 | -0.483 | -0.600 | 0.850 | 0.063 | -0.410 | 1.338 |

## Top eligible rows

No eligible row.

## Reading

The stable path is the product of coefficients. On the locked UMI>0 scores, a is positive and b is negative, so a×b is negative: higher TACSTD2 tracks higher CLDN4, and higher CLDN4 tracks lower T/NK. The cohort-stratified bootstrap interval for a×b stays below zero. Sobel p for that product is 0.0024. The marginal TACSTD2 coefficient is a small negative number, and the direct coefficient after CLDN4 is positive, so a×b/c sits above 1. That ratio is outside a partial-mediation proportion. 8.4% of bootstrap draws fall inside (0, 1).

CLDN4’s partial Spearman with T/NK, holding TACSTD2, stays near the marginal CLDN4 result and is negative in every cohort, with I² = 0. TACSTD2’s partial Spearman, holding CLDN4, is positive in the pool. GSE123902 is the only cohort whose own OLS shows a negative TACSTD2 coefficient that shrinks and stays negative (proportion 0.66).

No TACSTD2 summary in the 89-row grid is negatively correlated with T/NK in GSE205335. The most negative pooled ρ is −0.112, the locked percent. The largest same-sign attenuation anywhere in the grid is 0.33, for TACSTD2 percent versus T/NK with CLDN4 entered as malignant pseudobulk log1p(CP10k). That row’s pooled ρ is still −0.112, and the within-cohort OLS attenuation pattern is present in 1 of 4 cohorts.

The percent split points the same way. CLDN4-only (CLDN4 detected, TACSTD2 undetected) is negative versus T/NK in all four cohorts (I² = 0). Double-positive cells are negative in three cohorts and positive in GSE205335. TACSTD2-only cells are negative in one cohort. In the joint OLS all three compartment coefficients are negative; the TACSTD2-only coefficient is the least precise because that compartment is the smallest fraction.

Below each cohort’s CLDN4 median, TACSTD2 percent versus T/NK is negative in all four cohorts (pooled ρ −0.29, I² = 0, n = 34). Above the median the cohorts disagree and the pool is positive. The continuous TACSTD2 × CLDN4 product term is close to zero.

Taken together, the immune-cold association that is consistent across the four cohorts is CLDN4’s. TACSTD2 shares that direction through its correlation with CLDN4 (the indirect path). A residual TACSTD2 association in the same direction, after CLDN4, is not a 4-cohort result.

Cohort fixed effects absorb the mean difference between studies. They do not make the four studies exchangeable. GSE123902 and GSE189357 malignant cells are marker-gated (EPCAM or KRT8/18/19, and PTPRC = 0). GSE131907 and GSE205335 use the author malignant label. T/NK is the locked fraction: marker T/NK and not malignant, or the author T/NK label, over all cells in the unit.

A high mediation proportion on a weak TACSTD2–T/NK correlation is a large share of a small association. The table carries both numbers. The searched maximum is not a confirmatory p-value.

## Compartment split of the UMI > 0 fractions

On malignant cells, TACSTD2 % = (TACSTD2>0 and CLDN4>0) + (TACSTD2>0 and CLDN4=0). CLDN4 % = double-positive + CLDN4-only. The two identities hold to numerical error on all 65 units. This split asks which part of the TACSTD2-positive fraction carries the T/NK association.

| score | ρ | p | I² | cohorts < 0 | GSE123902 | GSE131907 | GSE205335 | GSE189357 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| pct_both_gt0 | -0.144 | 0.4542 | 45.8% | 3/4 | -0.352 | -0.248 | 0.307 | -0.483 |
| pct_tacstd2_only_gt0 | 0.154 | 0.2593 | 0.0% | 1/4 | -0.154 | 0.248 | 0.100 | 0.500 |
| pct_cldn4_only_gt0 | -0.350 | 0.0077 | 0.0% | 4/4 | -0.385 | -0.284 | -0.477 | -0.033 |
| tacstd2_pct_gt0 | -0.112 | 0.4603 | 15.4% | 3/4 | -0.324 | -0.151 | 0.208 | -0.483 |
| cldn4_pct_gt0 | -0.531 | 1.65e-05 | 0.0% | 4/4 | -0.659 | -0.522 | -0.435 | -0.600 |

Mean fractions across the 65 units: double-positive 0.449, TACSTD2-only 0.081, CLDN4-only 0.153.
OLS of T/NK fraction on the three compartments with cohort fixed effects (n = 65, R² = 0.392). Double-positive β = -0.3734 (p = 0.0024). TACSTD2-only β = -0.9102 (p = 0.0305). CLDN4-only β = -0.7771 (p = 6.95e-06). β is the change in T/NK fraction per 1.0 change in that compartment, so the TACSTD2-only slope applies to a compartment whose mean is 0.081.

## Caveats

Observational co-expression. CLDN4 may sit beside TACSTD2 on a shared epithelial program rather than on a path from TACSTD2 to immune composition. No intervention, no instrument, no spatial radius. Keratin adjustment is not in the grid; the locked CLDN4 result was already strongest without it. Logit uses a 10⁻⁴ clip at 0 and 1. Within-cohort quartile bins in GSE189357 are small. Four random-effect levels is a thin GLMM.

Reproduce: `bash methods/concordant4_tacstd2_cldn4_mediation/scripts/download.sh /tmp/geo_dl` then `python3 methods/concordant4_tacstd2_cldn4_mediation/scripts/extract_scores.py` then `python3 methods/concordant4_tacstd2_cldn4_mediation/scripts/analyze.py`. GEO matrices stay outside the repo. `results/unit_scores.tsv` is enough to rerun the grid.
