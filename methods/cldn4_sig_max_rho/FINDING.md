# Maximum |ρ|: locked 221-gene CLDN4-high signature vs CD8A / ImmuneScore

The 221 genes and their order are the PR 590 signature (`methods/cldn4_high_malignant_signature/tables/signature_genes.tsv`). Prefixes use that order. No gene was added or removed because of its correlation with CD8A or ImmuneScore. GSE10072, GSE11969, and GSE248378 stay closed. Bulk ρ is not a spatial exclusion result.

## Objective

Primary number: absolute DerSimonian–Laird meta-analytic Spearman versus CD8A on all QC-passing tumors. Studies are OncoSG, GSE273377 (discovery and validation inverse-variance combined first), GSE282774, and GSE233774 tumors. A spec is eligible only when all five strata have a finite correlation, so a cohort cannot be dropped to raise the meta |ρ|. The search is method (z-mean or ssGSEA), ssGSEA α in {0, 0.25, 0.75, 1}, prefix size 5 through 221, and purity mode {unadjusted, partial correlation on published PURITY or ESTIMATE StromalScore, gene-level residual on that same covariate}. The high-tumor median split is not in this objective.

The selected spec was chosen on these same cohorts, so its meta p-value is the p-value of a maximized |ρ|, not a single pre-specified test. The pre-specified 221-gene z-mean from PR 590 is the baseline row. The OncoSG correlation at the GEO-only winning spec was not used to pick that spec.

## Primary maximum

Maximum |meta ρ| versus CD8A is **-0.533** (fixed, p=3.05e-33, I²=0%, n sum=420, k=4). Spec: **ssGSEA, α=0.75, size 163, unadjusted**.

The locked baseline (z-mean, size 221, unadjusted) meta ρ is -0.409 (DL random, p=1.82e-05, I²=71%).

Per cohort, CD8A, at the winning spec:

OncoSG -0.567 (p=8.70e-16, n=169); GSE273377 discovery -0.519 (p=1.88e-08, n=103); GSE273377 validation -0.459 (p=0.000222, n=60); GSE282774 -0.560 (p=4.82e-06, n=58); GSE233774 tumor -0.461 (p=0.0104, n=30)

ImmuneScore at that same spec (not re-selected):

OncoSG -0.600 (p=7.08e-18, n=169); GSE273377 discovery -0.557 (p=9.59e-10, n=103); GSE273377 validation -0.482 (p=9.71e-05, n=60); GSE282774 -0.604 (p=5.23e-07, n=58); GSE233774 tumor -0.601 (p=0.00044, n=30)

ImmuneScore meta ρ at the CD8A spec: -0.574 (fixed, p=7.76e-40, I²=0%).

The eight largest eligible CD8A meta |ρ| values are ssGSEA α=0.75 size 163 none -0.533, ssGSEA α=0.75 size 156 none -0.533, ssGSEA α=0.75 size 157 none -0.533, ssGSEA α=0.75 size 159 none -0.532, ssGSEA α=1 size 163 none -0.532, ssGSEA α=0.75 size 162 none -0.532, ssGSEA α=0.75 size 155 none -0.532, ssGSEA α=0.75 size 160 none -0.531. The maximum is not a one-prefix spike.

## Purity at the winning score

Unadjusted correlation had the largest |meta ρ|. Partial correlation uses published PURITY on OncoSG and ESTIMATE StromalScore on GEO. Residual mode builds the score from gene-level linear residuals on that covariate and then correlates it with the raw endpoint.

unadjusted: meta ρ -0.533 (fixed, p=3.05e-33, I²=0%). OncoSG -0.567 (p=8.70e-16, n=169); GSE273377 discovery -0.519 (p=1.88e-08, n=103); GSE273377 validation -0.459 (p=0.000222, n=60); GSE282774 -0.560 (p=4.82e-06, n=58); GSE233774 tumor -0.461 (p=0.0104, n=30)

partial | purity/stroma: meta ρ -0.416 (fixed, p=5.08e-19, I²=0%). OncoSG -0.390 (p=1.76e-07, n=169); GSE273377 discovery -0.407 (p=2.17e-05, n=103); GSE273377 validation -0.466 (p=0.000202, n=60); GSE282774 -0.520 (p=3.44e-05, n=58); GSE233774 tumor -0.271 (p=0.156, n=30)

gene residual | purity/stroma: meta ρ -0.327 (fixed, p=6.82e-12, I²=0%). OncoSG -0.369 (p=8.02e-07, n=169); GSE273377 discovery -0.243 (p=0.0136, n=103); GSE273377 validation -0.432 (p=0.000572, n=60); GSE282774 -0.268 (p=0.0424, n=58); GSE233774 tumor -0.194 (p=0.305, n=30)


## ImmuneScore maximum

Repeating the objective for ImmuneScore gives meta ρ **-0.574** (fixed, p=7.76e-40, I²=0%, n sum=420). Spec: **ssGSEA, α=0.75, size 163, unadjusted**.

OncoSG -0.600 (p=7.08e-18, n=169); GSE273377 discovery -0.557 (p=9.59e-10, n=103); GSE273377 validation -0.482 (p=9.71e-05, n=60); GSE282774 -0.604 (p=5.23e-07, n=58); GSE233774 tumor -0.601 (p=0.00044, n=30)

## Held-out OncoSG

The GEO-only maximum (OncoSG not in the selection) is ssGSEA, α=1, size 109, unadjusted, GEO meta ρ -0.520 (fixed, p=2.94e-19, I²=0%, n sum=251).

Applied to OncoSG versus CD8A: ρ=-0.516 (p=7.23e-13, n=169). Versus ImmuneScore: ρ=-0.559 (p=2.67e-15, n=169).

## Largest |ρ| inside one cohort

Inside the primary search space (all tumors), the largest |ρ| is -0.690 in GSE282774 versus ImmuneScore (p=2.02e-09, n=58). Spec: ssGSEA α=0, size 22, purity=none. That cell was not required to be the same spec in the other cohorts.

Largest |ρ| for CD8A inside each cohort, same rule:

| cohort | n | ρ | p | spec |
|---|---:|---:|---:|---|
| OncoSG | 169 | -0.593 | 2.13e-17 | ssGSEA, α=0, size 164, none |
| GSE273377 discovery | 103 | -0.569 | 3.65e-10 | ssGSEA, α=1, size 207, none |
| GSE273377 validation | 60 | -0.494 | 6.95e-05 | ssGSEA, α=1, size 114, partial |
| GSE282774 | 58 | -0.641 | 5.95e-08 | ssGSEA, α=0, size 16, none |
| GSE233774 tumor | 30 | -0.558 | 0.00135 | ssGSEA, α=0.75, size 30, none |

## High-tumor subset

Not the primary objective. OncoSG keeps samples with published PURITY at or above the cohort median. GEO keeps samples with ESTIMATE StromalScore at or below the cohort median. z-means are recomputed inside the subset.

Maximum |meta ρ| versus CD8A on this subset is -0.450 (DL random, p=3.44e-09, I²=18%, n sum=213). Spec: ssGSEA, α=1, size 163, unadjusted.

OncoSG -0.351 (p=0.000858, n=87); GSE273377 discovery -0.562 (p=1.46e-05, n=52); GSE273377 validation -0.439 (p=0.0152, n=30); GSE282774 -0.595 (p=0.00066, n=29); GSE233774 tumor -0.225 (p=0.42, n=15)

## Reproduction of PR 590 cells

All 15 locked cells matched within the printed tolerance, including OncoSG z-mean size 221 versus CD8A at the full stored precision of the PR 590 summary.

Purity covariate: OncoSG uses the cBioPortal published PURITY column. GEO uses ESTIMATE stromal ssGSEA (α=0.25, package 1.0.13 stromal set, common-gene filter). The Affymetrix cosine purity formula is not used. Partial correlation is Pearson of rank residuals. Residual mode is an ordinary Spearman of a score built from gene-level linear residuals on that covariate; the covariate is not the immune endpoint.

Full grid: `tables/cohort_grid.tsv`. One row per spec: `tables/meta_grid.tsv`. Locked checks: `tables/locked_checks.tsv`.
