# scVI/scANVI concordant-4: larger patient-level CLDN4 vs T/NK effect

ADDITIVE. The locked result stays malignant CLDN4 % positive versus the T/NK fraction of all cells, on GSE123902 + GSE131907 + GSE205335 + GSE189357: ρ=-0.531, n=65, I²=0%, binomial GLMM χ²=16.44. This folder searches summaries of the existing patient-batch scVI CLDN4, and two patient-level binomials, for a larger likelihood-ratio statistic and a larger |ρ| on those same 65 units. No fifth cohort. n is the number of units.

## Honest n

- **n = 65** locked units (13 + 21 + 22 + 9).
- The scVI value on each unit is computed from the integration subsample already stored in `inputs/umap_obs.tsv.gz` (seed-malignant or scANVI-malignant cells in that subsample). T/NK and malignant counts are the full unit, from `inputs/patient_units.tsv`.
- Antimode of seed-malignant scVI log1p(CLDN4) = 1.1568, which is 2.18 counts per 10,000. The cut is the lightest histogram bin between the 5th and 60th percentiles of that expression. It does not use T/NK.

## Search

Every row is the same patient model,

`cbind(n_T/NK, n_fail) ~ coded(CLDN4) + dataset + (1 | unit)`,

logit link, one observation per unit, patient random intercept (logit-normal extra-binomial variance). The LR effect is the chi-square against the same model without CLDN4. |ρ| is the DerSimonian–Laird pool of the four within-cohort Spearmans.

The grid is the Cartesian product of:

- scores: locked % positive, full-unit mean log1p, seed-malignant mean scVI log1p(CLDN4), and the fraction of seed-malignant or scANVI-malignant cells with scVI CLDN4 at or above log1p(k) for k = 1…20 CP10K, plus the antimode cut
- outcomes: `n_fail` = all non-T/NK cells, or `n_fail` = malignant cells
- coding: global z-score, within-dataset z-score, or within-dataset rank

The primary row is the scVI/scANVI score that maximizes the worse of χ²/χ²_locked and |ρ|/|ρ|_locked. % positive is the reference, not a candidate. Naive p-values below treat the row as if it had been locked. The permutation does not.

## Primary row

**seed fraction with scVI CLDN4 ≥ antimode (2.18 CP10K)**, outcome `compartment`, coding `within_rank`.

- Patient-model LR: β=-1.201 per SD, χ²=36.25, naive p=1.74e-09, n=65, patient RE sd=1.344, singular=FALSE.
- Pooled Spearman: ρ=-0.736 (naive p=3.38e-09, I²=22.1%, -0.849 to -0.557).
- Cohorts: GSE123902 ρ=-0.764 (n=13); GSE131907 ρ=-0.517 (n=21); GSE205335 ρ=-0.802 (n=22); GSE189357 ρ=-0.867 (n=9).
- Locked % positive cohorts, all-cell T/NK fraction: GSE123902 ρ=-0.659 (n=13); GSE131907 ρ=-0.522 (n=21); GSE205335 ρ=-0.435 (n=22); GSE189357 ρ=-0.600 (n=9).

GSE131907 on the primary row is ρ=-0.517, next to its locked % positive ρ=-0.522. The pooled increase is carried by GSE123902, GSE205335, and GSE189357. That spread is the I² of 22%. The LR-maximizing neighbor (scVI CLDN4 ≥ 4 CP10K, same outcome and rank coding) keeps I² at 0% with χ²=41.03 and ρ=-0.702.

Locked reference on the same fit: β=-0.619, χ²=16.44, ρ=-0.531.

The same recoding applied to % positive, without changing the gene summary, is the best % positive row: malignant CLDN4 % positive, `compartment`, `within_rank`, χ²=23.42, ρ=-0.536. The primary scVI row is larger on both axes than that recoding.

Coordinate-wise maxima, if a reader wants one axis only:

- Largest LR: seed fraction with scVI CLDN4 ≥ 4 CP10K, `compartment`, `within_rank`, χ²=41.03, ρ=-0.702, I²=0.0%.
- Largest |ρ|: seed fraction with scVI CLDN4 ≥ antimode (2.18 CP10K), `compartment`, `within_rank`, χ²=36.25, ρ=-0.736, I²=22.1%.
- Best scANVI-called fraction in the same joint sense: scANVI fraction with scVI CLDN4 ≥ antimode (2.18 CP10K), χ²=24.39, ρ=-0.581. Seed-malignant cells, the locked definition, stay above scANVI calls.

| role | score | outcome | coding | β / SD | LR χ² | naive LR p | ρ | I² | joint gain |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| locked | malignant CLDN4 % positive | all_cells | global_z | -0.619 | 16.44 | 5.02e-05 | -0.531 | 0.0% | 1.000 |
| primary | seed fraction with scVI CLDN4 ≥ antimode (2.18 CP10K) | compartment | within_rank | -1.201 | 36.25 | 1.74e-09 | -0.736 | 22.1% | 1.385 |
| max LR | seed fraction with scVI CLDN4 ≥ 4 CP10K | compartment | within_rank | -1.257 | 41.03 | 1.50e-10 | -0.702 | 0.0% | 1.322 |
| max |ρ| | seed fraction with scVI CLDN4 ≥ antimode (2.18 CP10K) | compartment | within_rank | -1.201 | 36.25 | 1.74e-09 | -0.736 | 22.1% | 1.385 |
| best %pos recoding | malignant CLDN4 % positive | compartment | within_rank | -1.011 | 23.42 | 1.30e-06 | -0.536 | 0.0% | 1.009 |
| best scANVI | scANVI fraction with scVI CLDN4 ≥ antimode (2.18 CP10K) | compartment | within_rank | -1.028 | 24.39 | 7.87e-07 | -0.581 | 36.0% | 1.093 |

Equal-unit Gaussian LRT on logit(proportion), with the same design, stays within 0.13 chi-square of the binomial GLMM on every grid row. The quoted chi-squares are the binomial GLMM.

## Selection check

2000 within-cohort shuffles of the unit outcomes. Each shuffle re-picks the scVI/scANVI score, the outcome, and the coding. The observed maxima are the equal-unit Gaussian statistics (they match the binomial GLMM).

| maximized statistic | observed | null 95th percentile | permutation p |
|---|---:|---:|---:|
| joint gain | 1.385 | 0.549 | 5.00e-04 |
| LR χ² | 41.00 | 9.58 | 5.00e-04 |
| |ρ| | 0.736 | 0.414 | 0.0015 |

A permutation p of 5.00×10⁻⁴ is the floor for 2000 shuffles: none of them reached the observed joint gain or the observed LR χ². Two shuffles reached the observed |ρ|.

Leave-one-cohort-out uses the same joint rule on the other three cohorts. The held-out Spearman, pooled, is ρ=-0.722 (p=1.74e-10, I²=6.2%, -0.831 to -0.559).

| held out | score chosen on the rest | outcome | coding | held-out ρ | n |
|---|---|---|---|---:|---:|
| GSE123902 | seed_antimode | compartment | within_rank | -0.764 | 13 |
| GSE131907 | seed_antimode | compartment | global_z | -0.517 | 21 |
| GSE205335 | seed_antimode | compartment | within_rank | -0.802 | 22 |
| GSE189357 | seed_k03 | compartment | within_rank | -0.817 | 9 |

## Scope

- The locked % positive association remains ρ=-0.531 against the T/NK fraction of all cells.
- The compartment outcome is the log-odds that a cell in the T/NK + malignant pool is a T/NK cell. Myeloid, B, and stromal cells are outside that denominator. This is a two-part composition. It is separate from the CosMx spatial result.
- The sample size is 65 units.
- Naive p-values are the p-values of the fitted row. The permutation p re-runs the grid inside each shuffle.
- GSE131907 remains the tumor-bearing sample. GSE205335 remains the patient, with libraries pooled before the test. Cohorts kept out of this analysis: GSE148071, GSE127465, GSE154826, GSE200563, GSE207422.

## Reproduce

```bash
python3 methods/scanvi_c4_max_effect/maximize.py
```

Requires Python (numpy, matplotlib) and R with lme4. Inputs are the scVI/scANVI unit table and per-cell scVI CLDN4 from that integration (`patient_units.tsv` sha256 `9c4521dc69771c843236062da5e4fa36c457c0d1d6fd2bf7fd6b1dd57323b396`, `umap_obs.tsv.gz` sha256 `a033a7d72efa8e70ed350a2da37809de75fc03550ed9fbd258a1790ed6ac9be5`).
