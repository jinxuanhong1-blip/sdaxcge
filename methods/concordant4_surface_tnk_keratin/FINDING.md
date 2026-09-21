# Concordant-4: CLDN1 vs CLDN4 vs CLDN7 vs EPCAM against patient T/NK

Additive. Does not replace the locked CLDN4 result (malignant % positive vs T/NK,
DL ρ = −0.531, p = 1.65×10⁻⁵, I² = 0%, N = 65). Same four cohorts, same units,
same malignant gate, same T/NK fraction. No GSE148071, GSE127465, GSE154826,
GSE200563, or E-MTAB-13526. No mouse Harmony. Cell counts are not n.

Question: after partialling malignant keratin, does public scRNA single out CLDN4
from CLDN7 (and from CLDN1 / EPCAM) the way a CLDN4-specific barrier claim would
require? The locked bulk surface ranking does not. LUAD placed CLDN4 7th (+0.122),
with CLDN7, EPCAM, and MUC1 stronger. That bulk number is a different estimand.
This file asks the specificity question on patient T/NK.

## Gate and scores

- GSE123902 donor (n = 13), GSE189357 patient (n = 9):
  malignant = (EPCAM or KRT8 or KRT18 or KRT19) > 0 and PTPRC = 0.
  T/NK = (CD3D or CD3E or CD8A or NKG7 or GNLY or KLRD1) > 0 and not malignant.
  `frac_tnk = n_tnk / n_cells`.
- GSE131907 sample (n = 21): author malignant subtype, tumor-bearing
  origins tLung / tL/B / mLN / PE / mBrain, n_malignant ≥ 20.
  T/NK = author T lymphocytes or NK cells.
- GSE205335 patient (n = 22): author lineage. Tissue labels that
  start with Normal are out of the patient denominator. The four normal-only donors
  are out. P4001 stays in this correlation (n_malignant = 27).
- Primary malignant score = % of malignant cells with count > 0.
  Secondary = mean log1p(count).
- Keratin covariate = mean, over malignant cells, of the per-cell mean of
  log1p(KRT8), log1p(KRT18), log1p(KRT19). Sensitivity covariate drops KRT8
  and uses only KRT18 and KRT19.
- Partial Spearman residualizes the ranks of both the surface score and
  `frac_tnk` on the keratin ranks, inside each cohort. Pooling is
  DerSimonian–Laird on Fisher z. For one covariate the z variance is 1/(n−4),
  not 1/(n−3).
- Specificity test (pre-specified): on within-cohort rank-z scores, pooled
  partial correlation, 10000 within-cohort shuffles of `frac_tnk` (seed 1).
  Δ = ρ(CLDN4) − ρ(CLDN7). CLDN4 is pinned over CLDN7 only if its keratin-partial
  DL ρ is more negative, the permutation p for Δ is < 0.05 with Δ < 0, and the
  unique partial of CLDN4 given keratin and CLDN7 stays negative with permutation
  p < 0.05. Otherwise the public scRNA does not pin it.
- MUC1 is context. It is in the bulk sentence and in the extraction. It is not
  a fifth member of the pin rule.
- p-values are descriptive. The permutation p is the one used for the pin call.

Honest n = **65** (13 + 21 + 22 + 9).

## Calibration against PR #539

Recomputed malignant CLDN4 % positive matches the locked patient table:
max absolute difference 0.0000 percentage points. Unadjusted DL
ρ = -0.531 (p = 1.65e-05, I² = 0.0%,
N = 65). Locked reference ρ = −0.531.

| cohort | n | unadjusted ρ |
|---|---:|---:|
| GSE123902 | 13 | -0.659 |
| GSE131907 | 21 | -0.522 |
| GSE205335 | 22 | -0.435 |
| GSE189357 | 9 | -0.600 |

## 1. Keratin-partial DL Spearman, malignant % positive vs T/NK

Filled points in `results/figures/forest_partial_rho.png` are these partial ρ
values. Open points are the unadjusted ρ.

| gene | k | N | partial ρ | p | I² | 95% CI |
|---|---:|---:|---:|---:|---:|---|
| CLDN1 | 4 | 65 | 0.225 | 0.1087 | 0.0% | -0.051 to 0.469 |
| CLDN4 | 4 | 65 | -0.478 | 0.0002717 | 0.0% | -0.664 to -0.236 |
| CLDN7 | 4 | 65 | -0.352 | 0.1016 | 55.6% | -0.669 to 0.072 |
| EPCAM | 4 | 65 | -0.368 | 0.00688 | 0.0% | -0.582 to -0.106 |
| MUC1 (context) | 4 | 65 | -0.300 | 0.2131 | 63.5% | -0.662 to 0.176 |

Rank of the partial ρ, most negative first (MUC1 included): CLDN4 -0.478, EPCAM -0.368, CLDN7 -0.352, MUC1 -0.300, CLDN1 0.225.
Strongest exclusion associate on this table: **CLDN4**.
Among CLDN1 / CLDN4 / CLDN7 / EPCAM only, most negative first: CLDN4 -0.478, EPCAM -0.368, CLDN7 -0.352, CLDN1 0.225.
On the point estimate, CLDN4 is first of those four.

Unadjusted DL, same score, for scale:

| gene | k | N | ρ | p | I² | 95% CI |
|---|---:|---:|---:|---:|---:|---|
| CLDN1 | 4 | 65 | 0.162 | 0.2344 | 0.0% | -0.105 to 0.407 |
| CLDN4 | 4 | 65 | -0.531 | 1.65e-05 | 0.0% | -0.697 to -0.312 |
| CLDN7 | 4 | 65 | -0.409 | 0.08864 | 68.4% | -0.732 to 0.066 |
| EPCAM | 4 | 65 | -0.434 | 0.001314 | 8.1% | -0.634 to -0.179 |
| MUC1 (context) | 4 | 65 | -0.347 | 0.2043 | 74.8% | -0.726 to 0.195 |

CLDN4 % positive spans 0.3–96.3 (median 61.7).
CLDN7 spans 13.0–91.8 (median 51.1).
EPCAM spans 8.3–95.5 (median 66.9).
None of the four is a constant, so the partials are not a ceiling artifact.
The marker gate does use EPCAM and the keratins, so EPCAM and keratin variation
in GSE123902 and GSE189357 is variation inside an already epithelial gate.
That is why the author-label subset is reported below.

Cohort partial ρ for CLDN4 and CLDN7 (% positive, keratin):

| cohort | n | CLDN4 ρ | CLDN4 p | CLDN7 ρ | CLDN7 p |
|---|---:|---:|---:|---:|---:|
| GSE123902 | 13 | -0.292 | 0.3574 | -0.409 | 0.1866 |
| GSE131907 | 21 | -0.514 | 0.02036 | -0.658 | 0.001626 |
| GSE205335 | 22 | -0.484 | 0.0261 | 0.086 | 0.7109 |
| GSE189357 | 9 | -0.617 | 0.103 | -0.340 | 0.4094 |

## 2. Does this pin CLDN4 vs CLDN7?

Pooled within-cohort rank partial (keratin), 10000 shuffles:

| contrast | pooled ρ CLDN4 | pooled ρ other | Δ (CLDN4 − other) | perm p | unique ρ CLDN4 given other | unique perm p | unique ρ other given CLDN4 | unique perm p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| vs CLDN7 | -0.517 | -0.312 | -0.206 | 0.1662 | -0.465 | 0.0002 | -0.184 | 0.1539 |
| vs CLDN1 | -0.517 | 0.264 | -0.782 | 1.00e-04 | -0.541 | 1.00e-04 | 0.318 | 0.013 |
| vs EPCAM | -0.517 | -0.390 | -0.127 | 0.1784 | -0.369 | 0.0039 | -0.013 | 0.9184 |
| vs MUC1 (context) | -0.517 | -0.162 | -0.356 | 0.0121 | -0.499 | 1.00e-04 | 0.037 | 0.7768 |

DL unique partial (both sides; keratin + the other gene), CLDN4 vs CLDN7:
CLDN4 ρ = -0.422 (p = 0.002565);
CLDN7 ρ = -0.253 (p = 0.308).

Collinearity, DL Spearman of malignant % positive (keratin row is the mean score):
CLDN4–CLDN7 ρ = 0.663 (p = 0.003681, I² = 72.9%);
CLDN4–EPCAM ρ = 0.829;
CLDN4–keratin ρ = 0.539.

**Decision.** Public concordant-4 scRNA does **not** pin CLDN4 over CLDN7. The difference is Δ = -0.206 (permutation p = 0.1662). That is not the bulk failure mode. After keratin and CLDN7 are both partialled, CLDN4 still tracks lower T/NK (pooled unique ρ = -0.465, permutation p = 0.0002; DL ρ = -0.422, p = 0.002565), while CLDN7 given CLDN4 does not (pooled unique ρ = -0.184, permutation p = 0.1539). CLDN4 is not a relabeling of CLDN7. The margin is too small, and too cohort-dependent, to crown it.

The bulk surface ranking had CLDN7, EPCAM, and MUC1 stronger than CLDN4. This T/NK partial does not. CLDN4 has the most negative point estimate. EPCAM also excludes zero (-0.368, p = 0.00688, I² = 0.0%). CLDN7's meta interval crosses zero (I² = 55.6%). Inside GSE131907 the keratin-partial ρ is stronger for CLDN7 than for CLDN4, and dropping GSE205335 reverses the DL rank.

EPCAM is not separable from CLDN4 (Δ = -0.127, permutation p = 0.1784), and its unique partial collapses. CLDN1 has the opposite sign and is separable (Δ = -0.782, permutation p = 1.00e-04). MUC1 is weaker than CLDN4 here (Δ = -0.356, permutation p = 0.0121), the reverse of the bulk rank. What still pins CLDN4 is the private KD co-culture, not a stable public ranking of CLDN4 over CLDN7.

Pooled permutation p for each gene alone, keratin partial, % positive, n = 65:
CLDN1 0.264 (p = 0.0393);
CLDN4 -0.517 (p = 1.00e-04);
CLDN7 -0.312 (p = 0.0145);
EPCAM -0.390 (p = 0.0027);
MUC1 -0.162 (p = 0.2119).

## 3. Sensitivities

Author labels only (GSE131907 + GSE205335). The malignant gate does not use
EPCAM or keratin. Keratin-partial % positive:

| gene | k | N | partial ρ | p | I² |
|---|---:|---:|---:|---:|---:|
| CLDN1 | 2 | 43 | 0.242 | 0.1437 | 0.0% |
| CLDN4 | 2 | 43 | -0.499 | 0.001188 | 0.0% |
| CLDN7 | 2 | 43 | -0.336 | 0.4245 | 85.1% |
| EPCAM | 2 | 43 | -0.352 | 0.02934 | 0.0% |
| MUC1 | 2 | 43 | -0.002 | 0.9946 | 61.5% |

On that subset, CLDN4 vs CLDN7 Δ = -0.217,
perm p = 0.2772. Unique CLDN4 ρ = -0.470
(perm p = 0.0022); unique CLDN7 ρ = -0.206
(perm p = 0.1986).

Marker-gate cohorts only (GSE123902 + GSE189357), same partial. Read these as
within-gate, not as a clean epithelial definition:

| gene | k | N | partial ρ | p | I² |
|---|---:|---:|---:|---:|---:|
| CLDN1 | 2 | 22 | 0.182 | 0.4911 | 0.0% |
| CLDN4 | 2 | 22 | -0.422 | 0.09183 | 0.0% |
| CLDN7 | 2 | 22 | -0.385 | 0.1288 | 0.0% |
| EPCAM | 2 | 22 | -0.406 | 0.1072 | 0.0% |
| MUC1 | 2 | 22 | -0.633 | 0.005267 | 0.0% |

Mean log1p instead of % positive, keratin partial, all four cohorts:

| gene | N | partial ρ | p | I² |
|---|---:|---:|---:|---:|
| CLDN1 | 65 | 0.258 | 0.0649 | 0.0% |
| CLDN4 | 65 | -0.338 | 0.01388 | 0.0% |
| CLDN7 | 65 | -0.244 | 0.3766 | 71.5% |
| EPCAM | 65 | -0.210 | 0.1354 | 0.0% |
| MUC1 | 65 | -0.206 | 0.4175 | 66.0% |

KRT18+KRT19 only (KRT8 left out of the covariate), % positive:
CLDN4 partial ρ = -0.482 (p = 0.0002353, I² = 0.0%);
CLDN7 partial ρ = -0.327 (p = 0.1621).
Δ CLDN4−CLDN7 = -0.211, perm p = 0.1604.

Leave-one-cohort-out, keratin-partial % positive DL ρ:

| dropped | gene | N | ρ | p |
|---|---|---:|---:|---:|
| GSE123902 | CLDN1 | 52 | 0.270 | 0.07937 |
| GSE123902 | CLDN4 | 52 | -0.515 | 0.0003157 |
| GSE123902 | CLDN7 | 52 | -0.336 | 0.2579 |
| GSE123902 | EPCAM | 52 | -0.335 | 0.02739 |
| GSE131907 | CLDN1 | 44 | 0.215 | 0.2175 |
| GSE131907 | CLDN4 | 44 | -0.458 | 0.005157 |
| GSE131907 | CLDN7 | 44 | -0.128 | 0.4651 |
| GSE131907 | EPCAM | 44 | -0.419 | 0.0116 |
| GSE205335 | CLDN1 | 43 | 0.217 | 0.2201 |
| GSE205335 | CLDN4 | 43 | -0.474 | 0.004123 |
| GSE205335 | CLDN7 | 43 | -0.548 | 0.0006065 |
| GSE205335 | EPCAM | 43 | -0.331 | 0.05536 |
| GSE189357 | CLDN1 | 56 | 0.197 | 0.1855 |
| GSE189357 | CLDN4 | 56 | -0.460 | 0.00097 |
| GSE189357 | CLDN7 | 56 | -0.357 | 0.1885 |
| GSE189357 | EPCAM | 56 | -0.385 | 0.007103 |

Within-cohort quartile of % positive, Q4 vs Q1, stacked. Rank-biserial on
raw T/NK fraction, then on the within-cohort rank residual of T/NK after keratin.
Negative means Q4 (higher surface gene) has lower T/NK.

| gene | n_Q1/n_Q4 | r raw | p raw | r keratin residual | p residual |
|---|---|---:|---:|---:|---:|
| CLDN1 | 19/16 | 0.303 | 0.1319 | 0.342 | 0.08813 |
| CLDN4 | 19/16 | -0.724 | 0.0002879 | -0.664 | 0.0008751 |
| CLDN7 | 19/16 | -0.382 | 0.05691 | -0.263 | 0.1909 |
| EPCAM | 19/16 | -0.533 | 0.007685 | -0.467 | 0.01957 |
| MUC1 | 19/16 | -0.401 | 0.04514 | -0.217 | 0.2818 |

## What this does not say

The locked unadjusted CLDN4 association stands (calibration above). This
analysis asks whether that association is CLDN4-specific once keratin, CLDN7,
CLDN1, and EPCAM are allowed to compete. It is not a spatial exclusion test
and it is not the private knockdown co-culture. Bulk rank +0.122 and these
partial ρ values are not the same number.

## Reproduce

GEO files are not in the repo. Put them in `/tmp/geo_c4` (or `GEO_DIR`):

- `GSE123902_RAW.tar`
- `GSE189357_RAW.tar`
- `GSE131907_Lung_Cancer_cell_annotation.txt.gz`
- `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`
- `GSE205335_Lung_IO_CellIdentity.txt.gz`
- `GSE205335_family.soft.gz`
- `GSE205335_Lung_IO_UMI_matrix.rds.gz` (double-gzipped; `extract_gse205335.R` unpacks it)

```
python3 methods/concordant4_surface_tnk_keratin/analyze.py
```

Python does the patient-level scores and the partial / permutation tests.
`Rscript methods/concordant4_surface_tnk_keratin/extract_gse205335.R` reads the RDS.
