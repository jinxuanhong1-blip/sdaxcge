# Strict malignant T/NK recut: combinatorial CLDN4 and TACSTD2

ADDITIVE. Patient is the unit. Existing scRNA patient tables are reused;
GEO is not re-downloaded. **Keep only author-malignant, marker-malignant,
or DRMref.** All-epithelial rows are dropped (GSE241934 residual Epi,
GSE131907/GSE253013/GSE207422/GSE291670 epithelium columns, E-MTAB-13526).
Leftover GSE267108/GSE274595 malignant splits have n<4 and cannot enter.
inferCNV-like is not author/marker/DRMref and is dropped (CLDN4 also absent).

GSE207422 A3 TACSTD2 is taken as given (n=12, mean ρ=−0.490, p=0.106).
CLDN4 is from the same given DRMref table. p-values are descriptive
(many subsets). The ρ ≤ −0.35 rule is **both genes**, not a p-value gate.

Locked search order: (1) the 6 primary-grid mean-vs-T/NK members,
(2) other mean singles, (3) %pos singles, (4) k=2 aligned mean,
(5) larger subsets. First hit is reported with its honest n and p;
every other both-≤−0.35 cut is in `tables/combo_table.tsv`.

## First cut with both genes ρ ≤ −0.35

- **GSE291670 · marker_malig/tnk/mean (primary-grid member)** · k=1 · N=6 · TACSTD2 ρ=-0.543 p=0.266 · CLDN4 ρ=-0.829 p=0.0416

GSE253013 marker_malig / tnk / mean is a near-miss: TACSTD2 ρ=−0.717
p=0.0298, CLDN4 ρ=−0.333 p=0.381 (CLDN4 just above −0.35). The same
two marker-malignant cohorts pooled (mean) are N=15, TACSTD2 ρ=−0.666
p=0.016 I²=0%, CLDN4 ρ=−0.582 p=0.102 I²=29%.

## Strict full-pool row (one row, not the answer)

6-unit grid after dropping all-epithelial. Mean log1p vs T/NK.

| k | N | family | TACSTD2 ρ (p, I²) | CLDN4 ρ (p, I²) | cohorts |
|---:|---:|---|---|---|---|
| 6 | 95 | strict_primary/tnk/mean | -0.164 (p=0.332, I²=49%) | -0.260 (p=0.0195, I²=0%) | GSE207422+GSE205335+GSE291670+GSE253013+GSE131907+GSE325414 |
| 6 | 90 | strict_nsclc_swap/tnk/mean | -0.216 (p=0.111, I²=21%) | -0.231 (p=0.0459, I²=0%) | GSE207422+GSE205335_NSCLC+GSE291670+GSE253013+GSE131907+GSE325414 |

## Strict primary-grid members

| cohort | malig | T/NK | score | n | TACSTD2 ρ (p) | CLDN4 ρ (p) | both ρ<0 | both ρ≤−0.35 |
|---|---|---|---|---:|---|---|---|---|
| GSE207422 | author_DRMref | tnk | mean | 12 | -0.490 (0.106) | -0.091 (0.779) | yes | no |
| GSE205335 | author_malig | tnk | mean | 22 | +0.284 (0.2) | -0.200 (0.371) | no | no |
| GSE291670 | marker_malig | tnk | mean | 6 | -0.543 (0.266) | -0.829 (0.0416) | yes | yes |
| GSE253013 | marker_malig | tnk | mean | 9 | -0.717 (0.0298) | -0.333 (0.381) | yes | no |
| GSE131907 | author_malig | tnk | mean | 21 | +0.082 (0.724) | -0.396 (0.0755) | no | no |
| GSE325414 | author_malig | tnk | mean | 25 | -0.073 (0.728) | -0.119 (0.57) | yes | no |

## Cuts with both genes ρ ≤ −0.35 (honest n and p)

Singles first, then pooled subsets. Full combo table:
`tables/combo_table.tsv`.

| k | N | family | TACSTD2 ρ (p, I²) | CLDN4 ρ (p, I²) | cohorts |
|---:|---:|---|---|---|---|
| 1 | 9 | strict_primary/tnk/pct | -0.833 (p=0.00527, I²=0%) | -0.533 (p=0.139, I²=0%) | GSE253013 |
| 1 | 6 | strict_primary/tnk/mean | -0.543 (p=0.266, I²=0%) | -0.829 (p=0.0416, I²=0%) | GSE291670 |
| 1 | 6 | strict_primary/tnk/pct | -0.429 (p=0.397, I²=0%) | -0.886 (p=0.0188, I²=0%) | GSE291670 |
| 2 | 31 | strict_primary/tnk/pct | -0.425 (p=0.52, I²=89%) | -0.460 (p=0.0129, I²=0%) | GSE205335+GSE253013 |
| 2 | 30 | strict_primary/tnk/pct | -0.553 (p=0.232, I²=80%) | -0.525 (p=0.00428, I²=0%) | GSE253013+GSE131907 |
| 2 | 18 | strict_primary/tnk/mean | -0.503 (p=0.0551, I²=0%) | -0.490 (p=0.318, I²=63%) | GSE207422+GSE291670 |
| 2 | 15 | strict_primary/tnk/pct | -0.736 (p=0.00766, I²=9%) | -0.714 (p=0.0217, I²=23%) | GSE291670+GSE253013 |
| 2 | 15 | strict_primary/tnk/mean | -0.666 (p=0.016, I²=0%) | -0.582 (p=0.102, I²=29%) | GSE291670+GSE253013 |
| 3 | 37 | strict_primary/tnk/pct | -0.413 (p=0.367, I²=79%) | -0.550 (p=0.00373, I²=12%) | GSE205335+GSE291670+GSE253013 |
| 3 | 36 | strict_primary/tnk/pct | -0.507 (p=0.114, I²=60%) | -0.588 (p=0.000461, I²=0%) | GSE291670+GSE253013+GSE131907 |
| 3 | 36 | strict_primary/tnk/mean | -0.371 (p=0.267, I²=59%) | -0.453 (p=0.0112, I²=0%) | GSE291670+GSE253013+GSE131907 |
| 3 | 27 | strict_primary/tnk/mean | -0.585 (p=0.00452, I²=0%) | -0.378 (p=0.157, I²=26%) | GSE207422+GSE291670+GSE253013 |
| 3 | 22 | marker_malig/tnk/pct | -0.499 (p=0.194, I²=55%) | -0.551 (p=0.091, I²=41%) | GSE207422+GSE253013+GSE291670 |
| 3 | 22 | marker_malig/tnk/mean | -0.513 (p=0.0409, I²=0%) | -0.362 (p=0.314, I²=44%) | GSE207422+GSE253013+GSE291670 |
| 4 | 48 | strict_primary/tnk/mean | -0.379 (p=0.106, I²=46%) | -0.370 (p=0.0197, I²=0%) | GSE207422+GSE291670+GSE253013+GSE131907 |

Count: 15 enumerated cuts with both genes ρ ≤ −0.35.

## Single-cohort strict T/NK (n≥4, both genes present)

| cohort | malig | T/NK | score | n | TACSTD2 ρ (p) | CLDN4 ρ (p) | both ρ≤−0.35 |
|---|---|---|---|---:|---|---|---|
| GSE253013 | marker_malig | tnk | pct | 9 | -0.833 (0.00527) | -0.533 (0.139) | yes |
| GSE253013 | marker_malig | tnk | mean_cp10k | 9 | -0.783 (0.0125) | -0.550 (0.125) | yes |
| GSE253013 | marker_malig | tnk | mean | 9 | -0.717 (0.0298) | -0.333 (0.381) | no |
| GSE291670 | marker_malig | tnk | mean | 6 | -0.543 (0.266) | -0.829 (0.0416) | yes |
| GSE207422 | author_DRMref | tnk | mean | 12 | -0.490 (0.106) | -0.091 (0.779) | no |
| GSE291670 | marker_malig | tnk | pct | 6 | -0.429 (0.397) | -0.886 (0.0188) | yes |
| GSE291670 | marker_malig | tnk_umi | mean | 6 | -0.371 (0.468) | -0.543 (0.266) | yes |
| GSE205335_NSCLC | author_malig | tnk | pct | 17 | -0.248 (0.338) | -0.206 (0.428) | no |
| GSE131907 | author_malig | tnk | pct | 21 | -0.151 (0.515) | -0.522 (0.0152) | no |
| GSE205335_NSCLC | author_malig | tnk | mean | 17 | -0.113 (0.667) | -0.022 (0.933) | no |
| GSE325414 | author_malig | tnk | mean | 25 | -0.073 (0.728) | -0.119 (0.57) | no |

Both-negative singles: 11 of 19 paired strict T/NK combos with n≥4. Both ρ≤−0.35: 5.

## Primary-grid subset pools with both genes ρ < 0 (k≥2)

All 2⁶−1 = 63 nonempty subsets of the 6 strict units were enumerated.
Rows below: both-negative, N≥15, lowest TACSTD2 RE p first (top 15).

| k | N | family | TACSTD2 ρ (p, I²) | CLDN4 ρ (p, I²) | both ρ≤−0.35 | cohorts |
|---:|---:|---|---|---|---|---|
| 3 | 27 | strict_primary/tnk/mean | -0.585 (p=0.00452, I²=0%) | -0.378 (p=0.157, I²=26%) | yes | GSE207422+GSE291670+GSE253013 |
| 2 | 21 | strict_primary/tnk/mean | -0.593 (p=0.0083, I²=0%) | -0.191 (p=0.454, I²=0%) | no | GSE207422+GSE253013 |
| 2 | 15 | strict_primary/tnk/mean | -0.666 (p=0.016, I²=0%) | -0.582 (p=0.102, I²=29%) | yes | GSE291670+GSE253013 |
| 4 | 52 | strict_primary/tnk/mean | -0.382 (p=0.0413, I²=25%) | -0.234 (p=0.159, I²=8%) | no | GSE207422+GSE291670+GSE253013+GSE325414 |
| 2 | 18 | strict_primary/tnk/mean | -0.503 (p=0.0551, I²=0%) | -0.490 (p=0.318, I²=63%) | yes | GSE207422+GSE291670 |
| 3 | 46 | strict_primary/tnk/mean | -0.389 (p=0.0935, I²=47%) | -0.149 (p=0.363, I²=0%) | no | GSE207422+GSE253013+GSE325414 |
| 4 | 48 | strict_primary/tnk/mean | -0.379 (p=0.106, I²=46%) | -0.370 (p=0.0197, I²=0%) | yes | GSE207422+GSE291670+GSE253013+GSE131907 |
| 5 | 73 | strict_primary/tnk/mean | -0.267 (p=0.118, I²=36%) | -0.279 (p=0.029, I²=0%) | no | GSE207422+GSE291670+GSE253013+GSE131907+GSE325414 |
| 3 | 40 | strict_primary/tnk/mean | -0.394 (p=0.147, I²=44%) | -0.343 (p=0.167, I²=34%) | no | GSE291670+GSE253013+GSE325414 |
| 3 | 43 | strict_primary/tnk/mean | -0.238 (p=0.157, I²=0%) | -0.261 (p=0.272, I²=36%) | no | GSE207422+GSE291670+GSE325414 |
| 4 | 67 | strict_primary/tnk/mean | -0.251 (p=0.19, I²=48%) | -0.233 (p=0.0778, I²=0%) | no | GSE207422+GSE253013+GSE131907+GSE325414 |
| 3 | 42 | strict_primary/tnk/mean | -0.366 (p=0.2, I²=62%) | -0.306 (p=0.0691, I²=0%) | no | GSE207422+GSE253013+GSE131907 |
| 5 | 74 | strict_primary/tnk/mean | -0.249 (p=0.236, I²=56%) | -0.216 (p=0.0919, I²=0%) | no | GSE207422+GSE205335+GSE291670+GSE253013+GSE325414 |
| 4 | 49 | strict_primary/tnk/mean | -0.350 (p=0.245, I²=67%) | -0.272 (p=0.09, I²=0%) | no | GSE207422+GSE205335+GSE291670+GSE253013 |
| 3 | 36 | strict_primary/tnk/mean | -0.371 (p=0.267, I²=59%) | -0.453 (p=0.0112, I²=0%) | yes | GSE291670+GSE253013+GSE131907 |

Primary-grid both-negative k≥2 subsets: 50 of 57 enumerated. Both ρ≤−0.35: 6.

## Aligned-family cuts with both genes ρ ≤ −0.35

| k | N | family | TACSTD2 ρ (p, I²) | CLDN4 ρ (p, I²) | cohorts |
|---:|---:|---|---|---|---|
| 1 | 6 | marker_malig/tnk/mean | -0.543 (p=0.266, I²=0%) | -0.829 (p=0.0416, I²=0%) | GSE291670 |
| 1 | 9 | marker_malig/tnk/pct | -0.833 (p=0.00527, I²=0%) | -0.533 (p=0.139, I²=0%) | GSE253013 |
| 1 | 6 | marker_malig/tnk/pct | -0.429 (p=0.397, I²=0%) | -0.886 (p=0.0188, I²=0%) | GSE291670 |
| 1 | 6 | marker_malig/tnk_umi/mean | -0.371 (p=0.468, I²=0%) | -0.543 (p=0.266, I²=0%) | GSE291670 |
| 2 | 15 | marker_malig/tnk/pct | -0.736 (p=0.00766, I²=9%) | -0.714 (p=0.0217, I²=23%) | GSE253013+GSE291670 |
| 2 | 15 | marker_malig/tnk/mean | -0.666 (p=0.016, I²=0%) | -0.582 (p=0.102, I²=29%) | GSE253013+GSE291670 |
| 3 | 22 | marker_malig/tnk/pct | -0.499 (p=0.194, I²=55%) | -0.551 (p=0.091, I²=41%) | GSE207422+GSE253013+GSE291670 |
| 3 | 22 | marker_malig/tnk/mean | -0.513 (p=0.0409, I²=0%) | -0.362 (p=0.314, I²=44%) | GSE207422+GSE253013+GSE291670 |

## Dropped rows

| source | reason | n | detail |
|---|---|---:|---|
| GSE241934_IIT | all-epithelial | 11 | author residual Epi only; no author-malignant / marker-malignant / DRMref column |
| GSE241934_Real | all-epithelial | 24 | author residual Epi only |
| GSE131907 author_epi | all-epithelial | 36 | author epithelium columns dropped; author_malig n_mal>=20 is kept |
| GSE253013 author_epi / marker_epi | all-epithelial | 9 | author Epithelial and tumor epithelium columns dropped; marker_malig kept |
| GSE207422 marker_epi | all-epithelial | 12 | post all-epithelial extra dropped; marker_malig and DRMref kept |
| GSE291670 marker_epi | all-epithelial | 6 | epi_* columns dropped; marker_malig kept |
| E-MTAB-13526 | all-epithelial | 13 | author Cell types epithelium only; no malignant split |
| GSE267108 malignant | malignant n<4 | 1 | leftover stats: mal_TACSTD2 n=1; epithelial-only would be dropped anyway |
| GSE274595 malignant | malignant n<4 | 0 | leftover stats: mal_TACSTD2 n=0 |
| GSE207422 inferCNV-like | not author/marker/DRMref | 12 | wave2 inferCNV defs dropped; CLDN4 absent in that table |

## Extra figures

See `figures/`: strict primary forest, NSCLC swap, both-≤−0.35 singles,
marker pair, trio, subset bars, and patient scatters.

## Reproduce

```bash
python3 methods/strict_malig_tnk/combinatorial_search.py
```

