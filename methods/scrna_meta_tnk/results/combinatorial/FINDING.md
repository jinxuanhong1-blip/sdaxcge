# Combinatorial search: TACSTD2-high ↔ lower T/NK (CLDN4 same direction)

Patient is the unit. GSE207422 A3 is taken as given. Numbers are the computed
Spearman / DerSimonian–Laird / Stouffer values; this is a search across
definitions and cohort subsets, so the p-values are descriptive.

inferCNV-like exists only for GSE207422 TACSTD2 (wave2 table). CLDN4 is not
in that table and is not invented, so inferCNV rows cannot be scored as
“both genes same direction.”

## Full-pool row (one row, not the only row)

Primary 8-unit mixed-definition grid, mean log1p vs T/NK. Same members as
`results/meta_pooled.tsv`.

| k | N | family | TACSTD2 ρ (p, I²) | CLDN4 ρ (p, I²) | cohorts |
|---:|---:|---|---|---|---|
| 8 | 145 | primary_mixed/tnk/mean | -0.109 (0.378, 39%) | -0.137 (0.129, 0%) | GSE207422+GSE205335+GSE241934_IIT+GSE241934_Real+GSE291670+GSE253013+GSE131907+GSE325414 |
| 8 | 140 | primary_nsclc_swap/tnk/mean | -0.149 (0.169, 21%) | -0.113 (0.222, 0%) | GSE207422+GSE205335_NSCLC+GSE241934_IIT+GSE241934_Real+GSE291670+GSE253013+GSE131907+GSE325414 |

Highlighted both-negative recoveries from the same search (descriptive p):

- GSE207422+GSE291670+GSE253013 · primary_mixed/tnk/mean · k=3 · N=27 · TACSTD2 ρ=-0.585 p=0.00452 I²=0% · CLDN4 ρ=-0.378 p=0.157 I²=26%
- GSE207422+GSE241934_IIT+GSE291670+GSE253013+GSE325414 · primary_mixed/tnk/mean · k=5 · N=63 · TACSTD2 ρ=-0.306 p=0.0447 I²=11% · CLDN4 ρ=-0.197 p=0.166 I²=0%
- GSE253013+GSE291670 · marker/tnk/mean · k=2 · N=15 · TACSTD2 ρ=-0.666 p=0.016 I²=0% · CLDN4 ρ=-0.582 p=0.102 I²=29%
- GSE253013+GSE291670 · marker/tnk/pct · k=2 · N=15 · TACSTD2 ρ=-0.736 p=0.00766 I²=9% · CLDN4 ρ=-0.714 p=0.0217 I²=23%

## Primary-grid members (the 8 units behind the full-pool row)

| cohort | malig | T/NK | score | n | TACSTD2 ρ (p) | CLDN4 ρ (p) | both ρ<0 |
|---|---|---|---|---:|---|---|---|
| GSE207422 | author_DRMref | tnk | mean | 12 | -0.490 (0.106) | -0.091 (0.779) | yes |
| GSE205335 | author_malig | tnk | mean | 22 | +0.284 (0.2) | -0.200 (0.371) | no |
| GSE241934_IIT | author_epi | tnk | mean | 11 | -0.073 (0.832) | -0.064 (0.853) | yes |
| GSE241934_Real | author_epi | tnk | mean | 24 | -0.167 (0.436) | +0.178 (0.405) | no |
| GSE291670 | marker_malig | tnk | mean | 6 | -0.543 (0.266) | -0.829 (0.0416) | yes |
| GSE253013 | marker_malig | tnk | mean | 9 | -0.717 (0.0298) | -0.333 (0.381) | yes |
| GSE131907 | author_epi | tnk | mean | 36 | +0.168 (0.326) | -0.210 (0.218) | no |
| GSE325414 | author_malig | tnk | mean | 25 | -0.073 (0.728) | -0.119 (0.57) | yes |

## Primary-grid subset pools with both genes pooled ρ < 0 (k≥2)

All 2⁸−1 = 255 nonempty subsets of the 8 primary units were enumerated.
Rows below are both-negative subsets with N≥20, lowest TACSTD2 RE p first
(top 20).

| k | N | family | TACSTD2 ρ (p, I²) | CLDN4 ρ (p, I²) | cohorts |
|---:|---:|---|---|---|---|
| 3 | 27 | primary_mixed/tnk/mean | -0.585 (0.00452, 0%) | -0.378 (0.157, 26%) | GSE207422+GSE291670+GSE253013 |
| 2 | 21 | primary_mixed/tnk/mean | -0.593 (0.0083, 0%) | -0.191 (0.454, 0%) | GSE207422+GSE253013 |
| 4 | 51 | primary_mixed/tnk/mean | -0.380 (0.0126, 0%) | -0.184 (0.438, 45%) | GSE207422+GSE241934_Real+GSE291670+GSE253013 |
| 4 | 38 | primary_mixed/tnk/mean | -0.451 (0.0132, 0%) | -0.266 (0.177, 5%) | GSE207422+GSE241934_IIT+GSE291670+GSE253013 |
| 5 | 62 | primary_mixed/tnk/mean | -0.331 (0.0183, 0%) | -0.123 (0.494, 27%) | GSE207422+GSE241934_IIT+GSE241934_Real+GSE291670+GSE253013 |
| 5 | 76 | primary_mixed/tnk/mean | -0.289 (0.0314, 10%) | -0.124 (0.436, 28%) | GSE207422+GSE241934_Real+GSE291670+GSE253013+GSE325414 |
| 6 | 87 | primary_mixed/tnk/mean | -0.252 (0.0323, 0%) | -0.095 (0.463, 10%) | GSE207422+GSE241934_IIT+GSE241934_Real+GSE291670+GSE253013+GSE325414 |
| 3 | 32 | primary_mixed/tnk/mean | -0.442 (0.0389, 17%) | -0.147 (0.477, 0%) | GSE207422+GSE241934_IIT+GSE253013 |
| 4 | 52 | primary_mixed/tnk/mean | -0.382 (0.0413, 25%) | -0.234 (0.159, 8%) | GSE207422+GSE291670+GSE253013+GSE325414 |
| 5 | 63 | primary_mixed/tnk/mean | -0.306 (0.0447, 11%) | -0.197 (0.166, 0%) | GSE207422+GSE241934_IIT+GSE291670+GSE253013+GSE325414 |
| 5 | 81 | primary_mixed/tnk/mean | -0.244 (0.0552, 8%) | -0.034 (0.781, 0%) | GSE207422+GSE241934_IIT+GSE241934_Real+GSE253013+GSE325414 |
| 4 | 70 | primary_mixed/tnk/mean | -0.289 (0.061, 27%) | -0.030 (0.818, 0%) | GSE207422+GSE241934_Real+GSE253013+GSE325414 |
| 4 | 50 | primary_mixed/tnk/mean | -0.300 (0.0732, 7%) | -0.181 (0.454, 45%) | GSE241934_IIT+GSE241934_Real+GSE291670+GSE253013 |
| 3 | 39 | primary_mixed/tnk/mean | -0.404 (0.0745, 26%) | -0.294 (0.41, 64%) | GSE241934_Real+GSE291670+GSE253013 |
| 3 | 42 | primary_mixed/tnk/mean | -0.299 (0.0763, 0%) | -0.182 (0.56, 60%) | GSE207422+GSE241934_Real+GSE291670 |
| 3 | 26 | primary_mixed/tnk/mean | -0.441 (0.0814, 18%) | -0.379 (0.169, 27%) | GSE241934_IIT+GSE291670+GSE253013 |
| 4 | 57 | primary_mixed/tnk/mean | -0.304 (0.0882, 28%) | -0.134 (0.367, 0%) | GSE207422+GSE241934_IIT+GSE253013+GSE325414 |
| 4 | 53 | primary_mixed/tnk/mean | -0.257 (0.0927, 0%) | -0.106 (0.626, 40%) | GSE207422+GSE241934_IIT+GSE241934_Real+GSE291670 |
| 3 | 46 | primary_mixed/tnk/mean | -0.389 (0.0935, 47%) | -0.149 (0.363, 0%) | GSE207422+GSE253013+GSE325414 |
| 5 | 75 | primary_mixed/tnk/mean | -0.213 (0.0943, 0%) | -0.121 (0.454, 28%) | GSE241934_IIT+GSE241934_Real+GSE291670+GSE253013+GSE325414 |

Primary-grid both-negative subset pools (k≥2): 197 of 247 enumerated k≥2 subsets.

## Single-cohort combinations with both genes ρ < 0 (n≥4)

Every available (cohort × malignant def × T/NK def × score) that has both
genes. inferCNV is absent here because CLDN4 is missing.

| cohort | malig | T/NK | score | n | TACSTD2 ρ (p) | CLDN4 ρ (p) |
|---|---|---|---|---:|---|---|
| GSE253013 | marker_malig | tnk | pct | 9 | -0.833 (0.00527) | -0.533 (0.139) |
| GSE253013 | marker_malig | tnk | mean_cp10k | 9 | -0.783 (0.0125) | -0.550 (0.125) |
| GSE253013 | marker_malig | tnk | mean | 9 | -0.717 (0.0298) | -0.333 (0.381) |
| GSE291670 | marker_epi | tnk | mean | 6 | -0.600 (0.208) | -0.829 (0.0416) |
| GSE291670 | marker_malig | tnk | mean | 6 | -0.543 (0.266) | -0.829 (0.0416) |
| GSE253013 | marker_epi | tnk | pct | 9 | -0.517 (0.154) | -0.067 (0.865) |
| GSE207422 | author_DRMref | tnk | mean | 12 | -0.490 (0.106) | -0.091 (0.779) |
| GSE291670 | marker_malig | tnk | pct | 6 | -0.429 (0.397) | -0.886 (0.0188) |
| GSE291670 | marker_epi | tnk_umi | mean | 6 | -0.429 (0.397) | -0.543 (0.266) |
| GSE207422 | marker_malig | cd8 | mean | 7 | -0.393 (0.383) | -0.071 (0.879) |
| GSE291670 | marker_malig | tnk_umi | mean | 6 | -0.371 (0.468) | -0.543 (0.266) |
| GSE131907_tLung | author_epi | tnk | pct | 11 | -0.264 (0.433) | -0.400 (0.223) |
| GSE205335_NSCLC | author_malig | tnk | pct | 17 | -0.248 (0.338) | -0.206 (0.428) |
| GSE207422 | marker_epi | cd8 | mean | 12 | -0.196 (0.542) | -0.007 (0.983) |
| GSE131907 | author_malig | tnk | pct | 21 | -0.151 (0.515) | -0.522 (0.0152) |
| GSE241934_IIT | author_epi | cd8 | mean | 11 | -0.136 (0.689) | -0.109 (0.75) |
| GSE241934_Real | author_epi | tnk | pct | 24 | -0.123 (0.567) | -0.034 (0.875) |
| GSE253013 | marker_epi | tnk | mean | 9 | -0.117 (0.765) | -0.200 (0.606) |
| GSE205335_NSCLC | author_malig | tnk | mean | 17 | -0.113 (0.667) | -0.022 (0.933) |
| GSE131907 | author_malig | cd8 | pct | 21 | -0.078 (0.737) | -0.573 (0.00666) |
| GSE325414 | author_malig | tnk | mean | 25 | -0.073 (0.728) | -0.119 (0.57) |
| GSE241934_IIT | author_epi | tnk | mean | 11 | -0.073 (0.832) | -0.064 (0.853) |
| GSE207422 | marker_malig | cd8 | pct | 7 | -0.071 (0.879) | -0.179 (0.702) |
| GSE207422 | marker_malig_nsclc | cd8 | pct | 9 | -0.033 (0.932) | -0.133 (0.732) |
| GSE207422 | marker_malig_nsclc | cd8 | mean | 9 | -0.017 (0.966) | -0.133 (0.732) |
| GSE131907 | author_epi | tnk | pct | 36 | -0.009 (0.96) | -0.365 (0.0286) |

Count: 26 both-negative single-cohort combos (of 56 paired combos with n≥4).

## Aligned-family subset pools with both genes pooled ρ < 0 (k≥2)

Families hold malignant / T/NK / score buckets fixed, then enumerate cohort
subsets. Showing the 15 lowest TACSTD2 RE p among both-negative subsets with N≥15.

| k | N | family | TACSTD2 ρ (p, I²) | CLDN4 ρ (p, I²) | cohorts |
|---:|---:|---|---|---|---|
| 2 | 15 | marker/tnk/pct | -0.736 (0.00766, 9%) | -0.714 (0.0217, 23%) | GSE253013+GSE291670 |
| 2 | 15 | marker/tnk/mean | -0.666 (0.016, 0%) | -0.582 (0.102, 29%) | GSE253013+GSE291670 |
| 3 | 22 | marker/tnk/mean | -0.513 (0.0409, 0%) | -0.362 (0.314, 44%) | GSE207422+GSE253013+GSE291670 |
| 4 | 72 | author/tnk/mean | -0.174 (0.173, 0%) | -0.003 (0.981, 0%) | GSE207422+GSE241934_IIT+GSE241934_Real+GSE325414 |
| 2 | 23 | author/tnk/mean | -0.307 (0.19, 0%) | -0.078 (0.747, 0%) | GSE207422+GSE241934_IIT |
| 3 | 22 | marker/tnk/pct | -0.499 (0.194, 55%) | -0.551 (0.091, 41%) | GSE207422+GSE253013+GSE291670 |
| 2 | 16 | marker/tnk/mean | -0.475 (0.23, 44%) | -0.105 (0.738, 0%) | GSE207422+GSE253013 |
| 3 | 48 | author/tnk/mean | -0.178 (0.262, 0%) | -0.101 (0.525, 0%) | GSE207422+GSE241934_IIT+GSE325414 |
| 2 | 37 | author/tnk/mean | -0.229 (0.289, 27%) | -0.111 (0.535, 0%) | GSE207422+GSE325414 |
| 5 | 93 | author/tnk/mean | -0.116 (0.304, 0%) | -0.099 (0.382, 0%) | GSE131907+GSE207422+GSE241934_IIT+GSE241934_Real+GSE325414 |
| 4 | 82 | author/tnk/mean | -0.121 (0.31, 0%) | -0.104 (0.42, 14%) | GSE131907+GSE207422+GSE241934_Real+GSE325414 |
| 4 | 68 | author/tnk/mean | -0.133 (0.319, 0%) | -0.093 (0.523, 14%) | GSE131907+GSE207422+GSE241934_IIT+GSE241934_Real |
| 3 | 57 | author/tnk/mean | -0.148 (0.341, 13%) | -0.103 (0.596, 43%) | GSE131907+GSE207422+GSE241934_Real |
| 3 | 56 | author/tnk/pct | -0.136 (0.348, 0%) | -0.213 (0.289, 45%) | GSE131907+GSE241934_IIT+GSE241934_Real |
| 2 | 15 | marker_epi/tnk/mean | -0.300 (0.354, 0%) | -0.543 (0.208, 48%) | GSE253013+GSE291670 |

Aligned-family both-negative subset pools (k≥2): 57 of 102 enumerated k≥2 aligned pools.

## inferCNV-like (GSE207422 TACSTD2 only)

CLDN4 is not in the wave2 table. These rows are TACSTD2 vs T/NK only.

| malig | T/NK | score | n | TACSTD2 ρ (p) |
|---|---|---|---:|---|
| infercnv_p95 | tnk_drm | pct | 12 | -0.387 (0.214) |
| infercnv_p95 | tnk_drm | mean | 12 | -0.350 (0.265) |
| infercnv_p90 | tnk_drm | mean | 12 | -0.343 (0.276) |
| infercnv_p90 | tnk_drm | pct | 12 | -0.324 (0.304) |
| infercnv_gmm | tnk_drm | mean | 12 | -0.266 (0.404) |
| infercnv_gmm | tnk_drm | pct | 12 | -0.126 (0.697) |
| infercnv_gmm | tnk | mean | 12 | +0.000 (1) |
| infercnv_gmm | tnk_umi | mean | 12 | +0.035 (0.914) |
| infercnv_p90 | tnk | mean | 12 | +0.063 (0.846) |
| infercnv_p95 | tnk | mean | 12 | +0.077 (0.812) |
| infercnv_p90 | tnk_umi | mean | 12 | +0.098 (0.762) |
| infercnv_p95 | tnk_umi | mean | 12 | +0.119 (0.713) |
| infercnv_gmm | tnk | pct | 12 | +0.140 (0.665) |
| infercnv_p95 | tnk | pct | 12 | +0.155 (0.631) |
| infercnv_gmm | tnk_umi | pct | 12 | +0.168 (0.602) |
| infercnv_p90 | tnk | pct | 12 | +0.190 (0.554) |
| infercnv_p95 | tnk_umi | pct | 12 | +0.197 (0.539) |
| infercnv_p90 | tnk_umi | pct | 12 | +0.225 (0.481) |

Figures: `figures/supportive_singles_forest.png`, `figures/supportive_subset_bars.png`,
`figures/supportive_trio_GSE207422_GSE291670_GSE253013.png`,
`figures/supportive_marker_malig_*_GSE253013_GSE291670.png`,
`figures/primary_full_pool_members_forest.png`, `figures/primary_both_negative_members_forest.png`,
`figures/primary_nsclc_swap_members_forest.png`, `figures/scatter_GSE253013_*.png`,
`figures/scatter_GSE291670_*.png`, `figures/scatter_GSE207422_A3_given*.png`.

