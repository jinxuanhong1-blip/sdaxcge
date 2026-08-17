# Combinatorial search: malignant CLDN4 vs T/NK (and CXCL13+ / B / cyto-high T)

CLDN4-first. Patient is the unit. Numbers are computed Spearman /
DerSimonian–Laird / Stouffer values across definitions and cohort subsets;
p-values are descriptive.

TACSTD2-primary combinations in PR #279 / #271 / #274 are **taken as given**
and are not re-audited. GSE207422 A3 TACSTD2 is taken as given. CLDN4 on the
same given 12-patient table is one honest row (not an A3 re-cut).

GSE253013 uses the existing 9-patient tumor extract. The 9 GB GEO RDS was
not downloaded.

## Full-pool row (one row, not the answer)

Primary 8-unit mixed-definition grid, mean log1p vs T/NK. Same members as
PR #279. This is the merge-everything row.

| k | N | family | CLDN4 ρ (p, I²) | cohorts |
|---:|---:|---|---|---|
| 8 | 145 | primary_mixed | -0.137 (0.129, 0%) | GSE207422+GSE205335+GSE241934_IIT+GSE241934_Real+GSE291670+GSE253013+GSE131907+GSE325414 |

The search below is the answer: combinations where **CLDN4 ρ < 0** vs T/NK
or CXCL13+, with honest n and p.

## Combinations that recover CLDN4 ρ < 0 vs T/NK

Highlighted cuts (descriptive p; ranked by CLDN4, not TACSTD2):

- GSE131907+GSE205335 · author/tnk/pct · k=2 · N=43 · CLDN4 ρ=-0.479 p=0.00152 I²=0%
- GSE253013+GSE291670 · marker/tnk/pct · k=2 · N=15 · CLDN4 ρ=-0.714 p=0.0217 I²=23%
- GSE253013+GSE291670 · marker/tnk/mean · k=2 · N=15 · CLDN4 ρ=-0.582 p=0.102 I²=29%
- GSE291670 · primary_mixed · k=1 · N=6 · CLDN4 ρ=-0.829 p=0.0416 I²=0%
- GSE291670+GSE253013 · primary_mixed · k=2 · N=15 · CLDN4 ρ=-0.582 p=0.102 I²=29%
- GSE205335+GSE291670+GSE253013+GSE131907 · primary_mixed · k=4 · N=73 · CLDN4 ρ=-0.265 p=0.0342 I²=0%
- GSE207422+GSE291670+GSE253013 · primary_mixed · k=3 · N=27 · CLDN4 ρ=-0.378 p=0.157 I²=26%
- GSE148071+GSE207422+GSE253013 · tls/cxcl13pos/mean · k=3 · N=60 · CLDN4 ρ=-0.425 p=0.00121 I²=0%
- GSE148071+GSE207422 · tls/cxcl13pos/mean · k=2 · N=51 · CLDN4 ρ=-0.428 p=0.00216 I²=0%
- GSE207422+GSE253013 · tls/cxcl13_mean/mean · k=2 · N=24 · CLDN4 ρ=-0.552 p=0.00833 I²=0%
- seven primary members with single-cohort CLDN4 ρ<0 (drop GSE241934_Real) · k=7 · N=121 · CLDN4 ρ=-0.202 p=0.0403 I²=0%

## Primary-grid members (the 8 units behind the full-pool row)

| cohort | malig | immune | score | n | CLDN4 ρ (p) | CLDN4 ρ<0 |
|---|---|---|---|---:|---|---|
| GSE207422 | author_DRMref | tnk | mean | 12 | -0.091 (0.779) | yes |
| GSE205335 | author_malig | tnk | mean | 22 | -0.200 (0.371) | yes |
| GSE241934_IIT | author_epi | tnk | mean | 11 | -0.064 (0.853) | yes |
| GSE241934_Real | author_epi | tnk | mean | 24 | +0.178 (0.405) | no |
| GSE291670 | marker_malig | tnk | mean | 6 | -0.829 (0.0416) | yes |
| GSE253013 | marker_malig | tnk | mean | 9 | -0.333 (0.381) | yes |
| GSE131907 | author_epi | tnk | mean | 36 | -0.210 (0.218) | yes |
| GSE325414 | author_malig | tnk | mean | 25 | -0.119 (0.57) | yes |

## Primary-grid subset pools with CLDN4 pooled ρ < 0 (k≥2)

All 2⁸−1 = 255 nonempty subsets of the 8 primary units were enumerated.
Rows below are CLDN4-negative subsets with N≥20, lowest CLDN4 RE p first
(top 20). Full-pool is not repeated here.

| k | N | family | CLDN4 ρ (p, I²) | cohorts |
|---:|---:|---|---|---|
| 4 | 73 | primary_mixed | -0.265 (0.0342, 0%) | GSE205335+GSE291670+GSE253013+GSE131907 |
| 5 | 98 | primary_mixed | -0.227 (0.0353, 0%) | GSE205335+GSE291670+GSE253013+GSE131907+GSE325414 |
| 6 | 110 | primary_mixed | -0.214 (0.0371, 0%) | GSE207422+GSE205335+GSE291670+GSE253013+GSE131907+GSE325414 |
| 5 | 85 | primary_mixed | -0.243 (0.038, 0%) | GSE207422+GSE205335+GSE291670+GSE253013+GSE131907 |
| 6 | 109 | primary_mixed | -0.213 (0.0391, 0%) | GSE205335+GSE241934_IIT+GSE291670+GSE253013+GSE131907+GSE325414 |
| 5 | 84 | primary_mixed | -0.242 (0.0402, 0%) | GSE205335+GSE241934_IIT+GSE291670+GSE253013+GSE131907 |
| 7 | 121 | primary_mixed | -0.202 (0.0403, 0%) | GSE207422+GSE205335+GSE241934_IIT+GSE291670+GSE253013+GSE131907+GSE325414 |
| 6 | 96 | primary_mixed | -0.225 (0.0431, 0%) | GSE207422+GSE205335+GSE241934_IIT+GSE291670+GSE253013+GSE131907 |
| 4 | 89 | primary_mixed | -0.219 (0.052, 0%) | GSE205335+GSE291670+GSE131907+GSE325414 |
| 5 | 101 | primary_mixed | -0.205 (0.0534, 0%) | GSE207422+GSE205335+GSE291670+GSE131907+GSE325414 |
| 4 | 76 | primary_mixed | -0.234 (0.0562, 0%) | GSE207422+GSE205335+GSE291670+GSE131907 |
| 5 | 100 | primary_mixed | -0.204 (0.0562, 0%) | GSE205335+GSE241934_IIT+GSE291670+GSE131907+GSE325414 |
| 6 | 112 | primary_mixed | -0.194 (0.0574, 0%) | GSE207422+GSE205335+GSE241934_IIT+GSE291670+GSE131907+GSE325414 |
| 4 | 63 | primary_mixed | -0.259 (0.0588, 0%) | GSE207422+GSE291670+GSE253013+GSE131907 |
| 5 | 88 | primary_mixed | -0.217 (0.059, 0%) | GSE207422+GSE291670+GSE253013+GSE131907+GSE325414 |
| 4 | 75 | primary_mixed | -0.233 (0.0594, 0%) | GSE205335+GSE241934_IIT+GSE291670+GSE131907 |
| 4 | 76 | primary_mixed | -0.237 (0.0595, 3%) | GSE291670+GSE253013+GSE131907+GSE325414 |
| 4 | 62 | primary_mixed | -0.258 (0.0622, 0%) | GSE241934_IIT+GSE291670+GSE253013+GSE131907 |
| 5 | 87 | primary_mixed | -0.216 (0.0623, 0%) | GSE241934_IIT+GSE291670+GSE253013+GSE131907+GSE325414 |
| 5 | 87 | primary_mixed | -0.216 (0.0628, 0%) | GSE207422+GSE205335+GSE241934_IIT+GSE291670+GSE131907 |

Primary-grid CLDN4-negative subset pools (k≥2): 237 of 247 enumerated k≥2 subsets.

## Single-cohort CLDN4 ρ < 0 vs T/NK (n≥4)

Every available (cohort × malignant def × T/NK × score) with CLDN4 ρ<0.

| cohort | malig | immune | score | n | CLDN4 ρ (p) | source |
|---|---|---|---|---:|---|---|
| GSE131907 | author_malig | tnk | pct | 21 | -0.522 (0.0152) | tnk_extract |
| GSE291670 | marker_malig | tnk | pct | 6 | -0.886 (0.0188) | tnk_extract |
| GSE131907 | author_epi | tnk | pct | 36 | -0.365 (0.0286) | tnk_extract |
| GSE291670 | marker_malig | tnk | mean | 6 | -0.829 (0.0416) | tnk_extract |
| GSE291670 | marker_epi | tnk | mean | 6 | -0.829 (0.0416) | tnk_extract |
| GSE205335 | author_malig | tnk | pct | 22 | -0.435 (0.0429) | tnk_extract |
| GSE131907 | author_malig | tnk | mean | 21 | -0.396 (0.0755) | tnk_extract |
| GSE253013 | marker_malig | tnk | mean_cp10k | 9 | -0.550 (0.125) | tnk_extract |
| GSE253013 | marker_malig | tnk | pct | 9 | -0.533 (0.139) | tnk_extract |
| GSE131907 | author_epi | tnk | mean | 36 | -0.210 (0.218) | tnk_extract |
| GSE131907_tLung | author_epi | tnk | pct | 11 | -0.400 (0.223) | tnk_extract |
| GSE205335 | author_malig | tnk | mean | 22 | -0.200 (0.371) | tnk_extract |
| GSE253013 | marker_malig | tnk | mean | 9 | -0.333 (0.381) | tnk_extract |
| GSE267108 | marker_epi | tnk | mean | 8 | -0.357 (0.385) | leftover |
| GSE205335_NSCLC | author_malig | tnk | pct | 17 | -0.206 (0.428) | tnk_extract |
| GSE325414 | author_malig | tnk | mean | 25 | -0.119 (0.57) | tnk_extract |
| GSE253013 | marker_epi | tnk | mean | 9 | -0.200 (0.606) | tnk_extract |
| GSE207422 | author_DRMref | tnk | mean | 12 | -0.091 (0.779) | a3_given_table |
| GSE241934_IIT | author_epi | tnk | mean | 11 | -0.064 (0.853) | tnk_extract |
| GSE253013 | marker_epi | tnk | pct | 9 | -0.067 (0.865) | tnk_extract |
| GSE241934_Real | author_epi | tnk | pct | 24 | -0.034 (0.875) | tnk_extract |
| GSE131907_tLung | author_epi | tnk | mean | 11 | -0.045 (0.894) | tnk_extract |
| GSE205335_NSCLC | author_malig | tnk | mean | 17 | -0.022 (0.933) | tnk_extract |

Count: 23 CLDN4-negative T/NK singles (of 31 T/NK singles with n≥4).

## CLDN4 ρ < 0 vs CXCL13+ (n≥4)

CXCL13+ T fraction and CXCL13 mean from the PR #274 extract. The TACSTD2
UCell-vs-CXCL13+ grid in PR #271 is TACSTD2-primary and is not re-audited.

| cohort | malig | immune | score | n | CLDN4 ρ (p) |
|---|---|---|---|---:|---|
| GSE241934_IIT | author_epi | cxcl13_mean | mean | 11 | -0.727 (0.0112) |
| GSE148071 | malig_or_epi | cxcl13pos | mean | 36 | -0.396 (0.0169) |
| GSE241934_IIT | author_epi | cxcl13pos | mean | 11 | -0.655 (0.0289) |
| GSE253013 | malig_or_epi | cxcl13_mean | mean | 9 | -0.667 (0.0499) |
| GSE207422 | author_epi_tls | cxcl13pos | mean | 15 | -0.511 (0.0517) |
| GSE207422 | author_epi_tls | cxcl13_mean | mean | 15 | -0.486 (0.0664) |
| GSE148071 | malig_or_epi | cxcl13_mean | mean | 41 | -0.270 (0.0883) |
| GSE253013 | malig_or_epi | cxcl13pos | mean | 9 | -0.400 (0.286) |
| GSE131907 | malig_or_epi | cxcl13_mean | mean | 32 | -0.097 (0.597) |

CXCL13+ / CXCL13-mean subset pools with CLDN4 ρ<0, N≥15, lowest p first:

| k | N | family | CLDN4 ρ (p, I²) | cohorts |
|---:|---:|---|---|---|
| 3 | 60 | tls/cxcl13pos/mean | -0.425 (0.00121, 0%) | GSE148071+GSE207422+GSE253013 |
| 2 | 51 | tls/cxcl13pos/mean | -0.428 (0.00216, 0%) | GSE148071+GSE207422 |
| 3 | 65 | tls/cxcl13_mean/mean | -0.369 (0.00374, 0%) | GSE148071+GSE207422+GSE253013 |
| 2 | 24 | tls/cxcl13_mean/mean | -0.552 (0.00833, 0%) | GSE207422+GSE253013 |
| 2 | 45 | tls/cxcl13pos/mean | -0.396 (0.00882, 0%) | GSE148071+GSE253013 |
| 4 | 97 | tls/cxcl13_mean/mean | -0.289 (0.0117, 11%) | GSE131907+GSE148071+GSE207422+GSE253013 |
| 2 | 56 | tls/cxcl13_mean/mean | -0.325 (0.0171, 0%) | GSE148071+GSE207422 |
| 3 | 88 | tls/cxcl13_mean/mean | -0.244 (0.0267, 0%) | GSE131907+GSE148071+GSE207422 |
| 2 | 24 | tls/cxcl13pos/mean | -0.475 (0.0283, 0%) | GSE207422+GSE253013 |
| 3 | 82 | tls/cxcl13_mean/mean | -0.255 (0.0607, 21%) | GSE131907+GSE148071+GSE253013 |
| 4 | 92 | tls/cxcl13pos/mean | -0.282 (0.0609, 40%) | GSE131907+GSE148071+GSE207422+GSE253013 |
| 3 | 56 | tls/cxcl13_mean/mean | -0.356 (0.0731, 40%) | GSE131907+GSE207422+GSE253013 |

## CLDN4 ρ < 0 vs B (n≥4)

B fraction from the PR #274 extract plus B/plasma columns already on the
T/NK tables. PR #274 TACSTD2-vs-B RE meta is taken as given.

| cohort | malig | immune | score | n | CLDN4 ρ (p) |
|---|---|---|---|---:|---|
| GSE131907 | author_malig | b | pct | 21 | -0.600 (0.00404) |
| GSE131907 | author_epi | b | pct | 36 | -0.409 (0.0133) |
| GSE131907 | author_malig | b | mean | 21 | -0.465 (0.0337) |
| GSE241934_IIT | author_epi | b | mean | 11 | -0.609 (0.0467) |
| GSE131907 | author_epi | b | mean | 36 | -0.300 (0.0751) |
| GSE131907 | malig_or_epi | b | mean | 32 | -0.303 (0.0921) |
| GSE154826 | malig_or_epi | b | mean | 30 | -0.252 (0.179) |
| GSE207422 | author_epi_tls | b | mean | 15 | -0.239 (0.39) |
| GSE205335 | author_malig | b_plasma | pct | 22 | -0.159 (0.481) |
| GSE205335_NSCLC | author_malig | b_plasma | pct | 17 | -0.081 (0.758) |

B-fraction subset pools with CLDN4 ρ<0, N≥15, lowest p first:

| k | N | family | CLDN4 ρ (p, I²) | cohorts |
|---:|---:|---|---|---|
| 2 | 32 | author/b/mean | -0.513 (0.00388, 0%) | GSE131907+GSE241934_IIT |
| 3 | 77 | tls/b/mean | -0.272 (0.0215, 0%) | GSE131907+GSE154826+GSE207422 |
| 2 | 62 | tls/b/mean | -0.279 (0.0322, 0%) | GSE131907+GSE154826 |
| 4 | 86 | tls/b/mean | -0.230 (0.0441, 0%) | GSE131907+GSE154826+GSE207422+GSE253013 |
| 2 | 47 | tls/b/mean | -0.284 (0.0611, 0%) | GSE131907+GSE207422 |
| 3 | 71 | tls/b/mean | -0.228 (0.0677, 0%) | GSE131907+GSE154826+GSE253013 |
| 3 | 54 | author/b/mean | -0.343 (0.0897, 47%) | GSE131907+GSE205335+GSE241934_IIT |
| 2 | 43 | author/b/pct | -0.400 (0.112, 62%) | GSE131907+GSE205335 |
| 2 | 45 | tls/b/mean | -0.248 (0.113, 0%) | GSE154826+GSE207422 |
| 3 | 56 | tls/b/mean | -0.217 (0.131, 0%) | GSE131907+GSE207422+GSE253013 |

## CLDN4 vs cyto-high T / T/NK cytotoxicity

Continuous T/NK cytotoxicity score (GZMB/PRF1/GNLY/NKG7 mean log1p) from
the PR #260 extract and the given GSE207422 table. This is the public
cyto-high T score that already exists; PR #271 cyto-high T fractions are
TACSTD2-primary (UCell module) and are not re-audited.

| cohort | malig | immune | score | n | CLDN4 ρ (p) |
|---|---|---|---|---:|---|
| GSE205335 | author_malig | cyto | mean | 22 | -0.246 (0.271) |
| GSE241934_IIT | author_epi | cyto | mean | 11 | -0.245 (0.467) |
| GSE207422 | author_DRMref | cyto | mean | 12 | -0.091 (0.779) |
| GSE207422 | marker_malig | cyto | mean | 7 | +0.000 (1) |
| GSE241934_Real | author_epi | cyto | mean | 29 | +0.153 (0.428) |
| GSE233203 | marker_or_author_malig | cyto | mean | 6 | +0.714 (0.111) |
| GSE291670 | marker_malig | cyto | mean | 6 | +0.771 (0.0724) |

Cytotoxicity subset pools with CLDN4 ρ<0, N≥15:

| k | N | family | CLDN4 ρ (p, I²) | cohorts |
|---:|---:|---|---|---|
| 2 | 33 | author/cyto/mean | -0.246 (0.193, 0%) | GSE205335+GSE241934_IIT |
| 3 | 45 | author/cyto/mean | -0.208 (0.206, 0%) | GSE205335+GSE207422+GSE241934_IIT |
| 2 | 34 | author/cyto/mean | -0.197 (0.291, 0%) | GSE205335+GSE207422 |
| 2 | 23 | author/cyto/mean | -0.165 (0.493, 0%) | GSE207422+GSE241934_IIT |
| 4 | 74 | author/cyto/mean | -0.058 (0.65, 0%) | GSE205335+GSE207422+GSE241934_IIT+GSE241934_Real |
| 3 | 62 | author/cyto/mean | -0.057 (0.694, 8%) | GSE205335+GSE241934_IIT+GSE241934_Real |
| 3 | 63 | author/cyto/mean | -0.029 (0.831, 0%) | GSE205335+GSE207422+GSE241934_Real |
| 2 | 51 | author/cyto/mean | -0.031 (0.879, 45%) | GSE205335+GSE241934_Real |

## Extra leftover n (not in the primary 8-unit grid)

| cohort | malig | immune | n | CLDN4 ρ (p) | note |
|---|---|---|---:|---|---|
| GSE267108 | marker_epi | tnk | 8 | -0.357 (0.385) | leftover PR #280; treatment-naive LUAD; extra n |
| GSE274595 | marker_epi | tnk | 7 | +0.250 (0.589) | leftover PR #280; surgical, not ICI; extra n |
| EMTAB13526 | author_epi | tnk | 11 | +0.327 (0.326) | leftover PR #280 / #259; Cvejic atlas, not ICI; extra n |

## Aligned-family CLDN4-negative pools (k≥2, N≥15, lowest CLDN4 p)

Families hold malignant / immune / score buckets fixed. Top 15.

| k | N | family | CLDN4 ρ (p, I²) | cohorts |
|---:|---:|---|---|---|
| 3 | 60 | tls/cxcl13pos/mean | -0.425 (0.00121, 0%) | GSE148071+GSE207422+GSE253013 |
| 2 | 43 | author/tnk/pct | -0.479 (0.00152, 0%) | GSE131907+GSE205335 |
| 2 | 43 | author/cd8/pct | -0.472 (0.00181, 0%) | GSE131907+GSE205335 |
| 2 | 51 | tls/cxcl13pos/mean | -0.428 (0.00216, 0%) | GSE148071+GSE207422 |
| 3 | 65 | tls/cxcl13_mean/mean | -0.369 (0.00374, 0%) | GSE148071+GSE207422+GSE253013 |
| 2 | 32 | author/b/mean | -0.513 (0.00388, 0%) | GSE131907+GSE241934_IIT |
| 2 | 24 | tls/cxcl13_mean/mean | -0.552 (0.00833, 0%) | GSE207422+GSE253013 |
| 2 | 45 | tls/cxcl13pos/mean | -0.396 (0.00882, 0%) | GSE148071+GSE253013 |
| 3 | 54 | author/tnk/pct | -0.389 (0.011, 13%) | GSE131907+GSE205335+GSE241934_IIT |
| 4 | 97 | tls/cxcl13_mean/mean | -0.289 (0.0117, 11%) | GSE131907+GSE148071+GSE207422+GSE253013 |
| 2 | 56 | tls/cxcl13_mean/mean | -0.325 (0.0171, 0%) | GSE148071+GSE207422 |
| 3 | 77 | tls/b/mean | -0.272 (0.0215, 0%) | GSE131907+GSE154826+GSE207422 |
| 2 | 15 | marker/tnk/pct | -0.714 (0.0217, 23%) | GSE253013+GSE291670 |
| 3 | 88 | tls/cxcl13_mean/mean | -0.244 (0.0267, 0%) | GSE131907+GSE148071+GSE207422 |
| 2 | 24 | tls/cxcl13pos/mean | -0.475 (0.0283, 0%) | GSE207422+GSE253013 |

## Given TACSTD2-primary combos (not re-audited)

From PR #279 (TACSTD2-first): GSE207422+GSE291670+GSE253013 mean/T/NK
N=27 TACSTD2 ρ=−0.585 p=0.00452; marker-malignant %pos GSE253013+GSE291670
N=15 TACSTD2 ρ=−0.736 p=0.00766. From PR #271: GSE207422 UCell
{TACSTD2, CLDN4, EPCAM} vs CXCL13+ ρ=−0.657 p=0.020 (n=12). From PR #274:
TACSTD2 vs B RE ρ=−0.110 p=0.47 (k=7, n=162). Those rows stay there.

Figures: `figures/primary_full_pool_members_forest.png`,
`figures/primary_cldn4_negative_members_forest.png`,
`figures/cldn4_tnk_neg_singles_forest.png`, `figures/cldn4_tnk_subset_bars.png`,
`figures/best_marker_malig_*_GSE253013_GSE291670.png`,
`figures/best_author_pct_GSE131907_GSE205335.png`,
`figures/best_cxcl13pos_GSE148071_GSE207422_GSE253013.png`,
`figures/best_GSE131907_author_malig_pct.png`,
`figures/cldn4_cxcl13pos_neg_forest.png`, `figures/cldn4_B_neg_forest.png`,
`figures/cldn4_cyto_forest.png`, `figures/scatter_*.png`.

Reproduce: `python3 methods/scrna_cldn4_combo/combinatorial_search.py`

