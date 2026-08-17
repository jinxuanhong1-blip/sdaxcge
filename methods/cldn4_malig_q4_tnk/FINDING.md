# Malignant CLDN4-only: Q4 vs Q1 and combinatorial search vs T/NK / CXCL13+ / B

ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4 score. Patient is
the unit (GSE131907 T/NK extract is sample-level; that is stated on those
rows). Numbers are computed Spearman / DerSimonian–Laird / Stouffer and
Mann–Whitney Q4 vs Q1 (rank-biserial *r*). p-values are descriptive.

Existing processed malignant scores only (PR #279 T/NK, PR #274 TLS/B/
CXCL13+). GSE207422 A3 TACSTD2 is taken as given. CLDN4 on the same
locked 12-patient DRMref table is one honest row. GSE253013 uses the
existing 9-patient tumor extract; the 9 GB GEO RDS was not downloaded.

**Dropped all-epithelial:** GSE241934 IIT/Real residual Epi, GSE131907 author_epi (n=36), GSE253013/GSE291670 marker_epi, leftover epi extracts (GSE267108 / GSE274595 / E-MTAB-13526), TLS epithelial/gate rows (GSE154826 gate; GSE241934_IIT epithelial).

PR #290 mixed epithelial + malignant. This recut keeps malignant
definitions only. TACSTD2-primary combos in PR #279 / #271 / #274 stay
there and are not re-ranked as the CLDN4 answer.

## Full-pool row (one row, not the answer)

Primary **6-unit malignant-only** grid, mean log1p vs T/NK. PR #279's
8-unit mixed grid minus the two all-epithelial GSE241934 units, and
GSE131907 switched from author_epi (n=36) to author_malig (n=21).
Q4 vs Q1 full-pool uses only **poolable** members (n≥8 and both tails
≥3). GSE291670 (n=6, 2/2) and GSE253013 (n=9, 3/2) stay as flagged
thin singles and are not Fisher-z pooled — |r|=1 on n=2 vs n=2 is
complete separation, not a meta-analytic effect.

| analysis | k | N | CLDN4 effect (p, I²) | cohorts |
|---|---:|---:|---|---|
| Spearman ρ | 6 | 95 | -0.260 (0.0195, 0%) | GSE207422+GSE205335+GSE291670+GSE253013+GSE131907+GSE325414 |
| Q4 vs Q1 r | 4 | 42 compared | -0.405 (0.0188, 0%) | GSE207422+GSE205335+GSE131907+GSE325414 |

The search below is the answer: combinations where **CLDN4 is negative**
vs T/NK, CXCL13+, or B, with honest n and p.

## Combinations that recover CLDN4-negative (honest n / p)

Highlighted cuts. Ranked by CLDN4, not TACSTD2. Q4 vs Q1 uses
rank-biserial *r* on the quartile tails (n_compared = n_Q1 + n_Q4).

- spearman · GSE131907+GSE205335 · author/tnk/pct · k=2 · N=43 · CLDN4 ρ=-0.479 p=0.00152 I²=0%
- spearman · GSE253013+GSE291670 · marker/tnk/pct · k=2 · N=15 · CLDN4 ρ=-0.714 p=0.0217 I²=23%
- spearman · GSE253013+GSE291670 · marker/tnk/mean · k=2 · N=15 · CLDN4 ρ=-0.582 p=0.102 I²=29%
- spearman · GSE291670 · primary_malig_mixed · k=1 · N=6 · CLDN4 ρ=-0.829 p=0.0416 I²=0%
- spearman · GSE291670+GSE253013 · primary_malig_mixed · k=2 · N=15 · CLDN4 ρ=-0.582 p=0.102 I²=29%
- spearman · GSE205335+GSE291670+GSE253013+GSE131907 · primary_malig_mixed · k=4 · N=58 · CLDN4 ρ=-0.354 p=0.012 I²=0%
- spearman · GSE207422+GSE291670+GSE253013 · primary_malig_mixed · k=3 · N=27 · CLDN4 ρ=-0.378 p=0.157 I²=26%
- spearman · GSE148071+GSE207422+GSE253013 · tls_malig/cxcl13pos/mean · k=3 · N=45 · CLDN4 ρ=-0.371 p=0.0194 I²=0%
- spearman · GSE148071+GSE207422 · tls_malig/cxcl13pos/mean · k=2 · N=39 · CLDN4 ρ=-0.381 p=0.0212 I²=0%
- spearman · GSE207422+GSE253013 · tls_malig/cxcl13_mean/mean · k=2 · N=14 · CLDN4 ρ=-0.323 p=0.344 I²=0%
- spearman · GSE131907+GSE205335 · author/b/pct · k=2 · N=43 · CLDN4 ρ=-0.400 p=0.112 I²=62%
- q4q1 · GSE131907+GSE205335 · author/tnk/pct · k=2 · N=23 · n_Q1=12 n_Q4=11 · CLDN4 r=-0.705 p=0.000301 I²=0%
- q4q1 · GSE205335+GSE131907+GSE325414 · primary_malig_mixed_poolable · k=3 · N=36 · n_Q1=19 n_Q4=17 · CLDN4 r=-0.412 p=0.0228 I²=0%
- q4q1 · GSE207422+GSE205335+GSE131907+GSE325414 · primary_malig_mixed_poolable · k=4 · N=42 · n_Q1=22 n_Q4=20 · CLDN4 r=-0.405 p=0.0188 I²=0%
- q4q1 · GSE131907+GSE205335 · author/b/pct · k=2 · N=23 · n_Q1=12 n_Q4=11 · CLDN4 r=-0.611 p=0.24 I²=84%
- q4q1 single · GSE131907 · author_malig/b/pct · n=21 · n_Q1=6 n_Q4=5 · CLDN4 r=-0.867 p=0.0173
- q4q1 single · GSE131907 · author_malig/b/mean · n=21 · n_Q1=6 n_Q4=5 · CLDN4 r=-0.800 p=0.0303
- q4q1 single · GSE148071 · tls_malignant/cxcl13pos/mean · n=31 · n_Q1=8 n_Q4=8 · CLDN4 r=-0.609 p=0.0444
- q4q1 single · GSE205335 · author_malig/tnk/pct · n=22 · n_Q1=6 n_Q4=6 · CLDN4 r=-0.778 p=0.026

## Primary-grid members (malignant-only; the 6 units behind the full-pool row)

| cohort | malig | immune | score | n | unit | Spearman ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | ρ<0 |
|---|---|---|---|---:|---|---|---|---|
| GSE207422 | author_DRMref | tnk | mean | 12 | patient | -0.091 (0.779) | -0.333 (0.7; 3/3) | yes |
| GSE205335 | author_malig | tnk | mean | 22 | patient | -0.200 (0.371) | -0.556 (0.132; 6/6) | yes |
| GSE291670 | marker_malig | tnk | mean | 6 | patient | -0.829 (0.0416) | -1.000 (0.333; 2/2) thin | yes |
| GSE253013 | marker_malig | tnk | mean | 9 | patient | -0.333 (0.381) | -1.000 (0.2; 3/2) thin | yes |
| GSE131907 | author_malig | tnk | mean | 21 | sample | -0.396 (0.0755) | -0.533 (0.177; 6/5) | yes |
| GSE325414 | author_malig | tnk | mean | 25 | donor | -0.119 (0.57) | -0.143 (0.731; 7/6) | yes |

## Primary-grid Spearman subset pools with CLDN4 ρ < 0 (k≥2)

All 2^6−1 = 63 nonempty subsets of the 6 malignant-only primary units were enumerated.
Rows below are CLDN4-negative subsets with N≥20, lowest RE p first
(top 15). Full-pool is not repeated here.

| k | N | family | CLDN4 ρ (p, I²) | cohorts |
|---:|---:|---|---|---|
| 3 | 36 | primary_malig_mixed | -0.453 (0.0112, 0%) | GSE291670+GSE253013+GSE131907 |
| 4 | 58 | primary_malig_mixed | -0.354 (0.012, 0%) | GSE205335+GSE291670+GSE253013+GSE131907 |
| 5 | 70 | primary_malig_mixed | -0.314 (0.0161, 0%) | GSE207422+GSE205335+GSE291670+GSE253013+GSE131907 |
| 5 | 83 | primary_malig_mixed | -0.281 (0.0171, 0%) | GSE205335+GSE291670+GSE253013+GSE131907+GSE325414 |
| 4 | 48 | primary_malig_mixed | -0.370 (0.0197, 0%) | GSE207422+GSE291670+GSE253013+GSE131907 |
| 4 | 61 | primary_malig_mixed | -0.313 (0.0283, 5%) | GSE207422+GSE205335+GSE291670+GSE131907 |
| 5 | 86 | primary_malig_mixed | -0.254 (0.0289, 0%) | GSE207422+GSE205335+GSE291670+GSE131907+GSE325414 |
| 5 | 73 | primary_malig_mixed | -0.279 (0.029, 0%) | GSE207422+GSE291670+GSE253013+GSE131907+GSE325414 |
| 4 | 61 | primary_malig_mixed | -0.321 (0.0306, 9%) | GSE291670+GSE253013+GSE131907+GSE325414 |
| 3 | 49 | primary_malig_mixed | -0.376 (0.0356, 22%) | GSE205335+GSE291670+GSE131907 |
| 4 | 74 | primary_malig_mixed | -0.285 (0.0357, 14%) | GSE205335+GSE291670+GSE131907+GSE325414 |
| 3 | 52 | primary_malig_mixed | -0.304 (0.0398, 0%) | GSE205335+GSE253013+GSE131907 |
| 4 | 77 | primary_malig_mixed | -0.243 (0.0456, 0%) | GSE205335+GSE253013+GSE131907+GSE325414 |
| 4 | 64 | primary_malig_mixed | -0.268 (0.0473, 0%) | GSE207422+GSE205335+GSE253013+GSE131907 |
| 5 | 89 | primary_malig_mixed | -0.225 (0.0489, 0%) | GSE207422+GSE205335+GSE253013+GSE131907+GSE325414 |

Primary-grid CLDN4-negative Spearman subset pools (k≥2): 57 of 57 enumerated k≥2 subsets.

## Primary-grid Q4 vs Q1 subset pools with r < 0 (k≥2)

Poolable primary members only (drop thin GSE291670 / GSE253013).
Effect is pooled rank-biserial *r* on quartile tails. N is n_Q1 + n_Q4
(compared patients), not the full cohort n. Thin singles stay in the
member table above and are not pooled.

| k | N_compared | family | CLDN4 r (p, I²) | cohorts |
|---:|---:|---|---|---|
| 3 | 29 | primary_malig_mixed_poolable | -0.517 (0.0106, 0%) | GSE207422+GSE205335+GSE131907 |
| 2 | 23 | primary_malig_mixed_poolable | -0.545 (0.0117, 0%) | GSE205335+GSE131907 |
| 3 | 36 | primary_malig_mixed_poolable | -0.412 (0.0228, 0%) | GSE205335+GSE131907+GSE325414 |
| 2 | 18 | primary_malig_mixed_poolable | -0.505 (0.0539, 0%) | GSE207422+GSE205335 |
| 2 | 17 | primary_malig_mixed_poolable | -0.483 (0.0804, 0%) | GSE207422+GSE131907 |
| 3 | 31 | primary_malig_mixed_poolable | -0.353 (0.0836, 0%) | GSE207422+GSE205335+GSE325414 |
| 3 | 30 | primary_malig_mixed_poolable | -0.332 (0.114, 0%) | GSE207422+GSE131907+GSE325414 |
| 2 | 25 | primary_malig_mixed_poolable | -0.357 (0.121, 9%) | GSE205335+GSE325414 |
| 2 | 24 | primary_malig_mixed_poolable | -0.331 (0.144, 0%) | GSE131907+GSE325414 |
| 2 | 19 | primary_malig_mixed_poolable | -0.188 (0.492, 0%) | GSE207422+GSE325414 |

## Single-cohort malignant CLDN4 vs T/NK (n≥4)

Every available (cohort × malignant def × T/NK × score). Epithelial
defs are not listed. Q4 vs Q1 is omitted when n<6 or quartiles collapse.

| cohort | malig | score | n | unit | ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | source |
|---|---|---|---:|---|---|---|---|
| GSE131907 | author_malig | pct | 21 | sample | -0.522 (0.0152) | -0.600 (0.126; 6/5) | tnk_extract |
| GSE291670 | marker_malig | pct | 6 | patient | -0.886 (0.0188) | -1.000 (0.333; 2/2) thin | tnk_extract |
| GSE291670 | marker_malig | mean | 6 | patient | -0.829 (0.0416) | -1.000 (0.333; 2/2) thin | tnk_extract |
| GSE205335 | author_malig | pct | 22 | patient | -0.435 (0.0429) | -0.778 (0.026; 6/6) | tnk_extract |
| GSE131907 | author_malig | mean | 21 | sample | -0.396 (0.0755) | -0.533 (0.177; 6/5) | tnk_extract |
| GSE253013 | marker_malig | mean_cp10k | 9 | patient | -0.550 (0.125) | -1.000 (0.2; 3/2) thin | tnk_extract |
| GSE253013 | marker_malig | pct | 9 | patient | -0.533 (0.139) | -1.000 (0.2; 3/2) thin | tnk_extract |
| GSE205335 | author_malig | mean | 22 | patient | -0.200 (0.371) | -0.556 (0.132; 6/6) | tnk_extract |
| GSE253013 | marker_malig | mean | 9 | patient | -0.333 (0.381) | -1.000 (0.2; 3/2) thin | tnk_extract |
| GSE205335_NSCLC | author_malig | pct | 17 | patient | -0.206 (0.428) | -0.300 (0.556; 5/4) | tnk_extract |
| GSE325414 | author_malig | mean | 25 | donor | -0.119 (0.57) | -0.143 (0.731; 7/6) | tnk_extract |
| GSE207422 | marker_malig | mean | 7 | patient | +0.250 (0.589) | +0.000 (1; 2/2) thin | tnk_extract |
| GSE207422 | author_DRMref | mean | 12 | patient | -0.091 (0.779) | -0.333 (0.7; 3/3) | a3_given_table |
| GSE207422 | marker_malig_nsclc | mean | 9 | patient | +0.033 (0.932) | +0.000 (1; 3/2) thin | tnk_extract |
| GSE207422 | marker_malig_nsclc | pct | 9 | patient | +0.033 (0.932) | +0.000 (1; 3/2) thin | tnk_extract |
| GSE205335_NSCLC | author_malig | mean | 17 | patient | -0.022 (0.933) | -0.400 (0.413; 5/4) | tnk_extract |
| GSE207422 | marker_malig | pct | 7 | patient | +0.000 (1) | +0.000 (1; 2/2) thin | tnk_extract |

Count: 13 CLDN4-negative T/NK Spearman singles (of 17 malignant T/NK singles with n≥4).

## Malignant-compartment CLDN4 vs CXCL13+ (n≥4)

PR #274 extract restricted to `compartment=malignant`. Epithelial and
gate rows (including GSE154826) are dropped. The TACSTD2 UCell-vs-
CXCL13+ grid in PR #271 is TACSTD2-primary and is not re-audited.

| cohort | malig | immune | score | n | ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) |
|---|---|---|---|---:|---|---|
| GSE253013 | tls_malignant | cxcl13_mean | pct | 6 | -0.886 (0.0188) | -1.000 (0.333; 2/2) thin |
| GSE148071 | tls_malignant | cxcl13pos | mean | 31 | -0.409 (0.0225) | -0.609 (0.0444; 8/8) |
| GSE148071 | tls_malignant | cxcl13_mean | pct | 35 | -0.343 (0.0439) | -0.395 (0.169; 9/9) |
| GSE148071 | tls_malignant | cxcl13_mean | mean | 35 | -0.281 (0.102) | -0.407 (0.153; 9/9) |
| GSE253013 | tls_malignant | cxcl13pos | pct | 6 | -0.600 (0.208) | -1.000 (0.333; 2/2) thin |
| GSE253013 | tls_malignant | cxcl13_mean | mean | 6 | -0.600 (0.208) | -0.500 (0.667; 2/2) thin |
| GSE148071 | tls_malignant | cxcl13pos | pct | 31 | -0.222 (0.23) | -0.297 (0.341; 8/8) |
| GSE241934_Real | tls_malignant | cxcl13pos | mean | 6 | +0.486 (0.329) | +1.000 (0.333; 2/2) thin |
| GSE241934_Real | tls_malignant | cxcl13pos | pct | 6 | +0.486 (0.329) | +1.000 (0.333; 2/2) thin |
| GSE241934_Real | tls_malignant | cxcl13_mean | mean | 6 | +0.429 (0.397) | +0.500 (0.667; 2/2) thin |
| GSE241934_Real | tls_malignant | cxcl13_mean | pct | 6 | +0.429 (0.397) | +0.500 (0.667; 2/2) thin |
| GSE131907 | tls_malignant | cxcl13_mean | pct | 31 | -0.109 (0.56) | +0.031 (0.959; 8/8) |
| GSE207422 | tls_malignant | cxcl13pos | pct | 8 | -0.238 (0.57) | -1.000 (0.333; 2/2) thin |
| GSE207422 | tls_malignant | cxcl13pos | mean | 8 | -0.214 (0.61) | +0.000 (1; 2/2) thin |
| GSE207422 | tls_malignant | cxcl13_mean | pct | 8 | -0.214 (0.61) | -1.000 (0.333; 2/2) thin |
| GSE131907 | tls_malignant | cxcl13pos | mean | 31 | +0.094 (0.617) | +0.344 (0.279; 8/8) |
| GSE253013 | tls_malignant | cxcl13pos | mean | 6 | -0.257 (0.623) | +0.000 (1; 2/2) thin |
| GSE131907 | tls_malignant | cxcl13_mean | mean | 31 | -0.063 (0.735) | +0.188 (0.574; 8/8) |
| GSE207422 | tls_malignant | cxcl13_mean | mean | 8 | -0.119 (0.779) | +0.000 (1; 2/2) thin |
| GSE131907 | tls_malignant | cxcl13pos | pct | 31 | +0.049 (0.794) | +0.250 (0.442; 8/8) |

CXCL13+ / CXCL13-mean Spearman pools with CLDN4 ρ<0, N≥15, lowest p first:

| k | N | family | CLDN4 ρ (p, I²) | cohorts |
|---:|---:|---|---|---|
| 3 | 45 | tls_malig/cxcl13pos/mean | -0.371 (0.0194, 0%) | GSE148071+GSE207422+GSE253013 |
| 2 | 37 | tls_malig/cxcl13pos/mean | -0.395 (0.0201, 0%) | GSE148071+GSE253013 |
| 2 | 39 | tls_malig/cxcl13pos/mean | -0.381 (0.0212, 0%) | GSE148071+GSE207422 |
| 2 | 43 | tls_malig/cxcl13_mean/pct | -0.326 (0.0396, 0%) | GSE148071+GSE207422 |
| 4 | 51 | tls_malig/cxcl13pos/mean | -0.308 (0.0464, 0%) | GSE148071+GSE207422+GSE241934_Real+GSE253013 |
| 3 | 74 | tls_malig/cxcl13_mean/pct | -0.235 (0.0534, 0%) | GSE131907+GSE148071+GSE207422 |
| 2 | 41 | tls_malig/cxcl13_mean/mean | -0.312 (0.0559, 0%) | GSE148071+GSE253013 |
| 3 | 49 | tls_malig/cxcl13_mean/mean | -0.289 (0.0597, 0%) | GSE148071+GSE207422+GSE253013 |
| 4 | 80 | tls_malig/cxcl13_mean/pct | -0.318 (0.0601, 37%) | GSE131907+GSE148071+GSE207422+GSE253013 |
| 2 | 66 | tls_malig/cxcl13_mean/pct | -0.237 (0.0614, 0%) | GSE131907+GSE148071 |

CXCL13+ Q4 vs Q1 pools with r<0, N_compared≥8, lowest p first:

| k | N_compared | family | CLDN4 r (p, I²) | cohorts |
|---:|---:|---|---|---|
| 2 | 34 | tls_malig/cxcl13_mean/pct | -0.202 (0.361, 29%) | GSE131907+GSE148071 |
| 2 | 34 | tls_malig/cxcl13_mean/mean | -0.129 (0.677, 63%) | GSE131907+GSE148071 |
| 2 | 32 | tls_malig/cxcl13pos/mean | -0.173 (0.743, 86%) | GSE131907+GSE148071 |
| 2 | 32 | tls_malig/cxcl13pos/pct | -0.025 (0.928, 51%) | GSE131907+GSE148071 |

## Malignant CLDN4 vs B (n≥4)

B fraction from the PR #274 malignant-compartment extract plus B/plasma
columns already on the T/NK malignant tables. PR #274 TACSTD2-vs-B RE
meta is taken as given.

| cohort | malig | immune | score | n | ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) |
|---|---|---|---|---:|---|---|
| GSE131907 | author_malig | b | pct | 21 | -0.600 (0.00404) | -0.867 (0.0173; 6/5) |
| GSE207422 | tls_malignant | b | pct | 8 | -0.786 (0.0208) | -1.000 (0.333; 2/2) thin |
| GSE131907 | author_malig | b | mean | 21 | -0.465 (0.0337) | -0.800 (0.0303; 6/5) |
| GSE207422 | tls_malignant | b | mean | 8 | -0.714 (0.0465) | -1.000 (0.333; 2/2) thin |
| GSE131907 | tls_malignant | b | pct | 31 | -0.346 (0.0569) | -0.344 (0.279; 8/8) |
| GSE131907 | tls_malignant | b | mean | 31 | -0.321 (0.0783) | -0.500 (0.105; 8/8) |
| GSE148071 | tls_malignant | b | mean | 35 | +0.235 (0.175) | +0.383 (0.183; 9/9) |
| GSE253013 | tls_malignant | b | mean | 6 | +0.600 (0.208) | +1.000 (0.333; 2/2) thin |
| GSE205335 | author_malig | b_plasma | pct | 22 | -0.159 (0.481) | -0.111 (0.818; 6/6) |
| GSE241934_Real | tls_malignant | b | mean | 6 | +0.314 (0.544) | +0.500 (0.667; 2/2) thin |
| GSE241934_Real | tls_malignant | b | pct | 6 | +0.314 (0.544) | +0.500 (0.667; 2/2) thin |
| GSE148071 | tls_malignant | b | pct | 35 | +0.099 (0.57) | -0.099 (0.757; 9/9) |
| GSE253013 | tls_malignant | b | pct | 6 | +0.257 (0.623) | +0.000 (1; 2/2) thin |
| GSE205335_NSCLC | author_malig | b_plasma | pct | 17 | -0.081 (0.758) | +0.000 (1; 5/4) |
| GSE205335_NSCLC | author_malig | b_plasma | mean | 17 | +0.078 (0.765) | +0.200 (0.73; 5/4) |
| GSE205335 | author_malig | b_plasma | mean | 22 | +0.001 (0.998) | +0.167 (0.699; 6/6) |

B-fraction Spearman pools with CLDN4 ρ<0, N≥15, lowest p first:

| k | N | family | CLDN4 ρ (p, I²) | cohorts |
|---:|---:|---|---|---|
| 2 | 39 | tls_malig/b/mean | -0.437 (0.0518, 26%) | GSE131907+GSE207422 |
| 2 | 39 | tls_malig/b/pct | -0.532 (0.072, 52%) | GSE131907+GSE207422 |
| 2 | 43 | author/b/pct | -0.400 (0.112, 62%) | GSE131907+GSE205335 |
| 3 | 45 | tls_malig/b/pct | -0.405 (0.137, 44%) | GSE131907+GSE207422+GSE253013 |
| 3 | 45 | tls_malig/b/mean | -0.347 (0.145, 30%) | GSE131907+GSE207422+GSE241934_Real |
| 2 | 37 | tls_malig/b/pct | -0.280 (0.152, 5%) | GSE131907+GSE253013 |
| 3 | 45 | tls_malig/b/pct | -0.395 (0.167, 47%) | GSE131907+GSE207422+GSE241934_Real |
| 3 | 74 | tls_malig/b/pct | -0.308 (0.23, 73%) | GSE131907+GSE148071+GSE207422 |
| 4 | 51 | tls_malig/b/pct | -0.302 (0.237, 40%) | GSE131907+GSE207422+GSE241934_Real+GSE253013 |
| 3 | 43 | tls_malig/b/pct | -0.216 (0.256, 6%) | GSE131907+GSE241934_Real+GSE253013 |

B-fraction Q4 vs Q1 pools with r<0, N_compared≥8, lowest p first:

| k | N_compared | family | CLDN4 r (p, I²) | cohorts |
|---:|---:|---|---|---|
| 2 | 23 | author/b/pct | -0.611 (0.24, 84%) | GSE131907+GSE205335 |
| 2 | 34 | tls_malig/b/pct | -0.216 (0.246, 0%) | GSE131907+GSE148071 |
| 2 | 23 | author/b/mean | -0.430 (0.468, 85%) | GSE131907+GSE205335 |
| 2 | 34 | tls_malig/b/mean | -0.068 (0.887, 84%) | GSE131907+GSE148071 |

## Aligned-family CLDN4-negative Spearman pools (k≥2, N≥15, lowest p)

Families hold malignant / immune / score buckets fixed. Epithelial
families are not enumerated. Top 15.

| k | N | family | CLDN4 ρ (p, I²) | cohorts |
|---:|---:|---|---|---|
| 2 | 43 | author/tnk/pct | -0.479 (0.00152, 0%) | GSE131907+GSE205335 |
| 2 | 43 | author/cd8/pct | -0.472 (0.00181, 0%) | GSE131907+GSE205335 |
| 3 | 45 | tls_malig/cxcl13pos/mean | -0.371 (0.0194, 0%) | GSE148071+GSE207422+GSE253013 |
| 2 | 37 | tls_malig/cxcl13pos/mean | -0.395 (0.0201, 0%) | GSE148071+GSE253013 |
| 2 | 39 | tls_malig/cxcl13pos/mean | -0.381 (0.0212, 0%) | GSE148071+GSE207422 |
| 2 | 15 | marker/tnk/pct | -0.714 (0.0217, 23%) | GSE253013+GSE291670 |
| 2 | 43 | tls_malig/cxcl13_mean/pct | -0.326 (0.0396, 0%) | GSE148071+GSE207422 |
| 4 | 51 | tls_malig/cxcl13pos/mean | -0.308 (0.0464, 0%) | GSE148071+GSE207422+GSE241934_Real+GSE253013 |
| 2 | 39 | tls_malig/b/mean | -0.437 (0.0518, 26%) | GSE131907+GSE207422 |
| 3 | 74 | tls_malig/cxcl13_mean/pct | -0.235 (0.0534, 0%) | GSE131907+GSE148071+GSE207422 |
| 2 | 41 | tls_malig/cxcl13_mean/mean | -0.312 (0.0559, 0%) | GSE148071+GSE253013 |
| 3 | 49 | tls_malig/cxcl13_mean/mean | -0.289 (0.0597, 0%) | GSE148071+GSE207422+GSE253013 |
| 4 | 80 | tls_malig/cxcl13_mean/pct | -0.318 (0.0601, 37%) | GSE131907+GSE148071+GSE207422+GSE253013 |
| 2 | 43 | author/tnk/mean | -0.299 (0.0608, 0%) | GSE131907+GSE205335 |
| 2 | 66 | tls_malig/cxcl13_mean/pct | -0.237 (0.0614, 0%) | GSE131907+GSE148071 |

## Aligned-family CLDN4-negative Q4 vs Q1 pools (k≥2, N_compared≥8, lowest p)

| k | N_compared | family | CLDN4 r (p, I²) | cohorts |
|---:|---:|---|---|---|
| 2 | 23 | author/tnk/pct | -0.705 (0.000301, 0%) | GSE131907+GSE205335 |
| 2 | 23 | author/cd8/pct | -0.667 (0.000907, 0%) | GSE131907+GSE205335 |
| 3 | 29 | author/tnk/mean | -0.517 (0.0106, 0%) | GSE131907+GSE205335+GSE207422 |
| 2 | 23 | author/tnk/mean | -0.545 (0.0117, 0%) | GSE131907+GSE205335 |
| 4 | 42 | author/tnk/mean | -0.405 (0.0188, 0%) | GSE131907+GSE205335+GSE207422+GSE325414 |
| 3 | 36 | author/tnk/mean | -0.412 (0.0228, 0%) | GSE131907+GSE205335+GSE325414 |
| 2 | 18 | author/tnk/mean | -0.505 (0.0539, 0%) | GSE205335+GSE207422 |
| 2 | 23 | author/cd8/mean | -0.433 (0.0561, 0%) | GSE131907+GSE205335 |
| 2 | 17 | author/tnk/mean | -0.483 (0.0804, 0%) | GSE131907+GSE207422 |
| 3 | 31 | author/tnk/mean | -0.353 (0.0836, 0%) | GSE205335+GSE207422+GSE325414 |
| 3 | 30 | author/tnk/mean | -0.332 (0.114, 0%) | GSE131907+GSE207422+GSE325414 |
| 2 | 25 | author/tnk/mean | -0.357 (0.121, 9%) | GSE205335+GSE325414 |

## What was not done

- No dual-high TACSTD2×CLDN4 (or TACSTD2+CLDN4+EPCAM) score.
- No all-epithelial recut. GSE241934 has only residual Epi in the
  existing extract and is dropped, not imputed as malignant.
- GSE207422 wave2 inferCNV table has TACSTD2 only; CLDN4 was not invented.
- TACSTD2-primary combos in PR #279 / #271 / #274 are not re-audited.
- GSE253013 9 GB RDS was not downloaded.

Figures: `figures/primary_malig_tnk_*_forest.png`,
`figures/cldn4_tnk_*_forest.png`, `figures/cldn4_*_subset_bars.png`,
`figures/best_*.png`, `figures/q4q1_box_*.png`, `figures/scatter_*.png`.

Combo table: `tables/highlighted_combos.tsv`.

Reproduce: `python3 methods/cldn4_malig_q4_tnk/analyze.py`

