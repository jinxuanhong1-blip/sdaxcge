# Finding — CLDN4-only pair/triple enumeration vs same-patient T/NK

ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. No CellChat.
Bigger merge is **not** the answer. Every pair and every triple that
could actually be scored is a row. The full-pool is one row.
PR #320 GSE131907+GSE205335 Q4 (n=23, r=−0.705) is **one given row**
and was not re-audited.

Patient is the unit (GSE131907 T/NK extract is sample-level;
GSE123902 is donor-level). Public processed scRNA only (<2 GB).
p-values are descriptive (many subsets).

## Verdict — which combinations DIFFER

The listed set splits. **Negative and significant** pairs are those
built from GSE123902, GSE131907, GSE205335, and GSE189357.
**GSE148071 and GSE127465 do not** — they are null/positive as
singles (T/NK ρ ≈ +0.08 / +0.11) and they dilute every pair they enter.
GSE207422 is a weak negative (n=12, mean ρ=−0.09) and does not carry a pair.

Full-pool mean N=106 ρ=−0.25 is **weaker** than the best pairs
(GSE123902+GSE189357 ρ=−0.64; GSE123902+GSE131907 %pos ρ=−0.58;
given PR #320 GSE131907+GSE205335 Q4 r=−0.705). Bigger merge is not better.

High-end methods (not run here) should go to the differing pairs,
not the 7-cohort pool: GSE123902+GSE131907, GSE123902+GSE205335,
GSE123902+GSE189357, GSE131907+GSE189357, GSE189357+GSE205335,
and the given PR #320 pair. Triples that keep those members stay
negative; LOO that drops GSE123902 or GSE131907 often loses p<0.05.

## Scoreable vs skipped (listed accessions)

| accession | status | n | CLDN4 def | T/NK | note |
|---|---|---:|---|---|---|
| GSE123902 | scored | 13 | marker_malig (mean+%pos) | yes | Laughney 2020 tumor/met donors; marker-malignant; normals dropped |
| GSE127465 | scored | 7 | tisch_malig (mean+%pos) | yes | TISCH patient units, n_mal>=20 and n_tnk>=20 |
| GSE131907 | scored | 21 | author_malig (mean+%pos) | yes | author Malignant; tumor sites; n_mal>=20; sample-level (Kim 2020) |
| GSE148071 | scored | 22 | tisch_malig (mean+%pos) | yes | TISCH malignant (drop 3 epithelial-like); n_mal>=20 |
| GSE189357 | scored | 9 | marker_malig (mean+%pos) | yes | Zhu/Wang AIS–IAC; 9 patients; marker-malignant |
| GSE205335 | scored | 22 | author_malig (mean+%pos) | yes | author malignant, all subtypes |
| GSE207422 | scored | 12 | author_DRMref (mean only) | yes | locked A3 DRMref 12-patient table; CLDN4 mean only (%pos absent) |
| GSE179994 | skip | 0 | absent (T-only RDS) | T-only | no usable processed TME matrix (PR #420); not invented |
| GSE117570 | skip | 1 | present in 3/4 patients | 1 eligible | TISCH; Spearman-eligible n=1 (P2); skip |
| GSE146100 | skip | 1 | epithelial proxy | present | 1 patient / 3 nodules; patient is the unit; skip |
| GSE154826 | skip | 0 | gate only | no compact T/NK extract | CITE-seq; epithelium restricted; no <2GB malignant+T/NK matrix |
| GSE139555 | skip | 0 | absent (T-sorted) | T-sorted only | TISCH: 0 epithelial/malignant cells |
| GSE200563 | skip | 0 | spatial Visium | not scRNA | paired brain-met spatial RNA, not scRNA TME |
| GSE229353 | skip | 7 | absent (CD45+ only) | CD45+ immune | 7 NSCLC CD45+ libraries; no malignant epithelium |

Do not invent accessions. Combos use only the scored rows above.

## Full-pool row (one row, not the answer)

| score | k | N | ρ (p, I²) | Q4 vs Q1 r (p) |
|---|---:|---:|---|---|
| mean | 7 | 106 | -0.252 (0.028, I²=12%) | -0.396 (0.013) |
| pct | 6 | 94 | -0.391 (0.003, I²=27%) | -0.521 (0.002) |

## Given PR #320 row (not re-audited)

| combo | analysis | N | effect | p |
|---|---|---:|---|---|
| GSE131907+GSE205335 | Q4 vs Q1 author %pos | 23 (12/11) | r=-0.705 | 0.000301 |

## Singles (context; not the combo answer)

| cohort | score | n | ρ | p | Q4 r (p) |
|---|---|---:|---:|---:|---|
| GSE123902 | mean | 13 | -0.654 | 0.015 | n<16 |
| GSE189357 | mean | 9 | -0.517 | 0.154 | n<16 |
| GSE131907 | mean | 21 | -0.396 | 0.075 | -0.533 (0.177) |
| GSE205335 | mean | 22 | -0.200 | 0.371 | -0.556 (0.132) |
| GSE207422 | mean | 12 | -0.091 | 0.779 | n<16 |
| GSE148071 | mean | 22 | 0.080 | 0.725 | -0.056 (0.937) |
| GSE127465 | mean | 7 | 0.107 | 0.819 | n<16 |
| GSE123902 | pct | 13 | -0.659 | 0.014 | n<16 |
| GSE189357 | pct | 9 | -0.600 | 0.088 | n<16 |
| GSE131907 | pct | 21 | -0.522 | 0.015 | -0.600 (0.126) |
| GSE205335 | pct | 22 | -0.435 | 0.043 | -0.778 (0.026) |
| GSE148071 | pct | 22 | -0.015 | 0.946 | -0.056 (0.937) |
| GSE127465 | pct | 7 | 0.179 | 0.702 | n<16 |

## Every pair actually scored

Ranked by |ρ|. Q4 vs Q1 uses within-cohort CLDN4 ranks then MWU (n≥16).

| score | combo | k | N | ρ | p | I² | Q4 r | Q4 p | n_Q1/Q4 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| pct | GSE123902+GSE189357 | 2 | 22 | -0.638 | 0.003 | 0 | -1.000 | 0.003 | 7/5 |
| mean | GSE123902+GSE189357 | 2 | 22 | -0.606 | 0.005 | 0 | -1.000 | 0.003 | 7/5 |
| pct | GSE123902+GSE131907 | 2 | 34 | -0.575 | 0.001 | 0 | -0.700 | 0.012 | 10/8 |
| pct | GSE131907+GSE189357 | 2 | 30 | -0.542 | 0.003 | 0 | -0.619 | 0.042 | 9/7 |
| pct | GSE123902+GSE205335 | 2 | 35 | -0.522 | 0.002 | 0 | -0.802 | 0.005 | 9/9 |
| mean | GSE123902+GSE131907 | 2 | 34 | -0.500 | 0.004 | 0 | -0.650 | 0.021 | 10/8 |
| pct | GSE131907+GSE205335 | 2 | 43 | -0.479 | 0.002 | 0 | -0.702 | 0.006 | 11/11 |
| pct | GSE189357+GSE205335 | 2 | 31 | -0.478 | 0.009 | 0 | -0.750 | 0.010 | 8/8 |
| mean | GSE131907+GSE189357 | 2 | 30 | -0.428 | 0.025 | 0 | -0.524 | 0.091 | 9/7 |
| mean | GSE123902+GSE205335 | 2 | 35 | -0.423 | 0.115 | 54 | -0.630 | 0.027 | 9/9 |
| mean | GSE123902+GSE207422 | 2 | 25 | -0.417 | 0.198 | 56 | -0.714 | 0.035 | 7/6 |
| mean | GSE123902+GSE127465 | 2 | 20 | -0.398 | 0.334 | 56 | -0.520 | 0.222 | 5/5 |
| pct | GSE123902+GSE127465 | 2 | 20 | -0.365 | 0.425 | 63 | -0.520 | 0.222 | 5/5 |
| pct | GSE123902+GSE148071 | 2 | 35 | -0.357 | 0.335 | 75 | -0.333 | 0.251 | 9/9 |
| pct | GSE127465+GSE131907 | 2 | 28 | -0.316 | 0.360 | 47 | -0.429 | 0.209 | 7/7 |
| mean | GSE123902+GSE148071 | 2 | 35 | -0.313 | 0.452 | 79 | -0.358 | 0.216 | 9/9 |
| mean | GSE127465+GSE131907 | 2 | 28 | -0.312 | 0.129 | 0 | -0.388 | 0.259 | 7/7 |
| mean | GSE131907+GSE207422 | 2 | 33 | -0.300 | 0.108 | 0 | -0.444 | 0.139 | 9/8 |
| mean | GSE131907+GSE205335 | 2 | 43 | -0.299 | 0.061 | 0 | -0.521 | 0.042 | 11/11 |
| pct | GSE127465+GSE189357 | 2 | 16 | -0.295 | 0.484 | 45 | -0.250 | 0.686 | 4/4 |
| pct | GSE127465+GSE205335 | 2 | 29 | -0.287 | 0.300 | 28 | -0.571 | 0.072 | 8/7 |
| pct | GSE131907+GSE148071 | 2 | 43 | -0.286 | 0.296 | 66 | -0.405 | 0.115 | 11/11 |
| mean | GSE127465+GSE189357 | 2 | 16 | -0.285 | 0.380 | 10 | -0.250 | 0.686 | 4/4 |
| mean | GSE189357+GSE205335 | 2 | 31 | -0.284 | 0.145 | 0 | -0.562 | 0.065 | 8/8 |
| mean | GSE189357+GSE207422 | 2 | 21 | -0.276 | 0.272 | 0 | -0.667 | 0.082 | 6/5 |
| pct | GSE148071+GSE189357 | 2 | 31 | -0.264 | 0.411 | 52 | -0.250 | 0.442 | 8/8 |
| pct | GSE148071+GSE205335 | 2 | 44 | -0.236 | 0.286 | 48 | -0.350 | 0.176 | 12/10 |
| mean | GSE205335+GSE207422 | 2 | 34 | -0.166 | 0.376 | 0 | -0.481 | 0.093 | 9/9 |
| mean | GSE131907+GSE148071 | 2 | 43 | -0.165 | 0.504 | 57 | -0.240 | 0.358 | 11/11 |
| mean | GSE148071+GSE189357 | 2 | 31 | -0.157 | 0.614 | 48 | -0.250 | 0.442 | 8/8 |
| mean | GSE127465+GSE205335 | 2 | 29 | -0.148 | 0.474 | 0 | -0.357 | 0.281 | 8/7 |
| mean | GSE127465+GSE148071 | 2 | 29 | 0.084 | 0.685 | 0 | -0.143 | 0.694 | 8/7 |
| mean | GSE148071+GSE205335 | 2 | 44 | -0.062 | 0.704 | 0 | -0.300 | 0.249 | 12/10 |
| mean | GSE127465+GSE207422 | 2 | 19 | -0.030 | 0.914 | 0 | -0.200 | 0.690 | 5/5 |
| mean | GSE148071+GSE207422 | 2 | 34 | 0.025 | 0.895 | 0 | -0.160 | 0.596 | 9/9 |
| pct | GSE127465+GSE148071 | 2 | 29 | 0.019 | 0.928 | 0 | -0.036 | 0.955 | 8/7 |

## Every triple actually scored

| score | combo | k | N | ρ | p | I² | Q4 r | Q4 p | n_Q1/Q4 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| pct | GSE123902+GSE131907+GSE189357 | 3 | 43 | -0.580 | 0.000 | 0 | -0.723 | 0.004 | 13/10 |
| pct | GSE123902+GSE189357+GSE205335 | 3 | 44 | -0.536 | 0.000 | 0 | -0.802 | 0.002 | 11/11 |
| pct | GSE123902+GSE131907+GSE205335 | 3 | 56 | -0.522 | 0.000 | 0 | -0.735 | 0.001 | 14/14 |
| mean | GSE123902+GSE131907+GSE189357 | 3 | 43 | -0.503 | 0.001 | 0 | -0.662 | 0.008 | 13/10 |
| pct | GSE131907+GSE189357+GSE205335 | 3 | 52 | -0.497 | 0.000 | 0 | -0.692 | 0.003 | 13/13 |
| pct | GSE123902+GSE127465+GSE189357 | 3 | 29 | -0.487 | 0.052 | 29 | -0.619 | 0.042 | 9/7 |
| mean | GSE123902+GSE127465+GSE189357 | 3 | 29 | -0.483 | 0.029 | 12 | -0.619 | 0.042 | 9/7 |
| pct | GSE123902+GSE127465+GSE131907 | 3 | 41 | -0.483 | 0.016 | 27 | -0.550 | 0.032 | 12/10 |
| pct | GSE127465+GSE131907+GSE189357 | 3 | 37 | -0.450 | 0.018 | 9 | -0.455 | 0.095 | 11/9 |
| mean | GSE123902+GSE189357+GSE207422 | 3 | 34 | -0.448 | 0.026 | 14 | -0.800 | 0.003 | 10/8 |
| pct | GSE123902+GSE127465+GSE205335 | 3 | 42 | -0.439 | 0.028 | 26 | -0.669 | 0.009 | 11/11 |
| mean | GSE123902+GSE127465+GSE131907 | 3 | 41 | -0.431 | 0.020 | 15 | -0.517 | 0.044 | 12/10 |
| pct | GSE127465+GSE131907+GSE205335 | 3 | 50 | -0.424 | 0.004 | 0 | -0.609 | 0.009 | 13/13 |
| mean | GSE123902+GSE189357+GSE205335 | 3 | 44 | -0.420 | 0.017 | 14 | -0.653 | 0.010 | 11/11 |
| pct | GSE123902+GSE148071+GSE189357 | 3 | 44 | -0.418 | 0.109 | 58 | -0.438 | 0.088 | 11/11 |
| mean | GSE123902+GSE131907+GSE207422 | 3 | 46 | -0.411 | 0.013 | 12 | -0.580 | 0.018 | 13/11 |
| pct | GSE123902+GSE131907+GSE148071 | 3 | 56 | -0.407 | 0.064 | 59 | -0.490 | 0.029 | 14/14 |
| pct | GSE127465+GSE189357+GSE205335 | 3 | 38 | -0.400 | 0.022 | 0 | -0.600 | 0.026 | 10/10 |
| mean | GSE123902+GSE131907+GSE205335 | 3 | 56 | -0.391 | 0.007 | 9 | -0.571 | 0.011 | 14/14 |
| pct | GSE123902+GSE148071+GSE205335 | 3 | 57 | -0.370 | 0.073 | 54 | -0.456 | 0.043 | 15/13 |
| mean | GSE123902+GSE148071+GSE189357 | 3 | 44 | -0.362 | 0.203 | 64 | -0.455 | 0.076 | 11/11 |
| mean | GSE127465+GSE131907+GSE189357 | 3 | 37 | -0.360 | 0.046 | 0 | -0.394 | 0.149 | 11/9 |
| pct | GSE131907+GSE148071+GSE189357 | 3 | 52 | -0.359 | 0.087 | 47 | -0.444 | 0.058 | 13/13 |
| mean | GSE131907+GSE189357+GSE207422 | 3 | 42 | -0.343 | 0.040 | 0 | -0.500 | 0.052 | 12/10 |
| pct | GSE131907+GSE148071+GSE205335 | 3 | 65 | -0.337 | 0.042 | 40 | -0.475 | 0.023 | 17/15 |
| mean | GSE123902+GSE127465+GSE205335 | 3 | 42 | -0.333 | 0.138 | 36 | -0.521 | 0.042 | 11/11 |
| mean | GSE131907+GSE189357+GSE205335 | 3 | 52 | -0.332 | 0.024 | 0 | -0.503 | 0.031 | 13/13 |
| mean | GSE123902+GSE205335+GSE207422 | 3 | 47 | -0.328 | 0.086 | 30 | -0.569 | 0.019 | 12/12 |
| mean | GSE123902+GSE131907+GSE148071 | 3 | 56 | -0.327 | 0.161 | 62 | -0.378 | 0.094 | 14/14 |
| mean | GSE123902+GSE127465+GSE207422 | 3 | 32 | -0.312 | 0.242 | 39 | -0.469 | 0.130 | 8/8 |
| pct | GSE148071+GSE189357+GSE205335 | 3 | 53 | -0.312 | 0.093 | 33 | -0.405 | 0.085 | 14/12 |
| mean | GSE127465+GSE131907+GSE205335 | 3 | 50 | -0.261 | 0.087 | 0 | -0.456 | 0.051 | 13/13 |
| mean | GSE131907+GSE205335+GSE207422 | 3 | 55 | -0.260 | 0.072 | 0 | -0.480 | 0.033 | 14/14 |
| mean | GSE123902+GSE148071+GSE205335 | 3 | 57 | -0.253 | 0.260 | 59 | -0.405 | 0.072 | 15/13 |
| mean | GSE127465+GSE131907+GSE207422 | 3 | 40 | -0.250 | 0.154 | 0 | -0.360 | 0.186 | 10/10 |
| pct | GSE123902+GSE127465+GSE148071 | 3 | 42 | -0.239 | 0.405 | 57 | -0.289 | 0.264 | 11/11 |
| mean | GSE131907+GSE148071+GSE189357 | 3 | 52 | -0.238 | 0.224 | 37 | -0.290 | 0.218 | 13/13 |
| mean | GSE189357+GSE205335+GSE207422 | 3 | 43 | -0.234 | 0.164 | 0 | -0.521 | 0.042 | 11/11 |
| mean | GSE123902+GSE148071+GSE207422 | 3 | 47 | -0.234 | 0.367 | 60 | -0.361 | 0.141 | 12/12 |
| mean | GSE127465+GSE189357+GSE205335 | 3 | 38 | -0.232 | 0.203 | 0 | -0.440 | 0.104 | 10/10 |
| pct | GSE127465+GSE131907+GSE148071 | 3 | 50 | -0.213 | 0.348 | 46 | -0.349 | 0.137 | 13/13 |
| mean | GSE123902+GSE127465+GSE148071 | 3 | 42 | -0.212 | 0.487 | 62 | -0.289 | 0.264 | 11/11 |
| mean | GSE127465+GSE189357+GSE207422 | 3 | 28 | -0.198 | 0.381 | 0 | -0.347 | 0.318 | 7/7 |
| pct | GSE127465+GSE148071+GSE205335 | 3 | 51 | -0.187 | 0.301 | 22 | -0.321 | 0.173 | 14/12 |
| mean | GSE131907+GSE148071+GSE205335 | 3 | 65 | -0.175 | 0.218 | 14 | -0.357 | 0.089 | 17/15 |
| pct | GSE127465+GSE148071+GSE189357 | 3 | 38 | -0.149 | 0.505 | 21 | -0.200 | 0.473 | 10/10 |
| mean | GSE131907+GSE148071+GSE207422 | 3 | 55 | -0.148 | 0.355 | 14 | -0.265 | 0.241 | 14/14 |
| mean | GSE148071+GSE189357+GSE205335 | 3 | 53 | -0.134 | 0.389 | 5 | -0.345 | 0.143 | 14/12 |
| mean | GSE127465+GSE205335+GSE207422 | 3 | 41 | -0.132 | 0.452 | 0 | -0.345 | 0.193 | 11/10 |
| mean | GSE127465+GSE131907+GSE148071 | 3 | 50 | -0.130 | 0.477 | 22 | -0.219 | 0.356 | 13/13 |
| mean | GSE148071+GSE189357+GSE207422 | 3 | 43 | -0.080 | 0.639 | 0 | -0.289 | 0.264 | 11/11 |
| mean | GSE148071+GSE205335+GSE207422 | 3 | 56 | -0.067 | 0.644 | 0 | -0.303 | 0.182 | 15/13 |
| mean | GSE127465+GSE148071+GSE189357 | 3 | 38 | -0.055 | 0.773 | 3 | -0.180 | 0.521 | 10/10 |
| mean | GSE127465+GSE148071+GSE205335 | 3 | 51 | -0.046 | 0.768 | 0 | -0.274 | 0.247 | 14/12 |
| mean | GSE127465+GSE148071+GSE207422 | 3 | 41 | 0.035 | 0.842 | 0 | -0.182 | 0.504 | 11/10 |

## Leave-one-cohort-out (every triple)

Full LOO table: `tables/triple_loo.tsv` (one row per dropped member).
Below: LOO rows that **flip** the parent (parent ρ<0 p<0.05 becomes
p≥0.05 or ρ≥0 after the drop). That is the honest ‘which member carries it’ cut.

| parent triple | score | dropped | remaining N | ρ | p |
|---|---|---|---:|---:|---:|
| GSE123902+GSE127465+GSE131907 | mean | GSE123902 | 28 | -0.312 | 0.129 |
| GSE123902+GSE127465+GSE131907 | mean | GSE131907 | 20 | -0.398 | 0.334 |
| GSE123902+GSE127465+GSE131907 | pct | GSE123902 | 28 | -0.316 | 0.360 |
| GSE123902+GSE127465+GSE131907 | pct | GSE131907 | 20 | -0.365 | 0.425 |
| GSE123902+GSE127465+GSE189357 | mean | GSE123902 | 16 | -0.285 | 0.380 |
| GSE123902+GSE127465+GSE189357 | mean | GSE189357 | 20 | -0.398 | 0.334 |
| GSE123902+GSE127465+GSE205335 | pct | GSE123902 | 29 | -0.287 | 0.300 |
| GSE123902+GSE127465+GSE205335 | pct | GSE205335 | 20 | -0.365 | 0.425 |
| GSE123902+GSE131907+GSE205335 | mean | GSE123902 | 43 | -0.299 | 0.061 |
| GSE123902+GSE131907+GSE205335 | mean | GSE131907 | 35 | -0.423 | 0.115 |
| GSE123902+GSE131907+GSE207422 | mean | GSE123902 | 33 | -0.300 | 0.108 |
| GSE123902+GSE131907+GSE207422 | mean | GSE131907 | 25 | -0.417 | 0.198 |
| GSE123902+GSE189357+GSE205335 | mean | GSE123902 | 31 | -0.284 | 0.145 |
| GSE123902+GSE189357+GSE205335 | mean | GSE189357 | 35 | -0.423 | 0.115 |
| GSE123902+GSE189357+GSE207422 | mean | GSE123902 | 21 | -0.276 | 0.272 |
| GSE123902+GSE189357+GSE207422 | mean | GSE189357 | 25 | -0.417 | 0.198 |
| GSE127465+GSE131907+GSE189357 | mean | GSE131907 | 16 | -0.285 | 0.380 |
| GSE127465+GSE131907+GSE189357 | mean | GSE189357 | 28 | -0.312 | 0.129 |
| GSE127465+GSE131907+GSE189357 | pct | GSE131907 | 16 | -0.295 | 0.484 |
| GSE127465+GSE131907+GSE189357 | pct | GSE189357 | 28 | -0.316 | 0.360 |
| GSE127465+GSE131907+GSE205335 | pct | GSE131907 | 29 | -0.287 | 0.300 |
| GSE127465+GSE131907+GSE205335 | pct | GSE205335 | 28 | -0.316 | 0.360 |
| GSE127465+GSE189357+GSE205335 | pct | GSE189357 | 29 | -0.287 | 0.300 |
| GSE127465+GSE189357+GSE205335 | pct | GSE205335 | 16 | -0.295 | 0.484 |
| GSE131907+GSE148071+GSE205335 | pct | GSE131907 | 44 | -0.236 | 0.286 |
| GSE131907+GSE148071+GSE205335 | pct | GSE205335 | 43 | -0.286 | 0.296 |
| GSE131907+GSE189357+GSE205335 | mean | GSE131907 | 31 | -0.284 | 0.145 |
| GSE131907+GSE189357+GSE205335 | mean | GSE189357 | 43 | -0.299 | 0.061 |
| GSE131907+GSE189357+GSE207422 | mean | GSE131907 | 21 | -0.276 | 0.272 |
| GSE131907+GSE189357+GSE207422 | mean | GSE189357 | 33 | -0.300 | 0.108 |

## Which combinations DIFFER (high-end methods list)

A combo **differs** if Spearman ρ<0 and p<0.05, or Q4 vs Q1 r<0
and p<0.05 (n≥16). These — not the full-pool — are the ones that
should get high-end methods (CellChat / LIANA / Milo). **CellChat
was not run in this PR.**

| why | score | combo | N | ρ (p) | Q4 r (p) |
|---|---|---|---:|---|---|
| rho<0 p<0.05+Q4 drop | mean | GSE123902+GSE131907 | 34 | -0.500 (0.004) | -0.650 (0.021) |
| rho<0 p<0.05 | mean | GSE123902+GSE189357 | 22 | -0.606 (0.005) | -1.000 (0.003) |
| Q4 drop | mean | GSE123902+GSE205335 | 35 | -0.423 (0.115) | -0.630 (0.027) |
| Q4 drop | mean | GSE123902+GSE207422 | 25 | -0.417 (0.198) | -0.714 (0.035) |
| rho<0 p<0.05 | mean | GSE131907+GSE189357 | 30 | -0.428 (0.025) | -0.524 (0.091) |
| Q4 drop | mean | GSE131907+GSE205335 | 43 | -0.299 (0.061) | -0.521 (0.042) |
| rho<0 p<0.05+Q4 drop | pct | GSE123902+GSE131907 | 34 | -0.575 (0.001) | -0.700 (0.012) |
| rho<0 p<0.05 | pct | GSE123902+GSE189357 | 22 | -0.638 (0.003) | -1.000 (0.003) |
| rho<0 p<0.05+Q4 drop | pct | GSE123902+GSE205335 | 35 | -0.522 (0.002) | -0.802 (0.005) |
| rho<0 p<0.05+Q4 drop | pct | GSE131907+GSE189357 | 30 | -0.542 (0.003) | -0.619 (0.042) |
| rho<0 p<0.05+Q4 drop | pct | GSE131907+GSE205335 | 43 | -0.479 (0.002) | -0.702 (0.006) |
| rho<0 p<0.05+Q4 drop | pct | GSE189357+GSE205335 | 31 | -0.478 (0.009) | -0.750 (0.010) |
| rho<0 p<0.05+Q4 drop | mean | GSE123902+GSE127465+GSE131907 | 41 | -0.431 (0.020) | -0.517 (0.044) |
| rho<0 p<0.05+Q4 drop | mean | GSE123902+GSE127465+GSE189357 | 29 | -0.483 (0.029) | -0.619 (0.042) |
| Q4 drop | mean | GSE123902+GSE127465+GSE205335 | 42 | -0.333 (0.138) | -0.521 (0.042) |
| rho<0 p<0.05+Q4 drop | mean | GSE123902+GSE131907+GSE189357 | 43 | -0.503 (0.001) | -0.662 (0.008) |
| rho<0 p<0.05+Q4 drop | mean | GSE123902+GSE131907+GSE205335 | 56 | -0.391 (0.007) | -0.571 (0.011) |
| rho<0 p<0.05+Q4 drop | mean | GSE123902+GSE131907+GSE207422 | 46 | -0.411 (0.013) | -0.580 (0.018) |
| rho<0 p<0.05+Q4 drop | mean | GSE123902+GSE189357+GSE205335 | 44 | -0.420 (0.017) | -0.653 (0.010) |
| rho<0 p<0.05+Q4 drop | mean | GSE123902+GSE189357+GSE207422 | 34 | -0.448 (0.026) | -0.800 (0.003) |
| Q4 drop | mean | GSE123902+GSE205335+GSE207422 | 47 | -0.328 (0.086) | -0.569 (0.019) |
| rho<0 p<0.05 | mean | GSE127465+GSE131907+GSE189357 | 37 | -0.360 (0.046) | -0.394 (0.149) |
| rho<0 p<0.05+Q4 drop | mean | GSE131907+GSE189357+GSE205335 | 52 | -0.332 (0.024) | -0.503 (0.031) |
| rho<0 p<0.05 | mean | GSE131907+GSE189357+GSE207422 | 42 | -0.343 (0.040) | -0.500 (0.052) |
| Q4 drop | mean | GSE131907+GSE205335+GSE207422 | 55 | -0.260 (0.072) | -0.480 (0.033) |
| Q4 drop | mean | GSE189357+GSE205335+GSE207422 | 43 | -0.234 (0.164) | -0.521 (0.042) |
| rho<0 p<0.05+Q4 drop | pct | GSE123902+GSE127465+GSE131907 | 41 | -0.483 (0.016) | -0.550 (0.032) |
| Q4 drop | pct | GSE123902+GSE127465+GSE189357 | 29 | -0.487 (0.052) | -0.619 (0.042) |
| rho<0 p<0.05+Q4 drop | pct | GSE123902+GSE127465+GSE205335 | 42 | -0.439 (0.028) | -0.669 (0.009) |
| Q4 drop | pct | GSE123902+GSE131907+GSE148071 | 56 | -0.407 (0.064) | -0.490 (0.029) |
| rho<0 p<0.05+Q4 drop | pct | GSE123902+GSE131907+GSE189357 | 43 | -0.580 (0.000) | -0.723 (0.004) |
| rho<0 p<0.05+Q4 drop | pct | GSE123902+GSE131907+GSE205335 | 56 | -0.522 (0.000) | -0.735 (0.001) |
| Q4 drop | pct | GSE123902+GSE148071+GSE205335 | 57 | -0.370 (0.073) | -0.456 (0.043) |
| rho<0 p<0.05+Q4 drop | pct | GSE123902+GSE189357+GSE205335 | 44 | -0.536 (0.000) | -0.802 (0.002) |
| rho<0 p<0.05 | pct | GSE127465+GSE131907+GSE189357 | 37 | -0.450 (0.018) | -0.455 (0.095) |
| rho<0 p<0.05+Q4 drop | pct | GSE127465+GSE131907+GSE205335 | 50 | -0.424 (0.004) | -0.609 (0.009) |
| rho<0 p<0.05+Q4 drop | pct | GSE127465+GSE189357+GSE205335 | 38 | -0.400 (0.022) | -0.600 (0.026) |
| rho<0 p<0.05+Q4 drop | pct | GSE131907+GSE148071+GSE205335 | 65 | -0.337 (0.042) | -0.475 (0.023) |
| rho<0 p<0.05+Q4 drop | pct | GSE131907+GSE189357+GSE205335 | 52 | -0.497 (0.000) | -0.692 (0.003) |
| PR320_given Q4 drop | pct | GSE131907+GSE205335 | 23 | — (—) | -0.705 (0.000) |

## Methods (locked)

- CLDN4 only. TACSTD2 is never a gate.
- Spearman of malignant CLDN4 (mean log1p or TISCH mean; and %pos) vs same-unit T/NK fraction.
- Multi-cohort effect = DerSimonian–Laird random-effects on Fisher-z(ρ). I² reported.
- Q4 vs Q1: within-cohort CLDN4 ranks, then Mann–Whitney on T/NK; rank-biserial r. Only if n≥16.
- Honest n = patients (GSE131907 samples; GSE123902 donors). Not cells.
- GSE123902 / GSE189357: marker-malignant = (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0; T/NK = (CD3D|CD3E|CD8A|NKG7|GNLY|KLRD1)>0 and not malignant.
- Reproduce: `python3 methods/scrna_cldn4_combo_enum/analyze.py`

