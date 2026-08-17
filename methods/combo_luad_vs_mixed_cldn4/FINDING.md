# Combinatorial histology cut: LUAD-only vs mixed-all vs drop-SCLC (CLDN4-only)

ADDITIVE. **CLDN4-only.** Not a bigger merge. No dual-high TACSTD2×CLDN4
score. Three public processed series only: GSE131907 (LUAD), GSE148071
(NSCLC, no processed LUAD/LUSC field), GSE205335 (mixed; some SCLC).
Patient is the unit (GSE131907 T/NK extract is **sample-level**; that is
stated on those rows). Spearman / DerSimonian–Laird / Stouffer and
Mann–Whitney Q4 vs Q1 (rank-biserial *r*). p-values are descriptive.

Existing processed malignant scores only. Matrices were not re-downloaded.
GSE148071 locked TISCH continuous ρ is taken as given and recomputed here
only to confirm the file still matches (eligible malignant n=22,
mean ρ=+0.080; locked eligible-all n=25 mean ρ=+0.135).

## Three-row table (the answer)

Primary estimand: **malignant CLDN4 %pos vs T/NK fraction**.
Q4 vs Q1 is extra (quartile tails inside each cohort, then Fisher-z pooled;
N_compared = n_Q1 + n_Q4, not the full n). **Bigger n is not better.**
The mixed-all row is larger because it adds unlabeled NSCLC and SCLC;
that is a different population, not more of the same LUAD signal.

| cut | k | N | members | Spearman ρ (p, I²) | Q4 vs Q1 r (p, I²; tails) |
|---|---:|---:|---|---|---|
| LUAD-only | 2 | 35 | GSE131907+GSE205335 | -0.412 (0.0183, 0%) | -0.442 (0.0869, 0%; n_Q1=10/n_Q4=9) |
| mixed-all | 3 | 65 | GSE131907+GSE148071+GSE205335 | -0.337 (0.0418, 40%) | -0.533 (0.0429, 55%; n_Q1=18/n_Q4=17) |
| drop-SCLC | 3 | 61 | GSE131907+GSE148071+GSE205335 | -0.253 (0.132, 35%) | -0.301 (0.128, 0%; n_Q1=17/n_Q4=16) |

## Which cut differs

- **LUAD-only** (k=2, N=35): ρ=-0.412 p=0.0183 I²=0%. GSE131907 (LUAD series, n=21, ρ=−0.522) + GSE205335 ADC (n=14, ρ=−0.204).
The LUAD pool is carried by GSE131907. GSE205335 ADC alone is not a hit.
GSE148071 is **out** — the locked extract has no LUAD/LUSC column.
- **mixed-all** (k=3, N=65): ρ=-0.337 p=0.0418 I²=40%. Adds GSE148071 malignant-eligible NSCLC (n=22, ρ=−0.015) and the full
GSE205335 mix (n=22, ρ=−0.435, Q4 r=−0.778). Larger N, **weaker**
Spearman and I²=40%. The extra Q4 row looks stronger only because
GSE205335 SCLC tails sit in Q4; that is not a LUAD effect.
- **drop-SCLC** (k=3, N=61): ρ=-0.253 p=0.132 I²=35%. GSE205335 without SCLC is n=18, ρ=−0.181, Q4 r=−0.200. The mixed
negative **collapses**. SCLC is not free extra n.

Do not write that mixed-all is the preferred estimate because N is larger
or because its Q4 p is 0.043. LUAD-only is the LUAD cut (I²=0%).
mixed-all is the histology-mixed cut. drop-SCLC is the sensitivity
that shows which histology was doing the work.

## Honest n

| item | n | note |
|---|---:|---|
| GSE131907 samples in extract | 58 | full T/NK extract; not the test n |
| GSE131907 author_malig tumor-site | 21 | origins tLung/tL/B/mLN/PE/mBrain, n_malignant≥20; sample-level; LUAD series |
| GSE148071 GEO patients | 42 | Wu 2021; one sample each; not the test n |
| GSE148071 TISCH eligible | 25 | ≥20 scored epi/mal + ≥20 T/NK; locked given-ρ set |
| GSE148071 TISCH eligible malignant | 22 | drops P5/P35/P39 epithelial-like; used in mixed-all and drop-SCLC |
| GSE148071 LUAD-labeled in extract | 0 | no processed LUAD/LUSC field; out of LUAD-only |
| GSE205335 patients in extract | 22 | locked author-malignant table |
| GSE205335 ADC (LUAD) | 14 | LUAD-only member |
| GSE205335 SQ | 3 | in mixed-all and drop-SCLC; not LUAD |
| GSE205335 SCLC | 4 | in mixed-all only; dropped in drop-SCLC |
| GSE205335 NUT | 1 | not SCLC; kept in drop-SCLC; not LUAD |
| LUAD-only N | 35 | GSE131907 + GSE205335 ADC |
| mixed-all N | 65 | three series, all eligible |
| drop-SCLC N | 61 | mixed-all minus GSE205335 SCLC |

## Members behind the three rows (malignant CLDN4 %pos vs T/NK)

| cut | cohort | histology | n | unit | ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) |
|---|---|---|---:|---|---|---|
| LUAD-only | GSE131907 | LUAD | 21 | sample | -0.522 (0.0152) | -0.600 (0.126; 6/5) |
| LUAD-only | GSE205335 | LUAD | 14 | patient | -0.204 (0.483) | -0.125 (0.886; 4/4) |
| mixed-all | GSE131907 | LUAD | 21 | sample | -0.522 (0.0152) | -0.600 (0.126; 6/5) |
| mixed-all | GSE148071 | NSCLC_unlabeled | 22 | patient | -0.015 (0.946) | -0.056 (0.937; 6/6) |
| mixed-all | GSE205335 | LUAD;LUSC;NUT;SCLC | 22 | patient | -0.435 (0.0429) | -0.778 (0.026; 6/6) |
| drop-SCLC | GSE131907 | LUAD | 21 | sample | -0.522 (0.0152) | -0.600 (0.126; 6/5) |
| drop-SCLC | GSE148071 | NSCLC_unlabeled | 22 | patient | -0.015 (0.946) | -0.056 (0.937; 6/6) |
| drop-SCLC | GSE205335 | LUAD;LUSC;NUT | 18 | patient | -0.181 (0.473) | -0.200 (0.69; 5/5) |

## Mean-score sensitivity (not the three-row)

Same cuts, malignant CLDN4 **mean** vs T/NK. GSE148071 locked given
Spearman is this score. Do not swap this in as the answer because it
is closer to zero.

| cut | k | N | members | Spearman ρ (p, I²) | Q4 vs Q1 r (p, I²; tails) |
|---|---:|---:|---|---|---|
| LUAD-only | 2 | 35 | GSE131907+GSE205335 | -0.308 (0.0868, 0%) | -0.521 (0.0374, 0%; n_Q1=10/n_Q4=9) |
| mixed-all | 3 | 65 | GSE131907+GSE148071+GSE205335 | -0.175 (0.218, 14%) | -0.396 (0.0326, 0%; n_Q1=18/n_Q4=17) |
| drop-SCLC | 3 | 61 | GSE131907+GSE148071+GSE205335 | -0.113 (0.472, 23%) | -0.294 (0.138, 0%; n_Q1=17/n_Q4=16) |

## Methods (this slice)

- Predictor: malignant CLDN4 only. TACSTD2 is not a gate.
- GSE131907: author Malignant cells, tumor-site origins
  (tLung / tL/B / mLN / PE / mBrain), n_malignant ≥ 20. Sample-level.
- GSE205335: author malignant table as locked (n=22). ADC = LUAD.
  drop-SCLC drops `SCLC` only; NUT (n=1) and SQ stay in that cut.
- GSE148071: TISCH eligible (≥20 scored epithelial/malignant and ≥20 T/NK)
  **and** `epi_definition==Malignant`. Three epithelial-like eligible
  units (P5/P35/P39) are dropped. Not in LUAD-only.
- Quartiles: `pd.qcut(rank(method='average'), 4)` inside each cohort.
  Two-sided Mann–Whitney U, Q4 vs Q1. r = 2U/(n4 n1) − 1.
- Pool: DerSimonian–Laird on Fisher-z. Q4 pool uses n_compared and
  drops thin tails (n<8 or either tail <3, or |r|≈1).
- Thin flag and p-values are descriptive.

## What was not done

- No dual-high TACSTD2×CLDN4 score.
- No merge of extra series (not GSE207422 / GSE253013 / GSE291670 / …).
- No invented GSE148071 LUAD/LUSC labels.
- No claim that the larger mixed n is the better estimate.
- GSE131907 was not collapsed from sample to patient (extract is sample-level).

Figures: `figures/forest_three_cuts_spearman.png` (extra forest),
`figures/forest_three_cuts_q4q1.png`, `figures/forest_members_pct.png`,
`figures/q4q1_box_*.png`.

Tables: `tables/three_row.tsv`, `tables/honest_n.tsv`,
`tables/members.tsv`, `tables/cuts_all_scores.tsv`.

Reproduce: `python3 methods/combo_luad_vs_mixed_cldn4/analyze.py`

