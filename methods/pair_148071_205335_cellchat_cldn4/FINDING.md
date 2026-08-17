# FINDING — pairwise GSE148071 + GSE205335: malignant CLDN4 vs T/NK + CellChat

ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. **GSE131907 is not in this pair.**
Patient is the unit. Prior single-cohort CellChat folders (PR #348, #362) and the
Q4 meta (PR #320) are given and are not re-ranked. p-values are descriptive.

## Honest n

| Item | n | Note |
| --- | ---: | --- |
| GSE148071 deposited | **42** | Wu et al. 2021; stage III/IV NSCLC biopsies |
| GSE148071 eligible (≥25 epi **and** ≥25 T/NK) | **25** | marker-argmax epithelium = putative malignant |
| GSE148071 excluded | **17** | mostly epithelium with almost no T/NK |
| GSE205335 locked extract | **22** | author malignant; ≥20 mal and ≥20 T/NK (PR #279/#320) |
| Pair (eligible) | **47** | 25 + 22; not 42 + 26 |
| Q4 vs Q1 tails (within-cohort, stacked) | **13 vs 12** | n_compared=25, not 47 |

Malignant definitions are **not the same**: GSE148071 is marker-argmax epithelium
(no GEO labels); GSE205335 is author `Malignant cells`. CLDN4 %pos is not stacked
across platforms. Combo ρ is Fisher-z / DerSimonian–Laird of the two cohort rhos.
GSE205335 Q4 mixes SCLC with ADC; that mix is kept.

## Malignant CLDN4 vs same-patient T/NK

| cohort | score | n | Spearman ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | median T/NK Q1 | median T/NK Q4 | Δ |
|---|---|---:|---|---|---:|---:|---:|
| GSE148071 | %pos | 25 | +0.069 (0.742) | +0.000 (1; 7/6) | 0.080 | 0.074 | -0.006 |
| GSE205335 | %pos | 22 | -0.435 (0.0429) | -0.778 (0.026; 6/6) | 0.514 | 0.116 | -0.398 |
| GSE148071 | mean | 25 | +0.189 (0.365) | -0.048 (0.945; 7/6) | 0.080 | 0.092 | +0.012 |
| GSE205335 | mean | 22 | -0.200 (0.371) | -0.556 (0.132; 6/6) | 0.331 | 0.166 | -0.165 |

### Combo ρ (this pair only)

| score | k | N | DL ρ (p, I²) | 95% CI | stacked rank ρ (p) | stacked Q4 vs Q1 r (p; n_Q1/n_Q4) |
|---|---:|---:|---|---|---|---|
| **%pos (primary)** | 2 | 47 | **-0.190** (0.474, I²=66%) | [-0.615, +0.321] | -0.155 (0.298) | -0.308 (0.201; 13/12) |
| mean | 2 | 47 | +0.003 (0.987, I²=37%) | [-0.365, +0.371] | -0.024 (0.875) | -0.231 (0.341; 13/12) |

Primary cut is **%pos**. GSE148071 is near-null (ρ=+0.069, p=0.74; Q4 vs Q1 r=0).
GSE205335 is the negative arm (ρ=−0.435, p=0.043; r=−0.778, 6/6). The pair is
**heterogeneous** (I²≈66%) and the DL combo is **not** the GSE131907+GSE205335
row from PR #320 (that pair is a different agent). This pair does not recover a
significant CLDN4–T/NK anti-correlation.

### Quartile tails (CLDN4 %pos, within-cohort)

| tail | patient | n_mal | n_TNK | CLDN4 %pos | T/NK frac |
|---|---|---:|---:|---:|---:|
| Q4 | GSE148071 P8 (NSCLC_unlabeled, NA) | 422 | 66 | 90.3 | 0.049 |
| Q4 | GSE148071 P32 (NSCLC_unlabeled, NA) | 381 | 79 | 90.0 | 0.133 |
| Q4 | GSE148071 P22 (NSCLC_unlabeled, NA) | 153 | 231 | 86.3 | 0.359 |
| Q4 | GSE148071 P13 (NSCLC_unlabeled, NA) | 165 | 90 | 83.6 | 0.096 |
| Q4 | GSE148071 P12 (NSCLC_unlabeled, NA) | 288 | 46 | 81.9 | 0.052 |
| Q4 | GSE148071 P28 (NSCLC_unlabeled, NA) | 693 | 64 | 81.8 | 0.027 |
| Q1 | GSE148071 P11 (NSCLC_unlabeled, NA) | 106 | 46 | 2.8 | 0.080 |
| Q1 | GSE148071 P42 (NSCLC_unlabeled, NA) | 151 | 454 | 9.3 | 0.454 |
| Q1 | GSE148071 P4 (NSCLC_unlabeled, NA) | 112 | 304 | 25.9 | 0.067 |
| Q1 | GSE148071 P14 (NSCLC_unlabeled, NA) | 690 | 31 | 34.9 | 0.017 |
| Q1 | GSE148071 P23 (NSCLC_unlabeled, NA) | 1466 | 32 | 40.8 | 0.013 |
| Q1 | GSE148071 P24 (NSCLC_unlabeled, NA) | 88 | 64 | 47.7 | 0.166 |
| Q1 | GSE148071 P27 (NSCLC_unlabeled, NA) | 34 | 87 | 50.0 | 0.123 |
| Q4 | GSE205335 P1025 (SCLC, PD) | 980 | 128 | 92.9 | 0.105 |
| Q4 | GSE205335 P1115 (SCLC, PR) | 3769 | 271 | 88.9 | 0.062 |
| Q4 | GSE205335 P1089 (ADC, PD) | 1063 | 272 | 87.8 | 0.169 |
| Q4 | GSE205335 P1084 (ADC, NE) | 196 | 2350 | 86.7 | 0.665 |
| Q4 | GSE205335 P1037 (SQ, PR) | 3671 | 611 | 86.1 | 0.123 |
| Q4 | GSE205335 P1016 (SCLC, PR) | 5336 | 763 | 82.6 | 0.110 |
| Q1 | GSE205335 P1063 (ADC, NE) | 131 | 3144 | 35.1 | 0.732 |
| Q1 | GSE205335 P1119 (ADC, PD) | 1222 | 621 | 36.1 | 0.260 |
| Q1 | GSE205335 P1090 (SQ, PR) | 566 | 247 | 38.2 | 0.253 |
| Q1 | GSE205335 P1015 (ADC, PD) | 291 | 473 | 41.6 | 0.396 |
| Q1 | GSE205335 P1062 (ADC, PD) | 188 | 1934 | 45.2 | 0.633 |
| Q1 | GSE205335 P4001 (ADC, SD) | 27 | 1801 | 48.1 | 0.780 |

## CellChat-style ligands

Not scored in this write-up (matrix step skipped or failed).
Re-run without `--skip-cellchat` to fill `results/ligand_table.tsv`.

## What is not supported

- Treating this pair as GSE131907+GSE205335. That combo is a different folder.
- Dual-high TACSTD2×CLDN4. Groups are CLDN4 only.
- n = 42 + 26 as the communication n. Eligible n is 25 + 22; Q4 vs Q1 is the tails.
- Stacking raw CLDN4 %pos across Singleron vs 10x / author vs marker-argmax.
- Patient-level ICI response or MPR (GSE148071 has none; GSE205335 RECIST is not MPR).
- Cell-pooled high/low permutation as the pair test.

## Files

- `results/combo_rho.tsv` — DL combo ρ + stacked Q4 vs Q1
- `results/q4q1_tnk_singles.tsv` — per-cohort Spearman and Q4 vs Q1
- `results/ligand_table.tsv` — outgoing CLDN4-high → T/NK
- `results/patients_with_quartiles.tsv` — 47-patient table
- `figures/scatter_combo_cldn4_tnk.png` — extra scatter
- `figures/q4q1_tnk_pct.png` — extra Q4 vs Q1 box
- `figures/fig_extra_ligand_table.png` — extra ligand-table figure
- `figures/fig_n_per_patient.png` — extra honest-n bars
- `METHODS.md` — pair rule, quartiles, Hill probability

