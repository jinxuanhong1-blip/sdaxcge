# FINDING — CLDN4-only high-end CellChat on GSE123902 + GSE205335

ADDITIVE. **CLDN4 only. No dual-high.** Patient (GSE123902 donor) is the unit.
Do **not** add GSE148071. Do **not** pile the seven-cohort pool.

The pair %pos Spearman is **taken as given** from PR #459 and is **not
re-audited**:

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE123902+GSE205335 | %pos | 35 | −0.522 (0.002, 0%) | −0.802 (0.005; 9 vs 9) |

Singles (given, same PR): GSE123902 n=13 ρ=−0.659; GSE205335 n=22 ρ=−0.435.

This folder adds CellChat-style **outgoing CLDN4-high malignant → same-patient
T/NK**. Jin et al. 2021 Hill probability on CellChatDB v2 protein pairs
(10% truncated mean, Kh=0.5, expr_prop ≥ 0.10). CellChat R was not run.

Primary high-end: **within-patient** malignant CLDN4 Q4 vs Q1. Sensitivity:
median split. Extra: between-patient Q4 vs Q1 using the **given** %pos
scores (same vectors as the n=35 row).

## Honest n

Given combo n=35 is **not** the CellChat n. Both CLDN4-high and CLDN4-low
malignant arms plus T/NK must meet the cell floor (Q4 vs Q1: n_mal≥40
and ≥10/arm; median: ≥20 malignant and
≥10 T/NK; all-mal outgoing: n_mal≥20 and
n_tnk≥20).

Within-patient Q4 vs Q1 eligible: **n=34**. Median-split eligible: **n=34**.
All-malignant outgoing (between-patient extra): **n=35**.
Units with ≥1 Q4-eligible patient: GSE123902, GSE205335.

| unit | locked n | all-mal LR | Q4-eligible | median-eligible | note |
|---|---:|---:|---:|---:|---|
| GSE123902 | 13 | 13 | 13 | 12 | Laughney 2020; marker-malignant (EPCAM|KRT8/18/19>0 & PTPRC==0); donor-level; normals dropped |
| GSE205335 | 22 | 22 | 21 | 22 | author lineage.sub == Malignant cells; lineage.total == T/NK cells |

GSE123902 malignant = marker-malignant (EPCAM|KRT8|KRT18|KRT19 > 0 and
PTPRC == 0). GSE205335 malignant = author `Malignant cells`. TACSTD2 is
never a gate.

## Primary: within-patient Q4 vs Q1 ΔP (outgoing Mal → T/NK)

ΔP = P(CLDN4-high mal → same-patient T/NK) − P(CLDN4-low mal → T/NK).
Each patient is one delta. Unit means are random-effects pooled
(DerSimonian–Laird). p-values are descriptive.

| pair | class | k | n_patients | mean ΔP | p_meta | p_Wilcoxon | I² | units |
|---|---|---:|---:|---:|---|---|---:|---|
| APP_CD74 | other | 2 | 33 | +0.275 | 8.18e-14 | 2.33e-10 | 0 | GSE123902+GSE205335 |
| MDK_NCL | other | 2 | 34 | +0.216 | 1.04e-08 | 2.95e-08 | 0 | GSE123902+GSE205335 |
| LAMB2_CD44 | other | 2 | 25 | +0.122 | 1.10e-06 | 1.13e-06 | 0 | GSE123902+GSE205335 |
| CCL20_CCR6 | other | 2 | 8 | +0.006 | 4.93e-06 | 0.00781 | +43 | GSE123902+GSE205335 |
| THBS1_CD47 | other | 2 | 21 | +0.041 | 5.27e-06 | 1.91e-06 | 0 | GSE123902+GSE205335 |
| THBS3_CD47 | other | 2 | 17 | +0.034 | 1.78e-05 | 1.53e-05 | 0 | GSE123902+GSE205335 |
| COL4A4_CD44 | other | 2 | 16 | +0.063 | 1.79e-05 | 6.10e-05 | 0 | GSE123902+GSE205335 |
| HLA-E_CD94:NKG2C | inhibitory | 1 | 11 | +0.134 | 1.91e-05 | 9.77e-04 | 0 | GSE205335 |
| HLA-E_KLRK1 | inhibitory | 1 | 19 | +0.136 | 2.98e-05 | 0.00141 | 0 | GSE205335 |
| LGALS9_CD45 | inhibitory | 2 | 27 | +0.119 | 8.48e-05 | 0.00138 | 0 | GSE123902+GSE205335 |
| LGALS9_CD44 | inhibitory | 2 | 28 | +0.092 | 1.01e-04 | 0.00266 | 0 | GSE123902+GSE205335 |
| CDH1_KLRG1 | barrier|inhibitory | 2 | 18 | +0.074 | 1.06e-04 | 7.63e-06 | +62 | GSE123902+GSE205335 |


Primary table (`results/lr_table.tsv`) keeps pairs with **n_patients ≥ 8**.
Tiny-n rows (n=2–4) stay in `results/lr_meta_q4q1.tsv` and are not the claim.

## Sensitivity: within-patient median split

| pair | class | k | n_patients | mean ΔP | p_meta | p_Wilcoxon | I² | units |
|---|---|---:|---:|---:|---|---|---:|---|
| LAMB3_CD44 | other | 2 | 27 | +0.280 | 7.34e-15 | 1.49e-08 | 0 | GSE123902+GSE205335 |
| APP_CD74 | other | 2 | 33 | +0.230 | 2.51e-12 | 4.66e-10 | 0 | GSE123902+GSE205335 |
| CD55_ADGRE5 | other | 2 | 31 | +0.118 | 7.89e-09 | 9.31e-10 | +40 | GSE123902+GSE205335 |
| LAMA3_CD44 | other | 2 | 19 | +0.124 | 1.01e-08 | 1.91e-05 | 0 | GSE123902+GSE205335 |
| LAMC2_CD44 | other | 2 | 26 | +0.205 | 1.91e-08 | 1.23e-05 | +12 | GSE123902+GSE205335 |
| LAMB2_CD44 | other | 2 | 26 | +0.102 | 2.13e-07 | 1.64e-06 | 0 | GSE123902+GSE205335 |
| MDK_NCL | other | 2 | 34 | +0.180 | 2.52e-07 | 4.32e-08 | 0 | GSE123902+GSE205335 |
| ICAM1_SPN | other | 2 | 23 | +0.117 | 2.84e-07 | 7.15e-07 | +43 | GSE123902+GSE205335 |
| IGFBP3_TMEM219 | other | 2 | 19 | +0.123 | 6.87e-06 | 9.54e-05 | 0 | GSE123902+GSE205335 |
| THBS1_CD47 | other | 2 | 21 | +0.027 | 2.07e-05 | 4.10e-05 | 0 | GSE123902+GSE205335 |
| CDH1_KLRG1 | barrier|inhibitory | 2 | 19 | +0.064 | 2.43e-05 | 3.81e-06 | +30 | GSE123902+GSE205335 |
| COL4A4_CD44 | other | 2 | 16 | +0.081 | 2.99e-05 | 1.53e-04 | 0 | GSE123902+GSE205335 |


Full table: `results/lr_meta_median.tsv`.

## Extra: between-patient high-end (given %pos Q4 vs Q1)

Quartiles use the **given** malignant CLDN4 %pos scores from the locked
tables (the same vectors as the n=35 ρ=−0.522 row). Per-patient P is
scored from all malignant cells → same-patient T/NK. This is not a
re-audit of the Spearman.

| pair | class | k | n_Q1/n_Q4 | mean ΔP | p_meta | I² | units |
|---|---|---:|---|---:|---|---:|---|
| CD99_CD99 | other | 2 | 10/8 | -0.331 | 0.00207 | 0 | GSE123902+GSE205335 |
| CLEC2D_KLRB1 | other | 2 | 9/4 | -0.037 | 0.0037 | 0 | GSE123902+GSE205335 |
| HLA-E_CD8A | inhibitory | 2 | 9/6 | -0.312 | 0.0114 | 0 | GSE123902+GSE205335 |
| ICAM1_SPN | other | 2 | 9/4 | +0.043 | 0.0369 | 0 | GSE123902+GSE205335 |
| APP_SORL1 | other | 2 | 6/6 | +0.115 | 0.0411 | 0 | GSE123902+GSE205335 |
| HLA-E_CD8B | inhibitory | 2 | 9/5 | -0.277 | 0.0666 | +17 | GSE123902+GSE205335 |
| HLA-E_KLRK1 | inhibitory | 1 | 6/4 | -0.238 | 0.0688 | 0 | GSE205335 |
| HLA-B_CD8A | other | 2 | 9/6 | -0.174 | 0.153 | +19 | GSE123902+GSE205335 |
| LGALS9_CD45 | inhibitory | 2 | 9/4 | -0.151 | 0.219 | 0 | GSE123902+GSE205335 |
| LAMB2_CD44 | other | 2 | 5/4 | +0.048 | 0.222 | +26 | GSE123902+GSE205335 |


Full table: `results/lr_meta_between_q4q1.tsv`.

## Extra figures

- `figures/fig_given_combo_rho.png` — given %pos ρ (not re-audited)
- `figures/fig_honest_n.png` — locked vs CellChat-eligible n
- `figures/fig_n_per_patient.png` — per-unit malignant / T/NK floors
- `figures/fig_extra_ligand_table.png` — top within-patient Q4 ΔP pairs
- `figures/fig_forest_q4q1.png` — meta forest
- `figures/fig_patient_delta_top.png` — patient ΔP by unit for the top pair
- `figures/fig_forest_between_q4q1.png` — between-patient high-end extra
- `figures/fig_forest_median.png` — median-split sensitivity

## What is not claimed

- The n=35 ρ=−0.522 / Q4 r=−0.802 row is given. It is not recomputed here.
- TACSTD2 is not used to define high/low. This is not dual-high.
- GSE148071 is not added. The seven-cohort pool is not piled.
- Cell-pooled permutations are not the test. Patient is the unit.
- CellChat R visualizations were not generated.

## Reproduce

```bash
python3 methods/pair_123902_205335_hiend_cldn4/scripts/download.py
python3 methods/pair_123902_205335_hiend_cldn4/scripts/analyze.py
```

Hill constants: trim=0.1, Kh=0.5, expr_prop=0.1.
