# FINDING — CLDN4-only high-end CellChat / LIANA on GSE123902 + GSE131907

ADDITIVE. **CLDN4 only. No dual-high.** Patient/donor is the unit
(GSE123902 donor; GSE131907 sample). Do **not** add GSE148071. Do
**not** pile the seven-cohort pool.

The pair %pos Spearman is **taken as given** from PR #459 and is **not
re-audited**:

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE123902+GSE131907 | %pos | 34 | −0.575 (0.001, 0%) | −0.700 (0.012; 10 vs 8) |

Singles (given, same PR): GSE123902 n=13 ρ=−0.659; GSE131907 n=21 ρ=−0.522.

This folder adds **outgoing CLDN4-high malignant → same-unit T/NK**:

- CellChat-style: Jin et al. 2021 Hill probability on CellChatDB v2 protein
  pairs (10% truncated mean, Kh=0.5, expr_prop ≥ 0.10).
- LIANA-style: CellPhoneDB mean-of-means on the **same** complexes
  (Efremova 2020 / Garcia-Alonso 2022). Partner = complex mean; pair =
  0.5 × (L + R).

CellChat R and LIANA `mt.cellphonedb` permutations were not run.

Primary high-end: **within-unit** malignant CLDN4 Q4 vs Q1. Sensitivity:
median split. Extra: between-unit Q4 vs Q1 using the **given** %pos
scores (same vectors as the n=34 row).

## Honest n

Given combo n=34 is **not** the LR n. Both CLDN4-high and CLDN4-low
malignant arms plus T/NK must meet the cell floor (Q4 vs Q1: n_mal≥40
and ≥10/arm; median: ≥20 malignant and
≥10 T/NK; all-mal outgoing: n_mal≥20 and
n_tnk≥20).

Within-unit Q4 vs Q1 eligible: **n=34**. Median-split eligible: **n=31**.
All-malignant outgoing (between-unit extra): **n=34**.
Units with ≥1 Q4-eligible patient: GSE123902, GSE131907.

| unit | locked n | all-mal LR | Q4-eligible | median-eligible | note |
|---|---:|---:|---:|---:|---|
| GSE123902 | 13 | 13 | 13 | 12 | Laughney 2020; marker-malignant (EPCAM|KRT8/18/19>0 & PTPRC==0); donor-level; normals dropped |
| GSE131907 | 21 | 21 | 21 | 19 | author Cell_subtype == Malignant cells; Cell_type in {T lymphocytes, NK cells}; sample-level; n_mal>=20 |

GSE123902 malignant = marker-malignant (EPCAM|KRT8|KRT18|KRT19 > 0 and
PTPRC == 0). GSE131907 malignant = author `Malignant cells` (not tS*;
tLung primaries are almost unlabeled as malignant and are not in the
locked n=21). TACSTD2 is never a gate.

## Primary: CellChat-style within-unit Q4 vs Q1 ΔP (outgoing Mal → T/NK)

ΔP = P(CLDN4-high mal → same-unit T/NK) − P(CLDN4-low mal → T/NK).
Each patient/donor/sample is one delta. Unit means are random-effects
pooled (DerSimonian–Laird). p-values are descriptive.

| pair | class | k | n_patients | mean Δ | p_meta | p_Wilcoxon | I² | units |
|---|---|---:|---:|---:|---|---|---:|---|
| LGALS9_CD44 | inhibitory | 2 | 27 | +0.093 | 8.15e-11 | 6.66e-06 | 0 | GSE123902+GSE131907 |
| CD55_ADGRE5 | other | 2 | 31 | +0.110 | 5.70e-10 | 3.07e-08 | 0 | GSE123902+GSE131907 |
| LAMA5_CD44 | other | 2 | 31 | +0.135 | 3.60e-09 | 1.86e-09 | 0 | GSE123902+GSE131907 |
| ICAM1_SPN | other | 2 | 25 | +0.103 | 8.84e-09 | 5.96e-08 | 0 | GSE123902+GSE131907 |
| LGALS9_CD45 | inhibitory | 2 | 27 | +0.111 | 9.37e-09 | 1.14e-05 | 0 | GSE123902+GSE131907 |
| LAMC1_CD44 | other | 2 | 31 | +0.083 | 2.99e-08 | 1.57e-07 | 0 | GSE123902+GSE131907 |
| CCL20_CCR6 | other | 2 | 8 | +0.004 | 5.39e-08 | 0.00781 | +16 | GSE123902+GSE131907 |
| LAMB3_CD44 | other | 2 | 31 | +0.174 | 1.00e-07 | 2.56e-06 | 0 | GSE123902+GSE131907 |
| HLA-A_CD8A | other | 2 | 30 | +0.042 | 1.08e-07 | 1.68e-06 | 0 | GSE123902+GSE131907 |
| LAMC2_CD44 | other | 2 | 26 | +0.107 | 2.99e-07 | 2.98e-08 | 0 | GSE123902+GSE131907 |
| LAMB2_CD44 | other | 2 | 29 | +0.090 | 3.16e-07 | 9.31e-08 | 0 | GSE123902+GSE131907 |
| LGALS9_P4HB | inhibitory | 2 | 26 | +0.048 | 8.89e-06 | 1.91e-05 | +20 | GSE123902+GSE131907 |


Primary table (`results/lr_table.tsv`) keeps pairs with **n_patients ≥ 8**.
Tiny-n rows stay in `results/lr_meta_q4q1.tsv` and are not the claim.

## LIANA-style within-unit Q4 vs Q1 ΔS (outgoing Mal → T/NK)

Same patients, same detect gate, CellPhoneDB-style mean-of-means
S = 0.5 × (L + R). ΔS = S_high − S_low.

| pair | class | k | n_patients | mean Δ | p_meta | p_Wilcoxon | I² | units |
|---|---|---:|---:|---:|---|---|---:|---|
| ICAM1_ITGAL | other | 2 | 18 | +0.203 | 9.60e-18 | 7.63e-06 | 0 | GSE123902+GSE131907 |
| ICAM1_SPN | other | 2 | 25 | +0.212 | 2.23e-16 | 5.96e-08 | 0 | GSE123902+GSE131907 |
| HLA-A_CD8A | other | 2 | 30 | +0.205 | 5.00e-12 | 2.76e-06 | 0 | GSE123902+GSE131907 |
| CDH1_KLRG1 | barrier|inhibitory | 2 | 12 | +0.270 | 2.55e-11 | 4.88e-04 | 0 | GSE123902+GSE131907 |
| HLA-E_CD94:NKG2A | inhibitory | 1 | 8 | +0.233 | 1.54e-10 | 0.00781 | 0 | GSE131907 |
| CD55_ADGRE5 | other | 2 | 31 | +0.229 | 9.88e-10 | 9.31e-09 | 0 | GSE123902+GSE131907 |
| HLA-A_CD8B | other | 2 | 26 | +0.204 | 9.99e-10 | 2.69e-05 | 0 | GSE123902+GSE131907 |
| THBS3_CD47 | other | 2 | 15 | +0.143 | 3.70e-09 | 1.22e-04 | 0 | GSE123902+GSE131907 |
| LAMB2_CD44 | other | 2 | 29 | +0.253 | 2.33e-08 | 1.86e-08 | 0 | GSE123902+GSE131907 |
| HLA-E_KLRC1 | inhibitory | 1 | 9 | +0.210 | 5.84e-08 | 0.00391 | 0 | GSE131907 |
| COL1A1_CD44 | other | 2 | 29 | +0.247 | 6.86e-07 | 2.43e-05 | 0 | GSE123902+GSE131907 |
| LAMC1_CD44 | other | 2 | 31 | +0.232 | 8.53e-07 | 1.02e-07 | 0 | GSE123902+GSE131907 |


LIANA table (`results/lr_table_liana.tsv`) keeps pairs with **n_patients ≥ 8**.
Full meta: `results/lr_meta_q4q1_liana.tsv`.

## Sensitivity: within-unit median split (CellChat)

| pair | class | k | n_patients | mean Δ | p_meta | p_Wilcoxon | I² | units |
|---|---|---:|---:|---:|---|---|---:|---|
| LAMA5_CD44 | other | 2 | 30 | +0.108 | 1.77e-11 | 9.31e-09 | 0 | GSE123902+GSE131907 |
| CD55_ADGRE5 | other | 2 | 29 | +0.090 | 4.50e-10 | 3.73e-08 | 0 | GSE123902+GSE131907 |
| ICAM1_SPN | other | 2 | 24 | +0.079 | 9.18e-09 | 1.19e-07 | +3 | GSE123902+GSE131907 |
| LAMC1_CD44 | other | 2 | 31 | +0.054 | 8.20e-08 | 1.54e-05 | 0 | GSE123902+GSE131907 |
| LGALS9_CD44 | inhibitory | 2 | 26 | +0.064 | 1.07e-07 | 5.76e-05 | 0 | GSE123902+GSE131907 |
| LGALS9_CD45 | inhibitory | 2 | 26 | +0.081 | 4.50e-06 | 5.76e-05 | 0 | GSE123902+GSE131907 |
| LAMB2_CD44 | other | 2 | 28 | +0.072 | 5.01e-06 | 6.56e-07 | 0 | GSE123902+GSE131907 |
| HLA-A_CD8A | other | 2 | 27 | +0.028 | 8.15e-06 | 2.59e-05 | 0 | GSE123902+GSE131907 |
| LGALS9_P4HB | inhibitory | 2 | 25 | +0.029 | 1.44e-05 | 1.29e-04 | +20 | GSE123902+GSE131907 |
| APP_SORL1 | other | 2 | 19 | +0.048 | 8.34e-05 | 7.63e-06 | 0 | GSE123902+GSE131907 |
| NECTIN2_TIGIT | barrier|inhibitory | 1 | 9 | +0.083 | 2.85e-04 | 0.00391 | 0 | GSE123902 |
| HLA-E_CD94:NKG2A | inhibitory | 1 | 8 | +0.049 | 3.89e-04 | 0.00781 | 0 | GSE131907 |


Full table: `results/lr_meta_median.tsv`.

## Extra: between-unit high-end (given %pos Q4 vs Q1)

Quartiles use the **given** malignant CLDN4 %pos scores from the locked
tables (the same vectors as the n=34 ρ=−0.575 row). Per-unit P (or S)
is scored from all malignant cells → that unit’s T/NK. This is not a
re-audit of the Spearman.

CellChat-style:

| pair | class | k | n_Q1/n_Q4 | mean Δ | p_meta | I² | units |
|---|---|---:|---|---:|---|---:|---|
| SPP1_CD44 | other | 2 | 10/5 | -0.302 | 4.86e-04 | 0 | GSE123902+GSE131907 |
| HLA-DPB1_CD4 | other | 2 | 7/5 | +0.082 | 0.0101 | 0 | GSE123902+GSE131907 |
| HLA-DRB1_CD4 | other | 2 | 7/5 | +0.092 | 0.0144 | 0 | GSE123902+GSE131907 |
| APP_CD74 | other | 2 | 10/8 | +0.145 | 0.0146 | 0 | GSE123902+GSE131907 |
| HLA-F_CD8B | inhibitory | 2 | 7/6 | -0.179 | 0.0203 | +38 | GSE123902+GSE131907 |
| LAMB2_CD44 | other | 2 | 8/7 | +0.031 | 0.0263 | 0 | GSE123902+GSE131907 |
| CLEC2D_KLRB1 | other | 2 | 8/5 | -0.028 | 0.032 | +16 | GSE123902+GSE131907 |
| APP_SORL1 | other | 2 | 6/4 | +0.093 | 0.0343 | +2 | GSE123902+GSE131907 |
| HLA-F_CD8A | inhibitory | 2 | 7/7 | -0.238 | 0.0367 | +54 | GSE123902+GSE131907 |
| HLA-DRA_CD4 | other | 2 | 7/5 | +0.134 | 0.0532 | 0 | GSE123902+GSE131907 |


LIANA-style:

| pair | class | k | n_Q1/n_Q4 | mean Δ | p_meta | I² | units |
|---|---|---:|---|---:|---|---:|---|
| LAMB2_CD44 | other | 2 | 8/7 | -0.158 | 1.72e-05 | +33 | GSE123902+GSE131907 |
| APP_SORL1 | other | 2 | 6/4 | +0.196 | 0.00431 | 0 | GSE123902+GSE131907 |
| LGALS9_CD45 | inhibitory | 2 | 7/6 | -0.255 | 0.00696 | +44 | GSE123902+GSE131907 |
| HLA-DMB_CD4 | other | 2 | 6/4 | +0.091 | 0.00862 | 0 | GSE123902+GSE131907 |
| HLA-F_CD8A | inhibitory | 2 | 7/7 | -0.295 | 0.027 | +24 | GSE123902+GSE131907 |
| LAMC1_CD44 | other | 2 | 8/7 | -0.148 | 0.0333 | +58 | GSE123902+GSE131907 |
| HLA-F_CD8B | inhibitory | 2 | 7/6 | -0.292 | 0.0366 | +65 | GSE123902+GSE131907 |
| HLA-B_CD8B | other | 2 | 8/6 | -0.329 | 0.0555 | 0 | GSE123902+GSE131907 |
| HLA-DMA_CD4 | other | 2 | 7/5 | +0.235 | 0.0616 | 0 | GSE123902+GSE131907 |
| SPP1_CD44 | other | 2 | 10/5 | -0.217 | 0.0679 | 0 | GSE123902+GSE131907 |


Full tables: `results/lr_meta_between_q4q1.tsv`, `results/lr_meta_between_q4q1_liana.tsv`.

## Extra figures

- `figures/fig_given_combo_rho.png` — given %pos ρ (not re-audited)
- `figures/fig_honest_n.png` — locked vs LR-eligible n
- `figures/fig_n_per_patient.png` — per-unit malignant / T/NK floors
- `figures/fig_extra_ligand_table.png` — top within-unit Q4 ΔP pairs (CellChat)
- `figures/fig_liana_outgoing.png` — top within-unit Q4 ΔS pairs (LIANA)
- `figures/fig_forest_q4q1.png` — CellChat meta forest
- `figures/fig_forest_liana_q4q1.png` — LIANA meta forest
- `figures/fig_patient_delta_top.png` — patient ΔP by unit for the top CellChat pair
- `figures/fig_forest_between_q4q1.png` — between-unit high-end extra
- `figures/fig_forest_median.png` — median-split sensitivity

## What is not claimed

- The n=34 ρ=−0.575 / Q4 r=−0.700 row is given. It is not recomputed here.
- TACSTD2 is not used to define high/low. This is not dual-high.
- GSE148071 is not added (it dilutes). The seven-cohort pool is not piled.
- Cell-pooled permutations are not the test. Patient/donor/sample is the unit.
- CellChat R visualizations and LIANA `mt.cellphonedb` were not generated.
- GSE131907 tS1–tS3 tLung epithelium is not in the locked author-malignant n=21.

## Reproduce

```bash
python3 methods/pair_123902_131907_hiend_cldn4/scripts/download.py
python3 methods/pair_123902_131907_hiend_cldn4/scripts/analyze.py
```

Hill constants: trim=0.1, Kh=0.5, expr_prop=0.1.
LIANA-style: mean-of-means on the same CellChat complexes.
