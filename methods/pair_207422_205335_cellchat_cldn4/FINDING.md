# FINDING — pairwise GSE207422 + GSE205335 malignant CLDN4 vs T/NK + outgoing ligands

ADDITIVE. **CLDN4 only.** No dual-high. **Include 207422.** Both cohorts are
ICI-adjacent (GSE207422 neoadjuvant PD-1; GSE205335 palliative ICI biopsy/effusion).
Patient is the unit. Prior TACSTD2 A3 and the multi-cohort Q4 meta are given and
are not re-ranked. p-values are descriptive.

## Honest n

| cohort | malignant label | eligible patients | Q4 vs Q1 tails | notes |
|---|---|---:|---|---|
| GSE207422 | DRMref `Malignant cells` (not Hu CopyKAT) | **12** post-tx | **3 vs 3** (n_compared=6) | P06 has 15 malignant cells; kept on the locked mean table, dropped from LR (need ≥20). Pre-tx P01/P05/P08 are not in DRMref. |
| GSE205335 | author `Malignant cells` | **22** | **6 vs 6** (n_compared=12) | Four GEO patients with (near-)zero malignant cells already out. Q4 is SCLC-heavy. |
| **combo** | within-cohort quartiles, then pool tails | **34** | **9 vs 9** (n_compared=18) | Quartiles are **not** cut on the stacked 34. |

Combo Spearman is DerSimonian–Laird random-effects on Fisher-z(ρ), not a stacked
34-patient correlation (CLDN4 scales differ). Q4 vs Q1 combo is the same estimator
on rank-biserial *r* with n = n_Q1+n_Q4 per cohort.

## Combo rho — malignant CLDN4 vs same-patient T/NK

| set | score | n | Spearman ρ (p, I²) | Q4 vs Q1 r (p; n_Q1/n_Q4) |
|---|---|---:|---|---|
| GSE207422 | mean log1p(CP10k) | 12 | -0.091 (0.779) | -0.333 (0.7; 3/3) |
| GSE205335 | mean log1p(CP10k) | 22 | -0.200 (0.371) | -0.556 (0.132; 6/6) |
| **combo 207422+205335** | **mean (aligned)** | **34** | **-0.166 (0.376, I²=0%)** | **-0.505 (0.0539; n_compared=18)** |
| GSE205335 | %pos (stronger single) | 22 | -0.435 (0.0429) | -0.778 (0.026; 6/6) |
| GSE207422 | %pos (matrix+DRMref, not on locked A3 table) | 12 | -0.678 (0.0153) | -1.000 (0.1; 3/3) |
| **combo 207422+205335** | **%pos (primary ρ)** | **34** | **-0.524 (0.00207, I²=0%)** | -0.999 (0.184; n_compared=18; I²_r high) |

**Combo ρ (primary, %pos):** Fisher-z RE on the two %pos Spearmans is **ρ=-0.524 (p=0.00207, I²=0%, n=34)**. GSE207422 %pos is from the public UMI on DRMref malignant cells (not on the
locked A3 TACSTD2 table). Mean combo is weaker (ρ=−0.166, p=0.376, I²=0%).

Q4 vs Q1 **%pos** combo is **not** the tail claim: GSE207422 3 vs 3 has |r|=1
(P06 enters %pos Q4 with only 15 malignant cells) and the Fisher-z *r* pool
blows up (r≈−1, I²=99%). Mean-quartile tails are the
stable Q4 vs Q1 row (9 vs 9, all LR-eligible; r=−0.505, p=0.0539, I²=0%).
Ligand tables below use **mean** quartiles for that reason.

## CellChat-style outgoing CLDN4-high (Q4) → T/NK

Jin et al. 2021 Hill probability on CellChatDB v2 protein pairs. Outgoing =
same-patient malignant → T/NK. Test = Mann–Whitney on per-patient *P*
(detected ≥3 Q1 and ≥3 Q4). Quartiles are within-cohort
on **mean** CLDN4 (LR-eligible 9 vs 9; %pos Q4 would drop P06). CellChat R was not run.

Detect-gated outgoing rows: **74**. p<0.05: **8**. GSE207422 scored on DRMref compartments from the public UMI (1537 CellChat pairs). LIANA-style combo uses mean-of-means on the same CellChat complexes in both cohorts.

| pair | class | n_Q1/n_Q4 | cohorts | median P Q1 | median P Q4 | Δ | r | p |
|---|---|---|---:|---:|---:|---:|---:|---|
| MDK–NCL ** | other | 9/9 | 2 | 0.438 | 0.757 | +0.320 | +0.975 | 5.74e-04 |
| ICAM1–SPN ** | other | 8/7 | 2 | 0.067 | 0.299 | +0.232 | +0.857 | 0.00373 |
| ICAM1–ITGAL ** | other | 8/6 | 2 | 0.055 | 0.393 | +0.338 | +0.875 | 0.00466 |
| MDK–ITGA4_ITGB1 ** | other | 9/8 | 2 | 0.238 | 0.539 | +0.301 | +0.778 | 0.00551 |
| ICAM1–ITGAL_ITGB2 ** | other | 8/6 | 2 | 0.083 | 0.427 | +0.344 | +0.833 | 0.00799 |
| MICA–NKG2D_HCST ** | recruit | 4/4 | 2 | 0.010 | 0.083 | +0.072 | +1.000 | 0.0286 |
| LAMB3–CD44 ** | other | 9/7 | 2 | 0.154 | 0.325 | +0.171 | +0.651 | 0.0311 |
| GDF15–TGFBR2 ** | other | 4/5 | 2 | 0.017 | 0.187 | +0.170 | +0.900 | 0.0317 |
| NECTIN2–TIGIT | barrier|inhibitory | 9/7 | 2 | 0.144 | 0.367 | +0.223 | +0.556 | 0.0712 |
| THBS1–CD47 | other | 3/5 | 2 | 0.003 | 0.020 | +0.017 | +0.867 | 0.0714 |
| LAMC2–CD44 | other | 7/6 | 2 | 0.146 | 0.276 | +0.130 | +0.619 | 0.0734 |
| LAMA5–CD44 | other | 8/8 | 2 | 0.071 | 0.179 | +0.108 | +0.531 | 0.083 |
| CD55–ADGRE5 | other | 9/8 | 2 | 0.158 | 0.301 | +0.143 | +0.417 | 0.167 |
| LAMC1–CD44 | other | 6/8 | 2 | 0.055 | 0.130 | +0.075 | +0.458 | 0.181 |
| LAMC2–ITGA1_ITGB1 | other | 5/4 | 2 | 0.030 | 0.140 | +0.110 | +0.600 | 0.19 |

Stars = p<0.05. The rest of the table is the next pairs by p; they are not claimed.

## LIANA-style outgoing CLDN4-high (Q4) → T/NK

CellPhoneDB-style score (Efremova 2020 / Garcia-Alonso 2022): partner =
min subunit mean on log1p(CP10k); pair = mean of the two partner means.
GSE207422 uses that rule on DRMref compartments. GSE205335 reuses the given
CellChat truncated-mean complexes (geom-mean, 10% trim) as the partner means
— same patients, not a second matrix pass. Detect gate and MWU as above.

Detect-gated outgoing rows: **74**. p<0.05: **8**.

| pair | class | n_Q1/n_Q4 | cohorts | median S Q1 | median S Q4 | Δ | r | p |
|---|---|---|---:|---:|---:|---:|---:|---|
| MDK–NCL ** | other | 9/9 | 2 | 0.704 | 1.354 | +0.650 | +1.000 | 4.12e-04 |
| MDK–ITGA4_ITGB1 ** | other | 9/8 | 2 | 0.466 | 0.987 | +0.522 | +0.917 | 5.76e-04 |
| CDH1–KLRG1 ** | barrier|inhibitory | 6/6 | 2 | 0.197 | 0.313 | +0.116 | +0.889 | 0.00866 |
| NECTIN2–CD226 ** | barrier|inhibitory | 7/5 | 2 | 0.125 | 0.241 | +0.116 | +0.829 | 0.0177 |
| ICAM1–SPN ** | other | 8/7 | 2 | 0.292 | 0.478 | +0.185 | +0.679 | 0.0289 |
| CDH1–ITGAE_ITGB7 ** | barrier|inhibitory | 7/7 | 2 | 0.229 | 0.373 | +0.144 | +0.673 | 0.0379 |
| ICAM1–ITGAL ** | other | 8/6 | 2 | 0.254 | 0.602 | +0.348 | +0.667 | 0.0426 |
| ICAM1–ITGAL_ITGB2 ** | other | 8/6 | 2 | 0.392 | 0.671 | +0.278 | +0.667 | 0.0426 |
| GDF15–TGFBR2 | other | 4/5 | 2 | 0.201 | 0.441 | +0.241 | +0.800 | 0.0635 |
| CLEC2B–KLRB1 | other | 6/3 | 2 | 0.250 | 0.132 | -0.118 | -0.778 | 0.0952 |
| MIF–CD74_CD44 | other | 6/7 | 2 | 1.963 | 1.712 | -0.251 | -0.524 | 0.138 |
| APP–SORL1 | other | 6/7 | 2 | 0.308 | 0.475 | +0.167 | +0.476 | 0.181 |
| COL6A1–CD44 | other | 4/4 | 2 | 0.436 | 0.659 | +0.223 | +0.625 | 0.2 |
| THBS1–CD47 | other | 3/5 | 2 | 0.203 | 0.272 | +0.069 | +0.600 | 0.25 |
| PVR–TIGIT | barrier|inhibitory | 3/6 | 2 | 0.526 | 0.304 | -0.222 | -0.556 | 0.262 |


## What is not claimed

- Dual-high TACSTD2×CLDN4. TACSTD2 is not a gate.
- A stacked 34-patient Spearman (batch/scale mix). Combo ρ is Fisher-z RE.
- Cell-pooled truncated means as the test. Patient is the unit.
- Hu CopyKAT malignant IDs on GSE207422 (not public). DRMref is the given label.
- MPR/RECIST as the split. Quartiles are CLDN4, not response.
- CellChat R visualizations or a LIANA `mt.cellphonedb` permutation run.

## Files

- `results/combo_rho.tsv` — singles + Fisher-z combo
- `results/patients_with_quartiles.tsv` — 34-patient table
- `results/ligand_table_cellchat_outgoing.tsv` — combo CellChat-style LR
- `results/ligand_table_liana_outgoing.tsv` — combo LIANA-style LR
- `figures/fig_extra_ligand_table.png` — extra ligand-table figure
- `figures/fig_combo_rho_forest.png` — extra combo-ρ forest
- `METHODS.md` — quartiles, Hill P, CellPhoneDB mean-of-means

