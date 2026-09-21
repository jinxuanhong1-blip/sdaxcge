# FINDING — max CellPhoneDB / Connectome Δ, barrier ligands, concordant-4

ADDITIVE. **CLDN4 only. No dual-high.** Concordant four only
(GSE123902 + GSE131907 + GSE205335 + GSE189357).
Do **not** add GSE148071 / GSE127465 / GSE154826 / GSE200563 / E-MTAB-13526.
This is an expression ligand–receptor contrast. It is **not** a spatial exclusion test.

The family is F11R, NECTIN2, CDH1, and LGALS9. Magnitudes are the liana 1.10
CellPhoneDB `lr_means` and Connectome `expr_prod` (formula checked against
`liana.method.cellphonedb`). Positive Δ is CLDN4-high minus CLDN4-low.
Honest n is the patient. The selection rule is in the script docstring and in METHODS.md.
It was applied by code to the full grid. It was not edited after the grid was seen.

## Winner

Gate **pos_vs_neg**, expr_prop **0.10**, receiver **T/NK**.

| | |
|---|---|
| n patients | 53 (7+17+20+9) |
| mean CellPhoneDB Δ | +0.2533 |
| mean Connectome Δ | +0.1250 |
| patient fraction Δ>0 | 0.981 |
| Wilcoxon p (CellPhoneDB) | 3.36e-10 |
| BH q across the T/NK grid | 7.34e-10 |
| cohorts CellPhoneDB | GSE123902:+;GSE131907:+;GSE205335:+;GSE189357:+ |
| cohorts Connectome | GSE123902:+;GSE131907:+;GSE205335:+;GSE189357:+ |
| mean Δ on edges expressed in both arms | +0.1515 |
| mean fraction of detected edges with low arm at 0 | 0.403 |
| malignant CLDN4 log1p mean, high vs low | 1.835 vs 0.000 |

Reference row is the published quartile split at expr_prop=0.10: n=55, mean Δ +0.1951, patient fraction >0 0.964. Winner / reference = 1.298.

The largest CellPhoneDB Δ on the grid is pos_vs_neg at expr_prop 0.20: mean Δ +0.3873, n=26 with cohort counts 2/6/14/4. That row is not eligible. The locked rule requires every cohort to contribute at least 3 patients.

The same gate on myeloid was not used for selection. n=25 (1+10+13+1), CellPhoneDB mean Δ +0.2344 (4/4; GSE123902:+;GSE131907:+;GSE205335:+;GSE189357:+), Connectome +0.0834 (3/4; GSE123902:-;GSE131907:+;GSE205335:+;GSE189357:+). Myeloid Connectome is not 4/4, and two cohorts contribute one patient each.

## Ligands at the winner

| ligand | n | Δ CPDB | Δ Connectome | frac>0 | p | 123902 | 131907 | 205335 | 189357 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| F11R | 43 | +0.2668 | +0.1453 | 1.000 | 2.27e-13 | +0.2710 | +0.2518 | +0.2567 | +0.3035 |
| NECTIN2 | 52 | +0.2398 | +0.1504 | 0.981 | 5.15e-10 | +0.0822 | +0.2682 | +0.2781 | +0.2278 |
| CDH1 | 50 | +0.2950 | +0.1217 | 1.000 | 1.78e-15 | +0.2824 | +0.2896 | +0.2799 | +0.3757 |
| LGALS9 | 53 | +0.2217 | +0.0992 | 0.679 | 1.31e-05 | +0.0699 | +0.3505 | +0.1711 | +0.2093 |
| FAMILY | 53 | +0.2533 | +0.1250 | 0.981 | 3.36e-10 | +0.1709 | +0.2918 | +0.2435 | +0.2665 |

LGALS9 is the least consistent ligand: patient fraction Δ>0 is 0.679. Every cohort mean for LGALS9 is still positive. F11R and CDH1 are positive in every patient in whom that ligand is detected.

## Edges at the winner

Δ is the mean of patients in whom the receptor complex passes expr_prop.
`low0` is the fraction of those patient-edges whose low arm is 0 and high arm is >0.
BH is across these edges.

| edge | ligand | n | Δ CPDB | Δ Connectome | frac>0 | p | q | cohorts + | low0 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CDH1-KLRG1 | CDH1 | 32 | +0.3296 | +0.1426 | 1.000 | 4.66e-10 | 7.45e-10 | 4/4 | 0.469 |
| CDH1-integrin | CDH1 | 42 | +0.3130 | +0.1259 | 1.000 | 4.55e-13 | 3.64e-12 | 4/4 | 0.500 |
| LGALS9-CD45 | LGALS9 | 60 | +0.2949 | +0.1445 | 0.683 | 1.29e-06 | 1.49e-06 | 4/4 | 0.400 |
| F11R-LFA1 | F11R | 41 | +0.2777 | +0.1632 | 1.000 | 9.09e-13 | 3.64e-12 | 4/4 | 0.390 |
| NECTIN2-TIGIT | NECTIN2 | 52 | +0.2487 | +0.1541 | 1.000 | 3.50e-10 | 7.01e-10 | 4/4 | 0.404 |
| NECTIN2-CD96 | NECTIN2 | 58 | +0.2461 | +0.1533 | 0.966 | 6.03e-11 | 1.61e-10 | 4/4 | 0.414 |
| LGALS9-CD44 | LGALS9 | 61 | +0.2240 | +0.0990 | 0.689 | 1.30e-06 | 1.49e-06 | 4/4 | 0.393 |
| F11R-F11R | F11R | 7 | +0.1698 | +0.0496 | 1.000 | NA | NA | 2/4 | 0.000 |
| LGALS9-TIM3 | LGALS9 | 29 | +0.1148 | +0.0345 | 0.793 | 2.20e-04 | 2.20e-04 | 4/4 | 0.379 |
| CDH1-CDH1 | CDH1 | 4 | +0.1060 | +0.0416 | 1.000 | NA | NA | 3/4 | 0.000 |

## Grid (T/NK, highest CellPhoneDB Δ first, 12 rows)

Eligible means n≥40, every cohort n≥3, and all four cohort means >0 for both magnitudes.

| gate | expr_prop | n | 123902/131907/205335/189357 | Δ CPDB | Δ Connectome | frac>0 | CPDB cohorts | Conn cohorts | eligible | p | q_grid |
|---|---:|---:|---|---:|---:|---:|---:|---:|---|---:|---:|
| pos_vs_neg | 0.20 | 26 | 2/6/14/4 | +0.3873 | +0.2219 | 1.000 | 4/4 | 4/4 | no | 2.98e-08 | 3.58e-08 |
| decile | 0.20 | 25 | 2/5/14/4 | +0.3708 | +0.2321 | 1.000 | 4/4 | 4/4 | no | 5.96e-08 | 6.81e-08 |
| d15 | 0.20 | 26 | 2/6/14/4 | +0.3580 | +0.2174 | 1.000 | 4/4 | 4/4 | no | 2.98e-08 | 3.58e-08 |
| quintile | 0.20 | 26 | 2/6/14/4 | +0.3369 | +0.2016 | 1.000 | 4/4 | 4/4 | no | 2.98e-08 | 3.58e-08 |
| ventile | 0.20 | 21 | 2/4/11/4 | +0.3140 | +0.1963 | 0.952 | 3/4 | 3/4 | no | 8.39e-05 | 8.39e-05 |
| q4q1 | 0.20 | 27 | 3/6/14/4 | +0.3091 | +0.1842 | 0.963 | 4/4 | 4/4 | no | 2.09e-07 | 2.28e-07 |
| tertile | 0.20 | 27 | 3/6/14/4 | +0.2855 | +0.1656 | 0.963 | 4/4 | 4/4 | no | 2.83e-07 | 2.95e-07 |
| pos_vs_neg | 0.10 | 53 | 7/17/20/9 | +0.2533 | +0.1250 | 0.981 | 4/4 | 4/4 | yes | 3.36e-10 | 7.34e-10 |
| decile | 0.10 | 52 | 6/17/20/9 | +0.2386 | +0.1277 | 0.981 | 4/4 | 4/4 | yes | 3.71e-10 | 7.43e-10 |
| ventile | 0.10 | 45 | 4/15/17/9 | +0.2313 | +0.1168 | 0.978 | 4/4 | 4/4 | yes | 5.00e-12 | 6.00e-11 |
| d15 | 0.10 | 54 | 7/18/20/9 | +0.2132 | +0.1171 | 0.981 | 4/4 | 4/4 | yes | 4.68e-10 | 8.64e-10 |
| median | 0.20 | 27 | 3/6/14/4 | +0.2129 | +0.1223 | 1.000 | 4/4 | 4/4 | no | 1.49e-08 | 2.10e-08 |

Full grid: `results/tables/grid.tsv` (24 T/NK rows, 16 eligible).
Patient family scores: `results/tables/patient_family.tsv`.
Figure: `results/figures/max_effect_grid.png`.

## How the family score is built

Each edge uses the liana complex rule (minimum subunit). `lr_means` averages the ligand and receptor means and is 0 if either side is 0. Connectome magnitude is the product of those means. An edge counts for a patient only when the receptor passes expr_prop on T/NK. The ligand score is the mean of that ligand's counted edges. The family score is the unweighted mean of ligands with a counted edge. A patient is in the test only if at least 3 of the 4 ligands are counted. Missing ligands are not filled with zero.

The both-arms row keeps only edges that are nonzero on both the high and the low malignant arm, so that number is not produced by the low arm falling under expr_prop.

## What is not claimed

- This does not measure spatial exclusion, contact, or muzzling.
- TACSTD2 is not a gate. This is not dual-high.
- GSE148071 and the other non-concordant sets are not in this fit.
- The Wilcoxon p describes the selected row. The grid BH q is the multiplicity-adjusted figure for that search.
- Rank-aggregate ρ is not the maximized quantity. ρ is a rank, not an expression Δ.
- Cell-pooled tests are not the honest n.
- Patients in the winner table: 53.
