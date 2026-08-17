# FINDING — CLDN4-only high-end CellChat on the strict-malignant 6-unit combo

ADDITIVE. **CLDN4 only. No dual-high.** Patient is the unit.

The strict-malignant 6-unit mean-vs-T/NK Spearman is **taken as given** and is
**not re-audited**: n=95, ρ=−0.260, p=0.0195, I²=0%
(PR #312 / `methods/strict_malig_tnk`). Units:
GSE207422 (DRMref, 12) + GSE205335 (author mal, 22) + GSE291670 (marker, 6)
+ GSE253013 (marker, 9) + GSE131907 (author mal n_mal≥20, 21) + GSE325414
(author mal, 25).

This folder adds CellChat-style **outgoing malignant → T/NK** on the
**high-end** (within-patient CLDN4 Q4 vs Q1; median split as sensitivity)
and meta-analyzes **patient ΔP** across units. Jin et al. 2021 Hill
probability on CellChatDB v2 protein pairs. CellChat R was not run.

GSE253013’s only public processed matrix is a **9.3 GB RDS** and was
**not downloaded**. That unit contributes locked n=9 to the given Spearman
row and **0** patients to the LR meta.

**Verdict (patient ΔP, not a Spearman re-audit):** outgoing inhibitory
pairs are higher from the CLDN4-high malignant arm (LGALS9–CD45/CD44/P4HB;
HLA-E/F–CD8). NECTIN2–TIGIT is also up (k=3, n=50, ΔP=+0.135). CXCL16–CXCR6
is up, not down (k=4, n=36, ΔP=+0.033), so this is not a clean immune-cold
copy. CD274–PDCD1 is essentially null. Between-patient given-CLDN4 Q4 vs Q1
is led by MDK–NCL (k=4, 21/20, ΔP=+0.345).

## Honest n

Locked combo n=95 (given). CellChat-eligible patients are fewer
because both CLDN4-high and CLDN4-low malignant arms plus T/NK must meet
the cell floor (Q4 vs Q1: n_mal≥40 and ≥10/arm;
median: ≥20 malignant and ≥10 T/NK).

Within-patient Q4 vs Q1 eligible: **n=83**. Median-split eligible: **n=85**.
Units with ≥1 Q4-eligible patient: GSE131907, GSE205335, GSE207422, GSE291670, GSE325414.

| unit | locked n | Q4-eligible | median-eligible | note |
|---|---:|---:|---:|---|
| GSE207422 | 12 | 11 | 12 | DRMref barcodes not public; epithelial proxy on the same 12 locked samples |
| GSE205335 | 22 | 21 | 22 | author malignant / T/NK |
| GSE291670 | 6 | 6 | 6 | marker malignant; no author labels |
| GSE253013 | 9 | 0 | 0 | 9.3 GB RDS not downloaded; LR n=0 |
| GSE131907 | 21 | 21 | 21 | author malignant, locked n_mal≥20 |
| GSE325414 | 25 | 24 | 24 | author EpithelialCells_TumorCells |

GSE207422 DRMref barcodes are not public. CellChat on that unit uses
epithelial lineage as the malignant proxy on the same 12 locked samples.

## Primary: within-patient Q4 vs Q1 ΔP (outgoing Mal → T/NK)

ΔP = P(CLDN4-high mal → same-patient T/NK) − P(CLDN4-low mal → T/NK).
Each patient is one delta. Unit means are random-effects pooled
(DerSimonian–Laird). Tables below keep **k≥2 units and n≥15 patients**.
The full unfiltered meta is in `results/lr_meta_q4q1.tsv`. p-values are
descriptive.

| pair | class | k | n_patients | mean ΔP | p_meta | p_Wilcoxon | I² | units |
|---|---|---:|---:|---:|---|---|---:|---|
| LGALS9_CD45 | inhibitory | 4 | 66 | +0.125 | 4.75e-14 | 2.85e-07 | 0 | GSE131907+GSE205335+GSE207422+GSE325414 |
| LGALS9_CD44 | inhibitory | 4 | 67 | +0.105 | 1.06e-13 | 2.64e-07 | 0 | GSE131907+GSE205335+GSE207422+GSE325414 |
| LGALS9_P4HB | inhibitory | 4 | 65 | +0.047 | 4.11e-10 | 1.81e-07 | +9 | GSE131907+GSE205335+GSE207422+GSE325414 |
| HLA-F_CD8A | inhibitory | 4 | 66 | +0.066 | 2.59e-08 | 5.54e-08 | 0 | GSE131907+GSE205335+GSE207422+GSE325414 |
| COL6A2_CD44 | other | 4 | 34 | +0.120 | 1.22e-07 | 3.22e-07 | 0 | GSE131907+GSE205335+GSE207422+GSE325414 |
| HLA-E_KLRC1 | inhibitory | 3 | 15 | +0.037 | 1.03e-06 | 6.10e-05 | +23 | GSE131907+GSE205335+GSE207422 |
| BAG6_NCR3-PS | other | 2 | 16 | +0.019 | 1.06e-06 | 6.10e-05 | 0 | GSE131907+GSE205335 |
| THBS3_CD47 | other | 3 | 29 | +0.023 | 1.14e-06 | 1.86e-08 | +48 | GSE131907+GSE205335+GSE325414 |
| COL9A2_CD44 | other | 4 | 47 | +0.071 | 3.67e-06 | 5.53e-07 | +9 | GSE131907+GSE205335+GSE207422+GSE325414 |
| ICAM1_ITGAL | other | 4 | 61 | +0.078 | 4.19e-06 | 5.52e-08 | +32 | GSE131907+GSE205335+GSE207422+GSE325414 |
| HLA-E_KLRC2 | inhibitory | 2 | 16 | +0.074 | 6.08e-06 | 6.10e-05 | 0 | GSE131907+GSE205335 |
| HLA-E_CD94:NKG2C | inhibitory | 2 | 16 | +0.095 | 2.01e-05 | 6.10e-05 | +38 | GSE131907+GSE205335 |

Pre-specified barrier / recruit pairs that also pass the reportable gate:

| pair | class | k | n_patients | mean ΔP | p_meta | I² | units |
|---|---|---:|---:|---:|---|---:|---|
| NECTIN2_TIGIT | barrier\|inhibitory | 3 | 50 | +0.135 | 0.00124 | +88 | GSE205335+GSE207422+GSE325414 |
| CXCL16_CXCR6 | recruit | 4 | 36 | +0.033 | 0.00353 | +81 | GSE131907+GSE205335+GSE207422+GSE325414 |
| CDH1_KLRG1 | barrier\|inhibitory | 4 | 34 | +0.043 | 0.000999 | +86 | GSE131907+GSE205335+GSE207422+GSE325414 |
| CD274_PDCD1 | inhibitory | 4 | 15 | +0.000 | 0.720 | 0 | GSE131907+GSE205335+GSE207422+GSE325414 |

Full table: `results/lr_meta_q4q1.tsv`. Reportable slice: `results/lr_meta_q4q1_reportable.tsv` (80 pairs).

## Sensitivity: within-patient median split

| pair | class | k | n_patients | mean ΔP | p_meta | p_Wilcoxon | I² | units |
|---|---|---:|---:|---:|---|---|---:|---|
| LGALS9_CD44 | inhibitory | 4 | 72 | +0.057 | 1.56e-09 | 6.63e-07 | 0 | GSE131907+GSE205335+GSE207422+GSE325414 |
| LGALS9_P4HB | inhibitory | 4 | 70 | +0.028 | 2.32e-09 | 6.07e-07 | 0 | GSE131907+GSE205335+GSE207422+GSE325414 |
| LGALS9_CD45 | inhibitory | 4 | 71 | +0.071 | 1.99e-07 | 7.40e-07 | 0 | GSE131907+GSE205335+GSE207422+GSE325414 |
| BAG6_NCR3-PS | other | 2 | 16 | +0.017 | 4.14e-06 | 6.10e-05 | 0 | GSE131907+GSE205335 |
| COL1A2_CD44 | other | 4 | 19 | +0.032 | 7.65e-06 | 0.00141 | 0 | GSE131907+GSE205335+GSE207422+GSE325414 |
| COL6A2_CD44 | other | 4 | 34 | +0.087 | 1.35e-05 | 1.94e-06 | +9 | GSE131907+GSE205335+GSE207422+GSE325414 |
| HLA-F_CD8A | inhibitory | 4 | 69 | +0.036 | 2.08e-05 | 1.65e-05 | 0 | GSE131907+GSE205335+GSE207422+GSE325414 |
| ICAM1_ITGAL | other | 4 | 63 | +0.055 | 3.08e-05 | 1.98e-08 | +51 | GSE131907+GSE205335+GSE207422+GSE325414 |
| COL9A2_CD44 | other | 4 | 51 | +0.040 | 9.73e-05 | 1.38e-06 | +28 | GSE131907+GSE205335+GSE207422+GSE325414 |
| MDK_NCL | other | 4 | 80 | +0.071 | 1.58e-04 | 1.25e-08 | +48 | GSE131907+GSE205335+GSE207422+GSE325414 |
| HLA-E_CD94:NKG2C | inhibitory | 2 | 16 | +0.078 | 1.93e-04 | 3.05e-05 | +21 | GSE131907+GSE205335 |
| HLA-E_KLRC2 | inhibitory | 2 | 16 | +0.061 | 1.94e-04 | 3.05e-05 | 0 | GSE131907+GSE205335 |


Full table: `results/lr_meta_median.tsv`.

## Extra: between-patient high-end (given CLDN4 Q4 vs Q1)

Quartiles use the **given** malignant CLDN4 scores from the locked tables
(the same vectors as the n=95 ρ=−0.260 row). Per-patient P is scored from
all malignant cells → same-patient T/NK. This is not a re-audit of the
Spearman.

| pair | class | k | n_Q1/n_Q4 | mean ΔP | p_meta | I² | units |
|---|---|---:|---|---:|---|---:|---|
| MDK_NCL | other | 4 | 21/20 | +0.345 | 5.19e-11 | 0 | GSE131907+GSE205335+GSE207422+GSE325414 |
| COL1A1_CD44 | other | 3 | 12/14 | +0.304 | 3.94e-04 | 0 | GSE131907+GSE207422+GSE325414 |
| THBS1_CD47 | other | 3 | 10/8 | +0.118 | 4.77e-04 | 0 | GSE131907+GSE207422+GSE325414 |
| ICAM1_SPN | other | 4 | 16/16 | +0.198 | 0.00115 | 0 | GSE131907+GSE205335+GSE207422+GSE325414 |
| ICAM1_ITGAL | other | 3 | 14/11 | +0.218 | 0.0585 | +32 | GSE205335+GSE207422+GSE325414 |
| COL6A1_CD44 | other | 3 | 8/13 | -0.191 | 0.06 | 0 | GSE131907+GSE205335+GSE325414 |
| LAMB3_CD44 | other | 4 | 19/18 | +0.130 | 0.0721 | 0 | GSE131907+GSE205335+GSE207422+GSE325414 |
| SPP1_CD44 | other | 2 | 11/6 | +0.191 | 0.0766 | 0 | GSE131907+GSE325414 |
| LGALS9_CD44 | inhibitory | 4 | 17/18 | +0.114 | 0.0977 | 0 | GSE131907+GSE205335+GSE207422+GSE325414 |
| LGALS9_P4HB | inhibitory | 4 | 16/18 | +0.069 | 0.152 | 0 | GSE131907+GSE205335+GSE207422+GSE325414 |


Full table: `results/lr_meta_between_q4q1.tsv`.

## Extra figures

- `figures/fig_honest_n.png` — locked vs CellChat-eligible n
- `figures/fig_extra_ligand_table.png` — top within-patient Q4 ΔP pairs
- `figures/fig_extra_preclass.png` — barrier / inhibitory / recruit ΔP
- `figures/fig_forest_q4q1.png` — meta forest
- `figures/fig_patient_delta_top.png` — patient ΔP by unit for LGALS9–CD45
- `figures/fig_forest_between_q4q1.png` — between-patient high-end extra
- `figures/fig_forest_median.png` — median-split sensitivity

## What is not claimed

- The n=95 ρ=−0.260 row is given. It is not recomputed here.
- TACSTD2 is not used to define high/low. This is not dual-high.
- GSE253013 LR is empty because the 9.3 GB RDS was not downloaded.
- Cell-pooled permutations (sibling CellChat PRs) are not the test.
- CellChat R visualizations were not generated.

## Reproduce

```bash
python3 methods/sixunit95_cldn4_hiend/scripts/download.py
python3 methods/sixunit95_cldn4_hiend/scripts/analyze.py
```
