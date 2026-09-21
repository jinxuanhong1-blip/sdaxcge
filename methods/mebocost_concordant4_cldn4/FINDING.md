# FINDING — MEBOCOST concordant-four CLDN4-high vs low metabolite senders

ADDITIVE layer on top of the locked CellChat protein LR result (PR #540).
This does not replace barrier/inhibitory ligand probabilities.
CLDN4 only. No dual-high. Concordant four only
(GSE123902 + GSE131907 + GSE205335 + GSE189357).
Do **not** add GSE148071 / GSE127465 / CD45-only. This is **not** a full-pool.

Engine: **MEBOCOST** `infer_commu` (kaifuchenlab), `met_est='mebocost'`
(enzyme-expression mean). Not a Python reimplementation of the score.
COMPASS flux was not applied. Senders = malignant CLDN4-high vs CLDN4-low;
receiver = T/NK. Score = MEBOCOST `Commu_Score` (sender metabolite × receiver sensor).
Δ = score(high→TNK) − score(low→TNK). Honest n = patient / locked sample.

CINE: no installable metabolite-CCC package by that name. Unrelated hits
(CineMA cardiac MRI, comet infrared `cine`) were not used. SpatialDM and
MultiNicheNet were not run because MEBOCOST imported and executed.

Kynurenine–AHR is **not** in the MEBOCOST human sensor table
(`human_met_sensor_update_Oct21_2025.tsv`). It is not scored. The closest
DB row is kynurenic acid–GPR35, which was not pre-specified.

Primary split is malignant **Q4 vs Q1**. Extra: median and %pos.
Permutation draws inside MEBOCOST: n_shuffle=10, seed=12345.
The paired test uses `Commu_Score`, which does not depend on the shuffle count.
Within-sample FDRs are stored and are not the primary test.

## Honest n

| gate | n | note |
|---|---:|---|
| Inventory units loaded | 65 | one row per locked unit that was read |
| Both compartments (n_mal≥20, n_tnk≥20) | 65 | before the Q4 cell floor |
| Q4 vs Q1 family units | 64 | n_mal≥40 and ≥1 detected pre-specified pair |

One locked unit does not enter Q4: GSE205335 P4001 has 27 malignant cells
(floor is n_mal≥40). Every other inventory unit (64/64) has at least one
detected pre-specified pair on the Q4 split.

GSE123902 / GSE189357: epithelium marker-malignant
(EPCAM\|KRT8\|KRT18\|KRT19 > 0 and PTPRC == 0) and T/NK markers, not malignant.
GSE131907 / GSE205335: author malignant and T/NK labels. Normal-tissue
biopsies in GSE205335 are excluded. TACSTD2 is never a gate.

## Primary — immunosuppressive metabolite family (Q4 vs Q1)

Pre-specified. Expectation, stated before the run: CLDN4-high > CLDN4-low
outgoing metabolite signal toward T/NK. ON-thesis would agree. A miss does
not reopen the locked protein-LR barrier result.

Family: n=64, mean Δ=+0.0155, p_W=5.13e-11, observed=high>low, agrees=yes.

How to read the table:

- Inside one unit the T/NK sensor average is shared by the high and low
  senders, so the sign of Δ is the sign of the sender metabolite difference
  when the sensor is expressed.
- The family mean is carried by prostaglandin E2–PTGER4 and PTGER2
  (product enzymes in this DB include PTGES, PTGES2, PTGES3, CBR1, CBR3).
  Those two pairs are detected in nearly every Q4 unit.
- Adenosine–ADORA2B / ADORA2A are the same direction where the sensor is
  detected. ADORA2A is missing from 12 of 13 GSE123902 matrices,
  so that cohort contributes 0 to ADORA2A.
- D-lactic acid–HCAR1 is detected in a minority of units and the Δ is small.
  MEBOCOST estimates that metabolite from HAGH/HAGHL, not from LDHA.
  Do not read it as Warburg L-lactate.
- L-lactic acid–SLC16A1 is **not scored**. In
  `metabolite_associated_gene_reaction_HMDB_summary.tsv` every L-lactic acid
  reaction is substrate-direction, and MEBOCOST only emits a metabolite when
  a product-direction enzyme is present. SLC16A1 is in the matrices. The
  metabolite is not. This is a database limit, not a zero effect.

| pair | metabolite | expect | n | 123902 | 131907 | 205335 | 189357 | mean Δ | p_W | observed | agrees |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|
| Adenosine_ADORA2A | Adenosine | high>low | 27 | 0 | 7 | 12 | 8 | +0.0045 | 1.49e-08 | high>low | yes |
| Adenosine_ADORA2B | Adenosine | high>low | 34 | 9 | 10 | 7 | 8 | +0.0003 | 1.16e-10 | high>low | yes |
| PGE2_PTGER2 | Prostaglandin E2 | high>low | 64 | 13 | 21 | 21 | 9 | +0.0109 | 3.74e-11 | high>low | yes |
| PGE2_PTGER4 | Prostaglandin E2 | high>low | 63 | 13 | 20 | 21 | 9 | +0.0349 | 6.23e-11 | high>low | yes |
| DLactate_HCAR1 | D-Lactic acid | high>low | 18 | 4 | 3 | 3 | 8 | +0.0002 | 3.28e-04 | high>low | yes |
| LLactate_SLC16A1 | L-Lactic acid | high>low | 0 | 0 | 0 | 0 | 0 | NA | NA | not_scored | not_scored |
| FAMILY_immunosuppressive_metabolite | family | high>low | 64 | 13 | 21 | 21 | 9 | +0.0155 | 5.13e-11 | high>low | yes |

## Extra — median

| pair | metabolite | expect | n | 123902 | 131907 | 205335 | 189357 | mean Δ | p_W | observed | agrees |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|
| Adenosine_ADORA2A | Adenosine | high>low | 26 | 0 | 6 | 12 | 8 | +0.0035 | 5.96e-08 | high>low | yes |
| Adenosine_ADORA2B | Adenosine | high>low | 33 | 9 | 9 | 7 | 8 | +0.0003 | 4.66e-10 | high>low | yes |
| PGE2_PTGER2 | Prostaglandin E2 | high>low | 62 | 12 | 19 | 22 | 9 | +0.0094 | 1.31e-10 | high>low | yes |
| PGE2_PTGER4 | Prostaglandin E2 | high>low | 61 | 12 | 18 | 22 | 9 | +0.0283 | 1.92e-10 | high>low | yes |
| DLactate_HCAR1 | D-Lactic acid | high>low | 16 | 3 | 2 | 3 | 8 | +0.0001 | 3.05e-05 | high>low | yes |
| LLactate_SLC16A1 | L-Lactic acid | high>low | 0 | 0 | 0 | 0 | 0 | NA | NA | not_scored | not_scored |
| FAMILY_immunosuppressive_metabolite | family | high>low | 62 | 12 | 19 | 22 | 9 | +0.0123 | 1.81e-10 | high>low | yes |

## Extra — pctpos

| pair | metabolite | expect | n | 123902 | 131907 | 205335 | 189357 | mean Δ | p_W | observed | agrees |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|
| Adenosine_ADORA2A | Adenosine | high>low | 25 | 0 | 6 | 12 | 7 | +0.0045 | 6.56e-06 | high>low | yes |
| Adenosine_ADORA2B | Adenosine | high>low | 32 | 9 | 9 | 7 | 7 | +0.0004 | 4.66e-10 | high>low | yes |
| PGE2_PTGER2 | Prostaglandin E2 | high>low | 62 | 12 | 19 | 22 | 9 | +0.0168 | 3.90e-11 | high>low | yes |
| PGE2_PTGER4 | Prostaglandin E2 | high>low | 61 | 12 | 18 | 22 | 9 | +0.0443 | 8.60e-11 | high>low | yes |
| DLactate_HCAR1 | D-Lactic acid | high>low | 16 | 3 | 2 | 3 | 8 | +0.0003 | 3.05e-05 | high>low | yes |
| LLactate_SLC16A1 | L-Lactic acid | high>low | 0 | 0 | 0 | 0 | 0 | NA | NA | not_scored | not_scored |
| FAMILY_immunosuppressive_metabolite | family | high>low | 62 | 12 | 19 | 22 | 9 | +0.0210 | 9.07e-11 | high>low | yes |

## What is not claimed

- TACSTD2 is not used to define high/low. This is not dual-high.
- GSE148071, GSE127465, and CD45-only libraries are not added.
- This is not a metabolite discovery screen. `descriptive_all_pairs_q4q1.tsv`
  lists every MEBOCOST metabolite–sensor with a non-zero high or low score
  toward T/NK. Those rows are not a new claim.
- Kynurenine–AHR was not tested; it is absent from this sensor table.
- L-lactic acid was not tested; this database has no product-direction enzymes for it.
- D-lactic acid–HCAR1 is not an LDHA result.
- Compass / scFEA flux constraints were not applied.
- Cell-pooled tests are not reported. Honest n is the patient/sample.
- Visium same-spot correlation is not used, and this is not a spatial exclusion test.

## Versions

python 3.12.3
numpy 2.4.4
pandas 3.0.6
scipy 1.18.1
anndata 0.13.4
scikit-learn 1.9.1
mebocost file /home/ubuntu/.local/lib/python3.12/site-packages/mebocost/mebocost.py
MEBOCOST_ROOT /tmp/src/MEBOCOST
MEBOCOST git 3e0be6f
sensor human_met_sensor_update_Oct21_2025.tsv

