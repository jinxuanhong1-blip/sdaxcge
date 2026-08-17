# Finding — triple-merge CLDN4-only AUCell (GSE131907 + GSE148071 + GSE205335)

ADDITIVE. **CLDN4 only.** No TACSTD2 gate and no dual-high score.
Patient is the unit. A10 ELF3–CLDN4 is **given** and is not re-proved.
This slice is the **triple**. The 131907+205335-only SCENIC agent is not
redone; that pair appears only as a sensitivity row.
pySCENIC / cisTarget were not run. AUCell is Aibar recovery on a fixed
prior + background gene universe. p-values are descriptive.

## Verdict

Eligible malignant patients: **n=90** (GSE131907 31, GSE148071 37, GSE205335 22). Stratified Q4 vs Q1 (CLDN4 mean, quartiles **inside each dataset**) is **23 vs 24**. Do not write n=112 (GEO sum) or cell n.

Primary tests use **within-dataset ranks** of malignant CLDN4 mean vs
patient-mean AUCell, then a stratified Q4 vs Q1 on those same ranks.
Positive r / ρ = CLDN4-high patients have **higher** regulon activity.

| Regulon | n | Spearman ρ (p) | stratified Q4 vs Q1 r (p; n_Q4/n_Q1) |
|---|---:|---|---|
| IFN (STAT1/IRF + ISGs) | 90 | -0.059 (0.582) | -0.178 (0.302; 23 vs 24) |
| MHC / antigen presentation | 90 | -0.211 (0.0461) | -0.438 (0.0103; 23 vs 24) |
| TJ / barrier (CLDN4 held out) | 90 | +0.521 (1.37e-07) | +0.728 (1.98e-05; 23 vs 24) |
| ELF3 CollecTRI (A10 given) | 90 | -0.087 (0.414) | -0.170 (0.322; 23 vs 24) |
| housekeeping control | 90 | -0.099 (0.352) | -0.047 (0.79; 23 vs 24) |

**What holds.** TJ / barrier AUCell (CLDN4 held out) is higher in CLDN4-high
patients in the triple and in each dataset (triple ρ=+0.521, p=1.37e-07,
n=90; stratified r=+0.728, 23 vs 24). MHC / antigen-presentation AUCell
is **lower** in CLDN4-high patients in the triple (ρ=-0.211, p=0.0461;
stratified r=-0.438, p=0.0103). That MHC inverse is **GSE148071-driven**
(n=37, ρ=-0.432, p=0.00759). GSE131907 and GSE205335 are each null for MHC.

**What does not hold.** IFN AUCell does not track CLDN4 (ρ=-0.059, p=0.582).
The 131907+205335 pair (n=53) is also null for IFN and MHC — adding
GSE148071 is the increment, not a redo of that pair.

**A10 given.** ELF3 CollecTRI (CLDN4 was never in that prior) does not
track CLDN4 here (ρ=-0.087, p=0.414). This is supporting, not a new
ELF3–CLDN4 claim.

**Control.** Housekeeping AUCell is null (ρ=-0.099, p=0.352).

Do not upgrade a cell-level pattern to this patient n. Do not read
GSE205335 RECIST as MPR. GSE131907 has no ICI labels.

## Honest n

| item | n | note |
|---|---:|---|
| GEO patients GSE131907 (Kim 2020) | 44 | 58 samples; not the test n |
| GEO patients GSE148071 (Wu 2021) | 42 | one sample each; not the test n |
| GEO patients GSE205335 (Hu 2022/2024) | 26 | 33 samples; not the test n |
| GEO patients summed | 112 | do not write this as the test n |
| GSE131907 patients with any scored malignant cell | 32 | 1 patient has <20 malignant |
| GSE131907 eligible (≥20 malignant) | 31 | test n for this dataset |
| GSE131907 malignant cells scored | 31136 | cell n; not the test n |
| GSE131907 Q1/Q4 tails (CLDN4 mean) | 16 | Q1=8 Q4=8 |
| GSE148071 patients with any scored malignant cell | 42 | 5 patients have <20 malignant |
| GSE148071 eligible (≥20 malignant) | 37 | test n for this dataset |
| GSE148071 malignant cells scored | 48118 | cell n; not the test n |
| GSE148071 Q1/Q4 tails (CLDN4 mean) | 19 | Q1=10 Q4=9 |
| GSE205335 patients with any scored malignant cell | 22 | Normal* tissues already dropped |
| GSE205335 eligible (≥20 malignant) | 22 | test n for this dataset |
| GSE205335 malignant cells scored | 28512 | cell n; not the test n |
| GSE205335 Q1/Q4 tails (CLDN4 mean) | 12 | Q1=6 Q4=6 |
| triple eligible patients | 90 | primary pooled n |
| triple stratified Q1+Q4 (CLDN4 mean) | 47 | Q1=24 Q4=23; n_compared not 112 |
| pair 131907+205335 eligible (sensitivity only) | 53 | not the primary; other agent owns this pair |

Quartiles are assigned **inside each dataset** on eligible patients,
then tails are stacked. n_compared = n_Q1 + n_Q4, not the GEO n and
not the continuous n.

## Per-dataset AUCell (CLDN4 mean, raw scores)

| Dataset | Regulon | n | ρ (p) | Q4 vs Q1 r (p; n_Q4/n_Q1) |
|---|---|---:|---|---|
| GSE131907 | IFN (STAT1/IRF + ISGs) | 31 | -0.127 (0.495) | -0.188 (0.574; 8 vs 8) |
| GSE131907 | MHC / antigen presentation | 31 | -0.021 (0.909) | -0.156 (0.645; 8 vs 8) |
| GSE131907 | TJ / barrier (CLDN4 held out) | 31 | +0.448 (0.0114) | +0.781 (0.00699; 8 vs 8) |
| GSE148071 | IFN (STAT1/IRF + ISGs) | 37 | -0.085 (0.618) | -0.333 (0.236; 9 vs 10) |
| GSE148071 | MHC / antigen presentation | 37 | -0.432 (0.00759) | -0.711 (0.0101; 9 vs 10) |
| GSE148071 | TJ / barrier (CLDN4 held out) | 37 | +0.547 (4.64e-04) | +0.711 (0.0101; 9 vs 10) |
| GSE205335 | IFN (STAT1/IRF + ISGs) | 22 | +0.098 (0.665) | +0.111 (0.818; 6 vs 6) |
| GSE205335 | MHC / antigen presentation | 22 | -0.128 (0.57) | -0.222 (0.589; 6 vs 6) |
| GSE205335 | TJ / barrier (CLDN4 held out) | 22 | +0.606 (0.0028) | +0.944 (0.00433; 6 vs 6) |

## Pair 131907+205335 (sensitivity only — not this agent's claim)

Eligible pair n=53. IFN ρ=-0.032 (p=0.82). MHC ρ=-0.063 (p=0.65).
TJ still holds (ρ=+0.512, p=8.9e-05). The MHC inverse appears only after
GSE148071 is added. This pair is **not** the claim of this PR.

## Methods (this slice)

- Predictor: malignant **CLDN4 only** (mean log1p(CP10k) on GEO UMI;
  TISCH log-normalized values on GSE148071). TACSTD2 is recorded but
  never a gate.
- Malignant labels are **author / TISCH major-lineage**, not a new CNV call.
  GSE131907: epithelial `Malignant cells` / tS1 / tS2 / tS3 in tumor
  sites (tLung, tL/B, mLN, mBrain, PE). GSE148071: TISCH `Malignant`.
  GSE205335: `lineage.sub == Malignant cells`, Normal* tissues dropped.
- Eligible patient: ≥20 malignant cells after that filter. Multiple
  tumor samples from one patient are pooled.
- Regulons: IFN = STAT1/IRF + Hallmark-like ISGs; MHC = NLRC5/CIITA +
  HLA/TAP/immunoproteasome; TJ = ELF3/GRHL + junction genes **minus CLDN4**;
  ELF3 = CollecTRI prior (A10 given; CLDN4 was never in that prior);
  HK = housekeeping control. Gene lists: `resources/regulons.json`.
- AUCell: Aibar linear recovery, threshold 5% of the **fixed universe**
  (regulon members + background). Not full-transcriptome pySCENIC.
  Module score (mean of members) is written beside AUCell.
- Pool: rank CLDN4 and each regulon **within dataset**, then one Spearman.
  Q4 vs Q1: `pd.qcut` on average ranks inside each dataset; Mann–Whitney
  on stacked tails. Rank-biserial r = 2U/(n4 n1) − 1.
- GSE148071 uses TISCH2 `expression.h5` (already log-normalized).
  Ranks are within-dataset, so the scale is not mixed into one AUCell.

## What was not done

- No dual-high TACSTD2×CLDN4 score.
- No re-run of the 131907+205335-only SCENIC agent as the primary.
- No pySCENIC GRNBoost2 / cisTarget motif ranking.
- No ICI response / MPR test (GSE131907 and GSE148071 have none;
  GSE205335 RECIST is not used as MPR).
- No cell-level p-value as the claim (pseudoreplication).

Tables: `tables/regulon_table.tsv` (primary), `tables/patient_scores.tsv`,
`tables/honest_n.tsv`, `tables/contrasts_all.tsv`.

Figures: `figures/fig_q4q1_aucell.png`, `figures/fig_scatter_cldn4_aucell.png`.

Reproduce: `python3 methods/triple_scrna_scenic_cldn4/download.py && 
python3 methods/triple_scrna_scenic_cldn4/analyze.py`

