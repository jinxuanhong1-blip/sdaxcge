# FINDING — treatment cut: naive vs ICI-adjacent malignant CLDN4 vs T/NK

ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. **Not a bigger merge.**
Patient is the unit. Locked tables are taken as given (PR #279 / #290 / #429).
p-values are descriptive.

The answer is the **cut table**. The naive+ICI pile GSE131907+GSE205335
(N=43, ρ=−0.479, p=0.00152) is **not** the headline.

## Cut table — malignant CLDN4 vs same-patient T/NK

Primary score is **%pos** where the locked table has it. GSE207422 is
**mean only** (no CLDN4 %pos on the A3 given table). GSE179994 is **n=0**.

| cut | cohorts | score | k | n | ρ [95% CI] | p | I² | headline |
|---|---|---|---:|---:|---|---:|---:|---|
| treatment-naive core | GSE131907 | %pos | 1 | **21** | **−0.522 [−0.778, −0.117]** | **0.0152** | — | yes |
| treatment-naive inclusive | GSE131907+GSE148071 | %pos | 2 | **46** | **−0.242 [−0.708, +0.370]** | **0.446** | **76%** | yes (this cut differs) |
| ICI-adjacent core | GSE205335 | %pos | 1 | **22** | **−0.435 [−0.724, −0.017]** | **0.0429** | — | yes |
| ICI-adjacent inclusive | GSE205335+GSE207422 | mean | 2 | **34** | −0.166 [−0.491, +0.200] | 0.376 | 0% | yes (mean-matched ±) |
| ICI-adjacent ± | GSE179994 | — | 0 | **0** | — | — | — | unusable |
| naive+ICI pile | GSE131907+GSE205335 | %pos | 2 | 43 | −0.479 [−0.688, −0.196] | 0.00152 | 0% | **NO — do not headline** |

Machine table: `results/cut_table.tsv`. Extra forest: `figures/forest_cuts_CLDN4_tnk.png`.
Member forest: `figures/forest_members_CLDN4_tnk.png`.

## Honest n

| item | n | note |
|---|---:|---|
| GSE131907 author malignant, ≥20 mal and ≥20 T/NK | **21** | Kim 2020 treatment-naive; tLung tS* is out of this extract |
| GSE148071 eligible, ≥25 epi and ≥25 T/NK | **25 / 42** | Wu 2021 diagnostic biopsies; no ICI labels; epithelium is putative |
| GSE205335 author malignant | **22** | Hu 2022 neoadjuvant ICI; mixed ADC/SQ/SCLC/NUT |
| GSE205335 NSCLC (ADC+SQ) sensitivity | 17 | ρ=−0.206, p=0.428 — not the ICI cut |
| GSE207422 author DRMref | **12** | A3 given; CLDN4 mean only |
| GSE179994 | **0** | T-cell-only processed matrix (PR #416 / #420) |
| Naive inclusive | **46** | 21+25; not a cell stack |
| ICI core | **22** | GSE207422 cannot enter the %pos cut |
| ICI inclusive mean | **34** | 22+12 |
| Naive+ICI pile (not the answer) | 43 | 21+22 |

Do not write n=42 for GSE148071. Do not write n=36 for GSE179994 T-cell metadata.
Do not cite barcode counts as n.

## What the cut does

The two **cores** are the same sign: treatment-naive GSE131907 and
ICI-adjacent GSE205335 are both CLDN4-negative vs T/NK. They do not
differ from each other (Wald on Fisher-z z=−0.35, p=0.73).

The **± members** are the cut:

- Adding GSE148071 to naive makes the combo **null and heterogeneous**
  (I²=76%). GSE148071 itself is ρ=+0.069, p=0.74, n=25.
- Adding GSE207422 to ICI can only be done on **mean**, and that combo
  is null (ρ=−0.166, p=0.376, n=34). GSE207422 mean is ρ=−0.091, p=0.78,
  n=12. GSE179994 cannot be added.
- Gluing the two cores (GSE131907+GSE205335) recovers the strong
  PR #290 / #320 row. That is a **naive+ICI pile**. It is shown in grey
  on the forest so it is not quoted as the treatment-cut answer.

Q4 vs Q1 extra (within-cohort, not the combo): GSE131907 6 vs 5,
r=−0.600, p=0.126; GSE148071 7 vs 6, r=0, p=1; GSE205335 6 vs 6,
r=−0.778, p=0.026.

GSE205335 ADC+SQ only (n=17) is a sensitivity, not a replacement ICI
cut. The ICI-core %pos row includes SCLC/NUT.

## CellChat-style — only on the cut that differs

The cut that differs is **treatment-naive inclusive**
(GSE131907+GSE148071): null / I²=76% vs ICI-core negative.

CellChat R was not run. No new matrix merge. Outgoing Mal → T/NK Hill
probabilities are **reused** from the two naive members (PR #347, PR #348;
Jin et al. 2021, CellChatDB v2, cell-pooled). Table:
`results/cellchat_reuse.tsv`.

Shared higher-from-CLDN4-high pairs: **HLA-E–CD8A**, **LGALS9–CD44**,
**F11R/JAM1–ITGAL/ITGB2**. GSE148071 also has LGALS9–PTPRC, NECTIN2–TIGIT,
MDK; CXCL16–CXCR6 is slightly *higher* in high, not lower; CD274–PDCD1
is not detected.

So the naive inclusive T/NK ρ going to null is **not** because those
barrier / inhibitory pairs disappear on GSE148071. Ligand probabilities
are cell-pooled and are not a 46-patient mixed model. CellChat was **not**
run on the naive+ICI pile and **not** run on GSE205335 / GSE207422 for
this cut.

## What is not claimed

- No dual-high gate.
- No stacked 10x+Singleron Spearman as the combo ρ.
- No GSE179994 malignant CLDN4 (n=0).
- No re-download of GEO UMI for this cut.
- No claim that ICI-adjacent GSE205335 is LUAD-only (NSCLC sensitivity is null).
- The grey GSE131907+GSE205335 row is documentation, not the answer.

## Reproduce

```bash
python3 methods/combo_naive_vs_ici_cldn4/analyze.py
```
