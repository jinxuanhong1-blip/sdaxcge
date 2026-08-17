# FINDING — QUAD GSE207422+GSE131907+GSE148071+GSE205335 malignant CLDN4 vs T/NK + CellChat

ADDITIVE. **CLDN4 only.** No dual-high. Patient is the unit.
**GSE207422 is IN the merge** (not excluded). GSE207422-only CLDN4 vs T/NK
(flat) and PR #320 (GSE131907+GSE205335 Q4 r=−0.705) are given and are
not re-audited or re-ranked here. p-values are descriptive.

## Honest n

Eligible = ≥20 malignant and ≥20 T/NK cells in the same patient.
Malignant definitions stay cohort-native: author malignant/tS* (GSE131907
tumor-origin, pooled per patient), author malignant (GSE205335), marker
epithelial (GSE148071; GSE207422 post-treatment; CopyKAT IDs are not public).

| cohort | eligible n | malig def | dropped for floor |
|---|---:|---|---:|
| GSE207422 | 12 | marker_epithelial_post | 0 / 12 post-treatment |
| GSE131907 | 31 | author_malig_tS_tumor_origin | 1 / 32 tumor-origin patients |
| GSE148071 | 25 | marker_epithelial | 17 / 42 biopsies |
| GSE205335 | 22 | author_malig | 4 / 26 (near-zero malignant capture) |

QUAD pooled **n=90** patients. Q4 vs Q1 uses
**within-cohort** CLDN4 %pos quartiles, then pools the tails
(n_Q1=24, n_Q4=23, n_compared=47), not the mid quartiles.

## Combo rho (primary = pooled patients)

| analysis | k | N | Spearman ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | Δ T/NK |
|---|---:|---:|---|---|---:|
| QUAD pooled patients | 4 | 90 | -0.113 (0.29) | -0.366 (0.0325; 24/23) | -0.216 |
| GSE207422 only | 1 | 12 | -0.399 (0.199) | -0.778 (0.2; 3/3) | -0.415 |
| GSE131907 only | 1 | 31 | -0.286 (0.119) | -0.344 (0.279; 8/8) | -0.149 |
| GSE148071 only | 1 | 25 | +0.069 (0.742) | +0.000 (1; 7/6) | -0.006 |
| GSE205335 only | 1 | 22 | -0.418 (0.0526) | -0.778 (0.026; 6/6) | -0.398 |
| LOO drop GSE207422 | 3 | 78 | -0.188 (0.0999) | -0.348 (0.0586; 21/20) | -0.213 |
| LOO drop GSE131907 | 3 | 59 | -0.050 (0.704) | -0.358 (0.093; 16/15) | -0.138 |
| LOO drop GSE148071 | 3 | 65 | -0.297 (0.0164) | -0.536 (0.008; 17/17) | -0.273 |
| LOO drop GSE205335 | 3 | 68 | +0.055 (0.658) | -0.242 (0.228; 18/17) | -0.189 |

Fisher-z meta of the four within-cohort Spearmans: k=4 N=90 ρ=−0.239 p=0.0317 I²=9%.

The **pooled-patient Spearman is the honest combo rho** (n=90, ρ=−0.113,
p=0.29). It is weaker than Fisher-z because GSE148071 is near-null
(ρ=+0.069; T/NK already thin, median ~0.08 in both tails) and is not
hidden by weighting. Q4 vs Q1 on the within-cohort tails is
r=−0.366 (p=0.0325; 24 vs 23).

Leave-one-cohort-out is the same pooled-patient Spearman after dropping
that cohort. Drop GSE148071: n=65 ρ=−0.297 p=0.016. Drop GSE205335:
n=68 ρ=+0.055 p=0.66. Drop GSE207422: n=78 ρ=−0.188 p=0.10.
Including GSE207422 is the point of this merge; the given GSE207422-only
flat T/NK claim and PR #320 pair are not re-ranked. Per-cohort rows
are inventory, not a new single-cohort audit.

GSE205335 Q4 vs Q1 r=−0.778 (6/6) matches the given PR #320 single
and is not re-interpreted here.

### Within-cohort tails (CLDN4 %pos, Q4 / Q1)

| cohort | tail | n | median CLDN4 %pos | median T/NK |
|---|---|---:|---:|---:|
| GSE207422 | Q4 | 3 | 94.3 | 0.199 |
| GSE207422 | Q1 | 3 | 62.5 | 0.614 |
| GSE131907 | Q4 | 8 | 93.0 | 0.203 |
| GSE131907 | Q1 | 8 | 53.8 | 0.352 |
| GSE148071 | Q4 | 6 | 85.0 | 0.074 |
| GSE148071 | Q1 | 7 | 34.9 | 0.080 |
| GSE205335 | Q4 | 6 | 87.3 | 0.116 |
| GSE205335 | Q1 | 6 | 39.9 | 0.514 |

Pooled within-cohort tails: Q4 n=23 median T/NK=0.119; Q1 n=24 median T/NK=0.335.

## CellChat-style ligands (patients that pass floors)

Jin et al. 2021 Hill probability on CellChatDB v2 protein pairs
(10% truncated mean, \(K_h=0.5\), `expr_prop ≥ 0.10`). Outgoing =
malignant → **same-patient** T/NK. Incoming = T/NK → malignant.
Test = Mann–Whitney on per-patient *P* among within-cohort Q4 vs Q1
tails (detected in ≥3 Q1 and ≥3 Q4).
CellChat R and LIANA were not run. Cell-pooled stacked means are not the test.

Detect-gated direction×pair rows: 157. p<0.05: 7 (outgoing 6).

| direction | pair | class | n_Q1/n_Q4 | median P Q1 | median P Q4 | Δ | r | p |
|---|---|---|---|---:|---:|---:|---:|---|
| outgoing | MDK–NCL ** | other | 24/23 | 0.523 | 0.723 | +0.200 | +0.554 | 0.00117 |
| outgoing | APP–CD74 ** | other | 23/23 | 0.491 | 0.666 | +0.175 | +0.422 | 0.0147 |
| outgoing | NECTIN2–CD226 ** | barrier|inhibitory | 9/3 | 0.030 | 0.104 | +0.074 | +0.926 | 0.0182 |
| outgoing | ICAM1–SPN ** | other | 19/15 | 0.095 | 0.187 | +0.092 | +0.481 | 0.0183 |
| incoming | GRN–SORT1 ** | other | 7/9 | 0.001 | 0.009 | +0.008 | +0.651 | 0.0311 |
| outgoing | LGALS9–P4HB ** | inhibitory | 14/15 | 0.031 | 0.200 | +0.168 | +0.438 | 0.0471 |
| outgoing | F11R–ITGAL_ITGB2 ** | barrier | 17/10 | 0.155 | 0.427 | +0.272 | +0.471 | 0.0473 |
| outgoing | SPP1–CD44 | other | 14/9 | 0.585 | 0.089 | -0.496 | -0.492 | 0.0547 |
| outgoing | CDH1–KLRG1 | barrier|inhibitory | 12/9 | 0.021 | 0.059 | +0.038 | +0.500 | 0.0597 |
| outgoing | CLEC2B–KLRB1 | other | 17/13 | 0.031 | 0.090 | +0.059 | +0.394 | 0.0719 |
| outgoing | MDK–ITGA4_ITGB1 | other | 20/13 | 0.211 | 0.412 | +0.200 | +0.354 | 0.0937 |
| outgoing | LGALS9–CD44 | inhibitory | 15/17 | 0.101 | 0.329 | +0.228 | +0.349 | 0.0966 |
| outgoing | CD69–KLRB1 | other | 9/3 | 0.124 | 0.006 | -0.118 | -0.704 | 0.1 |
| outgoing | HLA-E–CD94:NKG2C | inhibitory | 3/3 | 0.521 | 0.134 | -0.386 | -1.000 | 0.1 |
| outgoing | HLA-E–KLRC2 | inhibitory | 3/3 | 0.343 | 0.036 | -0.307 | -1.000 | 0.1 |

Stars mark p<0.05. The rest of the table is the next detect-gated
pairs by p; they are not claimed. Full table:
`tables/ligand_table.tsv` / `results/ligand_table.tsv`.

## Files

- `tables/combo_rho.tsv` — QUAD / per-cohort / LOO Spearman + Q4 vs Q1
- `tables/ligand_table.tsv` — CellChat-style differential pairs
- `results/patients_eligible.tsv` — patient table + within/global quartiles
- `figures/scatter_quad_cldn4_tnk.png` — extra scatter
- `figures/q4q1_tnk_within.png` — extra Q4 vs Q1 box
- `figures/loo_rho_forest.png` — extra LOO forest
- `figures/fig_extra_ligand_table.png` — extra ligand-table figure
- `METHODS.md` — floors, labels, Hill probability

