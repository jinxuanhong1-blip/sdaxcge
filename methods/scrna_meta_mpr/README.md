# methods/scrna_meta_mpr

Public-only pool of lung neoadjuvant/ICI scRNA series that have MPR or RECIST **and** malignant/epithelial TACSTD2/CLDN4.

A3 (GSE207422) MPR labels are taken as given. This folder does not re-annotate that series.

## n_patients (primary MPR pool)

**53 patients** (36 NMPR, 17 MPR): GSE207422 *n*=12 + GSE241934 IIT *n*=11 + GSE241934 Real *n*=24 + GSE291670 *n*=6.

GSE205335 RECIST is a separate table (*n*=16 malignant-evaluable; 10 NR vs 6 R).

GSE243013 (*n*=243, public MPR) is immune-only and is not in the malignant table.

## Primary inverse-variance meta (Hedges *g*, NMPR − MPR)

| Gene | *k* | *n* | *g* (FE = RE) | 95% CI | *p* | *I²* |
|---|---|---|---|---|---|---|
| TACSTD2 | 4 | 53 | −0.12 | −0.71 to 0.47 | 0.68 | 0 |
| CLDN4 | 4 | 53 | −0.10 | −0.69 to 0.49 | 0.74 | 0 |

Rank-biserial IV: TACSTD2 *r* = −0.05 (*p* = 0.77); CLDN4 *r* = −0.10 (*p* = 0.58). Same 53 patients.

No cohort Mann–Whitney *p* is < 0.05. GSE291670 3 vs 3 cannot go below *p* = 0.10.

## Combinations (15 subsets)

Leave-one-in (*k*=1), pairs (*k*=2), leave-one-out (*k*=3), full (*k*=4). NMPR>MPR = Hedges *g* > 0. GSE243013 is not pooled.

| Gene | subsets with *g* > 0 | of 15 | *g* > 0 and *p* < 0.05 |
|---|---|---|---|
| TACSTD2 | GSE207422; GSE207422+IIT; GSE207422+Real; drop GSE291670 | 4 | 0 |
| CLDN4 | Real; GSE207422+Real; IIT+Real; Real+GSE291670; drop IIT | 5 | 0 |

Full four-cohort *g* is negative for both genes. Tables: `results/combinations.tsv`, `results/combinations_nmpr_gt_mpr.tsv`. Figure: `figures/combinations_g.png`.

## Per-cohort (malignant/epithelial mean, patient unit)

See `results/cohort_effects_compact.tsv` and the forest plots in `figures/`.

## RECIST (GSE205335 only)

| Gene | n NR vs R | mean NR vs R | Hedges *g* | MWU *p* |
|---|---|---|---|---|
| TACSTD2 | 10 vs 6 | 1.06 vs 1.00 | +0.08 | 0.96 |
| CLDN4 | 10 vs 6 | 1.17 vs 1.36 | −0.35 | 0.56 |

Three responders and one non-responder had zero captured malignant cells and drop from this table.

## Files

- `playbook.md` — inclusion, metric, estimators
- `scripts/run_meta.py` — assemble + Hedges *g* / rank-biserial + forest
- `results/leftover_inventory.tsv` — every leftover 2023–2026 hit
- `figures/forest_mpr_tacstd2.png`, `forest_mpr_cldn4.png`, `forest_recist_gse205335.png`, `combinations_g.png`
