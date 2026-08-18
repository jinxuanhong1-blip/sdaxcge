# Concordant-4 UMAP split by unit CLDN4 quartile

ADDITIVE visualization only. **CLDN4-only.**

This figure shows the already-reported inverse association
(higher malignant CLDN4 ↔ fewer T/NK) on the cell-level Harmony UMAP.
A single mixed UMAP cannot show it: T/NK and malignant occupy different
clusters. The same embedding is therefore split by the **unit's** CLDN4
quartile.

**This is not a new test.** PR #503 T/NK numbers were not re-audited:
n=65, %pos ρ=−0.531, stacked Q4 vs Q1 n=19/16 r=−0.724.
No new Spearman is quoted.

Datasets ONLY: **GSE123902 + GSE131907 + GSE205335 + GSE189357**.
Not GSE148071 / GSE127465 / GSE154826 / GSE207422. No dual-high.

## Honest n

- **n_units Q1 = 19**, **n_units Q4 = 16**
  (PR #503 locked labels: within-cohort malignant CLDN4 %pos rank then qcut).
- **n_cells shown Q1 / Q4** are written by `analyze.py` after the atlas run
  (QC + cap ≤350/unit, cells from those units).
- Inferential n remains the PR #503 units (n=65; tails 19/16).

## Quartile rule

Unit-level malignant CLDN4 **%pos**, ranked **within each cohort**, then
`qcut` on average-tie ranks → Q1/Q2/Q3/Q4. Same function as PR #503
`assign_quartiles`. The figure uses the **locked PR #503 unit table**
(`data/tnk_units_pr503.tsv`) so Q1/Q4 membership is the reported 19/16,
not a re-cut on the 350/unit atlas subsample.

Every cell of a unit inherits that unit's quartile.

## What the panels show

- **a** UMAP of cells from CLDN4-Q1 units only. T+NK teal; malignant by
  cell-level CLDN4 (red scale); other grey. Shared xy limits with b.
- **b** Same for CLDN4-Q4 units. Same colour scale.
- **c** T/NK count-density of Q1 vs Q4 on the same UMAP coordinates.
- **d** Unit-averaged compartment fractions (T, NK, malignant, other).
- **e** Schematic only: CLDN4-high units → T/NK down. No new number.

## What this is not

- Not a re-audit of PR #503 T/NK ρ or IFN/MHC DE.
- Not a new Spearman, MWU, or Q4 vs Q1 test.
- Not a dual-high TACSTD2∩CLDN4 object.
- Not GSE148071 / GSE127465 / GSE154826 / GSE207422.
- Not evidence that CLDN4 *causes* T/NK exclusion.
- Cell counts are not the inferential n.

Reproduce:

```bash
python3 methods/concordant4_atlas_umap_annotate/download.py
python3 methods/concordant4_umap_cldn4_tnk_split/analyze.py
```
