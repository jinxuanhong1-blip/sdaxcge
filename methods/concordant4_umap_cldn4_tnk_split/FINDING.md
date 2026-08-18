# Concordant-4 UMAP split by unit CLDN4 quartile

ADDITIVE visualization only. **CLDN4-only.**

This figure shows the already-reported inverse association
(higher malignant CLDN4 ↔ fewer T/NK) on the cell-level Harmony UMAP.
A single mixed UMAP cannot show it: T/NK and malignant occupy different
clusters. The same embedding is therefore split by the **unit's** CLDN4
quartile.

**This is not a new test.** PR #503 T/NK numbers were not re-audited:
n=65, %pos ρ=-0.531, stacked Q4 vs Q1 n=19/16
r=-0.724. No new Spearman is quoted.

Datasets ONLY: **GSE123902 + GSE131907 + GSE205335 + GSE189357**.
Not GSE148071 / GSE127465 / GSE154826 / GSE207422. No dual-high.

## Honest n

- **n_units Q1 = 19**, **n_units Q4 = 16**
  (PR #503 locked labels: within-cohort malignant CLDN4 %pos rank then qcut).
- **n_cells shown Q1 = 6650**, **n_cells shown Q4 = 5503**
  (atlas cells after QC + cap ≤350/unit, from those units).
- T/NK cells in the density panel: Q1=3185,
  Q4=1115 (descriptive cell counts, not the inferential n).
- Inferential n remains the PR #503 units (n=65; tails 19/16).

| dataset | Q1 | Q2 | Q3 | Q4 |
|---|---:|---:|---:|---:|
| GSE123902 | 4 | 3 | 3 | 3 |
| GSE131907 | 6 | 5 | 5 | 5 |
| GSE205335 | 6 | 5 | 5 | 6 |
| GSE189357 | 3 | 2 | 2 | 2 |

## Quartile rule

Unit-level malignant CLDN4 **%pos**, ranked **within each cohort**, then
`qcut` on average-tie ranks → Q1/Q2/Q3/Q4. Same function as PR #503
`assign_quartiles`. The figure uses the **locked PR #503 unit table**
(`data/tnk_units_pr503.tsv`) so Q1/Q4 membership is the reported 19/16,
not a re-cut on the 350/unit atlas subsample.

Atlas-capped malignant CLDN4 %pos, run through the same within-cohort
rank-then-qcut, matched the locked label for **36/65** units
that had ≥1 malignant cell. Mismatches are expected: the atlas %pos is a
capped subsample. The figure follows the locked labels.

Every cell of a unit inherits that unit's quartile.

## What the panels show

- **a** UMAP of cells from CLDN4-Q1 units only. T+NK teal; malignant by
  cell-level CLDN4 (red scale); other grey. Shared xy limits with b.
- **b** Same for CLDN4-Q4 units. Same colour scale. Q4 is the hotter
  malignant CLDN4 / thinner T/NK cloud.
- **c** T/NK count-density (KDE × n) of Q1 vs Q4 on the same UMAP
  coordinates and the same density levels — the “fewer T/NK” picture.
- **d** Unit-averaged compartment fractions (T, NK, malignant, other) in
  Q1 vs Q4. Points are units. n_units = 19 / 16. Not cell-pooled.
- **e** Schematic only: CLDN4-high units → T/NK down. No new number.

Unit-averaged fractions (descriptive; not a new test):

| compartment | Q1 mean | Q4 mean |
|---|---:|---:|
| T | 0.379 | 0.155 |
| NK | 0.100 | 0.047 |
| malignant | 0.136 | 0.611 |
| other | 0.385 | 0.187 |

## Pipeline (reused, not re-audited)

Same as `methods/concordant4_atlas_umap_annotate` (PR #507):

- QC: n_genes ≥ 200, n_counts ≥ 500, mitochondrial % < 20.
- Cap ≤350 cells / unit.
- Inner-join genes → HVG 2000 → PCA 30 → Harmony `batch=dataset`
  (theta=2.0). Sample is not a second Harmony key.
- Neighbors k=15 on `X_pca_harmony`; one UMAP; Leiden 0.6.
- Annotation: malignant / T / NK / myeloid / B / other
  (author malignant where present, else marker scores).

## What this is not

- Not a re-audit of PR #503 T/NK ρ or IFN/MHC DE.
- Not a new Spearman, MWU, or Q4 vs Q1 test.
- Not a dual-high TACSTD2∩CLDN4 object.
- Not GSE148071 / GSE127465 / GSE154826 / GSE207422.
- Not evidence that CLDN4 *causes* T/NK exclusion.
- Cell counts are not the inferential n.

## Files

- `results/figures/fig_umap_cldn4_tnk_split.png` / `.svg` / `.pdf`
- `results/tables/unit_quartiles.tsv` — locked Q + atlas composition
- `results/tables/composition_unit_avg.tsv` — panel d means
- `data/tnk_units_pr503.tsv` — locked PR #503 unit table

Reproduce:

```bash
python3 methods/concordant4_atlas_umap_annotate/download.py
python3 methods/concordant4_umap_cldn4_tnk_split/analyze.py
```
