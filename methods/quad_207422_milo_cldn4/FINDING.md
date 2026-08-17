# FINDING — QUAD Milo neighbourhoods vs malignant CLDN4

Additive **CLDN4-only** neighbourhood DA on the public merge of
GSE207422 + GSE131907 + GSE148071 + GSE205335. GSE207422 is **included**
by importing PR #323 (null, honest n=7). That 207422-only Milo is **not**
re-run. No TACSTD2 gate and no dual-high.

The independent unit is the **patient / sample**, not the cell and not the
overlapping neighbourhood. Graphs are **per dataset** (tLung and mBrain
are also split inside GSE131907). Harmony on ~400k cells from four
chemistries / label schemes was not used (15 GB RAM; site and treatment
are confounded with dataset). miloR / edgeR were not used.

## Verdict

| Test | Honest n | Result |
|---|---|---|
| Sample malignant CLDN4 vs T/NK, per-graph | see forest | Stouffer signed z=-0.840, p=4.01e-01, k=4 graphs |
| Same, within-dataset ranks pooled | **78** samples | ρ=-0.097, p=3.98e-01 |
| Nhood abundance vs malignant CLDN4 (SpatialFDR<0.1) | samples below, not cells | stacked table has **11,060** nhoods |

Do not cite cell count as *n*. Neighbourhoods are not independent across
or within a graph. Dataset-native CLDN4 scales differ (log1p-CP10k vs
hunt log2TPM fallback); the pooled test uses **within-dataset ranks**.

## Honest n (sample / patient is the unit)

| Dataset | Graph | Unit | n scored (≥10 malignant) | n in graph | Notes |
|---|---|---|---:|---:|---|
| GSE207422 | post_treatment | sample | **7** | 12 | imported PR #323; 5/12 post samples dropped (<10 malignant-like); not re-run |
| GSE131907 | tLung | sample | **10** | 11 | Milo nhood table pending; sample scores from hunt PR 230 (log2TPM, author labels) |
| GSE148071 | all_biopsies | patient | **39** | 42 | 39/42 with ≥10 malignant-like; marker-reconstructed lineage (no public CNA IDs) |
| GSE205335 | non_normal | patient | **22** | 22 | imported PR #385; 22 patients; 4 normal-only patients dropped; MPR n=0 |

Combined scored n across the four studies (sum of the rows above, not a
single mixed graph): **78**. GSE131907 sites are not collapsed.

## Sample-level malignant CLDN4 vs T/NK

| Graph | n | ρ | p | source |
|---|---:|---:|---:|---|
| GSE207422 post (imported #323) | 7 | -0.071 | 8.79e-01 | imported_PR323 |
| GSE131907 tLung (hunt fallback) | 10 | -0.188 | 6.03e-01 | imported_hunt_PR230 |
| GSE131907 other origins (hunt, not primary n) | 21 | -0.617 | 2.89e-03 | imported_hunt_PR230 |
| GSE148071 biopsies | 39 | 0.004 | 9.79e-01 | imported_PR383_scores_only |
| GSE205335 patients | 22 | -0.263 | 2.38e-01 | imported_PR385 |

Signed Stouffer (weights √n): z=-0.840, p=4.01e-01,
k=4 graphs. Within-dataset rank pool: n=78,
ρ=-0.097, p=3.98e-01.

## Neighbourhood DA vs malignant CLDN4 (graph per dataset)

| Dataset | Graph | nhoods | testable | P<0.05 | min SpatialFDR | SpatialFDR<0.1 | BH<0.1 | table |
|---|---|---:|---:|---:|---:|---:|---:|---|
| GSE207422 | post_treatment | 5347 | 20 | 2 | 0.368 | **0** | 0 | yes_imported |
| GSE131907 | tLung | None | None | None | NA | **None** | None | no |
| GSE148071 | all_biopsies | 6401 | 881 | 81 | 0.000 | **8** | 8 | summary_only |
| GSE205335 | non_normal | 5713 | 360 | 17 | 0.899 | **0** | 0 | yes |

DA model = sample-level Spearman of `prop[sample, nhood] = n_cells(sample in nhood) / n_cells(sample)`
versus that sample’s mean malignant CLDN4. SpatialFDR = miloR `graphSpatialFDR`
k-distance weights. This is **not** edgeR QLF and **not** a joint Harmony graph.

Nhood tables still missing from this merge (counts may come from a published summary): GSE131907 nhoods, GSE148071 nhoods.

## What this does not say

- It does not redo the GSE207422-only Milo (PR #323). Those 20 testable
  nhoods, 0 SpatialFDR<0.1, n=7 are imported and sit in the stacked table.
- It does not treat 90k–200k cells as *n*.
- It does not invent a shared ICI / MPR label. Only GSE207422 has MPR;
  GSE205335 has RECIST; GSE131907 is treatment-naive; GSE148071 is a
  diagnostic-biopsy atlas.
- Neighbourhood “next to” is kNN co-membership, not histology.
- Unrestricted nhood CLDN4 vs T/NK is partly lineage geometry.
- mRNA ≠ protein. UMI CLDN4 is not an IHC H-score.
- No dual-high TACSTD2+CLDN4 gate.

## Files

- `tables/nhoods.tsv` — stacked neighbourhood table (the merge deliverable)
- `tables/da_malignant_cldn4.tsv` — CLDN4 DA columns only
- `tables/sample_scores.tsv`, `tables/honest_n.tsv`, `tables/forest_sample_cldn4_vs_tnk.tsv`
- `tables/summary.json`
- `figures/fig_sample_cldn4_vs_tnk_by_dataset.png`
- `figures/fig_sample_ranks_cldn4_vs_tnk.png`
- `figures/fig_forest_sample_cldn4_vs_tnk.png`
- `figures/fig_nhood_volcano_faceted.png`
- `figures/fig_spatialfdr_counts.png`
- `figures/fig_interface_cldn4_vs_tnk.png`
- `figures/fig_honest_n.png`

```bash
pip install -r methods/quad_207422_milo_cldn4/requirements.txt
python3 methods/quad_207422_milo_cldn4/scripts/download_gse131907.py
python3 methods/quad_207422_milo_cldn4/scripts/download_gse148071.py
python3 methods/quad_207422_milo_cldn4/scripts/run_gse148071.py \
  --outdir methods/quad_207422_milo_cldn4/results/GSE148071 \
  --finding methods/quad_207422_milo_cldn4/results/GSE148071/FINDING.md
python3 methods/quad_207422_milo_cldn4/scripts/run_gse131907.py \
  --outdir methods/quad_207422_milo_cldn4/results/GSE131907 \
  --finding methods/quad_207422_milo_cldn4/results/GSE131907/FINDING.md
python3 methods/quad_207422_milo_cldn4/scripts/merge_quad.py
```
