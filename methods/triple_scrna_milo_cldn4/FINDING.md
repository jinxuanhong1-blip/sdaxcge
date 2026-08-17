# FINDING — triple-merge Milo vs malignant CLDN4

**No cohort-level neighbourhood DA vs malignant CLDN4.** At the honest sample-count floor (n_present ≥ 8), SpatialFDR < 0.1 is **0 / 0 / 0** on GSE131907 / GSE148071 / GSE205335. Stacked patient-level malignant CLDN4 vs T/NK is null (n = **71**, ρ = −0.066, p = 0.58). The 11 SpatialFDR < 0.1 neighbourhoods at the ≥5-sample floor are all |ρ| = 1 on 5–6 samples and are not a merge-level claim.

**Additive CLDN4-only.** Neighbourhood abundance versus **sample-level malignant CLDN4** on the public processed files from GSE131907 (Kim 2020; tLung epithelium+immune), GSE148071 (Wu 2021; 42 advanced-NSCLC biopsies) and GSE205335 (Ahn/Lee; non-normal ICI biopsies, patient-pooled). **No dual-high TACSTD2+CLDN4 gate. No GSE207422.**

The independent unit is the **patient / sample**, not the cell and not the overlapping neighbourhood. Do not write n = 213,453.

**Not miloR.** Graph + SpatialFDR follow Dann et al. 2022; DA is a sample-level Spearman (primary) or Welch t-test (median split). Joint Harmony of the three matrices was **not** run (15 GB RAM; ~370k cells). Each dataset has its own kNN graph with per-sample PCA centering.

## Verdict

| Dataset | Honest n (scored / in graph) | testable / nhoods | P<0.05 | min SpatialFDR | SpatialFDR<0.1 | sample CLDN4–T/NK ρ (p) |
|---|---|---|---:|---|---:|---|
| GSE131907 | **10** / 11 | 210 / 2998 | 17 | 0.000 | **3** | -0.418 (2.29e-01) |
| GSE148071 | **39** / 42 | 881 / 6401 | 81 | 0.000 | **8** | -0.127 (4.43e-01) |
| GSE205335 | **22** / 22 | 360 / 5713 | 17 | 0.899 | **0** | -0.263 (2.38e-01) |
| **stacked patients** | **71** | — | — | — | — | -0.066 (5.84e-01) |

GSE131907: 3 SpatialFDR<0.1 nhoods (all |ρ|=1 on 5 samples; T/NK max=31 — T-cell-only neighbourhoods at the floor). GSE148071: 8 SpatialFDR<0.1 nhoods (all |ρ|=1 on 5–6 samples; T/NK max=0). GSE205335: 0.

A SpatialFDR<0.1 hit that is |ρ|=1 on the 5-sample testability floor is **not** a cohort DA claim. At n_present ≥ 8 the three graphs are the numbers to quote.

## Honest n

| item | n | note |
|---|---:|---|
| Datasets | **3** | GSE131907 tLung + GSE148071 + GSE205335; GSE207422 excluded |
| Samples in the three graphs | **75** | the DA unit |
| Samples with malignant CLDN4 (≥10 malignant cells) | **71** | primary Spearman *n* |
| Cells in the three graphs | 213,453 | **do not cite as n** |
| Neighbourhoods (stacked) | 15112 | overlapping; not independent |
| Testable nhoods (stacked, still per-graph FDR) | 1451 | ≥5 samples present **inside that dataset** |

GSE131907 uses **tLung only** so brain mets / PE / mLN are not mixed into the lung-tumor graph. GSE205335 drops normal LN/brain/lung and pools biopsies by patient. GSE148071 is one biopsy per patient; malignant-like is a marker proxy (author CNA IDs are not public).

## SpatialFDR sensitivity (min samples present)

| Dataset | floor | testable | P<0.05 | min SpatialFDR | SpatialFDR<0.1 | \|ρ\|=1 |
|---|---|---:|---:|---|---:|---:|
| GSE131907 | ≥5 | 210 | 17 | 0.000 | **3** | 3 |
| GSE131907 | ≥8 | 9 | 0 | 0.834 | **0** | 0 |
| GSE131907 | ≥10 | 1 | 0 | 0.933 | **0** | 0 |
| GSE131907 | ≥15 | 0 | 0 | NA | **0** | 0 |
| GSE148071 | ≥5 | 881 | 81 | 0.000 | **8** | 8 |
| GSE148071 | ≥8 | 192 | 25 | 0.178 | **0** | 0 |
| GSE148071 | ≥10 | 48 | 3 | 0.218 | **0** | 0 |
| GSE148071 | ≥15 | 0 | 0 | NA | **0** | 0 |
| GSE205335 | ≥5 | 360 | 17 | 0.899 | **0** | 0 |
| GSE205335 | ≥8 | 65 | 3 | 0.899 | **0** | 0 |
| GSE205335 | ≥10 | 16 | 0 | 0.937 | **0** | 0 |
| GSE205335 | ≥15 | 0 | 0 | NA | **0** | 0 |

## Composition (transcriptional kNN, not histology)

Epithelium clusters with epithelium. Unrestricted nhood CLDN4 vs T/NK is lineage geometry. The sample-paired Wilcoxon is the same geometry at patient resolution. Disjoint-interface and sample-level Spearman are the honest composition tests.

| Dataset | interface n | iface ρ (p) | disjoint iface | disjoint ρ (p) | paired T/NK high vs low |
|---|---:|---|---|---|---|
| GSE131907 | 17 | -0.091 (7.29e-01) | 2 / 266 | NA (NA) | 11; med 0.015 vs 0.031; p=1.00e+00 |
| GSE148071 | 37 | -0.274 (1.01e-01) | 4 / 559 | -0.949 (5.13e-02) | 42; med 0.000 vs 0.000; p=1.67e-01 |
| GSE205335 | 362 | -0.120 (2.19e-02) | 32 / 553 | -0.120 (5.12e-01) | 22; med 0.002 vs 0.057; p=2.19e-04 |

Stacked sample-level malignant CLDN4 vs T/NK: n=71, ρ=-0.066, p=5.84e-01. CLDN4 scales are not Harmony-aligned across studies; this is a mixed-study rank correlation of independently computed scores.

## What this does not say

1. It does not invent an ICI / MPR contrast on GSE131907 or GSE148071 (GEO has none). GSE205335 RECIST is not re-tested here.
2. It does not treat 213,453 cells as *n*.
3. Neighbourhood “next to” is kNN co-membership, not histology.
4. miloR / edgeR QLF numbers are not claimed.
5. Dual-high TACSTD2+CLDN4 was not used. TACSTD2 is not a gate.
6. GSE207422 was not downloaded or analysed.
7. Joint Harmony was not feasible from the public processed files at this memory budget; per-dataset graphs are the designed fallback.
8. mRNA ≠ protein.

This is extra weight on whether **neighbourhood DA vs malignant CLDN4** appears when the three public lung scRNA matrices are analysed additively. It is not a GSE207422 MPR replicate and it is not a spatial-niche paper.

## Files

- `tables/nhoods.tsv` (stacked; dataset column)
- `tables/sample_scores.tsv`, `tables/honest_n.tsv`, `tables/one_row.tsv`, `tables/summary.json`
- per-dataset folders under `tables/` and `figures/`
- extra figures: `figures/fig_triple_volcano.png`, `fig_patient_cldn4_vs_tnk.png`, `fig_honest_n.png`, `fig_testable_nhoods.png`, `fig_min_spatialfdr.png`, `fig_interface_cldn4_vs_tnk.png`, `fig_sample_paired_tnk.png`, `fig_rho_hist_testable.png`

```bash
python3 methods/triple_scrna_milo_cldn4/scripts/download.py
python3 methods/triple_scrna_milo_cldn4/scripts/analyze.py
```
