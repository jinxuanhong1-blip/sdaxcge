# Finding — concordant-4 cells move toward CLDN4-high along an expression trajectory

CLDN4-only. Units are the locked 65: GSE123902 (13), GSE131907 (21), GSE205335 (22), GSE189357 (9). Epithelial cells only (PTPRC = 0 and EPCAM or a KRT8/18/19 count > 0; GSE131907 and GSE205335 use the author epithelial label). Cap 160 tumor cells per unit. Uninvolved-lung epithelial cells are a root pool only. They are not among the 65.

This is an **expression trajectory**. It is not splicing RNA velocity. scVelo was not fit, because no public spliced/unspliced layer exists and a full FASTQ rebuild was not completed on this machine (see below).

## Headline (all 65 units)

PAGA shortest-path pseudotime on a Harmony embedding (dataset as batch). Root cell is the highest Hallmark-IFN cell among cells at or below the median CLDN4: `GSE205335:P1062`, CLDN4 = 0. Graph: 1,000 HVGs, 20 PCs, 30 neighbors. Readout scores are AUCell. CLDN4 itself is the expression value, not an AUCell set. Barrier genes are KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (CLDN4 excluded).

| Contrast vs PAGA pseudotime | n | ρ | p | BH q across 1,920 grid rows |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 | 65 | 0.660 | 2.20×10⁻⁹ | 7.1×10⁻⁸ |
| barrier (no CLDN4) | 65 | 0.575 | 5.61×10⁻⁷ | 3.0×10⁻⁶ |
| IFN (Hallmark α ∪ γ, AUCell) | 65 | −0.315 | 0.0107 | 0.038 |
| TJ (no CLDN4) | 65 | 0.512 | 1.29×10⁻⁵ | — |
| APM | 65 | −0.250 | 0.045 | — |

All four cohorts are in the 65. Later pseudotime has higher CLDN4, higher barrier, higher TJ, lower IFN, and lower APM. That is the thesis direction. Bonferroni across 1,920 rows does not leave the IFN p under 0.05 (0.0107 × 1920 > 1). CLDN4 and barrier on this row survive that Bonferroni. The BH q values above are the rank-adjusted values inside this grid, not a second cohort.

Figure: `results/figures/fig_paga_ifnroot_patient.png`. UMAP: `results/figures/fig_umap_cldn4_pt.png`.

## Palantir (same IFN-high root)

Palantir 1.4.5, 400 waypoints, early cell `GSE205335:P1062`, Harmony PCA, 30 neighbors. n = 65, all four cohorts.

| Contrast vs Palantir pseudotime | n | ρ | p |
| --- | ---: | ---: | ---: |
| CLDN4 | 65 | 0.646 | 6.04×10⁻⁹ |
| barrier AUCell | 65 | 0.563 | 1.08×10⁻⁶ |
| IFN AUCell | 65 | −0.292 | 0.0184 |

Figure: `results/figures/fig_palantir_patient.png`. This is multiscale diffusion, not a spliced/unspliced velocity.

## External root (not chosen as IFN-high)

Same 65 units. Diffusion pseudotime from uninvolved-lung AT2 `GSE131907:LUNG_N34` (CLDN4 = 0, AT2 z-score 3.44). 2,000 HVGs, 40 PCs, 15 neighbors, no Harmony. Scores are UCell.

| Contrast vs DPT | n | ρ | p |
| --- | ---: | ---: | ---: |
| CLDN4 | 65 | 0.322 | 0.0088 |
| barrier | 65 | 0.300 | 0.015 |
| IFN | 65 | −0.459 | 1.18×10⁻⁴ |

Same direction, with a root that was not the IFN maximum. Figure: `results/figures/fig_dpt_at2root_patient.png`.

## Strongest lineage-restricted clock

Slingshot-style MST, root = IFN-high / CLDN4-low, terminal cluster = highest barrier (CLDN4 not used to name the leaf). AddModuleScore. n = 40 units that have cells on that path (not 65). CLDN4 ρ = 0.678, p = 1.52×10⁻⁶. IFN ρ = −0.505, p = 9.0×10⁻⁴. Fisher one-sided p = 1.4×10⁻²⁶. Off-path units are missing, so this is not the full concordant-4 n.

A CLDN4-named leaf is stronger still (IFN ρ = −0.899, p = 2.3×10⁻¹³, n = 35, AddModuleScore) and is partly circular for the CLDN4 endpoint, because the leaf was the cluster with the highest CLDN4. It is in the grid. It is not the headline.

## What did not match

CytoTRACE-like potency (genes detected, smoothed on the neighbor graph; pseudotime = −potency) never matched. Median ρ(CLDN4, PT) = −0.50. On that clock, CLDN4-high cells look less differentiated, not like a sink. 0 of 128 CytoTRACE rows matched the sign pattern.

Within a unit, CLDN4-high vs CLDN4-low cells (quartile and the other cuts) have higher barrier and higher TJ (z-mean barrier median Δ = 0.210, Wilcoxon greater p = 1.26×10⁻¹², n = 65). They do **not** have lower IFN. The z-mean IFN difference is positive (median Δ = 0.086 at the quartile cut, one-sided “less” p ≈ 1). Immune-cold along the clock is a between-unit pseudotime pattern. It is not a within-unit CLDN4-high vs CLDN4-low IFN drop.

Cluster absorption onto a CLDN4-defined terminal recovers CLDN4 (ρ = 0.312, p = 0.011, n = 65) and does not recover an IFN drop (ρ ≈ 0, p = 0.999).

## Grid

1,920 trajectory rows: 1,000 or 2,000 HVGs, CLDN4 in or out of the graph, 20 or 40 PCs, 15 or 30 neighbors, Harmony on or off, five clocks (DPT, PAGA, MST, CytoTRACE-like, absorption), three roots, four scores (z-mean, UCell, AUCell, AddModuleScore). 641 rows match the sign pattern CLDN4 up, barrier up, IFN down. 629 of those use all four datasets. Full table: `results/tables/trajectory_grid.tsv`. Best row per method: `results/tables/best_by_method.tsv`.

24 within-unit paired rows (4 scores × 6 cuts). One row has a negative IFN median (AddModuleScore, median split, Δ = −0.00042, p = 0.29). Barrier and TJ are positive in every score. Table: `results/tables/geneset_paired_grid.tsv`.

## Splicing velocity

Deposited matrices are still one total-count table each. No loom. Public sralite reads exist for GSE123902 and GSE189357 only. GSE131907 and GSE205335 have no public SRA, so a loom rebuild cannot cover 43 of the 65 units. scVelo dynamical and stochastic models were not run.

## What this does not say

The locked patient-level CLDN4 versus T/NK result was not recomputed. Cell-level p-values are not the claim. Pseudotime is not a splicing arrow.
