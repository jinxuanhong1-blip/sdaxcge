# CosMx-native density fields: tumor CLDN4 vs predicted CD8

He et al. 2022 CosMx SMI NSCLC has no matched Visium companion, so Tangram, RCTD, and cell2location were not run. Author cell types were mapped onto edge-corrected Gaussian density fields, and tumor CLDN4 was regressed on the CD8 field.

The pre-specified 50 µm mean Spearman is slightly negative and sits outside a within-FOV label null. The eight sections do not share that sign. Lung9 and Lung12 carry the mean. Lung5 and Lung6 are weakly positive, on a near-zero CD8 floor. A section-rank test is compatible with no shift (Wilcoxon p = 0.38). Bandwidth 100 µm is weaker and is not separated from its null.

This is a different estimand from the locked marker-based neighbor counts and contact odds. Those results are unchanged.

## Visium gate

The public He et al. *Nat Biotechnol* 2022 release (NanoString SMI flat files and the Giotto object; 8 sections from 5 tissues) is CosMx SMI only. The paper’s data-availability statement points at the CosMx dataset. There is no Visium, Visium HD, or other spot assay on these sections. Unrelated public Visium LUAD cohorts are separate tissues and were not used as a companion.

Tangram, RCTD, and cell2location map a single-cell reference onto spot counts. With no spot assay, that branch stops. The CosMx-native substitute is below: the author labels already exist at single-cell resolution, and the “predicted CD8” field is the Gaussian intensity of those labels.

## Data

| Item | Value |
|---|---|
| Source | Public Giotto object, He et al. 2022 (author cell types, CLDN4 counts, millimetre coordinates) |
| Cells | 771,236 author-QC cells, 233 FOVs, 8 sections, 5 patients |
| Author tumor | 302,389 (`tumor 5/6/9/12/13`) |
| Author CD8 | 16,070 (`T CD8 naive` + `T CD8 memory`) |
| CLDN4 | On the 960-plex (`data/cosmx_960_genes.csv`) and present as counts in the cell table |
| Coordinates | Giotto millimetres × 1000 = µm. Density is computed inside each FOV; FOV boxes do not overlap |
| Input digest | `d3add2b1f9e54083cd5e55547ef728891f7fb695d95afb4cffd68c7ef8e3a573` |

Tumor CLDN4 is zero-inflated (45.5% zeros overall; 75.5% in Lung6, where the section 75th percentile is 0).

## Method

Pre-specified before the association was computed.

- Gaussian kernel, bandwidths **50 µm (primary)** and 100 µm, pixel 2 µm, truncate 4σ.
- Edge correction divides by the kernel mass inside the FOV rectangle. Primary cells are those with kernel mass ≥ 0.80 (234,809 / 302,389 tumor cells at 50 µm). A mask of 0.20 is a sensitivity check.
- **Predicted CD8 density** = edge-corrected intensity of author CD8 cells, in cells per (100 µm)², evaluated at each tumor cell. The index tumor cell is not a CD8 cell.
- **Predicted CD8 fraction** = CD8 intensity / leave-one-out total intensity. This is the continuous analog of an RCTD/cell2location proportion. Leave-one-out removes the index cell from the tumor and total fields.
- The 16 super-types are a partition of the author labels. On the primary mask their intensities sum to the total field (median ratio 1.000).
- **Primary regression:** within each section, Spearman ρ of tumor CLDN4 vs predicted CD8 density. The test statistic is the unweighted mean of the 8 section ρ values.
- **Companion contrast:** CLDN4-high = count ≥ section 75th percentile and ≥ 1; CLDN4-undetected = count 0. Statistic = mean of the 8 differences in median predicted CD8.
- **Control:** partial Spearman of CLDN4 vs predicted CD8 given leave-one-out tumor density.
- **Null:** 999 permutations that shuffle CLDN4 inside each FOV and leave coordinates and cell types fixed (seed 20260921). Two-sided p = (1 + count of \|null\| ≥ \|observed\|) / 1000. The smallest attainable p is 0.001. Cell-level Spearman p-values are not used; cells in a FOV are spatially dependent.
- Kernel check against the continuous Gaussian: interior median relative error 0.003, correlation 1.000. A planted left–right CD8 pattern recovered ρ = −0.88.

## Primary result, 50 µm

Mean section Spearman ρ = **−0.042**. None of 999 within-FOV shuffles matched that distance from 0 (two-sided permutation p = **0.001**).

The median section ρ is **−0.009**. Four sections are negative and four are positive. Wilcoxon signed-rank on the eight ρ values, testing a shift of the section median, is p = **0.38**. The mean is outside the null because Lung9_Rep2 (ρ = −0.195) and Lung12 (ρ = −0.135) are large; Lung5 and Lung6 are small and positive.

| Section | Patient | Masked tumor | High / undetected | Median CD8 density high | Median CD8 density undetected | Δ density | Spearman ρ | Partial ρ given tumor density |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Lung5_Rep1 | Lung5 | 14,371 | 4,666 / 4,196 | 0.018 | 0.019 | −0.0002 | +0.005 | −0.025 |
| Lung5_Rep2 | Lung5 | 14,778 | 4,037 / 4,734 | 0.029 | 0.021 | +0.008 | +0.015 | −0.010 |
| Lung5_Rep3 | Lung5 | 12,730 | 3,298 / 4,713 | 0.023 | 0.012 | +0.011 | +0.055 | +0.026 |
| Lung6 | Lung6 | 52,162 | 13,469 / 38,693 | 0.009 | 0.008 | +0.001 | +0.009 | +0.024 |
| Lung9_Rep1 | Lung9 | 30,462 | 10,829 / 8,175 | 0.268 | 0.379 | −0.110 | −0.069 | +0.017 |
| Lung9_Rep2 | Lung9 | 75,155 | 28,055 / 28,841 | 0.089 | 0.277 | −0.187 | −0.195 | −0.146 |
| Lung12 | Lung12 | 14,926 | 6,039 / 8,887 | 0.143 | 0.254 | −0.111 | −0.135 | −0.075 |
| Lung13 | Lung13 | 20,225 | 5,573 / 4,730 | 1.836 | 1.889 | −0.053 | −0.023 | −0.001 |

Δ is median predicted CD8 density in CLDN4-high minus the median in CLDN4-undetected cells, in cells per (100 µm)². Lung6’s high arm is CLDN4 ≥ 1 because that section’s 75th percentile is 0.

Mean Δ density = **−0.055** (permutation p = 0.001). Wilcoxon p on the eight differences = **0.25**. Lung5_Rep1’s Δ of −0.0002 is a numerical tie, so five negative signs are not five real deficits.

Patient-level mean Δ (replicates averaged first):

| Patient | Mean Δ density | Mean section ρ |
|---|---:|---:|
| Lung5 | +0.006 | +0.025 |
| Lung6 | +0.001 | +0.009 |
| Lung9 | −0.149 | −0.132 |
| Lung12 | −0.111 | −0.135 |
| Lung13 | −0.053 | −0.023 |

Three of five patients are negative. The two-sided sign test is p = 1.0. The permutation p of the patient-mean Δ is 0.001 because Lung9 and Lung12 are large, which is the same concentration the section mean already shows.

Within-FOV Spearman (213 FOVs with ≥30 masked tumor cells): 131/213 (62%) are negative. Lung9_Rep1 18/20, Lung9_Rep2 40/45, Lung12 23/28. Lung5_Rep1 8/24, Lung6 11/29.

Author-CD8 density around tumor is near a floor in Lung5 and Lung6 (medians about 0.01–0.03 cells per (100 µm)²). Those sections have little CD8 mass for a kernel to move. Lung13’s median is about 1.8 and the CLDN4 contrast there is small.

The CD8-fraction field gives the same split: mean section ρ = −0.043 (permutation p = 0.001), 4/8 sections negative, Wilcoxon p = 0.38. CD8+NK density (author CD8 plus NK) has the same four negative sections; the Lung5 patient-level Δ for CD8+NK is positive (+0.067).

Edge-mask sensitivity (kernel mass ≥ 0.20, all 302,389 tumor cells): mean ρ = −0.036, and the same four sections are negative.

## Bandwidth 100 µm

Mean section ρ = **−0.017** (permutation p = **0.076**). Mean Δ density = −0.017 (permutation p = 0.29). Partial ρ given tumor density = −0.008 (p = 0.22). Wilcoxon p on the eight ρ values = 0.95. The 50 µm mean does not recur as a 100 µm result. Lung9_Rep2 (ρ = −0.198) and Lung12 (ρ = −0.062) are still the negative sections; Lung9_Rep1 is about zero (−0.011).

## Tumor-density control

Mean partial Spearman at 50 µm, CLDN4 vs CD8 density given local tumor density, is **−0.024** (permutation p = 0.001). The crude mean shrinks from −0.042 to −0.024.

The partial is not uniform either. Lung9_Rep2 stays negative (−0.146) and Lung12 stays negative (−0.075). Lung9_Rep1 flips from −0.069 to +0.017, so that section’s crude CD8 association travels with tumor density. Lung5_Rep3 and Lung6 stay weakly positive.

Across the type panel, mean Spearman of CLDN4 with each density field at 50 µm:

| Field | Mean ρ |
|---|---:|
| Endothelial | −0.069 |
| Epithelial (non-tumor) | −0.054 |
| CD8 | −0.042 |
| B | −0.042 |
| Neutrophil | −0.040 |
| Plasmablast | −0.039 |
| Treg | −0.037 |
| Mast | −0.035 |
| Fibroblast | −0.034 |
| NK | −0.033 |
| CD4 | −0.015 |
| pDC | −0.007 |
| Monocyte | +0.017 |
| Macrophage | +0.023 |
| mDC | +0.030 |
| Tumor (leave-one-out) | +0.051 |

In Lung9_Rep2 and Lung12, CLDN4 is positively associated with tumor density (ρ = +0.170 and +0.209) and negatively associated with CD8, NK, B, neutrophil, and endothelium. In those two sections, CLDN4-high tumor cells sit in tumor-denser neighborhoods that contain less of several non-tumor types. CD8 is one of those types. Endothelium’s mean ρ is more negative than CD8’s. In Lung5 the tumor-density correlation is negative (about −0.08 to −0.10) and the CD8 correlation is weakly positive: the geography runs the other way.

## Figures

- `results/cosmx_density_field/figures/section_paired_cd8.png` — per-section high−undetected difference at 50 µm.
- `results/cosmx_density_field/figures/spearman_forest_h50.png` — section ρ against the within-FOV null.
- `results/cosmx_density_field/figures/cldn4_vs_predicted_cd8_deciles.png` — mean CLDN4 by decile of predicted CD8 density.
- `results/cosmx_density_field/figures/type_density_heatmap.png` — CLDN4 vs every super-type field.
- `results/cosmx_density_field/figures/example_fov_cd8_density.png` — Lung9_Rep1 FOV 10, the FOV whose ρ (−0.053) equals the FOV median.
- `results/cosmx_density_field/figures/permutation_null_mean_rho.png` — null of the 8-section mean ρ.

Tables: `results/cosmx_density_field/tables/section_h50.csv`, `section_h100.csv`, `type_spearman_h50.csv`, `fov_spearman_h50.csv`, `summary.json`.

## What this analysis is

A CosMx-native map of author cell types onto density fields, plus a regression of tumor CLDN4 on predicted CD8. The regression mean at 50 µm is slightly negative relative to a within-FOV null, concentrated in Lung9 and Lung12, and absent in Lung5 and Lung6. It is heterogeneous across the five patients, sensitive to bandwidth, and shared with other non-tumor fields in the sections where it appears.

Neighbor counts, contact odds, and Ripley g(r) on marker-defined cytotoxic cells are separate estimands and are not recomputed here. Effector-state ratios (GZMB, PRF1, NKG7, IFNG) are not in this model. No ICI labels and no private 8-KL data were used.

## Reproduce

```bash
python3 scripts/test_density_kernel.py
python3 scripts/cosmx_density_field_cldn4.py
```

Dependencies: `scripts/requirements-cosmx-density.txt`. The run evaluates 233 FOVs and 2 × 999 permutations. The cell-level dump `tumor_cell_density.csv.gz` is written locally and gitignored.
