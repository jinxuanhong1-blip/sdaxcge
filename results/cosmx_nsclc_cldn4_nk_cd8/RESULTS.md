# CosMx NSCLC CLDN4-only: CD8+NK neighbor counts and mixing

ADDITIVE, **CLDN4-only**, official CosMx NSCLC 960-plex (He et al. 2022): **all 8 sections / 5 patients**. No private 8-KL. Primary readout is **CD8+NK neighbor COUNT** and **mixing** on a **20 / 40 / 60 / 80 µm** grid, with a **section-level paired test (n=8)**. Nearest-distance is not the primary endpoint. Tumor cells are **not** gated as `CD8A==0`.

CLDN4 is **present** on the 960-plex in all 8 sections.

## Design (retuned)

- **Index:** tumor / epithelial RNA-compartment cells, including those with CD8A>0.
- **CLDN4-high / low:** Q4 vs Q1 of log-normalized CLDN4 among tumor cells, per section.
- **Neighbors:** CD8+NK together (July PPT cytotoxic): non-tumor CD8A/B+ **or** NKG7/GNLY/KLRD1+ CD3− (KLRD1 is not on this panel; NKG7/GNLY used).
- **Radii:** 20, 40, 60, 80 µm. Within-FOV only. 0.18 µm/pixel.
- **Mixing:** Keren immune-side (CD8+NK–tumor contacts / CD8+NK–CD8+NK contacts) and homogeneous mixing, pooled across FOVs within each section.
- **Inference:** paired Wilcoxon signed-rank on the 8 section means (high vs low). Honest n = 8 sections (5 patients).

## Inventory

| Section | Patient | QC cells | FOVs | Tumor | Tumor CD8A+ kept | CLDN4-high | CLDN4-low | CD8+NK |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Lung5_Rep1 | P5 | 98,046 | 30 | 40,546 | 2,786 | 10,142 | 23,959 | 13,315 |
| Lung5_Rep2 | P5 | 100,437 | 29 | 40,601 | 3,702 | 10,209 | 23,488 | 17,181 |
| Lung5_Rep3 | P5 | 97,876 | 30 | 39,182 | 2,088 | 9,811 | 25,573 | 12,572 |
| Lung6 | P6 | 90,126 | 30 | 73,849 | 4,831 | 17,592 | 56,257 | 2,541 |
| Lung9_Rep1 | P9 | 87,617 | 20 | 46,604 | 3,878 | 12,735 | 19,996 | 7,838 |
| Lung9_Rep2 | P9 | 139,544 | 45 | 92,858 | 7,202 | 23,343 | 39,248 | 8,841 |
| Lung12 | P12 | 71,440 | 28 | 30,124 | 3,731 | 7,575 | 21,278 | 11,860 |
| Lung13 | P13 | 81,248 | 20 | 33,221 | 2,457 | 8,827 | 13,346 | 9,098 |

Tumor cells with CD8A>0 were **kept** in the tumor index (n=30,675 across sections). The CLDN4-low arm is larger than a strict quartile because CosMx CLDN4 is zero-inflated (Q1 = 0 in these sections); Q4 remains the top quartile.

## Section-level paired CD8+NK neighbor counts

| Radius | Median section mean (high) | Median section mean (low) | Median Δ (high−low) | Sections high<low | Wilcoxon p |
|---:|---:|---:|---:|---:|---:|
| 20 | 0.275 | 0.653 | -0.365 | 8/8 | 0.008 |
| 40 | 1.521 | 2.688 | -1.203 | 8/8 | 0.008 |
| 60 | 3.765 | 6.090 | -2.219 | 8/8 | 0.008 |
| 80 | 7.061 | 10.823 | -3.301 | 8/8 | 0.008 |

Per-section means:

| Section | Patient | Δ20 | Δ40 | Δ60 | Δ80 | mean high@40 | mean low@40 |
|---|---|---:|---:|---:|---:|---:|---:|
| Lung5_Rep1 | P5 | -0.389 | -1.252 | -2.322 | -3.562 | 1.550 | 2.802 |
| Lung5_Rep2 | P5 | -0.484 | -1.507 | -2.757 | -4.056 | 2.276 | 3.783 |
| Lung5_Rep3 | P5 | -0.367 | -1.223 | -2.312 | -3.559 | 1.334 | 2.557 |
| Lung6 | P6 | -0.021 | -0.055 | -0.088 | -0.126 | 0.279 | 0.334 |
| Lung9_Rep1 | P9 | -0.444 | -1.674 | -3.426 | -5.458 | 1.492 | 3.165 |
| Lung9_Rep2 | P9 | -0.201 | -0.703 | -1.374 | -2.136 | 1.047 | 1.750 |
| Lung12 | P12 | -0.363 | -1.183 | -2.125 | -3.043 | 1.774 | 2.957 |
| Lung13 | P13 | -0.335 | -0.885 | -1.197 | -1.106 | 1.689 | 2.574 |

All **8/8** sections have a lower mean CD8+NK neighbor count around CLDN4-high than around CLDN4-low at every radius. Wilcoxon p = 0.008 is the two-sided signed-rank minimum for n=8 with no sign ties. Lung6 is immune-poor overall (few CD8+NK) but the Δ still has the same sign.

Figures: `results/cosmx_nsclc_cldn4_nk_cd8/figures/ecdf_cd8nk_counts.png`, `section_paired_cd8nk_counts.png`.

## Section-level paired mixing (CLDN4-high/low tumor vs CD8+NK)

| Radius | Keren immune Δ (high−low) | Wilcoxon p | Homog Δ | Wilcoxon p |
|---:|---:|---:|---:|---:|
| 20 | -1.455 | 0.008 | -0.149 | 0.016 |
| 40 | -1.484 | 0.008 | -0.130 | 0.016 |
| 60 | -1.506 | 0.008 | -0.112 | 0.016 |
| 80 | -1.514 | 0.008 | -0.095 | 0.023 |

Per-section Keren mixing at 40 µm (immune-side):

| Section | vs CLDN4-high | vs CLDN4-low | Δ |
|---|---:|---:|---:|
| Lung5_Rep1 | 0.394 | 1.681 | -1.288 |
| Lung5_Rep2 | 0.388 | 1.485 | -1.096 |
| Lung5_Rep3 | 0.379 | 1.892 | -1.513 |
| Lung6 | 1.953 | 7.482 | -5.529 |
| Lung9_Rep1 | 1.052 | 3.504 | -2.452 |
| Lung9_Rep2 | 2.119 | 5.955 | -3.836 |
| Lung12 | 0.395 | 1.849 | -1.454 |
| Lung13 | 0.514 | 1.184 | -0.670 |

Figures: `section_paired_mixing.png`, `section_paired_mixing_grid.png`.

## FOV forest of Δ count (secondary)

At 40 µm, 232 FOVs; median Δ mean CD8+NK count = -0.663; Wilcoxon signed-rank p = 8.10e-35; fraction Δ<0 = 0.905.

Figure: `fov_forest_delta_cd8nk_count_40um.png`. This is **not** a nearest-distance forest.

## What this does not claim

- Not a nearest-CD8 / nearest-NK distance primary analysis.
- Not ICI response and not private 8-KL.
- Marker CD8+NK is not the paper’s 18-type map (no author cell-type column in the public flat files).
- Q4 vs Q1 is a contrast, not a biological threshold. n=8 sections, not n=8 patients.

## Reproduce

```bash
python3 scripts/download_cosmx_nsclc.py
python3 scripts/cosmx_nsclc_cldn4_nk_cd4_radii.py
```

Raw tarballs are gitignored. Tables live under `results/cosmx_nsclc_cldn4_nk_cd8/tables/`.
