# RESULTS — public single-cell spatial CLDN4 vs CD8 (additive)

CLDN4-only. No TACSTD2∩CLDN4 dual-high. No private 8-KL. No proxy gene when CLDN4 is absent.

Question: among **tumor cells** on public imaging-ST NSCLC cores, are **CLDN4-high** cells farther from the nearest **CD8A+** cell, and do they have fewer CD8A+ neighbors at 25 / 50 / 100 µm, than CLDN4-low tumor cells? **KRT8** and **EPCAM** are the same-split controls.

Confirmatory unit = **FOV / core**, not cell.

**ICON CosMx (n=69 FOVs): no CLDN4-specific CD8 distance or radius deficit.** Per-FOV median µm to nearest CD8A+ is 21.4 (CLDN4-high) vs 20.6 (CLDN4-low), paired Δ=−0.09 µm, p=0.87. CD8 counts at 25 / 50 / 100 µm are identical at the FOV median (1 / 3 / 13). The two LUAD TMAs have **opposite** small signs, and those signs track the KRT8 / EPCAM controls.

---

## Panel check (required before any distance)

| Dataset | Platform / panel | CLDN4 | CD8A | KRT8 | EPCAM | Action |
|---|---|---|---|---|---|---|
| Nat Commun 2025 paper list | CosMx Universal 1000 | **yes** | yes | yes | yes | analyze GSE299786 ICON |
| Nat Commun 2025 paper list | Xenium lung 289 + 50 custom | **no** (CLDN5 only) | yes | no | yes | record; no proxy |
| Nat Commun 2025 paper list | MERFISH IO 500 | **no** (CLDN5 only) | yes | no | yes | record; no proxy |
| GSE299786 matrix (LUAD TMA1) | CosMx 1000 | **yes** | yes | yes | yes | distances |
| GSE299886 matrix (LUAD TMA1) | MERFISH IO 500 | **no** | yes | no | yes | no distances |
| GSE300007 matrix (LUAD TMA2 UM) | Xenium | **no** | yes | no | yes | no distances |
| GSE311609 NSCLC Prime 5K | Xenium | **no** (CLDN1/5/7/18) | yes | no | yes | no distances |
| GSE311609 NSCLC custom IO | Xenium | **no** (no CLDN gene) | yes | yes | yes | no distances |
| GSE311609 NSCLC lung 289 | Xenium | **no** (CLDN5 only) | yes | no | yes | no distances |

GSE311609 is public. Compact `cell_feature_matrix.h5` + `cells.parquet` were downloaded per NSCLC panel (not the 179 GB RAW tar; breast skipped). CLDN4 is absent on all three NSCLC panels. No proxy.

GEO names LUAD TMA1/TMA2 are the paper ICON NSCLC cores. Paper ICON1 Xenium is not in GSE300007 (only LUAD TMA2 + MESO).

---

## Methods (locked)

- **Tumor cell:** `EPCAM≥1` or `KRT8≥1` or `KRT19≥1`, and `CD8A=0`.
- **CD8 T:** `CD8A≥1` (no author cell-type column in the public CosMx flat files).
- **CLDN4-high vs low:** within each FOV, among tumor cells. If the FOV median CLDN4 count is 0, high = `CLDN4≥1` and low = `CLDN4=0`. Otherwise a median split. ICON split rule: 55 FOVs pos-vs-zero, 14 median-split.
- **KRT8 / EPCAM controls:** the same within-FOV split on that marker, same tumor and CD8 definitions.
- **Geometry:** CosMx `CenterX/Y_global_px` × 0.12028 µm/px (NanoString; metadata has `Area` in pixels but no `Area.um2`). Neighborhoods are **within FOV only**.
- **QC:** `cell_ID>0`, RNA counts ≥ 20. FOV kept if ≥10 tumor cells per arm and ≥5 CD8A+ cells.
- **Confirmatory test:** paired Wilcoxon (Pratt) on per-FOV medians. ICON = LUAD TMA1 + TMA2 FOVs.
- Mesothelioma is secondary: CLDN4 is on the CosMx panel but rarely detected.

---

## GSE299786 CosMx — ICON / LUAD (primary)

QC cells: TMA1 37,419 (36 FOVs); TMA2 39,123 (38 FOVs). CLDN4 %pos in tumor: 30.6% (TMA1), 47.5% (TMA2). µm/px = 0.12028.

**n FOV (ICON, QC) = 69.** Tumor cells used: 13,917 CLDN4-high, 26,277 CLDN4-low. CD8A+ cells in those FOVs: 7,485.

| Test | CLDN4-high | CLDN4-low | high−low | n FOV | p |
|---|---:|---:|---:|---:|---:|
| Median µm to nearest CD8A+ | 21.4 | 20.6 | −0.09 | 69 | 0.87 |
| Median CD8A+ count in 25 µm | 1 | 1 | 0 | 69 | 0.22 |
| Median CD8A+ count in 50 µm | 3 | 3 | 0 | 69 | 0.57 |
| Median CD8A+ count in 100 µm | 13 | 13 | 0 | 69 | 0.11 |

FOVs where CLDN4-high tumor is farther from CD8: **33 / 69**. Closer: **36 / 69**.

Same FOV-level test for controls (ICON):

| Marker | median µm high | median µm low | Δ | n FOV | p | farther / closer |
|---|---:|---:|---:|---:|---:|---|
| CLDN4 | 21.4 | 20.6 | −0.09 | 69 | 0.87 | 33 / 36 |
| KRT8 | 21.3 | 20.5 | +0.09 | 68 | 0.66 | 37 / 31 |
| EPCAM | 21.2 | 21.2 | −0.57 | 69 | 0.15 | 25 / 44 |

### Per-TMA (do not pool the signs away)

| TMA | n FOV | CLDN4 Δ µm | p | KRT8 Δ µm | p | EPCAM Δ µm | p |
|---|---:|---:|---:|---:|---:|---:|---:|
| LUAD TMA1 | 32 | **+0.82** | 0.043 | +0.45 | 0.009 | −0.22 | 0.90 |
| LUAD TMA2 | 37 | **−0.87** | 0.019 | −0.47 | 0.075 | −0.66 | 0.042 |

TMA1 “CLDN4-high farther” is matched by **KRT8**. TMA2 “CLDN4-high closer” is matched by **EPCAM**. Neither is a CLDN4-only neighborhood effect.

### Distance ECDFs

![ICON FOV-level ECDF](figures/ecdf_ICON_fov_median_dist.png)

![LUAD TMA1 cell-level ECDF](figures/ecdf_LUAD_TMA1_cldn4.png)

![LUAD TMA2 cell-level ECDF](figures/ecdf_LUAD_TMA2_cldn4.png)

![ICON FOV box, nearest CD8](figures/box_ICON_fov_median_dist.png)

![ICON FOV box, CD8 counts at 25/50/100 µm](figures/box_ICON_cd8_radii.png)

### Example ICON FOV maps

![LUAD TMA1 FOV 17](figures/fov_LUAD_TMA1_FOV17.png)

![LUAD TMA2 FOV 6](figures/fov_LUAD_TMA2_FOV6.png)

Additional maps: `fov_LUAD_TMA1_FOV22.png`, `FOV43`, `FOV46`; `fov_LUAD_TMA2_FOV14.png`, `FOV20`, `FOV25`.

---

## Mesothelioma (secondary)

CLDN4 is on the CosMx panel, so MESO was scored. It is **not useful** as a CLDN4-high vs low tumor contrast: only 7.3% (TMA1) and 9.4% (TMA2) of RNA-epithelial cells are CLDN4+. Paired FOV Δ = −0.13 µm (p=0.81) and +0.37 µm (p=0.10). ECDFs: `ecdf_MESO_TMA1_cldn4.png`, `ecdf_MESO_TMA2_cldn4.png`.

---

## Datasets with CLDN4 absent (no distances)

- **GSE300007 Xenium** (lung 289 + 50 custom; LUAD TMA2 UM/MM). Features: CLDN5, CD8A, EPCAM, KRT19. **No CLDN4. No KRT8.**
- **GSE299886 MERFISH** (IO 500). Paper list and LUAD TMA1/TMA2 `cell_by_gene` headers: CLDN5, CD8A, EPCAM. **No CLDN4. No KRT8.**
- **GSE311609 Xenium** NSCLC compact matrices. Prime 5K: CLDN1/5/7/18, not CLDN4. Custom IO: no CLDN gene. Lung 289: CLDN5 only. Breast not used.

---

## What this does not claim

- It does not use private 8-KL spatial data.
- It does not impute CLDN4 from CLDN5, TACSTD2, or a TJ score.
- It does not treat cell-level CosMx p-values as the confirmatory n (TMA2 cell-level MW p=9.6×10⁻⁶ is the large-n artifact; FOV unit is p=0.019 and opposite to TMA1).
- It does not re-score prior Visium/GeoMx leftover ρ.

---

## Files

- `methods/public_spatial_cldn4/analyze.py`
- `results/public_spatial_cldn4/tables/panel_check.json`
- `results/public_spatial_cldn4/tables/cosmx_summary.json`
- `results/public_spatial_cldn4/tables/ICON_CLDN4_per_fov.csv`
- `results/public_spatial_cldn4/figures/ecdf_*.png`
- `results/public_spatial_cldn4/figures/fov_*.png`

Reproduce: `python3 methods/public_spatial_cldn4/analyze.py`
