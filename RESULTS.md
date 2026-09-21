# CosMx NSCLC: largest CLDN4 cytotoxic-neighbor effect that keeps the sign

Sensitivity layer on He et al. 2022 CosMx NSCLC (figshare 25976224, `cosmx_human_nsclc_clustered.h5ad`, 8 sections / 5 donors). The locked exclusion summary is unchanged: cytotoxic neighbor ratio **0.36 at 50 µm and 0.52 at 100 µm**, 8/8 sections, 5/5 donors, sign P = 0.031. Nothing below replaces that summary.

## Call

Largest ratio effect among specs with ≥7/8 sections and ≥4/5 donors of one sign: **detected_absent, cd8nk, 20 µm, section threshold**. Equal-weight section means 0.041 (high) / 0.071 (low), high/low ratio **0.573**, Δ **-0.030**. Sign: 7/8 sections and 4/5 donors (exclusion). |log ratio| = 0.557. One-sided sign P = 0.0352 (sections, out of 8) and 0.1875 (donors, out of 5). Two-sided Wilcoxon on usable section deltas P = 0.0156. FOVs with the same high<low direction: 156/217 (high>low 57/217).

Largest |Δ| under the same sign bar: **top30_vs_bottom, cd8, 50 µm, fov threshold**. Means 0.355 / 0.437, ratio 0.813, Δ **-0.082**. Sign: 7/8 sections and 4/5 donors (exclusion). Wilcoxon P = 0.0547.

Every eligible high/low ratio in this grid is below 1, and every eligible Δ is negative. The numerically largest ratio is 0.884 (top25_vs_rest, cd8, 50 µm, fov threshold) and the numerically smallest is 0.573 (detected_absent, cd8nk, 20 µm, section threshold). The least negative Δ is -0.001 (top20_vs_rest, nk, 10 µm, section threshold) and the most negative Δ is -0.082 (top30_vs_bottom, cd8, 50 µm, fov threshold). Eligible binary specs: 110 of 840.

Strongest fully concordant ratio (8/8 sections and 5/5 donors, low-arm mean ≥ 0.05): 0.634 at q4_q1, cd8nk, 20 µm, fov threshold, means 0.052 / 0.082, Δ -0.030. The call above is a larger fold and is not 8/8.

Neighborhoods whose equal-weight low-arm mean is below 0.05 do not get a ratio in the call. The smallest such high/low ratio with the sign bar is 0.226 (top20_vs_bottom, nk, 10 µm, fov threshold; means 0.0011 / 0.0048; 7/8 sections, 4/5 donors). That is a ratio of rare events, not the reported sensitivity.

Continuous CLDN4 (Spearman of log1p CP10k vs neighbor count), largest |mean ρ| with the same sign bar: **continuous, cd8, 40 µm, donor threshold**, mean section ρ = -0.043 (exclusion; 7/8 negative, 4/5 donors negative). This is not a high/low ratio and is not the count-sweep call.

## Selection rule

Eligible: ≥7 sections with Δ<0 and ≥4 donors with Δ<0, or the same counts with Δ>0. There are 8 sections and 5 donors; a section that cannot form both arms is not a supporting sign. Δ is the equal-weight mean of section (high − low) means. The ratio is the equal-weight mean of section high means divided by the equal-weight mean of section low means. It is left undefined when the low mean is below 0.05. Ratio effect is |log(ratio)|. Delta effect is |Δ|. Ties break on more sections in the majority direction, then section over FOV over donor. P-values describe the winning row. The grid was searched, so they are not a second locked test. Median split is in the grid as a calibration cut, not as a replacement for Q4/Q1.

## Definitions

- Malignant index: author `cell_type` matched to the section (`tumor 5/6/9/12/13`). CD8 labels: T CD8 memory, T CD8 naive. NK labels: NK. CD8+NK is their union. GZMB+: raw GZMB > 0 and cell_type does not start with 'tumor'.
- CLDN4 is log1p(count / n_counts × 10,000). Detected means raw count > 0. Thresholds are computed inside the named unit (section, FOV, or donor). Q4/Q1 keeps cells ≥ the 75th percentile versus cells ≤ the 25th percentile; overlap is removed. If the 75th percentile is 0, that group collapses to detected versus absent. Top X% versus bottom requires the upper quantile to sit strictly above the lower quantile. Top X% versus rest is cells at or above the upper quantile versus everyone below it. A section needs ≥30 cells in each arm. An FOV needs ≥5 in each arm, and a section needs ≥3 such FOVs.
- Counts are other cells of that cytotoxic class within the radius, in the same section, global centroids × 0.18 µm/px. The index cell is not in the reference set. Donor is not a spatial window: Lung5's three sections keep separate coordinate systems.
- Donor Δ averages the donor's usable sections with equal section weight. FOV signs use FOVs with ≥5 cells in each arm.

## Calibration row (not the lock)

Section-level median split, CD8+NK. This is the unstratified contrast that previously read about 0.80 at 50 µm on a closely related tumor definition. It is a pipeline check. It is not the locked 0.36 / 0.52 summary.

| Radius | High | Low | Ratio | Δ | Sections neg | Donors neg | Wilcoxon |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 0.005 | 0.010 | NA | -0.005 | 8/8 | 5/5 | 0.0078 |
| 50 | 0.415 | 0.512 | 0.810 | -0.097 | 4/8 | 3/5 | 0.3125 |
| 100 | 2.278 | 2.499 | 0.912 | -0.221 | 4/8 | 3/5 | 0.3125 |

## Section means for the ratio call

| Section | Donor | High | Low | Δ | Ratio | n high | n low | Usable |
|---|---|---:|---:|---:|---:|---:|---:|---|
| LUAD-5 R1 | Lung5 | 0.0136 | 0.0193 | -0.0057 | 0.706 | 12611 | 5226 | yes |
| LUAD-5 R2 | Lung5 | 0.0134 | 0.0150 | -0.0016 | 0.895 | 12377 | 5537 | yes |
| LUAD-5 R3 | Lung5 | 0.0158 | 0.0184 | -0.0026 | 0.858 | 9985 | 5908 | yes |
| LUSC-6 | Lung6 | 0.0122 | 0.0112 | 0.0010 | 1.086 | 16323 | 49873 | yes |
| LUAD-9 R1 | Lung9 | 0.0380 | 0.0999 | -0.0619 | 0.380 | 27834 | 11030 | yes |
| LUAD-9 R2 | Lung9 | 0.0306 | 0.0909 | -0.0603 | 0.336 | 57922 | 36954 | yes |
| LUAD-12 | Lung12 | 0.0252 | 0.0823 | -0.0571 | 0.306 | 7260 | 11007 | yes |
| LUAD-13 | Lung13 | 0.1779 | 0.2332 | -0.0553 | 0.763 | 19177 | 6853 | yes |

Detected versus absent does not use a quantile, so a donor-wide threshold labels the same cells as the section threshold. The donor-unit row matches this table.

## Section means for the Δ call

| Section | Donor | High | Low | Δ | Ratio | n high | n low | Usable |
|---|---|---:|---:|---:|---:|---:|---:|---|
| LUAD-5 R1 | Lung5 | 0.0831 | 0.0833 | -0.0003 | 0.997 | 5358 | 5945 | yes |
| LUAD-5 R2 | Lung5 | 0.0917 | 0.1168 | -0.0251 | 0.785 | 5382 | 5851 | yes |
| LUAD-5 R3 | Lung5 | 0.1190 | 0.1582 | -0.0392 | 0.752 | 4478 | 5488 | yes |
| LUSC-6 | Lung6 | 0.1524 | 0.1029 | 0.0494 | 1.480 | 3895 | 8439 | yes |
| LUAD-9 R1 | Lung9 | 0.4617 | 0.6985 | -0.2368 | 0.661 | 11669 | 12108 | yes |
| LUAD-9 R2 | Lung9 | 0.2184 | 0.3516 | -0.1332 | 0.621 | 28384 | 36713 | yes |
| LUAD-12 | Lung12 | 0.1723 | 0.3449 | -0.1726 | 0.500 | 3815 | 6705 | yes |
| LUAD-13 | Lung13 | 1.5431 | 1.6387 | -0.0956 | 0.942 | 7822 | 8186 | yes |

## Eligible specs with the largest |log ratio|

| Cut | Cytotoxic | µm | Unit | High | Low | Ratio | Δ | Section / donor signs | FOV high<low | Direction |
|---|---|---:|---|---:|---:|---:|---:|---|---|---|
| detected_absent | cd8nk | 20 | section | 0.041 | 0.071 | 0.573 | -0.030 | 7/8 neg, 1/8 pos; donors 4/5 neg, 1/5 pos | 156/217 | exclusion |
| detected_absent | cd8nk | 20 | donor | 0.041 | 0.071 | 0.573 | -0.030 | 7/8 neg, 1/8 pos; donors 4/5 neg, 1/5 pos | 156/217 | exclusion |
| detected_absent | cd8 | 20 | section | 0.033 | 0.053 | 0.618 | -0.020 | 7/8 neg, 1/8 pos; donors 4/5 neg, 1/5 pos | 146/217 | exclusion |
| detected_absent | cd8 | 20 | donor | 0.033 | 0.053 | 0.618 | -0.020 | 7/8 neg, 1/8 pos; donors 4/5 neg, 1/5 pos | 146/217 | exclusion |
| detected_absent | cd8nk | 25 | section | 0.073 | 0.117 | 0.623 | -0.044 | 7/8 neg, 1/8 pos; donors 4/5 neg, 1/5 pos | 162/217 | exclusion |
| detected_absent | cd8nk | 25 | donor | 0.073 | 0.117 | 0.623 | -0.044 | 7/8 neg, 1/8 pos; donors 4/5 neg, 1/5 pos | 162/217 | exclusion |
| q4_q1 | cd8nk | 20 | fov | 0.052 | 0.082 | 0.634 | -0.030 | 8/8 neg, 0/8 pos; donors 5/5 neg, 0/5 pos | 146/215 | exclusion |
| top30_vs_bottom | cd8 | 25 | section | 0.062 | 0.098 | 0.635 | -0.036 | 7/8 neg, 0/8 pos; donors 4/5 neg, 0/5 pos | 139/188 | exclusion |
| top30_vs_bottom | cd8 | 25 | donor | 0.063 | 0.098 | 0.637 | -0.036 | 7/8 neg, 0/8 pos; donors 4/5 neg, 0/5 pos | 137/188 | exclusion |
| detected_absent | cd8 | 20 | fov | 0.037 | 0.059 | 0.637 | -0.021 | 7/8 neg, 1/8 pos; donors 5/5 neg, 0/5 pos | 146/217 | exclusion |
| top20_vs_bottom | cd8nk | 20 | fov | 0.051 | 0.081 | 0.637 | -0.029 | 7/8 neg, 1/8 pos; donors 4/5 neg, 1/5 pos | 134/202 | exclusion |
| q4_q1 | cd8 | 25 | section | 0.056 | 0.087 | 0.646 | -0.031 | 7/8 neg, 1/8 pos; donors 4/5 neg, 1/5 pos | 146/217 | exclusion |

## Inventory

| Section | Donor | Cells | Malignant | CLDN4>0 | CD8 | NK | GZMB+ non-tumor | Median NN µm |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| LUAD-5 R1 | Lung5 | 98002 | 17837 | 12611 | 2020 | 1668 | 4757 | 8.31 |
| LUAD-5 R2 | Lung5 | 100335 | 17914 | 12377 | 1848 | 1748 | 6666 | 8.16 |
| LUAD-5 R3 | Lung5 | 97809 | 15893 | 9985 | 2294 | 1614 | 4295 | 8.34 |
| LUSC-6 | Lung6 | 89975 | 66196 | 16323 | 1230 | 184 | 834 | 9.80 |
| LUAD-9 R1 | Lung9 | 87606 | 38864 | 27834 | 2323 | 746 | 2268 | 8.13 |
| LUAD-9 R2 | Lung9 | 139504 | 94876 | 57922 | 1679 | 1182 | 1737 | 8.69 |
| LUAD-12 | Lung12 | 71304 | 18267 | 7260 | 1310 | 567 | 4043 | 8.89 |
| LUAD-13 | Lung13 | 81236 | 26030 | 19177 | 3251 | 105 | 3311 | 7.76 |

Malignant cells in the neighbor table: 295877. GZMB+ non-tumor cells: 27911. Largest shares: fibroblast 15%, neutrophil 13%, macrophage 9%. CD8 is 7% and NK is 6% of that set, so GZMB+ is not a purified cytotoxic class. Full GZMB+ counts are in `summary.json`.

## What this does not claim

- It does not replace 0.36 / 0.52, and it does not restate that pair as a new estimate from this grid.
- It does not say nearby effectors are muzzled. This sweep counts cells, not GZMB/PRF1/NKG7/IFNG inside the cells that are present.
- Sign tests and Wilcoxon p-values are descriptive. FOVs are nested in five donors, and the grid was searched.
- No private 8-KL. No ICI labels.

```bash
python3 scripts/download_cosmx_nsclc_h5ad.py
python3 scripts/cosmx_cytotoxic_neighbor_sweep.py
```
