# CosMx NSCLC: smallest TACSTD2-high immune and CD8+NK ratio at 8/8 and 5/5

He et al. 2022 CosMx 960-plex (figshare 25976224, `cosmx_human_nsclc_clustered.h5ad`; 295,877 patient-matched tumor cells, 8 sections, 5 patients). This note searches the TACSTD2 high-versus-low neighbor contrast. It does **not** replace the locked CLDN4 cytotoxic ratios 0.36 at 50 µm and 0.52 at 100 µm, and it does **not** say that effector cells which are present have lower GZMB, PRF1, NKG7, or IFNG. No private 8-KL.

## What was held fixed

Index cells are the author tumor label matched to the section (`tumor 5/6/9/12/13`). Centroids are global pixels × 0.18 µm. The index cell is excluded. Immune is the same 14 labels as the previous TACSTD2 note (B, plasmablast, CD4 naive/memory, CD8 naive/memory, Treg, NK, mDC, pDC, monocyte, macrophage, mast, neutrophil). CD8+NK is T CD8 memory, T CD8 naive, and NK. CD8 alone is scored and is not allowed to replace CD8+NK. The ratio is the unweighted mean of the eight section high means divided by the unweighted mean of the eight section low means. A patient value is the unweighted mean of that patient's sections.

Radii were 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 28, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100 µm. Cuts were detected, ge2, ge3, ge4, ge5, ge6, ge8, ge10, ge12, ge15, ge20, pos_q50, pos_q60, pos_q70, pos_q75, pos_q80, pos_q90, pos_q95, log_median, log_q4q1, log_q90_vs_q10, log_q80_vs_q20, log_q75_vs_q25, log_q70_vs_q30, log_q90_vs_rest, log_q80_vs_rest, log_q75_vs_rest, log_q70_vs_rest, log_pos_q75, log_pos_q90, log_pos_q95, rank10, rank20, rank30, fov_median. `detected` and `geK` compare raw TACSTD2 count > 0 or ≥ K with count = 0. `pos_q*` is the numpy quantile of positive raw counts versus count = 0. `log_*` uses log1p(count / n_counts × 10,000) inside the section. `rank10/20/30` takes that exact fraction of cells from each tail; ties are broken by a fixed per-section jitter (seed from the section name and the cut), not by file order. `fov_median` is the within-FOV log median. Estimands: `frac` (neighbor fraction, cells with no neighbor dropped), `count` (neighbor count, empty balls kept as zero), `frac0` (empty balls kept as fraction 0). A fraction spec needs every low-arm section mean ≥ 0.005 and every high-arm section mean > 0. A count spec needs every low-arm section mean ≥ 0.001, every low arm summing to ≥ 20 events, and every high arm summing to ≥ 5 events. `frac0` is eligible only when `frac` for the same masks is also 8/8 and 5/5 and clears the fraction floor. 9486 specs had both arms at n ≥ 30 in every section. 1149 of those were eligible; 1121 of the eligible specs are immune or CD8+NK.

## Calibration

These equal-weight means match the previous TACSTD2 note within 5×10⁻⁴. Ranking started only after that check.

| Cut | Class | Estimand | µm | High | Low | Ratio |
|---|---|---|---:|---:|---:|---:|
| log_median | immune | frac | 10 | 0.0390 | 0.0751 | 0.519 |
| log_median | immune | frac | 20 | 0.0614 | 0.0955 | 0.643 |
| log_median | cd8nk | count | 10 | 0.0043 | 0.0095 | 0.455 |
| log_median | cd8nk | count | 20 | 0.0394 | 0.0588 | 0.670 |
| detected | immune | count | 10 | 0.0349 | 0.0886 | 0.393 |
| log_median | immune | count | 50 | 5.3003 | 6.2335 | 0.850 |
| detected | immune | frac | 10 | 0.0382 | 0.0780 | 0.490 |
| log_q4q1 | immune | frac | 10 | 0.0431 | 0.0780 | 0.553 |

Where the log median is 0 (Lung6, Lung9, Lung12, Lung13), `log_median` is the detected-versus-absent split.

## Calls

Section and patient sign tests sit on the 8/8 and 5/5 floors (sign p = 0.0039 and 0.031; Wilcoxon p = 0.0078 and 0.0625) whenever every unit has the same sign, so those p-values do not choose the spec. FOV p-values are nominal: FOVs sit inside five patients, and the grid was searched.

| Call | Spec | High | Low | Ratio | Δ | Weakest section ratio | Pooled ratio | FOVs lower |
|---|---|---:|---:|---:|---:|---:|---:|---|
| Smallest immune ratio | `ge5`, immune, frac0, 9 µm | 0.0070 | 0.0399 | **0.176** | **-0.0329** | 0.3721 | 0.143 | 150/156 |
| Largest immune \|Δ\| | `fov_median`, immune, count, 60 µm | 8.4695 | 9.7959 | **0.865** | **-1.3264** | 0.9996 | 0.863 | 156/216 |
| Smallest CD8+NK ratio | `detected`, cd8nk, count, 11 µm | 0.0066 | 0.0133 | **0.496** | **-0.0067** | 0.8015 | 0.401 | 152/214 |
| Largest CD8+NK \|Δ\| | `log_median`, cd8nk, count, 40 µm | 0.2392 | 0.2937 | **0.814** | **-0.0546** | 0.9940 | 0.717 | 157/214 |

Per estimand, still immune or CD8+NK, still eligible:

| Class | Estimand | Smallest ratio | That Δ | Largest \|Δ\| | That ratio |
|---|---|---|---:|---|---:|
| immune | frac | `ge4`, immune, frac, 8 µm = 0.330 | -0.0485 | `ge4`, immune, frac, 8 µm = -0.0485 | 0.330 |
| immune | count | `ge3`, immune, count, 8 µm = 0.181 | -0.0274 | `fov_median`, immune, count, 60 µm = -1.3264 | 0.865 |
| immune | frac0 | `ge5`, immune, frac0, 9 µm = 0.176 | -0.0329 | `pos_q80`, immune, frac0, 15 µm = -0.0432 | 0.487 |
| cd8nk | frac | none eligible |  | none eligible |  |
| cd8nk | count | `detected`, cd8nk, count, 11 µm = 0.496 | -0.0067 | `log_median`, cd8nk, count, 40 µm = -0.0546 | 0.814 |
| cd8nk | frac0 | none eligible |  | none eligible |  |

CD8 alone, smallest eligible ratio 0.711 (`pos_q60`, cd8, count, 20 µm), largest |Δ| -0.0160 (`pos_q60`, cd8, count, 23 µm). Reported so a CD8-only cut cannot be mistaken for the CD8+NK call.

## How the calls sit next to the published contrasts

The immune ratio call is the zero-filled fraction: empty balls count as 0. On the same masks, the fraction among cells that have a neighbor is **0.380** (high 0.0287, low 0.0754, Δ -0.0467, eligible). Quietest contact-fraction high arm: n = 44, section mean 0.0054.

The neighbor-count ratio on those same masks is 0.149 (Δ -0.0484; quietest high arm 2.0 immune neighbors). That count spec does not clear the event floor, so it is not a call.

The published CD8+NK count ratio at 10 µm, log median, is 0.455 (Δ -0.0052), 8/8 and 5/5. It is not eligible: the quietest low arm has 17 events and the quietest high arm has 5, against floors of 20 and 5. The eligible CD8+NK ratio is larger than that published ratio. No cut in this grid is a smaller CD8+NK ratio that also clears the event floor.

Immune absolute-drop call, weakest section: LUAD-13 ratio 0.9996 (Δ -0.0091, high 23.4448, low 23.4540). The equal-weight Δ is not the drop in every section.
CD8+NK absolute-drop call, weakest section: LUAD-9 R1 ratio 0.9940 (Δ -0.0018, high 0.2961, low 0.2979). The equal-weight Δ is not the drop in every section.

## Smallest immune ratio

`ge5`, immune, frac0, 9 µm. Equal-weight high 0.0070, low 0.0399, ratio 0.176, Δ -0.0329. Sections 8/8, patients 5/5. Weakest section ratio 0.3721. Cell-weighted ratio 0.143. Smallest high-arm section mean 0.0018; smallest low-arm section mean 0.0057. FOVs with at least 8 cells in each arm: 150/156 lower, median Δ -0.0214, nominal Wilcoxon p = 1.36e-24. Mean degree 0.43 (high) and 0.93 (low).

| Section | Patient | n high | n low | High | Low | Ratio | Δ | High events | Low events | Degree high | Degree low |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LUAD-5 R1 | Lung5 | 1654 | 6018 | 0.0057 | 0.0266 | 0.216 | -0.0209 | 9.5 | 160.1 | 0.68 | 1.24 |
| LUAD-5 R2 | Lung5 | 1155 | 6562 | 0.0030 | 0.0211 | 0.144 | -0.0181 | 3.5 | 138.3 | 0.77 | 1.44 |
| LUAD-5 R3 | Lung5 | 1117 | 6504 | 0.0048 | 0.0244 | 0.196 | -0.0196 | 5.3 | 158.7 | 0.65 | 1.14 |
| LUSC-6 | Lung6 | 1108 | 41541 | 0.0018 | 0.0057 | 0.317 | -0.0039 | 2.0 | 236.5 | 0.22 | 0.48 |
| LUAD-9 R1 | Lung9 | 458 | 28431 | 0.0044 | 0.0292 | 0.150 | -0.0248 | 2.0 | 829.9 | 0.36 | 0.87 |
| LUAD-9 R2 | Lung9 | 986 | 69427 | 0.0030 | 0.0314 | 0.097 | -0.0284 | 3.0 | 2182.3 | 0.32 | 0.85 |
| LUAD-12 | Lung12 | 381 | 12348 | 0.0079 | 0.1121 | 0.070 | -0.1042 | 3.0 | 1383.9 | 0.14 | 0.57 |
| LUAD-13 | Lung13 | 214 | 16977 | 0.0257 | 0.0691 | 0.372 | -0.0434 | 5.5 | 1172.5 | 0.29 | 0.82 |

## Largest immune absolute drop

`fov_median`, immune, count, 60 µm. Equal-weight high 8.4695, low 9.7959, ratio 0.865, Δ -1.3264. Sections 8/8, patients 5/5. Weakest section ratio 0.9996. Cell-weighted ratio 0.863. Smallest high-arm section mean 1.7545; smallest low-arm section mean 2.0494. FOVs with at least 8 cells in each arm: 156/216 lower, median Δ -0.9013, nominal Wilcoxon p = 3.60e-12. Mean degree 61.00 (high) and 62.35 (low).

| Section | Patient | n high | n low | High | Low | Ratio | Δ | High events | Low events | Degree high | Degree low |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LUAD-5 R1 | Lung5 | 8890 | 8946 | 7.4987 | 8.1178 | 0.924 | -0.6192 | 66663.0 | 72622.0 | 61.23 | 62.73 |
| LUAD-5 R2 | Lung5 | 8933 | 8981 | 8.4388 | 8.7197 | 0.968 | -0.2809 | 75384.0 | 78312.0 | 64.88 | 66.50 |
| LUAD-5 R3 | Lung5 | 7747 | 8146 | 7.5672 | 8.0323 | 0.942 | -0.4651 | 58623.0 | 65431.0 | 57.80 | 59.32 |
| LUSC-6 | Lung6 | 24233 | 41963 | 1.7545 | 2.0494 | 0.856 | -0.2948 | 42518.0 | 85997.0 | 59.01 | 58.31 |
| LUAD-9 R1 | Lung9 | 10433 | 28431 | 6.1859 | 7.0608 | 0.876 | -0.8750 | 64537.0 | 200746.0 | 72.81 | 74.09 |
| LUAD-9 R2 | Lung9 | 25449 | 69427 | 3.4896 | 5.8194 | 0.600 | -2.3298 | 88806.0 | 404022.0 | 61.66 | 63.85 |
| LUAD-12 | Lung12 | 5870 | 12397 | 9.3765 | 15.1137 | 0.620 | -5.7372 | 55040.0 | 187365.0 | 41.70 | 44.26 |
| LUAD-13 | Lung13 | 9053 | 16977 | 23.4448 | 23.4540 | 1.000 | -0.0091 | 212246.0 | 398178.0 | 68.89 | 69.77 |

## Smallest CD8+NK ratio

`detected`, cd8nk, count, 11 µm. Equal-weight high 0.0066, low 0.0133, ratio 0.496, Δ -0.0067. Sections 8/8, patients 5/5. Weakest section ratio 0.8015. Cell-weighted ratio 0.401. Smallest high-arm section mean 0.0011; smallest low-arm section mean 0.0021. FOVs with at least 8 cells in each arm: 152/214 lower, median Δ -0.0044, nominal Wilcoxon p = 5.13e-24. Mean degree 1.55 (high) and 1.87 (low).

| Section | Patient | n high | n low | High | Low | Ratio | Δ | High events | Low events | Degree high | Degree low |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LUAD-5 R1 | Lung5 | 11819 | 6018 | 0.0017 | 0.0043 | 0.392 | -0.0026 | 20.0 | 26.0 | 2.04 | 2.33 |
| LUAD-5 R2 | Lung5 | 11352 | 6562 | 0.0011 | 0.0030 | 0.376 | -0.0019 | 13.0 | 20.0 | 2.22 | 2.59 |
| LUAD-5 R3 | Lung5 | 9389 | 6504 | 0.0011 | 0.0034 | 0.315 | -0.0023 | 10.0 | 22.0 | 1.94 | 2.19 |
| LUSC-6 | Lung6 | 24655 | 41541 | 0.0014 | 0.0021 | 0.666 | -0.0007 | 34.0 | 86.0 | 1.05 | 1.25 |
| LUAD-9 R1 | Lung9 | 10433 | 28431 | 0.0105 | 0.0132 | 0.802 | -0.0026 | 110.0 | 374.0 | 1.54 | 1.91 |
| LUAD-9 R2 | Lung9 | 25449 | 69427 | 0.0059 | 0.0158 | 0.370 | -0.0100 | 149.0 | 1099.0 | 1.45 | 1.80 |
| LUAD-12 | Lung12 | 5919 | 12348 | 0.0032 | 0.0194 | 0.166 | -0.0161 | 19.0 | 239.0 | 0.85 | 1.19 |
| LUAD-13 | Lung13 | 9053 | 16977 | 0.0277 | 0.0449 | 0.617 | -0.0172 | 251.0 | 763.0 | 1.35 | 1.68 |

## Largest CD8+NK absolute drop

`log_median`, cd8nk, count, 40 µm. Equal-weight high 0.2392, low 0.2937, ratio 0.814, Δ -0.0546. Sections 8/8, patients 5/5. Weakest section ratio 0.9940. Cell-weighted ratio 0.717. Smallest high-arm section mean 0.0733; smallest low-arm section mean 0.0755. FOVs with at least 8 cells in each arm: 157/214 lower, median Δ -0.0418, nominal Wilcoxon p = 9.88e-14. Mean degree 27.71 (high) and 28.63 (low).

| Section | Patient | n high | n low | High | Low | Ratio | Δ | High events | Low events | Degree high | Degree low |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LUAD-5 R1 | Lung5 | 8918 | 8919 | 0.0983 | 0.1182 | 0.832 | -0.0198 | 877.0 | 1054.0 | 28.82 | 29.67 |
| LUAD-5 R2 | Lung5 | 8957 | 8957 | 0.0979 | 0.1159 | 0.845 | -0.0180 | 877.0 | 1038.0 | 30.76 | 31.73 |
| LUAD-5 R3 | Lung5 | 7942 | 7951 | 0.1150 | 0.1187 | 0.968 | -0.0038 | 913.0 | 944.0 | 27.59 | 28.15 |
| LUSC-6 | Lung6 | 24655 | 41541 | 0.0733 | 0.0755 | 0.971 | -0.0022 | 1808.0 | 3136.0 | 26.33 | 26.22 |
| LUAD-9 R1 | Lung9 | 10433 | 28431 | 0.2961 | 0.2979 | 0.994 | -0.0018 | 3089.0 | 8469.0 | 32.30 | 33.44 |
| LUAD-9 R2 | Lung9 | 25449 | 69427 | 0.1509 | 0.2853 | 0.529 | -0.1344 | 3840.0 | 19807.0 | 27.97 | 29.40 |
| LUAD-12 | Lung12 | 5919 | 12348 | 0.1387 | 0.3400 | 0.408 | -0.2013 | 821.0 | 4198.0 | 18.27 | 19.99 |
| LUAD-13 | Lung13 | 9053 | 16977 | 0.9432 | 0.9985 | 0.945 | -0.0553 | 8539.0 | 16952.0 | 29.62 | 30.46 |

## Concordant specs that missed the floor

The smallest immune or CD8+NK ratio with a strict 8/8 and 5/5 sign, ignoring the floor and the frac0 contact rule, is **0.000** (`ge4`, immune, frac0, 5 µm; high 0.0000, low 0.0017, min low 0.0001, min high events 0.0, min low events 4.0). It is not a call. The floor is what keeps a near-empty arm from manufacturing a ratio.

The five smallest such non-calls:

| Cut | Class | Estimand | µm | Ratio | Δ | Min low | Min high events | Min low events |
|---|---|---|---:|---:|---:|---:|---:|---:|
| ge4 | immune | frac0 | 5 | 0.000 | -0.0017 | 0.0001 | 0.0 | 4.0 |
| ge4 | immune | count | 5 | 0.000 | -0.0018 | 0.0001 | 0.0 | 4.0 |
| ge5 | cd8nk | count | 7 | 0.000 | -0.0021 | 0.0002 | 0.0 | 2.0 |
| ge5 | cd8nk | frac0 | 7 | 0.000 | -0.0019 | 0.0001 | 0.0 | 0.8 |
| ge5 | immune | count | 5 | 0.000 | -0.0018 | 0.0001 | 0.0 | 4.0 |

## Ten smallest eligible immune or CD8+NK ratios

| Cut | Class | Estimand | µm | Ratio | Δ | Max section ratio | Pooled ratio | Min n high |
|---|---|---|---:|---:|---:|---:|---:|---:|
| ge5 | immune | frac0 | 9 | 0.176 | -0.0329 | 0.372 | 0.143 | 214 |
| ge3 | immune | count | 8 | 0.181 | -0.0274 | 0.423 | 0.155 | 1119 |
| pos_q80 | immune | count | 8 | 0.204 | -0.0267 | 0.380 | 0.185 | 1346 |
| ge3 | immune | count | 9 | 0.231 | -0.0437 | 0.472 | 0.195 | 1119 |
| pos_q95 | immune | count | 10 | 0.235 | -0.0678 | 0.463 | 0.181 | 381 |
| pos_q70 | immune | count | 8 | 0.235 | -0.0256 | 0.380 | 0.209 | 2727 |
| pos_q75 | immune | count | 8 | 0.235 | -0.0256 | 0.380 | 0.209 | 2727 |
| ge4 | immune | count | 10 | 0.235 | -0.0678 | 0.463 | 0.201 | 456 |
| ge2 | immune | count | 8 | 0.239 | -0.0255 | 0.380 | 0.204 | 2727 |
| pos_q90 | immune | count | 10 | 0.239 | -0.0674 | 0.516 | 0.207 | 716 |

## Ten largest eligible absolute drops, immune or CD8+NK

| Cut | Class | Estimand | µm | Δ | Ratio | Max section ratio |
|---|---|---|---:|---:|---:|---:|
| fov_median | immune | count | 60 | -1.3264 | 0.865 | 1.000 |
| detected | immune | count | 60 | -1.3264 | 0.865 | 1.000 |
| detected | immune | count | 55 | -1.2139 | 0.848 | 0.994 |
| fov_median | immune | count | 55 | -1.1914 | 0.850 | 0.994 |
| detected | immune | count | 50 | -1.0848 | 0.829 | 0.987 |
| fov_median | immune | count | 50 | -1.0497 | 0.833 | 0.987 |
| ge2 | immune | count | 45 | -1.0194 | 0.792 | 0.992 |
| pos_q60 | immune | count | 45 | -0.9858 | 0.799 | 0.992 |
| detected | immune | count | 45 | -0.9453 | 0.807 | 0.980 |
| log_median | immune | count | 50 | -0.9332 | 0.850 | 1.000 |

Count Δ grows with the area of the ball. A larger count Δ at a long radius is a larger absolute drop, not automatically a stronger fold. The fold call is the ratio table.

## What this does not claim

- Not a re-estimate of the locked CLDN4 50/100 µm cytotoxic ratios 0.36 / 0.52.
- Not muzzling. Neighbor counts are not GZMB, PRF1, NKG7, or IFNG levels inside the effector cells that remain.
- Not a smaller section-level p-value than the 8/8 Wilcoxon floor.
- FOV p-values are nominal.
- CD8 alone is not the CD8+NK result.
- No ICI labels. No private 8-KL. No Visium same-spot correlation written as exclusion.

```bash
python3 scripts/download_cosmx_nsclc_h5ad.py
python3 scripts/cosmx_tacstd2_max_effect.py
```

