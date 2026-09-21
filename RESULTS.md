# CosMx core vs margin: is CLDN4-high exclusion stronger in the core?

Public CosMx NSCLC (He et al. 2022; CellCharter object, figshare 25976224; 765,771 cells, 8 sections, 5 patients). CLDN4-only. No private 8-KL.

**Answer.** CLDN4-high tumor cells have fewer NK+CD8 neighbors in the core as well as at the margin. The absolute deficit is not larger in the core than at the margin. The high/low count ratio is lower in the core in 7/8 sections, and that fold difference is significant against a FOV label-shuffle null. Infiltration-depth AUC is lower for CLDN4-high in 7/8 sections and is not claimed: the paired test is p = 0.078, and LUSC-6 runs the other way.

This layer does not revise the locked unstratified cytotoxic ratios (0.36 at 50 µm, 0.52 at 100 µm; 8/8 sections, 5/5 patients, sign P = 0.031).

## Definitions (set before the cross-section test)

- Tumor: author `cell_type` `tumor 5/6/9/12/13` (epithelial is not tumor).
- Cytotoxic neighbors: `NK`, `T CD8 memory`, `T CD8 naive`.
- CLDN4-high: log-normalized `X` above the within-section tumor median. LUSC-6 and LUAD-12 have median 0, so high means CLDN4 > 0 (24% and 39% of tumor cells).
- Immune-rich regions: CellCharter niches `immune`, `lymphoid structure`, `macrophages`, `myeloid-enriched stroma`, `neutrophils`, `plasmablast-enriched stroma`. Generic `stroma` is a sensitivity.
- Distance: µm to the nearest cell in an immune-rich niche (0 if the tumor cell sits in one). Pixel size 0.18 µm.
- Margin: distance ≤ 50 µm. Core: distance ≥ 150 µm.
- Primary metric: mean cytotoxic-neighbor count at 50 µm. Ratio = high / low (reported when the low-arm mean is at least 0.05).
- Exclusion: high − low < 0. Stronger in the core, on the absolute scale: (core delta − margin delta) < 0.
- Tests: section-paired Wilcoxon (n = 8) and patient sign test (n = 5, sections averaged within patient). Null: 999 shuffles of the CLDN4 label inside each FOV. The shuffle null is not centered at 0, because a random split already produces a larger absolute gap where counts are higher. Two-sided permutation p-values are recentered on the null median.

Every section had at least 40 tumor cells in each CLDN4 arm in both strata. Lung6 is far from immune-rich niches (median distance 677 µm; 1,392 margin tumor cells). Lung13 is mostly margin (median 19 µm; 1,283 core tumor cells).

## Primary result, 50 µm

| Stratum | Section-mean count, high | Section-mean count, low | Section-mean Δ (high−low) | Section-mean ratio | Sections with high < low | Wilcoxon p | Patients |
|---|---:|---:|---:|---:|---:|---:|---:|
| Core ≥150 µm | 0.268 | 0.369 | **−0.101** | **0.763** | **7/8** | **0.023** | 4/5 (sign p = 0.19) |
| Margin ≤50 µm | 1.073 | 1.248 | **−0.175** | **0.840** | **7/8** | 0.078 | 4/5 (sign p = 0.19) |

The core deficit is real on the section test (7/8, p = 0.023). It is not larger than the margin deficit. Mean (core Δ − margin Δ) = **+0.074**. The core gap is the more negative one in only **2/8** sections (LUSC-6, LUAD-13). Section Wilcoxon p = 0.55. Patients 2/5 (sign p = 0.81). Relative to the FOV-shuffle null (null median of the interaction +0.030), the margin-larger absolute gap has one-sided p = 0.022 and recentered two-sided p = 0.064.

The fold ratio tells a different scale. Mean (core ratio − margin ratio) = **−0.077**, with the core ratio smaller in **7/8** sections (the exception is LUAD-5 R3). Section Wilcoxon p = 0.15. The shuffle null expects the core ratio to sit slightly above the margin ratio (null median +0.020) because a ratio is noisy when counts are low. Against that null, the lower core ratio has one-sided p = 0.002 and recentered two-sided p = **0.003**.

So: exclusion is present in the core; the absolute neighbor gap is not stronger there; the fold reduction is.

### Per section, 50 µm mean cytotoxic neighbors

| Section | Core high | Core low | Core ratio | Margin high | Margin low | Margin ratio |
|---|---:|---:|---:|---:|---:|---:|
| LUAD-5 R1 | 0.100 | 0.146 | 0.69 | 0.502 | 0.564 | 0.89 |
| LUAD-5 R2 | 0.099 | 0.131 | 0.75 | 0.485 | 0.636 | 0.76 |
| LUAD-5 R3 | 0.116 | 0.124 | 0.94 | 0.577 | 0.780 | 0.74 |
| LUSC-6 | 0.109 | 0.100 | 1.10 | 1.378 | 1.123 | 1.23 |
| LUAD-9 R1 | 0.220 | 0.275 | 0.80 | 2.011 | 2.270 | 0.89 |
| LUAD-9 R2 | 0.225 | 0.357 | 0.63 | 1.373 | 1.991 | 0.69 |
| LUAD-12 | 0.180 | 0.453 | 0.40 | 0.385 | 0.712 | 0.54 |
| LUAD-13 | 1.090 | 1.366 | 0.80 | 1.870 | 1.905 | 0.98 |

LUSC-6 is the section where CLDN4-high has more, not fewer, cytotoxic neighbors in both strata. It is the only squamous sample.

Patient-averaged 50 µm ratios (core, margin): Lung5 0.79, 0.80; Lung6 1.10, 1.23; Lung9 0.71, 0.79; Lung12 0.40, 0.54; Lung13 0.80, 0.98. The core ratio is lower in 4/5 patients. Lung5 is flat. The absolute-gap interaction favors the core only in Lung6 and Lung13.

Contact (any cytotoxic neighbor within 50 µm): Mantel–Haenszel OR 0.73 in the core and 0.75 at the margin. The section log-OR interaction is null (p = 0.95, 4/8). The contact deficit is shared; it is not stronger in the core.

Density-normalized cytotoxic fraction moves the same way as the counts (core Δ −0.003, 7/8, p = 0.016; margin Δ −0.005, 7/8, p = 0.078). The absolute fraction gap is not larger in the core (3/8, p = 0.55).

## 100 µm

| Stratum | Section-mean Δ | Section-mean ratio | Sections high < low | Wilcoxon p | Patients |
|---|---:|---:|---:|---:|---:|
| Core | −0.240 | 0.875 | 6/8 | 0.039 | 4/5 |
| Margin | −0.202 | 0.940 | 4/8 | 0.64 | 3/5 |

Raw interaction of absolute deltas: 4/8, Wilcoxon p = 0.84. The shuffle null, though, expects a much larger margin gap (null median of core Δ − margin Δ = +0.149) because 100 µm neighborhoods at the margin contain more cytotoxic cells. The observed interaction (−0.038) is on the core side of that null (one-sided p = 0.001, recentered two-sided p = 0.003). Core ratios are again lower in 7/8 sections (Wilcoxon p = 0.15; shuffle recentered p = 0.001).

## Sensitivity: within-section distance tertiles

Far tertile versus near tertile, 50 µm counts. This keeps every section balanced, including Lung6 and Lung13.

- Far: Δ = −0.100, ratio 0.707, **8/8** sections, Wilcoxon p = 0.0078, **5/5** patients, sign p = 0.031.
- Near: Δ = −0.167, ratio 0.808, 6/8, p = 0.078, 3/5 patients.
- Absolute interaction still favors the near tertile (3/8 far-stronger, p = 0.25).

The far tertile is where the exclusion sign is uniform across patients. The absolute gap remains larger in the near tertile.

Author niche labels, without a distance cut: tumor cells in `tumor interior` versus `tumor-stroma boundary`. Interior ratio 0.83 versus boundary 0.94. The core ratio is smaller in 7/8 sections (Wilcoxon p = 0.023). The absolute-delta interaction is not (5/8, p = 0.55). Adding generic stroma to the immune-rich set does not create an absolute core-stronger gap (2/8, p = 0.25).

## Infiltration-depth AUC — not claimed

AUC is the trapezoid of the section’s mean 50 µm cytotoxic-neighbor count in 20 µm bins of distance to an immune-rich niche, from 0 to 200 µm.

Section-mean AUC: CLDN4-high **90.2**, CLDN4-low **108.1**. High < low in **7/8** sections and 4/5 patients. Paired Wilcoxon p = **0.078**. LUSC-6 is the opposite section (AUC 115 versus 86) and is large enough to keep the rank test above 0.05. The section-mean curve itself is lower for CLDN4-high at every bin from 10 to 190 µm.

That direction matches exclusion. It is not significant on the pre-specified section test, and one section opposes it, so infiltration-depth AUC is not claimed.

## Unstratified reference

Same cells, no core/margin split. 50 µm section-mean ratio **0.80**, Δ −0.126, 7/8 sections, Wilcoxon p = 0.023, patients 4/5. At 100 µm the unstratified contrast is 4/8 (p = 0.31). These are not the locked 0.36 / 0.52 ratios and do not replace them.

## Figures

- `results/cosmx_core_margin/figures/paired_ratio_core_margin.png` — section-paired high/low ratios.
- `results/cosmx_core_margin/figures/section_delta_50um.png` — absolute 50 µm gaps.
- `results/cosmx_core_margin/figures/depth_curve_50um.png` — depth curve and section AUCs.
- `results/cosmx_core_margin/figures/map_core_margin.png` — LUAD-9 R2 tumor cells, margin in orange, core in purple.

Tables: `results/cosmx_core_margin/tables/`. Full tests: `results/cosmx_core_margin/stats.json`.

```bash
python3 scripts/cosmx_core_margin_exclusion.py
```

The h5ad is figshare file 46841842 and is not committed.
