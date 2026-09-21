# CosMx CLDN4 sweep — effect sizes for immune-cold and barrier ligands

This file is a parameter sweep. The pre-specified v1 analysis (tertiles, 40 µm proximity, SpatialDM `l=30` µm, COMMOT `dis_thr=50` µm) is unchanged and is still the primary CCC report in `RESULTS.md`. Nothing here replaces the locked contact odds ratio or the Ripley g(r) result.

Effect sizes below are computed from the same QC cells, epithelial tumor gate, and CD8 gate as v1. The CLDN4 split, the radius, and the FOV filter change. A negative CD8-fraction delta means CLDN4-high tumor cells less often have a CD8 cell inside the radius (immune-cold). A ligand fold above 1 means higher expression on the CLDN4-high arm.

Wilcoxon tests are on the 8 slide means. p-values on a cell that was chosen from this grid are descriptive of the sweep.

## Selection rules

- Eligible: 8 slides and at least 80 FOVs.
- Immune-cold: most negative median slide Δ in the fraction of tumor cells with a CD8 neighbor within the radius, with the high arm lower on at least 6 slides.
- Barrier ligands: largest mean log2 fold of interface expression for CDH1, ICAM1, and MIF together, each fold > 1, each high arm higher on at least 6 slides, and both medians at least 0.05 so a tiny denominator cannot win.
- CDH1 is also one of the epithelial markers used to call tumor. ICAM1 and MIF are not.

## Strongest immune-cold panel

CLDN4-high tumor cells are less often within the radius of a CD8 cell than CLDN4-low tumor cells in this setting (median slide difference -9.5 percentage points).

- Setting: split=tertile, radius=40 µm, FOV filter ≥10/10/5 (high/low/CD8), CLDN4 gap ≥0.50
- FOVs used: 184
- Levels: high 0.306, low 0.434, median slide Δ -0.095, fold 0.70
- Slides: 0/8 Δ>0, 8/8 Δ<0, Wilcoxon p=0.008; tissues 0/5 Δ>0 and 5/5 Δ<0

## Strongest barrier-ligand panel

Interface means are log1p median-library expression on tumor cells that have a CD8 cell within the radius. CDH1 high 0.152 vs low 0.089 (fold 1.70, Δ 0.050, 7/8 slides Δ>0, p=0.023, tissues 5/5) ICAM1 high 0.163 vs low 0.089 (fold 1.84, Δ 0.049, 7/8 slides Δ>0, p=0.023, tissues 5/5) MIF high 0.440 vs low 0.319 (fold 1.38, Δ 0.095, 8/8 slides Δ>0, p=0.008, tissues 5/5)

- Setting: split=tertile, radius=15 µm, FOV filter ≥40/40/15 (high/low/CD8), CLDN4 gap ≥0.00
- FOVs used: 187
- Levels: high 0.152, low 0.089, median slide Δ 0.050, fold 1.70
- Slides: 7/8 Δ>0, 1/8 Δ<0, Wilcoxon p=0.023; tissues 5/5 Δ>0 and 0/5 Δ<0

Cell-level Cohen's d on the same interface cells (median across slides of the per-FOV d) is a standardized effect, not a fold. CDH1 d=0.16 (7/8, p=0.023); ICAM1 d=0.13 (8/8, p=0.008); MIF d=0.24 (8/8, p=0.008). These d values are modest. The fold is a shift in the mean, and the two cell distributions still overlap.

Detection rate on those same interface cells: CDH1 23.9% vs 11.6% (+11.0 pp, 8/8, p=0.008, tissues 5/5); ICAM1 22.2% vs 10.8% (+9.6 pp, 8/8, p=0.008, tissues 5/5); MIF 48.4% vs 33.8% (+16.9 pp, 8/8, p=0.008, tissues 5/5).

## Immune-compartment neighborhood

The CD8-fraction rule above is the neighborhood of the CCC receiver. The broader immune-compartment call (the same epithelial / immune / stromal argmax used to define tumor) is a larger immune-cold contrast. It is reported here and is not a ligand–receptor score.

Strongest eligible row: split=tertile, radius=20 µm, FOV filter ≥10/10/5 (high/low/CD8), CLDN4 gap ≥0.50. High 0.272 vs low 0.492 (-23.5 percentage points, fold 0.55), 8/8 slides Δ<0, Wilcoxon p=0.008, tissues 5/5 Δ<0, 184 FOVs.

Without a CLDN4-gap filter and with the loosest cell floor (≥10/10/5), the strongest radius/split is quartile at 20 µm: high 0.317 vs low 0.515 (-22.0 percentage points), 8/8 slides, p=0.008, tissues 5/5, 225 FOVs.

The same tertile / 40 µm CD8 contrast with no gap filter is high 0.368 vs low 0.456 (-9.1 percentage points), 8/8, p=0.008, tissues 5/5, 225 FOVs. The gap filter changes this by a fraction of a percentage point.

## Same settings, the other ligands

These rows use the barrier-winning setting when one exists, so ligands that do not support the thesis stay visible. Interface log-mean, CLDN4-high minus CLDN4-low.

- CCL5: high 0.051 vs low 0.073, fold 0.70, Δ -0.016, 2/8 Δ>0, p=0.250, tissues 2/5 Δ>0
- CD274: high 0.058 vs low 0.061, fold 0.95, Δ -0.008, 2/8 Δ>0, p=0.312, tissues 2/5 Δ>0
- CDH1: high 0.152 vs low 0.089, fold 1.70, Δ 0.050, 7/8 Δ>0, p=0.023, tissues 5/5 Δ>0
- CXCL16: high 0.157 vs low 0.113, fold 1.39, Δ 0.049, 7/8 Δ>0, p=0.016, tissues 5/5 Δ>0
- CXCL9: high 0.063 vs low 0.057, fold 1.11, Δ 0.003, 4/8 Δ>0, p=0.844, tissues 2/5 Δ>0
- HLA-A: high 0.930 vs low 0.824, fold 1.13, Δ 0.089, 7/8 Δ>0, p=0.023, tissues 5/5 Δ>0
- ICAM1: high 0.163 vs low 0.089, fold 1.84, Δ 0.049, 7/8 Δ>0, p=0.023, tissues 5/5 Δ>0
- IFNG_on_cd8: high 0.046 vs low 0.043, fold 1.08, Δ 0.003, 6/8 Δ>0, p=0.383, tissues 3/5 Δ>0
- LGALS9: high 0.093 vs low 0.084, fold 1.10, Δ 0.004, 6/8 Δ>0, p=0.383, tissues 3/5 Δ>0
- MIF: high 0.440 vs low 0.319, fold 1.38, Δ 0.095, 8/8 Δ>0, p=0.008, tissues 5/5 Δ>0
- PDCD1LG2: high 0.040 vs low 0.032, fold 1.23, Δ 0.010, 6/8 Δ>0, p=0.078, tissues 5/5 Δ>0
- TGFB1: high 0.066 vs low 0.085, fold 0.78, Δ -0.006, 3/8 Δ>0, p=0.547, tissues 2/5 Δ>0

IFNG_on_cd8 is IFNG on CD8 cells sitting within the radius of that tumor class. A flat contrast is the not-muzzling check.

## How much of the grid points the same way

Eligible CD8-fraction rows: 168. Median Δ < 0 (CLDN4-high has fewer nearby CD8): 168. Of those, Wilcoxon p≤0.05: 147.

Base filter only (≥10/10/5, no CLDN4-gap filter), median slide Δ in CD8 fraction (negative = immune-cold):

| split | 15 µm | 20 µm | 25 µm | 40 µm | 50 µm | 80 µm | 100 µm |
|---|---:|---:|---:|---:|---:|---:|---:|
| decile | -0.042 | -0.052 | -0.055 | -0.065 | -0.058 | -0.045 | -0.025 |
| quartile | -0.043 | -0.059 | -0.074 | -0.091 | -0.087 | -0.050 | -0.030 |
| quintile | -0.042 | -0.056 | -0.074 | -0.088 | -0.085 | -0.053 | -0.038 |
| tertile | -0.038 | -0.055 | -0.069 | -0.091 | -0.082 | -0.043 | -0.027 |

Full grid: `sweep_table.tsv`. One row per split, radius, FOV filter, endpoint, and ligand.

## CCC tools at the selected settings

### immune_cold_tertile_40um_g50 (tertile, 40 µm)

- CD274|PDCD1 commot: high 0.002 vs low 0.003, fold 0.68, Δ -5.34e-04, 2/8 slides Δ>0, p=0.250, tissues 2/5 Δ>0
- CD274|PDCD1 commot_frac_pos: high 0.011 vs low 0.009, fold 1.23, Δ 2.64e-04, 4/8 slides Δ>0, p=0.383, tissues 3/5 Δ>0
- CD274|PDCD1 commot_supply: high 0.049 vs low 0.046, fold 1.07, Δ -0.004, 3/8 slides Δ>0, p=0.844, tissues 2/5 Δ>0
- CD274|PDCD1 expr: high 0.043 vs low 0.058, fold 0.74, Δ -0.004, 4/8 slides Δ>0, p=0.461, tissues 3/5 Δ>0
- CD274|PDCD1 sdm_prox: high 0.054 vs low 0.022, fold 2.41, Δ 0.051, 5/8 slides Δ>0, p=0.148, tissues 4/5 Δ>0
- CD274|PDCD1 squidpy: high 0.058 vs low 0.058, fold 0.99, Δ -0.004, 2/8 slides Δ>0, p=0.148, tissues 2/5 Δ>0
- CDH1|CDH1 commot: high 1.36e-04 vs low 1.63e-04, fold 0.83, Δ 4.35e-05, 5/8 slides Δ>0, p=0.547, tissues 3/5 Δ>0
- CDH1|CDH1 commot_frac_pos: high 0.041 vs low 0.021, fold 1.98, Δ 0.020, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- CDH1|CDH1 commot_supply: high 7.85e-04 vs low 7.80e-04, fold 1.01, Δ -9.82e-05, 3/8 slides Δ>0, p=0.844, tissues 2/5 Δ>0
- CDH1|CDH1 expr: high 0.258 vs low 0.136, fold 1.90, Δ 0.103, 7/8 slides Δ>0, p=0.016, tissues 5/5 Δ>0
- CDH1|CDH1 sdm_prox: high 0.061 vs low 0.040, fold 1.54, Δ 0.013, 6/8 slides Δ>0, p=0.078, tissues 4/5 Δ>0
- CDH1|CDH1 squidpy: high 0.149 vs low 0.095, fold 1.56, Δ 0.050, 7/8 slides Δ>0, p=0.016, tissues 5/5 Δ>0
- CDH1|ITGA2_ITGB1 commot: high 0.002 vs low 0.001, fold 1.98, Δ 8.62e-04, 7/8 slides Δ>0, p=0.016, tissues 5/5 Δ>0
- CDH1|ITGA2_ITGB1 commot_frac_pos: high 0.034 vs low 0.012, fold 2.72, Δ 0.020, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- CDH1|ITGA2_ITGB1 commot_supply: high 0.009 vs low 0.010, fold 0.85, Δ 0.002, 5/8 slides Δ>0, p=0.461, tissues 4/5 Δ>0
- CDH1|ITGA2_ITGB1 expr: high 0.258 vs low 0.136, fold 1.90, Δ 0.103, 7/8 slides Δ>0, p=0.016, tissues 5/5 Δ>0
- CDH1|ITGA2_ITGB1 sdm_prox: high 0.074 vs low 0.040, fold 1.82, Δ 0.044, 6/8 slides Δ>0, p=0.055, tissues 5/5 Δ>0
- CDH1|ITGA2_ITGB1 squidpy: high 0.141 vs low 0.081, fold 1.75, Δ 0.050, 7/8 slides Δ>0, p=0.016, tissues 5/5 Δ>0
- HLA-A|CD8A commot: high 0.016 vs low 0.021, fold 0.73, Δ -0.005, 0/8 slides Δ>0, p=0.008, tissues 0/5 Δ>0
- HLA-A|CD8A commot_frac_pos: high 0.557 vs low 0.440, fold 1.27, Δ 0.101, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- HLA-A|CD8A commot_supply: high 0.020 vs low 0.031, fold 0.67, Δ -0.010, 0/8 slides Δ>0, p=0.008, tissues 0/5 Δ>0
- HLA-A|CD8A expr: high 0.817 vs low 0.759, fold 1.08, Δ 0.093, 7/8 slides Δ>0, p=0.055, tissues 4/5 Δ>0
- HLA-A|CD8A sdm_prox: high 0.436 vs low 0.404, fold 1.08, Δ 0.084, 7/8 slides Δ>0, p=0.039, tissues 4/5 Δ>0
- HLA-A|CD8A squidpy: high 0.696 vs low 0.650, fold 1.07, Δ 0.046, 7/8 slides Δ>0, p=0.078, tissues 4/5 Δ>0
- ICAM1|ITGAL commot: high 0.013 vs low 0.007, fold 1.84, Δ 0.005, 7/8 slides Δ>0, p=0.016, tissues 4/5 Δ>0
- ICAM1|ITGAL commot_frac_pos: high 0.044 vs low 0.027, fold 1.65, Δ 0.017, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- ICAM1|ITGAL commot_supply: high 0.076 vs low 0.068, fold 1.11, Δ 0.004, 6/8 slides Δ>0, p=0.195, tissues 3/5 Δ>0
- ICAM1|ITGAL expr: high 0.131 vs low 0.098, fold 1.34, Δ 0.033, 7/8 slides Δ>0, p=0.016, tissues 5/5 Δ>0
- ICAM1|ITGAL sdm_prox: high 0.092 vs low 0.040, fold 2.30, Δ 0.063, 7/8 slides Δ>0, p=0.195, tissues 4/5 Δ>0
- ICAM1|ITGAL squidpy: high 0.143 vs low 0.111, fold 1.28, Δ 0.022, 7/8 slides Δ>0, p=0.016, tissues 5/5 Δ>0
- ICAM1|ITGAL_ITGB2 commot: high 0.004 vs low 0.002, fold 1.78, Δ 0.002, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- ICAM1|ITGAL_ITGB2 commot_frac_pos: high 0.017 vs low 0.009, fold 1.99, Δ 0.008, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- ICAM1|ITGAL_ITGB2 commot_supply: high 0.024 vs low 0.019, fold 1.26, Δ 0.005, 6/8 slides Δ>0, p=0.195, tissues 4/5 Δ>0
- ICAM1|ITGAL_ITGB2 expr: high 0.131 vs low 0.098, fold 1.34, Δ 0.033, 7/8 slides Δ>0, p=0.016, tissues 5/5 Δ>0
- ICAM1|ITGAL_ITGB2 sdm_prox: high 0.082 vs low 0.042, fold 1.96, Δ 0.037, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- ICAM1|ITGAL_ITGB2 squidpy: high 0.099 vs low 0.071, fold 1.40, Δ 0.022, 7/8 slides Δ>0, p=0.016, tissues 5/5 Δ>0
- IFNG|IFNGR1_IFNGR2 commot: high 0.002 vs low 0.001, fold 1.22, Δ -1.81e-04, 4/8 slides Δ>0, p=0.641, tissues 3/5 Δ>0
- IFNG|IFNGR1_IFNGR2 commot_frac_pos: high 0.008 vs low 0.005, fold 1.66, Δ 0.002, 5/8 slides Δ>0, p=0.148, tissues 4/5 Δ>0
- IFNG|IFNGR1_IFNGR2 commot_supply: high 0.041 vs low 0.045, fold 0.90, Δ -0.006, 4/8 slides Δ>0, p=0.641, tissues 3/5 Δ>0
- IFNG|IFNGR1_IFNGR2 expr: high 0.034 vs low 0.033, fold 1.01, Δ 9.26e-04, 6/8 slides Δ>0, p=0.250, tissues 4/5 Δ>0
- IFNG|IFNGR1_IFNGR2 sdm_prox: high 4.01e-04 vs low 0.004, fold 0.10, Δ 0.005, 4/8 slides Δ>0, p=0.742, tissues 3/5 Δ>0
- IFNG|IFNGR1_IFNGR2 squidpy: high 0.035 vs low 0.027, fold 1.30, Δ 0.006, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- LGALS9|HAVCR2 commot: high 0.003 vs low 0.003, fold 1.03, Δ 1.23e-04, 5/8 slides Δ>0, p=0.641, tissues 3/5 Δ>0
- LGALS9|HAVCR2 commot_frac_pos: high 0.024 vs low 0.017, fold 1.42, Δ 0.008, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- LGALS9|HAVCR2 commot_supply: high 0.033 vs low 0.037, fold 0.88, Δ -0.006, 3/8 slides Δ>0, p=0.844, tissues 2/5 Δ>0
- LGALS9|HAVCR2 expr: high 0.085 vs low 0.076, fold 1.12, Δ 0.012, 6/8 slides Δ>0, p=0.055, tissues 4/5 Δ>0
- LGALS9|HAVCR2 sdm_prox: high 0.043 vs low 0.024, fold 1.79, Δ -0.003, 4/8 slides Δ>0, p=0.742, tissues 3/5 Δ>0
- LGALS9|HAVCR2 squidpy: high 0.070 vs low 0.067, fold 1.05, Δ 0.006, 7/8 slides Δ>0, p=0.055, tissues 4/5 Δ>0
- MIF|CD74_CD44 commot: high 0.018 vs low 0.020, fold 0.89, Δ -5.40e-04, 3/8 slides Δ>0, p=0.547, tissues 2/5 Δ>0
- MIF|CD74_CD44 commot_frac_pos: high 0.184 vs low 0.137, fold 1.34, Δ 0.047, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- MIF|CD74_CD44 commot_supply: high 0.034 vs low 0.051, fold 0.67, Δ -0.012, 0/8 slides Δ>0, p=0.008, tissues 0/5 Δ>0
- MIF|CD74_CD44 expr: high 0.637 vs low 0.507, fold 1.25, Δ 0.148, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- MIF|CD74_CD44 sdm_prox: high 0.171 vs low 0.198, fold 0.87, Δ -0.006, 4/8 slides Δ>0, p=0.945, tissues 3/5 Δ>0
- MIF|CD74_CD44 squidpy: high 0.418 vs low 0.362, fold 1.15, Δ 0.072, 7/8 slides Δ>0, p=0.016, tissues 4/5 Δ>0
- PDCD1LG2|PDCD1 commot: high 0.001 vs low 0.002, fold 0.64, Δ -4.65e-04, 2/8 slides Δ>0, p=0.312, tissues 2/5 Δ>0
- PDCD1LG2|PDCD1 commot_frac_pos: high 0.007 vs low 0.006, fold 1.20, Δ 0.002, 5/8 slides Δ>0, p=0.109, tissues 4/5 Δ>0
- PDCD1LG2|PDCD1 commot_supply: high 0.057 vs low 0.067, fold 0.85, Δ 0.005, 4/8 slides Δ>0, p=0.945, tissues 3/5 Δ>0
- PDCD1LG2|PDCD1 expr: high 0.028 vs low 0.032, fold 0.87, Δ -0.005, 1/8 slides Δ>0, p=0.016, tissues 0/5 Δ>0
- PDCD1LG2|PDCD1 sdm_prox: high 0.017 vs low 0.050, fold 0.34, Δ -0.007, 4/8 slides Δ>0, p=0.461, tissues 3/5 Δ>0
- PDCD1LG2|PDCD1 squidpy: high 0.047 vs low 0.049, fold 0.96, Δ -0.002, 1/8 slides Δ>0, p=0.016, tissues 1/5 Δ>0
- TGFB1|TGFBR1_TGFBR2 commot: high 0.002 vs low 0.001, fold 1.25, Δ 2.77e-04, 7/8 slides Δ>0, p=0.039, tissues 4/5 Δ>0
- TGFB1|TGFBR1_TGFBR2 commot_frac_pos: high 0.009 vs low 0.007, fold 1.25, Δ 0.002, 7/8 slides Δ>0, p=0.078, tissues 4/5 Δ>0
- TGFB1|TGFBR1_TGFBR2 commot_supply: high 0.022 vs low 0.014, fold 1.56, Δ 0.008, 6/8 slides Δ>0, p=0.039, tissues 3/5 Δ>0
- TGFB1|TGFBR1_TGFBR2 expr: high 0.064 vs low 0.081, fold 0.79, Δ -0.016, 2/8 slides Δ>0, p=0.078, tissues 2/5 Δ>0
- TGFB1|TGFBR1_TGFBR2 sdm_prox: high 0.041 vs low 0.026, fold 1.61, Δ 0.024, 6/8 slides Δ>0, p=0.250, tissues 3/5 Δ>0
- TGFB1|TGFBR1_TGFBR2 squidpy: high 0.044 vs low 0.052, fold 0.85, Δ -0.008, 2/8 slides Δ>0, p=0.078, tissues 2/5 Δ>0

### barrier_tertile_15um_g0 (tertile, 15 µm)

- CD274|PDCD1 commot: high 0.003 vs low 0.003, fold 0.83, Δ -2.09e-04, 3/8 slides Δ>0, p=0.688, tissues 2/5 Δ>0
- CD274|PDCD1 commot_frac_pos: high 0.010 vs low 0.008, fold 1.19, Δ 4.90e-04, 5/8 slides Δ>0, p=0.297, tissues 3/5 Δ>0
- CD274|PDCD1 commot_supply: high 0.045 vs low 0.039, fold 1.15, Δ 0.004, 4/8 slides Δ>0, p=0.938, tissues 2/5 Δ>0
- CD274|PDCD1 expr: high 0.049 vs low 0.060, fold 0.81, Δ -0.008, 2/8 slides Δ>0, p=0.148, tissues 1/5 Δ>0
- CD274|PDCD1 sdm_prox: high 0.030 vs low 0.008, fold 3.59, Δ -0.007, 3/8 slides Δ>0, p=0.945, tissues 2/5 Δ>0
- CD274|PDCD1 squidpy: high 0.063 vs low 0.060, fold 1.04, Δ -0.002, 3/8 slides Δ>0, p=0.461, tissues 2/5 Δ>0
- CDH1|CDH1 commot: high 1.78e-05 vs low 5.37e-06, fold 3.31, Δ 9.55e-06, 7/8 slides Δ>0, p=0.016, tissues 4/5 Δ>0
- CDH1|CDH1 commot_frac_pos: high 0.028 vs low 0.011, fold 2.44, Δ 0.016, 7/8 slides Δ>0, p=0.016, tissues 4/5 Δ>0
- CDH1|CDH1 commot_supply: high 1.34e-04 vs low 5.11e-05, fold 2.63, Δ 7.49e-05, 5/8 slides Δ>0, p=0.078, tissues 3/5 Δ>0
- CDH1|CDH1 expr: high 0.154 vs low 0.102, fold 1.51, Δ 0.052, 7/8 slides Δ>0, p=0.023, tissues 5/5 Δ>0
- CDH1|CDH1 sdm_prox: high 0.037 vs low 0.027, fold 1.38, Δ 0.016, 4/8 slides Δ>0, p=0.383, tissues 2/5 Δ>0
- CDH1|CDH1 squidpy: high 0.103 vs low 0.074, fold 1.39, Δ 0.024, 7/8 slides Δ>0, p=0.016, tissues 5/5 Δ>0
- CDH1|ITGA2_ITGB1 commot: high 0.002 vs low 9.34e-04, fold 2.14, Δ 0.001, 7/8 slides Δ>0, p=0.016, tissues 4/5 Δ>0
- CDH1|ITGA2_ITGB1 commot_frac_pos: high 0.017 vs low 0.005, fold 3.32, Δ 0.012, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- CDH1|ITGA2_ITGB1 commot_supply: high 0.014 vs low 0.008, fold 1.70, Δ 0.003, 5/8 slides Δ>0, p=0.383, tissues 4/5 Δ>0
- CDH1|ITGA2_ITGB1 expr: high 0.154 vs low 0.102, fold 1.51, Δ 0.052, 7/8 slides Δ>0, p=0.023, tissues 5/5 Δ>0
- CDH1|ITGA2_ITGB1 sdm_prox: high 0.100 vs low 0.008, fold 12.24, Δ 0.092, 7/8 slides Δ>0, p=0.055, tissues 4/5 Δ>0
- CDH1|ITGA2_ITGB1 squidpy: high 0.085 vs low 0.061, fold 1.39, Δ 0.024, 7/8 slides Δ>0, p=0.016, tissues 5/5 Δ>0
- HLA-A|CD8A commot: high 0.034 vs low 0.033, fold 1.04, Δ 8.88e-04, 4/8 slides Δ>0, p=1.000, tissues 3/5 Δ>0
- HLA-A|CD8A commot_frac_pos: high 0.567 vs low 0.499, fold 1.14, Δ 0.101, 7/8 slides Δ>0, p=0.195, tissues 4/5 Δ>0
- HLA-A|CD8A commot_supply: high 0.041 vs low 0.044, fold 0.92, Δ -0.004, 4/8 slides Δ>0, p=1.000, tissues 3/5 Δ>0
- HLA-A|CD8A expr: high 0.922 vs low 0.801, fold 1.15, Δ 0.111, 7/8 slides Δ>0, p=0.195, tissues 4/5 Δ>0
- HLA-A|CD8A sdm_prox: high 0.534 vs low 0.351, fold 1.52, Δ 0.148, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- HLA-A|CD8A squidpy: high 0.756 vs low 0.684, fold 1.10, Δ 0.056, 7/8 slides Δ>0, p=0.195, tissues 4/5 Δ>0
- ICAM1|ITGAL commot: high 0.016 vs low 0.005, fold 3.13, Δ 0.008, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- ICAM1|ITGAL commot_frac_pos: high 0.031 vs low 0.013, fold 2.34, Δ 0.023, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- ICAM1|ITGAL commot_supply: high 0.076 vs low 0.058, fold 1.31, Δ 0.015, 7/8 slides Δ>0, p=0.016, tissues 4/5 Δ>0
- ICAM1|ITGAL expr: high 0.179 vs low 0.090, fold 1.98, Δ 0.083, 7/8 slides Δ>0, p=0.023, tissues 4/5 Δ>0
- ICAM1|ITGAL sdm_prox: high 0.049 vs low 0.010, fold 4.75, Δ 0.037, 5/8 slides Δ>0, p=0.461, tissues 3/5 Δ>0
- ICAM1|ITGAL squidpy: high 0.148 vs low 0.108, fold 1.37, Δ 0.032, 7/8 slides Δ>0, p=0.016, tissues 5/5 Δ>0
- ICAM1|ITGAL_ITGB2 commot: high 0.003 vs low 0.001, fold 2.79, Δ 0.001, 7/8 slides Δ>0, p=0.016, tissues 4/5 Δ>0
- ICAM1|ITGAL_ITGB2 commot_frac_pos: high 0.014 vs low 0.005, fold 2.51, Δ 0.008, 7/8 slides Δ>0, p=0.016, tissues 4/5 Δ>0
- ICAM1|ITGAL_ITGB2 commot_supply: high 0.024 vs low 0.019, fold 1.23, Δ 0.006, 6/8 slides Δ>0, p=0.156, tissues 3/5 Δ>0
- ICAM1|ITGAL_ITGB2 expr: high 0.179 vs low 0.090, fold 1.98, Δ 0.083, 7/8 slides Δ>0, p=0.023, tissues 4/5 Δ>0
- ICAM1|ITGAL_ITGB2 sdm_prox: high 0.029 vs low 0.008, fold 3.81, Δ 0.026, 4/8 slides Δ>0, p=0.742, tissues 3/5 Δ>0
- ICAM1|ITGAL_ITGB2 squidpy: high 0.106 vs low 0.066, fold 1.61, Δ 0.032, 7/8 slides Δ>0, p=0.016, tissues 5/5 Δ>0
- IFNG|IFNGR1_IFNGR2 commot: high 3.80e-04 vs low 5.16e-04, fold 0.74, Δ -1.45e-05, 3/8 slides Δ>0, p=0.812, tissues 2/5 Δ>0
- IFNG|IFNGR1_IFNGR2 commot_frac_pos: high 0.001 vs low 0.001, fold 1.03, Δ 6.60e-04, 4/8 slides Δ>0, p=0.375, tissues 3/5 Δ>0
- IFNG|IFNGR1_IFNGR2 commot_supply: high 0.008 vs low 0.010, fold 0.82, Δ -2.21e-04, 3/8 slides Δ>0, p=0.938, tissues 2/5 Δ>0
- IFNG|IFNGR1_IFNGR2 expr: high 0.046 vs low 0.040, fold 1.15, Δ 0.004, 5/8 slides Δ>0, p=0.461, tissues 3/5 Δ>0
- IFNG|IFNGR1_IFNGR2 sdm_prox: high -0.010 vs low -0.004, fold NA, Δ -0.003, 4/8 slides Δ>0, p=0.844, tissues 3/5 Δ>0
- IFNG|IFNGR1_IFNGR2 squidpy: high 0.032 vs low 0.027, fold 1.18, Δ 0.003, 6/8 slides Δ>0, p=0.078, tissues 4/5 Δ>0
- LGALS9|HAVCR2 commot: high 0.003 vs low 0.002, fold 1.38, Δ -3.29e-05, 3/8 slides Δ>0, p=1.000, tissues 2/5 Δ>0
- LGALS9|HAVCR2 commot_frac_pos: high 0.016 vs low 0.012, fold 1.34, Δ 0.003, 6/8 slides Δ>0, p=0.109, tissues 3/5 Δ>0
- LGALS9|HAVCR2 commot_supply: high 0.034 vs low 0.030, fold 1.14, Δ -0.009, 2/8 slides Δ>0, p=0.938, tissues 2/5 Δ>0
- LGALS9|HAVCR2 expr: high 0.082 vs low 0.083, fold 0.99, Δ 0.004, 5/8 slides Δ>0, p=0.547, tissues 3/5 Δ>0
- LGALS9|HAVCR2 sdm_prox: high 0.046 vs low 0.034, fold 1.35, Δ -0.009, 4/8 slides Δ>0, p=0.945, tissues 3/5 Δ>0
- LGALS9|HAVCR2 squidpy: high 0.071 vs low 0.074, fold 0.97, Δ -1.12e-04, 3/8 slides Δ>0, p=0.844, tissues 2/5 Δ>0
- MIF|CD74_CD44 commot: high 0.030 vs low 0.024, fold 1.25, Δ 0.003, 6/8 slides Δ>0, p=0.109, tissues 4/5 Δ>0
- MIF|CD74_CD44 commot_frac_pos: high 0.153 vs low 0.111, fold 1.38, Δ 0.047, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- MIF|CD74_CD44 commot_supply: high 0.061 vs low 0.077, fold 0.80, Δ -0.003, 2/8 slides Δ>0, p=0.547, tissues 2/5 Δ>0
- MIF|CD74_CD44 expr: high 0.449 vs low 0.319, fold 1.41, Δ 0.085, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- MIF|CD74_CD44 sdm_prox: high 0.243 vs low 0.166, fold 1.46, Δ 0.050, 5/8 slides Δ>0, p=0.148, tissues 4/5 Δ>0
- MIF|CD74_CD44 squidpy: high 0.354 vs low 0.291, fold 1.22, Δ 0.040, 8/8 slides Δ>0, p=0.008, tissues 5/5 Δ>0
- PDCD1LG2|PDCD1 commot: high 0.002 vs low 0.002, fold 0.73, Δ -4.89e-04, 2/8 slides Δ>0, p=0.578, tissues 1/5 Δ>0
- PDCD1LG2|PDCD1 commot_frac_pos: high 0.006 vs low 0.007, fold 0.79, Δ 9.18e-04, 4/8 slides Δ>0, p=1.000, tissues 2/5 Δ>0
- PDCD1LG2|PDCD1 commot_supply: high 0.064 vs low 0.085, fold 0.75, Δ -0.013, 3/7 slides Δ>0, p=1.000, tissues 1/5 Δ>0
- PDCD1LG2|PDCD1 expr: high 0.035 vs low 0.032, fold 1.10, Δ 0.004, 5/8 slides Δ>0, p=0.383, tissues 4/5 Δ>0
- PDCD1LG2|PDCD1 sdm_prox: high 0.033 vs low 0.060, fold 0.55, Δ -0.026, 3/8 slides Δ>0, p=0.688, tissues 1/5 Δ>0
- PDCD1LG2|PDCD1 squidpy: high 0.049 vs low 0.046, fold 1.07, Δ 0.001, 6/8 slides Δ>0, p=0.383, tissues 4/5 Δ>0
- TGFB1|TGFBR1_TGFBR2 commot: high 0.002 vs low 0.002, fold 1.25, Δ 6.65e-04, 6/8 slides Δ>0, p=0.461, tissues 3/5 Δ>0
- TGFB1|TGFBR1_TGFBR2 commot_frac_pos: high 0.011 vs low 0.006, fold 1.68, Δ 0.003, 7/8 slides Δ>0, p=0.195, tissues 4/5 Δ>0
- TGFB1|TGFBR1_TGFBR2 commot_supply: high 0.025 vs low 0.028, fold 0.90, Δ 0.021, 5/8 slides Δ>0, p=0.383, tissues 2/5 Δ>0
- TGFB1|TGFBR1_TGFBR2 expr: high 0.064 vs low 0.082, fold 0.79, Δ -0.007, 3/8 slides Δ>0, p=0.641, tissues 3/5 Δ>0
- TGFB1|TGFBR1_TGFBR2 sdm_prox: high 0.034 vs low 0.013, fold 2.70, Δ 0.026, 5/8 slides Δ>0, p=0.641, tissues 2/5 Δ>0
- TGFB1|TGFBR1_TGFBR2 squidpy: high 0.043 vs low 0.054, fold 0.80, Δ -0.009, 3/8 slides Δ>0, p=0.383, tissues 2/5 Δ>0

SpatialDM length scale is set so the RBF weight exp(−d² / 2l²) falls through 0.15 near the selected radius (`l = radius / sqrt(−2 ln 0.15)`). COMMOT `dis_thr` is that same radius. Squidpy `gr.ligrec` is restricted to the same radius. Caps, permutation count (1,000), and the CellChat pair list match v1.

A Moran I delta stays on a correlation scale. A COMMOT per-sender delta stays on a transport-mass scale and will look small even when the high/low fold is not. The fraction of senders with any positive transport is the transport score on a probability scale.

