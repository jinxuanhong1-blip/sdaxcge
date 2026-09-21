# CosMx NSCLC: MISTy multi-view importance of CLDN4 beyond keratin

Additive public layer on He et al. 2022 CosMx SMI NSCLC (figshare 25976224, `cosmx_human_nsclc_clustered.h5ad`). It asks whether CLDN4 in tumor cells predicts local CD8/NK density after keratin is in the model. It does not replace the locked CLDN4-high versus CLDN4-low neighbor result (exclusion, not muzzling). No private 8-KL matrices. No ICI labels.

## Cohort and definitions

Author labels: tumor index = tumor 5, tumor 6, tumor 9, tumor 12, tumor 13 (n=302,313 tumor cells). CD8/NK = T CD8 memory + T CD8 naive + NK. Normal `epithelial` is not in the index. Coordinates are pixels; 0.18 µm/px because each FOV spans about 5440 × 3620 px (0.98 × 0.65 mm), the CosMx FOV. Interior cells are those whose neighborhood ball sits inside the FOV union after closing ~20 µm tiling gaps.

Primary target: log1p(CD8+NK count within 50 µm), monotone with density at a fixed radius. Density in cells/mm² is stored in the cell table used to build the models. Keratin block: KRT8, KRT18, KRT19, KRT7, KRT5, KRT17, KRT6A, KRT6B, KRT6C. Expression is log1p of library-size normalized counts (scale = median total counts).

Views, in the MISTy layout (intraview / juxtaview / paraview, zone of indifference so juxta and para do not share neighbors):

- Intraview: CLDN4 and keratin in the index tumor cell.
- Juxtaview: mean CLDN4 and keratin in other tumor cells within 20 µm.
- Paraview: inverse-distance-weighted mean in other tumor cells at 20–100 µm.

Folds: leave-one-sample-out (8 sections) and 8-fold GroupKFold on FOV (`sample::fov`). Predictors and the outcome are centered within FOV inside each fold, so the fit is a within-field slope. Held-out R² is computed after centering the prediction and the outcome in the held-out sample (sample folds) or FOV (FOV folds).

Interior tumor cells at 50 µm: 296,521 (296,403 enter the models after dropping FOVs with fewer than 25 interior tumor cells). At 100 µm: 282,751. Sections: 8. Patients: 5. Within-section Spearman of CLDN4 versus the KRT8/18/19 score on interior cells ranges from 0.061 to 0.290, so CLDN4 is not a keratin duplicate.

## Primary result

Within each section, after FOV centering, partial Spearman of intraview CLDN4 versus 50 µm log1p(CD8+NK count), residualizing both on the keratin gene block: **7/8 sections negative**, median -0.029 (sample bootstrap 95% interval -0.080 to -0.007). Exact Wilcoxon signed-rank two-sided p=0.148, one-sided (less) p=0.074. The positive section is LUSC-6 (squamous; CLDN4 detected in 24.5% of its tumor cells). Lung9 and Lung12 carry the negative rank association. Lung5 is near zero.

| Section | Patient | n | Marginal Spearman | Partial Spearman \| KRT | Standardized beta |
|---|---|---:|---:|---:|---:|
| LUAD-5 R1 | Lung5 | 17998 | -0.001 | -0.007 | 0.018 |
| LUAD-5 R2 | Lung5 | 18533 | -0.011 | -0.023 | -0.013 |
| LUAD-5 R3 | Lung5 | 16006 | -0.024 | -0.018 | 0.001 |
| LUSC-6 | Lung6 | 65767 | 0.044 | 0.085 | 0.031 |
| LUAD-9 R1 | Lung9 | 38634 | -0.100 | -0.094 | -0.100 |
| LUAD-9 R2 | Lung9 | 95492 | -0.104 | -0.080 | -0.095 |
| LUAD-12 | Lung12 | 18176 | -0.125 | -0.060 | -0.063 |
| LUAD-13 | Lung13 | 25797 | -0.039 | -0.036 | -0.036 |

Standardized partial beta (SD of the FOV-centered outcome per SD of CLDN4, keratin genes in the same model): median -0.024, 5/8 sections negative, two-sided p=0.195. Lung5 replicate 1 and replicate 3 have a negative rank correlation and a beta indistinguishable from zero.

Marginal Spearman, same outcome, no keratin adjustment: median -0.032, 7/8 negative, two-sided p=0.078. Keratin adjustment shifts the median from -0.032 to -0.029.

Patient means of the partial Spearman (Lung5 and Lung9 replicates averaged): 4/5 patients negative, median -0.036, two-sided p=0.438.

Out-of-fold gain from adding the three CLDN4 views on top of the three keratin views (ΔR²):

- Sample folds, joint ridge: KRT R² -0.0045, KRT+CLDN4 R² 0.0173, ΔR² +0.0218.
- FOV folds, joint ridge: KRT R² 0.0322, KRT+CLDN4 R² 0.0502, ΔR² +0.0180.
- Sample folds, MISTy late fusion (ridge view models, ridge meta-model): KRT R² -0.0015, KRT+CLDN4 R² 0.0208, ΔR² +0.0222.
- FOV folds, MISTy late fusion: KRT R² 0.0323, KRT+CLDN4 R² 0.0488, ΔR² +0.0165.

Pooled sample-holdout ΔR² is +0.022. The median section's leave-one-section-out ΔR² is -0.0004. Lung9 R1 and R2 each gain about +0.04 R². Lung5 holdouts lose R² when CLDN4 is added. LUSC-6 stays worse than its own mean (R² -0.30 with keratin, -0.20 after CLDN4). FOV holdout estimates within-study prediction (+0.018 joint ridge, +0.017 MISTy). A 100 µm neighborhood can cross a FOV boundary, so adjacent fields are not fully sealed. Sample holdout does not use cells from the held-out section.

Sample-holdout fusion weights sum to 0.073 on the keratin views and 0.066 on the CLDN4 views (47.5% of the summed weight). The paraview is the largest scale for both blocks. A positive fusion weight means the meta-model uses that view's prediction; the biological direction is the partial correlation above, which is negative in 7/8 sections.

Local CD8/NK counts are sparse (section median count at 50 µm is 0 except LUAD-13), so absolute R² stays small. This layer reports the keratin-adjusted sign and the incremental R², not a reconstruction of the neighborhood.

## Sensitivities

Keratin control restricted to the mean of KRT8, KRT18, and KRT19 (the TCGA-style simple-keratin score), 50 µm: partial Spearman median -0.032, 7/8 negative, two-sided p=0.055.

Intraview CLDN4 residualized on all three keratin views (intra, juxta, para), 50 µm: partial Spearman median -0.031, 7/8 negative, two-sided p=0.195.

Within-section incremental R² of the CLDN4 view block after the keratin view block, 50 µm: median 0.0054. Median partial R² (share of leftover variance after keratin views) 0.0055. In-sample R² cannot be negative once CLDN4 is added, so the sign test is the partial correlation above, and the out-of-fold ΔR² is the predictive check.

Composition rather than area-density: CD8+NK fraction of all neighbors within 50 µm, partial Spearman | keratin genes: median -0.034, 7/8 negative, two-sided p=0.195.

Same density target after keratin genes plus log1p(local cell count): partial Spearman median -0.032, 7/8 negative, two-sided p=0.148.

100 µm log1p(CD8+NK count), partial Spearman | keratin genes: median -0.007, 4/8 negative, two-sided p=0.461.

100 µm sample-fold joint ridge: KRT R² -0.0277, KRT+CLDN4 R² -0.0020, ΔR² +0.0257.
100 µm FOV-fold joint ridge: KRT R² 0.0564, KRT+CLDN4 R² 0.0767, ΔR² +0.0203.

Section-level partial correlations are in `tables/sample_partial_cldn4.csv`. View weights are in `tables/misty_view_coefficients.csv`.

## What this does not say

Predictive contribution is not a causal effect of CLDN4, and it is not a ligand-receptor claim. CD8/NK cells are the author T CD8 and NK classes; effector transcripts inside those cells (GZMB, PRF1, IFNG) are not re-tested here. The locked reading of this object remains exclusion, not muzzling. Visium same-spot correlations are not used.

## How to rerun

```bash
python3 scripts/download_cosmx_nsclc_h5ad.py
python3 scripts/cosmx_misty_cldn4_beyond_krt.py
```

