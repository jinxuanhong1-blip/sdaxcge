# CosMx NSCLC: largest leak-free FOV-holdout ΔR² for CLDN4 beyond keratin

Additive layer on He et al. 2022 CosMx SMI NSCLC (figshare 25976224). Question: after keratin, how much does tumor CLDN4 improve out-of-fold prediction of the local immune neighborhood? No private 8-KL. No ICI labels. The locked neighbor-count contrast is unchanged: CLDN4-high tumor has fewer nearby CD8/NK cells (exclusion, not muzzling).

## Grid maximum

The score is min(pooled ΔR², median-section ΔR²). Eligible specifications keep all 8 sections and both deltas above 0. A within-section shuffle of CLDN4 has to stay below the pooled gain. The search that produced the number below is 2,800 unique fits (5,600 rows, because a filter of ≥400 cells per FOV matched the unfiltered grid: every FOV already had at least 400 cells). The reported gain is that maximum. The permutation p-value is for the chosen specification only.

**Pooled FOV-holdout ΔR² +0.164** (keratin R² 0.115 → keratin+CLDN4 R² 0.279). **Median section ΔR² +0.157**. Six of eight sections have a positive ΔR². The largest section gain is LUAD-5 R3 **+0.275**.

Specification:

- Bins of 300 µm that sit inside one FOV, at least 15 tumor cells (2,010 bins).
- Kernel: Gaussian σ = 80 µm, truncated at 3σ, using only cells in that FOV. Keratin and CLDN4 are tumor-cell means. The NK target is the same kernel on all cells.
- Baseline: that keratin mean (KRT8, KRT18, KRT19, KRT7, KRT5, KRT17).
- Added: the CLDN4 mean, and the fraction of tumor cells in the bin at or above the training FOVs’ 80th percentile of CLDN4.
- Target: NK fraction.
- Ridge fit inside each section. Folds: 5-fold GroupKFold on FOV.

A hard 300 µm bin with the same NK target and a mean plus training-fold p75 CLDN4 term, and no overlapping kernel, reaches pooled ΔR² +0.141 and median-section ΔR² +0.140 (keratin R² 0.079 → 0.220). Averaging the Gaussian model’s out-of-fold predictions to one row per FOV (212 FOVs) gives pooled ΔR² **+0.309** (keratin R² 0.179 → 0.488). The gain is still there when each FOV is one point.

Forty within-section shuffles of CLDN4, with keratin and the NK target fixed, have null mean −0.010 and null maximum +0.006. None reached +0.164 (p = 0.024 at this permutation count).

| Section | Bins | Keratin R² | Keratin+CLDN4 R² | ΔR² | Partial Spearman |
|---|---:|---:|---:|---:|---:|
| LUAD-5 R1 | 149 | 0.044 | 0.218 | +0.174 | +0.444 |
| LUAD-5 R2 | 150 | 0.086 | 0.229 | +0.144 | +0.378 |
| LUAD-5 R3 | 147 | −0.018 | 0.257 | +0.275 | +0.458 |
| LUSC-6 | 334 | 0.077 | 0.055 | −0.021 | +0.065 |
| LUAD-9 R1 | 211 | 0.053 | 0.266 | +0.213 | −0.349 |
| LUAD-9 R2 | 534 | 0.245 | 0.415 | +0.170 | −0.376 |
| LUAD-12 | 256 | 0.056 | 0.088 | +0.032 | −0.248 |
| LUAD-13 | 229 | −0.037 | −0.073 | −0.036 | +0.095 |

Partial Spearman is CLDN4 versus the NK fraction after a linear keratin residual, on the same Gaussian scale. The median section partial is **+0.080**. Three sections are negative (LUAD-9 R1, LUAD-9 R2, LUAD-12). Lung5, which carries the largest ΔR² values, is positive: higher CLDN4, higher NK fraction. LUSC-6 and LUAD-13 add no gain. Blue bars in the figure are sections where the ΔR² is positive and the partial is negative. Gold bars are positive ΔR² with a positive partial.

## Leave-one-section-out

The same features, with each section’s percentile taken from its own tumor cells and the slope fit on the other sections, give pooled ΔR² **−0.027** (keratin R² 0.021 → full R² −0.005). The pooled number recenters predictions inside the held-out section. The raw section R² values do not:

| Section | Keratin R² | Keratin+CLDN4 R² | ΔR² |
|---|---:|---:|---:|
| LUAD-5 R1 | −0.223 | −0.180 | +0.043 |
| LUAD-5 R2 | −0.201 | −0.089 | +0.112 |
| LUAD-5 R3 | −0.285 | −0.203 | +0.082 |
| LUSC-6 | −1.989 | −7.493 | −5.504 |
| LUAD-9 R1 | −0.022 | −0.101 | −0.079 |
| LUAD-9 R2 | −0.015 | −0.136 | −0.121 |
| LUAD-12 | 0.043 | −0.142 | −0.186 |
| LUAD-13 | −2.000 | −1.692 | +0.308 |

LUSC-6’s held-out R² falls to −7.5. Lung9 and Lung12, the sections whose within-section partials are negative, lose R² when the slope is learned elsewhere. The best within-section gain on this specification (LUAD-5 R3, +0.275) does not survive as a shared slope. The largest single-section ΔR² anywhere in the grid is +0.605 (400 µm, Gaussian σ = 160, NK, mean plus Q4), and that specification’s median-section ΔR² is only +0.066, so the score did not select it.

## What the sign does to the exclusion question

Specs whose partial Spearman is negative in at least 6 of 8 sections top out at pooled ΔR² **+0.017** and median-section ΔR² **+0.018** (CD8 fraction, 400 µm bins, Gaussian σ = 160, CLDN4 mean plus % ≥ p80; 4 of 8 section deltas positive). That is smaller than the 220 µm CD8+NK bin result below. The large NK gain and the exclusion neighbor-count result are different summaries.

## Earlier scales, same object

Cell-level log1p(CD8+NK count within 50 µm): leave-one-section-out joint ridge ΔR² **+0.022**, FOV-fold ΔR² **+0.018**, median section holdout about 0. Partial Spearman after the keratin genes was negative in 7/8 sections (median −0.029).

220 µm bins, CD8+NK fraction, keratin mean plus training-fold p75 CLDN4, FOV holdout inside the section: pooled ΔR² **+0.051** (keratin R² 0.106 → 0.157), median section ΔR² **+0.058**, 8/8 sections ≥ 0. On that same bin scale the median-section partial Spearman is **+0.058** (4/8 negative). Leave-one-section-out on those features is ΔR² **−0.046**.

## What this does not say

ΔR² is variance explained in held-out FOVs. It is unsigned, so a section can gain R² with either sign. The grid maximum was chosen on the same out-of-fold score it reports. Section holdout does not improve. Effector transcripts inside CD8/NK cells are not retested. The locked reading of the neighbor counts remains exclusion, not muzzling.

## Rerun

```bash
python3 scripts/download_cosmx_nsclc_h5ad.py
python3 scripts/cosmx_misty_max_effect.py
```

The 220 µm CD8+NK fit is `scripts/cosmx_misty_neighborhood_oof.py`. The cell-level MISTy is `scripts/cosmx_misty_cldn4_beyond_krt.py`.
