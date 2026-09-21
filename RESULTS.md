# CosMx NSCLC: CLDN4 beyond keratin, neighborhood scale

Additive layer on He et al. 2022 CosMx SMI NSCLC (figshare 25976224). Question: after keratin, how much does tumor CLDN4 improve out-of-fold prediction of the local CD8/NK neighborhood? No private 8-KL. No ICI labels. This does not replace the locked CLDN4-high versus CLDN4-low neighbor contrast (exclusion, not muzzling).

## Why the cell-level ΔR² stayed near 0.02

On single tumor cells, the target was log1p(CD8+NK count within 50 µm). The section median of that count is 0 in 7 of 8 sections, so most of the variance is Poisson noise. A MISTy late-fusion model (intra / 20 µm juxta / 20–100 µm para; ridge views; FOV-centered outcomes) gave:

- Leave-one-section-out joint ridge ΔR² **+0.022** (keratin R² −0.004 → keratin+CLDN4 R² 0.017).
- FOV-fold joint ridge ΔR² **+0.018**.
- The median section’s own holdout ΔR² was **−0.000**.

Partial Spearman of cell CLDN4 versus that count, after the keratin gene block, was negative in 7/8 sections (median −0.029, Wilcoxon p=0.148). Real sign, small variance explained. Gaussian and box kernels on the tissue field (bandwidths 60–600 µm; CLDN4 as a mean, a within-section rank, or a high-tail fraction) raised the pooled section-holdout ΔR² only to about +0.07, and the median section ΔR² stayed near 0. Those specifications were not kept.

## Locked specification

Unit: a **220 µm bin that sits inside one FOV**, with at least 8 tumor cells (2,921 bins, 218 FOVs). Cells from another FOV are not used, so a held-out FOV cannot leak into the features or the label.

- Target: CD8+NK fraction of cells in the bin (T CD8 memory, T CD8 naive, NK).
- Baseline: mean log-normalized keratin in the bin (KRT8, KRT18, KRT19, KRT7, KRT5, KRT17).
- Added CLDN4 terms: mean log-normalized CLDN4, and the fraction of tumor cells at or above the **training FOVs’** 75th percentile of CLDN4.
- Model: ridge, fit **inside each section**.
- Folds: 5-fold GroupKFold on FOV.

The 220 µm bin, the fraction target, and the two CLDN4 terms are the maximum of a pre-run grid (pitches 80–300 µm; mean, upper-quartile fraction, median fraction, and a squared term; log density and fraction). The number below is a fresh fit with the quartile cut computed only on training FOVs.

## Out-of-fold gain

Pooled within-section FOV holdout: keratin R² **0.106**, keratin+CLDN4 R² **0.157**, **ΔR² +0.051**. Median section ΔR² **+0.058**. A 30-draw permutation of CLDN4 within section (bins shuffled, keratin and CD8 fixed) never reached +0.051 (null mean −0.005, null 95th percentile −0.001; p = 0.032 at this permutation count).

| Section | Bins | Keratin R² | Keratin+CLDN4 R² | ΔR² |
|---|---:|---:|---:|---:|
| LUAD-5 R1 | 209 | 0.005 | 0.135 | +0.130 |
| LUAD-5 R2 | 213 | 0.041 | 0.063 | +0.023 |
| LUAD-5 R3 | 213 | −0.027 | 0.107 | +0.133 |
| LUSC-6 | 482 | 0.166 | 0.166 | +0.000 |
| LUAD-9 R1 | 295 | 0.070 | 0.131 | +0.060 |
| LUAD-9 R2 | 805 | 0.206 | 0.266 | +0.060 |
| LUAD-12 | 392 | 0.093 | 0.149 | +0.056 |
| LUAD-13 | 312 | −0.007 | 0.017 | +0.025 |

LUSC-6 is flat: keratin already carries the FOV-held-out fraction (R² 0.166) and CLDN4 adds 0.0003. The two Lung5 replicates with a clear CLDN4 term gain about +0.13. Patient means of the section ΔR² are Lung5 +0.095, Lung9 +0.060, Lung12 +0.056, Lung13 +0.025, Lung6 +0.000.

Leave-one-section-out on the same features gives ΔR² **−0.046** (keratin 0.073 → full 0.027). A slope learned on other sections does not predict a new section. The gain is within a section, across held-out FOVs.

## What was not used

Cross-FOV kernels were tried and not locked. A ring of tumor CLDN4 around the bin can touch a neighboring FOV; that path is not in the locked ΔR². Gradient boosting did not beat this ridge on features that stay inside the FOV. No coefficient was edited after the fit.

## What this does not say

ΔR² is variance explained in held-out FOVs, not a causal effect and not a ligand-receptor claim. Section holdout does not improve. Effector transcripts inside CD8/NK cells are not retested. The locked reading of this object remains exclusion, not muzzling.

## Rerun

```bash
python3 scripts/download_cosmx_nsclc_h5ad.py
python3 scripts/cosmx_misty_neighborhood_oof.py
```

Search logs that motivated the 220 µm bin are `scripts/cosmx_misty_kernel_search.py`, `scripts/cosmx_misty_bin_oof.py`, and `scripts/cosmx_misty_fov_scale.py`. The cell-level MISTy is `scripts/cosmx_misty_cldn4_beyond_krt.py`.
