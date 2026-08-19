# GSE307534 Visium: CLDN4 vs CD8 spatial statistics

Public dataset only (GEO [GSE307534](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307534)): Visium CytAssist of normal lung, AAH, AIS, MIA, and invasive LUAD. Analysis is **CLDN4-only**. No private gene sets were used.

## Question

Do CLDN4-high tumor/epithelial spots spatially exclude CD8?

## Data and gene check

- Sections with processed spots downloaded: 56.
- Sections analyzed (CLDN4 and CD8A present): 56.
- Sections skipped (both or one of CLDN4/CD8A absent, or load error): 0.
- CLDN4 present in 56 / 56 downloaded sections.
- CD8A present in 56 / 56 downloaded sections.
- NKG7 present in 56 / 56; GNLY present in 56 / 56.
- EPCAM present in 56 / 56; KRT8 present in 56 / 56.

## Methods (computed, not assumed)

- Filtered Space Ranger matrices and `tissue_positions.csv` from each GSM tar; H&E images were not used.
- Expression: log1p(10⁴ × counts / spot UMI). Duplicate gene symbols were summed.
- Epithelial-like spots: top quartile of the mean z-score of available EPCAM and KRT8.
- CLDN4-high / low among epithelial-like spots: top vs bottom quartile (Q4 vs Q1).
- CD8A-high: spots at or above the 75th percentile of CD8A (or CD8A>0 if that percentile is 0).
- Spatial weights: row-standardized 6-nearest neighbors in micron coordinates (55 µm / `spot_diameter_fullres`).
- Bivariate Moran I(x,y) = z_x' W z_y / z_x' z_x with mean-centered (not variance-standardized) x,y; |I| and |L| can exceed 1 when CLDN4 and CD8A have different variances. Lee's L as in Lee (2001) with row-standardized W.
- Residual spatial lag: OLS-residualize CLDN4 and CD8A on KRT8, then Spearman of residual CLDN4 vs W·residual CD8A.
- Getis-Ord Gi* with binary kNN including self; hotspot/coldspot |z|>1.96. Overlap = CLDN4 hot ∩ CD8A cold.
- Cross-K / L₁₂(r) between CLDN4 Q4 epithelial-like points and CD8A-high points; L(r)=√(K/π)−r; summary radius = 2 × median nearest-neighbor distance.
- Q4 vs Q1: mean distance (µm) to nearest CD8A-high spot, and mean CD8A sum in 6-NN.
- Tumor domain (when defined): largest 6-NN connected component of epithelial-like spots with ≥40 spots. Hop-depth from non-domain spots; margin = depth 1; CLDN4-high core = CLDN4 Q4 among domain spots and depth ≥ max(2, domain-depth Q3).
- Empirical p: 199 permutations of CD8A among epithelial-like spots (non-epithelial CD8A left in place). One-sided p is in the exclusion direction (negative association / larger Q4 distance / more Gi* overlap / lower core CD8A).
- Primary cut is invasive LUAD if those sections are available; precursor vs invasive is a secondary cut.
- No private 8-KL or other private signatures.

## Histology check (whether precursors are normal-like for CLDN4)

Median across sections of per-section median CLDN4 in epithelial-like spots:

- Normal: 0 (n=1 section)
- AAH: 0 (n=11 sections)
- AIS: 0 (n=14 sections)
- MIA: 0.0115 (n=4 sections)
- LUAD: 0.101 (n=26 sections)

AAH/AIS median CLDN4 in epithelial-like spots matches the single Normal section (0 on the log-norm scale). MIA is intermediate. Primary reporting set is therefore **invasive LUAD** (n=26). Precursor vs invasive remains the secondary cut.

## Meta results

### Primary: invasive LUAD

n sections = 26.

| statistic | n | median [Q1, Q3] | n<0 | n>0 | Wilcoxon p vs 0 | median epi-perm p (exclusion dir.) | n perm p<0.05 |
|---|---:|---|---:|---:|---:|---:|---:|
| Spearman CLDN4 vs CD8A (all spots) | 26 | 0.000992 [-0.108, 0.0937] | 13 | 13 | 0.727 | NA | NA |
| Spearman CLDN4 vs CD8A (epi-like spots) | 26 | 0.0398 [-0.17, 0.108] | 12 | 14 | 1 | 0.96 | 11 |
| Spearman CLDN4 vs CD8A+NKG7 | 26 | -0.0326 [-0.148, 0.0863] | 14 | 12 | 0.437 | NA | NA |
| Spearman after KRT8 residual | 26 | -0.0325 [-0.0715, 0.039] | 16 | 10 | 0.269 | NA | NA |
| CLDN4 residual vs spatial lag of CD8A residual | 26 | -0.0297 [-0.0497, 0.0475] | 17 | 9 | 0.653 | 0.215 | 12 |
| Bivariate Moran I(CLDN4, CD8A) | 26 | -0.0167 [-0.888, 0.632] | 14 | 12 | 0.784 | 0.258 | 13 |
| Lee's L(CLDN4, CD8A) | 26 | -0.0847 [-3.02, 1.41] | 14 | 12 | 0.689 | 0.27 | 13 |
| Δ nearest-µm (CLDN4 Q4 − Q1) to CD8A-high | 26 | -10.8 [-38.9, 22.9] | 16 | 10 | 0.34 | 0.728 | 9 |
| Δ kNN CD8A sum (Q4 − Q1) | 26 | -0.0635 [-0.307, 0.688] | 14 | 12 | 0.745 | 0.52 | 10 |
| Gi* overlap n (CLDN4 hot ∩ CD8A cold) | 26 | 106 [34, 516] | 0 | 20 | 8.86e-05 | 0.005 | 20 |
| Gi* overlap minus independence expectation | 26 | 0 [-84.2, 136] | 10 | 10 | 0.263 | NA | NA |
| Cross-L at 2× median NN | 26 | -12.3 [-37, 20.4] | 16 | 10 | 0.34 | 0.06 | 13 |
| Mean CD8A in CLDN4-high core − margin | 26 | -0.0469 [-0.221, 0.0517] | 16 | 10 | 0.0941 | 0.155 | 11 |
| Spearman CD8A vs hop-depth (tumor domain) | 26 | -0.0559 [-0.143, 0.0448] | 16 | 10 | 0.0382 | 0.025 | 14 |

### Secondary: precursor (AAH / AIS / MIA)

n sections = 29.

| statistic | n | median [Q1, Q3] | n<0 | n>0 | Wilcoxon p vs 0 | median epi-perm p (exclusion dir.) | n perm p<0.05 |
|---|---:|---|---:|---:|---:|---:|---:|
| Spearman CLDN4 vs CD8A (all spots) | 29 | 0.0701 [0.00779, 0.113] | 7 | 22 | 8.72e-05 | NA | NA |
| Spearman CLDN4 vs CD8A (epi-like spots) | 29 | 0.0565 [-0.0268, 0.0998] | 9 | 20 | 0.076 | 1 | 7 |
| Spearman CLDN4 vs CD8A+NKG7 | 29 | 0.0786 [0.00537, 0.114] | 6 | 23 | 0.00011 | NA | NA |
| Spearman after KRT8 residual | 29 | 0.0118 [-0.0267, 0.0515] | 13 | 16 | 0.169 | NA | NA |
| CLDN4 residual vs spatial lag of CD8A residual | 29 | -0.0257 [-0.0415, -0.00554] | 23 | 6 | 0.00198 | 0.205 | 11 |
| Bivariate Moran I(CLDN4, CD8A) | 29 | 0.359 [-0.231, 1.94] | 9 | 20 | 0.0386 | 0.575 | 10 |
| Lee's L(CLDN4, CD8A) | 29 | 0.422 [-1.11, 4.75] | 9 | 20 | 0.0798 | 0.365 | 10 |
| Δ nearest-µm (CLDN4 Q4 − Q1) to CD8A-high | 29 | -12.6 [-24.3, -0.0892] | 22 | 7 | 0.0689 | 0.99 | 7 |
| Δ kNN CD8A sum (Q4 − Q1) | 29 | 0.0504 [-0.0846, 0.483] | 11 | 18 | 0.247 | 0.795 | 9 |
| Gi* overlap n (CLDN4 hot ∩ CD8A cold) | 29 | 53 [19, 116] | 0 | 24 | 1.82e-05 | 0.005 | 21 |
| Gi* overlap minus independence expectation | 29 | -8.45 [-41.8, 0] | 17 | 7 | 0.511 | NA | NA |
| Cross-L at 2× median NN | 29 | 9 [-8.76, 20.1] | 13 | 16 | 0.23 | 0.785 | 8 |
| Mean CD8A in CLDN4-high core − margin | 29 | -0.0258 [-0.0981, 0.0513] | 17 | 12 | 0.19 | 0.31 | 5 |
| Spearman CD8A vs hop-depth (tumor domain) | 29 | -0.0416 [-0.106, 0.0243] | 19 | 10 | 0.0243 | 0.19 | 10 |

### All analyzed sections

n sections = 56.

| statistic | n | median [Q1, Q3] | n<0 | n>0 | Wilcoxon p vs 0 | median epi-perm p (exclusion dir.) | n perm p<0.05 |
|---|---:|---|---:|---:|---:|---:|---:|
| Spearman CLDN4 vs CD8A (all spots) | 56 | 0.0381 [-0.0301, 0.106] | 20 | 36 | 0.0353 | NA | NA |
| Spearman CLDN4 vs CD8A (epi-like spots) | 56 | 0.0408 [-0.0369, 0.103] | 22 | 34 | 0.26 | 0.985 | 19 |
| Spearman CLDN4 vs CD8A+NKG7 | 56 | 0.0328 [-0.0426, 0.107] | 20 | 36 | 0.103 | NA | NA |
| Spearman after KRT8 residual | 56 | -0.00672 [-0.0453, 0.0453] | 30 | 26 | 0.987 | NA | NA |
| CLDN4 residual vs spatial lag of CD8A residual | 56 | -0.0272 [-0.0478, 0.0248] | 41 | 15 | 0.0184 | 0.177 | 24 |
| Bivariate Moran I(CLDN4, CD8A) | 56 | 0.00516 [-0.365, 1.1] | 24 | 32 | 0.221 | 0.463 | 24 |
| Lee's L(CLDN4, CD8A) | 56 | 0.0365 [-2.17, 3] | 24 | 32 | 0.396 | 0.357 | 24 |
| Δ nearest-µm (CLDN4 Q4 − Q1) to CD8A-high | 56 | -11.6 [-28, 12.4] | 38 | 18 | 0.0574 | 0.975 | 17 |
| Δ kNN CD8A sum (Q4 − Q1) | 56 | 0.0262 [-0.245, 0.59] | 26 | 30 | 0.374 | 0.74 | 20 |
| Gi* overlap n (CLDN4 hot ∩ CD8A cold) | 56 | 73 [26.5, 194] | 0 | 45 | 5.17e-09 | 0.005 | 42 |
| Gi* overlap minus independence expectation | 56 | 0 [-57.3, 73] | 27 | 18 | 0.731 | NA | NA |
| Cross-L at 2× median NN | 56 | -2.71 [-21.5, 20.2] | 30 | 26 | 0.948 | 0.532 | 22 |
| Mean CD8A in CLDN4-high core − margin | 56 | -0.0386 [-0.147, 0.0521] | 34 | 22 | 0.0233 | 0.223 | 17 |
| Spearman CD8A vs hop-depth (tumor domain) | 56 | -0.0519 [-0.128, 0.0258] | 36 | 20 | 0.0016 | 0.165 | 25 |

## Observed direction on the primary set

- n=26 sections.
- Median Spearman CLDN4 vs CD8A = 0.000992 (13 negative, 13 positive; Wilcoxon p=0.727).
- Median bivariate Moran I = -0.0167; median Lee's L = -0.0847.
- Median residual CLDN4 vs lag residual CD8A ρ = -0.0297.
- Median Δ nearest-µm (Q4 − Q1 to CD8A-high) = -10.8 (positive = Q4 farther from CD8A-high).
- Median Gi* overlap count (CLDN4 hot ∩ CD8A cold) = 106. Median (overlap − independence expectation) = 0 (10 sections below, 10 above, 6 ties at 0; Wilcoxon p=0.263). Six LUAD sections (P1–P6) have no CD8A Gi* coldspots (|z|>1.96), so overlap is 0 there.
- Median cross-L at 2× NN = -12.3 (negative = fewer CD8A-high near CLDN4 Q4 than CSR).
- Tumor-domain CD8A (core − margin) median = -0.0469. Hop-depth Spearman of CD8A vs interior distance median = -0.0559 (Wilcoxon p=0.0382; 14/26 sections have epi-perm p<0.05).
- Section-level Spearman p-values treat spots as independent and are not used for inference; section-level Wilcoxon and epithelial CD8A permutations are the meta tests.

Wilcoxon p vs 0 on the raw Gi* overlap *count* is not an exclusion test (counts are ≥0). The relevant Gi* summary is overlap minus the independence expectation.

## How to read the exclusion-direction numbers

- Negative Spearman / bivariate Moran / Lee's L / residual lag ρ: lower CD8A where CLDN4 is higher.
- Positive Δ nearest-µm (Q4 − Q1): CLDN4-high epithelial-like spots are farther from CD8A-high spots than CLDN4-low.
- Negative Δ kNN CD8A (Q4 − Q1): fewer CD8A UMIs in the neighborhood of CLDN4-high epithelial-like spots.
- Positive Gi* overlap: CLDN4 hotspots coincide with CD8A coldspots.
- Negative cross-L at 2× NN: fewer CD8A-high points near CLDN4 Q4 points than a CSR / labeling null.
- Negative (core − margin) CD8A: lower CD8A in the CLDN4-high core than at the tumor-domain margin.

These are observed statistics and permutation p-values. They are not a clinical claim.

## Per-section LUAD table (primary set)

n=26 invasive LUAD sections. Full columns: `tables/luad_key_stats.csv` and `tables/section_stats.csv`.

| section | n | ρ CLDN4–CD8A (spot p) | epi-perm p | residual lag ρ (perm p) | Lee L (perm p) | Δnn µm Q4−Q1 (perm p) | ΔkNN CD8A | cross-L 2NN | Gi* overlap (obs−exp) | core−margin CD8A | hop-depth ρ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| P1_LUAD | 6603 | 0.12 (1.43e-22) | 0.94 | 0.076 (1) | 0.102 (0.625) | -36 (1) | 0.178 | 65 | 0 (0) | -0.0892 | 0.0447 |
| P2_LUAD | 10995 | 0.0383 (5.88e-05) | 0.98 | 0.149 (1) | -0.00182 (0.755) | -52 (1) | -0.0386 | 35 | 0 (0) | 0.028 | 0.0637 |
| P3_LUAD | 14336 | 0.0949 (4.84e-30) | 0.005 | -0.0287 (0.005) | 0.161 (0.045) | -11 (0.585) | -0.231 | 4.8 | 0 (0) | -0.0527 | -0.0611 |
| P4_LUAD | 12978 | -0.0105 (0.233) | 0.095 | 0.0654 (0.01) | -0.365 (0.005) | -0.53 (0.31) | -0.0885 | -9.1 | 0 (0) | 0.00216 | -0.0526 |
| P5_LUAD | 12084 | -0.0967 (1.73e-26) | 1 | 0.185 (1) | -0.758 (0.98) | -45 (1) | 0.259 | -16 | 0 (0) | 0.0234 | 0.0543 |
| P6_LUAD | 11724 | 0.0367 (7.15e-05) | 0.025 | -0.0184 (0.005) | -0.167 (0.005) | -10 (0.87) | -0.308 | 11 | 0 (0) | -0.0629 | -0.0662 |
| P7_LUAD | 13616 | 0.223 (1.74e-153) | 1 | 0.139 (1) | 12 (1) | -17 (1) | 0.702 | 41 | 32 (-118) | 0.13 | 0.15 |
| P7_LUAD-1 | 13074 | 0.0124 (0.155) | 1 | 0.07 (1) | 0.965 (1) | -0.22 (0.14) | 1.23 | -6.7 | 127 (-88) | 0.119 | -0.0919 |
| P8_LUAD | 12642 | -0.325 (8.66e-309) | 0.005 | -0.0474 (0.01) | -11.3 (0.005) | 302 (0.005) | -3.62 | -138 | 1145 (858) | -0.433 | -0.32 |
| P9_LUAD | 12127 | -0.0698 (1.4e-14) | 0.005 | -0.0712 (0.005) | -3.14 (0.005) | 57 (0.005) | -2.01 | -42 | 222 (53) | -0.337 | -0.38 |
| P10_LUAD | 13441 | -0.108 (3.58e-36) | 0.005 | -0.136 (0.005) | -0.634 (0.005) | 21 (0.005) | -0.259 | -16 | 713 (340) | -0.247 | -0.0679 |
| P11_LUAD | 11794 | 0.0671 (3.02e-13) | 1 | -0.0307 (0.005) | 1.28 (1) | -362 (1) | 4.37 | -39 | 82 (-240) | 0.304 | 0.239 |
| P12_LUAD | 14177 | -0.108 (8.49e-38) | 1 | -0.0398 (0.375) | -3.91 (1) | -40 (1) | 0.808 | -26 | 473 (126) | 0.0596 | 0.0644 |
| P13_LUAD | 4880 | 0.106 (1.46e-13) | 1 | -0.0136 (0.28) | 5.74 (1) | -57 (1) | 2.81 | 45 | 48 (-48) | 0.201 | 0.0111 |
| P14_LUAD | 11488 | 0.219 (5.79e-125) | 1 | 0.0386 (0.98) | 8.81 (0.985) | -13 (0.995) | 0.371 | 20 | 45 (-144) | -0.0745 | -0.0592 |
| P15_LUAD | 14231 | -0.152 (2.01e-74) | 0.005 | 0.03 (1) | -4.81 (0.005) | 91 (0.005) | -2.84 | -76 | 974 (694) | -0.359 | -0.237 |
| P16_LUAD | 7736 | -0.0198 (0.0817) | 0.005 | -0.0339 (0.005) | 0.0298 (0.005) | 5 (0.03) | -0.146 | 43 | 143 (-96) | -0.675 | -0.326 |
| P17_LUAD | 13515 | 0.0759 (9.76e-19) | 1 | -0.0552 (0.005) | 2.85 (1) | -76 (1) | 1.76 | -11 | 101 (-114) | -0.131 | -0.103 |
| P18_LUAD | 11971 | 0.09 (6.16e-23) | 0.03 | -0.00179 (0.53) | 1.97 (0.005) | 0.38 (0.485) | -0.206 | 20 | 87 (-72) | -0.142 | -0.157 |
| P19_LUAD | 12572 | -0.0817 (4.44e-20) | 1 | -0.0405 (0.97) | -2.01 (1) | -24 (1) | 0.645 | -14 | 531 (140) | -0.036 | 0.000291 |
| P20_LUAD | 13604 | -0.142 (2.88e-62) | 0.005 | -0.0504 (0.395) | -2.66 (0.005) | 28 (0.005) | -0.772 | -32 | 560 (305) | -0.0412 | -0.0475 |
| P21_LUAD | 13348 | -0.213 (1.05e-136) | 0.99 | -0.0426 (0.02) | -5.06 (0.495) | -29 (1) | 0.404 | -48 | 552 (244) | 0.0862 | 0.0448 |
| P22_LUAD | 12229 | 0.239 (2.34e-158) | 1 | 0.0505 (0.985) | 6.03 (1) | -44 (1) | 1.51 | 54 | 40 (-74) | 0.154 | 0.0641 |
| P23_LUAD | 10894 | -0.357 (<1e-300) | 0.995 | -0.0822 (0.15) | -7.34 (0.005) | 32 (0.005) | -0.304 | -71 | 698 (337) | -0.00829 | -0.0512 |
| P24_LUAD | 6276 | 0.166 (4.36e-40) | 0.005 | -0.0784 (0.02) | 1.45 (0.005) | 36 (0.005) | -0.984 | -30 | 112 (-95) | -0.455 | -0.36 |
| P25_LUAD | 3791 | -0.226 (3.29e-45) | 0.005 | -0.0632 (0.005) | -4.8 (0.005) | 23 (0.01) | -0.537 | -47 | 207 (76) | -0.31 | -0.205 |

P23_LUAD spot Spearman p underflowed to 0 in `scipy.stats.spearmanr`; shown as <1e-300.

## Figures

- `figures/summary_leeL_forest.png`, `figures/summary_biv_moran_forest.png`, `figures/summary_resid_lag_forest.png`
- `figures/summary_delta_nn_forest.png`, `figures/summary_gi_overlap.png`, `figures/summary_gi_overlap_minus_expected.png`
- `figures/summary_crossL_overlay.png`
- `figures/summary_by_histology_leeL.png`, `figures/summary_by_histology_delta_nn.png`, `figures/summary_cldn4_epi_by_histology.png`
- `figures/maps/*_cldn4_cd8a.png` (H&E-free)
- `figures/gi_star/*_gi_star.png`
- `figures/moran/*_moran_lag.png`
- `figures/kfunction/*_crossL.png`

## Software

numpy 2.4.4, pandas 3.0.5, scipy 1.18.0, matplotlib 3.11.1.

