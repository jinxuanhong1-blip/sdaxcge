# Donor-level spatial null: CosMx CLDN4 vs CD8+NK

Official NanoString CosMx NSCLC 960-plex (He et al. 2022), eight sections and five donors (Lung5 ×3, Lung6, Lung9 ×2, Lung12, Lung13). Additive to the locked radius and contact summaries. No private 8-KL data. Flat metadata has no author cell-type column, so tumor / CD8 / NK calls use the same RNA compartment score plus PanCK-high / CD45-low support as the radius script, and CD8+NK cells are required to be non-tumor.

## What is stronger than the sign-test floor

A one-sided sign test on eight sections cannot return a p-value below **0.003906** (1/256). On five donors the floor is **0.03125** (1/32). Wilcoxon on eight paired sections has the same two-sided floor, 1/128. Those tests use only the sign.

The p-value here is from a permutation of the CLDN4 mark. The statistic is the unweighted mean of five donor effects, so Lung5's three sections are one donor. Because the null has thousands of draws, the p-value can fall below both floors when the within-FOV pairing is extreme. The sign counts stay in the table as the reproducibility summary. They are not replaced.

This is not a p-value for drawing a new donor. The permutation null is: inside each FOV, tumor CLDN4 values are exchangeable with respect to a neighborhood that was computed once from fixed coordinates and fixed CD8/NK labels.

## How spatial labels are kept

- Neighbor counts and the 20 µm contact flag are computed once, inside the FOV, from centroid coordinates (`0.18 µm/px`).
- A permutation moves CLDN4 values among tumor cells that share a FOV. It does not move a value into another FOV, section, or donor, and it does not move a cell.
- Between-FOV geography is therefore left intact. A section-pooled high-versus-low gap that comes only from CLDN4-high cells sitting in immune-cold FOVs is invariant under this null and is not counted as evidence.
- The bootstrap resamples whole FOVs inside each section, then rebuilds the donor mean. Cells are not drawn independently of their neighbors.
- A second, coarser interval resamples the five donor effects. With five donors that interval cannot invent a p-value below 1/32, and it is not used for that claim.

## Pre-specified endpoints

All four use the within-FOV contrast, averaged with equal FOV weight inside a section, equal section weight inside a donor, and equal donor weight.

1. Median split of tumor CLDN4 inside the FOV: mean CD8+NK count at **50 µm**, high minus low.
2. Same split: CD8+NK **contact rate at 20 µm** (at least one CD8 or NK centroid within 20 µm), high minus low.
3. **Continuous** CLDN4: within-FOV Spearman versus the 50 µm CD8+NK count.
4. **Continuous** CLDN4: within-FOV Spearman versus the 20 µm contact indicator.

Exclusion is the negative direction. One-sided permutation p = (1 + number of null draws at least as small as observed) / (1 + n_perm). Radii 20 / 40 / 80 / 100 µm, the CD8-only sensitivity, the cell-weighted FOV mean, and the section-pooled median are secondary.

## Inventory

QC cells 766,334. Tumor 397,006. CD8 45,326. NK 44,287. CD8+NK 83,239. FOVs used 232/232 (≥20 tumor and ≥5 CD8+NK). Tumor cells inside those FOVs 397,006.

| Section | Donor | QC cells | Tumor | CD8 | NK | CD8+NK | FOVs used |
|---|---|---:|---:|---:|---:|---:|---:|
| Lung5_Rep1 | Lung5 | 98,046 | 40,547 | 7,010 | 7,222 | 13,313 | 30/30 |
| Lung5_Rep2 | Lung5 | 100,437 | 40,606 | 9,771 | 8,806 | 17,179 | 29/29 |
| Lung5_Rep3 | Lung5 | 97,876 | 39,185 | 6,133 | 7,253 | 12,571 | 30/30 |
| Lung6 | Lung6 | 90,126 | 73,852 | 1,299 | 1,390 | 2,541 | 30/30 |
| Lung9_Rep1 | Lung9 | 87,617 | 46,610 | 4,049 | 4,407 | 7,836 | 20/20 |
| Lung9_Rep2 | Lung9 | 139,544 | 92,862 | 4,985 | 4,628 | 8,840 | 45/45 |
| Lung12 | Lung12 | 71,440 | 30,124 | 7,097 | 5,844 | 11,860 | 28/28 |
| Lung13 | Lung13 | 81,248 | 33,220 | 4,982 | 4,737 | 9,099 | 20/20 |

## Primary results

| Endpoint | Donor-mean T | One-sided perm p | Two-sided perm p | FOV-block 95% CI | Sections neg | Donors neg | Below 1/256 | Below 1/32 |
|---|---:|---:|---:|---:|---:|---:|---|---|
| Median, 50 µm count | -0.976 | 2.00e-04 | 2.00e-04 | [-1.083, -0.873] | 7/8 | 4/5 | yes | yes |
| Median, 20 µm contact | -0.129 | 2.00e-04 | 2.00e-04 | [-0.138, -0.120] | 8/8 | 5/5 | yes | yes |
| Spearman, 50 µm count | -0.127 | 2.00e-04 | 2.00e-04 | [-0.139, -0.116] | 7/8 | 4/5 | yes | yes |
| Spearman, 20 µm contact | -0.126 | 2.00e-04 | 2.00e-04 | [-0.136, -0.117] | 8/8 | 5/5 | yes | yes |

Permutation draws: 4999. FOV-block bootstrap draws: 1999. Seed 25976224. The smallest one-sided permutation p this run can return is 1/5000.

At 50 µm the within-FOV sign counts are 7/8 sections (sign p = 0.0352) and 4/5 donors (sign p = 0.1875). Lung6, the immune-poor donor, is the exception at 50 µm (0.038; interval [-0.012, 0.102]). The permutation p is 2.00e-04 because the donor-mean magnitude is extreme: the most negative null draw was -0.067, against an observed -0.976. The within-FOV null mean is 0.0000.

The 20 µm contact contrast is negative in 8/8 sections (sign p = 0.0039) and 5/5 donors (sign p = 0.0312). Those sign tests sit on the 8-section floor (0.003906) and the 5-donor floor (0.03125). The permutation p is lower than both because it uses the size of the deficit.

## Radius and contact family

Within-FOV median contrast (high − low), equal-donor mean. Negative means CLDN4-high tumor has fewer CD8+NK neighbors or contacts.

| Outcome | T | perm p (one-sided) | FOV-block 95% CI | Section signs | Donor signs | Cell-weighted T | Cell-weighted p |
|---|---:|---:|---:|---:|---:|---:|---:|
| CD8+NK count, 20 µm | -0.237 | 2.00e-04 | [-0.259, -0.217] | 8/8 | 5/5 | -0.235 | 2.00e-04 |
| CD8+NK count, 40 µm | -0.717 | 2.00e-04 | [-0.789, -0.646] | 7/8 | 4/5 | -0.717 | 2.00e-04 |
| CD8+NK count, 50 µm | -0.976 | 2.00e-04 | [-1.083, -0.873] | 7/8 | 4/5 | -0.975 | 2.00e-04 |
| CD8+NK count, 80 µm | -1.684 | 2.00e-04 | [-1.920, -1.467] | 7/8 | 4/5 | -1.692 | 2.00e-04 |
| CD8+NK count, 100 µm | -2.044 | 2.00e-04 | [-2.350, -1.732] | 7/8 | 4/5 | -2.069 | 2.00e-04 |
| CD8+NK contact, 20 µm | -0.129 | 2.00e-04 | [-0.138, -0.120] | 8/8 | 5/5 | -0.129 | 2.00e-04 |
| CD8 count, 50 µm | -0.536 | 2.00e-04 | [-0.604, -0.467] | 7/8 | 4/5 | -0.536 | 2.00e-04 |
| CD8 contact, 20 µm | -0.086 | 2.00e-04 | [-0.094, -0.078] | 8/8 | 5/5 | -0.086 | 2.00e-04 |

## Continuous CLDN4

Within-FOV Spearman, then the same equal-weight average. A negative value means higher CLDN4 tracks fewer CD8+NK neighbors. This Spearman uses every tumor cell, including CLDN4-undetected cells. It is not a second median split.

| Outcome | Mean Spearman | perm p (one-sided) | FOV-block 95% CI | Section signs | Donor signs |
|---|---:|---:|---:|---:|---:|
| CD8+NK count, 20 µm | -0.130 | 2.00e-04 | [-0.139, -0.120] | 8/8 | 5/5 |
| CD8+NK count, 40 µm | -0.136 | 2.00e-04 | [-0.148, -0.124] | 7/8 | 4/5 |
| CD8+NK count, 50 µm | -0.127 | 2.00e-04 | [-0.139, -0.116] | 7/8 | 4/5 |
| CD8+NK count, 80 µm | -0.101 | 2.00e-04 | [-0.112, -0.089] | 7/8 | 4/5 |
| CD8+NK count, 100 µm | -0.083 | 2.00e-04 | [-0.096, -0.071] | 7/8 | 4/5 |
| CD8+NK contact, 20 µm | -0.126 | 2.00e-04 | [-0.136, -0.117] | 8/8 | 5/5 |

In 149 of 232 FOVs the tumor CLDN4 median is 0, so the within-FOV median split is detected versus undetected in those fields. About 60% of tumor cells in a typical FOV have no CLDN4 count. A five-bin cut of all cells is not a dose curve under that tie: the zeros share one rank. The figure instead shows the undetected cells and, among CLDN4 > 0, four equal-count quartiles (P1 low to P4 high).

The pre-specified Spearman above includes the undetected cells. A separate check, run after seeing how often the median is zero, restricts to CLDN4 > 0 and repeats the same donor-level permutation. Among expressors the association does **not** continue toward exclusion. The donor-mean contrasts are positive: higher CLDN4 among detected cells goes with more CD8+NK neighbors, not fewer. The one-sided exclusion p-value is 1. Lung13 is the donor that still goes the exclusion way; Lung5, Lung6, and Lung12 go the other way.

| Among CLDN4 > 0 | Donor-mean T | Exclusion perm p | FOV-block 95% CI | Donors negative |
|---|---:|---:|---:|---:|
| Spearman, 50 µm count | 0.054 | 1.0000 | [0.038, 0.072] | 2/5 |
| Spearman, 20 µm contact | 0.064 | 1.0000 | [0.051, 0.076] | 1/5 |
| Median split, 50 µm count | 0.317 | 1.0000 | [0.202, 0.419] | 2/5 |
| Median split, 20 µm contact | 0.039 | 1.0000 | [0.030, 0.050] | 2/5 |

Positive-only permutation draws: 1999. FOVs with at least 30 CLDN4-positive tumor cells: 230.

## Section-pooled median, same spatial null

This contrast pools tumor cells inside a section and splits on the section median. It can be large when CLDN4-high cells occupy different FOVs from CLDN4-low cells. The permutation still shuffles CLDN4 only inside FOVs, so that between-FOV arrangement stays in the null and the section-pooled null mean is not zero. At 50 µm that null mean is -0.324 and the observed donor-mean is -1.298. The excess, -0.974, matches the within-FOV donor-mean of -0.976.

| Outcome | Section-pooled T | perm p (one-sided) | FOV-block 95% CI | Section signs | Donor signs | Mean high | Mean low |
|---|---:|---:|---:|---:|---:|---:|---:|
| CD8+NK count, 20 µm | -0.289 | 2.00e-04 | [-0.313, -0.263] | 8/8 | 5/5 | 0.236 | 0.525 |
| CD8+NK count, 40 µm | -0.932 | 2.00e-04 | [-1.022, -0.830] | 8/8 | 5/5 | 1.291 | 2.222 |
| CD8+NK count, 50 µm | -1.298 | 2.00e-04 | [-1.432, -1.146] | 8/8 | 5/5 | 2.190 | 3.488 |
| CD8+NK count, 80 µm | -2.418 | 2.00e-04 | [-2.739, -2.060] | 8/8 | 5/5 | 6.459 | 8.877 |
| CD8+NK count, 100 µm | -3.125 | 2.00e-04 | [-3.605, -2.603] | 8/8 | 5/5 | 10.588 | 13.714 |
| CD8+NK contact, 20 µm | -0.160 | 2.00e-04 | [-0.170, -0.147] | 8/8 | 5/5 | 0.172 | 0.331 |

The locked public summary (50 / 100 µm cytotoxic ratio 0.36 / 0.52, 8/8 and 5/5, sign P = 0.031) is a different estimand and is not recomputed here. In this run the ratio of the section-pooled donor-mean counts is 2.190/3.488 = 0.63 at 50 µm and 10.588/13.714 = 0.77 at 100 µm. Those ratios are milder than 0.36 / 0.52. The sign pattern for the section-pooled contrast is still 8/8 sections and 5/5 donors. This does not overwrite the locked summary.

## Donor effects for the primary median contrasts

| Donor | 50 µm count Δ | 50 µm block CI | 20 µm contact Δ | Contact block CI |
|---|---:|---:|---:|---:|
| Lung5 | -1.449 | [-1.640, -1.257] | -0.183 | [-0.202, -0.164] |
| Lung6 | 0.038 | [-0.012, 0.102] | -0.006 | [-0.013, 0.003] |
| Lung9 | -0.853 | [-1.084, -0.651] | -0.098 | [-0.119, -0.078] |
| Lung12 | -1.710 | [-2.117, -1.325] | -0.203 | [-0.233, -0.176] |
| Lung13 | -0.904 | [-1.097, -0.727] | -0.152 | [-0.177, -0.129] |

## Reading rule

A permutation p below 1/256 means the donor-averaged within-FOV contrast is extreme under FOV-restricted label exchange. It does not mean the five-donor sign test has been given a smaller discrete floor. If the section-pooled gap is large and the within-FOV permutation p is not, the gap is carried by which FOV a cell sits in, and this null correctly refuses to call that a within-neighborhood CLDN4 effect.

## Methods notes

- CLDN4 is log-normalized inside each section: size factor = median total of the marker panel used for typing, divided by that cell's total, then log1p. The QC filter (≥20 counts and ≥5 genes) uses the full 960-plex total, excluding NegPrb columns.
- Tumor: epithelial RNA score above immune and stromal scores, or PanCK at/above the section median and CD45 below it, and not called immune. CD8: CD8A or CD8B count > 0 and not tumor. NK: NKG7 or GNLY count > 0, no CD3, and not tumor. KLRD1 is not on this 960-plex. CD8+NK is the union. Tumor cells with CD8A stay in the tumor index.
- Median split inside a FOV: CLDN4 strictly above the FOV median versus at or below it. A FOV needs at least 5 cells on each side to enter the median contrast. Spearman does not use that split.
- Contact radius 20 µm. Count radii 20, 40, 50, 80, 100 µm. Neighbors are counted inside the FOV only.
- Cell-weighted sensitivity weights FOVs by tumor-cell count, then still averages sections and donors equally.

## Figures

- `results/cosmx_donor_spatial_null/figures/null_primary.png`
- `results/cosmx_donor_spatial_null/figures/donor_forest.png`
- `results/cosmx_donor_spatial_null/figures/cldn4_quintiles.png`
- `results/cosmx_donor_spatial_null/figures/radius_profile.png`

