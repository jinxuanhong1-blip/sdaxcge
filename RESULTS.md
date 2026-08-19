# RESULTS: CosMx NSCLC — CLDN4-high tumor vs CD8+NK neighbors (all 8 sections)

## Paper sentence

In the official CosMx NSCLC FFPE 960-plex cohort (all 8 sections / 5 donors; 302,313 author-typed tumor cells, 15,955 CD8 T cells + 7,814 NK = 23,769 CD8+NK; 765,771 cells, 232 FOVs), CD8+NK neighbor counts at 40 µm around CLDN4-high vs CLDN4-low tumor cells were heterogeneous and not consistently reduced (section-paired median ∆[high−low] = +0.010 cells, two-sided Wilcoxon p = 0.945, one-sided high<low p = 0.578; 3/8 sections with fewer neighbors around high, 5/8 with more; donor-paired mean ∆ = -0.019, two-sided p = 0.812, n = 5 donors). Mixing-score ∆med = +0.0064 (two-sided p = 0.383); density-normalized CD8+NK fraction ∆med = +0.0005 (two-sided p = 0.742). Tumor cells were not gated on CD8A==0.

## Dataset and panel

- **Dataset:** official NanoString/Bruker CosMx SMI NSCLC FFPE (8 FFPE sections / 5 patients, 960-plex). https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/
- **Mirror used:** Zenodo 15487520 `cosmx_lung` (counts, coordinates, author `cell_type`).
- **Panel check:** CLDN4 is on the 960-gene matrix. Genes present: CLDN4, CD8A, CD8B, KRT8, EPCAM, NKG7.
- **Pixel size:** 0.18 µm/pixel (official SMI-ReadMe).
- **Tumor:** author types `tumor 5/6/9/12/13`. Tumor cells with CD8A>0 were **kept** (no CD8A==0 gate).
- **CD8+NK:** author `T CD8 naive` + `T CD8 memory` + `NK`.
- **CLDN4-high/low:** per-section median split of residual log1p(CLDN4) after OLS on log1p(KRT8), log1p(EPCAM), log1p(n_counts).
- **Units of inference:** section (n=8, paired) and donor (n=5, paired). FOV cells are not the test unit.
- No private 8-KL data.

- Loaded **765,771 cells** (5,465 official cells dropped for non-finite coordinates), **232 FOVs**, **all 8/8 whole-section CosMx samples** (not a TMA subset), **5 donors**.
- Author tumor = 302,313 (high 151,121 / low 151,192); CD8 = 15,955; NK = 7,814; CD8+NK = 23,769.
- Fraction of tumor cells with CD8A>0 (kept): 0.086.

## Primary: CD8+NK neighbor counts at 20 / 40 / 60 µm

Mean CD8+NK cells within radius r of each tumor cell, then averaged within section. Paired Wilcoxon signed-rank on the 8 section means (high vs low) and on the 5 donor means (section means averaged within donor).

| r | CLDN4-low mean | CLDN4-high mean | ∆med (high−low) | section p (high<low) | donor ∆mean | donor p | sections with high<low |
|---|---:|---:|---:|---:|---:|---:|---:|
| 20 µm | 0.058 | 0.054 | 0.002 | 0.578 | -0.007 | 0.500 | 3/8 |
| 40 µm | 0.293 | 0.286 | 0.010 | 0.578 | -0.019 | 0.406 | 3/8 |
| 60 µm | 0.748 | 0.748 | 0.028 | 0.578 | -0.028 | 0.500 | 4/8 |


### Per-section 40 µm (signed)

| Section | Donor | n tumor | n CD8+NK | low count | high count | ∆ (high−low) | low mix | high mix | perm p (∆count high<low) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Lung5_Rep1 | Lung5 | 18,600 | 3,688 | 0.111 | 0.134 | +0.023 | 0.014 | 0.024 | 1.000 |
| Lung5_Rep2 | Lung5 | 19,110 | 3,596 | 0.125 | 0.142 | +0.018 | 0.021 | 0.026 | 1.000 |
| Lung5_Rep3 | Lung5 | 16,555 | 3,908 | 0.113 | 0.159 | +0.046 | 0.015 | 0.032 | 1.000 |
| Lung6 | Lung6 | 66,981 | 1,414 | 0.056 | 0.105 | +0.049 | 0.003 | 0.011 | 1.000 |
| Lung9_Rep1 | Lung9 | 39,422 | 3,069 | 0.318 | 0.320 | +0.002 | 0.019 | 0.027 | 0.595 |
| Lung9_Rep2 | Lung9 | 96,822 | 2,861 | 0.283 | 0.248 | -0.035 | 0.020 | 0.020 | 0.005 |
| Lung12 | Lung12 | 18,512 | 1,877 | 0.299 | 0.255 | -0.044 | 0.059 | 0.051 | 0.005 |
| Lung13 | Lung13 | 26,311 | 3,356 | 1.039 | 0.927 | -0.112 | 0.066 | 0.056 | 0.005 |

Lung5 (3 sections) and Lung6 have **more** CD8+NK around CLDN4-high tumor; Lung12 and Lung13 have **fewer** (section label-permutation p = 0.005). The cohort-level paired test is therefore a mix of opposite donor directions, not a uniform exclusion.

Per-section values: `tables/section_neighbor_stats.csv`. Per-donor: `tables/donor_neighbor_stats.csv`.

Section-level high/low label permutation (n=199) for ∆count at 40 µm: median per-section p = 0.797; sections with perm p<0.05: 3/8.

## Mixing score and density-normalized CD8+NK fraction

- **Mixing score** at r: among neighbors that are tumor or CD8+NK, the fraction that are CD8+NK (immune–tumor mixing vs tumor–tumor self-aggregation).
- **Density-normalized fraction:** CD8+NK / all cells within r (accounts for local packing).
- **Enrichment:** density-normalized fraction / section CD8+NK prevalence.

- Mixing at 40 µm: low 0.0272 vs high 0.0310; ∆med = +0.0064; one-sided high<low p = 0.844; two-sided p = 0.383.
- Density-normalized CD8+NK fraction at 40 µm: low 0.0103 vs high 0.0105; ∆med = +0.0005; one-sided high<low p = 0.680; two-sided p = 0.742.
- Enrichment vs section prevalence at 40 µm: low 0.326 vs high 0.338.

## 1. CLDN4 in tumor vs other types

- Mean raw CLDN4: tumor 1.444 vs other 0.146.
- log1p CLDN4 MWU (tumor > other): p < 1e-300 (n_tumor=302,313, n_other=463,458).

## Secondary: nearest-µm (not primary)

Nearest CD8+NK distance is reported only as a negative control on why it is the wrong lead: section median 86.4 vs 86.6 µm (high−low ∆med = -0.62 µm; p = 0.727). Counts and mixing at 20/40/60 µm are the primary spatial readouts.

## Secondary: Ripley's K and pair-correlation g(r)

Bivariate K̂ and ĝ from CLDN4-high or low tumor cells to CD8+NK, border-corrected, averaged within section then across sections. CSR reference is πr² (K) or g=1.
- Mean K̂(40) high vs low: 2420.0 vs 2552.1 (CSR π·40² = 5026.5).
- FOVs contributing to K/g: sum of per-section Ripley FOVs = 216.

## Counts

| Level | n |
|---|---|
| Donors | 5 |
| Sections | 8 / 8 |
| FOVs | 232 |
| Cells | 765,771 |
| Author tumor | 302,313 |
| CLDN4-high / low tumor | 151,121 / 151,192 |
| Author CD8 | 15,955 |
| Author NK | 7,814 |
| CD8+NK | 23,769 |

## Figures

- `figures/fig01_maps.png` — example FOV maps from the 8-section cohort.
- `figures/fig02_neighbor_counts_20_40_60.png` — **primary** paired section CD8+NK counts.
- `figures/fig03_mixing_and_fraction.png` — mixing score and density-normalized fraction at 40 µm.
- `figures/fig04_donor_paired_40um.png` — donor-paired 40 µm counts.
- `figures/fig05_count_ecdf_40um.png` — tumor-cell count ECDF (supporting).
- `figures/fig06_cldn4_tumor_vs_other.png` — CLDN4 expression.
- `figures/fig07_nn_distance_secondary.png` — nearest-µm (not primary).
- `figures/fig08_K_hat.png` — Ripley's K (secondary).
- `figures/fig09_g_r.png` — pair-correlation g(r) (secondary).

## Methods

All 8 official CosMx NSCLC sections were used (Lung5 Rep1–3, Lung6, Lung9 Rep1–2, Lung12, Lung13). Coordinates were converted with 0.18 µm/pixel. Author cell types define tumor and CD8+NK; tumor cells were not excluded for CD8A expression. CLDN4-high vs low is a within-section median split of residual log1p(CLDN4) after KRT8, EPCAM, and library size. For each tumor cell, CD8+NK neighbors, all-cell neighbors, and tumor neighbors were counted at 20, 40, and 60 µm with a KD-tree. Mixing = n_CD8NK / (n_CD8NK + n_tumor_neighbors). Density-normalized fraction = n_CD8NK / n_all_neighbors. Section means were tested with paired Wilcoxon (high vs low). Donor means are unweighted averages of that donor's sections. High/low labels were also shuffled within section (199 times) to get a section-level permutation p for ∆count and ∆mix at 40 µm. Ripley's K and g(r) were computed per FOV and averaged; they are secondary to neighbor counts.

## Notes

- Additive public CLDN4-only spatial analysis. No private 8-KL.
- Nearest-µm is retained in the tables but is not the lead metric.
