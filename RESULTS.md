# RESULTS: cross-type Ripley K / L / g(r) for CLDN4-high tumor vs CD8 T cells

Official CosMx NSCLC FFPE 960-plex (He et al., *Nat Biotechnol* 2022; 8 samples / 5 patients).
This is **not** a Spearman-only proximity screen. All primary numbers below are from
cross-type point-process summaries with permutation envelopes.

## Data and panel check

- Source: NanoString/Bruker public CosMx NSCLC FFPE release (`All SMI Giotto object.tar.gz` from `nanostring-public-share`, 2021-10-28), the processed object that carries the **author cell-type calls** used in the CosMx data viewer / He et al. pipeline.
- Cells after author QC in that object: **771,236**.
- Samples (8): Lung12, Lung13, Lung5_Rep1, Lung5_Rep2, Lung5_Rep3, Lung6, Lung9_Rep1, Lung9_Rep2.
- Patients (5): Lung12, Lung13, Lung5, Lung6, Lung9.
- FOVs total: **233**. FOVs analysed (passing count filters): **159**. Skipped: **74**.
- **CLDN4 on panel: True** (960 RNA targets in `data/cosmx_960_genes.csv`; CLDN4 is present in the official `exprMat` gene list from SMI-ReadMe.html).
- Pixel size stated in the official SMI ReadMe: **0.18 µm/pixel**. Coordinates in the Giotto object are millimetres (FOV x-span ≈ 0.98 mm); analysis uses micrometres.
- No private KL / KP / 8-KL mouse data were used.

### Author cell types used

Tumor cells: all labels beginning with `tumor` (`tumor 5`, `tumor 6`, `tumor 9`, `tumor 12`, `tumor 13`).
CD8 T cells: `T CD8 memory` + `T CD8 naive`.
Epithelial (for λ only): `epithelial` plus all tumor labels.

Full author type counts:

```
cell_type
tumor 9         134099
fibroblast       93413
neutrophil       67726
tumor 6          66805
macrophage       53195
tumor 5          52286
plasmablast      38575
endothelial      38326
T CD4 naive      29327
B-cell           27512
tumor 13         26757
epithelial       24284
tumor 12         22442
T CD4 memory     16850
mDC              16796
pDC              13284
mast             11264
T CD8 naive      10565
monocyte          8296
NK                7875
Treg              6054
T CD8 memory      5505
```

### CLDN4-high definition (pre-specified)

Among author tumor cells, **CLDN4-high** = CLDN4 count ≥ sample-specific 75th percentile **and** CLDN4 ≥ 1.
Lung6 tumor CLDN4 75th percentile is 0, so CLDN4-high there reduces to CLDN4 ≥ 1 (detected).

| sample     |   n_tumor |   cldn4_mean |   cldn4_q75 |   n_high |
|:-----------|----------:|-------------:|------------:|---------:|
| Lung12     |     18512 |     0.860199 |           1 |     7288 |
| Lung13     |     26311 |     2.34085  |           3 |     9392 |
| Lung5_Rep1 |     18600 |     2.13785  |           3 |     5703 |
| Lung5_Rep2 |     19186 |     1.82357  |           3 |     5017 |
| Lung5_Rep3 |     16555 |     1.77868  |           2 |     6426 |
| Lung6      |     66981 |     0.373449 |           0 |    16402 |
| Lung9_Rep1 |     39422 |     2.16115  |           3 |    13001 |
| Lung9_Rep2 |     96822 |     1.49573  |           2 |    34641 |

CD8 T cells per sample:

| sample     |   n_cd8 |
|:-----------|--------:|
| Lung12     |    1310 |
| Lung13     |    3251 |
| Lung5_Rep1 |    2020 |
| Lung5_Rep2 |    1963 |
| Lung5_Rep3 |    2294 |
| Lung6      |    1230 |
| Lung9_Rep1 |    2323 |
| Lung9_Rep2 |    1679 |

## Estimators

- **Homogeneous cross-K / L / g(r)** between CLDN4-high tumor points and CD8 points, border (reduced-sample) edge correction, rectangular FOV window.
- **Inhomogeneous K / g(r)** using λ from an 80 µm Gaussian-smoothed intensity of **all epithelial/tumor cells** evaluated at CLDN4-high locations, and CD8 intensity at CD8 locations. This asks whether CD8 are depleted around CLDN4-high tumor cells *after accounting for tumor/epithelial density*.
- **Null 1 (primary): 399 label permutations** — randomly re-choose the same number of “high” cells among tumor cells in the FOV; CD8 positions fixed. This is the test of CLDN4-specificity vs “any tumor is dense”.
- **Null 2: 199 CSR envelopes** — CD8 relocated uniformly in the FOV rectangle; CLDN4-high fixed.
- Theoretical homogeneous CSR reference: K(r) = πr², g(r) = 1.
- FOV filters: ≥20 CLDN4-high tumor, ≥20 CD8, ≥40 tumor cells.
- r grid: 10–150 µm (5 µm). g(r) uses 5 µm rings. r* search window: 15–120 µm.

## Where exclusion is strongest

Two nested questions are separated:

1. **CSR / homogeneous g(r) vs πr²:** CD8 vs CLDN4-high tumor. g << 1 is the expected tumor-vs-stroma geometry (CD8 sit outside tumor nests). This is **not** a CLDN4-specific claim.
2. **Label permutation + inhomogeneous λ(tumor/epithelial):** does CD8 avoid *CLDN4-high* tumor cells more than a random tumor subset of the same size? That is Δg_inhom = g_obs − permutation median. Negative Δg is extra exclusion.

- **Meta r\*** (radius minimising mean Δg_inhom across FOVs, search 15–120 µm): **22.5 µm**.
- Mean g_inhom(r\*) across FOVs: **0.091** (SEM 0.008). This quantity is << 1 even under the label null because λ is tumor/epithelial density, not the CLDN4-high process intensity; **do not treat g_inhom < 1 as the CLDN4 test**.
- Mean Δg_inhom(r\*): **-0.043** (Wilcoxon signed-rank one-sided p = **1.102e-12**).
- Stouffer combined one-sided p (label permutation, inhomogeneous g at r\*): **1.979e-10**.
- FOVs with g_inhom(r\*) < 1 (CSR-like, not CLDN4-specific): **159/159**.
- FOVs with label-permutation p < 0.05 at r\*: **40/159**.
- FOVs with BH-FDR q < 0.05 on that p: **17/159**.
- Median per-FOV r of minimum *raw* g_inhom: **12.5 µm** (the 10–15 µm ring is dominated by cell-body exclusion; r* above is the meta Δg minimum).

The CLDN4-specific statement is therefore: at **22.5 µm**, CD8–CLDN4-high pairs are fewer than expected from random labeling of tumor cells (mean Δg_inhom < 0; combined p as above). The effect is **not uniform**: it is strongest in Lung9 and Lung12 and is **not detected** in the Lung5 serial sections (median p ≈ 0.5–0.7). Lung5_Rep3’s mean Δg is pulled by one sparse-ring FOV; its median p remains 0.53.

A CSR calibration of the homogeneous estimator on simulated uniform points recovered K/(πr²) ≈ 0.98 and g ≈ 0.99 (`scripts/test_csr_calibration.py`).

Interpretation is restricted to the envelopes. Causal barrier function of CLDN4 is not claimed.

## Per-FOV table (r*, g, p)

Full per-FOV numbers: `results/tables/per_fov_g_inhom.csv`.

FOVs with the most negative Δg_inhom(r*) among well-powered FOVs (n_high≥80 and n_CD8≥50):

| sample     | patient   |   fov |   n_high |   n_cd8 |   r_star_um |   g_inhom_rstar |   delta_g_inhom_rstar |   p_label_inhom |   p_label_hom |   p_csr_hom |   r_gmin_um |   g_inhom_min |
|:-----------|:----------|------:|---------:|--------:|------------:|----------------:|----------------------:|----------------:|--------------:|------------:|------------:|--------------:|
| Lung12     | Lung12    |    23 |      120 |     121 |        22.5 |       0.0917216 |            -0.199514  |          0.015  |        0.0025 |       0.005 |        12.5 |     0         |
| Lung12     | Lung12    |    27 |      198 |     144 |        22.5 |       0.0301671 |            -0.131731  |          0.01   |        0.01   |       0.005 |        12.5 |     0         |
| Lung12     | Lung12    |    24 |      210 |      69 |        22.5 |       0         |            -0.130636  |          0.0425 |        0.0425 |       0.005 |        12.5 |     0         |
| Lung12     | Lung12    |     8 |      438 |      58 |        22.5 |       0.025307  |            -0.129422  |          0.0025 |        0.005  |       0.005 |        22.5 |     0.025307  |
| Lung12     | Lung12    |    25 |      299 |      53 |        22.5 |       0.238431  |            -0.123677  |          0.1175 |        0.085  |       0.05  |        12.5 |     0.108584  |
| Lung9_Rep2 | Lung9     |    44 |      221 |      92 |        22.5 |       0.0104992 |            -0.109856  |          0.0025 |        0.0025 |       0.005 |        12.5 |     0         |
| Lung13     | Lung13    |     4 |      480 |      67 |        22.5 |       0.156297  |            -0.101452  |          0.0175 |        0.1825 |       0.12  |        12.5 |     0.117705  |
| Lung9_Rep2 | Lung9     |    38 |      195 |      78 |        22.5 |       0.0240149 |            -0.0976584 |          0.01   |        0.01   |       0.005 |        12.5 |     0         |
| Lung9_Rep2 | Lung9     |    39 |      711 |      53 |        22.5 |       0.106532  |            -0.0964016 |          0.0025 |        0.0075 |       0.005 |        12.5 |     0.0669766 |
| Lung9_Rep2 | Lung9     |     3 |     1189 |      52 |        22.5 |       0.123814  |            -0.0955816 |          0.01   |        0.0025 |       0.005 |        12.5 |     0.0339025 |
| Lung13     | Lung13    |     6 |      613 |     154 |        22.5 |       0.150791  |            -0.0902404 |          0.0025 |        0.045  |       0.005 |        12.5 |     0.131829  |
| Lung13     | Lung13    |     9 |      522 |     179 |        22.5 |       0.205674  |            -0.0877112 |          0.0075 |        0.01   |       0.005 |        22.5 |     0.205674  |

FOVs with the largest (least exclusive) Δg_inhom(r*) in that same subset:

| sample     | patient   |   fov |   n_high |   n_cd8 |   r_star_um |   g_inhom_rstar |   delta_g_inhom_rstar |   p_label_inhom |   p_label_hom |   p_csr_hom |   r_gmin_um |   g_inhom_min |
|:-----------|:----------|------:|---------:|--------:|------------:|----------------:|----------------------:|----------------:|--------------:|------------:|------------:|--------------:|
| Lung13     | Lung13    |     3 |      420 |     122 |        22.5 |       0.397349  |            0.0361548  |          0.7425 |        0.3675 |       0.27  |        17.5 |     0.240136  |
| Lung13     | Lung13    |     8 |      266 |     213 |        22.5 |       0.26482   |            0.0304498  |          0.775  |        0.78   |       0.11  |        12.5 |     0.155969  |
| Lung13     | Lung13    |     1 |      377 |      85 |        22.5 |       0.420597  |            0.0243545  |          0.6325 |        0.5525 |       0.365 |        52.5 |     0.199447  |
| Lung6      | Lung6     |    30 |       90 |      92 |        22.5 |       0.0159068 |            0.0159068  |          0.895  |        0.9625 |       0.05  |        12.5 |     0         |
| Lung13     | Lung13    |    19 |      635 |     120 |        22.5 |       0.33955   |            0.010455   |          0.5925 |        0.7775 |       0.115 |       147.5 |     0.251192  |
| Lung5_Rep1 | Lung5     |    12 |      235 |      61 |        22.5 |       0.0358815 |            0.0102588  |          0.675  |        0.49   |       0.005 |        12.5 |     0         |
| Lung9_Rep1 | Lung9     |    16 |      593 |      67 |        22.5 |       0.119395  |            0.00999103 |          0.65   |        0.1675 |       0.005 |        12.5 |     0.0364291 |
| Lung5_Rep3 | Lung5     |     7 |      528 |      73 |        22.5 |       0.0536587 |            0.00784781 |          0.5125 |        0.78   |       0.005 |        57.5 |     0.0169248 |

## Sample-level summary of g_inhom(r*)

| sample     |   n_fov |    mean_g |   median_g |   mean_delta_g |   frac_g_lt1 |   frac_p05 |   median_p |
|:-----------|--------:|----------:|-----------:|---------------:|-------------:|-----------:|-----------:|
| Lung12     |      25 | 0.103795  |  0.0917216 |   -0.080454    |            1 |     0.32   |    0.0925  |
| Lung13     |      20 | 0.266556  |  0.270341  |   -0.0256486   |            1 |     0.2    |    0.32625 |
| Lung5_Rep1 |      15 | 0.0163875 |  0         |   -0.000937563 |            1 |     0      |    0.675   |
| Lung5_Rep2 |      17 | 0.0274043 |  0         |    0.0130462   |            1 |     0      |    0.5325  |
| Lung5_Rep3 |      16 | 0.0247899 |  0.0102501 |   -0.107283    |            1 |     0.0625 |    0.52625 |
| Lung6      |      16 | 0.0563435 |  0.0395769 |    0.00267952  |            1 |     0      |    0.43    |
| Lung9_Rep1 |      18 | 0.0750274 |  0.0704343 |   -0.0402238   |            1 |     0.5    |    0.05125 |
| Lung9_Rep2 |      32 | 0.100825  |  0.0954826 |   -0.063884    |            1 |     0.5625 |    0.035   |

## Figures

- `results/figures/meta_g_inhom_across_fovs.png` — mean Δg_inhom(r) ± 1.96 SEM and FOV-level p/envelope fractions.
- `results/figures/meta_g_inhom_raw_mean.png` — mean raw g_inhom(r) (not the label-permutation null).
- `results/figures/g_inhom_rstar_by_sample.png` — per-FOV g_inhom(r*) by sample.
- `results/figures/kg_envelope_strongest_exclusion_Lung12_fov23.png` — K, L−r, and g vs 399 label-permutation envelopes (well-powered FOV with most negative Δg).
- `results/figures/kg_envelope_median_fov_Lung9_Rep1_fov8.png` — same for a median well-powered FOV.
- `results/figures/kg_envelope_weakest_exclusion_Lung13_fov3.png` — same for the least exclusive well-powered FOV.
- `results/figures/g_nulls_csr_vs_label_Lung12_fov23.png` — CSR vs label-permutation envelopes on homogeneous g(r).
- `results/figures/spatial_strongest_exclusion_Lung12_fov23.png` and `spatial_median_fov_Lung9_Rep1_fov8.png` — cell maps (CLDN4-high tumor vs CD8).

Curve objects (every FOV, every r) are in `results/tables/fov_curves.jsonl`.

## What this does *not* claim

- Causal barrier function of CLDN4.
- Results from any private 8-KL / KP mouse cohort.
- A Spearman correlation as the primary spatial test (none is reported as the headline).

