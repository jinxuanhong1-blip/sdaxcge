# CLDN4 tumor-geography: CD8 infiltration depth as a fence test

CLDN4-only. Official CosMx NSCLC (He et al., *Nat Biotechnol* 2022; NanoString / Bruker S3 flat files) and GEO GSE307534 Visium CytAssist LUAD. Invasive LUAD samples only (GEO titles labeled LUAD; AAH / AIS / MIA / Normal excluded). No private 8-KL signature.

This is a **spatial geometry test**, not a bulk CLDN4–CD8 correlation. The unit is a tumor domain. The x-axis is micrometers from the tumor–stroma margin. The null is a rotated/shifted CLDN4 field with the tumor mask and CD8 held fixed. A fence would be: CD8 high in immediate stroma, a drop across the first 50 µm of a CLDN4-high rim, and a lower 0–200 µm CD8 AUC inside CLDN4-high domains than inside CLDN4-low domains.

## Data

| Cohort | Source | Units | Filter | Tumor definition | Fence gene | CD8 |
|---|---|---|---|---|---|---|
| CosMx NSCLC | Official S3 flat files, 8 FFPE sections / 5 patients (Lung5 ×3, Lung6, Lung9 ×2, Lung12, Lung13) | 800,327 cells | drop `cell_ID==0` | EPCAM/KRT RNA GMM + PanCK protein | `CLDN4` only | CD8+ = CD8A+CD8B ≥ 1; density among **all** cells in each distance band |
| Visium LUAD | GEO GSE307534, 26 LUAD-titled slides | 296,126 spots | invasive LUAD label only | EPCAM/KRT RNA GMM | `CLDN4` only | CD8A/B log1p-normalized module |

CosMx scale is the official 0.18 µm/pixel (NanoString SMI-ReadMe). Visium coordinates use `spot_diameter_fullres` → 55 µm.

## Domain and core definitions

- **Tumor domain**: connected component of tumor cells/spots (link radius 25 µm CosMx / 150 µm Visium). Domains smaller than 80 cells (CosMx) or 12 spots (Visium) are discarded.
- **CLDN4-high core**: inside tumor, Getis-Ord Gi\* z ≥ 1.645 in a 50 µm (CosMx) / 150 µm (Visium) neighborhood **and** CLDN4 ≥ tumor median, unioned with local CLDN4 density peaks (neighborhood mean in the top 20%).
- **CLDN4-high domain**: ≥15% of its tumor units sit in a CLDN4-high core.
- Claudin-4 is the only fence molecule. No multi-gene tight-junction score.

## Metrics

- **Infiltration curve**: CD8 vs signed Euclidean distance-from-margin (distance transform of the rasterized tumor mask). 20 µm bins on CosMx, 50 µm bins on Visium. Every cell/spot in the band is counted, including infiltrating immune cells, not only EPCAM/KRT-high units.
- **AUC₀–₂₀₀**: trapezoidal area under the inward curve. Lower AUC = less CD8 inside.
- **ΔAUC**: AUC(CLDN4-low) − AUC(CLDN4-high). Positive would mean CLDN4-high domains exclude CD8 more.
- **Barrier index**: (CD8_stroma − CD8_rim) / (CD8_stroma + CD8_rim). CosMx rim = first 50 µm inside vs immediate 50 µm of stroma, restricted to the CLDN4-high rim when cores contact the margin. Visium rim = 100 µm (one spot ring); 50 µm is below Visium Nyquist. Range (−1, 1). **Positive = drop across the fence.**
- **Permutation**: rotate and shift the CLDN4 field (per FOV on CosMx; whole slide on Visium), recompute cores and the barrier, 199 (CosMx) / 99 (Visium) times. Empirical p = (1 + #{perm ≥ observed}) / (1 + n). Tests whether an exclusion-signed barrier is registered to CLDN4 geography.

## Results

### CosMx NSCLC — the geometrically decisive assay

Eight sections, **242** tumor domains (110 CLDN4-high, 132 CLDN4-low).

| Metric | CLDN4-high | CLDN4-low | Test |
|---|---|---|---|
| Barrier index | **−0.196 ± 0.031** (n=8) | **−0.163 ± 0.025** (n=8) | paired Wilcoxon (high > low): p = 0.95 |
| CD8+ fraction, immediate stroma | 0.100 ± 0.012 | — | — |
| CD8+ fraction, first 50 µm inside | 0.149 ± 0.018 | — | raw drop −0.048 ± 0.010 |
| AUC₀–₂₀₀ (cells·µm / mm²) | 191,330 ± 24,378 | 151,924 ± 18,619 | ΔAUC = **−39,405 ± 16,537**; Wilcoxon vs 0: p = 0.99 |
| Rotate/shift p (barrier, Stouffer) | **0.50** | — | not in the exclusion tail |

Every CosMx section has a **negative** barrier index: CD8 density *rises* from immediate stroma into the first 50 µm of tumor. The rise is at least as large at CLDN4-high rims as at CLDN4-low rims. 0–200 µm CD8 AUC is **higher**, not lower, inside CLDN4-high domains (7/8 sections have ΔAUC < 0).

| Sample | Cells | Domains (high/low) | Barrier high | Barrier low | ΔAUC | perm p |
|---|---:|---:|---:|---:|---:|---:|
| Lung12 | 73,997 | 15/11 | −0.236 | −0.179 | −50,268 | 0.71 |
| Lung13 | 82,843 | 18/31 | −0.146 | −0.156 | +32,574 | 0.14 |
| Lung5_Rep1 | 100,292 | 14/14 | −0.199 | −0.119 | −3,775 | 0.96 |
| Lung5_Rep2 | 106,660 | 14/12 | −0.085 | −0.076 | −78,756 | 0.84 |
| Lung5_Rep3 | 100,264 | 19/8 | −0.087 | −0.092 | −6,001 | 0.38 |
| Lung6 | 93,795 | 9/10 | −0.217 | −0.160 | −51,422 | 0.08 |
| Lung9_Rep1 | 91,972 | 6/17 | −0.339 | −0.260 | −116,444 | 0.91 |
| Lung9_Rep2 | 150,504 | 15/29 | −0.258 | −0.263 | −41,151 | 0.04 |

Lung9_Rep2 is the only section whose rotate/shift p falls below 0.05 (p = 0.04), and its absolute barrier is still negative. Combined CosMx permutation p = 0.50.

The paper infiltration-depth figure (`fig1_infiltration_depth_curves.png`) shows that CosMx CD8 density increases after the margin (x = 0) in both CLDN4 classes, with the CLDN4-high curve at or above the CLDN4-low curve through the first 100 µm.

### Visium GSE307534 invasive LUAD — slide-scale replicate

26 LUAD-titled slides, **419** domains (137 CLDN4-high, 282 CLDN4-low). Barrier uses a 100 µm rim (one spot).

- Barrier index, CLDN4-high: **0.060 ± 0.046** (n=26); CLDN4-low: **0.031 ± 0.041** (n=25); paired Wilcoxon p = 0.53.
- ΔAUC (low − high): **−0.295 ± 0.270** (n=23); Wilcoxon vs 0, p = 0.97.
- Rotate/shift Stouffer p for barrier: **0.96**.

Visium is mixed at the 100 µm ring (some slides exclusion-signed, some not) and does not show a lower 0–200 µm CD8 AUC in CLDN4-high domains. It cannot resolve a 50 µm CosMx rim.

## Interpretation

The fence test was run. On official CosMx NSCLC it does **not** return a CD8-excluding CLDN4 fence.

CD8 accumulates at the invasive front: first-50 µm tumor > immediate stroma in all eight sections. CLDN4-high cores are where that front is at least as CD8-rich as CLDN4-low cores, not a dead zone behind a claudin rim. Rotate/shift of the CLDN4 field does not produce an exclusion-signed excess. GSE307534 invasive LUAD agrees on the AUC sign (CLDN4-high is not CD8-poor over 0–200 µm) and does not support a registered barrier after permutation.

That is still tumor geography. The figure is an infiltration-depth curve with a signed margin, a barrier index, and a CLDN4-field null — the correct objects for a fence claim. The measurement rejects the fence on these two official lung spatial cohorts.

This is not a patient-level correlation and not a statement that CLDN4 is irrelevant to immunity. It is a statement about **where CD8 sits relative to CLDN4-high tumor cores**, in micrometers.

## What this does not claim

- Causality (CLDN4 as a physical barrier vs a marker of an invasive-front state).
- A multi-gene keratin / claudin program. Only `CLDN4` calls cores.
- Private 8-KL or any unpublished gene set.
- Precursor lesions in GSE307534 (AAH / AIS / MIA / Normal were excluded).
- That a 50 µm CosMx rim can be read out of 100 µm Visium spots.

## Outputs

- `results/figures/fig1_infiltration_depth_curves.png` — paper infiltration-depth figure (CosMx + Visium)
- `results/figures/fig2_barrier_index_permutation.png` — barrier index and rotate/shift null
- `results/figures/fig3_spatial_fence_maps.png` — example CosMx fields (Lung12, Lung13)
- `results/tables/sample_metrics.csv`, `domain_metrics.csv`, `infiltration_curves.csv`

## Methods notes

Signed distance is a Euclidean distance transform of a rasterized tumor mask (8 µm CosMx / 25 µm Visium), not a kNN-to-stroma proxy, so holes and concave margins are handled. Gi\* uses binary weights and the Ord–Getis z-score including the self-neighbor. Permutation is a rigid rotation plus a random shift of the CLDN4 scalar field; CD8 and the tumor mask stay put. Code: `analysis/cldn4_tumor_geography.py`.
