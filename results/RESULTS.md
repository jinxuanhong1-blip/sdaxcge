# Barrier index and CD8 infiltration-depth AUC at CLDN4-high vs CLDN4-low tumor domains

CLDN4-only. Primary endpoints are the **barrier index** and the **0–200 µm CD8 infiltration-depth AUC**, compared between CLDN4-high and CLDN4-low tumor domains. This is not a nearest-µm cell–cell distance and not a bulk correlation.

Official CosMx NSCLC: **all 8 sections** (Lung5_Rep1/2/3, Lung6, Lung9_Rep1/2, Lung12, Lung13; He et al. 2022; NanoString/Bruker S3). GEO GSE307534 Visium LUAD: invasive (LUAD-titled) slides only. Tumor = EPCAM/KRT RNA GMM plus PanCK protein on CosMx. **Tumor is not gated as CD8A==0.** Infiltrating CD8 inside the epithelial mask stays in the CD8 numerator.

## Primary endpoints — CosMx, all 8 official sections

| Endpoint | CLDN4-high | CLDN4-low | Comparison |
|---|---|---|---|
| **Barrier index** (section means) | **−0.196 ± 0.031** (n=8) | **−0.163 ± 0.025** (n=8) | paired Wilcoxon high>low p = 0.95 |
| **CD8 AUC₀–₂₀₀** (section means) | **191,330 ± 24,378** | **151,924 ± 18,619** | paired Wilcoxon high>low p = **0.020** |
| ΔAUC (low − high) | **−39,405 ± 16,537** | — | Wilcoxon vs 0 p = **0.039** |
| Barrier, domain-level | −0.171 ± 0.017 (n=110) | −0.222 ± 0.016 (n=132) | MWU p = 0.050 |
| AUC₀–₂₀₀, domain-level | 96,526 ± 5,457 (n=110) | 90,080 ± 3,426 (n=132) | MWU p = 0.89 |
| Rotate/shift p (barrier, Stouffer) | 0.50 | — | CLDN4 field vs tumor outline |

Barrier index = (CD8_stroma − CD8_rim) / (CD8_stroma + CD8_rim) on the first 50 µm inside the tumor versus immediate stroma. Positive would be a drop across a CLDN4 rim. AUC₀–₂₀₀ is the trapezoid of CD8 density versus inward distance-from-margin (signed distance transform), 0–200 µm.

On all eight CosMx sections the section-level barrier is negative (rim CD8 > stroma CD8). Section-level infiltration-depth AUC is **higher** in CLDN4-high domains than in CLDN4-low domains (7/8 sections; paired p = 0.020). Domain-level AUC does not separate the classes (p = 0.89); domain-level barrier is slightly less negative in CLDN4-high domains (p = 0.050).

| Sample | Cells | Domains high/low | Barrier high | Barrier low | AUC high | AUC low | ΔAUC | perm p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Lung5_Rep1 | 100,292 | 14/14 | −0.199 | −0.119 | 127,407 | 123,631 | −3,775 | 0.96 |
| Lung5_Rep2 | 106,660 | 14/12 | −0.085 | −0.076 | 195,108 | 116,352 | −78,756 | 0.84 |
| Lung5_Rep3 | 100,264 | 19/8 | −0.087 | −0.092 | 86,667 | 80,666 | −6,001 | 0.38 |
| Lung6 | 93,795 | 9/10 | −0.217 | −0.160 | 154,021 | 102,599 | −51,422 | 0.08 |
| Lung9_Rep1 | 91,972 | 6/17 | −0.339 | −0.260 | 289,443 | 172,998 | −116,444 | 0.91 |
| Lung9_Rep2 | 150,504 | 15/29 | −0.258 | −0.263 | 229,368 | 188,217 | −41,151 | 0.04 |
| Lung12 | 73,997 | 15/11 | −0.236 | −0.179 | 266,641 | 216,373 | −50,268 | 0.71 |
| Lung13 | 82,843 | 18/31 | −0.146 | −0.156 | 181,983 | 214,557 | +32,574 | 0.14 |

800,327 cells; 242 tumor domains (110 CLDN4-high, 132 CLDN4-low).

## Visium GSE307534 invasive LUAD (supporting)

26 LUAD-titled slides. Barrier uses a 100 µm rim (one spot ring). Same two endpoints.

- Barrier index: CLDN4-high **0.060 ± 0.046** (n=26) vs CLDN4-low **0.031 ± 0.041** (n=25); paired Wilcoxon p = 0.53.
- ΔAUC (low − high): **−0.295 ± 0.270** (n=23); Wilcoxon vs 0 p = 0.97.
- Domain-level barrier MWU p = 0.66; domain-level AUC MWU p = 0.27.
- Rotate/shift Stouffer p = 0.96.

## Definitions (CLDN4-only)

- **Tumor domain**: connected EPCAM/KRT-high (± PanCK) component. CosMx link 25 µm, min 80 cells. **No CD8A==0 gate.**
- **CLDN4-high core**: Gi\* z ≥ 1.645 and CLDN4 ≥ tumor median, or local CLDN4 density peak (top 20%), inside tumor. Fence gene is `CLDN4` only.
- **CLDN4-high domain**: ≥15% of tumor units in a CLDN4-high core.
- **CD8**: CosMx CD8A+CD8B ≥ 1 among all cells in the distance band (density per mm²). Visium CD8A/B log1p module.
- **Permutation**: rotate/shift the CLDN4 field; tumor mask and CD8 stay put. Empirical p = (1 + #{perm ≥ observed}) / (1 + n).

No nearest-µm pairwise distance is used as an endpoint. No private 8-KL set. GSE307534 precursors (AAH/AIS/MIA/Normal) are excluded.

## Figure

`results/figures/fig1_barrier_and_auc.png` is the paper panel: CosMx barrier (A), CosMx AUC₀–₂₀₀ (B), Visium barrier (C), CosMx 0–200 µm depth curve (D). Supporting depth curves: `fig2_infiltration_depth_curves.png`. Spatial maps: `fig3_spatial_fence_maps.png`.

## Methods notes

Signed distance is a Euclidean distance transform of a rasterized tumor mask (8 µm CosMx / 25 µm Visium). Gi\* uses binary weights and the Ord–Getis z-score including self. Code: `analysis/cldn4_tumor_geography.py`. Tables: `results/tables/sample_metrics.csv`, `domain_metrics.csv`, `infiltration_curves.csv`.
