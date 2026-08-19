# RESULTS — CLDN4-only high-end spatial autocorrelation vs CD8A

**Scope.** Open Visium LUAD/NSCLC sections that could be downloaded quickly.
CLDN4-only (no TACSTD2). No private 8-KL. No GSE307534 re-download (9.4 GB RAW / image-heavy per-sample tars).
Series without both **CLDN4** and **CD8A** are skipped. No claim is written if the test does not support it.

## Datasets run

| Series | Sections (QC OK) | Histology | Notes |
|---|---:|---|---|
| 10x_NSCLC_demo | 2 | LUSC_FFPE, NEC_FFPE | Public 10x NSCLC demos (SCC 6.5 mm; neuroendocrine 11 mm). Not LUAD. |
| GSE189487 | 6 | LUAD_early | Early LUAD Visium with tissue_positions. |
| GSE273378 | 16 | LUAD_stageI | Stage I LUAD Visium (16 sections). Matrices + coords only; images not kept. |
| GSE277206 | 2 | LUAD_early_FFPE | CytAssist FFPE early LUAD (never-smoker progression). H5 only; array coords from 10x 11 mm barcode map. |
| Zenodo_13337961 | 2 | LUAD_lepidic, LUAD_solid | Two LUAD FFPE CytAssist sections (lepidic / solid). Matrices only; array coords from 6.5 mm barcode map. |

QC-passed sections: **28**. LUAD-labelled: **26**.

## Methods (short)

- Spots: Space Ranger filtered matrices; QC min 100 counts and 50 genes.
- Expression: log1p counts-per-10k.
- Weights: k=6 KNN on Visium array (row, col); row-standardized for Moran / Lee / SLX / SEM; binary for Gi*.
- Bivariate Moran's I: `esda.Moran_BV` (CLDN4 vs CD8A), 299 permutations.
- Lee's L: Lee 2001 spatial Pearson, L = Z'(V'V)Z / 1'(V'V)1 off-diagonal, 299 permutations.
- Partial Moran / Lee: residualize CLDN4 and CD8A on KRT8 (OLS), then repeat.
- Gi*: `esda.G_Local(..., star=True)` analytic z; hot/cold |z| > 1.96.
- SLX: OLS `CD8A ~ CLDN4 + KRT8 + W·CLDN4 + W·KRT8`.
- Spatial error: `spreg.GM_Error` `CD8A ~ CLDN4 + KRT8`.
- Meta: inverse-variance fixed + DerSimonian–Laird random-effects on section I; Wilcoxon signed-rank on I.

## Headline (honest)

Bivariate Moran's I (CLDN4 vs CD8A) is negative in **20/28** sections (median I = **-0.010**).
Wilcoxon signed-rank on I (alternative: I < 0): W=118.0, p=0.0267.
Random-effects meta I (all QC sections) = **-0.019** (SE 0.011, p=0.077, I²=96.8%, n=28).
LUAD-only RE meta I = **-0.023** (p=0.021, n=26).
Heterogeneity is high (I² ≈ 96%). The strongest LUAD anti-association is Zenodo solid (**I = −0.243**); the 10x neuroendocrine demo is a **positive** outlier (**I = +0.103**) and is excluded from the LUAD-only meta.

After residualizing both genes on **KRT8**, partial Moran I is negative in **19/28** sections (median **-0.008**).
Wilcoxon on partial I (alternative: I < 0): W=134.0, p=0.06.
**No claim of CLDN4-specific CD8 exclusion after epithelium control.** Raw bivariate I can be negative because CLDN4 marks epithelial/tumor spots; partial Moran after KRT8 is small or not consistently negative. That is composition, not a proven barrier niche.

## Lee's L

Median Lee's L = **-0.010** (20/28 negative). Median partial L (KRT8 residuals) = **-0.004**.

## Gi* overlap (CLDN4-hot ∩ CD8-cold)

Analytic CD8-cold (`Gi* z < −1.96`) is **essentially empty** on these Visium sections: CD8A is zero-inflated, so local Gi* cannot go far below the already-low background (median min Gi*_CD8A z = **-0.92** when present). That is a real assay limit, not a plotting bug.
Analytic overlap (CLDN4 z>1.96 ∩ CD8 z<−1.96): median fraction **0.000**.

Maps therefore mark **CD8-low as the section 10th percentile of Gi*_CD8A** (rank-based cold) and CLDN4-hot as `z > 1.96`.
Rank-based overlap fraction median = **0.029** (median expected under independence 0.036). Fisher exact one-sided enrichment p < 0.05 in **6/28** sections. A Wilcoxon test of (observed − expected) overlap is not one-sided significant across sections — **no claim of systematic CLDN4-hot / CD8-low Gi* coincidence**.
Maps: `figures/gi_overlay_<section>.png` (CLDN4, CD8A, Gi* category overlay). Forest: `figures/forest_bivariate_I.png`. Gallery: `figures/gallery_key_sections.png`.

## SLX / spatial error

SLX same-spot CLDN4 coefficient: median **-0.005** (4/28 p<0.05).
SLX spatial-lag W·CLDN4 coefficient: median **-0.010** (8/28 p<0.05).
SEM (spatial error) CLDN4 coefficient: median **-0.004**; median λ = **+0.104**.

## What this is not

- Not a private 8-sample KL GEMM analysis.
- Not a GSE307534-only redo.
- Not TACSTD2 / TROP2.
- Not a claim that CLDN4 *causes* CD8 exclusion. Gi* overlap is a hotspot coincidence test.
- 10x demos are NSCLC but not LUAD (SCC; neuroendocrine) and are tagged `is_luad=false`.

## Per-section table

| section_id | series | n_spots | moran_I | moran_I_p | lee_L | partial_moran_I | partial_moran_I_p | frac_overlap | overlap_OR | slx_CLDN4_coef | slx_W_CLDN4_coef | sem_CLDN4_coef | sem_lambda |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE277206_V0801 | GSE277206 | 10697 | -0.0252 | 0.0033 | -0.0240 | -0.0302 | 0.0033 | 0.0634 | 1.4397 | -0.0120 | -0.0755 | -0.0137 | 0.1917 |
| GSE277206_V0901 | GSE277206 | 4774 | -0.0071 | 0.2133 | -0.0046 | -0.0060 | 0.2067 | 0.0333 | 0.4259 | 0.0124 | -0.0083 | 0.0107 | 0.0716 |
| GSE189487_TD1 | GSE189487 | 4094 | -0.0035 | 0.3467 | 0.0014 | -0.0077 | 0.1900 | 0.0601 | 0.9284 | -0.0047 | -0.0107 | -0.0034 | 0.0295 |
| GSE189487_TD2 | GSE189487 | 4990 | -0.0021 | 0.3600 | -0.0047 | -0.0026 | 0.3467 | 0.0182 | 1.1826 | -0.0053 | -0.0032 | -0.0038 | 0.1030 |
| GSE189487_TD3 | GSE189487 | 1139 | -0.0481 | 0.0033 | -0.0340 | -0.0426 | 0.0033 | 0.0702 | 2.5146 | -0.0131 | -0.0733 | -0.0209 | -0.0030 |
| GSE189487_TD5 | GSE189487 | 1700 | -0.1321 | 0.0033 | -0.1290 | -0.0506 | 0.0033 | 0.0447 | 2.0899 | -0.0110 | -0.0877 | -0.0386 | 0.0920 |
| GSE189487_TD6 | GSE189487 | 3911 | -0.0099 | 0.0833 | -0.0086 | -0.0092 | 0.1033 | 0.0251 | 1.3692 | -0.0064 | -0.0060 | -0.0066 | -0.0052 |
| GSE189487_TD8 | GSE189487 | 1760 | 0.0169 | 0.1367 | 0.0241 | 0.0140 | 0.1967 | 0.0170 | 0.4521 | 0.0097 | 0.0134 | 0.0124 | 0.0452 |
| GSE273378_LM_SD_1216_1 | GSE273378 | 2401 | 0.0338 | 0.0067 | 0.0337 | 0.0268 | 0.0167 | 0.0521 | 0.6351 | -0.0051 | 0.0320 | 0.0027 | 0.1324 |
| GSE273378_LM_SD_16 | GSE273378 | 2755 | -0.0118 | 0.1367 | -0.0108 | -0.0108 | 0.1200 | 0.0338 | 1.2355 | -0.0018 | -0.0111 | -0.0040 | 0.0308 |
| GSE273378_LM_SD_11 | GSE273378 | 2724 | -0.0559 | 0.0033 | -0.0568 | -0.0106 | 0.1133 | 0.0551 | 1.9157 | -0.0045 | -0.0200 | -0.0065 | 0.1165 |
| GSE273378_LM_SD_2 | GSE273378 | 1728 | -0.0529 | 0.0033 | -0.0571 | 0.0074 | 0.2500 | 0.0532 | 1.4611 | -0.0071 | 0.0224 | -0.0090 | 0.1265 |
| GSE273378_LM_SD_3 | GSE273378 | 3902 | -0.0095 | 0.0967 | -0.0060 | -0.0082 | 0.1833 | 0.0062 | 0.6104 | -0.0013 | -0.0181 | -0.0040 | 0.0953 |
| GSE273378_LM_SD_4 | GSE273378 | 3524 | 0.0255 | 0.0200 | 0.0291 | -0.0066 | 0.2233 | 0.0105 | 0.5821 | -0.0145 | -0.0144 | -0.0087 | 0.1610 |
| GSE273378_LM_SD_5 | GSE273378 | 3230 | -0.0352 | 0.0033 | -0.0361 | -0.0209 | 0.0133 | 0.0477 | 1.0051 | 0.0011 | -0.0213 | -0.0038 | 0.0477 |
| GSE273378_LM_SD_6 | GSE273378 | 3194 | 0.0273 | 0.0033 | 0.0331 | 0.0224 | 0.0033 | 0.0216 | 0.3864 | 0.0094 | 0.1109 | 0.0159 | 0.1048 |
| GSE273378_LM_SD_7 | GSE273378 | 2822 | 0.0048 | 0.3933 | 0.0167 | -0.0088 | 0.2100 | 0.0188 | 0.3968 | 0.0031 | -0.0135 | 0.0055 | 0.1108 |
| GSE273378_LM_SD_1216_8 | GSE273378 | 2836 | -0.0059 | 0.2533 | -0.0071 | -0.0115 | 0.1167 | 0.0409 | 1.1600 | 0.0060 | -0.0210 | 0.0046 | 0.0741 |
| GSE273378_LM_SD_9 | GSE273378 | 2801 | -0.0112 | 0.2000 | -0.0112 | -0.0110 | 0.1767 | 0.0618 | 0.5694 | -0.0131 | -0.0042 | -0.0083 | 0.1234 |
| GSE273378_LM_SD_10 | GSE273378 | 2620 | -0.0160 | 0.0867 | -0.0130 | 0.0058 | 0.3067 | 0.0363 | 0.7050 | 0.0067 | 0.0194 | 0.0060 | 0.0274 |
| GSE273378_LM_SD_1216_12 | GSE273378 | 3604 | -0.0073 | 0.1967 | -0.0106 | 0.0004 | 0.4700 | 0.0386 | 0.6325 | -0.0088 | 0.0106 | -0.0091 | 0.0014 |
| GSE273378_LM_SD_13 | GSE273378 | 2912 | 0.0151 | 0.1533 | 0.0144 | 0.0122 | 0.0733 | 0.0237 | 0.4100 | -0.0071 | 0.0173 | 0.0012 | 0.1210 |
| GSE273378_LM_SD_1216_14 | GSE273378 | 1915 | -0.0163 | 0.1267 | -0.0199 | -0.0094 | 0.1767 | 0.0198 | 1.1358 | -0.0298 | 0.0086 | -0.0293 | 0.1194 |
| GSE273378_LM_SD_15 | GSE273378 | 1744 | 0.0081 | 0.2400 | -0.0031 | 0.0101 | 0.1867 | 0.0195 | 0.9159 | -0.0032 | 0.0039 | -0.0028 | 0.0605 |
| 10x_NSCLC_SCC_FFPE | 10x_NSCLC_demo | 3849 | -0.0317 | 0.0033 | -0.0227 | -0.0071 | 0.1967 | 0.0099 | 0.6361 | 0.0180 | -0.0148 | 0.0129 | 0.1316 |
| 10x_NSCLC_NEC_11mm | 10x_NSCLC_demo | 6195 | 0.1029 | 0.0033 | 0.0963 | 0.1054 | 0.0033 | 0.0098 | 0.4939 | 0.0170 | 0.1795 | 0.0568 | 0.1441 |
| Zenodo13337961_lepidic | Zenodo_13337961 | 4992 | -0.0376 | 0.0033 | -0.0251 | -0.0169 | 0.0200 | 0.0066 | 0.1311 | -0.0268 | -0.0116 | -0.0179 | 0.3843 |
| Zenodo13337961_solid | Zenodo_13337961 | 3701 | -0.2427 | 0.0033 | -0.2298 | -0.0843 | 0.0033 | 0.0113 | 2.1419 | -0.0150 | -0.1439 | -0.0441 | 0.3482 |

## Reproduce

```bash
bash methods/spatial_autocorr_cldn4_cd8a/download.sh
python3 methods/spatial_autocorr_cldn4_cd8a/analyze.py
```

Outputs live under `methods/spatial_autocorr_cldn4_cd8a/{tables,figures}`.
