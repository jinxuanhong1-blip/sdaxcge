# RESULTS: Visium spatial lag (CLDN4 vs neighbor CD8A)

**CLDN4-only. Public Visium LUAD/NSCLC. No private KL. No claim-failed. No fabrication.**

Primary readout is **spatial lag**, not same-spot Spearman and not nearest-spot µm.
A 55 µm Visium spot mixes tumor and T cells; within-spot ρ and nearest-spot distance
are the wrong exclusion test and are not reported as the punch.

- Sections attempted: **87**
- Sections kept (QC + CLDN4/CD8A contrast): **71**
- Sections dropped (no contrast / load / QC): **16**
- Kept by source: GSE307534 n=46, GSE273378 n=12, GSE189487 n=6, GSE300676 n=6, 10x n=1

## Methods (pre-specified)

- Counts: Space Ranger filtered matrix; in-tissue barcodes.
- QC: ≥200 genes and ≥100 UMI.
- Normalization: log1p(CPM to 10,000).
- Epithelial-like: mean of available `KRT8/KRT18/KRT19/EPCAM` ≥ section median.
- Neighbors: Visium hex ring-1 (self excluded). Sensitivity: Euclidean 80–150 µm annulus.
- Lag-ρ: Spearman(index CLDN4, mean neighbor CD8A) on epithelial-like spots with ≥1 neighbor.
- Q4 vs Q1: among those spots, equal-sized bottom/top 25% of CLDN4; Δ = mean(neighbor CD8A | Q4) − mean(neighbor CD8A | Q1).
- Morisita–Horn: 200 µm grid counts of CLDN4-Q4 epi vs CD8A-high spots (q75, or CD8A>0 if q75=0).
- KRT8 residual: OLS residual of index CLDN4 on KRT8, then Spearman vs *neighbor* CD8A.
- Contrast drop: CD8A+ spots <50 or <5% of QC, **or** CLDN4+ among epi <10% or IQR=0, **or** <80 epi spots with a neighbor, **or** neighbor-CD8A SD=0.
- Optional interface: epi-like spots with ≥1 non-epi ring-1 neighbor.
- Cross-section: Wilcoxon signed-rank vs 0 on section lag-ρ and on section Δ.
- Unit: **section** (one capture area).

## Primary: section-level spatial lag

- n sections = **71**
- Lag-ρ (ring-1): median +0.040 (IQR -0.011 to +0.096)
- Signs: **25 negative / 46 positive** / 0 zero
- Wilcoxon signed-rank on section lag-ρ vs 0: p=2.74e-05

- Q4−Q1 Δ mean neighbor CD8A: median 0.000 (IQR -0.009 to 0.001)
- Signs of Δ (negative = Q4 has lower neighbor CD8A): **28 negative / 43 positive** / 0 zero
- Wilcoxon signed-rank on section Δ vs 0: p=0.409

Negative lag-ρ / negative Δ is the exclusion direction. The pooled section
vector is **not** exclusion-signed for lag-ρ (more positive than negative).
Q4−Q1 Δ is consistent with zero.

### By source (same pre-specified metrics; one series dominates n)

| source | n | lag-ρ median | lag n_neg / n_pos | lag Wilcoxon p | Δ median | Δ n_neg / n_pos | Δ Wilcoxon p |
|---|---:|---:|---|---:|---:|---|---:|
| 10x | 1 | +0.040 | 0 / 1 | NA | -0.001 | 1 / 0 | NA |
| GSE189487 | 6 | -0.020 | 5 / 1 | 0.438 | -0.011 | 4 / 2 | 0.219 |
| GSE273378 | 12 | -0.011 | 9 / 3 | 0.424 | -0.012 | 8 / 4 | 0.064 |
| GSE300676 | 6 | -0.054 | 4 / 2 | 0.156 | -0.018 | 5 / 1 | 0.062 |
| GSE307534 | 46 | +0.085 | 7 / 39 | 2.77e-08 | 0.000 | 10 / 36 | 0.019 |

## Mixing (Morisita–Horn)

- Morisita–Horn (200 µm grid; 0=segregated, 1=mixed): median 0.387 (IQR 0.332 to 0.437; n=71)
- Wilcoxon of (Morisita−0.5) vs 0 (negative = more segregated than 0.5): p=1.20e-12 (65 below 0.5 / 6 above)

## KRT8 residual on the lag (not same-spot)

- Lag-ρ of (CLDN4 | KRT8 residual) vs ring-1 CD8A: median +0.014 (IQR -0.018 to +0.042)
- Signs: **26 negative / 45 positive**; Wilcoxon p=0.043

## Sensitivity

- 80–150 µm annulus lag-ρ: median +0.045; 24 neg / 47 pos; Wilcoxon p=2.95e-05
- Interface-only lag-ρ: n=71, median +0.039; 17 neg / 54 pos; Wilcoxon p=4.88e-07
- Interface Q4−Q1 Δ neighbor CD8A: median 0.000; 20 neg / 51 pos; Wilcoxon p=0.056

Forest: `methods/visium_spatial_lag_cldn4/results/figures/forest_lag_and_delta.png`.

## Per-section lag table

| section | dataset | histo | n_qc | n_epi+nbr | lag-ρ ring-1 (95% CI) | p | mean nbCD8A Q4 / Q1 | Δ | Morisita | resid lag-ρ |
|---|---|---|---:|---:|---|---:|---|---:|---:|---:|
| 10x_LUSC_FFPE | 10x_CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma | LUSC | 3830 | 1915 | +0.040 (-0.005 to +0.084) | 0.082 | 0.090 / 0.091 | -0.001 | 0.511 | +0.023 |
| GSE189487_TD1_IAC | GSE189487 | LUAD_IAC | 4092 | 2046 | -0.036 (-0.079 to +0.008) | 0.105 | 0.092 / 0.110 | -0.017 | 0.334 | -0.034 |
| GSE189487_TD2_IAC | GSE189487 | LUAD_IAC | 4409 | 2204 | -0.010 (-0.051 to +0.032) | 0.652 | 0.062 / 0.067 | -0.005 | 0.320 | -0.013 |
| GSE189487_TD3_MIA | GSE189487 | LUAD_MIA | 1137 | 567 | -0.029 (-0.111 to +0.053) | 0.485 | 0.063 / 0.081 | -0.018 | 0.211 | -0.027 |
| GSE189487_TD5_AIS | GSE189487 | LUAD_AIS | 1699 | 848 | -0.037 (-0.104 to +0.030) | 0.280 | 0.109 / 0.128 | -0.019 | 0.300 | -0.058 |
| GSE189487_TD6_MIA | GSE189487 | LUAD_MIA | 3909 | 1954 | -0.010 (-0.054 to +0.035) | 0.663 | 0.035 / 0.032 | 0.003 | 0.203 | -0.010 |
| GSE189487_TD8_AIS | GSE189487 | LUAD_AIS | 1760 | 880 | +0.050 (-0.016 to +0.116) | 0.138 | 0.077 / 0.069 | 0.008 | 0.484 | +0.061 |
| GSE273378_LM_SD_1216_1 | GSE273378 | LUAD_stageI | 2400 | 1199 | +0.064 (+0.007 to +0.120) | 0.028 | 0.095 / 0.080 | 0.015 | 0.371 | +0.066 |
| GSE273378_LM_SD_11 | GSE273378 | LUAD_stageI | 2724 | 1361 | -0.011 (-0.064 to +0.042) | 0.676 | 0.029 / 0.030 | -0.001 | 0.194 | -0.009 |
| GSE273378_LM_SD_2 | GSE273378 | LUAD_stageI | 1728 | 864 | -0.035 (-0.101 to +0.032) | 0.304 | 0.066 / 0.085 | -0.019 | 0.363 | +0.018 |
| GSE273378_LM_SD_3 | GSE273378 | LUAD_stageI | 3853 | 1927 | -0.019 (-0.064 to +0.026) | 0.401 | 0.157 / 0.172 | -0.015 | 0.475 | -0.020 |
| GSE273378_LM_SD_4 | GSE273378 | LUAD_stageI | 3500 | 1749 | -0.023 (-0.070 to +0.023) | 0.327 | 0.314 / 0.337 | -0.023 | 0.368 | -0.047 |
| GSE273378_LM_SD_7 | GSE273378 | LUAD_stageI | 2753 | 1369 | +0.028 (-0.025 to +0.081) | 0.299 | 0.173 / 0.169 | 0.005 | 0.408 | +0.029 |
| GSE273378_LM_SD_1216_8 | GSE273378 | LUAD_stageI | 2756 | 1377 | -0.087 (-0.139 to -0.034) | 0.001 | 0.133 / 0.163 | -0.030 | 0.301 | -0.081 |
| GSE273378_LM_SD_9 | GSE273378 | LUAD_stageI | 2801 | 1399 | -0.003 (-0.055 to +0.050) | 0.917 | 0.049 / 0.074 | -0.025 | 0.333 | -0.007 |
| GSE273378_LM_SD_10 | GSE273378 | LUAD_stageI | 2620 | 1309 | -0.008 (-0.062 to +0.047) | 0.784 | 0.078 / 0.088 | -0.010 | 0.416 | +0.018 |
| GSE273378_LM_SD_13 | GSE273378 | LUAD_stageI | 2912 | 1453 | -0.010 (-0.062 to +0.041) | 0.691 | 0.140 / 0.161 | -0.021 | 0.436 | +0.027 |
| GSE273378_LM_SD_1216_14 | GSE273378 | LUAD_stageI | 1915 | 955 | +0.033 (-0.031 to +0.096) | 0.312 | 0.136 / 0.122 | 0.013 | 0.423 | +0.004 |
| GSE273378_LM_SD_15 | GSE273378 | LUAD_stageI | 1741 | 870 | -0.024 (-0.090 to +0.043) | 0.481 | 0.048 / 0.048 | 0.000 | 0.342 | -0.020 |
| GSE300676_CRC_mPAP3_B | GSE300676 | LUAD_micropapillary | 3979 | 1990 | -0.046 (-0.090 to -0.002) | 0.041 | 0.049 / 0.064 | -0.015 | 0.297 | +0.014 |
| GSE300676_CRC_mPAP4_A | GSE300676 | LUAD_micropapillary | 3878 | 1935 | +0.000 (-0.044 to +0.045) | 0.983 | 0.162 / 0.178 | -0.016 | 0.331 | +0.006 |
| GSE300676_CRC_mPAP4_B | GSE300676 | LUAD_micropapillary | 3507 | 1754 | -0.071 (-0.118 to -0.025) | 0.003 | 0.153 / 0.183 | -0.030 | 0.442 | -0.078 |
| GSE300676_CRC_mPAP1_A | GSE300676 | LUAD_micropapillary | 4508 | 2254 | +0.034 (-0.007 to +0.075) | 0.104 | 0.114 / 0.109 | 0.005 | 0.418 | -0.035 |
| GSE300676_CRC_mPAP1_B | GSE300676 | LUAD_micropapillary | 4463 | 2229 | -0.063 (-0.104 to -0.022) | 0.003 | 0.204 / 0.237 | -0.033 | 0.326 | -0.058 |
| GSE300676_CRC_mPAP2_B | GSE300676 | LUAD_micropapillary | 3897 | 1949 | -0.065 (-0.109 to -0.020) | 0.004 | 0.053 / 0.073 | -0.020 | 0.310 | -0.071 |
| GSE307534_P1_AAH | GSE307534 | AAH_precursor | 9911 | 4956 | +0.034 (+0.006 to +0.061) | 0.018 | 0.137 / 0.145 | -0.008 | 0.414 | +0.060 |
| GSE307534_P1_LUAD | GSE307534 | LUAD_invasive | 6119 | 3059 | +0.098 (+0.063 to +0.133) | 5.87e-08 | 0.153 / 0.140 | 0.013 | 0.438 | +0.106 |
| GSE307534_P3_AIS | GSE307534 | AIS_precursor | 9902 | 4951 | -0.005 (-0.033 to +0.023) | 0.710 | 0.093 / 0.101 | -0.007 | 0.398 | +0.001 |
| GSE307534_P3_LUAD | GSE307534 | LUAD_invasive | 14172 | 7086 | -0.055 (-0.078 to -0.032) | 3.58e-06 | 0.068 / 0.090 | -0.022 | 0.340 | -0.066 |
| GSE307534_P4_AAH | GSE307534 | AAH_precursor | 12208 | 6102 | +0.032 (+0.007 to +0.057) | 0.012 | 0.099 / 0.101 | -0.002 | 0.380 | +0.029 |
| GSE307534_P4_LUAD | GSE307534 | LUAD_invasive | 12975 | 6488 | -0.067 (-0.091 to -0.043) | 6.74e-08 | 0.036 / 0.059 | -0.023 | 0.297 | -0.038 |
| GSE307534_P5_AIS | GSE307534 | AIS_precursor | 10149 | 5072 | +0.068 (+0.041 to +0.096) | 1.10e-06 | 0.063 / 0.065 | -0.003 | 0.316 | +0.074 |
| GSE307534_P5_LUAD | GSE307534 | LUAD_invasive | 12083 | 6039 | +0.163 (+0.139 to +0.188) | 2.21e-37 | 0.110 / 0.067 | 0.043 | 0.268 | +0.206 |
| GSE307534_P6_AAH | GSE307534 | AAH_precursor | 11197 | 5599 | +0.045 (+0.019 to +0.071) | 8.26e-04 | 0.002 / 0.001 | 0.000 | 0.362 | +0.026 |
| GSE307534_P6_LUAD | GSE307534 | LUAD_invasive | 11710 | 5855 | -0.048 (-0.074 to -0.023) | 2.28e-04 | 0.110 / 0.138 | -0.029 | 0.359 | -0.072 |
| GSE307534_P7_LUAD | GSE307534 | LUAD_invasive | 13615 | 6807 | +0.129 (+0.105 to +0.152) | 1.35e-26 | 0.002 / 0.002 | 0.001 | 0.482 | +0.107 |
| GSE307534_P7_LUAD-1 | GSE307534 | LUAD_invasive | 12976 | 6488 | +0.155 (+0.131 to +0.178) | 4.58e-36 | 0.003 / 0.002 | 0.001 | 0.447 | +0.114 |
| GSE307534_P9_AIS | GSE307534 | AIS_precursor | 10525 | 5263 | +0.057 (+0.030 to +0.084) | 3.83e-05 | 0.001 / 0.001 | 0.000 | 0.364 | +0.028 |
| GSE307534_P9_LUAD | GSE307534 | LUAD_invasive | 12067 | 6034 | -0.041 (-0.067 to -0.016) | 0.001 | 0.001 / 0.001 | -0.000 | 0.237 | -0.028 |
| GSE307534_P10_MIA | GSE307534 | MIA_precursor | 11973 | 5987 | +0.099 (+0.073 to +0.124) | 2.06e-14 | 0.002 / 0.001 | 0.000 | 0.370 | +0.063 |
| GSE307534_P10_LUAD | GSE307534 | LUAD_invasive | 13439 | 6720 | -0.002 (-0.025 to +0.022) | 0.900 | 0.006 / 0.006 | -0.000 | 0.367 | -0.083 |
| GSE307534_P11_AAH | GSE307534 | AAH_precursor | 10873 | 5437 | +0.137 (+0.111 to +0.163) | 3.41e-24 | 0.002 / 0.001 | 0.001 | 0.432 | +0.043 |
| GSE307534_P11_LUAD | GSE307534 | LUAD_invasive | 11782 | 5891 | +0.351 (+0.328 to +0.373) | 4.68e-170 | 0.010 / 0.004 | 0.006 | 0.470 | +0.012 |
| GSE307534_P12_AIS | GSE307534 | AIS_precursor | 12167 | 6084 | +0.163 (+0.138 to +0.187) | 2.74e-37 | 0.003 / 0.002 | 0.001 | 0.504 | +0.067 |
| GSE307534_P12_LUAD | GSE307534 | LUAD_invasive | 13824 | 6911 | +0.124 (+0.101 to +0.148) | 3.11e-25 | 0.002 / 0.001 | 0.000 | 0.395 | +0.036 |
| GSE307534_P13_MIA | GSE307534 | MIA_precursor | 13517 | 6759 | +0.170 (+0.147 to +0.193) | 3.59e-45 | 0.002 / 0.001 | 0.001 | 0.448 | +0.001 |
| GSE307534_P13_LUAD | GSE307534 | LUAD_invasive | 4860 | 2430 | +0.190 (+0.151 to +0.228) | 4.41e-21 | 0.003 / 0.002 | 0.002 | 0.504 | -0.016 |
| GSE307534_P14_AIS | GSE307534 | AIS_precursor | 12329 | 6165 | +0.072 (+0.047 to +0.097) | 1.54e-08 | 0.001 / 0.001 | 0.000 | 0.281 | +0.014 |
| GSE307534_P14_LUAD | GSE307534 | LUAD_invasive | 11477 | 5739 | +0.149 (+0.124 to +0.174) | 7.91e-30 | 0.002 / 0.001 | 0.001 | 0.405 | +0.023 |
| GSE307534_P15_MIA | GSE307534 | MIA_precursor | 13173 | 6587 | +0.092 (+0.068 to +0.116) | 6.24e-14 | 0.002 / 0.001 | 0.000 | 0.424 | +0.002 |
| GSE307534_P15_LUAD | GSE307534 | LUAD_invasive | 14207 | 7104 | +0.080 (+0.057 to +0.103) | 1.59e-11 | 0.002 / 0.001 | 0.001 | 0.213 | +0.086 |
| GSE307534_P16_AIS | GSE307534 | AIS_precursor | 12701 | 6351 | +0.276 (+0.253 to +0.298) | 4.48e-111 | 0.002 / 0.001 | 0.001 | 0.527 | +0.088 |
| GSE307534_P16_LUAD | GSE307534 | LUAD_invasive | 7596 | 3798 | +0.025 (-0.007 to +0.056) | 0.129 | 0.003 / 0.003 | 0.000 | 0.335 | +0.008 |
| GSE307534_P17_AIS | GSE307534 | AIS_precursor | 13474 | 6736 | +0.151 (+0.128 to +0.174) | 1.18e-35 | 0.001 / 0.001 | 0.000 | 0.387 | +0.069 |
| GSE307534_P17_LUAD | GSE307534 | LUAD_invasive | 13510 | 6755 | +0.032 (+0.008 to +0.056) | 0.008 | 0.004 / 0.004 | 0.000 | 0.426 | +0.008 |
| GSE307534_P18_LUAD | GSE307534 | LUAD_invasive | 11970 | 5985 | +0.087 (+0.062 to +0.112) | 1.45e-11 | 0.004 / 0.003 | 0.001 | 0.431 | +0.044 |
| GSE307534_P19_AIS | GSE307534 | AIS_precursor | 13882 | 6942 | +0.095 (+0.072 to +0.118) | 2.25e-15 | 0.005 / 0.003 | 0.001 | 0.412 | +0.019 |
| GSE307534_P19_LUAD | GSE307534 | LUAD_invasive | 12561 | 6281 | +0.120 (+0.096 to +0.145) | 9.76e-22 | 0.003 / 0.002 | 0.001 | 0.432 | +0.055 |
| GSE307534_P20_AAH | GSE307534 | AAH_precursor | 12299 | 6150 | +0.090 (+0.065 to +0.115) | 1.29e-12 | 0.009 / 0.008 | 0.002 | 0.522 | +0.031 |
| GSE307534_P20_LUAD | GSE307534 | LUAD_invasive | 13602 | 6801 | +0.144 (+0.121 to +0.167) | 8.04e-33 | 0.010 / 0.007 | 0.003 | 0.352 | +0.025 |
| GSE307534_P21_AIS | GSE307534 | AIS_precursor | 13369 | 6685 | +0.042 (+0.018 to +0.066) | 6.30e-04 | 0.003 / 0.002 | 0.001 | 0.375 | -0.013 |
| GSE307534_P21_AIS-1 | GSE307534 | AIS_precursor | 12904 | 6452 | +0.067 (+0.042 to +0.091) | 8.07e-08 | 0.003 / 0.003 | 0.000 | 0.479 | +0.042 |
| GSE307534_P21_LUAD | GSE307534 | LUAD_invasive | 13312 | 6656 | +0.045 (+0.021 to +0.069) | 2.24e-04 | 0.003 / 0.002 | 0.000 | 0.276 | -0.020 |
| GSE307534_P22_AAH | GSE307534 | AAH_precursor | 7939 | 3970 | +0.101 (+0.070 to +0.132) | 1.77e-10 | 0.002 / 0.001 | 0.001 | 0.432 | +0.050 |
| GSE307534_P22_AIS | GSE307534 | AIS_precursor | 12816 | 6408 | +0.089 (+0.065 to +0.114) | 7.68e-13 | 0.002 / 0.002 | 0.000 | 0.446 | +0.040 |
| GSE307534_P22_LUAD | GSE307534 | LUAD_invasive | 12226 | 6113 | +0.158 (+0.134 to +0.183) | 1.21e-35 | 0.008 / 0.005 | 0.003 | 0.554 | +0.061 |
| GSE307534_P23_AIS-1 | GSE307534 | AIS_precursor | 10811 | 5406 | +0.035 (+0.008 to +0.062) | 0.010 | 0.002 / 0.002 | 0.000 | 0.367 | -0.009 |
| GSE307534_P24_AAH | GSE307534 | AAH_precursor | 6598 | 3299 | +0.082 (+0.048 to +0.116) | 2.48e-06 | 0.002 / 0.002 | 0.000 | 0.476 | +0.032 |
| GSE307534_P24_LUAD | GSE307534 | LUAD_invasive | 6275 | 3138 | +0.045 (+0.010 to +0.080) | 0.011 | 0.006 / 0.005 | 0.001 | 0.418 | +0.007 |
| GSE307534_P25_AAH | GSE307534 | AAH_precursor | 9252 | 4626 | +0.026 (-0.003 to +0.055) | 0.075 | 0.003 / 0.003 | 0.000 | 0.461 | -0.002 |
| GSE307534_P25_LUAD | GSE307534 | LUAD_invasive | 3791 | 1896 | -0.085 (-0.129 to -0.040) | 2.17e-04 | 0.009 / 0.012 | -0.003 | 0.294 | -0.074 |

## Dropped / not used

- `10x_NEC_11mm_FFPE` (10x_CytAssist_11mm_FFPE_Human_Lung_Cancer): no_contrast:no_CLDN4_contrast
- `GSE273378_LM_SD_16` (GSE273378): no_contrast:no_CD8A_contrast
- `GSE273378_LM_SD_5` (GSE273378): no_contrast:no_CD8A_contrast
- `GSE273378_LM_SD_6` (GSE273378): no_contrast:no_CLDN4_contrast
- `GSE273378_LM_SD_1216_12` (GSE273378): no_contrast:no_CD8A_contrast
- `GSE300676_CRC_mPAP3_A` (GSE300676): no_contrast:no_CD8A_contrast
- `GSE300676_CRC_mPAP2_A` (GSE300676): no_contrast:no_CD8A_contrast
- `GSE307534_P2_AAH` (GSE307534): no_contrast:no_CD8A_contrast
- `GSE307534_P2_LUAD` (GSE307534): no_contrast:no_CD8A_contrast
- `GSE307534_P4_AAH-1` (GSE307534): no_contrast:no_CD8A_contrast,no_CLDN4_contrast
- `GSE307534_P8_AIS` (GSE307534): no_contrast:no_CD8A_contrast
- `GSE307534_P8_LUAD` (GSE307534): no_contrast:no_CD8A_contrast
- `GSE307534_P9_AAH` (GSE307534): no_contrast:no_CD8A_contrast,no_CLDN4_contrast
- `GSE307534_P18_MIA` (GSE307534): no_contrast:no_CD8A_contrast
- `GSE307534_P23_AIS` (GSE307534): no_contrast:no_CD8A_contrast
- `GSE307534_P23_LUAD` (GSE307534): no_contrast:no_CD8A_contrast
- Private 8-KL: not used.
- Same-spot Spearman and nearest-spot µm: computed nowhere as a result; they are the discarded readout.
- GSE300676 sample IDs are `CRC_mPAP*` as deposited; the GEO series is LUAD micropapillary Visium.
- GSE307534 normal lung sections were not requested as LUAD/NSCLC and were not downloaded.

