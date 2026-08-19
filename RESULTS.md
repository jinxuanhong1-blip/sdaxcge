# RESULTS: public Visium CLDN4 vs CD8A (additive)

**CLDN4-only. Public data only. No private 8-KL. No claim language.**

- Sections attempted: **32**
- Sections with numbers (QC pass + CLDN4/CD8A present): **32**
- By source: 10x n=2, GSE189487 n=6, GSE273378 n=16, GSE300676 n=8

## Methods (pre-specified)

- Counts: Space Ranger filtered matrix (in-tissue barcodes).
- QC: ≥200 genes and ≥100 UMI per spot.
- Normalization: log1p(CPM to 10,000).
- Epithelial-like: mean of available `KRT8/KRT18/KRT19/EPCAM` ≥ section median.
- CD8A-high: CD8A ≥ 75th percentile of QC spots when that cutoff is >0; otherwise CD8A > 0. Distance excludes self.
- Distance: Euclidean µm from pixel coordinates (55 µm / `spot_diameter_fullres`, or nearest-neighbor calibrated to 100 µm).
- Neighbor CD8A: mean CD8A of Visium hex ring-1 neighbors.
- KRT8 residual: OLS residual of log1p-CPM CLDN4 on KRT8.
- Q4 vs Q1: Mann–Whitney U, two-sided, among epithelial-like spots.
- Unit: **section** (one capture area).

## Inventory

| section | dataset | histology | n_qc | n_epi | n_Q1 | n_Q4 | n_CD8A-high | status |
|---|---|---|---:|---:|---:|---:|---:|---|
| 10x_LUSC_FFPE | 10x_CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma | LUSC | 3830 | 1915 | 478 | 478 | 712 | ok |
| 10x_NEC_11mm_FFPE | 10x_CytAssist_11mm_FFPE_Human_Lung_Cancer | lung_neuroendocrine | 6195 | 3098 | 774 | 774 | 974 | ok |
| GSE189487_TD1_IAC | GSE189487 | LUAD_IAC | 4092 | 2046 | 511 | 511 | 541 | ok |
| GSE189487_TD2_IAC | GSE189487 | LUAD_IAC | 4409 | 2205 | 551 | 551 | 327 | ok |
| GSE189487_TD3_MIA | GSE189487 | LUAD_MIA | 1137 | 569 | 142 | 142 | 109 | ok |
| GSE189487_TD5_AIS | GSE189487 | LUAD_AIS | 1699 | 850 | 212 | 212 | 314 | ok |
| GSE189487_TD6_MIA | GSE189487 | LUAD_MIA | 3909 | 1955 | 488 | 488 | 196 | ok |
| GSE189487_TD8_AIS | GSE189487 | LUAD_AIS | 1760 | 880 | 220 | 220 | 257 | ok |
| GSE273378_GSM8427428_LM_SD_1216_1 | GSE273378 | LUAD | 2400 | 1200 | 300 | 300 | 322 | ok |
| GSE273378_GSM8427429_LM_SD_16 | GSE273378 | LUAD | 2755 | 1378 | 344 | 344 | 54 | ok |
| GSE273378_GSM8427430_LM_SD_11 | GSE273378 | LUAD | 2724 | 1362 | 340 | 340 | 143 | ok |
| GSE273378_GSM8427431_LM_SD_2 | GSE273378 | LUAD | 1728 | 864 | 216 | 216 | 232 | ok |
| GSE273378_GSM8427432_LM_SD_3 | GSE273378 | LUAD | 3853 | 1927 | 481 | 481 | 853 | ok |
| GSE273378_GSM8427433_LM_SD_4 | GSE273378 | LUAD | 3500 | 1750 | 437 | 437 | 875 | ok |
| GSE273378_GSM8427434_LM_SD_5 | GSE273378 | LUAD | 3230 | 1615 | 403 | 403 | 154 | ok |
| GSE273378_GSM8427435_LM_SD_6 | GSE273378 | LUAD | 2530 | 1265 | 316 | 316 | 322 | ok |
| GSE273378_GSM8427436_LM_SD_7 | GSE273378 | LUAD | 2753 | 1377 | 344 | 344 | 689 | ok |
| GSE273378_GSM8427437_LM_SD_1216_8 | GSE273378 | LUAD | 2756 | 1378 | 344 | 344 | 377 | ok |
| GSE273378_GSM8427438_LM_SD_9 | GSE273378 | LUAD | 2801 | 1401 | 350 | 350 | 290 | ok |
| GSE273378_GSM8427439_LM_SD_10 | GSE273378 | LUAD | 2620 | 1310 | 327 | 327 | 379 | ok |
| GSE273378_GSM8427440_LM_SD_1216_12 | GSE273378 | LUAD | 3604 | 1802 | 450 | 450 | 94 | ok |
| GSE273378_GSM8427441_LM_SD_13 | GSE273378 | LUAD | 2912 | 1456 | 364 | 364 | 577 | ok |
| GSE273378_GSM8427442_LM_SD_1216_14 | GSE273378 | LUAD | 1915 | 958 | 239 | 239 | 357 | ok |
| GSE273378_GSM8427443_LM_SD_15 | GSE273378 | LUAD | 1741 | 871 | 217 | 217 | 134 | ok |
| GSE300676_GSM9066288_CRC_mPAP3_A | GSE300676 | LUAD | 4639 | 2320 | 580 | 580 | 125 | ok |
| GSE300676_GSM9066289_CRC_mPAP3_B | GSE300676 | LUAD | 3979 | 1990 | 497 | 497 | 655 | ok |
| GSE300676_GSM9066290_CRC_mPAP4_A | GSE300676 | LUAD | 3878 | 1939 | 484 | 484 | 970 | ok |
| GSE300676_GSM9066291_CRC_mPAP4_B | GSE300676 | LUAD | 3507 | 1754 | 438 | 438 | 877 | ok |
| GSE300676_GSM9066292_CRC_mPAP1_A | GSE300676 | LUAD | 4508 | 2254 | 563 | 563 | 1127 | ok |
| GSE300676_GSM9066293_CRC_mPAP1_B | GSE300676 | LUAD | 4463 | 2232 | 558 | 558 | 1116 | ok |
| GSE300676_GSM9066294_CRC_mPAP2_A | GSE300676 | LUAD | 4725 | 2363 | 590 | 590 | 142 | ok |
| GSE300676_GSM9066295_CRC_mPAP2_B | GSE300676 | LUAD | 3897 | 1949 | 487 | 487 | 485 | ok |

## Per-section numbers

| section | n_qc | ρ CLDN4–CD8A | p | ρ (epi) | p_epi | med dist Q4 / Q1 µm | Δdist | p_dist | med nbCD8A Q4 / Q1 | Δnb | p_nb | ρ KRT8-resid–CD8A | p_resid |
|---|---:|---:|---:|---:|---:|---|---:|---:|---|---:|---:|---:|---:|
| 10x_LUSC_FFPE | 3830 | +0.064 | 7.00e-05 | +0.047 | 0.038 | 78.8 / 78.8 | -0.0 | 0.008 | 0.068 / 0.047 | 0.021 | 0.067 | +0.065 | 6.27e-05 |
| 10x_NEC_11mm_FFPE | 6195 | +0.059 | 3.05e-06 | +0.033 | 0.065 | 91.9 / 91.9 | 0.0 | 0.876 | 0.021 / 0.020 | 0.000 | 0.573 | +0.053 | 2.97e-05 |
| GSE189487_TD1_IAC | 4092 | +0.001 | 0.931 | -0.003 | 0.903 | 100.3 / 100.3 | 0.0 | 0.163 | 0.073 / 0.092 | -0.019 | 0.027 | -0.005 | 0.755 |
| GSE189487_TD2_IAC | 4409 | -0.011 | 0.472 | +0.022 | 0.306 | 173.5 / 174.0 | -0.5 | 0.444 | 0.000 / 0.000 | 0.000 | 0.570 | -0.011 | 0.465 |
| GSE189487_TD3_MIA | 1137 | -0.036 | 0.220 | -0.034 | 0.422 | 174.0 / 173.5 | 0.5 | 0.321 | 0.000 / 0.000 | 0.000 | 0.509 | -0.038 | 0.201 |
| GSE189487_TD5_AIS | 1699 | -0.097 | 6.13e-05 | -0.036 | 0.291 | 100.5 / 100.4 | 0.1 | 0.676 | 0.090 / 0.104 | -0.014 | 0.260 | -0.057 | 0.020 |
| GSE189487_TD6_MIA | 3909 | -0.025 | 0.114 | -0.009 | 0.687 | 200.7 / 200.7 | 0.0 | 0.374 | 0.000 / 0.000 | 0.000 | 0.846 | -0.021 | 0.192 |
| GSE189487_TD8_AIS | 1760 | +0.077 | 0.001 | +0.078 | 0.020 | 100.2 / 100.5 | -0.4 | 0.011 | 0.059 / 0.054 | 0.005 | 0.195 | +0.074 | 0.002 |
| GSE273378_GSM8427428_LM_SD_1216_1 | 2400 | +0.047 | 0.020 | +0.006 | 0.842 | 85.1 / 146.6 | -61.5 | 6.88e-06 | 0.059 / 0.000 | 0.059 | 0.024 | +0.040 | 0.053 |
| GSE273378_GSM8427429_LM_SD_16 | 2755 | -0.008 | 0.681 | -0.014 | 0.593 | 305.2 / 293.4 | 11.8 | 0.022 | 0.000 / 0.000 | 0.000 | 0.474 | -0.006 | 0.766 |
| GSE273378_GSM8427430_LM_SD_11 | 2724 | -0.044 | 0.022 | -0.018 | 0.499 | 224.3 / 224.1 | 0.3 | 0.661 | 0.000 / 0.000 | 0.000 | 0.812 | +0.000 | 0.991 |
| GSE273378_GSM8427431_LM_SD_2 | 1728 | -0.047 | 0.049 | -0.062 | 0.070 | 146.5 / 146.7 | -0.2 | 0.248 | 0.000 / 0.000 | 0.000 | 0.473 | -0.019 | 0.422 |
| GSE273378_GSM8427432_LM_SD_3 | 3853 | +0.007 | 0.646 | +0.014 | 0.549 | 84.8 / 84.8 | -0.0 | 0.916 | 0.141 / 0.151 | -0.010 | 0.287 | +0.011 | 0.499 |
| GSE273378_GSM8427433_LM_SD_4 | 3500 | +0.029 | 0.087 | +0.001 | 0.977 | 84.8 / 84.8 | 0.0 | 0.012 | 0.300 / 0.315 | -0.015 | 0.239 | +0.001 | 0.932 |
| GSE273378_GSM8427434_LM_SD_5 | 3230 | -0.005 | 0.785 | +0.010 | 0.698 | 223.8 / 224.2 | -0.4 | 0.103 | 0.000 / 0.000 | 0.000 | 0.393 | -0.011 | 0.536 |
| GSE273378_GSM8427435_LM_SD_6 | 2530 | +0.054 | 0.007 | +0.077 | 0.006 | 85.0 / 115.9 | -30.8 | 0.135 | 0.258 / 0.091 | 0.167 | 0.757 | +0.002 | 0.923 |
| GSE273378_GSM8427436_LM_SD_7 | 2753 | +0.072 | 1.70e-04 | +0.040 | 0.142 | 84.8 / 84.8 | -0.0 | 3.13e-04 | 0.147 / 0.142 | 0.005 | 0.326 | +0.069 | 2.75e-04 |
| GSE273378_GSM8427437_LM_SD_1216_8 | 2756 | +0.020 | 0.301 | +0.018 | 0.512 | 146.5 / 85.1 | 61.5 | 0.053 | 0.000 / 0.109 | -0.109 | 0.023 | +0.009 | 0.626 |
| GSE273378_GSM8427438_LM_SD_9 | 2801 | +0.009 | 0.632 | +0.006 | 0.835 | 146.7 / 147.1 | -0.4 | 2.06e-04 | 0.000 / 0.000 | 0.000 | 0.781 | +0.007 | 0.709 |
| GSE273378_GSM8427439_LM_SD_10 | 2620 | +0.026 | 0.182 | +0.015 | 0.579 | 85.0 / 146.5 | -61.4 | 0.018 | 0.055 / 0.000 | 0.055 | 0.833 | +0.041 | 0.034 |
| GSE273378_GSM8427440_LM_SD_1216_12 | 3604 | -0.018 | 0.291 | +0.001 | 0.970 | 224.4 / 225.0 | -0.6 | 0.081 | 0.000 / 0.000 | 0.000 | 0.156 | -0.014 | 0.415 |
| GSE273378_GSM8427441_LM_SD_13 | 2912 | +0.026 | 0.159 | +0.000 | 0.996 | 84.8 / 84.8 | -0.0 | 0.016 | 0.110 / 0.139 | -0.028 | 0.248 | +0.009 | 0.616 |
| GSE273378_GSM8427442_LM_SD_1216_14 | 1915 | -0.047 | 0.038 | -0.032 | 0.324 | 84.8 / 84.8 | -0.0 | 0.291 | 0.117 / 0.106 | 0.011 | 0.393 | -0.039 | 0.091 |
| GSE273378_GSM8427443_LM_SD_15 | 1741 | -0.004 | 0.859 | -0.033 | 0.335 | 169.2 / 146.9 | 22.3 | 0.062 | 0.000 / 0.000 | 0.000 | 0.479 | -0.002 | 0.919 |
| GSE300676_GSM9066288_CRC_mPAP3_A | 4639 | -0.034 | 0.019 | +0.005 | 0.812 | 293.6 / 254.7 | 38.9 | 0.736 | 0.000 / 0.000 | 0.000 | 0.868 | -0.032 | 0.029 |
| GSE300676_GSM9066289_CRC_mPAP3_B | 3979 | -0.014 | 0.361 | +0.021 | 0.342 | 85.2 / 85.3 | -0.1 | 0.354 | 0.037 / 0.034 | 0.002 | 0.033 | +0.010 | 0.533 |
| GSE300676_GSM9066290_CRC_mPAP4_A | 3878 | +0.010 | 0.526 | +0.025 | 0.276 | 84.9 / 84.9 | -0.0 | 0.291 | 0.149 / 0.153 | -0.003 | 0.757 | +0.012 | 0.441 |
| GSE300676_GSM9066291_CRC_mPAP4_B | 3507 | -0.014 | 0.404 | -0.011 | 0.632 | 84.8 / 84.8 | -0.0 | 0.634 | 0.131 / 0.163 | -0.031 | 0.005 | -0.016 | 0.348 |
| GSE300676_GSM9066292_CRC_mPAP1_A | 4508 | +0.017 | 0.248 | +0.070 | 8.37e-04 | 84.6 / 84.6 | 0.0 | 0.450 | 0.093 / 0.087 | 0.006 | 0.218 | +0.025 | 0.098 |
| GSE300676_GSM9066293_CRC_mPAP1_B | 4463 | +0.005 | 0.759 | -0.023 | 0.280 | 84.9 / 84.8 | 0.1 | 0.007 | 0.182 / 0.215 | -0.034 | 0.009 | -0.024 | 0.109 |
| GSE300676_GSM9066294_CRC_mPAP2_A | 4725 | -0.043 | 0.003 | -0.043 | 0.035 | 305.8 / 294.6 | 11.2 | 0.392 | 0.000 / 0.000 | 0.000 | 0.356 | -0.034 | 0.018 |
| GSE300676_GSM9066295_CRC_mPAP2_B | 3897 | -0.052 | 0.001 | -0.044 | 0.051 | 146.7 / 146.7 | 0.1 | 0.091 | 0.000 / 0.000 | 0.000 | 0.014 | -0.047 | 0.003 |

## Cross-section summary (honest n = number of QC-pass sections)

- n sections = **32**
- Median Spearman CLDN4 vs CD8A (all QC spots): -0.001 (IQR -0.028 to +0.026; n=32)
- Signs of ρ CLDN4–CD8A: 16 positive / 16 negative / 0 zero; two-sided sign-test p=1.000
- Median Spearman CLDN4 vs CD8A (epithelial-like): +0.003 (IQR -0.019 to +0.021; n=32)
- Median Δ distance Q4−Q1 (µm; positive = Q4 farther from CD8A-high): -0.000 (IQR -0.264 to +0.083; n=32)
- Signs of Δ distance: 13 positive / 17 negative / 2 zero; two-sided sign-test p=0.585
- Median Δ neighbor CD8A Q4−Q1 (negative = Q4 has lower neighbor CD8A): +0.000 (IQR -0.005 to +0.003; n=32)
- Signs of Δ neighbor CD8A: 10 positive / 9 negative / 13 zero; two-sided sign-test p=1.000
- Median Spearman KRT8-residual CLDN4 vs CD8A: -0.001 (IQR -0.020 to +0.011; n=32)
- Signs of KRT8-residual ρ: 16 positive / 16 negative / 0 zero; two-sided sign-test p=1.000

Maps: `methods/visium_10x_nsclc_cldn4/results/figures/map_<section>.png`. Summary: `summary_per_section.png`.

## Skipped / not used

- Private 8-KL: not used.
- JGAS / HUM0394 / EGA controlled: skipped.
- Visium HD 10x LUAD (fixed frozen) binned_outputs.tar.gz is 8.9 GB; not downloaded here (separate HD hunt).
- GSE307534 (56 Visium LUAD/precursor, 9.4 GB) and GSE277206: left to dedicated hunts.
- GSE303162: mouse LLC Visium, not human.
- GSE300676 sample names are `CRC_mPAP*` as deposited; the GEO series is LUAD micropapillary Visium.

