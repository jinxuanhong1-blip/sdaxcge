# GSE207422 public-UMI definition grid

ADDITIVE extra. Public GEO UMI only. A3 slide numbers taken as given (malignant TACSTD2 higher in NMPR than MPR; per-patient Spearman vs T/NK ρ ≈ −0.40 to −0.50). This slice does not re-litigate the slide. Author CopyKAT / epithelium RDS barcodes are **not public** and are not in the numeric grid.

Matrix: **92,330** cells × **24,292** genes (user PPT ~90,652; deposited matrix is 92,330). Unit = **12 post-treatment patients** (MPR n=4 including pCR P06; NMPR n=8). Extra figure: `fig_def_grid.png`. Gold boxes = ρ ≤ −0.35 or NMPR>MPR.

## Cell counts by definition

| Axis | Definition | n cells |
|---|---|---|
| lineage | Hu-mean epithelial | 10,669 |
| lineage | Hu-%pos epithelial | 10,667 |
| lineage | stromal (fibroblast+endothelial) | 1,303 |
| immune | T | 34,719 |
| immune | NK | 3,702 |
| immune | CD8 only | 15,999 |
| immune | CXCL13+ T/NK | 8,517 |
| immune | cytotoxicity-high T | 8,680 |
| malignant | hu_markers_mean | 6,401 |
| malignant | hu_markers_pctpos | 4,928 |
| malignant | epcam_krt_pos | 10,607 |
| malignant | epcam_krt_mean | 5,147 |
| malignant | copykat_stromal_p95 | 6,073 |
| malignant | infercnv_stromal_p95 | 1,520 |
| malignant | drmref_public | 2,051 |
| malignant | drmref_like | 7,005 |
| malignant | author CopyKAT | not public |

UCell module = {TACSTD2, CLDN4, EPCAM}, maxRank=1500. CopyKAT-like / inferCNV-like use stromal reference; they are window-smoothed expression scores, not the R packages.

## Highlight: ρ ≤ −0.35 vs T+NK

Requested flag: TACSTD2 vs T/NK fraction, ρ ≤ −0.35. n=12 keeps all four MPR samples; n<12 drops MPR samples with zero malignant cells.

| Malignant | Score | ρ | p | n |
|---|---|---|---|---|
| EPCAM+KRT (pos) | %pos UMI>0 | **-0.420** | 0.175 | 12 |
| Hu markers (mean) | mean log1p(CP10k) | **-0.467** | 0.205 | 9 |
| Hu markers (mean) | %pos UMI>0 | **-0.393** | 0.295 | 9 |
| DRMref-like rebuild | mean log1p(CP10k) | **-0.350** | 0.356 | 9 |

## Highlight: ρ ≤ −0.35 (all five immune defs)

33 combination(s). Strongest n=12 rows are EPCAM+KRT (pos) or DRMref public UCell vs CXCL13+ / cytotoxicity-high T.

| Malignant | Immune | Score | ρ | p | n |
|---|---|---|---|---|---|
| EPCAM+KRT (pos) | CXCL13+ | UCell module | **-0.657** | 0.020 | 12 |
| DRMref public labels | CXCL13+ | UCell module | **-0.643** | 0.024 | 12 |
| EPCAM+KRT (pos) | cyto-high T | UCell module | **-0.622** | 0.031 | 12 |
| DRMref public labels | cyto-high T | UCell module | **-0.622** | 0.031 | 12 |
| Hu markers (mean) | CXCL13+ | mean log1p(CP10k) | **-0.583** | 0.099 | 9 |
| Hu markers (mean) | CD8 only | mean log1p(CP10k) | **-0.567** | 0.112 | 9 |
| Hu markers (mean) | CD8+NK | mean log1p(CP10k) | **-0.483** | 0.187 | 9 |
| Hu markers (mean) | CD8+NK | %pos UMI>0 | **-0.477** | 0.194 | 9 |
| Hu markers (mean) | T+NK | mean log1p(CP10k) | **-0.467** | 0.205 | 9 |
| DRMref-like rebuild | CD8 only | mean log1p(CP10k) | **-0.467** | 0.205 | 9 |
| DRMref-like rebuild | CXCL13+ | mean log1p(CP10k) | **-0.467** | 0.205 | 9 |
| DRMref public labels | CD8+NK | UCell module | **-0.462** | 0.131 | 12 |
| Hu markers (mean) | CXCL13+ | UCell module | **-0.450** | 0.224 | 9 |
| CopyKAT-like stromal p95 | CXCL13+ | UCell module | **-0.448** | 0.145 | 12 |
| Hu markers (mean) | CD8 only | %pos UMI>0 | **-0.444** | 0.232 | 9 |
| DRMref public labels | CD8 only | UCell module | **-0.441** | 0.152 | 12 |
| Hu markers (mean) | CXCL13+ | %pos UMI>0 | **-0.435** | 0.242 | 9 |
| Hu markers (mean) | cyto-high T | mean log1p(CP10k) | **-0.433** | 0.244 | 9 |
| EPCAM+KRT (pos) | CD8 only | UCell module | **-0.427** | 0.167 | 12 |
| EPCAM+KRT (pos) | T+NK | %pos UMI>0 | **-0.420** | 0.175 | 12 |
| EPCAM+KRT (pos) | CD8+NK | UCell module | **-0.413** | 0.183 | 12 |
| Hu markers (%pos) | CXCL13+ | mean log1p(CP10k) | **-0.405** | 0.320 | 8 |
| Hu markers (mean) | CD8 only | UCell module | **-0.400** | 0.286 | 9 |
| DRMref-like rebuild | CXCL13+ | UCell module | **-0.400** | 0.286 | 9 |
| Hu markers (mean) | T+NK | %pos UMI>0 | **-0.393** | 0.295 | 9 |
| Hu markers (%pos) | CD8 only | mean log1p(CP10k) | **-0.381** | 0.352 | 8 |
| EPCAM+KRT (pos) | CXCL13+ | %pos UMI>0 | **-0.378** | 0.226 | 12 |
| EPCAM+KRT (pos) | CD8 only | %pos UMI>0 | **-0.378** | 0.226 | 12 |
| EPCAM+KRT (pos) | CD8+NK | %pos UMI>0 | **-0.371** | 0.236 | 12 |
| inferCNV-like stromal p95 | CXCL13+ | UCell module | **-0.364** | 0.245 | 12 |
| DRMref-like rebuild | T+NK | mean log1p(CP10k) | **-0.350** | 0.356 | 9 |
| DRMref-like rebuild | cyto-high T | mean log1p(CP10k) | **-0.350** | 0.356 | 9 |
| DRMref-like rebuild | cyto-high T | UCell module | **-0.350** | 0.356 | 9 |

## Highlight: NMPR>MPR (Δ > 0)

20 combination(s). Exact Wilcoxon. None of these is a retuning of the slide.

| Malignant | Score | mean NMPR | mean MPR | Δ | U | p | n (NMPR vs MPR) |
|---|---|---|---|---|---|---|---|
| Hu markers (mean) | mean log1p(CP10k) | 2.016 | 0.380 | 1.636 | 14.0 | 0.056 | 7 vs 2 |
| Hu markers (mean) | UCell module | 0.713 | 0.453 | 0.260 | 14.0 | 0.056 | 7 vs 2 |
| DRMref-like rebuild | UCell module | 0.669 | 0.523 | 0.145 | 14.0 | 0.056 | 7 vs 2 |
| Hu markers (%pos) | mean log1p(CP10k) | 1.988 | 0.362 | 1.627 | 12.0 | 0.071 | 6 vs 2 |
| Hu markers (%pos) | UCell module | 0.698 | 0.456 | 0.242 | 12.0 | 0.071 | 6 vs 2 |
| EPCAM+KRT (pos) | mean log1p(CP10k) | 1.612 | 1.258 | 0.354 | 26.0 | 0.109 | 8 vs 4 |
| DRMref public labels | mean log1p(CP10k) | 1.571 | 1.120 | 0.450 | 25.0 | 0.154 | 8 vs 4 |
| CopyKAT-like stromal p95 | mean log1p(CP10k) | 1.863 | 1.453 | 0.410 | 25.0 | 0.154 | 8 vs 4 |
| Hu markers (mean) | %pos UMI>0 | 0.849 | 0.361 | 0.488 | 12.0 | 0.194 | 7 vs 2 |
| DRMref public labels | UCell module | 0.651 | 0.584 | 0.067 | 24.0 | 0.214 | 8 vs 4 |
| Hu markers (%pos) | %pos UMI>0 | 0.790 | 0.319 | 0.472 | 10.0 | 0.286 | 6 vs 2 |
| CopyKAT-like stromal p95 | %pos UMI>0 | 0.975 | 0.952 | 0.023 | 22.0 | 0.364 | 8 vs 4 |
| EPCAM+KRT (pos) | UCell module | 0.673 | 0.623 | 0.050 | 22.0 | 0.368 | 8 vs 4 |
| CopyKAT-like stromal p95 | UCell module | 0.702 | 0.657 | 0.045 | 20.0 | 0.570 | 8 vs 4 |
| inferCNV-like stromal p95 | UCell module | 0.725 | 0.705 | 0.020 | 19.0 | 0.683 | 8 vs 4 |
| DRMref public labels | %pos UMI>0 | 0.751 | 0.660 | 0.091 | 19.0 | 0.683 | 8 vs 4 |
| EPCAM+KRT (mean) | mean log1p(CP10k) | 1.847 | 1.735 | 0.112 | 18.0 | 0.808 | 8 vs 4 |
| EPCAM+KRT (mean) | %pos UMI>0 | 0.847 | 0.821 | 0.026 | 17.5 | 0.842 | 8 vs 4 |
| inferCNV-like stromal p95 | mean log1p(CP10k) | 1.417 | 1.369 | 0.048 | 17.0 | 0.933 | 8 vs 4 |
| DRMref-like rebuild | mean log1p(CP10k) | 1.646 | 1.537 | 0.109 | 7.0 | 1.000 | 7 vs 2 |

## Full Spearman grid (every ρ / p / n)

Post-treatment patients only. Samples with zero malignant cells under that definition drop out (n < 12).

### Score: mean log1p(CP10k)

| Malignant | T+NK ρ (p, n) | CD8 only ρ (p, n) | CD8+NK ρ (p, n) | CXCL13+ ρ (p, n) | cyto-high T ρ (p, n) |
|---|---|---|---|---|---|
| Hu markers (mean) | **-0.467** (0.205, n=9) | **-0.567** (0.112, n=9) | **-0.483** (0.187, n=9) | **-0.583** (0.099, n=9) | **-0.433** (0.244, n=9) |
| Hu markers (%pos) | -0.238 (0.570, n=8) | **-0.381** (0.352, n=8) | -0.310 (0.456, n=8) | **-0.405** (0.320, n=8) | -0.333 (0.420, n=8) |
| EPCAM+KRT (pos) | -0.119 (0.713, n=12) | -0.245 (0.443, n=12) | -0.189 (0.557, n=12) | -0.294 (0.354, n=12) | -0.084 (0.795, n=12) |
| EPCAM+KRT (mean) | +0.077 (0.812, n=12) | +0.035 (0.914, n=12) | +0.091 (0.779, n=12) | -0.224 (0.484, n=12) | -0.007 (0.983, n=12) |
| CopyKAT-like stromal p95 | +0.077 (0.812, n=12) | -0.028 (0.931, n=12) | +0.014 (0.966, n=12) | -0.133 (0.681, n=12) | +0.063 (0.846, n=12) |
| inferCNV-like stromal p95 | +0.049 (0.880, n=12) | +0.049 (0.880, n=12) | +0.112 (0.729, n=12) | -0.098 (0.762, n=12) | +0.112 (0.729, n=12) |
| DRMref public labels | -0.126 (0.697, n=12) | -0.196 (0.542, n=12) | -0.140 (0.665, n=12) | -0.350 (0.265, n=12) | -0.168 (0.602, n=12) |
| DRMref-like rebuild | **-0.350** (0.356, n=9) | **-0.467** (0.205, n=9) | -0.283 (0.460, n=9) | **-0.467** (0.205, n=9) | **-0.350** (0.356, n=9) |

### Score: %pos UMI>0

| Malignant | T+NK ρ (p, n) | CD8 only ρ (p, n) | CD8+NK ρ (p, n) | CXCL13+ ρ (p, n) | cyto-high T ρ (p, n) |
|---|---|---|---|---|---|
| Hu markers (mean) | **-0.393** (0.295, n=9) | **-0.444** (0.232, n=9) | **-0.477** (0.194, n=9) | **-0.435** (0.242, n=9) | -0.301 (0.431, n=9) |
| Hu markers (%pos) | -0.190 (0.651, n=8) | -0.262 (0.531, n=8) | -0.310 (0.456, n=8) | -0.238 (0.570, n=8) | -0.238 (0.570, n=8) |
| EPCAM+KRT (pos) | **-0.420** (0.175, n=12) | **-0.378** (0.226, n=12) | **-0.371** (0.236, n=12) | **-0.378** (0.226, n=12) | -0.133 (0.681, n=12) |
| EPCAM+KRT (mean) | -0.046 (0.888, n=12) | -0.018 (0.957, n=12) | -0.067 (0.837, n=12) | -0.245 (0.442, n=12) | -0.102 (0.753, n=12) |
| CopyKAT-like stromal p95 | +0.155 (0.631, n=12) | +0.120 (0.711, n=12) | +0.148 (0.646, n=12) | +0.000 (1.000, n=12) | +0.225 (0.481, n=12) |
| inferCNV-like stromal p95 | -0.207 (0.519, n=12) | -0.144 (0.656, n=12) | -0.098 (0.762, n=12) | -0.294 (0.353, n=12) | -0.084 (0.795, n=12) |
| DRMref public labels | -0.189 (0.557, n=12) | -0.175 (0.587, n=12) | -0.175 (0.587, n=12) | -0.343 (0.276, n=12) | -0.091 (0.779, n=12) |
| DRMref-like rebuild | -0.100 (0.798, n=9) | -0.200 (0.606, n=9) | -0.050 (0.898, n=9) | -0.117 (0.765, n=9) | -0.067 (0.865, n=9) |

### Score: UCell module

| Malignant | T+NK ρ (p, n) | CD8 only ρ (p, n) | CD8+NK ρ (p, n) | CXCL13+ ρ (p, n) | cyto-high T ρ (p, n) |
|---|---|---|---|---|---|
| Hu markers (mean) | -0.233 (0.546, n=9) | **-0.400** (0.286, n=9) | -0.283 (0.460, n=9) | **-0.450** (0.224, n=9) | -0.267 (0.488, n=9) |
| Hu markers (%pos) | -0.095 (0.823, n=8) | -0.286 (0.493, n=8) | -0.190 (0.651, n=8) | -0.333 (0.420, n=8) | -0.238 (0.570, n=8) |
| EPCAM+KRT (pos) | -0.315 (0.319, n=12) | **-0.427** (0.167, n=12) | **-0.413** (0.183, n=12) | **-0.657** (0.020, n=12) | **-0.622** (0.031, n=12) |
| EPCAM+KRT (mean) | +0.084 (0.795, n=12) | +0.112 (0.729, n=12) | +0.098 (0.762, n=12) | -0.245 (0.443, n=12) | -0.091 (0.779, n=12) |
| CopyKAT-like stromal p95 | -0.168 (0.602, n=12) | -0.161 (0.618, n=12) | -0.154 (0.633, n=12) | **-0.448** (0.145, n=12) | -0.266 (0.404, n=12) |
| inferCNV-like stromal p95 | -0.035 (0.914, n=12) | -0.014 (0.966, n=12) | -0.028 (0.931, n=12) | **-0.364** (0.245, n=12) | -0.224 (0.484, n=12) |
| DRMref public labels | -0.294 (0.354, n=12) | **-0.441** (0.152, n=12) | **-0.462** (0.131, n=12) | **-0.643** (0.024, n=12) | **-0.622** (0.031, n=12) |
| DRMref-like rebuild | -0.083 (0.831, n=9) | -0.267 (0.488, n=9) | -0.217 (0.576, n=9) | **-0.400** (0.286, n=9) | **-0.350** (0.356, n=9) |

## Full NMPR vs MPR grid (every Δ / p / n)

| Malignant | Score | mean NMPR | mean MPR | Δ | U | p | p one-sided NMPR>MPR | n |
|---|---|---|---|---|---|---|---|---|
| Hu markers (mean) | mean log1p(CP10k) | 2.016 | 0.380 | **1.636** | 14.0 | 0.056 | 0.028 | 7 vs 2 |
| Hu markers (mean) | %pos UMI>0 | 0.849 | 0.361 | **0.488** | 12.0 | 0.194 | 0.083 | 7 vs 2 |
| Hu markers (mean) | UCell module | 0.713 | 0.453 | **0.260** | 14.0 | 0.056 | 0.028 | 7 vs 2 |
| Hu markers (%pos) | mean log1p(CP10k) | 1.988 | 0.362 | **1.627** | 12.0 | 0.071 | 0.036 | 6 vs 2 |
| Hu markers (%pos) | %pos UMI>0 | 0.790 | 0.319 | **0.472** | 10.0 | 0.286 | 0.107 | 6 vs 2 |
| Hu markers (%pos) | UCell module | 0.698 | 0.456 | **0.242** | 12.0 | 0.071 | 0.036 | 6 vs 2 |
| EPCAM+KRT (pos) | mean log1p(CP10k) | 1.612 | 1.258 | **0.354** | 26.0 | 0.109 | 0.059 | 8 vs 4 |
| EPCAM+KRT (pos) | %pos UMI>0 | 0.830 | 0.840 | -0.010 | 15.0 | 0.933 | 0.568 | 8 vs 4 |
| EPCAM+KRT (pos) | UCell module | 0.673 | 0.623 | **0.050** | 22.0 | 0.368 | 0.170 | 8 vs 4 |
| EPCAM+KRT (mean) | mean log1p(CP10k) | 1.847 | 1.735 | **0.112** | 18.0 | 0.808 | 0.366 | 8 vs 4 |
| EPCAM+KRT (mean) | %pos UMI>0 | 0.847 | 0.821 | **0.026** | 17.5 | 0.842 | 0.418 | 8 vs 4 |
| EPCAM+KRT (mean) | UCell module | 0.728 | 0.737 | -0.009 | 14.0 | 0.808 | 0.558 | 8 vs 4 |
| CopyKAT-like stromal p95 | mean log1p(CP10k) | 1.863 | 1.453 | **0.410** | 25.0 | 0.154 | 0.075 | 8 vs 4 |
| CopyKAT-like stromal p95 | %pos UMI>0 | 0.975 | 0.952 | **0.023** | 22.0 | 0.364 | 0.125 | 8 vs 4 |
| CopyKAT-like stromal p95 | UCell module | 0.702 | 0.657 | **0.045** | 20.0 | 0.570 | 0.255 | 8 vs 4 |
| inferCNV-like stromal p95 | mean log1p(CP10k) | 1.417 | 1.369 | **0.048** | 17.0 | 0.933 | 0.438 | 8 vs 4 |
| inferCNV-like stromal p95 | %pos UMI>0 | 0.682 | 0.877 | -0.195 | 9.5 | 0.301 | 0.869 | 8 vs 4 |
| inferCNV-like stromal p95 | UCell module | 0.725 | 0.705 | **0.020** | 19.0 | 0.683 | 0.396 | 8 vs 4 |
| DRMref public labels | mean log1p(CP10k) | 1.571 | 1.120 | **0.450** | 25.0 | 0.154 | 0.055 | 8 vs 4 |
| DRMref public labels | %pos UMI>0 | 0.751 | 0.660 | **0.091** | 19.0 | 0.683 | 0.204 | 8 vs 4 |
| DRMref public labels | UCell module | 0.651 | 0.584 | **0.067** | 24.0 | 0.214 | 0.117 | 8 vs 4 |
| DRMref-like rebuild | mean log1p(CP10k) | 1.646 | 1.537 | **0.109** | 7.0 | 1.000 | 0.500 | 7 vs 2 |
| DRMref-like rebuild | %pos UMI>0 | 0.762 | 0.796 | -0.034 | 6.0 | 0.889 | 0.611 | 7 vs 2 |
| DRMref-like rebuild | UCell module | 0.669 | 0.523 | **0.145** | 14.0 | 0.056 | 0.028 | 7 vs 2 |

## Per-sample malignant n and TACSTD2 (post-treatment)

Malignant cell counts:

| patient | response | hu_markers_mean_n | hu_markers_pctpos_n | epcam_krt_pos_n | epcam_krt_mean_n | copykat_stromal_p95_n | infercnv_stromal_p95_n | drmref_public_n | drmref_like_n |
|---|---|---|---|---|---|---|---|---|---|
| P02 | NMPR | 0 | 0 | 92 | 16 | 13 | 16 | 100 | 2 |
| P03 | MPR | 655 | 392 | 978 | 34 | 494 | 194 | 440 | 372 |
| P04 | NMPR | 33 | 30 | 136 | 35 | 35 | 25 | 71 | 51 |
| P06 | MPR | 0 | 0 | 66 | 4 | 19 | 5 | 15 | 0 |
| P07 | NMPR | 2996 | 2270 | 5315 | 3002 | 3609 | 484 | 516 | 4014 |
| P09 | NMPR | 205 | 135 | 237 | 139 | 196 | 116 | 69 | 236 |
| P10 | NMPR | 186 | 143 | 168 | 12 | 75 | 35 | 216 | 144 |
| P11 | MPR | 0 | 0 | 433 | 23 | 274 | 118 | 54 | 0 |
| P12 | NMPR | 362 | 267 | 409 | 212 | 90 | 97 | 181 | 383 |
| P13 | NMPR | 1 | 0 | 35 | 17 | 21 | 3 | 25 | 0 |
| P14 | MPR | 1 | 1 | 162 | 6 | 118 | 51 | 87 | 2 |
| P15 | NMPR | 1 | 1 | 486 | 121 | 169 | 32 | 277 | 4 |

Immune fractions:

| patient | response | frac_T_NK | frac_CD8_only | frac_CD8_NK | frac_CXCL13_pos | frac_cyto_high_T |
|---|---|---|---|---|---|---|
| P02 | NMPR | 0.542 | 0.169 | 0.188 | 0.037 | 0.017 |
| P03 | MPR | 0.546 | 0.292 | 0.331 | 0.247 | 0.220 |
| P04 | NMPR | 0.669 | 0.282 | 0.331 | 0.123 | 0.131 |
| P06 | MPR | 0.629 | 0.389 | 0.464 | 0.163 | 0.134 |
| P07 | NMPR | 0.119 | 0.039 | 0.056 | 0.013 | 0.015 |
| P09 | NMPR | 0.677 | 0.359 | 0.401 | 0.339 | 0.327 |
| P10 | NMPR | 0.614 | 0.345 | 0.426 | 0.212 | 0.268 |
| P11 | MPR | 0.450 | 0.161 | 0.243 | 0.030 | 0.050 |
| P12 | NMPR | 0.222 | 0.099 | 0.115 | 0.063 | 0.045 |
| P13 | NMPR | 0.199 | 0.065 | 0.076 | 0.014 | 0.043 |
| P14 | MPR | 0.516 | 0.155 | 0.196 | 0.046 | 0.030 |
| P15 | NMPR | 0.189 | 0.054 | 0.110 | 0.005 | 0.017 |

## Notes

- Author CopyKAT IDs remain unpublished. CopyKAT-like and inferCNV-like rows are window-smoothed expression CNV vs stromal cells on the public UMI.
- DRMref public labels are third-party marker annotations (30,877 cells), not Hu CopyKAT. DRMref-like is a 12-type marker rebuild.
- n=12 is small. A true ρ = −0.45 has two-sided Spearman p ≈ 0.14 at n=12. Every p is reported as computed.
- Highlight threshold ρ ≤ −0.35 is the requested flag for this extra figure, not a new claim cutoff.

