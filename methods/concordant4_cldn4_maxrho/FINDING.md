# Concordant-4 CLDN4 vs T/NK — cut / %pos / keratin / histology sweep

Exploratory maximum on the locked four cohorts only: GSE123902, GSE131907, GSE205335, GSE189357.
No GSE148071, GSE127465, GSE207422, GSE154826, GSE200563, or E-MTAB-13526.
Patient / donor / sample is the unit. p-values are descriptive. This maximum is the best row of a pre-specified grid, not a confirmatory p-value.

## Rule

A panel is 4-cohort consistent when all four cohorts remain, each Spearman is negative, and the stacked within-cohort Q4 vs Q1 Cliff delta on raw T/NK fraction is negative.
Unadjusted panels need n≥6 per cohort. A one-keratin partial needs n≥7. The three-keratin joint partial needs n≥9.
Pooled ρ is DerSimonian–Laird on Fisher z. Partials use variance 1/(n−3−k).
Best |ρ| below is the sign-consistent maximum. A second row requires I²=0, which is the locked pool's consistency standard. Joint ranking is |ρ|+|Cliff δ|.

## Locked reproduction

Locked published panel is malignant CLDN4 %pos (UMI>0), no keratin adjustment, all histologies, n_mal≥20.
Recomputed: N=65, ρ=-0.531 (p=1.65e-05, I²=0.0%, -0.697 to -0.312), Cliff δ=-0.724 (n_Q1/n_Q4=19/16, p=2.88e-04).
The published point was ρ=−0.531. Cohort ρ were −0.659 / −0.522 / −0.435 / −0.600.

Recomputed cohort ρ: GSE123902 -0.659 (n=13), GSE131907 -0.522 (n=21), GSE205335 -0.435 (n=22), GSE189357 -0.600 (n=9).

The published mean row was ρ=−0.403. Recomputed with the scale actually stored in each locked table (log1p of raw UMI in GSE123902, GSE131907, and GSE189357; log1p CP10K in GSE205335): ρ=-0.403 (p=0.002, I²=0.0%, Cliff δ=-0.592).
That mixed scale is a reproduction check only. The search below uses one scale in all four cohorts.

## Best |ρ| panel

Grid size 6435. Sign-consistent panels 1100.
Rows with both a larger |ρ| and a larger |Cliff δ| than the locked %pos panel: 0.

The joint maximum is the locked panel. Malignant CLDN4 %pos (UMI>0), no keratin adjustment, all histologies, n_mal≥20: ρ=-0.531, Cliff δ=-0.724, I²=0.0%, N=65.

Largest |ρ| with I²=0:

`pct_gt0|partial|KRT18_cp10k|all|min20`. N=65, ρ=-0.533 (p=3.21e-05, I²=0.0%, -0.703 to -0.304), Cliff δ=-0.697 (n_Q1/n_Q4=19/16, p=4.77e-04).
Cohort ρ: GSE123902 -0.547 (n=13), GSE131907 -0.488 (n=21), GSE205335 -0.506 (n=22), GSE189357 -0.712 (n=9).
The |ρ| gain over the locked panel is in the third decimal. |Cliff δ| is smaller.

Largest |ρ| with only the sign constraint (I² not required). This pool is heterogeneous:

`pct_ge2|partial|KRT8_raw|all|min100`. N=60, ρ=-0.552 (p=0.023, I²=66.5%, -0.821 to -0.084), Cliff δ=-0.498 (n_Q1/n_Q4=17/15, p=0.017).

| cohort | n | ρ | p |
|---|---:|---:|---:|
| GSE123902 | 11 | -0.107 | 0.754 |
| GSE131907 | 19 | -0.491 | 0.033 |
| GSE205335 | 21 | -0.357 | 0.112 |
| GSE189357 | 9 | -0.940 | 1.64e-04 |

Best mean or upper-percentile score in the sign-consistent grid: `mean_log1p_cp10k|partial|KRT5_cp10k|adc|min20` ρ=-0.480, Cliff δ=-0.446, I²=0.0%, N=52.
Best row whose histology filter is not the full set: `pct_gt0|none|none|drop_ais|min20` ρ=-0.527, Cliff δ=-0.722, I²=0.0%, N=62.
Neither is larger in |ρ| than the locked %pos panel.

## Cliff and joint maxima

Largest |Cliff δ| is a different panel: `pct_gt0|partial|KRT_SIMPLE_raw|all|min100` ρ=-0.491, Cliff δ=-0.757, N=60.
Largest |ρ|+|Cliff δ| is `pct_gt0|none|none|all|min20` ρ=-0.531, Cliff δ=-0.724, N=65.

## Top consistent panels by |ρ|

| panel | N | ρ | I² | Cliff δ |
|---|---:|---:|---:|---:|
| `pct_ge2|partial|KRT8_raw|all|min100` | 60 | -0.552 | 66.5% | -0.498 |
| `pct_ge3|partial|KRT8_raw|all|min100` | 60 | -0.542 | 68.4% | -0.576 |
| `pct_ge2|semipartial|KRT8_raw|all|min100` | 60 | -0.542 | 66.3% | -0.498 |
| `pct_gt0|partial|KRT18_cp10k|all|min20` | 65 | -0.533 | 0.0% | -0.697 |
| `pct_ge3|semipartial|KRT8_raw|all|min100` | 60 | -0.532 | 67.9% | -0.576 |
| `pct_ge2|partial|KRT8_raw|all|min20` | 65 | -0.531 | 67.8% | -0.513 |
| `pct_gt0|none|none|all|min20` | 65 | -0.531 | 0.0% | -0.724 |
| `pct_gt0|none|none|drop_ais|min20` | 62 | -0.527 | 0.0% | -0.722 |
| `pct_ge3|partial|KRT8_raw|all|min20` | 65 | -0.524 | 68.1% | -0.441 |
| `pct_gt0|partial|KRT5_raw|all|min20` | 60 | -0.524 | 0.0% | -0.686 |
| `pct_gt0|partial|KRT5_cp10k|all|min20` | 60 | -0.523 | 0.0% | -0.686 |
| `pct_ge2|partial|KRT8_raw|all|min50` | 63 | -0.523 | 69.0% | -0.474 |

## What the grid was

Cuts: malignant CLDN4 % of cells with UMI>0, ≥2, ≥3, ≥5, ≥10, above the cohort median, and above the cohort 75th percentile. Inclusion floors n_malignant ≥20, ≥50, ≥100.
Means: mean log1p(UMI), mean log1p(CP10K), mean among detected cells, and the 90th percentile, each on raw UMI and on CP10K.
Keratin partials: Spearman partial and semipartial on KRT8, KRT18, KRT19, KRT7, KRT5, KRT6A, the mean of KRT8/18/19, and the joint KRT8+KRT18+KRT19 residual. Each keratin covariate is fit on raw log1p and on log1p CP10K. Quartiles for Cliff δ use the residual CLDN4 score; the outcome stays raw T/NK fraction.
Histology: all units; GSE205335 ADC only; GSE205335 ADC+SQ; drop GSE189357 AIS; ADC-only plus drop AIS (invasive adenocarcinoma spectrum). GSE123902 and GSE131907 are lung adenocarcinoma in the GEO records, so those filters do not drop them. GSE189357 histological type is the GEO field (AIS / MIA / IAC).

## Caveats

GSE123902 and GSE189357 malignant cells are marker-gated (EPCAM or KRT8/18/19, and PTPRC==0), so a keratin covariate is not independent of the gate. GSE131907 and GSE205335 use author malignant labels.
The searched |ρ| is the maximum of the grid. It is not a replacement of the pre-specified locked %pos estimate unless the locked row itself is the maximum.
Not causal. Not a cell-level correlation. Not a spatial exclusion claim.

Reproduce: `python3 methods/concordant4_cldn4_maxrho/scripts/extract_scores.py` then `python3 methods/concordant4_cldn4_maxrho/scripts/sweep.py`. GEO matrices are read from `/tmp/geo_dl` and are not stored in this repo.
