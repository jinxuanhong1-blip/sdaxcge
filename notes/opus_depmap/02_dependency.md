# 02 · Dependency of TACSTD2 and CLDN4 in lung lines

Release DepMap Public 24Q4. Lung models with CRISPR screens: **126** (NSCLC 98, SCLC 25); non-lung comparison set: 1052 models.

Chronos scale reminder: 0 = median non-essential gene, -1 = median common essential gene. In these lung lines the mean gene effect of the 726 non-essential control genes is -0.0160 with a between-gene SD of 0.1162; the 1242 common essential controls average -1.077. Those two numbers set the scale for everything below.

## Headline numbers

| gene | lung mean (SD) | lung min | dependent lung lines (p>0.5) | dependent lines pan-cancer | mean in non-essential SD units |
| --- | --- | --- | --- | --- | --- |
| TACSTD2 | -0.0287 (0.121) | -0.408 | 0 / 126 | 1 / 1178 | -0.11 |
| CLDN4 | 0.0342 (0.132) | -0.294 | 0 / 126 | 0 / 1178 | +0.43 |

Full statistics, including the benchmark genes, are in `results/opus_depmap/tables/dependency_summary.csv`; genome-wide calibration is in `dependency_genomewide_calibration.csv`; the expression-versus-own-dependency test is in `dependency_vs_expression.csv`.

Figure: `results/opus_depmap/figures/fig1_dependency_distributions.png`.
