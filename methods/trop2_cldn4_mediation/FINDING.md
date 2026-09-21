# TROP2–immune and CLDN4

**Call: null.** An immune-cold TACSTD2 (TROP2) association that depends partly on CLDN4 is not supported. Of the twelve primary rows, seven are null, four are partial, and one is opposite. None is support.

The locked CLDN4 results stand. Concordant-4 malignant CLDN4 percent detected versus the T/NK fraction is ρ = −0.531 (p = 1.6×10⁻⁵, I² = 0%, N = 65). Given TACSTD2, that partial is ρ = −0.543 (p = 2.0×10⁻⁵). On CosMx, the section-mean Spearman of continuous CLDN4 versus the CD8+NK neighbor count at 50 µm is ρ = −0.037, the same figure as the prior continuous sweep.

## Evidence matrix

| Dataset | Immune readout | n | TACSTD2–immune | CLDN4–immune | TACSTD2 given CLDN4 | CLDN4 given TACSTD2 | % attenuation | ACME (95% interval) | Call |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| He2022 CosMx | CD8+NK neighbors, 50 µm | 5 donors | -0.047 (p=0.004) | -0.054 (p=0.086) | -0.037 (p<0.001) | -0.045 (p=0.144) | — | -0.010 (-0.023, +0.004) | null |
| He2022 CosMx | CD8+NK neighbors, 100 µm | 5 donors | -0.039 (p=0.073) | -0.048 (p=0.203) | -0.029 (p=0.011) | -0.040 (p=0.270) | — | -0.007 (-0.025, +0.009) | null |
| GSE123902 | malignant % detected | 13 | -0.324 (p=0.280) | -0.659 (p=0.014) | -0.174 (p=0.588) | -0.623 (p=0.031) | +46% | -0.226 (-0.528, +0.902) | partial |
| GSE123902 | malignant mean log1p(UMI) | 13 | -0.209 (p=0.494) | -0.654 (p=0.015) | -0.020 (p=0.950) | -0.634 (p=0.027) | +90% | -0.109 (-0.738, +0.420) | partial |
| GSE131907 | malignant % detected | 21 | -0.151 (p=0.515) | -0.522 (p=0.015) | +0.316 (p=0.174) | -0.575 (p=0.008) | +310% | -0.409 (-1.581, +0.009) | null |
| GSE131907 | malignant mean log1p(UMI) | 21 | +0.082 (p=0.724) | -0.396 (p=0.075) | +0.422 (p=0.064) | -0.550 (p=0.012) | -415% | -0.356 (-0.779, -0.042) | null |
| GSE205335 | malignant % detected | 22 | +0.208 (p=0.352) | -0.435 (p=0.043) | +0.345 (p=0.126) | -0.503 (p=0.020) | -65% | -0.016 (-0.288, +0.255) | null |
| GSE205335 | malignant mean log1p(UMI) | 22 | +0.292 (p=0.187) | -0.170 (p=0.450) | +0.402 (p=0.071) | -0.332 (p=0.142) | -38% | -0.100 (-0.376, +0.291) | null |
| GSE189357 | malignant % detected | 9 | -0.483 (p=0.187) | -0.600 (p=0.088) | +0.063 (p=0.882) | -0.410 (p=0.313) | +113% | -0.988 (-8.898, -0.271) | partial |
| GSE189357 | malignant mean log1p(UMI) | 9 | -0.400 (p=0.286) | -0.517 (p=0.154) | +0.112 (p=0.792) | -0.372 (p=0.364) | +128% | -0.916 (-4.746, +2.846) | partial |
| concordant-4 | malignant % detected | 65 (4 cohorts) | -0.112 (p=0.460) | -0.531 (p=1.6e-05) | +0.216 (p=0.124) | -0.543 (p=2.0e-05) | +294% | -0.076 (-0.311, +0.159) | null |
| concordant-4 | malignant mean log1p(UMI) | 65 (4 cohorts) | +0.041 (p=0.777) | -0.394 (p=0.002) | +0.310 (p=0.025) | -0.476 (p=2.9e-04) | — | -0.204 (-0.423, +0.016) | opposite |

Coefficients are Spearman ρ. CosMx values are the mean of five within-donor coefficients; the p-value and the ACME interval resample donors. Concordant-4 cohort rows are within-cohort tests. The concordant-4 row is DerSimonian–Laird. ACME is outcome standard deviations per TACSTD2 standard deviation. Attenuation is blank when |raw ρ| < 0.05. Attenuation above 100%, and negative attenuation, are printed and are not mediated proportions.

Figure: `results/figures/evidence_matrix.png`.

## Reading the rows

On CosMx the within-donor correlations are about −0.04 at both radii. At 50 µm the TACSTD2 mean is detectably negative (ρ = −0.047; donor-bootstrap interval −0.101 to −0.006) and the ACME interval includes 0. TACSTD2 is negative in 5/8 sections and 4/5 donors. That is too small for CLDN4 to be carrying a TROP2–immune association.

On the locked concordant-4 percent-detected scale, TACSTD2 versus T/NK is ρ = −0.112 (p = 0.46; 95% CI −0.388 to +0.183). Given CLDN4 the partial is +0.216 (p = 0.12). ACME = −0.076 (−0.311 to +0.159). The +294% figure is the partial crossing zero on a null total effect.

The mean log1p meta-analysis is the opposite call. Raw TACSTD2 ρ = +0.041. Given CLDN4, ρ = +0.310 (p = 0.025; 95% CI +0.040 to +0.537): higher TACSTD2 goes with a higher T/NK fraction once CLDN4 is held fixed.

GSE123902 (n = 13) is partial on both readouts. The point estimates shrink by 46% and 90%, the ACME is negative, and every interval includes 0. GSE189357 (n = 9) is partial on the point estimates. Its percent-detected ACME interval runs from −8.90 to −0.27. Leaving one patient out keeps the ACME negative and moves it between −0.75 and −5.64. GSE131907 mean log1p has an ACME interval entirely below 0 (−0.779 to −0.042) while raw TACSTD2 ρ = +0.082, so the indirect path and the total association disagree. GSE205335 TACSTD2 point estimates are positive (+0.208 and +0.292), with intervals that include 0.

Stacking all 65 patients without a cohort term, percent-detected TACSTD2 ρ = −0.229 (p = 0.067) and the ACME interval lies below 0, while the partial flips to +0.171. The prespecified within-cohort pool does not show that pattern.

Across the eight CosMx sections, attenuation of the section-mean TACSTD2 Spearman is 3% at 50 µm and 5% at 100 µm. Section-mean CLDN4 versus the 100 µm neighbor count is ρ = +0.738 (p = 0.037). Leaving out LUSC-6, LUAD-9 R1, or LUAD-13 moves that p-value to 0.15 (`results/tables/cosmx_between_section_loo.tsv`).

The locked CosMx high-versus-low cytotoxic ratios, 0.36 at 50 µm and 0.52 at 100 µm, remain the exclusion summary. This matrix is the continuous partial correlation and the product-method ACME. Definitions and the call rule are in `METHODS.md`.
