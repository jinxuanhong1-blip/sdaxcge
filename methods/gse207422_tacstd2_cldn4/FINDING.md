# GSE207422: TACSTD2 and CLDN4, scored separately

Public neoadjuvant PD-1 plus chemotherapy NSCLC scRNA-seq (Hu et al., Genome Medicine 2023, PMID 36869384).
Slide 5–6 use this cohort for a TROP2-high state: higher tumor fraction (stated P=0.047), lower CD8 T/NK (stated P=0.015), and lower TROP2 in MPR than in NMPR.
This note scores TACSTD2, CLDN4, and dual-high with the same rules. CLDN4 is not required to match TACSTD2.

## Data

- GEO processed UMI matrix re-downloaded (184,001,817 bytes, sha256 `aba15960fc7ee6a2443511bce5177e4d71b964131b6e98597e7a85d0a213ba36`).
- 24,292 genes × 92,330 cells. The slide text says 90,652 cells; the deposited matrix has 92,330, which matches the paper. No cells were dropped to force the slide count.
- Median genes/cell = 1253 (paper 1,256).
- Sample metadata: 15 patients, one sample each. 12 post-treatment surgery samples (MPR n=4, including pCR P06; NMPR n=8) and 3 pre-treatment biopsies.
- GEO does not deposit author barcode labels or CopyKAT calls. Lineages are reconstructed from canonical markers. TACSTD2 and CLDN4 are not used to call a cell.
- Epithelial: highest lineage score is epithelial. A3-malignant: epithelial and zero UMI for SFTPA2, AGER, SCGB1A1, SCGB3A1, TPPP3.
- A3 leaves fewer than 10 malignant cells in 5 of 12 post-treatment samples (P06, P11, P13, P14, P15), including 3 of 4 MPR samples (P06 pCR, P11, P14). An MPR contrast on the A3 subset is 1 MPR vs 6 NMPR. The cohort that still contains every post-treatment sample is the epithelial mean (n=12).

## Where each gene is expressed (post-treatment)

Median percent of cells with UMI > 0, across post-treatment samples that have at least 5 cells in the compartment.

- TACSTD2 in A3-malignant: 66.9% (n samples=7)
- TACSTD2 in epithelial: 80.7% (n samples=12)
- TACSTD2 in CD8/NK: 1.3% (n samples=12)
- TACSTD2 in T/NK lineage: 0.5% (n samples=12)
- CLDN4 in A3-malignant: 70.0% (n samples=7)
- CLDN4 in epithelial: 85.5% (n samples=12)
- CLDN4 in CD8/NK: 1.4% (n samples=12)
- CLDN4 in T/NK lineage: 0.7% (n samples=12)

Both genes are epithelial-restricted (about 1% of CD8/NK cells positive, under 1% of T/NK-lineage cells). A whole-sample average therefore moves with epithelial fraction. That all-cell contrast is reported below and is not an independent test of tumor content. The expression values in the main tests are means inside epithelial cells, or inside A3-malignant cells.

## Post-treatment samples that still contain the MPR cases (epithelial mean, n=12)

Median cuts on the epithelial mean log1p(CP10k): TACSTD2 = 1.550, CLDN4 = 1.529. High = at or above that median. Dual-high = both.

- **TACSTD2 high vs low.** epithelial fraction: high median 0.043 (n=6) vs low median 0.055 (n=6), Δ=-0.012, exact P=0.699; CD8 T + NK fraction: high median 0.222 (n=6) vs low median 0.197 (n=6), Δ=0.024, exact P=0.699; all T + NK fraction: high median 0.403 (n=6) vs low median 0.527 (n=6), Δ=-0.124, exact P=0.937; MPR rate: high 0/6 vs low 4/6, Fisher P=0.061.
- **CLDN4 high vs low.** epithelial fraction: high median 0.038 (n=6) vs low median 0.050 (n=6), Δ=-0.013, exact P=0.818; CD8 T + NK fraction: high median 0.197 (n=6) vs low median 0.249 (n=6), Δ=-0.052, exact P=0.818; all T + NK fraction: high median 0.494 (n=6) vs low median 0.528 (n=6), Δ=-0.034, exact P=0.937; MPR rate: high 2/6 vs low 2/6, Fisher P=1.000.
- **Dual-high vs not** (n dual-high = 3). epithelial fraction: high median 0.056 (n=3) vs low median 0.034 (n=9), Δ=0.022, exact P=1.000; CD8 T + NK fraction: high median 0.078 (n=3) vs low median 0.222 (n=9), Δ=-0.143, exact P=0.282; all T + NK fraction: high median 0.199 (n=3) vs low median 0.540 (n=9), Δ=-0.341, exact P=0.600; MPR rate: high 0/3 vs low 4/9, Fisher P=0.491.

Continuous Spearman, same 12 samples:

- TACSTD2: epithelial fraction ρ=-0.07 P=0.829 (n=12); CD8 T + NK fraction ρ=-0.18 P=0.572 (n=12); all T + NK fraction ρ=-0.10 P=0.762 (n=12); residual tumor fraction ρ=0.48 P=0.117 (n=12).
- CLDN4: epithelial fraction ρ=-0.06 P=0.846 (n=12); CD8 T + NK fraction ρ=-0.10 P=0.746 (n=12); all T + NK fraction ρ=0.03 P=0.914 (n=12); residual tumor fraction ρ=0.12 P=0.720 (n=12).

Expression itself, NMPR versus MPR (pCR counted as MPR). A positive NMPR−MPR difference is higher expression in NMPR, which is the direction stated for TROP2 on the slide.

- TACSTD2: NMPR median 1.695 (n=8) vs MPR median 1.329 (n=4), Δ(NMPR−MPR)=0.366, exact P=0.109.
- CLDN4: NMPR median 1.577 (n=8) vs MPR median 1.515 (n=4), Δ(NMPR−MPR)=0.062, exact P=1.000.

The TACSTD2 MPR values sit inside the NMPR range (see `figures/fig_epithelial_mpr.png`). The median split puts all 4 MPR samples on the TACSTD2-low side; the continuous test still overlaps. CLDN4 MPR values cover the same span as NMPR.

Group membership:

- TACSTD2 high: n=6, MPR=0, NMPR=6, Adeno=2, Squamous=4 (P04,P07,P09,P10,P13,P15)
- TACSTD2 low: n=6, MPR=4, NMPR=2, Adeno=4, Squamous=2 (P02,P03,P06,P11,P12,P14)
- CLDN4 high: n=6, MPR=2, NMPR=4, Adeno=3, Squamous=3 (P02,P06,P07,P09,P11,P13)
- CLDN4 low: n=6, MPR=2, NMPR=4, Adeno=3, Squamous=3 (P03,P04,P10,P12,P14,P15)
- dual-high not dual-high: n=9, MPR=4, NMPR=5, Adeno=6, Squamous=3 (P02,P03,P04,P06,P10,P11,P12,P14,P15)
- dual-high dual-high: n=3, MPR=0, NMPR=3, Adeno=0, Squamous=3 (P07,P09,P13)

Whole-sample mean versus epithelial fraction, same 12 samples: TACSTD2 ρ=0.91 (P=4.19e-05); CLDN4 ρ=0.97 (P=3.88e-07). That correlation is the epithelial restriction above. It is not an independent test of tumor content, and it is not used as the tumor-fraction result.

## A3-malignant subset (not an MPR test)

Samples with ≥10 A3-malignant cells: n=7. Cuts: TACSTD2 = 1.473, CLDN4 = 1.354. One MPR sample remains (P03).

- **TACSTD2.** A3-malignant fraction: high median 0.028 (n=4) vs low median 0.053 (n=3), Δ=-0.025, exact P=1.000; CD8 T + NK fraction: high median 0.377 (n=4) vs low median 0.172 (n=3), Δ=0.205, exact P=0.400; all T + NK fraction: high median 0.636 (n=4) vs low median 0.540 (n=3), Δ=0.097, exact P=0.400; MPR rate: high 0/4 vs low 1/3, Fisher P=0.429.
- **CLDN4.** A3-malignant fraction: high median 0.019 (n=4) vs low median 0.053 (n=3), Δ=-0.034, exact P=0.629; CD8 T + NK fraction: high median 0.252 (n=4) vs low median 0.325 (n=3), Δ=-0.073, exact P=0.857; all T + NK fraction: high median 0.603 (n=4) vs low median 0.542 (n=3), Δ=0.061, exact P=0.857; MPR rate: high 0/4 vs low 1/3, Fisher P=0.429.
- **Dual-high** (n=3). A3-malignant fraction: high median 0.033 (n=3) vs low median 0.039 (n=4), Δ=-0.006, exact P=0.857; CD8 T + NK fraction: high median 0.332 (n=3) vs low median 0.248 (n=4), Δ=0.084, exact P=1.000; all T + NK fraction: high median 0.666 (n=3) vs low median 0.541 (n=4), Δ=0.125, exact P=0.629; MPR rate: high 0/3 vs low 1/4, Fisher P=1.000.

- TACSTD2 expression: NMPR median 1.653 (n=6) vs MPR median 0.734 (n=1), Δ(NMPR−MPR)=0.919, exact P=0.571.
- CLDN4 expression: NMPR median 1.527 (n=6) vs MPR median 1.022 (n=1), Δ(NMPR−MPR)=0.506, exact P=0.857.

Spearman on this subset:

- TACSTD2: A3-malignant fraction ρ=0.32 P=0.482 (n=7); CD8 T + NK fraction ρ=0.00 P=1.000 (n=7); all T + NK fraction ρ=0.11 P=0.819 (n=7); residual tumor fraction ρ=0.32 P=0.478 (n=7).
- CLDN4: A3-malignant fraction ρ=-0.32 P=0.482 (n=7); CD8 T + NK fraction ρ=-0.07 P=0.879 (n=7); all T + NK fraction ρ=0.14 P=0.760 (n=7); residual tumor fraction ρ=0.47 P=0.289 (n=7).

Membership:

- TACSTD2 high: n=4, MPR=0, NMPR=4, Adeno=1, Squamous=3 (P04,P07,P09,P10)
- TACSTD2 low: n=3, MPR=1, NMPR=2, Adeno=2, Squamous=1 (P02,P03,P12)
- CLDN4 high: n=4, MPR=0, NMPR=4, Adeno=2, Squamous=2 (P02,P04,P07,P09)
- CLDN4 low: n=3, MPR=1, NMPR=2, Adeno=1, Squamous=2 (P03,P10,P12)
- dual-high not dual-high: n=4, MPR=1, NMPR=3, Adeno=2, Squamous=2 (P02,P03,P10,P12)
- dual-high dual-high: n=3, MPR=0, NMPR=3, Adeno=1, Squamous=2 (P04,P07,P09)

## How the two genes compare

On all 12 post-treatment samples, a higher epithelial TACSTD2 does not bring a higher epithelial fraction (high−low Δ=-0.012, P=0.699) or a lower CD8 T/NK fraction (Δ=0.024, P=0.699; Spearman ρ=-0.18, P=0.572). The same is true for CLDN4 (tumor-fraction Δ=-0.013, P=0.818; CD8/NK Δ=-0.052, P=0.818; Spearman ρ=-0.10, P=0.746). CLDN4 versus the lineage T/NK fraction is ρ=0.03 (P=0.914); TACSTD2 versus that fraction is ρ=-0.10 (P=0.762). The MPR contrast is where the genes differ. Epithelial TACSTD2 is higher in NMPR than in MPR (Δ=0.366, P=0.109), and the TACSTD2-high half contains no MPR sample. Epithelial CLDN4 does not shift between NMPR and MPR (Δ=0.062, P=1.000), and the CLDN4 median split is 2 MPR versus 2 MPR. Dual-high is reported as its own row; it is not used to rewrite the CLDN4 estimate.

## Other cuts

Same tests under other pre-specified definitions. Not used to pick a CLDN4 result.

- response_known_ge10_A3 | TACSTD2 | A3_malignant_mean_log1p_cp10k | A3-malignant fraction: high 0.033 vs low 0.028 (n=5 vs 4), Δ=0.004, P=0.413
- response_known_ge10_A3 | TACSTD2 | A3_malignant_mean_log1p_cp10k | CD8 T + NK fraction: high 0.332 vs low 0.146 (n=5 vs 4), Δ=0.186, P=0.730
- response_known_ge10_A3 | TACSTD2 | A3_malignant_mean_log1p_cp10k | all T + NK fraction: high 0.607 vs low 0.541 (n=5 vs 4), Δ=0.065, P=0.556
- response_known_ge10_A3 | TACSTD2 | A3_malignant_mean_log1p_cp10k | MPR rate high 0/5 vs low 1/4, Fisher P=0.444
- response_known_ge10_A3 | CLDN4 | A3_malignant_mean_log1p_cp10k | A3-malignant fraction: high 0.006 vs low 0.056 (n=5 vs 4), Δ=-0.049, P=0.190
- response_known_ge10_A3 | CLDN4 | A3_malignant_mean_log1p_cp10k | CD8 T + NK fraction: high 0.172 vs low 0.223 (n=5 vs 4), Δ=-0.051, P=0.905
- response_known_ge10_A3 | CLDN4 | A3_malignant_mean_log1p_cp10k | all T + NK fraction: high 0.543 vs low 0.402 (n=5 vs 4), Δ=0.141, P=0.556
- response_known_ge10_A3 | CLDN4 | A3_malignant_mean_log1p_cp10k | MPR rate high 0/5 vs low 1/4, Fisher P=0.444
- response_known_ge10_A3 | dual-high | A3_malignant_mean_log1p_cp10k | A3-malignant fraction: high 0.033 vs low 0.039 (n=3 vs 6), Δ=-0.006, P=0.905
- response_known_ge10_A3 | dual-high | A3_malignant_mean_log1p_cp10k | CD8 T + NK fraction: high 0.332 vs low 0.146 (n=3 vs 6), Δ=0.186, P=0.905
- response_known_ge10_A3 | dual-high | A3_malignant_mean_log1p_cp10k | all T + NK fraction: high 0.666 vs low 0.541 (n=3 vs 6), Δ=0.125, P=0.548
- response_known_ge10_A3 | dual-high | A3_malignant_mean_log1p_cp10k | MPR rate high 0/3 vs low 1/6, Fisher P=1.000
- response_known_ge10_A3 | TACSTD2 | A3_malignant_mean_log1p_cp10k | NMPR vs MPR expression: medians 1.653 vs 0.734 (n=8 vs 1), P=0.444
- response_known_ge10_A3 | CLDN4 | A3_malignant_mean_log1p_cp10k | NMPR vs MPR expression: medians 1.389 vs 1.022 (n=8 vs 1), P=0.889
- post_ge10_A3_PTPRCneg | TACSTD2 | A3_PTPRCneg_mean_log1p_cp10k | PTPRC-negative malignant fraction: high 0.021 vs low 0.047 (n=3 vs 3), Δ=-0.026, P=1.000
- post_ge10_A3_PTPRCneg | TACSTD2 | A3_PTPRCneg_mean_log1p_cp10k | CD8 T + NK fraction: high 0.332 vs low 0.325 (n=3 vs 3), Δ=0.007, P=1.000
- post_ge10_A3_PTPRCneg | TACSTD2 | A3_PTPRCneg_mean_log1p_cp10k | all T + NK fraction: high 0.666 vs low 0.542 (n=3 vs 3), Δ=0.124, P=0.700
- post_ge10_A3_PTPRCneg | TACSTD2 | A3_PTPRCneg_mean_log1p_cp10k | MPR rate high 0/3 vs low 1/3, Fisher P=1.000
- post_ge10_A3_PTPRCneg | CLDN4 | A3_PTPRCneg_mean_log1p_cp10k | PTPRC-negative malignant fraction: high 0.021 vs low 0.047 (n=3 vs 3), Δ=-0.026, P=1.000
- post_ge10_A3_PTPRCneg | CLDN4 | A3_PTPRCneg_mean_log1p_cp10k | CD8 T + NK fraction: high 0.332 vs low 0.325 (n=3 vs 3), Δ=0.007, P=1.000
- post_ge10_A3_PTPRCneg | CLDN4 | A3_PTPRCneg_mean_log1p_cp10k | all T + NK fraction: high 0.666 vs low 0.542 (n=3 vs 3), Δ=0.124, P=0.700
- post_ge10_A3_PTPRCneg | CLDN4 | A3_PTPRCneg_mean_log1p_cp10k | MPR rate high 0/3 vs low 1/3, Fisher P=1.000
- post_ge10_A3_PTPRCneg | dual-high | A3_PTPRCneg_mean_log1p_cp10k | PTPRC-negative malignant fraction: high 0.021 vs low 0.047 (n=3 vs 3), Δ=-0.026, P=1.000
- post_ge10_A3_PTPRCneg | dual-high | A3_PTPRCneg_mean_log1p_cp10k | CD8 T + NK fraction: high 0.332 vs low 0.325 (n=3 vs 3), Δ=0.007, P=1.000
- post_ge10_A3_PTPRCneg | dual-high | A3_PTPRCneg_mean_log1p_cp10k | all T + NK fraction: high 0.666 vs low 0.542 (n=3 vs 3), Δ=0.124, P=0.700
- post_ge10_A3_PTPRCneg | dual-high | A3_PTPRCneg_mean_log1p_cp10k | MPR rate high 0/3 vs low 1/3, Fisher P=1.000
- post_ge10_A3_PTPRCneg | TACSTD2 | A3_PTPRCneg_mean_log1p_cp10k | NMPR vs MPR expression: medians 1.724 vs 0.708 (n=5 vs 1), P=0.333
- post_ge10_A3_PTPRCneg | CLDN4 | A3_PTPRCneg_mean_log1p_cp10k | NMPR vs MPR expression: medians 1.236 vs 0.995 (n=5 vs 1), P=1.000
- post_all12_allcell_mean | TACSTD2 | all_cell_mean_log1p_cp10k_NOT_independent_of_tumor_fraction | A3-malignant fraction: high 0.028 vs low 0.001 (n=6 vs 6), Δ=0.027, P=0.394
- post_all12_allcell_mean | TACSTD2 | all_cell_mean_log1p_cp10k_NOT_independent_of_tumor_fraction | CD8 T + NK fraction: high 0.273 vs low 0.172 (n=6 vs 6), Δ=0.101, P=0.937
- post_all12_allcell_mean | TACSTD2 | all_cell_mean_log1p_cp10k_NOT_independent_of_tumor_fraction | all T + NK fraction: high 0.495 vs low 0.527 (n=6 vs 6), Δ=-0.031, P=0.818
- post_all12_allcell_mean | TACSTD2 | all_cell_mean_log1p_cp10k_NOT_independent_of_tumor_fraction | MPR rate high 2/6 vs low 2/6, Fisher P=1.000
- post_all12_allcell_mean | CLDN4 | all_cell_mean_log1p_cp10k_NOT_independent_of_tumor_fraction | A3-malignant fraction: high 0.043 vs low 0.001 (n=6 vs 6), Δ=0.042, P=0.240
- post_all12_allcell_mean | CLDN4 | all_cell_mean_log1p_cp10k_NOT_independent_of_tumor_fraction | CD8 T + NK fraction: high 0.171 vs low 0.253 (n=6 vs 6), Δ=-0.082, P=0.394
- post_all12_allcell_mean | CLDN4 | all_cell_mean_log1p_cp10k_NOT_independent_of_tumor_fraction | all T + NK fraction: high 0.334 vs low 0.573 (n=6 vs 6), Δ=-0.239, P=0.310
- post_all12_allcell_mean | CLDN4 | all_cell_mean_log1p_cp10k_NOT_independent_of_tumor_fraction | MPR rate high 2/6 vs low 2/6, Fisher P=1.000
- post_all12_allcell_mean | dual-high | all_cell_mean_log1p_cp10k_NOT_independent_of_tumor_fraction | A3-malignant fraction: high 0.033 vs low 0.002 (n=5 vs 7), Δ=0.031, P=0.530
- post_all12_allcell_mean | dual-high | all_cell_mean_log1p_cp10k_NOT_independent_of_tumor_fraction | CD8 T + NK fraction: high 0.222 vs low 0.173 (n=5 vs 7), Δ=0.049, P=0.639
- post_all12_allcell_mean | dual-high | all_cell_mean_log1p_cp10k_NOT_independent_of_tumor_fraction | all T + NK fraction: high 0.448 vs low 0.540 (n=5 vs 7), Δ=-0.091, P=0.530
- post_all12_allcell_mean | dual-high | all_cell_mean_log1p_cp10k_NOT_independent_of_tumor_fraction | MPR rate high 2/5 vs low 2/7, Fisher P=1.000
- post_all12_allcell_mean | TACSTD2 | all_cell_mean_log1p_cp10k_NOT_independent_of_tumor_fraction | NMPR vs MPR expression: medians 0.097 vs 0.096 (n=8 vs 4), P=0.933
- post_all12_allcell_mean | CLDN4 | all_cell_mean_log1p_cp10k_NOT_independent_of_tumor_fraction | NMPR vs MPR expression: medians 0.075 vs 0.118 (n=8 vs 4), P=0.808

All-cell means versus epithelial fraction are in `tables/spearman.tsv` under cohort `post_all12_allcell_mean`. Full test rows: `tables/tests.tsv`.

## Caveats

- Barcode labels are reconstructed. This is not a cell-for-cell replay of the slide, and the slide's stated P values are not copied forward as results.
- n=12 post-treatment samples. A median split is one cut; the Spearman column does not depend on that cut. With 4 MPR samples, a Fisher P near 0.06 is a small-sample count.
- P07 (squamous, NMPR) is an epithelial-fraction outlier (about half the cells). Medians, not means, are the high-vs-low summaries.
- CD8/NK is marker-defined: CD3+ and CD8A UMI > 0, or CD3-negative with NCR1, KLRF1, or GNLY together with NKG7. The all-T/NK column is the lineage-score compartment used in the earlier CLDN4-only run.
- Adenocarcinoma versus squamous is unbalanced inside some splits. Counts are listed; there is no multivariable model at this n.
- Dual-high is the intersection of the two median groups. A dual-high count is not a CLDN4 result.

