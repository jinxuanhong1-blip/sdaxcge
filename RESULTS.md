# CosMx He2022: TACSTD2-local immune neighbors and CLDN4

Question: when a tumor cell is TACSTD2-high, are immune and CD8 neighbors fewer at 10, 20, 50, and 100 µm, and does that contrast shrink once CLDN4 is stratified, residualized, or entered in a partial field model?

Data: He et al. 2022 CosMx NSCLC, figshare 25976224, 295877 author-matched tumor cells, 8 sections, 5 patients. Scale 0.18 µm/px, checked by a median nearest-neighbor distance between 3 and 40 µm in every section. The locked CLDN4 cytotoxic ratios 0.36 at 50 µm and 0.52 at 100 µm, 8/8 and 5/5, sign P = 0.031, are not recomputed and are not replaced. Nearby effectors are not scored for GZMB, PRF1, NKG7, or IFNG in this file.

## Answer

At 10 and 20 µm, tumor cells above the section median of TACSTD2 have a lower immune-neighbor fraction in 8/8 sections and 5/5 patients. Equal-weight ratios are 0.519 (0.039 vs 0.075) at 10 µm and 0.643 (0.061 vs 0.096) at 20 µm. CD8 neighbor counts are lower in the same 8/8 and 5/5, with ratios 0.502 (0.0036 vs 0.0072) and 0.726 (0.032 vs 0.044).

That short-range immune-fraction contrast among CLDN4-low tumor cells is 0.507 at 10 µm (8/8 and 5/5) and 0.587 at 20 µm (8/8 and 5/5). Among CLDN4-high cells it is weaker (7/8 at 10 µm, ratio 0.718; 4/8 at 20 µm, ratio 0.887). The median section-level rank attenuation of TACSTD2, after CLDN4 is in the same model, is 0.124 at 10 µm and 0.103 at 20 µm. Lung5 sections sit near zero or below (the TACSTD2 coefficient does not shrink), while LUAD-9 R1 is 0.430 at 10 µm. The 10–20 µm immune-fraction pattern is not accounted for by the cell's own CLDN4.

At 50 µm the immune-fraction ratio is 0.848, lower in 7/8 sections and 5/5 patients. The section that is higher is LUAD-5 R3. At 100 µm the ratio is 0.957, 4/8 and 3/5. Higher sections: LUAD-5 R1, LUAD-5 R2, LUAD-5 R3, LUAD-13.
Counted as lower at 50 µm with a section ratio above 0.99: LUAD-13 ratio 1.000.

At 10 µm the four median quadrants, as equal-weight section means of immune fraction, are both low 0.091, TACSTD2-high / CLDN4-low 0.047, TACSTD2-low / CLDN4-high 0.047, and both high 0.033. TACSTD2-high / CLDN4-low is below both-low in 8/8 sections. Each marker is associated with a colder 10 µm neighborhood when the other marker is low.

The Gaussian field does not reproduce that cell-level 8/8. At 10 µm, LUSC-6, LUAD-12, and LUAD-13 keep no FOV with at least 40 evaluation points after the kernel-mass and local-tumor-count filters, so those sections are absent from the 10 µm field mean. At 20 µm, where all 8 sections pass, the partial Spearman of the TACSTD2 field given the CLDN4 field is negative in 5/8 sections (mean -0.022). Section coefficients change sign. The smoothed field is not an 8/8 TACSTD2 exclusion result.

## How high and low were fixed

Primary split: within each section, tumor cells above the median log-normalized TACSTD2 versus cells at or below that median. Companion splits, fixed before the summaries: raw count detected versus absent, and Q4 versus Q1. Where the log-normalized median is 0, the median split is the detected-versus-absent split. That happens in the sections marked below. It is reported as the same contrast twice, not as extra confirmation.

| Section | Patient | Tumor cells | TACSTD2 > 0 | CLDN4 > 0 | Median is detected | Spearman | NN µm |
|---|---|---|---|---|---|---|---|
| LUAD-5 R1 | Lung5 | 17837 | 0.663 | 0.707 | no | 0.210 | 8.31 |
| LUAD-5 R2 | Lung5 | 17914 | 0.634 | 0.691 | no | 0.175 | 8.16 |
| LUAD-5 R3 | Lung5 | 15893 | 0.591 | 0.628 | no | 0.206 | 8.34 |
| LUSC-6 | Lung6 | 66196 | 0.372 | 0.247 | yes | 0.171 | 9.80 |
| LUAD-9 R1 | Lung9 | 38864 | 0.268 | 0.716 | yes | 0.167 | 8.13 |
| LUAD-9 R2 | Lung9 | 94876 | 0.268 | 0.611 | yes | 0.184 | 8.69 |
| LUAD-12 | Lung12 | 18267 | 0.324 | 0.397 | yes | 0.257 | 8.89 |
| LUAD-13 | Lung13 | 26030 | 0.348 | 0.737 | yes | 0.093 | 7.76 |

Patient is the equal-weight mean of its usable sections. A section that cannot fill both arms (fewer than 30 cells with a finite score) counts against 8/8 and against 5/5. The sign test is one-sided in the exclusion direction (high mean lower than low mean). Its floor is 0.0039 on 8 sections and 0.031 on 5 patients. Wilcoxon p-values are two-sided on the defined units.

## Marginal TACSTD2 high versus low

Equal-weight mean of usable section means. Ratio is high / low. Tallies are high < low.

### Primary: within-section median

| Readout | µm | High | Low | Ratio | Sections | Patients | Sign P (5) | Wilcoxon P (8) |
|---|---|---|---|---|---|---|---|---|
| CD8 count | 10 | 0.0036 | 0.0072 | 0.502 | 8/8 | 5/5 | 0.031 | 0.008 |
| CD8 count | 20 | 0.032 | 0.044 | 0.726 | 8/8 | 5/5 | 0.031 | 0.008 |
| CD8 count | 50 | 0.334 | 0.371 | 0.902 | 6/8 | 5/5 | 0.031 | 0.148 |
| CD8 count | 100 | 1.768 | 1.824 | 0.970 | 5/8 | 3/5 | 0.500 | 0.742 |
| Immune fraction | 10 | 0.039 | 0.075 | 0.519 | 8/8 | 5/5 | 0.031 | 0.008 |
| Immune fraction | 20 | 0.061 | 0.096 | 0.643 | 8/8 | 5/5 | 0.031 | 0.008 |
| Immune fraction | 50 | 0.125 | 0.147 | 0.848 | 7/8 | 5/5 | 0.031 | 0.023 |
| Immune fraction | 100 | 0.185 | 0.193 | 0.957 | 4/8 | 3/5 | 0.500 | 1.000 |

### Companion: detected versus absent

| Readout | µm | High | Low | Ratio | Sections | Patients | Sign P (5) | Wilcoxon P (8) |
|---|---|---|---|---|---|---|---|---|
| CD8 count | 10 | 0.0036 | 0.0072 | 0.505 | 8/8 | 5/5 | 0.031 | 0.008 |
| CD8 count | 20 | 0.032 | 0.045 | 0.726 | 8/8 | 5/5 | 0.031 | 0.008 |
| CD8 count | 50 | 0.335 | 0.369 | 0.907 | 6/8 | 4/5 | 0.188 | 0.312 |
| CD8 count | 100 | 1.771 | 1.818 | 0.974 | 4/8 | 3/5 | 0.500 | 0.945 |
| Immune fraction | 10 | 0.038 | 0.078 | 0.490 | 8/8 | 5/5 | 0.031 | 0.008 |
| Immune fraction | 20 | 0.060 | 0.099 | 0.606 | 8/8 | 5/5 | 0.031 | 0.008 |
| Immune fraction | 50 | 0.123 | 0.150 | 0.823 | 8/8 | 5/5 | 0.031 | 0.008 |
| Immune fraction | 100 | 0.183 | 0.195 | 0.938 | 5/8 | 3/5 | 0.500 | 0.742 |

### Companion: Q4 versus Q1

| Readout | µm | High | Low | Ratio | Sections | Patients | Sign P (5) | Wilcoxon P (8) |
|---|---|---|---|---|---|---|---|---|
| CD8 count | 10 | 0.0042 | 0.0072 | 0.591 | 7/8 | 4/5 | 0.188 | 0.016 |
| CD8 count | 20 | 0.033 | 0.045 | 0.742 | 8/8 | 5/5 | 0.031 | 0.008 |
| CD8 count | 50 | 0.333 | 0.369 | 0.900 | 6/8 | 4/5 | 0.188 | 0.250 |
| CD8 count | 100 | 1.765 | 1.818 | 0.971 | 4/8 | 2/5 | 0.812 | 0.844 |
| Immune fraction | 10 | 0.043 | 0.078 | 0.553 | 8/8 | 5/5 | 0.031 | 0.008 |
| Immune fraction | 20 | 0.068 | 0.099 | 0.683 | 8/8 | 5/5 | 0.031 | 0.008 |
| Immune fraction | 50 | 0.132 | 0.150 | 0.884 | 5/8 | 4/5 | 0.188 | 0.547 |
| Immune fraction | 100 | 0.193 | 0.195 | 0.992 | 4/8 | 3/5 | 0.500 | 0.945 |

## Same-pipeline CLDN4 calibration

Within-section median of CLDN4, CD8+NK neighbor count, same trees and same section weights. This is a calibration of the code path. It is not the locked 0.36 / 0.52 summary.

| µm | CLDN4 high | CLDN4 low | Ratio | Sections | Patients | Δ |
|---|---|---|---|---|---|---|
| 10 | 0.0050 | 0.0102 | 0.494 | 8/8 | 5/5 | -0.0051 |
| 20 | 0.040 | 0.063 | 0.642 | 5/8 | 3/5 | -0.023 |
| 50 | 0.415 | 0.512 | 0.810 | 4/8 | 3/5 | -0.097 |
| 100 | 2.278 | 2.499 | 0.912 | 4/8 | 3/5 | -0.221 |

## Stratified by CLDN4

Inside each section, tumor cells are split at the CLDN4 median. Inside each stratum, TACSTD2 is split at that stratum's own median. The readout is still high versus low TACSTD2.

| Readout | µm | Stratum | TACSTD2 high | TACSTD2 low | Ratio | Sections | Patients | Usable |
|---|---|---|---|---|---|---|---|---|
| CD8 count | 10 | CLDN4 high | 0.0028 | 0.0052 | 0.538 | 6/8 | 4/5 | 8 |
| CD8 count | 10 | CLDN4 low | 0.0048 | 0.0085 | 0.567 | 8/8 | 5/5 | 8 |
| CD8 count | 20 | CLDN4 high | 0.028 | 0.034 | 0.819 | 4/8 | 4/5 | 8 |
| CD8 count | 20 | CLDN4 low | 0.038 | 0.052 | 0.746 | 8/8 | 5/5 | 8 |
| CD8 count | 50 | CLDN4 high | 0.318 | 0.322 | 0.987 | 4/8 | 2/5 | 8 |
| CD8 count | 50 | CLDN4 low | 0.362 | 0.404 | 0.897 | 7/8 | 5/5 | 8 |
| CD8 count | 100 | CLDN4 high | 1.712 | 1.699 | 1.008 | 3/8 | 1/5 | 8 |
| CD8 count | 100 | CLDN4 low | 1.873 | 1.908 | 0.982 | 5/8 | 3/5 | 8 |
| Immune fraction | 10 | CLDN4 high | 0.033 | 0.046 | 0.718 | 7/8 | 4/5 | 8 |
| Immune fraction | 10 | CLDN4 low | 0.047 | 0.092 | 0.507 | 8/8 | 5/5 | 8 |
| Immune fraction | 20 | CLDN4 high | 0.061 | 0.069 | 0.887 | 4/8 | 3/5 | 8 |
| Immune fraction | 20 | CLDN4 low | 0.066 | 0.113 | 0.587 | 8/8 | 5/5 | 8 |
| Immune fraction | 50 | CLDN4 high | 0.128 | 0.128 | 1.005 | 2/8 | 2/5 | 8 |
| Immune fraction | 50 | CLDN4 low | 0.125 | 0.159 | 0.785 | 7/8 | 4/5 | 8 |
| Immune fraction | 100 | CLDN4 high | 0.191 | 0.182 | 1.053 | 2/8 | 2/5 | 8 |
| Immune fraction | 100 | CLDN4 low | 0.181 | 0.200 | 0.903 | 7/8 | 4/5 | 8 |

## Quadrants

Section-wide medians, crossed. `TACSTD2-high / CLDN4-low` is the cell that is high for TROP2 without being high for CLDN4. Means are equal-weight averages of usable sections. A quadrant with fewer than 30 finite cells in a section is unused in that section.

| Readout | µm | Both high | TAC high, CLDN4 low | TAC low, CLDN4 high | Both low | HL<LL | HH<LH | HH<LL |
|---|---|---|---|---|---|---|---|---|
| CD8 count | 10 | 0.0028 | 0.0048 | 0.0052 | 0.0085 | 8/8 | 6/8 | 7/8 |
| CD8 count | 20 | 0.028 | 0.038 | 0.035 | 0.051 | 8/8 | 4/8 | 7/8 |
| CD8 count | 50 | 0.318 | 0.362 | 0.323 | 0.404 | 7/8 | 4/8 | 6/8 |
| CD8 count | 100 | 1.710 | 1.871 | 1.700 | 1.907 | 5/8 | 3/8 | 5/8 |
| Immune fraction | 10 | 0.033 | 0.047 | 0.047 | 0.091 | 8/8 | 7/8 | 8/8 |
| Immune fraction | 20 | 0.059 | 0.066 | 0.070 | 0.111 | 8/8 | 6/8 | 7/8 |
| Immune fraction | 50 | 0.127 | 0.125 | 0.129 | 0.158 | 7/8 | 2/8 | 5/8 |
| Immune fraction | 100 | 0.189 | 0.181 | 0.182 | 0.199 | 7/8 | 2/8 | 5/8 |

## Residual and mediation-like attenuation

Within each section the outcome, TACSTD2, and CLDN4 are rank-transformed and standardized. Attenuation is 1 minus the joint standardized coefficient of TACSTD2 divided by its simple coefficient. The product-of-paths proportion matches that number in this linear specification; both are stored, and the report quotes the attenuation. A simple coefficient smaller than 0.01 in absolute value leaves the attenuation undefined. The linear (log-normalized, not rank) attenuation is the companion. Residual Δ is the median-split TACSTD2 contrast after the outcome is regressed on CLDN4. Negative means TACSTD2-high cells still have the lower neighborhood. These are observational associations. They are not a causal mediation estimate, and overlapping neighborhoods are not treated as independent cells.

| Readout | µm | Rank β total | Rank β with CLDN4 | TACSTD2 attenuation | Defined | CLDN4 attenuation | Residual Δ | Res. sections | Res. patients | Index partial ρ | Partial < 0 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CD8 count | 10 | -0.0096 | -0.0095 | 0.169 | 4/8 | 0.085 | -2.85e-04 | 7/8 | 5/5 | -0.0093 | 6/8 |
| CD8 count | 20 | -0.0151 | -0.0139 | 0.175 | 4/8 | 0.082 | -0.0022 | 7/8 | 5/5 | -0.0136 | 7/8 |
| CD8 count | 50 | -0.0092 | -0.0075 | 0.218 | 5/8 | 0.053 | -0.0092 | 6/8 | 4/5 | -0.0074 | 5/8 |
| CD8 count | 100 | -0.0015 | -0.0055 | 0.043 | 5/8 | -0.0068 | -0.029 | 5/8 | 3/5 | -0.0054 | 5/8 |
| Immune fraction | 10 | -0.051 | -0.049 | 0.124 | 8/8 | 0.113 | -0.0139 | 8/8 | 5/5 | -0.048 | 8/8 |
| Immune fraction | 20 | -0.039 | -0.043 | 0.103 | 8/8 | 0.0128 | -0.0140 | 8/8 | 5/5 | -0.042 | 8/8 |
| Immune fraction | 50 | -0.0049 | -0.0187 | 0.217 | 3/8 | -0.0178 | -0.0072 | 7/8 | 4/5 | -0.0184 | 5/8 |
| Immune fraction | 100 | 0.0174 | 0.0144 | 0.262 | 7/8 | 0.049 | 0.0053 | 3/8 | 3/5 | 0.0141 | 2/8 |

## Partial field models

Each section is rasterized. Gaussian bandwidths are 10, 20, 50, and 100 µm (pixel = bandwidth/5, truncate 4). The TACSTD2 field and the CLDN4 field are kernel-weighted means on tumor cells. The immune field is the kernel-weighted immune share of all cells. The CD8 field is kernel-weighted CD8 density (cells per µm²). Evaluation points are up to 250 tumor cells per FOV, kept where kernel mass is at least 0.80 and the local tumor disk holds at least 3 cells. A FOV needs 40 points. A section needs 3 such FOVs. Reported values are unweighted means of FOV Spearman coefficients, then of sections. Partial ρ residualizes both ranks on the CLDN4 field. β is the standardized coefficient in immune ~ TACSTD2 field + CLDN4 field.

| Field | µm | Marginal ρ | Sections | Patients | Partial ρ | Sections | Patients | β TACSTD2 | β CLDN4 | Field TAC–CLDN4 ρ |
|---|---|---|---|---|---|---|---|---|---|---|
| CD8 density | 10 | -0.095 | 5/8 | 2/5 | -0.063 | 4/8 | 2/5 | -0.038 | -0.061 | 0.286 |
| CD8 density | 20 | -0.0077 | 5/8 | 4/5 | 0.0130 | 3/8 | 3/5 | 0.021 | -0.053 | 0.398 |
| CD8 density | 50 | 0.0018 | 3/8 | 2/5 | 0.031 | 2/8 | 2/5 | 0.0034 | -0.056 | 0.477 |
| CD8 density | 100 | 0.051 | 2/8 | 1/5 | 0.044 | 3/8 | 2/5 | 0.070 | -0.036 | 0.504 |
| Immune fraction | 10 | -0.061 | 4/8 | 2/5 | -0.023 | 3/8 | 2/5 | 0.0192 | -0.101 | 0.286 |
| Immune fraction | 20 | -0.052 | 4/8 | 3/5 | -0.022 | 5/8 | 2/5 | -0.0100 | -0.074 | 0.398 |
| Immune fraction | 50 | -0.028 | 2/8 | 2/5 | 0.0141 | 4/8 | 3/5 | 0.034 | -0.082 | 0.477 |
| Immune fraction | 100 | 0.036 | 2/8 | 1/5 | 0.037 | 4/8 | 3/5 | 0.047 | -0.028 | 0.504 |

## Reading the dependence test

A radius is called concordant only at 8/8 sections and 5/5 patients. It is called largely shared with CLDN4 only when that concordance holds, neither CLDN4 stratum is itself 8/8 and 5/5, and the median rank attenuation is above 0.5. It is called not accounted for by CLDN4 only when concordance holds, at least one stratum is itself 8/8 and 5/5, and the median rank attenuation is below 0.25. Every other concordant radius is mixed. A radius that is not 8/8 and 5/5 is reported as not concordant. Companion cuts do not get this label.

- **Immune fraction, 10 µm.** Median-split ratio 0.519 (0.039 vs 0.075), 8/8 sections, 5/5 patients. Strata: cldn4_high 7/8 and 4/5, ratio 0.718; cldn4_low 8/8 and 5/5, ratio 0.507. Rank attenuation of TACSTD2 0.124 (8/8 defined); rank attenuation of CLDN4 0.113. The median split stays lower inside a CLDN4 stratum at 8/8 and 5/5, and the rank attenuation is below one quarter.
- **Immune fraction, 20 µm.** Median-split ratio 0.643 (0.061 vs 0.096), 8/8 sections, 5/5 patients. Strata: cldn4_high 4/8 and 3/5, ratio 0.887; cldn4_low 8/8 and 5/5, ratio 0.587. Rank attenuation of TACSTD2 0.103 (8/8 defined); rank attenuation of CLDN4 0.0128. The median split stays lower inside a CLDN4 stratum at 8/8 and 5/5, and the rank attenuation is below one quarter.
- **Immune fraction, 50 µm.** Median-split ratio 0.848 (0.125 vs 0.147), 7/8 sections, 5/5 patients. Strata: cldn4_high 2/8 and 2/5, ratio 1.005; cldn4_low 7/8 and 4/5, ratio 0.785. Rank attenuation of TACSTD2 0.217 (3/8 defined); rank attenuation of CLDN4 -0.0178. The pre-specified median split is not lower in all 8 sections and all 5 patients.
- **Immune fraction, 100 µm.** Median-split ratio 0.957 (0.185 vs 0.193), 4/8 sections, 3/5 patients. Strata: cldn4_high 2/8 and 2/5, ratio 1.053; cldn4_low 7/8 and 4/5, ratio 0.903. Rank attenuation of TACSTD2 0.262 (7/8 defined); rank attenuation of CLDN4 0.049. The pre-specified median split is not lower in all 8 sections and all 5 patients.
- **CD8 count, 10 µm.** Median-split ratio 0.502 (0.0036 vs 0.0072), 8/8 sections, 5/5 patients. Strata: cldn4_high 6/8 and 4/5, ratio 0.538; cldn4_low 8/8 and 5/5, ratio 0.567. Rank attenuation of TACSTD2 0.169 (4/8 defined); rank attenuation of CLDN4 0.085. The median split is lower in 8/8 and 5/5, and the CLDN4 adjustment is only partial.
- **CD8 count, 20 µm.** Median-split ratio 0.726 (0.032 vs 0.044), 8/8 sections, 5/5 patients. Strata: cldn4_high 4/8 and 4/5, ratio 0.819; cldn4_low 8/8 and 5/5, ratio 0.746. Rank attenuation of TACSTD2 0.175 (4/8 defined); rank attenuation of CLDN4 0.082. The median split is lower in 8/8 and 5/5, and the CLDN4 adjustment is only partial.
- **CD8 count, 50 µm.** Median-split ratio 0.902 (0.334 vs 0.371), 6/8 sections, 5/5 patients. Strata: cldn4_high 4/8 and 2/5, ratio 0.987; cldn4_low 7/8 and 5/5, ratio 0.897. Rank attenuation of TACSTD2 0.218 (5/8 defined); rank attenuation of CLDN4 0.053. The pre-specified median split is not lower in all 8 sections and all 5 patients.
- **CD8 count, 100 µm.** Median-split ratio 0.970 (1.768 vs 1.824), 5/8 sections, 3/5 patients. Strata: cldn4_high 3/8 and 1/5, ratio 1.008; cldn4_low 5/8 and 3/5, ratio 0.982. Rank attenuation of TACSTD2 0.043 (5/8 defined); rank attenuation of CLDN4 -0.0068. The pre-specified median split is not lower in all 8 sections and all 5 patients.

## Section ratios for the median split

Ratio is the section high mean divided by the section low mean. A ratio below 1 is the exclusion direction.

### Immune fraction

| Section | Patient | 10 | 20 | 50 | 100 |
|---|---|---|---|---|---|
| 5-R1 | Lung5 | 0.571 | 0.758 | 0.985 | 1.077 |
| 5-R2 | Lung5 | 0.622 | 0.773 | 0.974 | 1.103 |
| 5-R3 | Lung5 | 0.531 | 0.842 | 1.008 | 1.096 |
| 6 | Lung6 | 0.661 | 0.669 | 0.808 | 0.877 |
| 9-R1 | Lung9 | 0.550 | 0.637 | 0.842 | 0.960 |
| 9-R2 | Lung9 | 0.341 | 0.388 | 0.562 | 0.710 |
| 12 | Lung12 | 0.400 | 0.412 | 0.620 | 0.774 |
| 13 | Lung13 | 0.780 | 0.914 | 1.000 | 1.019 |

### CD8 count

| Section | Patient | 10 | 20 | 50 | 100 |
|---|---|---|---|---|---|
| 5-R1 | Lung5 | 0.0000 | 0.386 | 0.851 | 0.931 |
| 5-R2 | Lung5 | 0.0000 | 0.781 | 1.026 | 1.132 |
| 5-R3 | Lung5 | 0.334 | 0.930 | 0.909 | 0.984 |
| 6 | Lung6 | 0.842 | 0.831 | 0.982 | 0.954 |
| 9-R1 | Lung9 | 0.851 | 0.946 | 1.120 | 1.130 |
| 9-R2 | Lung9 | 0.415 | 0.505 | 0.624 | 0.755 |
| 12 | Lung12 | 0.113 | 0.280 | 0.499 | 0.714 |
| 13 | Lung13 | 0.555 | 0.829 | 0.976 | 1.005 |

CD8 counts at 10 µm are sparse. High-arm and low-arm event totals (mean count times cells in the arm) are:

| Section | High-arm CD8 events | Low-arm CD8 events |
|---|---|---|
| LUAD-5 R1 | 0.0 | 3.0 |
| LUAD-5 R2 | 0.0 | 2.0 |
| LUAD-5 R3 | 1.0 | 3.0 |
| LUSC-6 | 15.0 | 30.0 |
| LUAD-9 R1 | 64.0 | 205.0 |
| LUAD-9 R2 | 74.0 | 487.0 |
| LUAD-12 | 6.0 | 111.0 |
| LUAD-13 | 163.0 | 551.0 |

Sections with fewer than 2 CD8-neighbor events in the high arm at 10 µm: LUAD-5 R1 0 vs 3; LUAD-5 R2 0 vs 2; LUAD-5 R3 1 vs 3. They still count toward 8/8. The immune-fraction contrast is not carried by these CD8 events.

## Empty-neighborhood check

Immune fraction is undefined when a tumor cell has no other cell inside the radius, and those cells are left out of the fraction mean. After that was visible, the share of tumor cells with at least one immune neighbor was computed on every tumor cell, including isolated ones. It is a check on the fraction, not a second primary endpoint.

| Section | 10 µm high/low | 20 µm high/low |
|---|---|---|
| 5-R1 | 0.021/0.038 | 0.165/0.201 |
| 5-R2 | 0.023/0.034 | 0.170/0.206 |
| 5-R3 | 0.0193/0.035 | 0.169/0.193 |
| 6 | 0.0062/0.0112 | 0.068/0.087 |
| 9-R1 | 0.027/0.054 | 0.154/0.194 |
| 9-R2 | 0.0169/0.056 | 0.099/0.202 |
| 12 | 0.057/0.183 | 0.261/0.508 |
| 13 | 0.084/0.134 | 0.517/0.550 |

Any-immune contact is lower for TACSTD2-high cells in 8/8 sections at 10 µm and 8/8 at 20 µm.

## Field coefficients by section

Mean partial ρ in the summary table averages usable sections only. A section with fewer than 3 usable FOVs is omitted from that mean and counts against 8/8.

| Section | µm | FOVs used | Marginal ρ | Partial ρ | β TACSTD2 | β CLDN4 |
|---|---|---|---|---|---|---|
| 5-R1 | 10 | 10 | -0.033 | 0.0122 | 0.0090 | -0.125 |
| 5-R2 | 10 | 11 | -0.046 | -0.0153 | 0.068 | -0.034 |
| 5-R3 | 10 | 8 | -0.166 | -0.085 | -0.095 | -0.080 |
| 6 | 10 | 0 | unused | unused | unused | unused |
| 9-R1 | 10 | 5 | 0.060 | 0.078 | 0.182 | -0.179 |
| 9-R2 | 10 | 7 | -0.121 | -0.107 | -0.068 | -0.088 |
| 12 | 10 | 0 | unused | unused | unused | unused |
| 13 | 10 | 0 | unused | unused | unused | unused |
| 5-R1 | 20 | 22 | -0.049 | -0.110 | -0.111 | 0.199 |
| 5-R2 | 20 | 21 | 0.0040 | -0.0190 | -0.0031 | 0.120 |
| 5-R3 | 20 | 21 | -0.052 | -0.054 | -0.0122 | 0.040 |
| 6 | 20 | 29 | 0.023 | 0.0061 | -0.024 | 0.078 |
| 9-R1 | 20 | 19 | 0.039 | 0.159 | 0.199 | -0.410 |
| 9-R2 | 20 | 44 | -0.165 | -0.023 | 0.028 | -0.351 |
| 12 | 20 | 23 | -0.241 | -0.177 | -0.208 | -0.154 |
| 13 | 20 | 20 | 0.022 | 0.040 | 0.051 | -0.117 |
| 5-R1 | 50 | 22 | 0.053 | -0.055 | -0.097 | 0.349 |
| 5-R2 | 50 | 22 | 0.090 | 0.024 | 0.136 | 0.133 |
| 5-R3 | 50 | 21 | 0.083 | -0.0192 | -0.0036 | 0.343 |
| 6 | 50 | 29 | 0.0112 | -0.029 | -0.055 | 0.136 |
| 9-R1 | 50 | 19 | 0.091 | 0.251 | 0.306 | -0.560 |
| 9-R2 | 50 | 45 | -0.180 | 0.076 | 0.158 | -0.605 |
| 12 | 50 | 28 | -0.471 | -0.252 | -0.280 | -0.356 |
| 13 | 50 | 20 | 0.097 | 0.117 | 0.110 | -0.101 |
| 5-R1 | 100 | 22 | 0.132 | -0.134 | -0.222 | 0.579 |
| 5-R2 | 100 | 22 | 0.183 | 0.056 | 0.106 | 0.328 |
| 5-R3 | 100 | 22 | 0.224 | -0.0014 | -0.0123 | 0.478 |
| 6 | 100 | 29 | 0.038 | -0.021 | -0.052 | 0.132 |
| 9-R1 | 100 | 19 | 0.162 | 0.310 | 0.367 | -0.581 |
| 9-R2 | 100 | 45 | -0.123 | 0.118 | 0.223 | -0.579 |
| 12 | 100 | 28 | -0.455 | -0.207 | -0.235 | -0.380 |
| 13 | 100 | 20 | 0.130 | 0.174 | 0.201 | -0.201 |

## What this file does not say

It does not replace the locked CLDN4 exclusion summary. It does not say that effector cells next to TACSTD2-high tumor cells are muzzled. It does not use an ICI label. It does not pool private KL tumors with this public object. FOV coefficients are descriptive because cells inside a FOV share neighbors.

```bash
python3 scripts/download_cosmx_nsclc_h5ad.py
python3 scripts/cosmx_tacstd2_cldn4_exclusion.py
```
