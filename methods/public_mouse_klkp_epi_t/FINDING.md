# Public mouse K/KP/KL lung scRNA: epithelium vs T, with NHEJ/STING/IFN

Processed GEO matrices only. Private 8-KL matrices were not used and were not merged.
Mice are not pooled across accessions. Scores are means on the log-normalized matrix described in `tables/qc_scale.tsv`.
Epithelial means require at least 20 epithelial cells in that mouse. T fraction is T cells / QC cells in that mouse.

## Cell compartments

Marker rule: T if Cd3d or Cd3e is detected and the T-marker mean is at least the epithelial-marker mean; epithelium if at least two of Epcam/Krt8/Krt18/Sftpc/Sftpb/Nkx2-1 are detected, Cd3 is absent, and the epithelial mean exceeds the T mean. A doublet requires Cd3, Epcam, and a second epithelial marker. Ambient Sftpc alone does not make a doublet. GSE127465 uses deposited Major cell type for T and NK, and the marker rule for epithelium.

| dataset | n_mice | n_cells | n_epi | n_T |
| --- | --- | --- | --- | --- |
| GSE127465 | 4 | 15915 | 15 | 1979 |
| GSE136246 | 7 | 15158 | 4371 | 510 |
| GSE149813 | 2 | 9384 | 9372 | 1 |
| GSE154977 | 4 | 11017 | 11009 | 0 |
| GSE154989 | 30 | 3891 | 3847 | 0 |
| GSE165641 | 2 | 6594 | 1043 | 243 |
| GSE179501 | 4 | 20122 | 1564 | 4401 |
| GSE179502 | 6 | 12643 | 12478 | 1 |
| GSE180963 | 2 | 14260 | 2177 | 5424 |
| GSE267321 | 6 | 7956 | 319 | 1634 |

## Threshold sweep

Inside each accession, epithelial Cldn4 is cut at 0 and at the pooled epithelial quantiles 0.50, 0.60, 0.70, 0.75, 0.80, 0.90, and 0.95. The mouse feature is the fraction of that mouse's epithelial cells above the cut. The outcome is the T-cell fraction. Spearman is reported only for n>=4 mice with an epithelial score. The pre-specified primary cut is the 75th percentile. The strongest grid point is the minimum Spearman among eligible cuts. Eligibility requires a mixed lung design (not CD45-sorted, not epithelium-sorted, not a tumor-cell plate, not the subcutaneous non-malignant series), at least four mice, and a non-zero T count. The minimum over the grid is a grid result, not a locked single test. When Cldn4 is zero in most epithelial cells, several nominal quantiles collapse to the same cut (Cldn4 > 0).

Strongest grid point: GSE136246 threshold gt0 (value 0), n=7, Spearman rho=-0.3571, two-sided p=0.4316, mean T fraction 0.0340.
Same rho in that accession for: gt0=0, q0.50=0, q0.60=0, q0.70=0, q0.75=0, q0.80=0, q0.90=0, q0.95=0.796746.
Primary q=0.75 in that accession: value 0, rho=-0.3571, p=0.4316.
Other eligible accession GSE179501 gt0: n=4, rho=0.4000, p=0.6.
Read the grid minimum together with the other eligible rows. A single negative rho on a zero-inflated Cldn4 cut is not a stable Cldn4-high / T-low result.

| dataset | threshold_name | threshold_value | n_corr | rho_frac_T | p_frac_T | eligible_for_strongest | reason |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE154989 | gt0 | 0 | 28 |  |  | False | T compartment too thin |
| GSE154989 | q0.50 | 0.0198552 | 28 |  |  | False | T compartment too thin |
| GSE154989 | q0.60 | 0.445674 | 28 |  |  | False | T compartment too thin |
| GSE154989 | q0.70 | 0.972537 | 28 |  |  | False | T compartment too thin |
| GSE154989 | q0.75 | 1.26604 | 28 |  |  | False | T compartment too thin |
| GSE154989 | q0.80 | 1.58637 | 28 |  |  | False | T compartment too thin |
| GSE154989 | q0.90 | 2.35332 | 28 |  |  | False | T compartment too thin |
| GSE154989 | q0.95 | 2.77168 | 28 |  |  | False | T compartment too thin |
| GSE154977 | gt0 | 0 | 4 |  |  | False | T compartment too thin |
| GSE154977 | q0.50 | 0.197909 | 4 |  |  | False | T compartment too thin |
| GSE154977 | q0.60 | 0.440212 | 4 |  |  | False | T compartment too thin |
| GSE154977 | q0.70 | 0.710056 | 4 |  |  | False | T compartment too thin |
| GSE154977 | q0.75 | 0.883222 | 4 |  |  | False | T compartment too thin |
| GSE154977 | q0.80 | 1.07367 | 4 |  |  | False | T compartment too thin |
| GSE154977 | q0.90 | 1.56064 | 4 |  |  | False | T compartment too thin |
| GSE154977 | q0.95 | 1.9255 | 4 |  |  | False | T compartment too thin |
| GSE179501 | gt0 | 0 | 4 | 0.4 | 0.6 | True |  |
| GSE179501 | q0.50 | 0 | 4 | 0.4 | 0.6 | True |  |
| GSE179501 | q0.60 | 0 | 4 | 0.4 | 0.6 | True |  |
| GSE179501 | q0.70 | 0 | 4 | 0.4 | 0.6 | True |  |
| GSE179501 | q0.75 | 0 | 4 | 0.4 | 0.6 | True |  |
| GSE179501 | q0.80 | 0 | 4 | 0.4 | 0.6 | True |  |
| GSE179501 | q0.90 | 0 | 4 | 0.4 | 0.6 | True |  |
| GSE179501 | q0.95 | 0.865899 | 4 | 0.4 | 0.6 | True |  |
| GSE179502 | gt0 | 0 | 6 |  |  | False | T compartment too thin |
| GSE179502 | q0.50 | 0 | 6 |  |  | False | T compartment too thin |
| GSE179502 | q0.60 | 0 | 6 |  |  | False | T compartment too thin |
| GSE179502 | q0.70 | 0 | 6 |  |  | False | T compartment too thin |
| GSE179502 | q0.75 | 0 | 6 |  |  | False | T compartment too thin |
| GSE179502 | q0.80 | 0.155227 | 6 |  |  | False | T compartment too thin |
| GSE179502 | q0.90 | 0.843092 | 6 |  |  | False | T compartment too thin |
| GSE179502 | q0.95 | 1.5111 | 6 |  |  | False | T compartment too thin |
| GSE165641 | gt0 | 0 | 2 |  |  | False | n<4 |
| GSE165641 | q0.50 | 0 | 2 |  |  | False | n<4 |
| GSE165641 | q0.60 | 0 | 2 |  |  | False | n<4 |
| GSE165641 | q0.70 | 0 | 2 |  |  | False | n<4 |
| GSE165641 | q0.75 | 0 | 2 |  |  | False | n<4 |
| GSE165641 | q0.80 | 0.539805 | 2 |  |  | False | n<4 |
| GSE165641 | q0.90 | 1.5873 | 2 |  |  | False | n<4 |
| GSE165641 | q0.95 | 2.03884 | 2 |  |  | False | n<4 |
| GSE180963 | gt0 | 0 | 2 |  |  | False | n<4 |
| GSE180963 | q0.50 | 0 | 2 |  |  | False | n<4 |
| GSE180963 | q0.60 | 0 | 2 |  |  | False | n<4 |
| GSE180963 | q0.70 | 0 | 2 |  |  | False | n<4 |
| GSE180963 | q0.75 | 0 | 2 |  |  | False | n<4 |
| GSE180963 | q0.80 | 0 | 2 |  |  | False | n<4 |
| GSE180963 | q0.90 | 0 | 2 |  |  | False | n<4 |
| GSE180963 | q0.95 | 0 | 2 |  |  | False | n<4 |
| GSE149813 | gt0 | 0 | 2 |  |  | False | n<4 |
| GSE149813 | q0.50 | 0 | 2 |  |  | False | n<4 |
| GSE149813 | q0.60 | 0 | 2 |  |  | False | n<4 |
| GSE149813 | q0.70 | 0 | 2 |  |  | False | n<4 |
| GSE149813 | q0.75 | 0 | 2 |  |  | False | n<4 |
| GSE149813 | q0.80 | 0 | 2 |  |  | False | n<4 |
| GSE149813 | q0.90 | 0 | 2 |  |  | False | n<4 |
| GSE149813 | q0.95 | 0.768902 | 2 |  |  | False | n<4 |
| GSE267321 | gt0 | 0 | 6 | 0.676123 | 0.140357 | False |  |
| GSE267321 | q0.50 | 0 | 6 | 0.676123 | 0.140357 | False |  |
| GSE267321 | q0.60 | 0 | 6 | 0.676123 | 0.140357 | False |  |
| GSE267321 | q0.70 | 0 | 6 | 0.676123 | 0.140357 | False |  |
| GSE267321 | q0.75 | 0 | 6 | 0.676123 | 0.140357 | False |  |
| GSE267321 | q0.80 | 0 | 6 | 0.676123 | 0.140357 | False |  |
| GSE267321 | q0.90 | 0 | 6 | 0.676123 | 0.140357 | False |  |
| GSE267321 | q0.95 | 0 | 6 | 0.676123 | 0.140357 | False |  |
| GSE136246 | gt0 | 0 | 7 | -0.357143 | 0.431611 | True |  |
| GSE136246 | q0.50 | 0 | 7 | -0.357143 | 0.431611 | True |  |
| GSE136246 | q0.60 | 0 | 7 | -0.357143 | 0.431611 | True |  |
| GSE136246 | q0.70 | 0 | 7 | -0.357143 | 0.431611 | True |  |
| GSE136246 | q0.75 | 0 | 7 | -0.357143 | 0.431611 | True |  |
| GSE136246 | q0.80 | 0 | 7 | -0.357143 | 0.431611 | True |  |
| GSE136246 | q0.90 | 0 | 7 | -0.357143 | 0.431611 | True |  |
| GSE136246 | q0.95 | 0.796746 | 7 | -0.357143 | 0.431611 | True |  |

## Genotype contrasts on epithelial means

Mann-Whitney p is filled only when both sides have at least 3 mice. GSE154989 contrasts K vs KP on the plate tumor-cell series (no T-cell fraction). The all-stage contrast mixes stages that are not balanced across genotypes; within-stage rows are the ones to read, and none of those rows has n>=3 on both sides. GSE180963 is one K library and one KL library.

| dataset | subset | contrast | score | n_left | n_right | median_left | median_right | U | p | delta_median |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE154989 | all_stages | KP_minus_K | epi_mean_Cldn4 | 9 | 15 | 0.423231 | 0.86235 | 103 | 0.0368884 | 0.439119 |
| GSE154989 | all_stages | KP_minus_K | epi_mean_Tacstd2 | 9 | 15 | 0.520446 | 0.379711 | 56 | 0.51188 | -0.140734 |
| GSE154989 | all_stages | KP_minus_K | epi_mean_NHEJ | 9 | 15 | 0.0975027 | 0.1129 | 113 | 0.00729036 | 0.0153969 |
| GSE154989 | all_stages | KP_minus_K | epi_mean_STING | 9 | 15 | 0.0726137 | 0.0730306 | 61 | 0.720515 | 0.000416882 |
| GSE154989 | all_stages | KP_minus_K | epi_mean_IFN | 9 | 15 | 0.024249 | 0.0203956 | 43 | 0.152406 | -0.00385336 |
| GSE154989 | 01_T_early_ND | KP_minus_K | epi_mean_Cldn4 | 0 | 0 |  |  |  |  |  |
| GSE154989 | 01_T_early_ND | KP_minus_K | epi_mean_Tacstd2 | 0 | 0 |  |  |  |  |  |
| GSE154989 | 01_T_early_ND | KP_minus_K | epi_mean_NHEJ | 0 | 0 |  |  |  |  |  |
| GSE154989 | 01_T_early_ND | KP_minus_K | epi_mean_STING | 0 | 0 |  |  |  |  |  |
| GSE154989 | 01_T_early_ND | KP_minus_K | epi_mean_IFN | 0 | 0 |  |  |  |  |  |
| GSE154989 | 02_KorKP_early_ND | KP_minus_K | epi_mean_Cldn4 | 3 | 1 | 0.145135 | 0.250465 |  |  | 0.10533 |
| GSE154989 | 02_KorKP_early_ND | KP_minus_K | epi_mean_Tacstd2 | 3 | 1 | 0.519281 | 0.470574 |  |  | -0.0487071 |
| GSE154989 | 02_KorKP_early_ND | KP_minus_K | epi_mean_NHEJ | 3 | 1 | 0.0975027 | 0.106121 |  |  | 0.00861813 |
| GSE154989 | 02_KorKP_early_ND | KP_minus_K | epi_mean_STING | 3 | 1 | 0.0726137 | 0.118892 |  |  | 0.046278 |
| GSE154989 | 02_KorKP_early_ND | KP_minus_K | epi_mean_IFN | 3 | 1 | 0.0276872 | 0.0203956 |  |  | -0.0072916 |
| GSE154989 | 04_K_12w_ND | KP_minus_K | epi_mean_Cldn4 | 3 | 0 | 0.423231 |  |  |  |  |
| GSE154989 | 04_K_12w_ND | KP_minus_K | epi_mean_Tacstd2 | 3 | 0 | 0.581598 |  |  |  |  |
| GSE154989 | 04_K_12w_ND | KP_minus_K | epi_mean_NHEJ | 3 | 0 | 0.0864197 |  |  |  |  |
| GSE154989 | 04_K_12w_ND | KP_minus_K | epi_mean_STING | 3 | 0 | 0.0750375 |  |  |  |  |
| GSE154989 | 04_K_12w_ND | KP_minus_K | epi_mean_IFN | 3 | 0 | 0.024249 |  |  |  |  |
| GSE154989 | 05_K_30w_ND | KP_minus_K | epi_mean_Cldn4 | 3 | 0 | 0.451802 |  |  |  |  |
| GSE154989 | 05_K_30w_ND | KP_minus_K | epi_mean_Tacstd2 | 3 | 0 | 0.520446 |  |  |  |  |
| GSE154989 | 05_K_30w_ND | KP_minus_K | epi_mean_NHEJ | 3 | 0 | 0.108125 |  |  |  |  |
| GSE154989 | 05_K_30w_ND | KP_minus_K | epi_mean_STING | 3 | 0 | 0.0696538 |  |  |  |  |
| GSE154989 | 05_K_30w_ND | KP_minus_K | epi_mean_IFN | 3 | 0 | 0.0220683 |  |  |  |  |
| GSE154989 | 06_KP_12w_ND | KP_minus_K | epi_mean_Cldn4 | 0 | 6 |  | 0.906663 |  |  |  |
| GSE154989 | 06_KP_12w_ND | KP_minus_K | epi_mean_Tacstd2 | 0 | 6 |  | 0.551143 |  |  |  |
| GSE154989 | 06_KP_12w_ND | KP_minus_K | epi_mean_NHEJ | 0 | 6 |  | 0.117104 |  |  |  |
| GSE154989 | 06_KP_12w_ND | KP_minus_K | epi_mean_STING | 0 | 6 |  | 0.0752139 |  |  |  |
| GSE154989 | 06_KP_12w_ND | KP_minus_K | epi_mean_IFN | 0 | 6 |  | 0.0228855 |  |  |  |
| GSE154989 | 07_KP_20w_ND | KP_minus_K | epi_mean_Cldn4 | 0 | 5 |  | 0.565399 |  |  |  |
| GSE154989 | 07_KP_20w_ND | KP_minus_K | epi_mean_Tacstd2 | 0 | 5 |  | 0.193152 |  |  |  |
| GSE154989 | 07_KP_20w_ND | KP_minus_K | epi_mean_NHEJ | 0 | 5 |  | 0.0998615 |  |  |  |
| GSE154989 | 07_KP_20w_ND | KP_minus_K | epi_mean_STING | 0 | 5 |  | 0.0659448 |  |  |  |
| GSE154989 | 07_KP_20w_ND | KP_minus_K | epi_mean_IFN | 0 | 5 |  | 0.0179345 |  |  |  |
| GSE154989 | 08_KP_30w_ND | KP_minus_K | epi_mean_Cldn4 | 0 | 3 |  | 0.976445 |  |  |  |
| GSE154989 | 08_KP_30w_ND | KP_minus_K | epi_mean_Tacstd2 | 0 | 3 |  | 0.379711 |  |  |  |
| GSE154989 | 08_KP_30w_ND | KP_minus_K | epi_mean_NHEJ | 0 | 3 |  | 0.121652 |  |  |  |
| GSE154989 | 08_KP_30w_ND | KP_minus_K | epi_mean_STING | 0 | 3 |  | 0.0677737 |  |  |  |
| GSE154989 | 08_KP_30w_ND | KP_minus_K | epi_mean_IFN | 0 | 3 |  | 0.0142699 |  |  |  |
| GSE179501 | all_mice | Restored_minus_Non-Restored | epi_mean_Cldn4 | 2 | 2 | 0.0891247 | 0.101714 |  |  | 0.0125896 |
| GSE179501 | all_mice | Restored_minus_Non-Restored | epi_mean_Tacstd2 | 2 | 2 | 0.419254 | 0.315181 |  |  | -0.104073 |
| GSE179501 | all_mice | Restored_minus_Non-Restored | epi_mean_NHEJ | 2 | 2 | 0.046897 | 0.0514238 |  |  | 0.00452681 |
| GSE179501 | all_mice | Restored_minus_Non-Restored | epi_mean_STING | 2 | 2 | 0.104994 | 0.107563 |  |  | 0.00256908 |
| GSE179501 | all_mice | Restored_minus_Non-Restored | epi_mean_IFN | 2 | 2 | 0.0181198 | 0.0319369 |  |  | 0.0138172 |
| GSE180963 | all | KL_minus_K | epi_mean_Cldn4 | 1 | 1 | 0.00068244 | 0.0173158 |  |  | 0.0166333 |
| GSE180963 | all | KL_minus_K | epi_mean_Tacstd2 | 1 | 1 | 0.0214211 | 0.144694 |  |  | 0.123273 |
| GSE180963 | all | KL_minus_K | epi_mean_NHEJ | 1 | 1 | 0.077699 | 0.0691365 |  |  | -0.00856254 |
| GSE180963 | all | KL_minus_K | epi_mean_STING | 1 | 1 | 0.101692 | 0.132605 |  |  | 0.0309134 |
| GSE180963 | all | KL_minus_K | epi_mean_IFN | 1 | 1 | 0.0714614 | 0.166445 |  |  | 0.0949837 |
| GSE267321 | all | KK_minus_K | epi_mean_Cldn4 | 2 | 2 | 0.12037 | 0.00595238 |  |  | -0.114418 |
| GSE267321 | all | KK_minus_K | epi_mean_Tacstd2 | 2 | 2 | 0.0185185 | 0.0238095 |  |  | 0.00529101 |
| GSE267321 | all | KK_minus_K | epi_mean_NHEJ | 2 | 2 | 0.0593915 | 0.0544643 |  |  | -0.00492725 |
| GSE267321 | all | KK_minus_K | epi_mean_STING | 2 | 2 | 0.155234 | 0.276786 |  |  | 0.121551 |
| GSE267321 | all | KK_minus_K | epi_mean_IFN | 2 | 2 | 0.772323 | 1.03274 |  |  | 0.260415 |
| GSE267321 | all | KLK_minus_K | epi_mean_Cldn4 | 2 | 2 | 0.12037 | 0 |  |  | -0.12037 |
| GSE267321 | all | KLK_minus_K | epi_mean_Tacstd2 | 2 | 2 | 0.0185185 | 0.104651 |  |  | 0.0861326 |
| GSE267321 | all | KLK_minus_K | epi_mean_NHEJ | 2 | 2 | 0.0593915 | 0.0501762 |  |  | -0.00921535 |
| GSE267321 | all | KLK_minus_K | epi_mean_STING | 2 | 2 | 0.155234 | 0.244362 |  |  | 0.0891279 |
| GSE267321 | all | KLK_minus_K | epi_mean_IFN | 2 | 2 | 0.772323 | 1.86934 |  |  | 1.09702 |
| GSE267321 | all | KLK_minus_KK | epi_mean_Cldn4 | 2 | 2 | 0.00595238 | 0 |  |  | -0.00595238 |
| GSE267321 | all | KLK_minus_KK | epi_mean_Tacstd2 | 2 | 2 | 0.0238095 | 0.104651 |  |  | 0.0808416 |
| GSE267321 | all | KLK_minus_KK | epi_mean_NHEJ | 2 | 2 | 0.0544643 | 0.0501762 |  |  | -0.00428811 |
| GSE267321 | all | KLK_minus_KK | epi_mean_STING | 2 | 2 | 0.276786 | 0.244362 |  |  | -0.0324235 |
| GSE267321 | all | KLK_minus_KK | epi_mean_IFN | 2 | 2 | 1.03274 | 1.86934 |  |  | 0.836607 |

## Design notes

- GSE165641: two KL lung libraries, mixed cells. GSM5047303 characteristics literally say Lkb2fl/fl; the series title is Lkb1. Genotype group is KL for both. n=2, so the sweep is descriptive.
- GSE180963: one KrasG12D/+ library and one KrasG12D/+;Lkb1fl/fl library. The GEO growth protocol also describes lenti targeting. n=2.
- GSE154977: four KP 30-week 10x libraries (two ND, two Cis72). Sweep eligibility is revoked if T cells are absent after the marker rule.
- GSE179502: six Lkb1-XTR mice, sorted neoplastic cells. Epithelial scores only.
- GSE179501: four Lkb1-XTR mice, total viable cells (the unsorted sister series). This is the KL mixed-lung series with mouse n=4.
- GSE154989: K and KP plate tumor cells. Mouse id strips the tumor suffix (`_T#`). T-cell sweep is not eligible.
- GSE267321: LKR13 subcutaneous K / KK / KLK, deposited as non-malignant cells. Not lung, not eligible for the lung T-cell sweep.
- GSE127465: mouse arm is CD45-positive cells from two healthy and two tumor-bearing lungs. Author T/NK labels are used. Epithelial Cldn4 is not an epithelial compartment here.
- GSE136246: Laughney KP lung, RBC-depleted, seven library tokens parsed from cell names. Treatment is the library token (NegCtrl/VHL/PTC), not relabeled.
- GSE149813: two K mice, YFP-sorted epithelium at 7 weeks. Extra YFP-positive epithelial means are in the mouse table.
- Found and not scored: GSE277777 (KP tumor-cell objects; combined h5ad is 8.0 GB) and GSE319598 (KcP eGFP+ tumor cells; stroma is pooled across mice, so per-mouse T is not in the deposit).

QC scale modes are in `tables/qc_scale.tsv`. Label means are in `tables/label_means.tsv`.
