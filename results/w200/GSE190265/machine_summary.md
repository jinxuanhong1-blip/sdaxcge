# Machine summary (numbers only)

GSE190265 TPM genes in file: 33976 (last=ZZZ3); missing=[]
GSE190266 TPM genes in file: 16383 (last=MTMR14); missing=['TACSTD2']

GSE190265 match: {'n_clin': 43, 'n_expr': 43, 'n_overlap': 43, 'n_DCB': 14, 'n_NDB': 29, 'n_excluded': 0}
GSE190266 match: {'n_clin': 70, 'n_expr': 70, 'n_overlap': 70, 'n_DCB': 17, 'n_NDB': 52, 'n_excluded': 1}

## per-cohort
| cohort | gene | n_DCB | n_NDB | median_TPM_DCB | median_TPM_NDB | median_log2_DCB | median_log2_NDB | delta_median_log2 | AUC_DCB_gt_NDB | rank_biserial | wilcoxon_p | spearman_rho_vs_PFS | spearman_p | logrank_high_vs_low_Z | logrank_high_vs_low_p | pfs_time_capped | note | wilcoxon_p_BH |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE190265 | TACSTD2 | 14 | 29 | 72.5354 | 61.1058 | 6.1974 | 5.9567 | 0.2407 | 0.4581 | -0.0837 | 0.6689 | -0.0315 | 0.8411 | -0.1526 | 0.8787 | False |  | 0.9277 |
| GSE190265 | CLDN4 | 14 | 29 | 22.4481 | 17.1151 | 4.5420 | 4.1791 | 0.3629 | 0.4901 | -0.0197 | 0.9277 | 0.0454 | 0.7726 | -0.4270 | 0.6694 | False |  | 0.9277 |
| GSE190266 | TACSTD2 | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | gene absent from deposited TPM matrix | NA |
| GSE190266 | CLDN4 | 17 | 52 | 55.8765 | 21.7273 | 5.8298 | 4.5018 | 1.3280 | 0.6584 | 0.3167 | 0.0516 | 0.1973 | 0.1043 | -1.0211 | 0.3072 | True |  | 0.1548 |

## van Elteren
| gene | n_cohorts | van_elteren_Z | van_elteren_p |
| --- | --- | --- | --- |
| TACSTD2 | 1 | -0.4406 | 0.6595 |
| CLDN4 | 2 | 1.6709 | 0.0947 |

## pooled z
| gene | n_DCB | n_NDB | n_cohorts | median_z_DCB | median_z_NDB | AUC_DCB_gt_NDB | rank_biserial | wilcoxon_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TACSTD2 | 14 | 29 | 1 | 0.2423 | 0.1164 | 0.4581 | -0.0837 | 0.6689 |
| CLDN4 | 31 | 81 | 2 | 0.2381 | 0.0095 | 0.5671 | 0.1342 | 0.2744 |

## histology x DCB
| cohort | histology | n_DCB | n_NDB |
| --- | --- | --- | --- |
| GSE190265 | non-squamous | 4 | 10 |
| GSE190265 | squamous | 5 | 6 |
| GSE190265 | unknown | 5 | 13 |
| GSE190266 | non-squamous | 15 | 39 |
| GSE190266 | squamous | 1 | 12 |
| GSE190266 | unknown | 1 | 1 |

## expression vs histology
| cohort | gene | n_nonsquamous | n_squamous | median_log2_nonsquamous | median_log2_squamous | wilcoxon_p |
| --- | --- | --- | --- | --- | --- | --- |
| GSE190265 | TACSTD2 | 14 | 11 | 6.0239 | 6.5441 | 0.3112 |
| GSE190265 | CLDN4 | 14 | 11 | 4.2722 | 3.9574 | 0.1187 |
| GSE190266 | CLDN4 | 54 | 13 | 4.7944 | 4.3870 | 0.3526 |
