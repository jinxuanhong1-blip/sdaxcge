# Finding — CPTAC LUAD / LSCC protein: CLDN4 and TACSTD2 vs MHC-I, IFN, and CD8A

Public CPTAC TMT freeze v1.2 tumor protein. **LUAD** (Gillette *Cell* 2020, n=110) and **LSCC** (Satpathy *Cell* 2021, n=108) kept separate. Pairwise-complete Spearman ρ, two-sided p, 2,000-resample bootstrap 95% CI (seed `20260921`). Honest n. Do not write n=110 / n=108 for a CLDN4 test.

This page is protein–protein. It does not re-fit CLDN4 protein vs ImmuneScore, GEP18, or CD8A **RNA** (PR #315, #289, #245) and it does not replace the CLDN4 vs CD274 protein page (PR #356). TACSTD2, CD8A **protein**, and the CLDN4 missingness account are the addition.

## Treatment-naive

The freeze phenotype table has **no treatment, neoadjuvant, or ICI column** (LUAD 183 columns, LSCC 193 columns; deconvolution, hallmark ssGSEA, CNV, purity, TMB). Therapy-like column names found: **none**. Hallmark interferon columns in that file are **RNA ssGSEA**, not protein, and are not tested here.

Treatment-naive is the published cohort definition, not a label in this matrix: both studies are prospectively collected, previously untreated surgical resections. There is no ICI arm, no on-treatment biopsy, and no response label. A correlation here cannot be written as acquired resistance or post-ICI change.

## Direction under test

Private PDX proteome numbers are **not** in this repository and are not restated. The human-protein direction that would support an inverse PDX IFN/MHC result is pre-specified: **higher CLDN4 or TACSTD2 protein, lower MHC-I or IFN-core protein**. A pair is called **supports inverse** only when all three hold: bootstrap 95% CI entirely below 0, Benjamini–Hochberg q < 0.05 inside the primary family (2 cohorts × 2 predictors × MHC-I, IFN-core, and CD8A), and the WES-purity partial CI also entirely below 0. CD8A uses the same numeric rule but is an immune-cell protein, not an IFN/MHC score. Null and positive pairs are not support.

## Verdict

| Cohort | Pair | n | ρ | 95% CI | p | q | partial ρ (WES) | Call |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| LUAD | CLDN4 vs MHC-I | 79 | +0.049 | -0.187 to +0.290 | 0.665 | 0.788 | +0.071 | null; does not support |
| LUAD | CLDN4 vs IFN-core | 79 | -0.192 | -0.421 to +0.047 | 0.0894 | 0.268 | -0.154 | null; does not support |
| LUAD | CLDN4 vs CD8A | 79 | -0.098 | -0.320 to +0.148 | 0.391 | 0.586 | -0.044 | null; does not support |
| LUAD | TACSTD2 vs MHC-I | 110 | -0.083 | -0.282 to +0.120 | 0.387 | 0.586 | -0.039 | null; does not support |
| LUAD | TACSTD2 vs IFN-core | 110 | -0.026 | -0.224 to +0.182 | 0.788 | 0.788 | +0.024 | null; does not support |
| LUAD | TACSTD2 vs CD8A | 110 | -0.038 | -0.230 to +0.163 | 0.691 | 0.788 | +0.010 | null; does not support |
| LSCC | CLDN4 vs MHC-I | 78 | -0.410 | -0.574 to -0.211 | 2.0e-04 | 0.00117 | -0.320 | supports inverse |
| LSCC | CLDN4 vs IFN-core | 78 | -0.175 | -0.388 to +0.075 | 0.125 | 0.299 | -0.113 | null; does not support |
| LSCC | CLDN4 vs CD8A | 78 | -0.444 | -0.628 to -0.225 | 4.6e-05 | 5.6e-04 | -0.376 | supports inverse |
| LSCC | TACSTD2 vs MHC-I | 108 | +0.029 | -0.171 to +0.229 | 0.764 | 0.788 | +0.070 | null; does not support |
| LSCC | TACSTD2 vs IFN-core | 108 | +0.230 | +0.053 to +0.408 | 0.0167 | 0.0666 | +0.252 | positive CI; primary-family q ≥ 0.05 |
| LSCC | TACSTD2 vs CD8A | 108 | -0.103 | -0.290 to +0.089 | 0.287 | 0.574 | -0.089 | null; does not support |

**MHC-I, squamous CLDN4 only.** LSCC CLDN4 protein vs the HLA-A/B/C score ρ=-0.410, n=78, 95% CI -0.574 to -0.211, p=2.0e-04, q=0.00117, WES partial ρ=-0.320 (CI -0.523 to -0.083, n=77). Call: supports inverse. The partial CI stays below 0, so this is not removed by WES purity. LUAD CLDN4 vs MHC-I is null (ρ=+0.049, n=79, 95% CI -0.187 to +0.290, p=0.665, q=0.788, WES partial ρ=+0.071 (CI -0.168 to +0.299, n=77)). TACSTD2 vs MHC-I is null in both histologies.

**IFN-core score does not support the inverse direction.** LSCC CLDN4 ρ=-0.175, n=78, 95% CI -0.388 to +0.075, p=0.125, q=0.299, WES partial ρ=-0.113 (CI -0.346 to +0.127, n=77). LUAD CLDN4 ρ=-0.192, n=79, 95% CI -0.421 to +0.047, p=0.0894, q=0.268, WES partial ρ=-0.154 (CI -0.379 to +0.079, n=77). Both confidence intervals include 0. LSCC CLDN4 is inverse at HLA-A and HLA-C and at JAK1, JAK2, IFNGR1, STAT1, IRF9, TAP1, TAP2, TAPBP, PSMB8, PSMB9, and NLRC5 (gene table; each of those CIs lies below 0). HLA-B is the weaker heavy chain (its CI crosses 0). ISG15, MX1, OAS1, and IFI35 are null in that same cohort and pull the 10-gene mean-z score across 0. Do not promote the JAK rows into an IFN-core or IFN-ligand result. IFNG, IFNA1, and IFNB1 protein are absent.

**TACSTD2 is not the CLDN4 result.** It does not track MHC-I or CD8A. LSCC TACSTD2 vs IFN-core is the opposite sign: ρ=+0.230, n=108, 95% CI +0.053 to +0.408, p=0.0167, q=0.0666, WES partial ρ=+0.252 (CI +0.058 to +0.426, n=107). Call: positive CI; primary-family q ≥ 0.05. The partial CI is also entirely above 0. That pair is not inverse support, and it does not clear q < 0.05 in the 12-test primary family. Do not write TACSTD2 protein as a human-protein copy of the LSCC CLDN4–MHC-I inverse.

**CD8A protein, squamous CLDN4 only.** LSCC CLDN4 vs CD8A ρ=-0.444, n=78, 95% CI -0.628 to -0.225, p=4.6e-05, q=5.6e-04, WES partial ρ=-0.376 (CI -0.584 to -0.137, n=77). Call: supports inverse. CD8A protein is quantified in every LSCC tumor (108/108); n=78 is CLDN4 completeness, not CD8A dropout. This is an immune-cell protein association, not an IFN-core result. LUAD CLDN4 vs CD8A is null, and TACSTD2 vs CD8A is null in both histologies. CD8B is absent in LUAD and quantified in 10/108 LSCC tumors, so it is not an endpoint.

Human-protein support for an inverse PDX IFN/MHC direction, on the pre-specified rule, is **LSCC CLDN4 vs MHC-I protein only**. The IFN-core score does not support it. LUAD does not support it. TACSTD2 does not support it. The additional inverse that clears the same bar is LSCC CLDN4 vs CD8A protein.

## CLDN4 missingness — completeness

CLDN4 is on the TMT matrix (`ENSG00000189143.9` in both cohorts). Missing means the abundance cell is empty after numeric parse, not a join failure and not an exact zero. Exact zeros are real log2 ratios equal to the reference and are counted separately. Values are not imputed. Complete-case Spearman uses only quantified CLDN4 tumors.

| Cohort | Protein | Row | Quantified | Missing | % missing | Exact zeros |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| LUAD | CLDN4 | yes | 79 | 31 | 28.2% | 0 |
| LUAD | TACSTD2 | yes | 110 | 0 | 0.0% | 0 |
| LUAD | CD8A | yes | 110 | 0 | 0.0% | 0 |
| LUAD | CD8B | absent | 0 | 110 | 100.0% | 0 |
| LUAD | HLA-A | yes | 110 | 0 | 0.0% | 0 |
| LUAD | HLA-B | yes | 110 | 0 | 0.0% | 0 |
| LUAD | HLA-C | yes | 110 | 0 | 0.0% | 0 |
| LUAD | B2M | absent | 0 | 110 | 100.0% | 0 |
| LUAD | IFNG | absent | 0 | 110 | 100.0% | 0 |
| LUAD | IFNA1 | absent | 0 | 110 | 100.0% | 0 |
| LUAD | IFNB1 | absent | 0 | 110 | 100.0% | 0 |
| LSCC | CLDN4 | yes | 78 | 30 | 27.8% | 0 |
| LSCC | TACSTD2 | yes | 108 | 0 | 0.0% | 0 |
| LSCC | CD8A | yes | 108 | 0 | 0.0% | 0 |
| LSCC | CD8B | yes | 10 | 98 | 90.7% | 0 |
| LSCC | HLA-A | yes | 108 | 0 | 0.0% | 0 |
| LSCC | HLA-B | yes | 108 | 0 | 0.0% | 0 |
| LSCC | HLA-C | yes | 108 | 0 | 0.0% | 0 |
| LSCC | B2M | absent | 0 | 108 | 100.0% | 0 |
| LSCC | IFNG | absent | 0 | 108 | 100.0% | 0 |
| LSCC | IFNA1 | absent | 0 | 108 | 100.0% | 0 |
| LSCC | IFNB1 | absent | 0 | 108 | 100.0% | 0 |

LUAD: CLDN4 quantified **79 / 110** (28.2% missing; range 21.43 to 24.96, median 23.07). TACSTD2 quantified **110 / 110** (0.0% missing). Both quantified 79; CLDN4 only 0; TACSTD2 only 31; both missing 0. LSCC: CLDN4 quantified **78 / 108** (27.8% missing; range 21.53 to 24.73, median 23.34). TACSTD2 quantified **108 / 108** (0.0% missing). Both quantified 78; CLDN4 only 0; TACSTD2 only 30; both missing 0.

The CLDN4 Spearman n is 79 (LUAD) and 78 (LSCC). TACSTD2 and CD8A are quantified in every tumor, so those predictors and that endpoint do not add further dropout. HLA-A/B/C are also complete. IFNG, IFNA1, IFNB1, and B2M are a different kind of gap: the row is absent, not a within-row NA. There is nothing to correlate.

### Do the CLDN4-missing tumors differ?

Mann–Whitney, two-sided, quantified CLDN4 vs missing CLDN4, on tumors that have the endpoint. Delta is median(missing) − median(quantified). A higher MHC-I or IFN score in the missing group would mean the dropped tumors sit toward the high-IFN side and the complete-case ρ is not the full cohort. Minimum n per side = 8.

| Cohort | Endpoint | n quantified vs missing | Median in | Median out | Δ (out−in) | MW p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| LUAD | TACSTD2 | 79 vs 31 | 26.509 | 26.640 | +0.131 | 0.403 |
| LUAD | MHC1_protein | 79 vs 31 | 0.005 | -0.140 | -0.145 | 0.285 |
| LUAD | IFN_core_protein | 79 vs 31 | -0.037 | 0.138 | +0.176 | 0.832 |
| LUAD | CD8A | 79 vs 31 | 21.810 | 21.618 | -0.192 | 0.0687 |
| LUAD | HLA-A | 79 vs 31 | 27.860 | 27.733 | -0.128 | 0.453 |
| LUAD | HLA-B | 79 vs 31 | 26.988 | 26.905 | -0.083 | 0.656 |
| LUAD | HLA-C | 79 vs 31 | 26.132 | 25.964 | -0.167 | 0.232 |
| LUAD | WES_purity | 77 vs 31 | 0.310 | 0.340 | +0.030 | 0.703 |
| LSCC | TACSTD2 | 78 vs 30 | 26.877 | 26.847 | -0.030 | 0.223 |
| LSCC | MHC1_protein | 78 vs 30 | -0.118 | 0.137 | +0.255 | 0.0527 |
| LSCC | IFN_core_protein | 78 vs 30 | -0.068 | 0.094 | +0.162 | 0.157 |
| LSCC | CD8A | 78 vs 30 | 22.381 | 22.687 | +0.307 | 0.0284 |
| LSCC | HLA-A | 78 vs 30 | 27.526 | 27.572 | +0.046 | 0.378 |
| LSCC | HLA-B | 78 vs 30 | 26.916 | 27.021 | +0.105 | 0.112 |
| LSCC | HLA-C | 78 vs 30 | 26.259 | 26.313 | +0.054 | 0.201 |
| LSCC | WES_purity | 77 vs 30 | 0.370 | 0.405 | +0.035 | 0.983 |

CLDN4 dropout is not a purity split and not a TACSTD2 split. WES purity Δ=+0.030, p=0.703 (77 vs 31) in LUAD and Δ=+0.035, p=0.983 (77 vs 30) in LSCC. TACSTD2 Δ=+0.131, p=0.403 (79 vs 31) in LUAD and Δ=-0.030, p=0.223 (78 vs 30) in LSCC. In LUAD the missing tumors are not a distinct MHC-I, IFN-core, or CD8A group (MHC-I Δ=-0.145, p=0.285 (79 vs 31); IFN-core Δ=+0.176, p=0.832 (79 vs 31); CD8A Δ=-0.192, p=0.0687 (79 vs 31)). In LSCC the 30 CLDN4-missing tumors have higher CD8A (Δ=+0.307, p=0.0284 (78 vs 30)) and a same-direction MHC-I shift that does not clear 0.05 (Δ=+0.255, p=0.0527 (78 vs 30)). IFN-core does not differ (Δ=+0.162, p=0.157 (78 vs 30)). Those 30 tumors sit toward higher CD8A, the same direction as the complete-case inverse, not toward a hidden CLDN4-high / IFN-high state. They are still excluded from the Spearman: CLDN4 is not imputed, and the LSCC CLDN4 n stays 78.

TACSTD2 has no missing tumors, so a TACSTD2 quantified-vs-missing contrast does not exist. Full contrast table: `tables/missingness_contrast.tsv`.

## Primary correlations

MHC-I score = mean of per-gene z-scores for HLA-A, HLA-B, and HLA-C (B2M is absent, so the score is 3/4; a tumor needs all three HLA rows). IFN-core = mean of z-scores for the signaling+ISG members with a usable row (LUAD 10/10: STAT1, STAT2, IRF1, IRF9, JAK1, JAK2, ISG15, MX1, OAS1, IFI35; LSCC 10/10: STAT1, STAT2, IRF1, IRF9, JAK1, JAK2, ISG15, MX1, OAS1, IFI35). A tumor needs at least 4 of those members. CD8A is the single protein, not a score. Partial ρ residualizes ranks on WES purity.

WES purity coverage: LUAD 108/110, LSCC 107/108.

Forest: `figures/fig_forest.png`. Scatters: `figures/fig_scatter_cldn4.png`, `figures/fig_scatter_tacstd2.png`.

## Gene-level proteins (supporting, not the support call)

Locked list. Cells are pairwise n, ρ, p. Absent or pairwise n<8 → no ρ. These p-values are **not** in the primary-family FDR. Pairs with n=10 (LSCC CD8B, LSCC IFNAR2) are shown and not interpreted.

### IFN ligands

| Protein | LUAD CLDN4 | LSCC CLDN4 | LUAD TACSTD2 | LSCC TACSTD2 |
| --- | ---: | ---: | ---: | ---: |
| IFNG | absent | absent | absent | absent |
| IFNA1 | absent | absent | absent | absent |
| IFNB1 | absent | absent | absent | absent |

### IFN receptors

| Protein | LUAD CLDN4 | LSCC CLDN4 | LUAD TACSTD2 | LSCC TACSTD2 |
| --- | ---: | ---: | ---: | ---: |
| IFNAR1 | 61, -0.156, 0.229 | 45, -0.422, 0.00386 | 88, +0.003, 0.977 | 60, +0.029, 0.824 |
| IFNAR2 | absent | n=0 | absent | 10, -0.733, 0.0158 |
| IFNGR1 | 79, -0.007, 0.948 | 78, -0.346, 0.00194 | 110, +0.198, 0.038 | 108, +0.234, 0.0148 |
| IFNGR2 | absent | absent | absent | absent |

### IFN signaling

| Protein | LUAD CLDN4 | LSCC CLDN4 | LUAD TACSTD2 | LSCC TACSTD2 |
| --- | ---: | ---: | ---: | ---: |
| STAT1 | 79, -0.218, 0.0541 | 78, -0.254, 0.0247 | 110, -0.125, 0.193 | 108, +0.147, 0.129 |
| STAT2 | 79, -0.136, 0.231 | 78, -0.094, 0.413 | 110, -0.027, 0.78 | 108, +0.100, 0.302 |
| IRF1 | 67, -0.107, 0.388 | 65, -0.219, 0.0791 | 85, -0.144, 0.188 | 90, +0.138, 0.195 |
| IRF9 | 79, -0.283, 0.0115 | 78, -0.281, 0.0128 | 110, -0.180, 0.0599 | 108, +0.263, 0.00597 |
| JAK1 | 79, -0.036, 0.753 | 78, -0.369, 8.8e-04 | 110, +0.256, 0.0069 | 108, -0.002, 0.986 |
| JAK2 | 79, -0.119, 0.297 | 78, -0.419, 1.3e-04 | 110, +0.018, 0.852 | 108, +0.043, 0.658 |

### ISGs

| Protein | LUAD CLDN4 | LSCC CLDN4 | LUAD TACSTD2 | LSCC TACSTD2 |
| --- | ---: | ---: | ---: | ---: |
| ISG15 | 79, -0.055, 0.633 | 78, +0.010, 0.93 | 110, +0.051, 0.599 | 108, +0.182, 0.06 |
| MX1 | 79, -0.037, 0.746 | 78, -0.036, 0.757 | 110, +0.126, 0.19 | 108, +0.239, 0.0126 |
| OAS1 | 79, -0.183, 0.107 | 78, +0.056, 0.624 | 110, -0.058, 0.549 | 108, +0.189, 0.0498 |
| IFI35 | 79, -0.128, 0.26 | 78, +0.015, 0.9 | 110, +0.004, 0.965 | 108, +0.249, 0.00947 |

### MHC-I

| Protein | LUAD CLDN4 | LSCC CLDN4 | LUAD TACSTD2 | LSCC TACSTD2 |
| --- | ---: | ---: | ---: | ---: |
| HLA-A | 79, -0.101, 0.374 | 78, -0.332, 0.00295 | 110, -0.003, 0.979 | 108, -0.064, 0.513 |
| HLA-B | 79, +0.039, 0.734 | 78, -0.227, 0.0454 | 110, -0.008, 0.938 | 108, +0.065, 0.501 |
| HLA-C | 79, +0.071, 0.536 | 78, -0.344, 0.00203 | 110, -0.112, 0.242 | 108, +0.014, 0.882 |
| B2M | absent | absent | absent | absent |

### Antigen processing (extra)

| Protein | LUAD CLDN4 | LSCC CLDN4 | LUAD TACSTD2 | LSCC TACSTD2 |
| --- | ---: | ---: | ---: | ---: |
| TAP1 | 79, -0.107, 0.348 | 78, -0.329, 0.00331 | 110, -0.053, 0.584 | 108, +0.121, 0.214 |
| TAP2 | 79, +0.071, 0.532 | 78, -0.286, 0.0112 | 110, -0.030, 0.754 | 108, +0.134, 0.165 |
| TAPBP | 79, -0.103, 0.367 | 78, -0.333, 0.00288 | 110, -0.093, 0.334 | 108, +0.122, 0.209 |
| PSMB8 | 79, -0.117, 0.306 | 78, -0.323, 0.0039 | 110, -0.260, 0.00612 | 108, +0.090, 0.353 |
| PSMB9 | 79, -0.054, 0.634 | 78, -0.285, 0.0113 | 110, -0.230, 0.0158 | 108, +0.042, 0.668 |
| NLRC5 | 75, +0.083, 0.479 | 78, -0.315, 0.005 | 106, +0.042, 0.669 | 103, +0.165, 0.0951 |

### CD8

| Protein | LUAD CLDN4 | LSCC CLDN4 | LUAD TACSTD2 | LSCC TACSTD2 |
| --- | ---: | ---: | ---: | ---: |
| CD8A | 79, -0.098, 0.391 | 78, -0.444, 4.6e-05 | 110, -0.038, 0.691 | 108, -0.103, 0.287 |
| CD8B | absent | n=5 | absent | 10, +0.055, 0.881 |

## CLDN4 vs TACSTD2 protein

| Cohort | n | ρ | 95% CI | p |
| --- | ---: | ---: | ---: | ---: |
| LUAD | 79 | +0.275 | +0.058 to +0.477 | 0.0143 |
| LSCC | 78 | +0.081 | -0.154 to +0.318 | 0.483 |

Coexpression is context. It is not an IFN/MHC endpoint and it is outside the primary-family FDR.

## What this does not claim

- It does not put a private PDX number into the public table. Support is sign-and-uncertainty on CPTAC protein only.
- It does not test ICI, PD-L1 IHC, or RNA. CD274 protein was null in both histologies on PR #356 and is not re-fit.
- It does not treat IFN **ligand** protein as measured. IFNG, IFNA1, and IFNB1 rows are absent. B2M is absent.
- It does not impute the CLDN4-missing tumors, and it does not write n=110 or n=108 as the CLDN4 n.
- It does not pool LUAD and LSCC.
- Gene-level HLA or JAK rows are not a second discovery set. The support call is the primary scores.
- Hallmark IFN-γ / IFN-α columns in the phenotype file are RNA signatures and are unused.
- WGS purity is stored on the phenotype table and is not a second primary adjustment. WES is the pre-specified residual, matching PR #289.

## Methods

- **Source:** open S3 `cptac-pancancer-data / data_freeze_v1.2_reorganized/{LUAD,LSCC}/`. Filenames from the LinkedOmics CPTAC-pancan index.
- **Protein file:** `{LUAD,LSCC}_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt`. Already log2 vs reference; not logged again.
- **Phenotype file:** `{LUAD,LSCC}_phenotype.txt`, joined on `idx` = protein column name (110/110 and 108/108). Used for the treatment-column audit and `WES_purity`.
- **IDs:** Ensembl prefix (`ENSG` ± version). Map in `genes.py`.
- **Predictors:** CLDN4 and TACSTD2 protein.
- **Primary endpoints:** MHC-I score (HLA-A/B/C), IFN-core score (STAT1, STAT2, IRF1, IRF9, JAK1, JAK2, ISG15, MX1, OAS1, IFI35, dropping any row that is absent or has fewer than 8 quantified tumors), CD8A protein.
- **Gene-level list:** IFN ligands IFNG, IFNA1, IFNB1; receptors IFNAR1, IFNAR2, IFNGR1, IFNGR2; the IFN-core members; MHC-I HLA-A, HLA-B, HLA-C, B2M; antigen processing TAP1, TAP2, TAPBP, PSMB8, PSMB9, NLRC5; CD8B.
- **Missingness:** empty cell = missing. Exact zero counted and kept as a value. Mann–Whitney compares endpoint distributions in predictor-quantified vs predictor-missing tumors.
- **Partial:** rank-transform predictor, endpoint, and WES purity; residualize the two ranks on the purity rank; Pearson of residuals; bootstrap of that partial ρ.
- **FDR:** Benjamini–Hochberg across the primary family only.
- **Not done:** no RNA endpoints, no ImmuneScore re-fit, no imputation, no ICI model, no LUAD+LSCC pool, no post-hoc gene-set search.

## Outputs

- `tables/coverage.tsv` — row present, quantified, missing, % missing, exact zeros, range
- `tables/joint_missingness.tsv` — CLDN4 × TACSTD2 detection
- `tables/missingness_contrast.tsv` — quantified vs missing
- `tables/correlations.tsv` — n, ρ, p, CI, partial, q, call
- `tables/sample_scores.tsv`
- `tables/summary.json`
- `figures/fig_missingness.png`
- `figures/fig_forest.png`
- `figures/fig_scatter_cldn4.png`
- `figures/fig_scatter_tacstd2.png`

```bash
python3 -m pip install -r methods/cptac_ifn_mhc_protein/requirements.txt
python3 methods/cptac_ifn_mhc_protein/download.py --outdir data/cptac_ifn_mhc_protein
python3 methods/cptac_ifn_mhc_protein/analyze.py --data data/cptac_ifn_mhc_protein --outdir methods/cptac_ifn_mhc_protein
```
