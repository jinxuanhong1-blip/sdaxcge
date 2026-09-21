# Finding — GTEx + TCGA pan-cancer: CLDN4 vs PRKDC / LIG4 / STING1 / HLA-A/B/C

**Additive public bulk RNA.** UCSC Xena Toil `TcgaTargetGtex_rsem_gene_tpm`, log2(RSEM TPM+0.001), GENCODE v23, one pipeline for TCGA and GTEx. STING1 is **TMEM173** on this freeze (`ENSG00000184584.12`). This folder does not re-audit CosMx exclusion, concordant-4, the keratin partials, or TISMO.

Primary question, pre-specified: within-cohort Spearman of **CLDN4** vs **PRKDC**, **LIG4**, **STING1**, **HLA-A**, **HLA-B**, and **HLA-C** in **TCGA-LUAD primary**, **TCGA-LUSC primary**, and **GTEx lung**. Those 18 tests are one BH family. Other TCGA primaries and other GTEx tissues are consistency context with their own BH families. Cohorts are not pooled into one correlation.

## Honest n

| Cohort | Samples before donor collapse | Donors (analysis n) | Donors with >1 aliquot |
|---|---:|---:|---:|
| TCGA-LUAD primary | 513 | **513** | 0 |
| TCGA-LUSC primary | 498 | **498** | 0 |
| GTEx lung | 288 | **266** | 14 |
| TCGA-LUAD solid-tissue normal | 59 | 59 | 0 |
| TCGA-LUSC solid-tissue normal | 50 | 50 | 0 |
| TCGA-LUAD HiSeqV2 primary (sensitivity) | 515 | 515 | 0 |
| TCGA-LUSC HiSeqV2 primary (sensitivity) | 502 | 502 | 0 |

Analysis n is donors after averaging replicate aliquots. No imputation. Pan-cancer context uses 32 TCGA solid-primary cohorts plus LAML (blood primary) and SKCM metastatic (most SKCM RNA is metastatic; primary SKCM is its own solid row). GTEx context uses 46 normal tissues with at least 30 donors. TARGET samples are not included. GTEx cell lines are not included.

## Lung primary result

BH q is across these 18 tests only. A q on this table is not the pan-cancer q.

| Cohort | Endpoint | n | ρ | 95% CI | p | q |
|---|---|---:|---:|---|---:|---:|
| TCGA-LUAD | PRKDC | 513 | +0.171 | [+0.086, +0.254] | 9.50e-05 | 4.27e-04 |
| TCGA-LUAD | LIG4 | 513 | +0.051 | [-0.036, +0.137] | 0.248 | 0.410 |
| TCGA-LUAD | STING1 | 513 | +0.316 | [+0.236, +0.392] | 2.36e-13 | 4.24e-12 |
| TCGA-LUAD | HLA-A | 513 | +0.156 | [+0.071, +0.240] | 3.81e-04 | 0.001 |
| TCGA-LUAD | HLA-B | 513 | +0.101 | [+0.015, +0.186] | 0.021 | 0.055 |
| TCGA-LUAD | HLA-C | 513 | +0.099 | [+0.012, +0.184] | 0.025 | 0.057 |
| TCGA-LUSC | PRKDC | 498 | -0.028 | [-0.116, +0.060] | 0.533 | 0.639 |
| TCGA-LUSC | LIG4 | 498 | +0.029 | [-0.059, +0.117] | 0.513 | 0.639 |
| TCGA-LUSC | STING1 | 498 | +0.035 | [-0.053, +0.123] | 0.435 | 0.602 |
| TCGA-LUSC | HLA-A | 498 | +0.018 | [-0.070, +0.106] | 0.689 | 0.730 |
| TCGA-LUSC | HLA-B | 498 | -0.026 | [-0.113, +0.062] | 0.570 | 0.641 |
| TCGA-LUSC | HLA-C | 498 | -0.006 | [-0.094, +0.082] | 0.889 | 0.889 |
| GTEx-Lung | PRKDC | 266 | -0.067 | [-0.186, +0.053] | 0.273 | 0.410 |
| GTEx-Lung | LIG4 | 266 | -0.070 | [-0.189, +0.050] | 0.252 | 0.410 |
| GTEx-Lung | STING1 | 266 | -0.131 | [-0.248, -0.011] | 0.032 | 0.065 |
| GTEx-Lung | HLA-A | 266 | -0.238 | [-0.348, -0.121] | 8.85e-05 | 4.27e-04 |
| GTEx-Lung | HLA-B | 266 | -0.352 | [-0.453, -0.242] | 3.45e-09 | 3.11e-08 |
| GTEx-Lung | HLA-C | 266 | -0.162 | [-0.276, -0.042] | 0.008 | 0.025 |

## What holds

The three lung cohorts are not one correlation. LUSC is null. LUAD and GTEx lung point different ways, and the LUAD NHEJ estimates change sign on the HiSeqV2 freeze.

**LUSC primary, Toil** (n=498). The largest |ρ| among the six is 0.035. Every lung-family q is above 0.05.

**LUAD STING1 is the tumor association that is positive on both RNA freezes.** Toil ρ=+0.316 (n=513, p=2.36e-13, q=4.24e-12). It stays after PTPRC (+0.358) and after EPCAM (+0.357). STING1 itself tracks PTPRC (ρ=+0.298). The PTPRC partial stays positive and is slightly larger than the unadjusted ρ. HiSeqV2 LUAD STING1 ρ=+0.193 (n=515, q=6.23e-05 in that separate 12). LUSC STING1 is null on both freezes (Toil +0.035, HiSeqV2 -0.005).

**LUAD PRKDC on Toil is small, and the HiSeqV2 freeze is negative.** Toil ρ=+0.171 (q=4.27e-04). PTPRC leaves it (+0.169); EPCAM shrinks it to +0.062 (EPCAM-partial q=0.244). HiSeqV2 LUAD PRKDC is -0.126 (q=0.013). LUSC Toil PRKDC is -0.028; HiSeqV2 LUSC PRKDC is -0.223.

**LIG4 is null in the Toil lung 18.** LUAD +0.051, LUSC +0.029, GTEx lung -0.070, all lung-family q above 0.05. HiSeqV2 LUAD LIG4 is -0.168 (q=5.13e-04). The Toil lung estimates for LIG4 sit near zero; the HiSeqV2 LUAD estimate is negative.

**GTEx lung HLA is negative on Toil, and most of that is shared with PTPRC.** HLA-A -0.238, HLA-B -0.352, HLA-C -0.162, n=266, all three lung-family q below 0.05. In the same donors CLDN4 vs PTPRC is -0.366, and HLA-B vs PTPRC is +0.659. After PTPRC the partials are HLA-A -0.031, HLA-B -0.159 (partial q=0.029), HLA-C +0.017. Read that as bulk composition in normal lung: CLDN4 against leukocyte HLA. The LUAD tumor HLA estimates are small and positive. A bulk ρ is not a spatial neighborhood.

**LUAD HLA on Toil is a small positive. HiSeqV2 HLA is null, and that null matches the earlier HiSeqV2 table.** Toil LUAD HLA-A +0.156 (q=0.001), HLA-B +0.101 (q=0.055), HLA-C +0.099 (q=0.057). HiSeqV2 LUAD n=515: HLA-A +0.059 (p=0.180), HLA-B +0.035 (p=0.424), HLA-C +0.023 (p=0.596). Quote both freezes for HLA. The HiSeqV2 unadjusted result is null.

No lung primary cohort has more than 10% of donors at the log2(TPM+0.001) floor for CLDN4 or for any of the six endpoints.

## Composition sensitivity (lung 18, separate BH)

Partials are rank-residual Pearson on the same Toil matrix. **PTPRC** is a leukocyte marker. **EPCAM** is an epithelial marker. Neither covariate is inside the endpoint definitions. Residualizing RNA does not equal adjusting for DNA purity. An EPCAM partial that shrinks is shared epithelial variance, not by itself proof that the unadjusted association was only composition. q below is BH within the 18 partials of that covariate, not the unadjusted 18.

CLDN4 versus the two covariates (context for the partials):

| Cohort | n | CLDN4 vs PTPRC ρ | p | CLDN4 vs EPCAM ρ | p |
|---|---:|---:|---:|---:|---:|
| TCGA-LUAD | 513 | -0.082 | 0.065 | +0.382 | 2.82e-19 |
| TCGA-LUSC | 498 | -0.052 | 0.250 | +0.293 | 2.66e-11 |
| GTEx-Lung | 266 | -0.366 | 7.32e-10 | +0.763 | 6.23e-52 |

**PRKDC.** TCGA-LUAD unadj +0.171 → PTPRC +0.169 (q=5.61e-04), EPCAM +0.062 (q=0.244), both +0.060; TCGA-LUSC unadj -0.028 → PTPRC -0.037 (q=0.607), EPCAM -0.078 (q=0.163), both -0.077; GTEx-Lung unadj -0.067 → PTPRC -0.055 (q=0.607), EPCAM -0.252 (q=1.22e-04), both -0.239.

**LIG4.** TCGA-LUAD unadj +0.051 → PTPRC +0.080 (q=0.185), EPCAM +0.007 (q=0.872), both -0.002; TCGA-LUSC unadj +0.029 → PTPRC +0.043 (q=0.607), EPCAM +0.027 (q=0.585), both +0.025; GTEx-Lung unadj -0.070 → PTPRC +0.023 (q=0.850), EPCAM -0.308 (q=2.73e-06), both -0.257.

**STING1.** TCGA-LUAD unadj +0.316 → PTPRC +0.358 (q=1.24e-15), EPCAM +0.357 (q=1.53e-15), both +0.366; TCGA-LUSC unadj +0.035 → PTPRC +0.061 (q=0.400), EPCAM +0.164 (q=7.24e-04), both +0.171; GTEx-Lung unadj -0.131 → PTPRC +0.008 (q=0.895), EPCAM +0.088 (q=0.244), both +0.176.

**HLA-A.** TCGA-LUAD unadj +0.156 → PTPRC +0.194 (q=8.42e-05), EPCAM +0.184 (q=1.21e-04), both +0.186; TCGA-LUSC unadj +0.018 → PTPRC +0.050 (q=0.540), EPCAM +0.084 (q=0.138), both +0.089; GTEx-Lung unadj -0.238 → PTPRC -0.031 (q=0.796), EPCAM -0.061 (q=0.361), both +0.077.

**HLA-B.** TCGA-LUAD unadj +0.101 → PTPRC +0.169 (q=5.61e-04), EPCAM +0.187 (q=1.20e-04), both +0.201; TCGA-LUSC unadj -0.026 → PTPRC +0.009 (q=0.895), EPCAM +0.056 (q=0.275), both +0.061; GTEx-Lung unadj -0.352 → PTPRC -0.159 (q=0.029), EPCAM -0.093 (q=0.239), both +0.060.

**HLA-C.** TCGA-LUAD unadj +0.099 → PTPRC +0.146 (q=0.003), EPCAM +0.146 (q=0.002), both +0.148; TCGA-LUSC unadj -0.006 → PTPRC +0.022 (q=0.796), EPCAM +0.061 (q=0.245), both +0.063; GTEx-Lung unadj -0.162 → PTPRC +0.017 (q=0.886), EPCAM -0.072 (q=0.288), both +0.035.

## HiSeqV2 lung freeze (sensitivity)

Legacy UCSC Xena HiSeqV2 log2(norm_count+1), primary tumors (`-01`), donor-averaged. This is the freeze earlier LUAD/LUSC HLA tables used. It is not the Toil TPM matrix. BH here is a separate 12-test family (2 cohorts × 6 genes).

| Cohort | Endpoint | n | ρ | 95% CI | p | q |
|---|---|---:|---:|---|---:|---:|
| TCGA-LUAD-HiSeqV2 | PRKDC | 515 | -0.126 | [-0.210, -0.040] | 0.004 | 0.013 |
| TCGA-LUAD-HiSeqV2 | LIG4 | 515 | -0.168 | [-0.251, -0.083] | 1.28e-04 | 5.13e-04 |
| TCGA-LUAD-HiSeqV2 | STING1 | 515 | +0.193 | [+0.108, +0.275] | 1.04e-05 | 6.23e-05 |
| TCGA-LUAD-HiSeqV2 | HLA-A | 515 | +0.059 | [-0.027, +0.145] | 0.180 | 0.359 |
| TCGA-LUAD-HiSeqV2 | HLA-B | 515 | +0.035 | [-0.051, +0.121] | 0.424 | 0.706 |
| TCGA-LUAD-HiSeqV2 | HLA-C | 515 | +0.023 | [-0.063, +0.110] | 0.596 | 0.715 |
| TCGA-LUSC-HiSeqV2 | PRKDC | 502 | -0.223 | [-0.305, -0.139] | 4.24e-07 | 5.09e-06 |
| TCGA-LUSC-HiSeqV2 | LIG4 | 502 | -0.086 | [-0.172, +0.002] | 0.055 | 0.133 |
| TCGA-LUSC-HiSeqV2 | STING1 | 502 | -0.005 | [-0.092, +0.083] | 0.915 | 0.915 |
| TCGA-LUSC-HiSeqV2 | HLA-A | 502 | +0.032 | [-0.055, +0.119] | 0.471 | 0.706 |
| TCGA-LUSC-HiSeqV2 | HLA-B | 502 | +0.006 | [-0.081, +0.094] | 0.885 | 0.915 |
| TCGA-LUSC-HiSeqV2 | HLA-C | 502 | +0.025 | [-0.062, +0.113] | 0.570 | 0.715 |

Same genes, two quantifications. Toil is log2(TPM+0.001). HiSeqV2 is log2(norm_count+1).

| Endpoint | Toil LUAD | HiSeqV2 LUAD | Toil LUSC | HiSeqV2 LUSC |
|---|---:|---:|---:|---:|
| PRKDC | +0.171 | -0.126 | -0.028 | -0.223 |
| LIG4 | +0.051 | -0.168 | +0.029 | -0.086 |
| STING1 | +0.316 | +0.193 | +0.035 | -0.005 |
| HLA-A | +0.156 | +0.059 | +0.018 | +0.032 |
| HLA-B | +0.101 | +0.035 | -0.026 | +0.006 |
| HLA-C | +0.099 | +0.023 | -0.006 | +0.025 |

STING1 in LUAD is positive on both. PRKDC and LIG4 in LUAD change sign. HLA-A/B/C in LUAD are null on HiSeqV2 and only HLA-A clears the Toil 18-test BH. LUSC HLA is null on both. Report both freezes.

HiSeqV2 LUAD HLA-A/B/C match the earlier unadjusted LUAD HiSeqV2 estimates (+0.059 / +0.035 / +0.023) within 0.001. That earlier table joined ESTIMATE and used the same HiSeqV2 freeze.

## TCGA pan-cancer context

Each solid-primary cohort is its own Spearman (n≥30 to enter BH). Exploratory BH is across TCGA cohorts × the six genes, and it includes LUAD, LUSC, LAML blood, and SKCM metastatic as separate rows. Median and IQR below are **solid primaries only** (not LAML, not SKCM metastatic). This is a description of heterogeneity, not a single pan-cancer p-value.

| Endpoint | Cohorts | Median ρ | IQR | Positive | Exploratory q<0.05 and + | Exploratory q<0.05 and − |
|---|---:|---:|---|---:|---:|---:|
| PRKDC | 32 | +0.165 | [+0.028, +0.256] | 25/32 | 16 | 1 |
| LIG4 | 32 | +0.119 | [-0.014, +0.224] | 23/32 | 14 | 1 |
| STING1 | 32 | +0.178 | [+0.032, +0.278] | 26/32 | 15 | 0 |
| HLA-A | 32 | +0.105 | [+0.000, +0.170] | 24/32 | 11 | 2 |
| HLA-B | 32 | +0.079 | [-0.056, +0.202] | 21/32 | 10 | 3 |
| HLA-C | 32 | +0.131 | [-0.018, +0.203] | 23/32 | 12 | 1 |

Solid-primary cohort count in that table: 32. Full per-cohort ρ is in `tables/correlations.tsv` and `figures/heatmap_tcga.png`.

## GTEx tissue context

Normal tissues with at least 30 donors. Cell lines are excluded. Brain subregions stay separate. BH is within this GTEx family. Lung is included here as one tissue; its primary q is the lung-family q above, not this one.

| Endpoint | Cohorts | Median ρ | IQR | Positive | Exploratory q<0.05 and + | Exploratory q<0.05 and − |
|---|---:|---:|---|---:|---:|---:|
| PRKDC | 46 | +0.077 | [-0.067, +0.352] | 28/46 | 18 | 6 |
| LIG4 | 46 | +0.087 | [-0.092, +0.324] | 29/46 | 17 | 7 |
| STING1 | 46 | +0.170 | [-0.078, +0.315] | 31/46 | 23 | 8 |
| HLA-A | 46 | +0.187 | [+0.058, +0.319] | 39/46 | 23 | 3 |
| HLA-B | 46 | +0.184 | [+0.096, +0.339] | 38/46 | 23 | 3 |
| HLA-C | 46 | +0.179 | [+0.069, +0.258] | 41/46 | 20 | 3 |

Full per-tissue ρ is in `tables/correlations.tsv` and `figures/heatmap_gtex.png`.

## Lung adjacent normal

TCGA solid-tissue normal, not pooled with GTEx and not part of the 18. These donors overlap the tumor cohorts as paired normals; the correlation is within the normal samples, not a tumor-minus-normal delta. BH is a separate 12-test family.

| Cohort | Endpoint | n | ρ | 95% CI | p | q |
|---|---|---:|---:|---|---:|---:|
| TCGA-LUAD-normal | PRKDC | 59 | +0.072 | [-0.188, +0.322] | 0.589 | 0.785 |
| TCGA-LUAD-normal | LIG4 | 59 | +0.388 | [+0.147, +0.586] | 0.002 | 0.028 |
| TCGA-LUAD-normal | STING1 | 59 | +0.215 | [-0.044, +0.446] | 0.103 | 0.247 |
| TCGA-LUAD-normal | HLA-A | 59 | -0.130 | [-0.374, +0.130] | 0.325 | 0.558 |
| TCGA-LUAD-normal | HLA-B | 59 | -0.140 | [-0.382, +0.121] | 0.291 | 0.558 |
| TCGA-LUAD-normal | HLA-C | 59 | +0.009 | [-0.248, +0.264] | 0.946 | 0.966 |
| TCGA-LUSC-normal | PRKDC | 50 | -0.247 | [-0.491, +0.034] | 0.084 | 0.247 |
| TCGA-LUSC-normal | LIG4 | 50 | +0.241 | [-0.040, +0.487] | 0.091 | 0.247 |
| TCGA-LUSC-normal | STING1 | 50 | +0.281 | [+0.003, +0.519] | 0.048 | 0.247 |
| TCGA-LUSC-normal | HLA-A | 50 | +0.006 | [-0.273, +0.284] | 0.966 | 0.966 |
| TCGA-LUSC-normal | HLA-B | 50 | +0.041 | [-0.240, +0.316] | 0.777 | 0.932 |
| TCGA-LUSC-normal | HLA-C | 50 | +0.121 | [-0.163, +0.386] | 0.402 | 0.603 |

The only adjacent test with q<0.05 in that 12 is LUAD normal vs LIG4 (ρ=+0.388, n=59, q=0.028). Tumor LIG4 is null (LUAD +0.051, LUSC +0.029). The normal correlation is not a substitute for the tumor result.

## Companions (not in any FDR family)

NHEJ_mean is the unweighted mean of PRKDC and LIG4 on the log2 TPM scale. MHC_I is the unweighted mean of HLA-A, HLA-B, and HLA-C. Neither is a gene-set enrichment.

| Cohort | Companion | n | ρ | p |
|---|---|---:|---:|---:|
| TCGA-LUAD | NHEJ_mean | 513 | +0.156 | 3.74e-04 |
| TCGA-LUAD | MHC_I | 513 | +0.129 | 0.004 |
| TCGA-LUSC | NHEJ_mean | 498 | -0.016 | 0.717 |
| TCGA-LUSC | MHC_I | 498 | -0.009 | 0.834 |
| GTEx-Lung | NHEJ_mean | 266 | -0.134 | 0.029 |
| GTEx-Lung | MHC_I | 266 | -0.286 | 2.09e-06 |

## What this measurement is

These are within-cohort rank correlations on bulk RNA. A positive CLDN4–PRKDC or CLDN4–LIG4 ρ means the two transcripts move together across donors in that tissue or tumor type. A CLDN4–HLA or CLDN4–STING1 ρ mixes tumor-cell and immune-cell RNA. It is not a tumor-cell program, not a spatial neighborhood, and not an ICI endpoint. Negative HLA ρ is not spatial exclusion. Positive NHEJ ρ is not a repair mechanism.

## Methods

- **Toil matrix:** `TcgaTargetGtex_rsem_gene_tpm`, log2(TPM+0.001). Phenotype `TcgaTargetGTEX_phenotype.txt.gz`. Gene map `gencode.v23.annotation.gene.probemap` (Ensembl ID with version). Zeros are stored at log2(0.001) ≈ −9.966; Spearman treats that floor as a tie.
- **Cohorts:** TCGA `_sample_type == Primary Tumor` split by `detailed_category`, plus LAML blood primary and SKCM metastatic as labeled extra rows. GTEx `_sample_type == Normal Tissue` split by `detailed_category`. Lung adjacent: TCGA `_sample_type == Solid Tissue Normal` for LUAD and LUSC only.
- **Donor collapse:** TCGA patient = first three barcode fields; GTEx donor = first two. Mean of aliquots. Spearman is unchanged by a monotone transform, so log2 TPM and TPM agree.
- **Unadjusted:** Spearman, two-sided. **CI:** Fisher z, variance 1/(n−3). Tests with n<10 are not computed. BH uses cohorts with n≥30.
- **Partial:** rank-residual Pearson given PTPRC, EPCAM, or both. t degrees of freedom n−k−2. Fisher z variance 1/(n−k−3). Algebraic first-order partial Spearman was the formula in the earlier ImmuneScore lung table; rank-residual Pearson is the sensitivity used here because it extends to two covariates. On a single covariate the two agree closely when ties are mild.
- **FDR families, kept separate:** (1) lung Toil 18, (2) TCGA exploratory cohorts × 6, (3) GTEx tissues × 6, (4) lung adjacent 12, (5) HiSeqV2 lung 12, (6) lung PTPRC partials 18, (7) lung EPCAM partials 18. Companions are outside every family.
- **HiSeqV2:** `TCGA.LUAD.sampleMap/HiSeqV2` and `TCGA.LUSC.sampleMap/HiSeqV2`, primary `01` only, collapsed to the 15-character barcode by mean, then to the patient.

## What is not done

- No pooled GTEx+TCGA correlation. No STAR-TPM pan-cancer re-run. No ESTIMATE or leukocyte-fraction residual (PTPRC/EPCAM are the pre-specified RNA covariates). No ICI outcome, protein, spatial statistic, or single-cell pseudobulk. No claim that this ranks CLDN4 against other surface genes.

## Files

- `tables/correlations.tsv` — every cohort, endpoint, partial, and q
- `tables/sign_summary.tsv` — median ρ and sign counts
- `tables/counts.tsv` — honest n per cohort
- `tables/samples.tsv` — donor-level matrix the correlations were computed from
- `tables/provenance.json` — URLs, sha256, Ensembl IDs
- `figures/forest_lung.png` — lung unadjusted and PTPRC-partial ρ
- `figures/heatmap_tcga.png` — TCGA per-cohort ρ
- `figures/heatmap_gtex.png` — GTEx per-tissue ρ
- Reproduce: `python3 methods/gtex_tcga_cldn4_nhej/analyze.py`
