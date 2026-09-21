# TACSTD2-high malignant cells are not the IFN-low arm

On the locked concordant-4 malignant pseudobulks, TACSTD2-high is IFN-enriched. Holding CLDN4 constant does not turn that enrichment into a downregulation, and the CLDN4-low half does not show a stable IFN decrease. The path TACSTD2 → CLDN4 → lower IFN is the indirect piece. On the pseudobulk scale its interval includes zero. On the percent-positive scale the indirect piece is negative and the TACSTD2 slope that remains is positive.

Observational. Not a knockdown. The locked T/NK correlation is not re-fit. Cohorts are GSE123902, GSE131907, GSE205335, and GSE189357. Expression n = 64. P4001 is in the percent-positive vector and absent from the UMI sum.

Positive NES means the set sits at the TACSTD2-high end of the ranking. The interval is a patient bootstrap stratified by cohort (B = 300; inner gene-set permutations = 200). The GSEA nominal p conditions on the ranking. Full-data NES uses 1,000 permutations.

## Calibration

The CLDN4 %pos Q4 vs Q1 family score, TMM fit on all 64 units, reproduces the PR #503 numbers.

| family | genes | Q1 / Q4 | logFC | p |
|---|---:|---|---:|---:|
| IFN Hallmark union | 221 | 18 / 16 | −0.584 | 0.00243 |
| MHC-I/APM | 21 | 18 / 16 | −0.779 | 0.00988 |

CLDN4-high is the IFN-low contrast. The TACSTD2 contrasts below are a different exposure on that same matrix.

## Pseudobulk NES

Exposure and covariate are malignant log2(TMM-CPM+1). Residualized rows are the TACSTD2 coefficient after CLDN4 and cohort. The low half is the within-cohort CLDN4 rank at or below the median (n = 34), with no extra CLDN4 term.

| contrast | set | NES | patient bootstrap 95% |
|---|---|---:|---|
| continuous, n = 64 | IFN-γ | +2.32 | +0.74 to +2.96 |
| continuous | IFN-α | +2.60 | +0.65 to +3.01 |
| continuous | MHC-I/APM | +1.55 | −1.61 to +2.35 |
| continuous given CLDN4 | IFN-γ | +2.36 | +1.04 to +3.18 |
| continuous given CLDN4 | IFN-α | +2.66 | +0.63 to +3.26 |
| continuous given CLDN4 | MHC-I/APM | +2.18 | −1.27 to +2.54 |
| CLDN4-low half | IFN-γ | +1.92 | −1.01 to +2.84 |
| CLDN4-low half | IFN-α | +2.04 | −1.25 to +2.97 |
| CLDN4-low half | MHC-I/APM | +1.34 | −1.93 to +2.40 |
| Q4 vs Q1, 19 vs 15 | IFN-γ | +2.31 | −1.73 to +3.36 |

The continuous IFN-γ and IFN-α intervals stay above zero after CLDN4 is in the model. The paired change, residualized NES minus total NES, is +0.15 for IFN-γ (−0.69 to +1.10) and +0.23 for IFN-α (−0.54 to +1.03). MHC-I intervals cover both signs. The Q4 contrast has the same positive point estimate and a patient interval that covers both signs; 15 versus 19 is a thin tail.

Leave-one-cohort-out keeps the continuous IFN-γ NES positive (lowest +1.73 after dropping GSE205335). In that same drop, the CLDN4-low MHC-I NES changes sign (−0.59). IFN-γ does not.

The Q4 IFN-γ leading edge includes CXCL9, CCL2, CD74, and HLA-DQA1 along with IFIT and GBP genes. It is the Hallmark list as ranked, including genes that are also immune-lineage transcripts.

## Percent-positive NES

Cell-level malignant percent positive, same 64 units. This is the sensitivity scale.

| contrast | set | NES | patient bootstrap 95% |
|---|---|---:|---|
| continuous | IFN-γ | +1.86 | −1.80 to +2.90 |
| continuous given CLDN4 | IFN-γ | +3.05 | +1.19 to +3.45 |
| continuous given CLDN4 | IFN-α | +2.89 | +0.92 to +3.35 |
| CLDN4-low half | IFN-γ | −0.85 | −2.79 to +2.04 |
| CLDN4-low half | MHC-I/APM | −1.13 | −2.70 to +1.86 |

On this scale the total continuous IFN-γ interval covers both signs. After CLDN4 is held constant, the IFN-γ and IFN-α intervals stay positive. The paired IFN-γ difference is +1.35 (+0.30 to +3.54). The low-half intervals cover both signs. A negative point estimate in that half is not a stable downregulation.

## Mediation TACSTD2 → CLDN4 → module score

Cohort-adjusted OLS. The outcome is the mean log2(TMM-CPM+1) of the gene set. a is CLDN4 on TACSTD2. b is the module on CLDN4 given TACSTD2. c is the total TACSTD2 slope. c′ is the slope after CLDN4. The product a×b uses a patient bootstrap (B = 4,000).

Pseudobulk slopes are log2 score per log2 TACSTD2. Percent-positive slopes are log2 score per percentage point. VIF of TACSTD2 on CLDN4 plus cohort is 1.74 (pseudobulk) and 1.62 (percent positive).

| scale | outcome | a (p) | b (p) | c (p) | c′ (p) | a×b (95%) |
|---|---|---|---|---|---|---|
| pseudobulk | IFN-γ | +0.458 (8.5×10⁻⁸) | −0.067 (0.27) | +0.077 (0.033) | +0.108 (0.020) | −0.031 (−0.099 to +0.013) |
| pseudobulk | IFN-α | +0.458 | −0.072 (0.33) | +0.078 (0.072) | +0.111 (0.046) | −0.033 (−0.107 to +0.022) |
| pseudobulk | MHC-I | +0.458 | −0.099 (0.29) | +0.055 (0.31) | +0.100 (0.15) | −0.045 (−0.128 to +0.044) |
| percent positive | IFN-γ | +0.491 (1.5×10⁻⁵) | −0.0131 (4.7×10⁻⁴) | +0.0045 (0.16) | +0.0109 (0.0017) | −0.0064 (−0.0116 to −0.0023) |
| percent positive | IFN-α | +0.491 | −0.0135 (0.0029) | +0.0055 (0.14) | +0.0121 (0.0042) | −0.0066 (−0.0130 to −0.0019) |
| percent positive | MHC-I | +0.491 | −0.0169 (0.0032) | +0.0014 (0.77) | +0.0096 (0.066) | −0.0083 (−0.0154 to −0.0019) |

TACSTD2 and CLDN4 move together (pseudobulk Spearman 0.485, p = 4.8×10⁻⁵; percent-positive Spearman 0.617). The indirect product has the CLDN4 sign, toward lower IFN. On the pseudobulk scale that product’s interval includes zero, and the total IFN-γ slope is positive. On the percent-positive scale the product’s interval lies below zero, and the direct IFN-γ slope is positive. The two pieces pull in opposite directions. Their sum is not a TACSTD2-high IFN decrease.

Inside the CLDN4-low half, pseudobulk TACSTD2 still varies (SD 2.37, versus 2.04 overall) and still tracks CLDN4 (Spearman 0.415). The half is not a TACSTD2-restricted slice.

Cohort-adjusted Spearman of pseudobulk TACSTD2 with the IFN-γ score is +0.214 (p = 0.090). After CLDN4 as well, it is +0.246 (p = 0.050). The four-cohort DerSimonian–Laird Spearman is +0.229 (I² = 0, 95% −0.039 to +0.466).

## What this does not claim

The patient bootstrap, not the nominal GSEA p, is the interval for a NES. Gene-set permutation p-values here are small for several positive IFN rankings whose patient intervals still cover zero. GSE189357 contributes two patients to each pseudobulk tail. No cohort was dropped to chase a sign.
