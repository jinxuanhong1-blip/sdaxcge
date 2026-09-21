# Protocol (frozen before TCGA association results)

Concordant-4 scRNA reference → TCGA-LUAD bulk deconvolution. This file is the
method. Numbers from the run live in `FINDING.md` and `results/`.

## Question

In TCGA-LUAD primary tumors, how does bulk **CLDN4** associate with an
**estimated CD8 T-cell fraction** built from concordant-4 single-cell
signatures, and how does that association compare with the **keratin-adjusted
bulk** association of CLDN4 with CD8A?

## What the estimators are

- **ν-SVR** is the Newman 2015 CIBERSORT fractions engine (linear kernel,
  ν ∈ {0.25, 0.5, 0.75}, lowest RMSE, negative weights set to 0, remainder
  rescaled to 1). It is the core of CIBERSORTx fraction mode. It is **not**
  the Stanford CIBERSORTx S-mode / B-mode container, which is licensed and is
  not run here.
- **Bayes MAP** is a Dirichlet-multinomial MAP on fixed scRNA reference
  profiles (uniform Dirichlet prior α = 1, simplex via softmax, L-BFGS-B).
  Gene weights are linear TPM. Reference profiles are not updated. This is
  **BayesPrism-style** (mixture + Dirichlet prior + fixed scRNA prior). It is
  **not** the Chu / Danko BayesPrism Gibbs sampler and not InstaPrism.
- **NNLS** is a pre-specified robustness estimator (non-negative least
  squares, same signature, rescaled to 1).

## Reference (concordant-4 only)

GSE123902, GSE131907, GSE205335, GSE189357. No GSE148071, GSE127465,
GSE154826, GSE200563, E-MTAB-13526. No private KL matrices. No mouse data.

Cell classes, fixed order: malignant, cd8, cd4, nk, b, myeloid, endothelial,
fibroblast.

| Dataset | Cells used | Labels |
|---|---|---|
| GSE131907 | `Sample_Origin == tLung` only | Author subtype. Malignant = tS1/tS2/tS3. CD8 = Exhausted CD8+ T, CD8 low T, Cytotoxic CD8+ T, Naive CD8+ T, and only if `Cell_type` is T lymphocytes. CD4 = CD4+ Th, Naive CD4+ T, Treg, Exhausted Tfh (T lymphocytes). NK = `Cell_type` NK cells and subtype NK. B / myeloid / endothelial / fibroblast = author `Cell_type`. Mixed CD8/CD4 and epithelial NA are dropped. |
| GSE205335 | All barcodes that match the UMI matrix | Author `lineage.sub`: Malignant cells, CD8+ T cells, CD4+ T cells, NK cells, B/Plasma cells, Myeloid cells. Endothelial and fibroblast from `lineage.total` (their subtype is NA). Normal epithelial, non-malignant, mast, oligodendrocytes dropped. Cohort includes LN and liver mets (GEO design); they stay in the reference. |
| GSE123902 | Primary tumor and metastasis dense matrices; adjacent normal dropped | Marker gate below. 13 tumor donors, normal-only donor excluded. |
| GSE189357 | All 9 resected samples (AIS/MIA/IAC) | Marker gate below. |

Marker gate (exclusive), used only when there is no author label. A count > 0
is detected.

1. T cell if CD3D or CD3E is detected. CD8 if T and (CD8A or CD8B) and not CD4. CD4 if T and CD4 and not CD8A and not CD8B. If both CD4 and CD8 are detected, assign the side with the larger marker count; ties are dropped. Other T cells are dropped.
2. Else NK if NKG7 or GNLY is detected.
3. Else B if MS4A1 or CD79A.
4. Else myeloid if LYZ, CD68, or CD14.
5. Else endothelial if PECAM1 or VWF.
6. Else fibroblast if COL1A1 or DCN.
7. Else malignant if (EPCAM or KRT19) and PTPRC is undetected.
8. Else drop.

No cell is called malignant because it expresses CLDN4.

Profile for a class in a dataset = pseudobulk CPM
`1e6 * sum(UMI) / sum(library UMI)` over cells of that class. Datasets with
fewer than 30 cells of a class do not contribute that class. The concordant-4
profile is the **unweighted mean** of dataset CPMs (a dataset with a missing
gene is skipped for that gene, not filled with zero). A gene must be present
in at least two datasets to enter the matrix.

## Signature features

Within the concordant-4 CPM matrix, each class contributes up to 120 genes
with CPM ≥ 1 in that class and log2((CPM+1)/(next-best CPM+1)) ≥ 0.585
(~1.5-fold). Union across classes. Always keep lineage markers when they are
in the matrix: CD8A, CD8B, CD3D, CD3E, CD4, FOXP3, NKG7, GZMB, MS4A1, CD79A,
LYZ, CD68, KRT8, KRT18, KRT19, KRT7, EPCAM, PECAM1, COL1A1, COL1A2.

**Held out of every signature (query genes):** CLDN4, TACSTD2.

Before deconvolution, also drop MT-*, RPL*, RPS*, HBA*, HBB*, HBD* and genes
whose max linear TPM in the bulk cohort exceeds 50× the max reference CPM
(BayesPrism-style outlier screen). The same gene list is used for ν-SVR,
Bayes MAP, and NNLS.

**Sensitivity signature:** the same rules, and also hold out KRT8, KRT18,
KRT19, KRT7, EPCAM so the fraction estimate cannot read the keratin/EPCAM
program that the bulk adjustment uses.

**Sensitivity reference:** malignant CPM averaged only from author-labeled
malignant cells (GSE131907 tS1–tS3 and GSE205335). Immune and stromal classes
stay as in the primary reference. Feature rules are the primary rules.

A signature is used only if, in the reference CPM, CD8A is highest in cd8,
MS4A1 or CD79A is highest in b, and EPCAM or KRT19 is highest in malignant.

## Bulk cohort

TCGA-LUAD STAR TPM from the UCSC Xena GDC hub
(`TCGA-LUAD.star_tpm.tsv.gz`), log2(TPM+1) → linear TPM = 2^x − 1.
Ensembl IDs stripped of version and mapped with the GSE189357 10x feature
table (gene symbols, upper case). Duplicate symbols are summed in linear TPM.

Samples: primary solid tumor (barcode sample type `01`). One sample per
patient, preferring vial A (`01A` < `01B`). Patients missing CLDN4, CD8A,
KRT8, KRT18, or KRT19 are dropped.

This is treatment-naive TCGA. It is not an ICI endpoint.

## Association tests

Expression covariates are log2(TPM+1). Fractions are the estimator output
(sum to 1). Spearman is two-sided. Partial Spearman residualizes midranks of
the two variables on midranks of the covariates (intercept included) and
applies a t test with df = n − 2 − k.

Fisher z interval: SE = 1/√(n−3) for Spearman and 1/√(n−3−k) for partial
Spearman.

Primary endpoint: Spearman of CLDN4 vs ν-SVR CD8 fraction.

Keratin-adjusted bulk comparator, same matrix, same patients: partial
Spearman of CLDN4 vs CD8A given KRT8, KRT18, and KRT19.

Also reported, same patients:

- CLDN4 vs ν-SVR CD8 given the ν-SVR malignant fraction
- CLDN4 vs ν-SVR CD8 given KRT8+KRT18+KRT19
- The same three contrasts for Bayes MAP CD8
- CLDN4 vs NNLS CD8 (unadjusted)
- CLDN4 vs ν-SVR (CD8+NK) fraction
- Partial Spearman of CLDN4 vs KRT8 given only KRT18 and KRT19
  (specificity of CLDN4 among keratins; LUAD only)

Benjamini–Hochberg q is computed across those biological tests. The primary
p-value is the unadjusted Spearman p. Controls are not in the BH family:

- Positive: CD8A vs ν-SVR CD8 (expect strong positive)
- Negative: ACTB vs ν-SVR CD8
- Sanity: KRT8 vs ν-SVR malignant fraction (expect positive)
- Concordance: ν-SVR CD8 vs Bayes MAP CD8

Fit QC: per-sample Pearson r of observed vs reconstructed linear TPM on the
signature genes. Median r is reported. Fractions are still computed if fit
is modest; the finding states the median r.

## Non-goals

No seven-cohort keratin re-analysis. No claim that this replaces the locked
keratin table. No spatial language. No merging with private 8KL data.
