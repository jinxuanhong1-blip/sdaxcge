# Methods

## Cohorts

Public processed single-cell matrices, concordant four only.

| Accession | What was downloaded | Malignant call | T/NK call | Unit |
|---|---|---|---|---|
| GSE123902 | GEO RAW dense CSV | EPCAM or KRT8/18/19 detected and PTPRC = 0 | CD3D, CD3E, CD8A, NKG7, GNLY, or KLRD1 detected, and not malignant | Donor. Primary tumor preferred over metastasis. Normal lung dropped. |
| GSE131907 | GEO raw UMI text + author annotation | Cell_subtype in {Malignant cells, tS1, tS2, tS3} | Cell_type in {T lymphocytes, NK cells} | Sample id from the locked 21-sample inventory (mLN, tL/B, mBrain with malignant cells). |
| GSE205335 | GEO RDS UMI (double-gzip) + CellIdentity | lineage.sub = Malignant cells | lineage.total = T/NK cells | Patient. Multiple GSMs pooled. Tissue starting with "Normal" dropped. |
| GSE189357 | GEO 10x MTX | Same marker rule as GSE123902 | Same marker rule as GSE123902 | Patient (TD1–TD9, tumor). |

Not used: GSE148071, GSE127465, GSE154826, GSE200563, E-MTAB-13526, GSE207422. TACSTD2 is never a gate.

## Normalization and the within-patient split

Counts are `log1p(1e4 * UMI / library size)`.

Primary split, inside each unit: rank malignant cells by CLDN4 log expression with ties broken by first occurrence (R `rank(..., ties.method="first")`). High = rank above the 75th percentile. Low = rank at or below the 25th percentile. A unit is eligible for the quartile contrast when it has at least 40 malignant cells and 20 T/NK cells, the quartile split is defined, mean CLDN4 is higher in the high arm than the low arm, and at least 25% of malignant cells have a CLDN4 UMI. The 25% floor keeps a zero-tie from filling the high arm. Units below that floor remain in the percent-positive split when both positive and negative cells exist.

Extra splits on the same eligible units: above vs at-or-below the unit median, and CLDN4-positive vs CLDN4-zero cells.

Sender Δ for a ligand is mean log expression in the high arm minus mean log expression in the low arm. The family score for a unit is the mean of those ligand Δs. The test is a Wilcoxon signed-rank against 0 across units. Cohort means are also combined with DerSimonian–Laird random effects; I² is reported from that meta-analysis.

## Receiver geneset

Empirical activity uses author-annotated T/NK only (GSE131907 and GSE205335). Marker-gated T/NK in GSE123902 and GSE189357 is not used to choose the geneset. Those two cohorts still contribute sender splits.

For each author cohort, Spearman correlation of T/NK mean log expression with malignant CLDN4 percent positive. The two correlations are Fisher-z combined with DerSimonian–Laird. A gene must have the same sign in both cohorts.

Benjamini–Hochberg FDR ≤ 0.10 is the primary cutoff. If fewer than 15 genes pass, the unadjusted p ≤ 0.01 set is used and labeled as a fallback. If that is still smaller than 15, the top 100 genes by absolute meta z are used and labeled as a fallback. In this run the p ≤ 0.01 fallback was used. Genes also have to sit in the T/NK expressed background before activity is scored, which dropped the up-set from 69 genes to 12.

Lung epithelial, secretory, and ciliated genes are removed before the cutoff: EPCAM, TACSTD2, CDH1, NKX2-1, FOXA2, GATA6, ELF3, SOX2, HNF1B, KRT5/7/8/17/18/19, CLDN3/4/7, MUC1/4/5B/16, CEACAM5/6, WFDC2, SCGB1A1, SCGB3A1/2, BPIFA1/B1, SFTPA1/A2, SFTPB/C/D, SFTA2/3, AGER, AQP1/4, FOXJ1, TMC5, ARMC3, EMP2, RAB25, GRB7, NAPSA, AGR2. An all-four correlation table is written and is not the activity geneset.

## NicheNet activity

Prior: human NicheNet v2 ligand–target matrix and ligand–receptor network, Zenodo record 7074291 (`ligand_target_matrix_nsga2r_final.rds`, `lr_network_human_21122021.rds`).

Code: the functions `predict_ligand_activities`, `convert_gene_list_settings_evaluation`, `convert_settings_ligand_prediction`, `get_single_ligand_importances`, `evaluate_target_prediction`, and `classification_evaluation_continuous_pred` sourced from saeyslab/nichenetr commit `66f90d5eeafef280b2b2f339b3fd70ffec1781dd`. AUPR uses ROCR precision-recall curves and `caTools::trapz`. AUPR-corrected is AUPR minus the fraction of background genes that sit in the geneset. That is the rank used in the NicheNet vignette.

Background gene: measured in a cohort, and detected in at least 10% of T/NK cells in at least 10% of eligible units of the cohorts where it was measured.

Potential ligand: a column of the ligand–target matrix, detected in at least 10% of malignant cells in at least 10% of eligible units, with at least one receptor in the v2 ligand–receptor network detected in T/NK at that same threshold.

Activity is scored for the empirical up and down genesets and for three a priori receiver programs (IFN, cytotoxicity, exhaustion). A continuous Pearson of regulatory potential against the meta-z vector is stored as a companion and is not the official binary activity.

The ligand-level contrast of core barrier (F11R, NECTIN2, CDH1, LGALS9) versus core IFN/recruit (CXCL9, CXCL10, CCL5, IFNG) is an exact reassignment of those labels. It is reported separately from the patient-level Wilcoxon on sender Δ.

A joint priority score on the empirical-up geneset is the sum of z-scored AUPR-corrected, z-scored sender Δ, and z-scored best-receptor detection. It is a differential-NicheNet-style summary, not a fit of the `multinichenetr` package.

## Pre-specified families

- Core barrier / inhibitory: F11R, NECTIN2, CDH1, LGALS9. Expect higher in CLDN4-high malignant cells.
- Extended barrier: PVR, NECTIN1, NECTIN3, CD274, PDCD1LG2, CEACAM1, TGFB1, CD47, HLA-E, MIF, CD24, LGALS3.
- Core IFN / recruit: CXCL9, CXCL10, CCL5, IFNG. Expect higher in CLDN4-low malignant cells.
- Extended IFN / recruit: CXCL11, CXCL16, CCL2, CCL4, CXCL12, ICAM1, TNFSF9, IL15, IL18.
- MHC-I, reported on its own: HLA-A, HLA-B, HLA-C.

PVRL2 is read as NECTIN2 and JAM1 as F11R when a matrix still uses the old symbol.
