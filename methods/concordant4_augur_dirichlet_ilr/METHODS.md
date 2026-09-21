# Methods

Additive layer on the locked concordant-4 (GSE123902, GSE131907, GSE205335, GSE189357). The malignant CLDN4 percent-positive score and the within-cohort quartile are copied from the locked patient table (n = 65). They are not recomputed, and no fifth cohort is added.

The unit is the donor (GSE123902), tumor-bearing sample (GSE131907), or patient (GSE205335, GSE189357). Cell counts are full-unit counts. The Harmony cap of 350 cells per unit is not used.

## Lineages

GSE131907 uses author `Cell_type` / `Cell_subtype`. T = T lymphocytes, NK = NK cells, myeloid = Myeloid cells, malignant = subtype Malignant cells. T + NK equals the locked T/NK count on every sample.

GSE205335 uses author `lineage.total` / `lineage.sub`. Samples whose tissue starts with Normal are dropped, matching the locked denominator. NK = lineage.sub NK cells. T = the rest of lineage.total T/NK cells. Myeloid = lineage.total Myeloid cells. Malignant = lineage.sub Malignant cells. T + NK equals the locked T/NK count on every patient.

GSE123902 and GSE189357 have no author lineage table. Cells are called from raw counts with a mutually exclusive marker gate:

- malignant: (EPCAM or KRT8 or KRT18 or KRT19) > 0 and PTPRC = 0
- T: not malignant, and (CD3D or CD3E or CD3G or CD4 or CD8A) > 0
- NK: not malignant or T, and (NKG7 or GNLY or KLRD1 or NCAM1 or NCR1) > 0
- myeloid: not malignant, T, or NK, and (LYZ or CD68 or CD14 or FCGR3A or CD163 or C1QA or MARCO or AIF1) > 0

The narrower locked T/NK gate, (CD3D or CD3E or CD8A or NKG7 or GNLY or KLRD1) > 0 and not malignant, is recomputed only as a check. It matches the locked count on all 22 marker-cohort units. The T and NK counts used in the models also include CD4 / CD3G T cells and NCAM1 / NCR1 NK cells that the locked gate left in the remainder. That is why T + NK is larger than the locked T/NK count on those two cohorts only.

Rest = all cells minus T minus NK minus myeloid. It contains the malignant cells.

## ILR

Counts get a half-cell pseudocount, then close. No unit had a zero T, NK, or myeloid count, so the pseudocount does not create a part.

The 4-part simplex is (T, NK, myeloid, Rest). The sequential binary partition is pre-specified:

1. T/NK/myeloid versus Rest. Scale sqrt(3/4).
2. T/NK versus myeloid. Scale sqrt(2/3).
3. T versus NK. Scale sqrt(1/2).

These three coordinates are an orthonormal basis. Balances 2 and 3 do not involve Rest, so they are identical to the two balances of the renormalized (T, NK, myeloid) simplex. The within-immune Dirichlet model is not redundant with them, because its multinomial total is T + NK + myeloid rather than all cells.

Each balance is regressed on cohort (reference GSE131907) plus the within-cohort z-score of malignant CLDN4 percent-positive. The primary p-value is a 1999-draw permutation that shuffles the z-score inside each cohort. HC1 intervals are descriptive. A joint test uses Wilks' lambda on the balance block, with the same permutation scheme.

Spearman associations are pooled across the four cohorts with DerSimonian–Laird on the Fisher z scale. Quartiles are the locked within-cohort labels (Q1 n = 19, Q4 n = 16), not recut.

## Dirichlet-multinomial

Counts are raw (no pseudocount). The mean is a softmax with Rest as the reference part in the 4-part model and myeloid as the reference in the within-immune model. Precision is a single shared parameter. Coefficients are maximum likelihood. The reported permutation p-value shuffles the CLDN4 z-score inside cohort (299 draws). The displayed fraction change is the mean predicted proportion at +0.5 SD minus the mean at −0.5 SD, with cohort indicators left as observed.

Centered log-ratio regressions are the symmetric companion: each part versus the geometric mean of all four parts, same design, 1999 permutations, Benjamini–Hochberg across the four parts.

## Augur-style priority

Cells from Q1 and Q4 units are subsampled (up to 25 per cell type per unit). Features are a shared NSCLC activity panel, log-normalized to library size. CLDN4 is removed from every classifier. Malignant cells also drop EPCAM, KRT7/8/18/19, CLDN3/7, and TACSTD2, so the malignant score is not the CLDN4 program itself.

Inside each training fold of StratifiedGroupKFold (folds grouped by unit), genes are residualized by cohort mean and the 40 highest-variance genes are kept. The classifier is L2-penalized logistic regression (C = 0.2). A random forest was fit first; its out-of-fold scores sat below 0.5 when Q1/Q4 labels were shuffled inside cohort, so it is not the reported scale. That forest table is `results/tables/augur_priority_forest_sensitivity.tsv`. The primary score is the AUC of the unit-mean out-of-fold probability. The cell-level AUC is the classical Augur number and is secondary, because cells from one unit are not independent. The null shuffles Q1/Q4 labels inside cohort, moving every cell of a unit together.

GSE205335 expression comes from the double-gzip GEO RDS. The other three cohorts are read from the GEO count matrices. The 350-cell Harmony object is not the input.
