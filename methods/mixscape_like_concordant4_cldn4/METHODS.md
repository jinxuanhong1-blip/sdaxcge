# Methods

## Question

In the locked concordant-4 cohorts, does a Mixscape-style local residual
that treats CLDN4-low malignant cells as KD-like resemble Hallmark
interferon, and does that residual change across CLDN4 quartiles?

## Not CRISPR

Papalexi et al. and Seurat `CalcPerturbSig` subtract, for each cell that
carries a targeting guide, the average expression of its k nearest
non-targeting cells. `RunMixscape` then fits a mixture to separate
perturbed cells from escapees.

This analysis copies only the local subtraction. CLDN4 expression
quartiles stand in for guide classes. There is no guide identity, no
non-targeting gRNA, and no escapee posterior. Cells are not knockouts.

## Cohorts and cells

Same 65 units as the concordant-4 Seurat integration:

| cohort | unit | malignant definition |
|---|---|---|
| GSE123902 | donor | one tumor sample per donor (PRIMARY preferred only when it sorts after METASTASIS by the locked rule: tissue order METASTASIS then PRIMARY, first row kept); marker malignant |
| GSE131907 | sample | author `Cell_subtype == Malignant cells`; the 21 tumor-bearing samples already in the locked table |
| GSE205335 | patient | author `lineage.sub == Malignant cells`; Normal tissue dropped |
| GSE189357 | patient | marker malignant |

Marker malignant: (EPCAM > 0 or KRT8 > 0 or KRT18 > 0 or KRT19 > 0) and PTPRC == 0, on raw counts.

QC, applied after the malignant gate: ≥200 genes, ≥500 UMIs, mitochondrial percent < 20. Mitochondrial genes are those with symbol prefix `MT-`.

## Quartiles

Within a unit, on CLDN4 UMI counts among QC-pass malignant cells:

- Q1: count ≤ empirical 25th percentile
- Q2: above the 25th and ≤ the 50th
- Q3: above the 50th and ≤ the 75th
- Q4: above the 75th

Ties stay in the lower bin. If most cells are CLDN4-negative, Q2 and sometimes Q3 are empty. Those bins are left empty. Cells are not randomly assigned across a tie.

Q4 is the NT-like pool. The lowest occupied bin with at least 20 cells (Q1 if it qualifies) is KD-like.

A unit is skipped when Q4 has fewer than 15 cells, when k would fall below 10, or when CLDN4 does not vary.

## Local signature

Counts are log1p-normalized to 10,000 UMIs per cell (Seurat `data` slot analog).

PCA uses the top 2,000 variable genes by variance of that log-normalized matrix, after dropping:

- CLDN4
- Hallmark interferon alpha and interferon gamma
- the custom MHC-I / antigen-processing list used in the concordant-4 Seurat run
- Hallmark spermatogenesis (IFN genes removed)
- detection-matched control genes chosen for those readouts

Genes must be detected in at least 5% of the unit’s QC malignant cells to enter the variable-gene list. Features are z-scored across cells and clipped to [−10, 10]. The first 15 PCs are used (fewer if the unit is smaller).

For each malignant cell, k = 20 nearest Q4 cells in that PCA (Euclidean). If Q4 has fewer than 21 cells, k = n_Q4 − 1. A cell is not its own neighbor.

The residual of a gene set is the mean, across genes in the set, of log-normalized expression minus the mean of the k neighbors. That is the Mixscape perturbation signature restricted to the set.

## Scores

- **IFN residual:** Hallmark IFNα ∪ IFNγ (symbols from the concordant-4 gene-set file; CLDN4 is not in the set).
- **IFN alpha and IFN gamma:** the two hallmarks separately.
- **MHC-I/APM:** the custom list with genes that also sit in the IFN union removed, so it is not a second copy of IFN.
- **Matched control:** for each IFN gene present, the non-readout gene with the closest detection rate, without replacement. The competitive delta is IFN residual minus matched residual.
- **Spermatogenesis:** same construction, as a hallmark that is not the IFN claim.
- **Unmatched IFN:** mean log-normalized IFN without neighbor subtraction. Reported so the local residual can be compared with a plain Q1 versus Q4 difference. It is not the Mixscape estimand.

The KD-like gene signature used for cosine similarity is the mean residual, in KD-like cells, of the variable genes plus the IFN genes. Cosine with the IFN indicator is compared with 499 permutations of gene labels inside the unit.

## Tests

The unit is the independent observation.

Primary dose test: Spearman of CLDN4 count versus IFN residual among non-Q4 cells in the unit (Q4 is the self-neighborhood null and is left out). Units are equally weighted. Inference is a one-sample t-test on Fisher z and a Wilcoxon signed-rank test across units. A DerSimonian–Laird meta-analysis combines the four cohort means, with weights from the empirical standard error of each cohort’s Fisher z.

Quartile dose: equal-unit mean residual by quartile. Paired Wilcoxon tests compare unit-level Q1 versus Q3 (shared NT-like pool) and Q1 versus Q4 (KD-like versus the self-null).

Similarity: Wilcoxon signed-rank of the KD-like IFN − matched delta against zero, and the same test for spermatogenesis.

p-values are descriptive.

## Software

Python 3 with numpy, scipy, scikit-learn, matplotlib, and pandas (GSE123902 CSV only). GSE205335 is read in R (`Matrix`) from the author RDS and written out as Matrix Market, because that matrix is distributed as an RDS. Patient and tissue labels for that cohort come from `data/GSE205335_sample_map.tsv`, parsed from the public GEO family soft. The signature math is the Python implementation above, not an R `CalcPerturbSig` call. A built-in simulation (CLDN4-low cells carry IFN counts; spermatogenesis does not) checks the sign before any cohort is scored.
