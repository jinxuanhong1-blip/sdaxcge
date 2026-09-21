# pySCENIC on GSE131907 malignant cells

Observational GRN on author-malignant cells from GSE131907 (Kim et al. 2020, PMID 32385277). The contrast is CLDN4-high versus CLDN4-low. It is not a knockdown.

## Why pySCENIC, not SCENIC+

SCENIC+ needs matched scATAC. GSE131907 on GEO is 10x scRNA-seq only. pySCENIC 0.12.1 is the method that matches this matrix: GRNBoost2, cisTarget, AUCell.

Concordant-4 (GSE123902, GSE131907, GSE205335, GSE189357) was not merged. The host has 15 GiB RAM. Two cisTarget ranking databases plus one malignant matrix already use that budget.

## Cells

Author `Cell_subtype` in {Malignant cells, tS1, tS2, tS3}. Empty subtype labels are not called malignant. Cells with UMI < 200 are dropped. The public raw UMI matrix is the input (the 2.86 GB log2TPM text file is not used).

## GRN

- Genes detected in ≥1% of QC malignant cells form the AUCell ranking universe.
- GRN features: top 2000 variable genes on log1p(CP10k), plus Lambert/Aerts TFs detected in ≥5% of QC malignant cells, plus the comparator program genes.
- Cells: a sample-stratified subsample (target 3000, seed 131907) drawn from samples with ≥20 QC malignant cells. The draw is not restricted to CLDN4-high cells.
- Expression for GRNBoost2 and for the TF–target correlation inside `modules_from_adjacencies`: log1p(CP10k), library size from all genes.
- GRNBoost2 early-stopping window stays at the arboreto default (25). Seed 131907.
- cisTarget: hg38 v10 cluster rankings, 500 bp up / 100 bp down, and ±10 kb. Motif annotation `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`. Rank threshold 1500. NES ≥ 3.0. Activating modules only (pySCENIC default).

`np.object = object` is set before importing pySCENIC, because NumPy 2 removed that alias and pySCENIC 0.12.1 still references it. arboreto 0.1.6 always builds an unused meta dataframe; dask 2026 rejects that empty list, so the call is skipped. pandas 3 drops grouping columns inside `DataFrameGroupBy.apply`, so the pySCENIC “top regulators per target” step is done with `nlargest` per target instead. That selects the same rows. The GradientBoostingRegressor settings are the arboreto GRNBoost2 defaults (`learning_rate=0.01`, `n_estimators=5000`, `max_features=0.1`, `subsample=0.9`, early-stop window 25).

## AUCell contrast

AUCell (pySCENIC, auc threshold 0.05, weights on) is run on every QC malignant cell. Within a cell, ranks of raw counts and of log1p(CP10k) are the same.

Primary signatures:

- each cisTarget regulon after CLDN4 is removed (CLDN4 is the split gene)
- the same regulon as returned, kept only to show whether CLDN4 membership changes the sign
- program sets that are **not** cisTarget regulons: MSigDB Hallmark IFN-α, IFN-γ, their union, a fixed MHC-I/APM list, a fixed tight-junction list, a keratin list, and a ribosome control. CLDN4 is removed from each.

High versus low: within each sample, median of malignant CLDN4 log1p(CP10k). If that median is 0, high = CLDN4 > 0. A sample is paired only when both tails have ≥20 cells. The patient value is the cell-weighted mean of AUCell in that patient's paired samples. The test is a two-sided Wilcoxon signed-rank on the patient deltas. Benjamini–Hochberg FDR is computed inside the cisTarget family and inside the program family separately. A tLung-only rerun recomputes the split on primary-tumor cells. A between-patient Spearman of CLDN4+ fraction versus mean AUCell is a different question and is labeled as such.

## Sign labels (not fit on this matrix)

- KD-like expectation, carried from earlier AUCell work: IFN/MHC Δ < 0 (lower in CLDN4-high); TJ/keratin Δ > 0.
- GSE207704 public CLDN4 CRISPR: IFN down after loss, which predicts observational Δ > 0.
- GSE50927 public whole-lung Cldn4 KO, n=1: IFN/MHC up after loss, which predicts observational Δ < 0.

SAME / OPPOSITE is the comparison of the observational median Δ to that reference. Program AUCell is not described as a cisTarget regulon.

The primary IFN board is STAT1, STAT2, and IRF1–IRF9 cisTarget regulons, plus Hallmark IFN-α, IFN-γ, and their union scored by program AUCell. A cisTarget regulon that only overlaps Hallmark IFN genes (hypergeometric k ≥ 5 and FDR < 0.05) is written in a separate table. That overlap does not make the transcription factor an IFN regulator, and those rows are not added to the IFN SAME / OPPOSITE counts. Their observational signs are still printed.

## What is not claimed

Motif recovery is not ChIP. Patient n is not cell n. This atlas has no ICI label. ELF3–CLDN4 is not re-discovered here.
