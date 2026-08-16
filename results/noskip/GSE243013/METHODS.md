# GSE243013 TACSTD2 / CLDN4 vs MPR

## Why this dataset was not skipped

GEO series [GSE243013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243013) (Liu et al., *Cell* 2025, PMID 40147443) has a public processed count matrix:

- `GSE243013_NSCLC_immune_scRNA_counts.mtx.gz` (Matrix Market, gzipped, ~6.6 GiB)
- `GSE243013_genes.csv.gz`
- `GSE243013_barcodes.csv.gz`
- `GSE243013_NSCLC_immune_scRNA_metadata.csv.gz` (includes `pathological_response`)

The GEO series matrix table itself is empty (`!Sample_data_row_count` is 0 for every GSM). The usable processed expression object is the supplementary MTX, not the series matrix.

`GSE243013_RAW.tar` is listed as scTCR-seq / sample TCR archives, not a gene-expression matrix. `GSE243013_NMF_all_group_5.csv.gz` is TIME-subtype labels only.

## Assay actually measured

Post-neoadjuvant chemo-immunotherapy surgical tumors, 10x 5' scRNA of **CD45+ immune cells**. Major types in the metadata are T/NK, B, and myeloid only. This is **not** a tumor-epithelial or bulk tumor matrix.

TACSTD2 and CLDN4 are epithelial genes. Any signal here is immune-compartment expression or ambient/misassigned transcripts, not cancer-cell TACSTD2/CLDN4. Results must be read with that limit.

## Matrix geometry (verified from file header + gene list)

```
%%MatrixMarket matrix coordinate real general
1254749 31831 2010550708
```

Rows = cells (same order as `barcodes.csv` / metadata `cellID`; verified identical).  
Columns = genes (order of `genes.csv`).  
TACSTD2 = column 1057. CLDN4 = column 11797. Each symbol occurs once.

## Response labels

Metadata field `pathological_response` values observed: `pCR`, `MPR`, `non-MPR`, and one `unknowm` (typo in the deposited file).

Paper definition (from the *Cell* text): residual viable tumor ≤10% = MPR; pCR is a subset of MPR. GEO splits pCR out as its own label. Primary contrast here:

- **MPR-any** = `MPR` or `pCR`
- **non-MPR** = `non-MPR`
- patient `P433` (`unknowm`) excluded from tests

GEO metadata has more `sampleID`s than the paper's n=234. Counts are taken from the deposited metadata, not forced to 234.

`pathological_response_rate` is mostly numeric 1−RVT, with a few non-numeric strings (`>90%`, `<40%`). Labels were used for grouping; rates were not re-derived except as a consistency check on numeric rows (no MPR/non-MPR mismatches among numeric rates).

## Scores per patient

From the two extracted MTX columns, joined to metadata `sampleID` and `total_counts`:

- detection fraction: cells with count > 0 / n cells
- mean raw count
- mean CPM: mean of (count / cell `total_counts` × 1e6)

## Statistics

Two-sided Mann–Whitney U, normal approximation with tie correction and continuity correction. Effect size: Cliff's delta and AUC = 0.5 + 0.5 × delta. No multiple-testing claim beyond reporting all pre-specified gene × metric × contrast tests.

Script: `scripts/analyze_gse243013_tacstd2_cldn4_mpr.py`.
