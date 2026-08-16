# Methods (this folder)

1. Confirmed GEO scRNA metadata is 15 samples, not 92,330 cells.
2. Downloaded BMC Additional files 1 and 3; confirmed they are clinical / gene-list tables.
3. Listed `Junjie-Hu/NSCLC-immunotherapy` via GitHub API (scripts only; no clone).
4. Downloaded DRMref Seurat RDS files (not committed). Extracted `meta.data` with R `SeuratObject` 5.4.0.
5. Streamed GEO UMI matrix for the `TACSTD2` row; joined on barcode (30,877 / 30,877 match).
6. Per post-treatment sample: mean `log1p(1e4 * TACSTD2_UMI / nCount_RNA)` in cells labeled `Malignant cells`.
7. MPR includes pCR (P06). Wilcoxon exact two-sided by enumerating all ways to assign 4 of 12 samples to MPR. Spearman uses midranks and a Student-t p-value on n−2 df.
8. T/NK fraction = (CD8+ T + CD4+ T + NK) / annotated cells in that sample.

Primary pre-specified tests: Wilcoxon on mean log1p(CP10k) and Spearman of that mean vs T/NK fraction. Other metrics are secondary.
