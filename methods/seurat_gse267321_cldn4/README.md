# Seurat GSE267321 LKR13 K / KK / KLK (Cldn4-only)

Additive public-mouse slice. Primary engine is **R + Seurat** (`CreateSeuratObject`) on the GEO normalized CSV. See [FINDING.md](FINDING.md).

```bash
Rscript methods/seurat_gse267321_cldn4/scripts/analyze_gse267321_seurat.R
```

Downloads `GSE267321_Normalized_expression_matrix_02122026.csv.gz` into `/tmp/gse267321/` (not committed). Writes `tables/genotype_table.tsv` and regenerates `FINDING.md`.

No Python-only primary. No private KL object.
