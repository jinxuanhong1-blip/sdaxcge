# seurat_gse180963_cldn4

Public [GSE180963](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180963) K vs KL lung GEMM scRNA. **R + Seurat** (`CreateSeuratObject`). Cldn4-only. Honest n = **2 mice**. See `FINDING.md`.

```bash
bash methods/seurat_gse180963_cldn4/scripts/install_r.sh
bash methods/seurat_gse180963_cldn4/scripts/download.sh /tmp/gse180963
Rscript methods/seurat_gse180963_cldn4/scripts/analyze.R --data /tmp/gse180963 --out methods/seurat_gse180963_cldn4
```

