# seurat_integrate_165641_180963_cldn4

Integrate public [GSE165641](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165641) + [GSE180963](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180963) KL/K GEMM 10x. **R + Seurat** (`CreateSeuratObject` per dataset) then **Harmony** via `IntegrateLayers`. Cldn4-only. Honest n = mice after merge. See `FINDING.md` after `analyze.R` runs.

```bash
bash methods/seurat_integrate_165641_180963_cldn4/scripts/install_r.sh
bash methods/seurat_integrate_165641_180963_cldn4/scripts/download.sh /tmp/geo/work
Rscript methods/seurat_integrate_165641_180963_cldn4/scripts/analyze.R --data /tmp/geo/work --out methods/seurat_integrate_165641_180963_cldn4
```
