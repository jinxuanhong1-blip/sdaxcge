# GSE154977 KP 30w 10x — Cldn4-only (Seurat)

ADDITIVE public mouse. Primary analysis is **R + Seurat** (not Python).

```bash
/home/ubuntu/micromamba/envs/seurat/bin/Rscript methods/seurat_gse154977_cldn4/analyze.R
```

Creates `objects/gse154977_kp30w_seurat.rds`, mouse-level tables, DimPlot/VlnPlot, and `FINDING.md` inputs under `tables/` / `figures/`.

Honest n = biological mice (4), not 11,017 cells.
