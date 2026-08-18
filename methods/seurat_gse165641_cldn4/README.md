# GSE165641 KL GEMM scRNA — Cldn4-only (Seurat)

Public 10x of **two** KrasG12D/+; Lkb1fl/fl (KL) tumor-bearing mice
([GSE165641](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165641);
Wang / Zhong, *Adv Sci* 2021, PMID 34369094).

Primary is **R + Seurat** `CreateSeuratObject` on Cell Ranger filtered MTX.
Honest **n = 2 mice**. No private 8-KL. No Python primary.

```bash
# once: Cell Ranger tars under /tmp/geo_gse165641
Rscript methods/seurat_gse165641_cldn4/analyze_gse165641.R \
  /tmp/geo_gse165641 methods/seurat_gse165641_cldn4
```

See `FINDING.md`.
