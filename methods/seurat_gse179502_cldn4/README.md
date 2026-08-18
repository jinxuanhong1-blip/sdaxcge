# GSE179502 Seurat Cldn4-only (Lkb1 restore)

ADDITIVE. **Cldn4 only.** No Tacstd2∩Cldn4 dual-high. No private 8 KL.
Public processed mtx + features + barcodes only.

```
Rscript methods/seurat_gse179502_cldn4/analyze.R
```

Requires R + Seurat. If Seurat cannot be installed, stop — do not score in Python.

This run: Seurat 5.5.1 `CreateSeuratObject` on the public GEO mtx (16,017 → 12,643 QC cells; honest n = 6 mice). See `FINDING.md`.
