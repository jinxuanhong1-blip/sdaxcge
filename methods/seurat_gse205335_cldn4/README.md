# Seurat GSE205335 — CLDN4-only malignant vs T/NK + IFN/MHC/TJ

ADDITIVE. Public processed GSE205335 only. Patient is the unit.
No dual-high. No GSE148071. Primary analysis is R + Seurat (`CreateSeuratObject`).

```bash
Rscript methods/seurat_gse205335_cldn4/scripts/download.R /tmp/gse205335
GSE205335_DATA=/tmp/gse205335 SEURAT_GSE205335_ROOT=methods/seurat_gse205335_cldn4 \
  Rscript methods/seurat_gse205335_cldn4/scripts/analyze.R
```

See `FINDING.md` for the patient-level result.
