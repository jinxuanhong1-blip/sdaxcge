# Seurat concordant-4 CLDN4-only (R primary)

ADDITIVE. **CLDN4-only.** Honest unit = patient / donor / sample.
Datasets ONLY: GSE123902 + GSE131907 + GSE205335 + GSE189357.
Not GSE148071, GSE127465, GSE207422, GSE154826. No dual-high.

Primary is **R + Seurat 5 + Harmony** (`harmony::RunHarmony` on the Seurat object).
Python is used only to stream the GSE131907 genes×cells text into a sparse
atlas subset (the full dense table does not fit in 16 GB).

```bash
Rscript methods/seurat_concordant4_cldn4/analyze.R
```

Writes `FINDING.md`, DimPlots, the patient/donor/sample table, and
`results/objects/seurat_harmony_integrated.rds`.
