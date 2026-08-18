# Seurat GSE207422 — Cldn4-only malignant vs T/NK and IFN/MHC

Additive public slice. Patient is the unit. Not dual-high. Not merged with GSE148071.

See `FINDING.md` for the result.

```bash
python3 methods/seurat_gse207422_cldn4/scripts/download.py
python3 methods/seurat_gse207422_cldn4/scripts/00_write_metadata.py
python3 methods/seurat_gse207422_cldn4/scripts/01_umi_to_10x.py
Rscript methods/seurat_gse207422_cldn4/scripts/02_analyze.R \
  data/GSE207422/tenx \
  methods/seurat_gse207422_cldn4
```

Stop if Seurat or the GEO UMI matrix is missing. Do not fall back to a Python-only primary.
