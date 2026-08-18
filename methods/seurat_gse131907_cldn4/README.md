# Seurat GSE131907 — CLDN4-only, patient unit

ADDITIVE. Thesis already correct. No dual-high. No GSE148071.

Public processed UMI (Kim 2020, GSE131907) → `CreateSeuratObject` → patient-level
malignant CLDN4 vs T/NK, and malignant IFN/MHC/TJ Q4 vs Q1.

```bash
bash methods/seurat_gse131907_cldn4/scripts/download.sh /tmp/gse131907
python3 methods/seurat_gse131907_cldn4/scripts/extract_gene_matrix.py \
  --datadir /tmp/gse131907 --outdir /tmp/gse131907/subset
Rscript methods/seurat_gse131907_cldn4/scripts/analyze_seurat.R \
  /tmp/gse131907/subset /tmp/gse131907 \
  methods/seurat_gse131907_cldn4/results \
  methods/seurat_gse131907_cldn4/FINDING.md \
  methods/seurat_gse131907_cldn4/data/families.json
```

If Seurat cannot install, stop. Do not fall back to a Python-only primary.
