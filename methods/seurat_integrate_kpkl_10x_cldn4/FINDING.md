# FINDING — Seurat/Harmony integrate KP+K+KL 10x (Cldn4-only)

**ADDITIVE. Public mouse. Cldn4-only. No dual-high.** Thesis is already correct and is not rewritten.

This stub is replaced by live numbers after:

```bash
bash methods/seurat_integrate_kpkl_10x_cldn4/scripts/install_r.sh
bash methods/seurat_integrate_kpkl_10x_cldn4/scripts/download.sh /tmp/kpkl_10x
Rscript methods/seurat_integrate_kpkl_10x_cldn4/scripts/analyze.R \
  --data /tmp/kpkl_10x \
  --out methods/seurat_integrate_kpkl_10x_cldn4
```

Requested trio (keep only if a processed matrix exists): GSE154977 (KP) + GSE180963 (K/KL) + GSE165641 (KL). Drop any series with no matrix. Do not add GSE179502, GSE154989, GSE267321, GSE127465, GSE50927, human, or private 8 KL.

Honest unit = **mouse**. Dataset is a covariate. Within-genotype and leave-one-dataset-out must be reported so a genotype mix is not sold as a Cldn4 effect.
