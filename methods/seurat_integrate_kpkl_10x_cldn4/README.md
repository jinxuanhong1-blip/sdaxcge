# seurat_integrate_kpkl_10x_cldn4

ADDITIVE public mouse. **R + Seurat / Harmony** integration of processed 10x GEMM lung series:

- GSE154977 (KP 30w) — include only if processed matrix exists
- GSE180963 (K / KL) — include only if processed matrix exists
- GSE165641 (KL) — include only if processed matrix exists

Cldn4-only. Honest unit = **mouse**. Dataset is a covariate. Within-genotype and leave-one-dataset-out are reported so a genotype mix is not sold as a Cldn4 effect.

Do **not** add GSE179502, GSE154989, GSE267321, GSE127465, GSE50927, human, or private 8 KL. No dual-high. No Python-only primary. Drop any series with no matrix.

```bash
bash methods/seurat_integrate_kpkl_10x_cldn4/scripts/install_r.sh
bash methods/seurat_integrate_kpkl_10x_cldn4/scripts/download.sh /tmp/kpkl_10x
Rscript methods/seurat_integrate_kpkl_10x_cldn4/scripts/analyze.R \
  --data /tmp/kpkl_10x \
  --out methods/seurat_integrate_kpkl_10x_cldn4
```

See `FINDING.md`.
