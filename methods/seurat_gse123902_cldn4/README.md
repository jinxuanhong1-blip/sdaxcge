# Seurat GSE123902 — CLDN4-only, patient/donor unit

ADDITIVE human LUAD/NSCLC. **CLDN4 only.** Public processed GSE123902 (Laughney et al., *Nat Med* 2020). Primary engine is **R + Seurat**. `CreateSeuratObject` is called on a 10x-style MTX written from the GEO dense UMI CSVs. Honest unit = patient/donor (LX ID). No TACSTD2∩CLDN4 dual-high. No GSE148071 merge.

## Reproduce

```bash
bash methods/seurat_gse123902_cldn4/scripts/download.sh
Rscript methods/seurat_gse123902_cldn4/scripts/run_seurat_gse123902.R
```

Requires R ≥ 4.3 and Seurat ≥ 5. If Seurat cannot be installed, stop (do not fall back to a Python-only primary).
