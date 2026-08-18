# Seurat / Harmony pair GSE123902 + GSE205335 (CLDN4-only)

ADDITIVE. **CLDN4 only.** Public human GSE123902 + GSE205335 only — the pair
that previously differed. Seurat v5 + Harmony. Honest unit = patient.
Primary: CLDN4 vs T/NK; malignant IFN/MHC Q4 vs Q1.

No GSE148071. No GSE127465. No dual-high TACSTD2∩CLDN4. No Python-only
primary. Concordant-4 four-way merge is a different agent.

See `FINDING.md`.

```bash
bash methods/seurat_integrate_123902_205335_cldn4/scripts/00_download.sh
Rscript methods/seurat_integrate_123902_205335_cldn4/scripts/01_seurat_harmony_integrate.R
```
