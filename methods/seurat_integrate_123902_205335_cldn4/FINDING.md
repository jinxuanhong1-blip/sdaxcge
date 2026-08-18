# FINDING — Seurat / Harmony pair GSE123902 + GSE205335, CLDN4-only

ADDITIVE. **CLDN4 only.** Thesis already correct. Public human GSE123902 +
GSE205335 only (the pair that previously differed). No GSE148071. No
GSE127465. No dual-high. No Python-only primary. Concordant-4 four-way
merge is a different agent.

Primary engine: **R + Seurat + Harmony**. Honest unit = **patient**.

_Pending `01_seurat_harmony_integrate.R`. Do not quote numbers from this stub._

See `METHODS.md`. Reproduce:

```bash
bash methods/seurat_integrate_123902_205335_cldn4/scripts/00_download.sh
Rscript methods/seurat_integrate_123902_205335_cldn4/scripts/01_seurat_harmony_integrate.R
```
