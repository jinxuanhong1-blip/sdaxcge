# Concordant-4 cell-level atlas (methods / identity)

ADDITIVE. **CLDN4-only.** How the four public scRNA sets were taken in,
integrated, clustered, and annotated.

Datasets: **GSE123902 + GSE131907 + GSE205335 + GSE189357**.
Not GSE148071 / GSE127465 / GSE154826 / GSE207422.

This figure is methods/identity. It does **not** re-run or re-audit the
T/NK ρ or IFN/MHC numbers in PR #503 (`n=65`).

```bash
python3 methods/concordant4_atlas_umap_annotate/download.py
python3 methods/concordant4_atlas_umap_annotate/analyze.py
```
