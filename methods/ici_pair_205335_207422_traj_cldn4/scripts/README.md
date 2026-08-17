# Scripts

1. `download.py` — GEO processed UMI + identity/metadata only (<2 GB).
2. `extract.py` — per-dataset epithelium h5ad. No joint object.
3. `analyze.py` — PAGA + DPT per dataset, then stacked patient tables.
4. `gene_sets.py` — locked AT2 / barrier (CLDN4 out) / A3 leftover genes.

Root is leftover AT2-like epithelium. Never CLDN4-high. Harmony is not called.
