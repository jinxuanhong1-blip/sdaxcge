# GSE189357 — REAL Slingshot/PAGA, CLDN4 only

ADDITIVE public slice. Tumor epithelium from Zhu et al. 2022 (GSE189357).
**CLDN4 only.** No dual-high. Patient is the unit. **n may be 9 — say so.**

1. `scripts/install_tools.sh` — scanpy + Bioconductor **slingshot** (required).
2. `scripts/download.py` — GEO processed 10x tar only.
3. `scripts/extract_epithelium.py` — marker-malignant tumor epithelium → h5ad.
4. `scripts/analyze.py` — Harmony + Leiden + PAGA + REAL Slingshot; writes the lineage table and FINDING.md.

Done when `results/tables/slingshot_lineages.tsv` exists.
