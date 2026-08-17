# Scripts

1. `install_tools.sh` — Python scanpy stack + R/Bioconductor slingshot.
2. `download.py` — public GEO processed files only (no FASTQ, no GSE189487).
3. `extract_epithelium.py` — marker-malignant tumor epithelium (TD1–TD9) → h5ad.
4. `run_slingshot.R` — REAL Slingshot; start cluster is never CLDN4-high.
5. `analyze.py` — Harmony + Leiden + PAGA + Slingshot; writes lineage table and FINDING.md.
6. `gene_sets.py` / `ifn_genes.json` — locked AT2 / barrier (CLDN4 excluded) / Hallmark IFN.
