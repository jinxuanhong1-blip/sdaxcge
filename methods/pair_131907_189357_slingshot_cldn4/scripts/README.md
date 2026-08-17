# Scripts

1. `download.py` — public GEO processed files only (GSE131907 UMI + annotation; GSE189357 RAW.tar). No GSE148071. No FASTQ.
2. `extract_epithelium.py` — GSE131907 author epithelium; GSE189357 marker epithelium. Cap ≤350/unit. Protect nLung AT2.
3. `run_slingshot.R` — real Bioconductor Slingshot. Start cluster is AT2-rooted, never CLDN4-high.
4. `analyze.py` — Harmony + Leiden + PAGA + Slingshot. Writes the lineage table.
5. `run_all.py` — the three Python steps in order.

Requires `Rscript` and `slingshot` in `R_LIBS_USER` (default `~/R/library`).
