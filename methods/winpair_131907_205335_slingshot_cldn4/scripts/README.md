# Scripts

1. `download.py` — public GEO processed files only (no FASTQ, no GSE207422, no 2.86 GB TPM).
2. `extract_epithelium.py` — author epithelium from both accessions → joint h5ad.
3. `analyze.py` — Harmony + Leiden + PAGA + AT2-rooted DPT; writes the sample-level table and FINDING.md.
4. `gene_sets.py` — locked AT2 / barrier (CLDN4 excluded) / comparator sets.
