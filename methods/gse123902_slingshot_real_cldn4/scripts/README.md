# Scripts

1. `install_tools.sh` — scanpy/PAGA plus R/Bioconductor slingshot when apt/R allow.
2. `download.py` — GEO `GSE123902_RAW.tar` (90.4 MB). Skips the 36.5 GB author H5.
3. `extract_epithelium.py` — marker epithelium from dense UMI CSVs → h5ad.
4. `analyze.py` — Leiden + PAGA + real Slingshot; writes `lineage.tsv` and FINDING.md.
5. `run_slingshot.R` — Bioconductor slingshot on PCA + Leiden if R is present.
6. `gene_sets.py` — locked AT2 / barrier (CLDN4 excluded) / IFN / comparator sets.
