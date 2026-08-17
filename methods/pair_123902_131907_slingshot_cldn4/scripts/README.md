# Scripts

1. `install_r_slingshot.sh` — R + Bioconductor slingshot (required; no DPT-only stop).
2. `download.py` — public GSE123902 RAW.tar + GSE131907 UMI/annotation.
3. `extract_epithelium.py` — joint epithelial AnnData (cap ≤350/unit).
4. `analyze.py` — Harmony, Leiden, PAGA, REAL Slingshot, sample-level tables.
5. `run_slingshot.R` — Street 2018 slingshot on Harmony/PCA + Leiden.
6. `gene_sets.py` — locked CLDN4 / barrier (no CLDN4) / IFN sets.
