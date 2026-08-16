# scripts/fable_tacstd2_kdko

Pipeline for the TROP2/TACSTD2 loss-of-function meta-analysis. Run in order.

| Script | Purpose |
|---|---|
| `01_search_geo.py` | Exhaustive GEO DataSets search → `results/.../geo_search_*` |
| `02_download.sh` | Download processed matrices / DESeq2 tables (<2 GB) → `results/.../data/` |
| `gene_panels.py` | CLDN4 / junction / IFN-immune gene panels (shared) |
| `stats_utils.py` | Welch t-test, BH-FDR, Mann-Whitney, Cohen's d, competitive & self-contained set tests |
| `analysis_common.py` | Shared engine: mygene mapping, gene-level collapse, per-gene + set-level tests, output writing |
| `03_analyze_gse334497.py` | Mouse 4T1 Trop2 KO vs WT (normalized counts) |
| `04_analyze_gse289287.py` | Human T-47D Trop-2 KO + DSG2 KO (author DESeq2) |
| `05_analyze_gse245459.py` | Human ovarian shTACSTD2 vs shNC (FPKM) |
| `06_analyze_gse15212.py` | Human colorectal siTACSTD2 vs neg-ctrl (Agilent array) |
| `07_synthesize.py` | Cross-dataset `SYNTHESIS_*` tables |

Requires network access to NCBI E-utilities/FTP and mygene.info, and Python packages
`pandas numpy scipy statsmodels requests` (GEOparse optional).
