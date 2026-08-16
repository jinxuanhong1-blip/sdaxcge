# scripts

| Script | Role |
|---|---|
| `00_download.py` | GEO files under the 2 GB budget; catalogs the 9.3 GB RDS |
| `01_extract.py` | stream marker panel from UMI txt / MTX (never loads the full matrix) |
| `02_build_composition.py` | patient-level counts + malignant TACSTD2 + median split |
| `03_fit_grid.py` | cohort × annotation × contrast DM-GLM / ALR / naive |
| `04_figures.py` | heatmap, honest-n bars, descriptive scatter |
| `selftest_dm.py` | planted TNK-down on synthetic DM data |
| `run_all.py` | 00→04 |
