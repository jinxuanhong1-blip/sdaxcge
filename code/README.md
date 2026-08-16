# Claim B5 analysis code

1. `03_extract_cldn4.py` — pull CLDN4 (`ENSG00000189143`) and clinical labels from ORCESTRA/PredictIO per-study TSVs (`data/raw/icb/`, Zenodo 10.5281/zenodo.7199344).
2. `03b_imvigor210_official.py` — replace the broken constant-floor ORCESTRA IMvigor210 CLDN4 vector with official CoreBiologies counts.
3. `04_meta_analysis.py` — per-cohort ORs, DerSimonian–Laird meta, sensitivities, forest plots → `results/claim_B5/`.

`01_prepare_data.py` / `02_qc_expression.py` were an earlier attempt on the merged Zenodo matrix (`7459023`). That merge **drops CLDN4 entirely** and is not used for the reported meta.
