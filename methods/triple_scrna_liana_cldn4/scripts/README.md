# Scripts

1. `00_download.py` — GEO processed files for GSE131907, GSE148071, GSE205335 + CellPhoneDB v5.
2. `01_build_pairs.py` — parse CellPhoneDB interactors; write `resources/cellphonedb_v5_lr_pairs.tsv`.
3. `02_extract.py` — stream lineage + LR genes; keep malignant and T/NK cells.
4. `03_run_ccc.py` — global CLDN4 median, patient-level paired scores, figures, FINDING.md.
