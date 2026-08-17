# Scripts

1. `download.py` — public GEO UMI + annotation for GSE131907 and GSE205335 only.
2. `extract_malignant.py` — author-malignant cells, cap ≤200 / unit, joint h5ad.
3. `analyze.py` — CytoTRACE2 (or documented Gulati 2020 CytoTRACE), patient-level table, figures, FINDING.md.

`run_all.sh` creates a numpy<2 venv (required by `cytotrace2-py`) and runs the three steps.
