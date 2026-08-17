1. `00_download.py` — GSE131907 + GSE205335 processed files and NicheNet-v2 priors. Not GSE207422.
2. `01_convert_prior.py` — RDS → TSV/parquet (no R).
3. `02_extract.py` — patient-level gene summaries; CLDN4-only split.
4. `03_analyze.py` — MultiNicheNet-style activity, paired DE, figures, FINDING.md.
