# Rerun B2_CCLE_all

```bash
python3 -m pip install pandas numpy scipy matplotlib statsmodels openpyxl
python3 scripts/B2_CCLE_all/00_download.py --workdir /tmp/ccle_data --outdir results/w200/B2_CCLE_all
python3 scripts/B2_CCLE_all/01_analyze.py --workdir /tmp/ccle_data --outdir results/w200/B2_CCLE_all
```

Raw DepMap RNA (~461 MB) is hashed and not retained. Extracted TACSTD2/CLDN columns and gene-means live under `/tmp/ccle_data/extracted/`.
Outputs: `results/w200/B2_CCLE_all/{tables,figures,WRITEUP.md,catalog.tsv,repro/}`.
