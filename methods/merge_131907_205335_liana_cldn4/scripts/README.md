# Scripts

```bash
python3 methods/merge_131907_205335_liana_cldn4/scripts/00_download.py
python3 methods/merge_131907_205335_liana_cldn4/scripts/01_analyze.py
```

`00_download.py` writes processed GEO files under `/tmp/gse131907` and
`/tmp/gse205335`. `01_analyze.py` scores CellPhoneDB-style pairs and writes
`results/lr_table.tsv` plus `FINDING.md`.
