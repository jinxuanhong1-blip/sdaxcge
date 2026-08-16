# w200 / GEO_2020 — leftover 2020 human lung ICI × TACSTD2 / CLDN4

Run in order (Python 3; pandas / numpy / scipy / matplotlib):

```
python3 scripts/w200/GEO_2020/00_prior_coverage.py
python3 scripts/w200/GEO_2020/01_search_geo.py
python3 scripts/w200/GEO_2020/02_fetch_metadata.py
python3 scripts/w200/GEO_2020/03_triage.py
python3 scripts/w200/GEO_2020/04_probe_leftovers.py
python3 scripts/w200/GEO_2020/05_inspect_matrices.py
python3 scripts/w200/GEO_2020/06_extract_and_test.py
```

`00_prior_coverage.py` expects `/tmp/inv/branch_gse.tsv` and `/tmp/inv/path_gse.tsv`
(git-grep of every remote branch). If those files are absent the coverage table
is empty and triage still runs.

Do not commit `downloads/platforms/` (GPL SOFT files are huge). The committed
two-gene extract for GSE141479 is
`results/w200/GEO_2020/downloads/GSE141479/GSE141479_TACSTD2_CLDN4_rows.txt`.
The full series matrix is
`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141479/matrix/GSE141479_series_matrix.txt.gz`.
