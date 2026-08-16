# Downloaded GEO processed data (not tracked in git)

This directory holds the **open processed supplementary + series-matrix files**
downloaded for the 23 in-scope series (human + lung + immune-checkpoint), one
sub-directory per accession (`GSE#####/`).

- Total on disk at generation time: **~3.2 GB across 79 files** (0 failures).
- Only files **< 2 GB each** were stored; `GSE72094_RAW.tar` (2.1 GB) was skipped
  per the "< 2 GB" rule (its expression is still available via the series matrix).
- The exact file list, listed vs stored byte sizes, status, and source URLs are in
  [`../download_manifest.tsv`](../download_manifest.tsv).

## Reproduce

```bash
python3 scripts/fable_geo_2015_2018/01_search.py     # GEO search -> candidates
python3 scripts/fable_geo_2015_2018/02_verify.py     # verify every accession
python3 scripts/fable_geo_2015_2018/03_download.py   # re-download into this dir
```

Nothing here is fabricated: every file comes from
`https://ftp.ncbi.nlm.nih.gov/geo/series/...` for a real, verified accession.
