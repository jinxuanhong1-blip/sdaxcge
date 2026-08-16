# GEO 2010–2014 leftover lung ICI / TACSTD2 / CLDN4

Parallel leftover slice. All outputs live under `results/w200/GEO_2010_2014/`
and `notes/geo_2010_2014/`. **No invented IDs**: every accession is a live
NCBI E-utilities hit, re-verified against its GEO SOFT record.

## Reproduce

From the repo root (public NCBI endpoints, no API key):

```bash
python3 scripts/geo_2010_2014/01_search.py
python3 scripts/geo_2010_2014/02_verify.py
python3 scripts/geo_2010_2014/03_download.py
python3 scripts/geo_2010_2014/04_build_clinical.py
python3 scripts/geo_2010_2014/05_target_gene_analysis.py
python3 scripts/geo_2010_2014/06_curate.py
```

`data/` is git-ignored (re-downloadable). Processed files ≥ 2 GB are skipped
and recorded in `download_manifest.tsv`.
