# GEO 2017–2018 leftover slice (TACSTD2 / CLDN4 × lung ICI)

Run from the repo root:

```
python3 scripts/w200/geo_2017_2018/01_search_geo.py
python3 scripts/w200/geo_2017_2018/02_build_tables.py
```

`01` hits NCBI E-utilities (`gds`, GSE, *Homo sapiens*, PDAT 2017-01-01..2018-12-31).
`02` does not re-download; it reads the matrices already under
`results/w200/GEO_2017_2018/downloads/` and writes the triage tables.

No association test is run. The two real lung ICI series in this window
(GSE93157, GSE110390) do not deposit TACSTD2 or CLDN4.
