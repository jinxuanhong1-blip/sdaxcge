# w200 GEO 2019 leftover (TACSTD2 / CLDN4 × lung ICI)

Calendar-2019 GEO slice (`PDAT` 2019-01-01 … 2019-12-31). Skip series already
analyzed in `cursor/fable-geo-2019-2021-c71e`. Do not invent accessions or
clinical labels.

```
python3 scripts/w200_geo_2019/01_search_geo.py
python3 scripts/w200_geo_2019/02_fetch_metadata.py
python3 scripts/w200_geo_2019/03_triage_and_soft.py
python3 scripts/w200_geo_2019/04_build_triage.py
```

Outputs live under `results/w200/GEO_2019/` and `notes/w200_geo_2019/`.
Raw series-matrix downloads are local-only (gitignored).
