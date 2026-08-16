# Pipeline — GEO 2019–2021 human lung ICI × TACSTD2/CLDN4

Run in order (Python 3; `pip install pandas numpy scipy lifelines statsmodels matplotlib`):

| step | script | output |
|---|---|---|
| 1 | `01_search_geo.py` | `results/.../search_uids.json` (90 real GSE UIDs) |
| 2 | `02_fetch_metadata.py` | `results/.../candidates_metadata.json` |
| 3 | `03_download.py` | `results/.../downloads/<GSE>/` (processed matrices + series_matrix; 2 GB/file cap) |
| 4 | `04_build_clinical.py` | `results/.../clinical/<GSE>_clinical.csv` |
| 5 | `05_analyze.py` | `results/.../tables/*`, `results/.../figures/*` |
| 6 | `06_triage.py` | `results/.../tables/triage_all_candidates.csv` |

All GEO accessions are real NCBI identifiers. No IDs were invented. Outputs are
confined to `notes/fable_geo_2019_2021/`, `scripts/fable_geo_2019_2021/`,
`results/fable_geo_2019_2021/`. See `notes/.../WRITEUP.md` (zh+en) and
`notes/.../verification_log.md`.
