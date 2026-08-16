# results/w200/GEO_2010_2014

Leftover GEO window **2010-01-01 .. 2014-12-31**, lung × ICI / immunotherapy,
TACSTD2 / CLDN4.

## Honest verdict

**No leftover human lung PD-1/PD-L1/CTLA-4 ICI treatment-response cohort.**
See `verdict.json` and `notes/geo_2010_2014/WRITEUP.md`.

**No invented IDs.** The 9 accessions in `search_candidates.tsv` are live
NCBI `gds` hits, re-verified in SOFT (`verified_series.tsv`).

## Key files

| file | what |
|---|---|
| `search_manifest.json` | exact queries + UID counts (4 / 9 / 0) |
| `search_candidates.tsv` | 9 real GSE |
| `verified_series.tsv` | SOFT + leftover_reason |
| `leftover_shortlist.tsv` | same 9, classified |
| `download_manifest.tsv` | 6 processed matrices downloaded |
| `clinical/` | per-sample characteristics |
| `analysis/gene_availability_audit.tsv` | TACSTD2/CLDN4/CD274/PDCD1 probes |
| `analysis/analysis_results.json` | summaries, ρ, Mann-Whitney |
| `analysis/GSE35640_tacstd2_cldn4.tsv` | per-sample leftover vaccine series |
| `analysis/GSE27556-GPL570_tacstd2_cldn4.tsv` | per-sample lung epitope series |
| `master_summary.tsv` / `catalog.tsv` | joined table |
| `verdict.json` | machine-readable verdict |

`data/` is git-ignored; re-run `03_download.py`.
