# B5 leftover HNSCC ICI CLDN4

Catalog leftover public HNSCC ICI cohorts and test CLDN4 only where the gene is measured and ICI labels are public. No tuning.

## Run

```bash
python3 results/w200/B5_HNSC_extra/analyze.py
```

Python 3 stdlib only. Writes `summary.json`, `cohort_catalog.tsv`, `marker_statistics.csv`, `survival_statistics.csv`, `per_sample_expression.csv`, `prat_panel_audit.tsv`.

## Leftovers

| GEO | Why leftover | CLDN4? | Labels? | Used |
|---|---|---|---|---|
| GSE93157 Prat HNSCC | prior w200 unused HNSCC slice | no | yes (all NR) | panel audit only |
| GSE159067 Foy | public R/M HNSCC ICI | no | yes | catalog |
| GSE179730 Liu | neoadjuvant nivo OCSCC | yes (sparse) | Table S2 | expression + RFS/OS |
| GSE190575 ALPHA | HNSCC ICI+afatinib | no | mixed combo | catalog |
| GSE212549 NIVACTOR train | R/M HNSCC ICI | yes | no public response | catalog |
| GSE212550 NIVACTOR test | R/M HNSCC ICI | yes | LTS vs STS | expression |
| EGAD50000002506 | ICI HNSCC RNA-seq | unknown | EGA | inaccessible |

## Locks

- Liu primary: Responder+Stable vs Progressor (paper clinical-benefit grouping).
- Liu sensitivity: pathologic Responder vs rest.
- Liu survival: last-observation text; RFS event = AWD or DOD; OS event = DOD; split = CLDN4 detectable vs zero (median is 0).
- NIVACTOR: deposited LTS vs STS only. Train set not analyzed (no labels).
- No threshold scan. No model fitting.

## Honest

`NO_REPRODUCIBLE_CLDN4_ICI_SIGNAL`. See `RESULTS.md`.
