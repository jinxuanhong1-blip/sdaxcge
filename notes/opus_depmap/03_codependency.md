# 03 · Co-dependency

A co-dependency scan is only meaningful when the query gene's Chronos profile is reproducible. The order of evidence below is deliberate.

## Technical floor

Lines screened at least twice: **32** (lung lines). Median between-screen correlation of the gene effect value across those lines:

| class | genes | median r |
| --- | --- | --- |
| common essential | 397 | 0.234 |
| non-essential control | 360 | 0.134 |
| target | 2 | 0.386 |

* **TACSTD2**: between-screen r = 0.430 (p = 0.0142, 32 lines).
* **CLDN4**: between-screen r = 0.343 (p = 0.0543, 32 lines).

## Scan statistics

| query | cohort | lines | genes | profile SD | max abs r | permutation null 95th pct | perm p | genes at FDR<5% |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TACSTD2 | pan-cancer | 1178 | 14873 | 0.126 | 0.253 | 0.187 | 0.005 | 3328 |
| TACSTD2 | lung | 126 | 12348 | 0.121 | 0.539 | 0.526 | 0.040 | 2 |
| CLDN4 | pan-cancer | 1178 | 14873 | 0.127 | 0.400 | 0.186 | 0.005 | 5630 |
| CLDN4 | lung | 126 | 12348 | 0.132 | 0.506 | 0.526 | 0.095 | 165 |
| EGFR | pan-cancer | 1178 | 14873 | 0.335 | 0.458 | nan | nan | 2155 |
| EGFR | lung | 126 | 12348 | 0.403 | 0.535 | nan | nan | 13 |
| CTNNB1 | pan-cancer | 1178 | 14873 | 0.366 | 0.681 | nan | nan | 1567 |
| CTNNB1 | lung | 126 | 12348 | 0.248 | 0.482 | nan | nan | 4 |
| MYC | pan-cancer | 1178 | 14873 | 0.659 | 0.415 | nan | nan | 2183 |
| MYC | lung | 126 | 12348 | 0.688 | 0.447 | nan | nan | 31 |

`EGFR`, `CTNNB1` and `MYC` are positive controls run through the identical code path; their scans are expected to return strong, interpretable partners.

Tables: `codependency_scan_stats.csv`, `codependency_top_hits.csv`, `codependency_top200_targets.csv`, `codependency_technical_floor.csv`. Figure: `fig2_codependency.png`.
