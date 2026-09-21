# HTAN open lung scRNA: CLDN4 vs DNA damage / STING

Portal snapshot date: 2026-09-21. Database `htan_2026_932`, table `files`, filter `organType` has Lung and `assayName = scRNA-seq`.

- `lung_scrna_files.tsv` — one row per portal file (2,390 files).
- `inventory_summary.tsv` — counts by atlas, level, download source, format.
- Phase-2 database `htan2_2026_926` had no lung RNA rows on that date.

Anonymous Synapse file download (`repo-prod` `/entity/{id}/file`) returned HTTP 403. The two combined Level-4 h5ads used in the score analysis are the CELLxGENE copies of `syn23626795` (Chan combined) and `syn51033592` (Glasner human LUAD). Per-sample BU and HTAPP matrices stayed on Synapse and were not scored.

Refresh the summary with `python3 scripts/htan_cldn4_sting/00_summarize_inventory.py`.
