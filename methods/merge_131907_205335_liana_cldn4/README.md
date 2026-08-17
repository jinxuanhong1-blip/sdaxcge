# Merged GSE131907 + GSE205335 CLDN4-only LIANA / CellPhoneDB

Outgoing **CLDN4-high malignant → same-patient T/NK** on the PR #320
slice that already differs (Q4 vs Q1 n=23 r=−0.705 is given; not
re-audited). No dual-high. No GSE207422.

```bash
python3 methods/merge_131907_205335_liana_cldn4/scripts/00_download.py
python3 methods/merge_131907_205335_liana_cldn4/scripts/01_analyze.py
```

Writeup: `FINDING.md`. Primary table: `results/lr_table.tsv`.
