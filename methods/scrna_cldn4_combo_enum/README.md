# CLDN4-only pair/triple enumeration vs same-patient T/NK

Additive. CLDN4 only. No dual-high. No CellChat.

Bigger merge is not the answer. This folder enumerates **every pair and every triple** that can actually be scored from the listed public processed scRNA cohorts.

Write-up: [`FINDING.md`](FINDING.md). Complete table: [`tables/combo_table.tsv`](tables/combo_table.tsv). Extra forests: [`figures/`](figures/).

```bash
python3 methods/scrna_cldn4_combo_enum/analyze.py
```
