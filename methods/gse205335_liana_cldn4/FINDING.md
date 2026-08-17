# GSE205335 — LIANA/CellPhoneDB-style LR from CLDN4-high malignant to T/NK

**Status:** pipeline is on the branch. `results/lr_table_cldn4_outgoing_tnk.tsv`
is written by `scripts/02_run_ccc.py` after the public processed UMI is scored.
This note will be replaced with computed numbers (not placeholders).

**CLDN4 only.** TACSTD2 is not a gate. Dual-high is not run. Not GSE207422
(PR #344). Patient is the unit.

```bash
cd methods/gse205335_liana_cldn4
python3 scripts/00_download.py
python3 scripts/01_build_pairs.py
python3 scripts/02_run_ccc.py
```
