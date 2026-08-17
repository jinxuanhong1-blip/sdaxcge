# QUAD merge — LIANA/CellPhoneDB-style LR from CLDN4-high malignant to T/NK

**Status:** pipeline is on the branch. `results/lr_table_cldn4_outgoing_tnk.tsv`
is written by `scripts/02_run_ccc.py` after the four public objects are scored.
This note will be replaced with computed numbers (not placeholders).

**CLDN4 only.** TACSTD2 is not a gate. Dual-high is not run. GSE207422 is
included in the merge because the user asked for it. This is not a re-audit
of the GSE207422-only LIANA (PR #344). Patient is the unit.

```bash
cd methods/quad_207422_liana_cldn4
python3 scripts/00_download.py
python3 scripts/01_build_pairs.py
python3 scripts/02_run_ccc.py
```
