# QUAD merge CLDN4-only LIANA / CellPhoneDB

Additive patient-level ligand–receptor scores from **CLDN4-high malignant
cells to T/NK** on the public merge of GSE207422 + GSE131907 + GSE148071 +
GSE205335.

CLDN4 only. TACSTD2 is not a gate. Dual-high is not run. This is not a
re-audit of the GSE207422-only LIANA (PR #344).

```bash
cd methods/quad_207422_liana_cldn4
python3 scripts/00_download.py
python3 scripts/01_build_pairs.py
python3 scripts/02_run_ccc.py
```

Done criterion: `results/lr_table_cldn4_outgoing_tnk.tsv` and `FINDING.md`.
