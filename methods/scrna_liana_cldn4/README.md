# CLDN4-only LIANA / CellPhoneDB-style LR on GSE207422

Additive slice. **CLDN4 only** — TACSTD2 is not a gate.

Primary question: do CLDN4-high malignant-like cells send weaker T-recruit / MHC-I
signals to T/NK than CLDN4-low cells on the public GSE207422 UMI?

Read `FINDING.md`. The ligand–receptor table is
`results/lr_table_cldn4_outgoing_tnk.tsv`.

```bash
python3 scripts/00_download.py
python3 scripts/01_build_pairs.py
python3 scripts/02_run_ccc.py
```
