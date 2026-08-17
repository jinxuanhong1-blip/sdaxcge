# GSE205335 — CLDN4-only LIANA / CellPhoneDB-style LR

Additive slice. **CLDN4 only** — TACSTD2 is not a gate. Dual-high is not run.
Not GSE207422 (PR #344).

Primary question: do CLDN4-high author-malignant cells send weaker T-recruit /
MHC-I signals to same-patient T/NK than CLDN4-low cells on the public
GSE205335 processed UMI?

Read `FINDING.md`. The ligand–receptor table is
`results/lr_table_cldn4_outgoing_tnk.tsv`.

```bash
cd methods/gse205335_liana_cldn4
python3 scripts/00_download.py
python3 scripts/01_build_pairs.py
python3 scripts/02_run_ccc.py
```
