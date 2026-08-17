# Triple-merge CLDN4-only LIANA / CellPhoneDB

Additive patient-level ligand–receptor scoring on the public merge
**GSE131907 + GSE148071 + GSE205335**.

- **CLDN4 only.** TACSTD2 is not a gate. Dual-high is not run.
- **Not** GSE207422. **Not** the 131907+205335-only pair merge.
- Inferential unit = patient. Cells are not replicates.
- Primary method = documented CellPhoneDB-style mean-of-means on log1p(CP10k).
- LIANA `mt.cellphonedb` is secondary if importable. CellChat is not run.

```bash
cd methods/triple_scrna_liana_cldn4
python3 scripts/00_download.py
python3 scripts/01_build_pairs.py
python3 scripts/02_extract.py
python3 scripts/03_run_ccc.py
```

Primary table: `results/lr_table_cldn4_outgoing_tnk.tsv`.
Note: `FINDING.md`.
