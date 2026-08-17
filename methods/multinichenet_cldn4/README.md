# Multi-sample NicheNet — CLDN4-only malignant → same-patient T/NK

Reusable Python pipeline for **multi-sample / MultiNicheNet-style** ligand
activity. **CLDN4 only.** No dual-high TACSTD2×CLDN4. Patient is the unit.

This run: PR #320 winning public merge **GSE131907 + GSE205335**, plus
**GSE207422** as a third cohort.

R `nichenetr` / `multinichenetr` are not required. Scoring uses the published
NicheNet-v2 ligand–target prior (Zenodo 7074291) in Python.

```bash
python3 scripts/00_download.py
python3 scripts/01_convert_prior.py
python3 scripts/02_extract_panels.py
python3 scripts/03_analyze.py
```

Writeup: [`FINDING.md`](FINDING.md). Ligand table: [`results/ligand_activity.tsv`](results/ligand_activity.tsv).
