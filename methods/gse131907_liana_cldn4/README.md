# GSE131907 LIANA/LR — CLDN4-high malignant → T/NK

Additive public test: ligand–receptor scores from **CLDN4-high malignant**
cells to **T/NK** on Kim et al. GSE131907, with **honest n** (samples /
patients, not cells).

- `FINDING.md` — verdict, n table, LR table
- `METHODS.md` — labels, score, gates
- `results/lr_table.tsv` — pooled pair table (done criterion)

```bash
python3 -m pip install -r methods/gse131907_liana_cldn4/requirements.txt
python3 methods/gse131907_liana_cldn4/scripts/00_download.py
python3 methods/gse131907_liana_cldn4/scripts/01_analyze.py
```
