# merge_131907_205335_cellchat_cldn4

ADDITIVE **CLDN4-only** CellChat-style ligand–receptor on the merged public
UMI slice **GSE131907** (author malignant / tS*) + **GSE205335** (author
malignant). Patient is the unit. GSE207422 is not run.

```bash
python3 methods/merge_131907_205335_cellchat_cldn4/scripts/download.py
python3 methods/merge_131907_205335_cellchat_cldn4/scripts/analyze.py
```

Primary deliverable: `results/ligand_table.tsv` (n / Δ / p) and `FINDING.md`.
