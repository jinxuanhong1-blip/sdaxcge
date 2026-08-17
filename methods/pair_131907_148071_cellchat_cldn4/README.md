# Pair GSE131907 + GSE148071 — CLDN4 vs T/NK + CellChat

ADDITIVE **CLDN4-only**. New pairwise merge of **GSE131907 + GSE148071**.
**Not +GSE205335.** No dual-high TACSTD2×CLDN4.

1. Patient-level malignant CLDN4 vs same-patient T/NK (honest n, Q4 vs Q1).
2. CellChat-style outgoing CLDN4-high → T/NK on the Q4 vs Q1 patients that pass the cell floors.

Writeup: `FINDING.md`. Combo ρ: `results/combo_rho.tsv`. LR table: `results/ligand_table.tsv`.

```bash
python3 methods/pair_131907_148071_cellchat_cldn4/download.py --out /tmp/pair_131907_148071
python3 methods/pair_131907_148071_cellchat_cldn4/analyze.py
# patient-level combo only:
python3 methods/pair_131907_148071_cellchat_cldn4/analyze.py --skip-cellchat
```
