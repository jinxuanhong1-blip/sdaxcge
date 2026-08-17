# QUAD CellChat-style CLDN4 (GSE207422 + GSE131907 + GSE148071 + GSE205335)

ADDITIVE. **CLDN4 only.** Patient is the unit. GSE207422 is **in** the merge.

```bash
python3 methods/quad_207422_131907_148071_205335_cellchat_cldn4/download.py --out /tmp/quad_cldn4_geo
python3 methods/quad_207422_131907_148071_205335_cellchat_cldn4/analyze.py
```

Write-up: `FINDING.md`. Combo rho: `tables/combo_rho.tsv`. LR table: `tables/ligand_table.tsv`.
