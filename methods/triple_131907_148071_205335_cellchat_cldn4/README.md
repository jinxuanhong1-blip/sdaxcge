# Triple merge: GSE131907 + GSE148071 + GSE205335, CLDN4-only CellChat

ADDITIVE. Patient is the unit. Scores the **triple** malignant-CLDN4 vs
same-patient T/NK table first, then Jin 2021 Hill outgoing ligands on
patients that pass floors.

```bash
python3 methods/triple_131907_148071_205335_cellchat_cldn4/analyze.py
python3 methods/triple_131907_148071_205335_cellchat_cldn4/analyze.py --skip-cellchat
```

Does **not** redo the 131907+205335-only CellChat. No GSE207422. No dual-high.
