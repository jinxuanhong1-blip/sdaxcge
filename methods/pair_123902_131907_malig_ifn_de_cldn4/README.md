# Pair GSE123902+GSE131907 — malignant-cell-intrinsic CLDN4 DE

ADDITIVE **CLDN4-only**. Tumor-cell-intrinsic program DE.
Not CellChat. Not T/NK infiltrate. No dual-high. No GSE148071.

Writeup: `FINDING.md`.

```bash
# optional rebuild of GSE123902 marker-malignant UMI-sum
bash methods/pair_123902_131907_malig_ifn_de_cldn4/download.sh
python3 methods/pair_123902_131907_malig_ifn_de_cldn4/build_gse123902_pseudobulk.py
python3 methods/pair_123902_131907_malig_ifn_de_cldn4/analyze.py
```
