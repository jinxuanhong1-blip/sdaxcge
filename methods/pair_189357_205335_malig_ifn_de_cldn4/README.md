# Pair GSE189357+GSE205335: CLDN4-only malignant IFN/MHC/TJ DE

ADDITIVE **CLDN4-only**. Not CellChat. No dual-high. No GSE148071.
Given PR #459 cut: %pos n=31 ρ=−0.478, Q4 r=−0.750 (T/NK ρ not re-audited).
Patient is the unit.

```bash
python3 methods/pair_189357_205335_malig_ifn_de_cldn4/download.py --skip-189357
python3 methods/pair_189357_205335_malig_ifn_de_cldn4/build_malignant_pseudobulk.py
python3 methods/pair_189357_205335_malig_ifn_de_cldn4/analyze.py
```

Writeup: `FINDING.md`.
Headline family DE table: `tables/family_de.tsv`.
