# Pair GSE123902+GSE189357: CLDN4-only malignant IFN/MHC/TJ DE

ADDITIVE **CLDN4-only**. Not CellChat. No dual-high.
Given PR #459 cut: %pos n=22 ρ=−0.638 (T/NK ρ not re-audited).
Patient/donor is the unit. Q4 vs Q1 tails are thin (7/5).

```bash
python3 methods/pair_123902_189357_malig_ifn_de_cldn4/download.py
python3 methods/pair_123902_189357_malig_ifn_de_cldn4/build_malignant_pseudobulk.py
python3 methods/pair_123902_189357_malig_ifn_de_cldn4/analyze.py
```

Writeup: `FINDING.md`.
Headline family DE table: `tables/de_q4q1_combined_families.tsv`.
