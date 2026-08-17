# Triple GSE131907+GSE205335+GSE189357: CLDN4-only malignant IFN/MHC/TJ DE

ADDITIVE **CLDN4-only**. Not infiltrate. No dual-high. No GSE148071.

The pair GSE123902+GSE189357 %pos n=22 ρ=−0.638 (PR #459) already has a
malignant IFN DE (PR #469). This extra is the remaining 189357-axis triple.

Given PR #459 triple cut: %pos n=52 ρ=−0.497 (T/NK ρ not re-audited).
Patient/sample is the unit. P4001 (27 malignant cells) is out of malignant DE.

```bash
python3 methods/pair_or_triple_189357_malig_ifn_de_cldn4/analyze.py
```

Writeup: `FINDING.md`.
Headline family DE table: `tables/family_de.tsv`.
