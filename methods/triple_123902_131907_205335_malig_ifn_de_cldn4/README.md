# Triple GSE123902+GSE131907+GSE205335 — malignant-cell-intrinsic CLDN4 DE

ADDITIVE **CLDN4-only**. Tumor-cell-intrinsic IFN/MHC/TJ DE on the PR #459
triple that **differs**. Not the 7-pool. Not +GSE148071. Not CellChat.
No dual-high.

Writeup: `FINDING.md`. Headline: `tables/family_de.tsv`.

```bash
# optional rebuild of GSE123902 marker-malignant UMI-sum
bash methods/triple_123902_131907_205335_malig_ifn_de_cldn4/download.sh
python3 methods/triple_123902_131907_205335_malig_ifn_de_cldn4/build_gse123902_pseudobulk.py
python3 methods/triple_123902_131907_205335_malig_ifn_de_cldn4/analyze.py
```
