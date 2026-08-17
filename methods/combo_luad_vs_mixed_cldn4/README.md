# LUAD-only vs mixed-all vs drop-SCLC (CLDN4-only)

Additive combinatorial **histology cut**, not a bigger merge.

Three public processed series: GSE131907 (LUAD), GSE148071 (NSCLC, unlabeled
in the extract), GSE205335 (mixed, some SCLC).

Malignant CLDN4 vs T/NK. Honest n. Q4 vs Q1 extra. No dual-high.

```bash
python3 methods/combo_luad_vs_mixed_cldn4/analyze.py
```

Answer: `FINDING.md` three-row table. Extra forest:
`figures/forest_three_cuts_spearman.png`.
