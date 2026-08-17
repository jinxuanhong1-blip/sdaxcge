# Combo notes used (not re-audited)

The marker-malignant combo that already differs is taken from PR #290 / #279:

- **Marker-malignant %pos vs T/NK: GSE253013 + GSE291670**
- **k=2 · N=15 · CLDN4 ρ=−0.714 (p=0.0217, I²=23%)**

`samples_n15.tsv` lists the 15 units already scored in those notes:

- GSE253013: 9 **tumor** patients (MRC001, MRC002, MRC003, MRC004, MRC006, MRC007, MRC008, MRC009, MRC010). Adjacent-normal (ANT) rows are not in the combo.
- GSE291670: 6 post-treatment tumors (MPR-1/2/3, Non-MPR-1/2/3).

This folder does **not** re-compute that Spearman. Adjacent-normal GSE253013 libraries and any extra GEO samples are not added. n stays 15.
