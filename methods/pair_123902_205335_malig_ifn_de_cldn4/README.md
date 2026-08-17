# Pair GSE123902+GSE205335 malignant IFN/MHC/TJ/keratin DE (CLDN4-only)

ADDITIVE. Tumor-cell-intrinsic program DE. Not CellChat. Not T/NK infiltrate.
No dual-high. No GSE148071. Patient/donor is the unit.

Writeup: `FINDING.md`. Headline family table: `tables/family_de.tsv` (n / logFC / p).

```bash
python3 methods/pair_123902_205335_malig_ifn_de_cldn4/scripts/download_gse123902.py
python3 methods/pair_123902_205335_malig_ifn_de_cldn4/scripts/build_gse123902.py
python3 methods/pair_123902_205335_malig_ifn_de_cldn4/analyze.py
```
