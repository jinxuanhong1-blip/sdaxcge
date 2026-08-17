# GSE131907 Milo neighbourhoods vs malignant CLDN4

Additive to prior GSE131907 epithelial TACSTD2 work and to the GSE207422
TACSTD2 Milo folder. This folder tests **neighbourhood differential
abundance versus sample-level malignant CLDN4** on the Kim et al. 2020
LUAD atlas (GSE131907).

GSE207422 Milo is out of scope here (different agent / folder).

**Not miloR.** The environment has no R / Bioconductor / edgeR. Scripts
reimplement Milo neighbourhoods + SpatialFDR and use sample-level Spearman
or Welch tests on neighbourhood proportions. Do not cite the p-values as
miloR output.

The independent unit is the **sample**, not the cell.

```bash
pip install -r methods/gse131907_milo_cldn4/requirements.txt
python3 methods/gse131907_milo_cldn4/scripts/download.py
python3 methods/gse131907_milo_cldn4/scripts/knn_nhood.py
python3 methods/gse131907_milo_cldn4/scripts/run_gse131907.py
```

`FINDING.md` is written by the runner from the numbers that finished.
