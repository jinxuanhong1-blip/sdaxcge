# align_tj — TROP2 / tight-junction public alignment

See [WRITEUP.md](WRITEUP.md) for the honest overlap call.

```
python3 scripts/01_tcga_prep_diff_coexpr.py
python3 scripts/02_tcga_gsea.py
python3 scripts/03_build_genesets.py
python3 scripts/04_gse207422_sc.py
python3 scripts/05_coxpresdb_human_mouse.py
python3 scripts/06_gse131907_malig.py
python3 scripts/07_tabula_muris_lung.py
python3 scripts/08_gse207422_bulk.py
python3 scripts/09_intersection.py
```

Inputs download into `data/` (gitignored). Outputs go to `tables/`.
