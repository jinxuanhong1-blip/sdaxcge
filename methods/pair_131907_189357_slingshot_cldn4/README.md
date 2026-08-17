# Pair GSE131907+GSE189357 CLDN4 Slingshot/PAGA

ADDITIVE **CLDN4-only** real Slingshot + PAGA on the PR #459 pair that already
differs. Root is GSE131907 nLung AT2, **not** CLDN4-high. No dual-high.
No GSE148071.

```bash
pip install -r methods/pair_131907_189357_slingshot_cldn4/requirements.txt
# R 4.3 + Bioconductor slingshot in ~/R/library
python3 methods/pair_131907_189357_slingshot_cldn4/scripts/run_all.py
```

Done when `results/tables/lineage_table.tsv` exists.
