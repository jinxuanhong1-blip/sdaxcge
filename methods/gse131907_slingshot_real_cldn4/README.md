# GSE131907 REAL Slingshot, CLDN4-only

Additive public slice. Epithelium / malignant only (Kim 2020).
Not a redo of PR #325 DPT. Root = nLung AT2, never CLDN4-high.
Sample is the unit. No dual-high gate.

```bash
pip install -r methods/gse131907_slingshot_real_cldn4/requirements.txt
python3 methods/gse131907_slingshot_real_cldn4/scripts/download.py \
  --out /tmp/gse131907_slingshot_data
python3 methods/gse131907_slingshot_real_cldn4/scripts/extract_epithelium.py \
  --data /tmp/gse131907_slingshot_data \
  --out /tmp/gse131907_slingshot_data/epithelium_malignant.h5ad
python3 methods/gse131907_slingshot_real_cldn4/scripts/analyze.py \
  --input /tmp/gse131907_slingshot_data/epithelium_malignant.h5ad \
  --outdir methods/gse131907_slingshot_real_cldn4/results \
  --finding methods/gse131907_slingshot_real_cldn4/FINDING.md
```

Done when `results/tables/slingshot_lineages.tsv` and
`results/tables/sample_level_spearman.tsv` exist.
