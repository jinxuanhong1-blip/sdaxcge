# Pair GSE123902+GSE131907 — CLDN4-only REAL Slingshot / PAGA

ADDITIVE. PR #459 pair that differs. No dual-high. No GSE148071.

Root = GSE131907 nLung author AT2, never CLDN4-high.
Unit = GSE123902 donor + GSE131907 sample.
Primary clock = Slingshot. PAGA is geometry.

```bash
bash methods/pair_123902_131907_slingshot_cldn4/scripts/install_r_slingshot.sh
pip install -r methods/pair_123902_131907_slingshot_cldn4/requirements.txt
python3 methods/pair_123902_131907_slingshot_cldn4/scripts/download.py
python3 methods/pair_123902_131907_slingshot_cldn4/scripts/extract_epithelium.py
python3 methods/pair_123902_131907_slingshot_cldn4/scripts/analyze.py \
  --input /tmp/pair_123902_131907/epithelium.h5ad \
  --outdir methods/pair_123902_131907_slingshot_cldn4/results \
  --finding methods/pair_123902_131907_slingshot_cldn4/FINDING.md
```

Done when `results/tables/slingshot_lineages.tsv` and
`results/tables/sample_level_spearman.tsv` exist.
