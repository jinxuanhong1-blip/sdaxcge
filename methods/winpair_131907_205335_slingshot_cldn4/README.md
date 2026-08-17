# Winning-pair CLDN4 trajectory (GSE131907 + GSE205335 epithelium)

ADDITIVE CLDN4-only Slingshot/Palantir-style trajectory on the winning pair.
Does not redo GSE131907-only PAGA (PR #325). Does not add GSE207422.

If Slingshot R is missing, documented AT2-rooted diffusion pseudotime is the clock.

```bash
pip install -r methods/winpair_131907_205335_slingshot_cldn4/requirements.txt
python3 methods/winpair_131907_205335_slingshot_cldn4/scripts/download.py --out /tmp/winpair_131907_205335
python3 methods/winpair_131907_205335_slingshot_cldn4/scripts/extract_epithelium.py \
  --data /tmp/winpair_131907_205335 \
  --out /tmp/winpair_131907_205335/epithelium.h5ad
python3 methods/winpair_131907_205335_slingshot_cldn4/scripts/analyze.py \
  --input /tmp/winpair_131907_205335/epithelium.h5ad \
  --outdir methods/winpair_131907_205335_slingshot_cldn4/results \
  --finding methods/winpair_131907_205335_slingshot_cldn4/FINDING.md
```

Done when `results/tables/sample_level_spearman.tsv` exists.
