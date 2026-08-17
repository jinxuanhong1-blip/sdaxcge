# Pair GSE123902+GSE189357 — CLDN4-only REAL Slingshot/PAGA

ADDITIVE. PR #459 pair that **differs** (n=22, %pos ρ=−0.638). Tails **7/5 thin**.
Root is GSE123902 NORMAL / not CLDN4-high. No dual-high gate.

Done when `results/tables/slingshot_lineages.tsv` exists.

```bash
pip install -r methods/pair_123902_189357_slingshot_cldn4/requirements.txt
bash methods/pair_123902_189357_slingshot_cldn4/scripts/install_r_slingshot.sh
python3 methods/pair_123902_189357_slingshot_cldn4/scripts/download.py \
  --out /tmp/geo_pair_123902_189357
python3 methods/pair_123902_189357_slingshot_cldn4/scripts/extract.py \
  --tars /tmp/geo_pair_123902_189357 \
  --out /tmp/geo_pair_123902_189357/epithelium.h5ad
python3 methods/pair_123902_189357_slingshot_cldn4/scripts/analyze.py \
  --input /tmp/geo_pair_123902_189357/epithelium.h5ad \
  --outdir methods/pair_123902_189357_slingshot_cldn4/results \
  --finding methods/pair_123902_189357_slingshot_cldn4/FINDING.md
```
