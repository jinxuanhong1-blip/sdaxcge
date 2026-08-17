# Winning-pair real Slingshot (CLDN4 only)

ADDITIVE CLDN4-only Slingshot + PAGA on GSE131907 + GSE205335 malignant/epithelial
cells. Root is AT2-like or lowest-CLDN4, never CLDN4-high.

PR #320 T/NK r=−0.705 is given and is not re-audited. No dual-high. No GSE148071.

See `FINDING.md` and `METHODS.md`.

```bash
pip install -r methods/winpair_slingshot_real_cldn4/requirements.txt
python3 methods/winpair_slingshot_real_cldn4/scripts/download.py --out /tmp/winpair_slingshot_real
python3 methods/winpair_slingshot_real_cldn4/scripts/extract.py --data /tmp/winpair_slingshot_real --out /tmp/winpair_slingshot_real/epithelium.h5ad
python3 methods/winpair_slingshot_real_cldn4/scripts/analyze.py \
  --input /tmp/winpair_slingshot_real/epithelium.h5ad \
  --outdir methods/winpair_slingshot_real_cldn4/results \
  --finding methods/winpair_slingshot_real_cldn4/FINDING.md
```
