# GSE127465 REAL Slingshot/PAGA — CLDN4 only

Additive public slice. Tumor epithelium from Zilionis inDrops (GSE127465) if CLDN4+ malignant cells exist. **n=7 is thin.** Root is not CLDN4-high. No dual-high gate.

Done when `results/tables/lineages.tsv` exists, or `FINDING.md` is a stop note with honest n=0.

```bash
bash methods/gse127465_slingshot_real_cldn4/scripts/install_tools.sh
python3 methods/gse127465_slingshot_real_cldn4/scripts/download.py --out /tmp/gse127465_slingshot
python3 methods/gse127465_slingshot_real_cldn4/scripts/extract_epithelium.py \
  --data /tmp/gse127465_slingshot \
  --out /tmp/gse127465_slingshot/epithelium.h5ad
python3 methods/gse127465_slingshot_real_cldn4/scripts/analyze.py \
  --input /tmp/gse127465_slingshot/epithelium.h5ad \
  --outdir methods/gse127465_slingshot_real_cldn4/results \
  --finding methods/gse127465_slingshot_real_cldn4/FINDING.md
```
