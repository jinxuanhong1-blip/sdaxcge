# EXTRA series — GSE291670

Second public neoadjuvant ICI lung scRNA (not a GSE207422 audit).

```bash
# data/GSE291670/GSE291670_RAW.tar from GEO (122 MB), then:
tar -xf data/GSE291670/GSE291670_RAW.tar -C data/GSE291670
python3 scripts/extra/GSE291670/analyze.py \
  --datadir data/GSE291670 \
  --outdir results/extra/GSE291670
```
