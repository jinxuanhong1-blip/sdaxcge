# GSE205335 — Milo / neighbourhood DA vs malignant CLDN4

Additive public-only analysis. Patient is the independent unit. miloR is not used.

```bash
pip install -r methods/gse205335_milo_cldn4/requirements.txt
python3 methods/gse205335_milo_cldn4/scripts/analyze.py
```

GEO processed UMI (~500 MB gzip, double-gzipped RDS) is downloaded to `/tmp/GSE205335` and is not committed. Controlled EGA raw (`EGAD00001008703`) is not accessed.

See `FINDING.md` for honest *n* and SpatialFDR.
