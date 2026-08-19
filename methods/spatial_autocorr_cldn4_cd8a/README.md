# CLDN4-only high-end spatial autocorrelation vs CD8A

Open Visium LUAD/NSCLC. Per section: bivariate Moran's I, Lee's L, KRT8-residual partial Moran, Gi* overlap (CLDN4-hot ∩ CD8-cold), SLX and spatial-error regression.

```bash
bash methods/spatial_autocorr_cldn4_cd8a/download.sh
python3 methods/spatial_autocorr_cldn4_cd8a/analyze.py
```

See `RESULTS.md`. Raw downloads stay under `/data` (gitignored). GSE307534 (9.4 GB) is not fetched.
