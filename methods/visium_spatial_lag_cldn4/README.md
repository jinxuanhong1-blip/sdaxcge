# Public Visium spatial-lag CLDN4 vs neighbor CD8A

Retune of the Visium exclusion test. Same-spot Spearman and nearest-spot µm
are not the punch (55 µm spots mix tumor + T). Primary metrics:

1. Lag-ρ: index epithelial-like CLDN4 vs ring-1 / 80–150 µm neighbor CD8A
2. Q4 vs Q1 mean neighbor CD8A (paired by section)
3. Morisita–Horn mixing
4. KRT8 residual on the lag
5. Drop no-contrast sections
6. Optional tumor–stroma interface

```bash
python3 methods/visium_spatial_lag_cldn4/scripts/download_public_visium.py
python3 methods/visium_spatial_lag_cldn4/scripts/analyze_spatial_lag.py
```

See `METHODS.md` and `RESULTS.md`.
