# CosMx NSCLC CLDN4–CD8 spatial point-process analysis

Cross-type Ripley's K, L, and pair-correlation g(r) between **CLDN4-high tumor cells** and **CD8 T cells** on the official CosMx NSCLC FFPE 960-plex dataset (8 samples / 5 patients; He et al. 2022).

See **[RESULTS.md](RESULTS.md)** for numbers, envelopes, and per-FOV meta-analysis.

## Run

```bash
python scripts/test_csr_calibration.py
python scripts/run_fov_analysis.py
```

The analysis-ready table `data/cosmx_nsclc_cells.csv.gz` is exported from NanoString's public `All SMI Giotto object.tar.gz` (author cell types + CLDN4 counts + millimetre coordinates).
