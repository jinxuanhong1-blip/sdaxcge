# Results · Milo neighborhood DA (GSE207422)

Filled after `python3 methods/scrna_milo/02_run_milo.py`. See `methods/scrna_milo/README.md` for the test and the SpatialFDR definition.

Primary objects:

- `summary.json` — n per arm, n neighborhoods, SpatialFDR counts
- `sample_table.tsv` — patient-level malignant TACSTD2 split and MPR labels
- `nhoods_{pca,harmony}_{tacstd2_high_vs_low,mpr_vs_nmpr}.tsv` — all neighborhoods
- `kept_*.tsv` — T/NK-depleted interface nhoods next to TACSTD2-high epithelium
- `fig_*.png` — volcano (SpatialFDR) and interface scatter
