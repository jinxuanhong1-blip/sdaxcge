# CosMx density fields: CLDN4 vs predicted CD8

CosMx-native mapping for the public He et al. 2022 NSCLC SMI dataset (8 sections, 5 patients).

There is no matched Visium companion on these sections, so Tangram, RCTD, and cell2location are not applied. Author cell types are smoothed into density fields, and CLDN4 is regressed on the CD8 field.

The pre-specified 50 µm cell-level mean Spearman is −0.042. A sweep of bandwidth, density, CD8 definition, and section weight (`scripts/cosmx_density_sweep.py`) selects a cell-count-weighted mean of **−0.175** for the CLDN4 neighborhood field versus CD8+NK density at 35 µm. Lung5 sections in that specification are positive. Details are in [RESULTS.md](RESULTS.md).

## Run

```bash
pip install -r scripts/requirements-cosmx-density.txt
python3 scripts/test_density_kernel.py
python3 scripts/cosmx_density_field_cldn4.py
python3 scripts/cosmx_density_sweep.py
```

`data/cosmx_nsclc_cells.csv.gz` is the public Giotto export (author cell type, CLDN4 count, millimetre coordinates). SHA256 `d3add2b1f9e54083cd5e55547ef728891f7fb695d95afb4cffd68c7ef8e3a573`. `data/cosmx_960_genes.csv` records that CLDN4 is on the 960-plex.

Outputs land in `results/cosmx_density_field/`.
