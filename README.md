# CLDN4 spatial neighbor-graph analysis

Public-only squidpy/knn spatial graphs for **CLDN4-high tumor/epithelium** versus CD8 / NK / Treg / macrophages.

- Official CosMx NSCLC FFPE (`Lung5_Rep2`, He et al. 2022), **50 µm** radius
- Visium LUAD **GSE307534 GSM9226169 P1_LUAD** (invasive), **1–2 hop** hex graph

See [RESULTS.md](RESULTS.md) for enrichment heatmaps, co-occurrence-vs-distance curves, composition permutation tests, and MISTy-style CD8 variance partitioning.

```bash
python3 scripts/run_spatial_graph_cldn4.py
```
