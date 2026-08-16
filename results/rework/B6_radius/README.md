# B6_radius — tumor-only CLDN4 vs immune neighborhoods at hex rings 1/2/3

Self-contained rework of claim **B6** after PR67 (`fable_spatial`) reported Visium CLDN4 vs immune-neighborhood **|median partial ρ| ≤ 0.06**.

This folder is the only write target. Do not invent results. Thresholds are pre-specified and are not tuned to enlarge |ρ|.

## Rework rules

1. **Tumor-spot only** (E-MTAB-13530 `P*_T*` sections; epithelial-high spots).
2. **Radii = Visium hex rings 1 / 2 / 3** (self excluded).
3. **Immune-rich spots are excluded from the CLDN4 (and TACSTD2) score** (immune score ≥ section 75th percentile).
4. **GeoMx (GSE271689) = tumor / PanCK compartment only** (no stromal pairing; no hex rings exist).

## Run

```bash
pip install -r results/rework/B6_radius/requirements.txt
python3 results/rework/B6_radius/download.py
python3 results/rework/B6_radius/analyze.py
```

Processed inputs cache under `/tmp/b6_radius_data` (not committed). Outputs (tables, figures, `summary.json`, `WRITEUP.md`) stay in this folder.
