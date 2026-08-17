# Locked patient-level extracts (given combo; do not re-audit ρ)

Copied from PR #459 `methods/scrna_cldn4_combo_enum/data/`.

- `GSE123902_marker_units.tsv` — Laughney 2020 dense CSV; marker-malignant; donor unit.
- `GSE189357_marker_units.tsv` — Zhu/Wang AIS–IAC 10x MTX; marker-malignant; patient unit.

These lock the given Spearman members (n=13 + n=9 = 22). CellChat uses the
public processed matrices downloaded by `scripts/download.py` (not these TSVs).
