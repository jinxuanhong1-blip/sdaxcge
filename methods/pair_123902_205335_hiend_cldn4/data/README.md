# Locked patient-level extracts (not re-audited)

From PR #459 `methods/scrna_cldn4_combo_enum/data/`:

- `GSE123902_marker_units.tsv` — Laughney 2020 marker-malignant donor table.
  Eligible tumor/met donors (n=13) after dropping normals.
- `GSE205335_patients.tsv` — author-malignant patients (n=22).

The given pair %pos Spearman uses `mal_CLDN4_pct` (GSE123902) and
`mal_CLDN4_pct_pos` (GSE205335) vs same-unit `frac_tnk`.
This folder does not re-score those vectors.
