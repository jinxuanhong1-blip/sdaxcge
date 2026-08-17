# Locked patient/sample extracts (not re-audited)

From PR #459 `methods/scrna_cldn4_combo_enum/data/`:

- `GSE123902_marker_units.tsv` — Laughney 2020 marker-malignant donor table.
  Eligible tumor/met donors (n=13) after dropping normals.
- `GSE131907_samples.tsv` — Kim 2020 author-malignant sample table
  (PR #279 / #320 extract). Eligible tumor-origin samples with
  `n_malignant ≥ 20` (n=21). Sample is the unit.

The given pair %pos Spearman uses `mal_CLDN4_pct` vs same-unit `frac_tnk`.
This folder does not re-score those vectors. GSE148071 is not here.
