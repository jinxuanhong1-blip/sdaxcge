# Locked inputs

`locked_patient_units.tsv` is the patient/sample table from the concordant-4 Seurat analysis (PR #539). Columns used here are `dataset`, `unit_id`, and `frac_tnk` (plus malignant counts for audit).

That table is not re-derived. The Spearman of malignant CLDN4 %pos versus T/NK (n=65, ρ=−0.531) is not recomputed in this folder.

Raw GEO matrices are downloaded by `scripts/download.py` and are not committed.
