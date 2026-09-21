# Inputs

Malignant UMI sums and the locked CLDN4 percent-positive tables are the concordant-4 matrices from PR #503 (`cursor/concordant4-cldn4-tnk-muscat-f40c`):

- `GSE123902_malignant_counts.tsv.gz` and `GSE123902_marker_units.tsv`
- `GSE131907_malignant_counts.tsv.gz` and `GSE131907_samples.tsv`
- `GSE205335_malignant_counts.tsv.gz` and `GSE205335_patients.tsv`
- `GSE189357_malignant_counts.tsv.gz` and `GSE189357_marker_units.tsv`
- `a8_sets.json` (Hallmark IFN-α, Hallmark IFN-γ, custom MHC-I/APM)

`patient_scores_pr638.tsv` is the cell-level malignant percent-positive table from PR #638. It is the TACSTD2 %pos source. It is not a re-count of cells.

P4001 is in the GSE205335 phenotype vector and is absent from the UMI sum. Expression analyses use the 64 units that have a count column.

These files are derived patient-level sums already used in the concordant-4 PRs. They are not a new GEO download and they are not the private 8-KL matrices.
