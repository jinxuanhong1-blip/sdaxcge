# Inputs

- `GSE123902_marker_units.tsv` / `GSE205335_patients.tsv` — PR #459 locked malignant CLDN4 %pos (T/NK ρ taken as given, not re-scored here).
- `GSE123902_malignant_counts.tsv.gz` — donor UMI-sum of marker-malignant cells (PRIMARY preferred; normals dropped). Built by `scripts/build_gse123902.py`.
- `GSE205335_malignant_counts.tsv.gz` — patient UMI-sum of author malignant cells (same matrix as the winning-pair muscat extra; P4001 n_mal=27 absent).
- `a8_sets.json` — Hallmark IFN, custom MHC-I/APM, KEGG/GO tight junction, keratinization / KRT epithelial.
