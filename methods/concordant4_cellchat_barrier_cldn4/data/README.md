Locked unit tables from prior CLDN4-only extracts (PR #459 members).

- `GSE123902_marker_units.tsv` — donor-level; PRIMARY preferred over METASTASIS.
- `GSE131907_samples.tsv` — author malignant + T/NK; lock n_malignant>0 and n_tnk>0.
- `GSE205335_patients.tsv` — author malignant + T/NK; lock n_malignant>0 and n_tnk>0.
- `GSE189357_marker_units.tsv` — TD1–TD9 marker-malignant.

Do not add GSE148071 / GSE127465 / CD45-only.
