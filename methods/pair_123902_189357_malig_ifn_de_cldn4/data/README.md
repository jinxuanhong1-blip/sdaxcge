# Inputs

Locked PR #459 pair units (not re-audited):

- `GSE123902_marker_units.tsv` — Laughney 2020 donor-level marker-malignant CLDN4 (normals dropped in analysis).
- `GSE189357_marker_units.tsv` — Zhu/Wang AIS–IAC patient-level marker-malignant CLDN4.

Gene families:

- `a8_sets.json` — Hallmark IFN-α/γ, custom MHC-I/APM, KEGG/GO tight junction (same file as the winning-pair DE extra).

Built by `build_malignant_pseudobulk.py` from public GEO processed matrices (<2 GB):

- `GSE123902_malignant_counts.tsv.gz` — donor UMI-sum of marker-malignant cells (tumor/met only).
- `GSE189357_malignant_counts.tsv.gz` — patient UMI-sum of marker-malignant cells.
- `*_malignant_meta.tsv` — cells summed per unit.

Marker-malignant (same gate as PR #459): `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`.
Patient/donor is the unit. No dual-high TACSTD2×CLDN4.
