# Inputs

Locked PR #459 pair units (not re-audited):

- `GSE189357_marker_units.tsv` — Zhu/Wang AIS–IAC patient-level marker-malignant CLDN4.
- `GSE205335_patients.tsv` — Hui/Zhang ICI scRNA patient-level author-malignant CLDN4 (n=22; 4 extra GEO patients with 0 malignant already out).

Gene families:

- `a8_sets.json` — Hallmark IFN-α/γ, custom MHC-I/APM, KEGG/GO tight junction.

Built by `build_malignant_pseudobulk.py` from public GEO processed matrices (no FASTQ):

- `GSE189357_malignant_counts.tsv.gz` — patient UMI-sum of marker-malignant cells.
- `GSE205335_malignant_counts.tsv.gz` — patient UMI-sum of author `Malignant cells` (n=22 locked).
- `*_malignant_meta.tsv` — cells summed per unit.

GSE189357 marker-malignant (same gate as PR #459): `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`.
GSE205335 malignant = author `lineage.sub == Malignant cells`.
Patient is the unit. No dual-high TACSTD2×CLDN4. T/NK is not scored here.
