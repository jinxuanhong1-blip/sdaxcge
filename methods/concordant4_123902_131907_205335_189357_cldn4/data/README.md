# Locked concordant-4 inputs

Public processed-matrix extracts only. Not a mega-merge. Not GSE148071 /
GSE127465 / GSE207422 / GSE154826 / CD45+ sets.

## Patient / donor / sample units (T/NK)

Locked PR #459 tables (not re-derived):

- `GSE123902_marker_units.tsv` — Laughney 2020 donors; marker-malignant.
- `GSE131907_samples.tsv` — Kim 2020; **sample** is the unit; author malignant.
- `GSE205335_patients.tsv` — advanced ICI; author malignant; RECIST not required.
- `GSE189357_marker_units.tsv` — Zhu/Wang AIS–IAC; patient unit; marker-malignant.

## Malignant UMI-sum matrices (tumor-cell-intrinsic DE)

Existing patient-pseudobulk counts from the pair/triple IFN-DE extras:

- `GSE123902_malignant_counts.tsv.gz` / `GSE123902_malignant_meta.tsv`
- `GSE189357_malignant_counts.tsv.gz` / `GSE189357_malignant_meta.tsv`
- `GSE131907_malignant_counts.tsv.gz`
- `GSE205335_malignant_counts.tsv.gz` / `GSE205335_malignant_meta.tsv` (P4001 absent)

`a8_sets.json` — Hallmark IFNα/γ, custom MHC-I/APM, KEGG/GO tight junction.
