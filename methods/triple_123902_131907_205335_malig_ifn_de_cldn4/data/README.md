# Inputs

Locked **CLDN4-only** units from PR #459 (GSE123902+GSE131907+GSE205335 %pos n=56).
Not the 7-pool. Not +GSE148071. No dual-high. No T/NK infiltrate re-audit.

- `GSE123902_marker_units.tsv` — Laughney 2020 donor/sample marker-malignant CLDN4 %pos (PR #459).
- `GSE131907_samples.tsv` — Kim 2020 sample-level author-malignant CLDN4 %pos (PR #459).
- `GSE205335_patients.tsv` — Hu 2023 patient-level author-malignant CLDN4 %pos (PR #459 / #320).
- `GSE123902_malignant_counts.tsv.gz` — marker-malignant UMI-sum, 13 tumor/met donors (`build_gse123902_pseudobulk.py`).
- `GSE131907_malignant_counts.tsv.gz` — author malignant UMI-sum, 21 tumor samples (same matrix as PR #456 / #472).
- `GSE205335_malignant_counts.tsv.gz` — author malignant UMI-sum, 21 patients (same matrix as PR #456). P4001 (27 malignant cells) is in the label table, not in this matrix.
- `GSE205335_malignant_meta.tsv` — audit of the 21 patients in the count matrix.
- `a8_sets.json` — Hallmark IFN, custom MHC-I/APM, KEGG/GO tight junction, KRT_EPITHELIAL.

Patient/donor is the unit (GSE131907: sample; GSE123902: donor; GSE205335: patient).
