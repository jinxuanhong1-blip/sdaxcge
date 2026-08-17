# Inputs

Locked **CLDN4-only** units from PR #459 (GSE123902+GSE131907 %pos n=34).
No dual-high. No GSE148071. No T/NK infiltrate re-audit.

- `GSE123902_marker_units.tsv` — Laughney 2020 donor/sample marker-malignant CLDN4 %pos (PR #459).
- `GSE131907_samples.tsv` — Kim 2020 sample-level author-malignant CLDN4 %pos (PR #459).
- `GSE131907_malignant_counts.tsv.gz` — author malignant UMI-sum, 21 tumor samples (same matrix as PR #456).
- `GSE123902_malignant_counts.tsv.gz` — marker-malignant UMI-sum, 13 tumor/met donors (`build_gse123902_pseudobulk.py`).
- `a8_sets.json` — Hallmark IFN, custom MHC-I/APM, KEGG/GO tight junction, KRT_EPITHELIAL.

Patient/donor is the unit (GSE131907: sample; GSE123902: donor).
