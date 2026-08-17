# Inputs

Locked PR #459 units (T/NK ρ not re-audited):

- `GSE131907_samples.tsv` — Kim 2020 author-malignant sample-level CLDN4 (PR #279 / #320 / #459).
- `GSE205335_patients.tsv` — Hu 2022 author-malignant patient-level CLDN4 (PR #279 / #320 / #459).
- `GSE189357_marker_units.tsv` — Zhu/Wang AIS–IAC marker-malignant patient-level CLDN4 (PR #459).

Gene families:

- `a8_sets.json` — Hallmark IFN-α/γ, custom MHC-I/APM, KEGG/GO tight junction (same file as the pair DE extras).

Malignant UMI-sum matrices (reused; not rebuilt here):

- `GSE131907_malignant_counts.tsv.gz` — author-malignant sample UMI-sum (PR #456 / #472).
- `GSE205335_malignant_counts.tsv.gz` — author-malignant patient UMI-sum (PR #456). P4001 (27 cells) is absent (n_mal≥30).
- `GSE189357_malignant_counts.tsv.gz` — marker-malignant patient UMI-sum (PR #469). Gate: `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`.
- `*_malignant_meta.tsv` — cells summed per unit.

No GSE148071. No dual-high TACSTD2×CLDN4. Patient/sample is the unit.
