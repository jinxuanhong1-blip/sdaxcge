# Locked patient-level extracts (PR #459)

Not re-audited:

- `GSE189357_marker_units.tsv` — Zhu/Wang AIS–IAC marker-malignant patients (TD1–TD9).
- `GSE205335_patients.tsv` — Hu et al. author-malignant patients (n=22).

Combo Spearman / Q4 vs Q1 on these tables is **given**
(PR #459: %pos n=31 ρ=−0.478, Q4 r=−0.750) and is not re-ranked.

UMI matrices are downloaded at run time (`scripts/download.py`) and are not stored in git.
