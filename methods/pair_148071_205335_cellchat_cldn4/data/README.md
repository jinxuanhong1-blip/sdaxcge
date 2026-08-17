# Locked extracts (not re-ranked)

- `GSE148071_per_sample.tsv` — 42-biopsy epithelial / T/NK counts and CLDN4
  scores from the GSE148071 CellChat folder (PR #348). Eligibility is
  ≥25 marker-argmax epithelial **and** ≥25 T/NK. Epithelium is
  **putative** malignant (no GEO labels / CopyKAT).
- `GSE205335_patients.tsv` — locked 22-patient author-malignant extract
  (PR #279 / #320 / #362). Eligibility was ≥20 author-malignant **and**
  ≥20 T/NK.
- `GSE205335_per_patient_lr.tsv.gz` — same-patient CellChat-style Hill
  probabilities from PR #362. Reused so the 500 MB GEO RDS is not
  downloaded again.

GSE131907 is **not** in this pair. No dual-high TACSTD2×CLDN4 score.
The 172 MB GSE148071 RAW tar is downloaded at run time and is not stored
in git.
