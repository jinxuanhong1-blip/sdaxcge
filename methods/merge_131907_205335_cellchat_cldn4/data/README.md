# Locked extracts (given; not re-audited)

- `GSE131907_samples.tsv` — PR #320 / #279 T/NK extract (sample-level).
  The companion between-unit table uses tumor-origin rows with
  `n_malignant >= 20` (author `Malignant cells` only; n=21).
- `GSE205335_patients.tsv` — PR #320 / #279 locked 22-patient table
  (≥20 author-malignant and ≥20 T/NK).

Primary paired CellChat uses GEO UMIs + author labels, including GSE131907
tS1/tS2/tS3, and does not treat these TSVs as the malignant definition.
