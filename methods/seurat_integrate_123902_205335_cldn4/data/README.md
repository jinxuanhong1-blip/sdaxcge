# Inputs (public GEO + locked pair tables)

ADDITIVE. **GSE123902 + GSE205335 only.** Not GSE148071. Not GSE127465.
Not the concordant-4 four-way merge.

- `GSE123902_marker_units.tsv` — PR #459 per-library marker-malignant counts
  (eligibility / filename map). Seurat lineage is recomputed; this is not the
  primary score table.
- `GSE205335_patients.tsv` — PR #459 / #320 patient-level author-malignant
  CLDN4 and T/NK (eligibility + locked full-sample `frac_tnk`).
- `GSE205335_gsm_sample_metadata.csv` — GSM → patient / `orig.ident` / tissue.

GEO files are downloaded to `/tmp/geo_pair_123902_205335` (not committed):

- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/
