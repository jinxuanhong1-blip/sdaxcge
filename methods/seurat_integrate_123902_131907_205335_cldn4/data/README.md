# Inputs (public GEO + locked PR #459 tables)

ADDITIVE. **CLDN4-only.** The triple that differs is given (PR #459).
Do **not** add GSE148071 or GSE189357.

- `GSE123902_marker_units.tsv` — Laughney 2020 donor-level marker-malignant
  (PRIMARY preferred over METASTASIS; normals dropped). Eligible n=13.
- `GSE131907_samples.tsv` — Kim 2020 author-malignant sample-level extract
  (PR #279 / #320 / #459). Eligible n=21 with n_malignant ≥ 20.
- `GSE205335_patients.tsv` — Hu 2023 author-malignant patients
  (PR #279 / #320 / #459). Eligible n=22.
- `GSE205335_gsm_sample_metadata.csv` — GSM → patient / `orig.ident`.
- `a8_sets.json` — Hallmark IFN, custom MHC-I (same freeze as the muscat extras).

GEO files are downloaded to `/tmp/triple_geo` (not committed):

- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/
