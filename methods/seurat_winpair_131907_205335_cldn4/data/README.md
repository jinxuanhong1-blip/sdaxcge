# Inputs (public GEO + locked PR #320 tables)

- `GSE131907_samples.tsv` / `GSE205335_patients.tsv` — PR #320 author-malignant CLDN4 %pos and T/NK fractions (winning pair).
- `GSE205335_gsm_sample_metadata.csv` — GSM → patient / `orig.ident`.
- `a8_sets.json` — Hallmark IFN, custom MHC-I (same freeze as the muscat/pseudobulk extras).

GEO files are downloaded to `/tmp/winpair_geo` (not committed):

- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/
