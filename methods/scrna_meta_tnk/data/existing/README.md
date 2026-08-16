# Per-patient tables (processed GEO)

These tables are the patient-level inputs for `ingest_cohorts.py`. They come from public processed GEO analyses already in this repository. GSE207422 is the A3-given DRMref table and is not re-audited here.

| File | Cohort | Provenance branch / note |
|---|---|---|
| `GSE207422_drmref_patients.tsv` | GSE207422 | A3 given (DRMref 12 post-tx patients) |
| `GSE205335_patients.tsv` | GSE205335 | Author malignant labels |
| `GSE131907_samples.tsv` | GSE131907 | Author cell types; tumor-site filter applied at ingest |
| `GSE253013_patients.tsv` | GSE253013 | Tumor + ANT; ingest keeps Tumor + eligible_malig |
| `GSE291670_patients.tsv` | GSE291670 | Marker malignant; lineage T/NK |
| `GSE325414_donors.csv` | GSE325414 | 2026 leftover; author malignant / T/NK; donor unit |
