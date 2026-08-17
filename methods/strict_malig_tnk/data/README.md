# Reused scRNA patient tables

Copied from existing repo extracts. This analysis does **not** re-download
GEO or re-annotate cells.

Kept for scoring (author-malignant, marker-malignant, or DRMref):

- `GSE207422_drmref_patients.tsv` — A3 given DRMref table (PR #279)
- `GSE207422_marker_fable.tsv` / `GSE207422_marker_nsclc.tsv` — marker-malignant extras
- `GSE205335_patients.tsv` — author malignant
- `GSE291670_patients.tsv` — marker malignant
- `GSE253013_patients.tsv` — marker malignant-like (tumor)
- `GSE131907_samples.tsv` — author malignant subtype (epithelial columns ignored)
- `GSE325414_donors.csv` — author malignant

`leftover/` is the drop audit only:

- GSE241934 IIT/Real — author residual epithelium (all-epithelial)
- E-MTAB-13526 — author epithelium only
- GSE267108 / GSE274595 — malignant split empty or n<4
- leftover `stats.tsv` — documents those n
