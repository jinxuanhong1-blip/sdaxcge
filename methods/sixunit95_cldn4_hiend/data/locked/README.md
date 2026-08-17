# Locked six-unit patient tables

Copied from `methods/strict_malig_tnk/data/` (PR #312). These are the
strict-malignant primary-grid members whose pooled CLDN4 vs T/NK Spearman
is taken as given: **n=95, ρ=−0.260**. This folder does **not** re-audit
that row.

| unit | file | malignant def | n |
|---|---|---|---:|
| GSE207422 | `GSE207422_drmref_patients.tsv` | author DRMref | 12 |
| GSE205335 | `GSE205335_patients.tsv` | author malignant | 22 |
| GSE291670 | `GSE291670_patients.tsv` | marker malignant | 6 |
| GSE253013 | `GSE253013_patients.tsv` | marker malignant-like, Tumor | 9 |
| GSE131907 | `GSE131907_samples.tsv` | author malignant, n_mal≥20 | 21 |
| GSE325414 | `GSE325414_donors.csv` | author malignant | 25 |

GSE253013’s only public matrix is the 9.3 GB RDS and is **not** downloaded.
