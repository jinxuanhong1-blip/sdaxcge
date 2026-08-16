# Harvested sample-level tables

These are **not** new GEO downloads except GSE241934. Sibling agents already computed per-sample TACSTD2/CLDN4 and T/NK on public matrices. This folder copies those tables so the histology split is additive and auditable.

| File | Source branch / extract | Histology field |
| --- | --- | --- |
| `gse207422_sample_table.tsv` | `cursor/nsclc-ici-tacstd2-cldn4-2443` | `Pathology` |
| `gse205335_patient_table.tsv` | same | `cancer_subtype` |
| `gse131907_sample_table.tsv` | same | series is LUAD |
| `gse253013_per_patient.tsv` | `cursor/gse253013-luad-io-scrna-4916` | series is LUAD |
| `gse291670_per_sample.tsv` | `cursor/a3-extra-gse291670-15cd` | none on GEO (kept for catalog) |
| `gse241934_per_sample.tsv` | `scripts/extract_gse241934.py` on GEO MTX | `Histology` |

Do not hand-edit ρ. Recompute with `scripts/01_compute_cohort_effects.py`.
