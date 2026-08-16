# TISCH2 NSCLC pool — epithelial TACSTD2/CLDN4 vs T/NK fraction

**Slice:** `methods/tisch_nsclc_pool/` only. Public TISCH2 h5 + CellMetainfo.
Numbers below are written from `tables/` and `summary.json`.

## Verdict (honest)

- **TACSTD2 mean vs T/NK (Fisher-z pool):** n_datasets=4, n_units=67, ρ=-0.016, p=0.908, 95% CI [-0.273, +0.244], I²=4%.
- **CLDN4 mean vs T/NK (Fisher-z pool):** n_datasets=4, n_units=67, ρ=+0.026, p=0.847, 95% CI [-0.234, +0.282], I²=0%.
- **ICI labels:** TISCH gallery marks GSE151537, GSE146100, and GSE176021_aPD1 as Immunotherapy. Downloadable CellMetainfo has **no Response/RECIST/MPR column**. GSE151537 and GSE176021 are T-sorted (no epithelium). GSE146100 is 1 patient / 3 nodules. This pool is **not** an ICI-response test.

## Per-dataset Spearman (eligible tumor-like units)

| Dataset | Unit | Contrast | n | ρ | p | note |
|---|---|---|---:|---:|---:|---|
| NSCLC_GSE117570 | sample | CLDN4_mean_vs_fracTNK | 2 | — | — | n=2<5; Spearman not computed |
| NSCLC_GSE117570 | sample | TACSTD2_mean_vs_fracTNK | 2 | — | — | n=2<5; Spearman not computed |
| NSCLC_GSE127465 | patient | CLDN4_mean_vs_fracTNK | 7 | +0.107 | 0.819 |  |
| NSCLC_GSE127465 | patient | TACSTD2_mean_vs_fracTNK | 7 | +0.321 | 0.482 |  |
| NSCLC_GSE131907 | sample | CLDN4_mean_vs_fracTNK | 27 | +0.038 | 0.851 |  |
| NSCLC_GSE131907 | sample | TACSTD2_mean_vs_fracTNK | 27 | +0.156 | 0.438 |  |
| NSCLC_GSE143423 | sample | CLDN4_mean_vs_fracTNK | 2 | — | — | n=2<5; Spearman not computed |
| NSCLC_GSE143423 | sample | TACSTD2_mean_vs_fracTNK | 2 | — | — | n=2<5; Spearman not computed |
| NSCLC_GSE146100 | sample | CLDN4_mean_vs_fracTNK | 3 | — | — | n=3<5; Spearman not computed |
| NSCLC_GSE146100 | sample | TACSTD2_mean_vs_fracTNK | 3 | — | — | n=3<5; Spearman not computed |
| NSCLC_GSE148071 | sample | CLDN4_mean_vs_fracTNK | 25 | +0.135 | 0.519 |  |
| NSCLC_GSE148071 | sample | TACSTD2_mean_vs_fracTNK | 25 | -0.138 | 0.512 |  |
| NSCLC_GSE149655 | sample | CLDN4_mean_vs_fracTNK | 2 | — | — | n=2<5; Spearman not computed |
| NSCLC_GSE149655 | sample | TACSTD2_mean_vs_fracTNK | 2 | — | — | n=2<5; Spearman not computed |
| NSCLC_GSE150660 | sample | CLDN4_mean_vs_fracTNK | 2 | — | — | n=2<5; Spearman not computed |
| NSCLC_GSE150660 | sample | TACSTD2_mean_vs_fracTNK | 2 | — | — | n=2<5; Spearman not computed |
| NSCLC_GSE153935 | sample | CLDN4_mean_vs_fracTNK | 5 | +0.200 | 0.747 |  |
| NSCLC_GSE153935 | sample | TACSTD2_mean_vs_fracTNK | 5 | -0.600 | 0.285 |  |
| NSCLC_GSE162498 | sample | CLDN4_mean_vs_fracTNK | 8 | -0.524 | 0.183 |  |
| NSCLC_GSE162498 | sample | TACSTD2_mean_vs_fracTNK | 8 | -0.524 | 0.183 |  |

## ICI / skip audit

| Dataset | Gallery | CellMetainfo ICI cols | Eligible units | In pool | Skip / ICI note |
|---|---|---|---:|---|---|
| NSCLC_GSE151537 | Immunotherapy | none | 0 | no | T-sorted only (no epithelial/malignant cells) |
| NSCLC_GSE146100 | Immunotherapy | none | 3 | no | eligible samples n=3 (need >=5 for Spearman / pool) |
| NSCLC_GSE117570 | None | none | 2 | no | eligible samples n=2 (need >=5 for Spearman / pool) |
| NSCLC_GSE127465 | None | none | 7 | yes | no ICI labels in TISCH gallery or CellMetainfo |
| NSCLC_GSE131907 | None | none | 27 | yes | no ICI labels in TISCH gallery or CellMetainfo |
| NSCLC_GSE148071 | None | none | 25 | yes | no ICI labels in TISCH gallery or CellMetainfo |
| NSCLC_EMTAB6149 | None | none | 0 | no | CellMetainfo has no Patient/Sample column |
| NSCLC_GSE127471 | None | none | 0 | no | PBMC only; 0 epithelial/malignant cells |
| NSCLC_GSE139555 | None | none | 0 | no | T/immune-sorted; 0 epithelial/malignant cells |
| NSCLC_GSE143423 | None | none | 2 | no | eligible samples n=2 (need >=5 for Spearman / pool) |
| NSCLC_GSE99254 | None | none | 0 | no | T-sorted; 0 epithelial/malignant cells |
| NSCLC_GSE149655 | None | none | 2 | no | eligible samples n=2 (need >=5 for Spearman / pool) |
| NSCLC_GSE150660 | None | none | 2 | no | eligible samples n=2 (need >=5 for Spearman / pool) |
| NSCLC_GSE153935 | None | none | 5 | yes | no ICI labels in TISCH gallery or CellMetainfo |
| NSCLC_GSE162498 | None | none | 8 | yes | no ICI labels in TISCH gallery or CellMetainfo |
| NSCLC_GSE176021_aPD1 | Immunotherapy | none | 0 | no | T-sorted (~817k cells); 0 epithelial/malignant; CellMetainfo has no Response column |
| NSCLC_GSE179373 | None | Treatment | 0 | no | T-sorted; 0 epithelial/malignant cells |

## Coverage

- Catalog datasets: 17
- Datasets with Spearman n≥5: 5
- Datasets in Fisher-z (n≥6): 4 (GSE153935 n=5 is Spearman-only)
- Eligible sample/patient units written: 83
- Generated: 2026-08-16T21:18:38.640023+00:00

## Caveats

- TISCH values are MAESTRO `log2(TPM/10+1)`, not raw UMI.
- GSE131907 has no TISCH Malignant call; tumor-tissue epithelial cells are a proxy.
- GSE131907 TISCH `Source` is swapped for several LUNG_T/LUNG_N pairs; author sample names are used.
- T/NK fraction is composition after dissociation and TISCH annotation, not spatial infiltration.
- Sample-level tests are the unit. Cell-level p-values are not reported as primary.
- GSE146100 nodule-level n=3 is listed in `per_unit_metrics.tsv` and excluded from Spearman.
- Fisher-z is the primary pool. Percentile-rank-within-dataset is secondary.

