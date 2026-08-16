# Results · naive LUAD/NSCLC Harmony joint object

Extra atlas *n*. Not an ICI / MPR / response test.

Joint tumor cells = 435,519. Harmony subsample = 187,168 cells (≤2,500 / donor),
63-gene panel, 20 PCs, batch = dataset.

Eligible tumor donors (≥20 epithelial and ≥20 T/NK): **n = 76**
(GSE131907 sites 34, GSE253013 9, GSE148071 26, GSE127465 7).
B/TLS-like contrasts use donors with ≥10 B cells (**n = 71**).

## Primary joint Spearman (donor unit)

| Contrast | n | ρ | p | BH-FDR |
|---|---:|---:|---:|---:|
| Epithelial TACSTD2 vs T/NK | 76 | +0.15 | 0.18 | 0.61 |
| Epithelial CLDN4 vs T/NK | 76 | +0.12 | 0.30 | 0.76 |
| Malignant-like TACSTD2 vs T/NK | 73 | −0.33 | 0.0038 | 0.038 |
| Malignant-like CLDN4 vs T/NK | 73 | −0.27 | 0.020 | 0.099 |
| Epithelial TACSTD2 vs B/TLS-like | 71 | −0.0004 | 1.00 | 1.00 |
| Epithelial CLDN4 vs B/TLS-like | 71 | +0.058 | 0.63 | 0.95 |
| Malignant-like TACSTD2 vs B/TLS-like | 68 | −0.040 | 0.75 | 0.95 |
| Malignant-like CLDN4 vs B/TLS-like | 68 | +0.023 | 0.85 | 0.95 |

Dataset-residualized ranks (malignant-like vs T/NK): TACSTD2 n=73 ρ=−0.25
p=0.031; CLDN4 n=73 ρ=−0.26 p=0.028.

Malignant-like TACSTD2 vs T/NK by atlas: GSE253013 n=9 ρ=−0.72 p=0.030;
GSE131907 n=31 ρ=−0.32 p=0.080; GSE148071 n=26 ρ=−0.15 p=0.45;
GSE127465 n=7 ρ=+0.43 p=0.34.

All-epithelial (not malignant-like) contrasts vs T/NK and all B/TLS-like
contrasts are near zero in the joint object.

Leader/GSE154826 was not merged (CD45+ CITE-seq; epithelium is doublet-gate
leak). Full table: `association_statistics.tsv`.
