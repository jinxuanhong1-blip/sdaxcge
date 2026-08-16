# Extra atlas n · Harmony joint naive LUAD/NSCLC scRNA

Public treatment-naïve / diagnostic tumor scRNA. Not an ICI test.

Harmony joint object: GSE131907 tumor sites (Kim 2020) + GSE253013 9 LUAD
(Xiang/Sze 2024) + GSE148071 (Wu 2021) + GSE127465 tumor (Zilionis 2019).
435,519 tumor cells; Harmony on a 63-gene lineage panel (batch = dataset).
Leader/GSE154826 omitted (CD45+ CITE-seq; epithelium is doublet-gate leak).

Unit = tumor donor. Eligible: ≥20 epithelial and ≥20 T/NK cells.

| Contrast | n | ρ | p |
|---|---:|---:|---:|
| Epithelial TACSTD2 vs T/NK | 76 | +0.15 | 0.18 |
| Epithelial CLDN4 vs T/NK | 76 | +0.12 | 0.30 |
| Malignant-like TACSTD2 vs T/NK | 73 | −0.33 | 0.0038 |
| Malignant-like CLDN4 vs T/NK | 73 | −0.27 | 0.020 |
| Epithelial TACSTD2 vs B/TLS-like | 71 | −0.0004 | 1.00 |
| Epithelial CLDN4 vs B/TLS-like | 71 | +0.058 | 0.63 |

Malignant-like TACSTD2 vs T/NK remains inverse after dataset-residualized
ranks (n=73, ρ=−0.25, p=0.031). Per-atlas: GSE253013 n=9 ρ=−0.72 p=0.030;
GSE131907 n=31 ρ=−0.32 p=0.080; GSE148071 n=26 ρ=−0.15 p=0.45;
GSE127465 n=7 ρ=+0.43 p=0.34. B/TLS-like is a B/plasma + chemokine proxy,
not histology.
