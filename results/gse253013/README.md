# GSE253013 · extra LUAD scRNA (TACSTD2 / CLDN4 vs T/NK)

Additive figure only. **GSE207422 A3 is not re-analyzed.**

| Item | Value |
|---|---|
| Accession | [GSE253013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253013) |
| Paper | Sze/Xiang, *Cancer Res* 2024, PMID 38335304 |
| Tissue | Lung tumor + adjacent non-tumor (not blood) |
| Patients / cells / GSM | 9 treatment-naïve LUAD / 256,379 / 89 |
| ICI / MPR / R labels | **None public** (`treatment` = None) |
| TACSTD2 / CLDN4 | Both present |
| Analysis unit | **Patient** (not cell, not 10x lane) |

## Primary (tumor, marker malignant-like)

- TACSTD2 mean log1p vs T/NK: **n=9, ρ=−0.72, p=0.030**
- TACSTD2 log1p(CP10k) vs T/NK: n=9, ρ=−0.78, p=0.013
- TACSTD2 %pos vs T/NK: n=9, ρ=−0.83, p=0.0053
- All-epithelial TACSTD2 vs T/NK: n=9, ρ=−0.12, p=0.77
- CLDN4 mean log1p vs T/NK: **n=9, ρ=−0.33, p=0.38**
- Response contrast: **not possible**

## Caveats

- n=9 is small. Three tumors have <50 malignant-like cells (MRC004=16, MRC007=18, MRC006=30). n_malig≥50: TACSTD2 ρ=−0.60, p=0.28, n=5.
- Author protocol sorted CD45+/CD45− in 6/9 tumors; T/NK fraction is sampling-biased.
- Malignant-like is a marker proxy, not CopyKAT.
- Author `Epithelial` TACSTD2 vs T-cell fraction: ρ=−0.58, p=0.099, n=9.

Paper text: `paper/extra_gse253013_luad_scrna.md`. Scripts: `methods/`.
