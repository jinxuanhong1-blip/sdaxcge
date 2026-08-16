# Extra — two-way public lung scRNA (timepoint × response)

Among public lung scRNA series, **GSE207422** is the only deposit that jointly labels timepoint (pre/post) and response (MPR/RECIST) on a matrix that includes malignant epithelial cells. Other public ICI lung scRNA sets are post-only (GSE291670, GSE243013, GSE146100), pre-only for the response core (GSE205335), T-cell-only (GSE179994), or CD45-sorted (GSE229353). No second series was stacked.

GSE207422 is cross-sectional, not paired: 3 pre-treatment biopsies and 12 post-treatment resections are different patients. The MPR 2×2 of malignant TACSTD2 / CLDN4 (patient mean log1p CP10k; pCR counted as MPR) is:

|  | MPR | NMPR |
|---|---|---|
| pre | empty (n=0) | TACSTD2 1.75 / CLDN4 1.23 (n=2) |
| post | TACSTD2 1.33 / CLDN4 1.55 (n=4) | TACSTD2 1.69 / CLDN4 1.52 (n=8) |

Post NMPR vs MPR: TACSTD2 Mann–Whitney p=0.21; CLDN4 p=0.93. NMPR post vs pre: TACSTD2 p=0.89; CLDN4 p=0.40. Pre-MPR is missing, so the interaction is not estimable. RECIST (PR vs SD) fills all four cells only if n=1 cells are kept; post SD vs PR is p=1.00 (TACSTD2) and p=0.93 (CLDN4).
