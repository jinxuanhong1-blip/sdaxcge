# Paper extras — TACSTD2 / CLDN4 in leftover public lung IO RNA

Source: `results/rework/A2_wave2/extra/`. User A2 durvalumab result taken as given.

## TACSTD2 / CLDN4 vs native response

| Series | Gene | Endpoint | n_pos / n_neg | median_pos | median_neg | Cliff δ | p | After ESTIMATE residual p |
|---|---|---|---|---|---|---|---|---|
| GSE207422 | TACSTD2 | MPR | 9/15 | 4.999 | 6.222 | -0.244 | 0.34 | 0.952 |
| GSE207422 | CLDN4 | MPR | 9/15 | 5.854 | 6.272 | -0.289 | 0.257 | 0.858 |
| GSE207422 | TACSTD2 | ORR | 17/7 | 6.052 | 6.222 | 0.092 | 0.757 | 0.664 |
| GSE207422 | CLDN4 | ORR | 17/7 | 5.874 | 6.493 | -0.277 | 0.318 | 0.13 |
| GSE329813 | TACSTD2 | MPR | 11/11 | 2.247 | 2.444 | -0.686 | 0.0071 | 0.0878 |
| GSE329813 | CLDN4 | MPR | NE | — | — | — | NE | NE |
| GSE126044 | TACSTD2 | ORR | 5/11 | 6.168 | 6.288 | -0.273 | 0.441 | 0.913 |
| GSE126044 | CLDN4 | ORR | 5/11 | 2.542 | 3.771 | -0.527 | 0.115 | 0.661 |
| GSE166449 | TACSTD2 | ORR | 7/15 | 3.346 | 2.436 | 0.238 | 0.407 | 0.63 |
| GSE166449 | CLDN4 | ORR | 7/15 | 1.310 | 1.511 | 0.029 | 0.945 | 0.68 |
| GSE135222 | TACSTD2 | DCB | 7/20 | 7.341 | 7.865 | -0.143 | 0.607 | 0.314 |
| GSE135222 | CLDN4 | DCB | 7/20 | 7.689 | 7.126 | 0.114 | 0.685 | 0.498 |
| GSE248378 | TACSTD2 | recurrence | 9/20 | 6.340 | 5.424 | 0.456 | 0.0562 | 0.465 |
| GSE248378 | CLDN4 | recurrence | 9/20 | 6.292 | 5.321 | 0.244 | 0.311 | 0.588 |

Positive Cliff δ means the gene is **higher** in the positive class (MPR / responder / DCB / recurrence). Signs are not forced onto the user A2 coefficient.

## Continuous PFS (where deposited)

| Series | Gene | n / events | unadj HR | unadj p | ESTIMATE-residual HR | residual p |
|---|---|---|---|---|---|---|
| GSE135222 | TACSTD2 | 27 / 21 | 1.035 | 0.776 | 1.161 | 0.317 |
| GSE135222 | CLDN4 | 27 / 21 | 1.050 | 0.608 | 1.163 | 0.242 |
| GSE248378 | TACSTD2 | 29 / 9 | 2.236 | 0.016 | 1.675 | 0.174 |
| GSE248378 | CLDN4 | 29 / 9 | 1.410 | 0.17 | 1.216 | 0.43 |

## TACSTD2 / CLDN4 vs CD8A and GEP18 after purity residual

| Series | Gene | Immune | n | raw ρ | raw p | partial ρ given ESTIMATE | partial p |
|---|---|---|---|---|---|---|---|
| GSE207422 | TACSTD2 | CD8A | 24 | -0.501 | 0.0127 | -0.332 | 0.121 |
| GSE207422 | TACSTD2 | GEP18 | 24 | -0.524 | 0.00853 | -0.259 | 0.233 |
| GSE207422 | CLDN4 | CD8A | 24 | -0.361 | 0.0832 | -0.238 | 0.274 |
| GSE207422 | CLDN4 | GEP18 | 24 | -0.459 | 0.024 | -0.325 | 0.13 |
| GSE126044 | TACSTD2 | CD8A | 16 | -0.218 | 0.418 | -0.046 | 0.87 |
| GSE126044 | TACSTD2 | GEP18 | 16 | -0.050 | 0.854 | 0.247 | 0.375 |
| GSE126044 | CLDN4 | CD8A | 16 | -0.456 | 0.0759 | -0.125 | 0.657 |
| GSE126044 | CLDN4 | GEP18 | 16 | -0.335 | 0.204 | 0.094 | 0.738 |
| GSE166449 | TACSTD2 | CD8A | 22 | -0.078 | 0.73 | -0.151 | 0.515 |
| GSE166449 | TACSTD2 | GEP18 | 22 | 0.129 | 0.566 | 0.178 | 0.44 |
| GSE166449 | CLDN4 | CD8A | 22 | -0.246 | 0.269 | -0.104 | 0.653 |
| GSE166449 | CLDN4 | GEP18 | 22 | -0.118 | 0.601 | 0.095 | 0.681 |
| GSE135222 | TACSTD2 | CD8A | 27 | -0.012 | 0.954 | -0.267 | 0.187 |
| GSE135222 | TACSTD2 | GEP18 | 27 | 0.168 | 0.403 | -0.047 | 0.821 |
| GSE135222 | CLDN4 | CD8A | 27 | 0.123 | 0.542 | -0.156 | 0.448 |
| GSE135222 | CLDN4 | GEP18 | 27 | 0.280 | 0.157 | 0.041 | 0.843 |
| GSE248378 | TACSTD2 | CD8A | 29 | -0.706 | 1.89e-05 | -0.577 | 0.0013 |
| GSE248378 | TACSTD2 | GEP18 | 29 | -0.630 | 0.000253 | -0.404 | 0.033 |
| GSE248378 | CLDN4 | CD8A | 29 | -0.464 | 0.0112 | -0.332 | 0.0843 |
| GSE248378 | CLDN4 | GEP18 | 29 | -0.382 | 0.0407 | -0.172 | 0.383 |
| GSE329813 | TACSTD2 | CD8A | 22 | -0.604 | 0.00294 | -0.562 | 0.00797 |
| GSE329813 | TACSTD2 | GEP18 | 22 | -0.378 | 0.083 | -0.200 | 0.384 |
| GSE329813 | CLDN4 | CD8A | 0.0 | NE | NE | NE | NE |
| GSE329813 | CLDN4 | GEP18 | 0.0 | NE | NE | NE | NE |

Partial Spearman is the A2-matched estimator (Pearson on ranks, n−3 df). GSE248378 is the open post-durvalumab series; the other rows are leftover public anti-PD-1 / chemo-IO lung RNA.

## ESTIMATE / GEP coverage

| Series | n | ImmuneSignature overlap | StromalSignature overlap | GEP18 genes present | TumorPurity OOB |
|---|---|---|---|---|---|
| GSE207422 | 24 | 141 | 139 | 18/18 | 0 |
| GSE126044 | 16 | 141 | 137 | 18/18 | 6 |
| GSE166449 | 22 | 141 | 138 | 18/18 | 12 |
| GSE135222 | 27 | 140 | 135 | 18/18 | 0 |
| GSE329813 | 22 | 76 | 58 | 17/18 | 0 |
| GSE248378 | 29 | 138 | 133 | 18/18 | 6 |

GSE329813 is a GeoMx panel (~1.8k genes). ESTIMATE overlap is incomplete; the purity residual is still computed on the genes that are present and is labelled as such.

Figures: `extra_TACSTD2_response_boxplots.png`, `extra_response_forest_delta.png`,
`extra_TACSTD2_CD8_GEP_partial_forest.png`, `extra_CLDN4_CD8_GEP_partial_forest.png`.

