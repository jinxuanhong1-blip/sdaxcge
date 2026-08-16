# Extra figure — T/NK cytotoxicity and exhaustion vs malignant TACSTD2 / CLDN4

Public lung ICI / neoadjuvant scRNA with both a malignant (or author epithelial) compartment and T/NK.
Patient is the unit. Scores are the mean of log1p(CP10k) over the named genes in the named compartment.

- Cytotoxicity in T/NK: GZMB, PRF1, GNLY, NKG7
- Exhaustion in T/NK: PDCD1, HAVCR2, LAG3, TIGIT, TOX
- Malignant TACSTD2 and CLDN4: mean log1p(CP10k) in malignant / author-Epi cells
- Keep a patient if ≥10 malignant and ≥10 T/NK cells
- Spearman on the patient table; Fisher-z IVW fixed-effect and DerSimonian–Laird random-effect meta (n≥5 / cohort)

## Cohorts

| Cohort | Setting | Labels | Patients in Spearman |
|---|---|---|---|
| GSE205335 | palliative ICI atlas (Park/Ahn/Lee); author malignant + T/NK | see METHODS.md | 22 |
| GSE207422 | neoadjuvant PD-1 + chemo (Hu 2023); post-tx; marker malignant-like | see METHODS.md | 7 |
| GSE233203 | leftover: pre-ABCP pleural effusion; marker compartments | see METHODS.md | 6 |
| GSE241934_IIT | NEOTIDE EGFR-mut neoadjuvant sintilimab + chemo; author Epi + T/NK | see METHODS.md | 11 |
| GSE241934_RWC | real-world neoadjuvant PD-1 + chemo; author Epi + T/NK | see METHODS.md | 29 |
| GSE291670 | neoadjuvant anlotinib + camrelizumab; marker malignant-like | see METHODS.md | 6 |

## Patient-level Spearman (honest n / ρ / p)

| Contrast | Cohort | n | ρ | p |
|---|---|---:|---:|---:|
| CLDN4 vs T/NK cytotoxicity | GSE205335 | 22 | -0.246 | 0.271 |
| CLDN4 vs T/NK cytotoxicity | GSE207422 | 7 | 0.000 | 1 |
| CLDN4 vs T/NK cytotoxicity | GSE233203 | 6 | 0.714 | 0.111 |
| CLDN4 vs T/NK cytotoxicity | GSE241934_IIT | 11 | -0.245 | 0.467 |
| CLDN4 vs T/NK cytotoxicity | GSE241934_RWC | 29 | 0.153 | 0.428 |
| CLDN4 vs T/NK cytotoxicity | GSE291670 | 6 | 0.771 | 0.0724 |
| CLDN4 vs T/NK exhaustion | GSE205335 | 22 | -0.193 | 0.391 |
| CLDN4 vs T/NK exhaustion | GSE207422 | 7 | 0.000 | 1 |
| CLDN4 vs T/NK exhaustion | GSE233203 | 6 | -0.257 | 0.623 |
| CLDN4 vs T/NK exhaustion | GSE241934_IIT | 11 | -0.536 | 0.089 |
| CLDN4 vs T/NK exhaustion | GSE241934_RWC | 29 | 0.018 | 0.927 |
| CLDN4 vs T/NK exhaustion | GSE291670 | 6 | -0.086 | 0.872 |
| TACSTD2 vs T/NK cytotoxicity | GSE205335 | 22 | 0.378 | 0.083 |
| TACSTD2 vs T/NK cytotoxicity | GSE207422 | 7 | -0.179 | 0.702 |
| TACSTD2 vs T/NK cytotoxicity | GSE233203 | 6 | 0.371 | 0.468 |
| TACSTD2 vs T/NK cytotoxicity | GSE241934_IIT | 11 | -0.382 | 0.247 |
| TACSTD2 vs T/NK cytotoxicity | GSE241934_RWC | 29 | -0.254 | 0.184 |
| TACSTD2 vs T/NK cytotoxicity | GSE291670 | 6 | 0.943 | 0.0048 |
| TACSTD2 vs T/NK exhaustion | GSE205335 | 22 | 0.228 | 0.308 |
| TACSTD2 vs T/NK exhaustion | GSE207422 | 7 | -0.143 | 0.76 |
| TACSTD2 vs T/NK exhaustion | GSE233203 | 6 | -0.657 | 0.156 |
| TACSTD2 vs T/NK exhaustion | GSE241934_IIT | 11 | -0.536 | 0.089 |
| TACSTD2 vs T/NK exhaustion | GSE241934_RWC | 29 | -0.314 | 0.0974 |
| TACSTD2 vs T/NK exhaustion | GSE291670 | 6 | 0.314 | 0.544 |

## Random-effect meta (Fisher z, n≥5)

| Contrast | k | N | ρ_RE | 95% CI | p_RE | I² |
|---|---:|---:|---:|---|---:|---:|
| TACSTD2 vs T/NK cytotoxicity | 6 | 81 | 0.176 | -0.313 to 0.591 | 0.488 | 68% |
| TACSTD2 vs T/NK exhaustion | 6 | 81 | -0.182 | -0.468 to 0.139 | 0.266 | 29% |
| CLDN4 vs T/NK cytotoxicity | 6 | 81 | 0.088 | -0.250 to 0.406 | 0.616 | 35% |
| CLDN4 vs T/NK exhaustion | 6 | 81 | -0.143 | -0.372 to 0.102 | 0.252 | 0% |

## Notes on n

- GSE207422: 12 post-tx samples; **7** have ≥10 malignant-like cells (P02/P06/P11/P13/P14 drop).
- GSE241934 RWC: 34 tumors; **29** have ≥10 author-Epi cells.
- GSE233203 leftover: 7 PE samples; **6** have ≥10 T/NK (NCCLu_397 has 3 T/NK).
- GSE291670 TACSTD2 vs cytotoxicity is the only cohort-level p<0.05 (**n=6, ρ=0.943, p=0.0048**). That point is why I²=68% on that contrast; dropping it gives RE ρ=−0.03, p=0.88 (N=75).
- Named-only (drop leftover): TACSTD2–cyto RE ρ=0.16, p=0.58, N=75. Named neoadjuvant-only TACSTD2–exhaustion RE ρ=−0.30, p=0.045, N=53, I²=0% (post-hoc slice; not the primary).

GSE146100 is a leftover with both compartments but n=1 patient and is not in the meta.
GSE243013 / GSE266035 / T-sorted series lack a malignant compartment and are not scored.
