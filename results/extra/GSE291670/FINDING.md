# EXTRA — GSE291670 malignant TACSTD2 / CLDN4 (not GSE207422)

Independent second public neoadjuvant ICI lung scRNA series. The GSE207422 CopyKAT slide is taken as given and is not re-audited here.

**Series:** Xia et al., *J Transl Med* 2025. GEO [GSE291670](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291670). Six post-treatment NSCLC tumors after neoadjuvant **anlotinib + camrelizumab**. GEO titles = MPR-1/2/3 vs Non-MPR-1/2/3.

**Why not the named alternatives:** GSE243013 (closest to “GSE253013”) is CD45+ immune-only (0 malignant cells). GSE266035 is blood T cells from one patient. GSE131907 tLung is treatment-naive with no MPR.

## n

| | Count |
|---|---:|
| Patients | 6 (MPR 3, NMPR 3) |
| QC cells | 28,847 |
| Epithelial | 8,448 |
| Malignant (epithelial, not high normal-lung score) | 6,335 |
| Lineage T/NK | 1,755 |

Author barcode labels are not on GEO. Malignant = marker epithelium minus a high normal-lung program. All six samples have ≥353 malignant cells.

## Extra figure (`fig_extra_malignant_tacstd2_cldn4.png`)

Patient is the unit. Exact Wilcoxon enumerates C(6,3)=20; the smallest two-sided p at 3 vs 3 is 0.10.

| Panel | Test | n | Result | p |
|---|---|---|---|---|
| A | malignant TACSTD2 mean log1p(CP10k), NMPR vs MPR | 3 vs 3 | NMPR 0.099 vs MPR 0.226; **Δ=−0.13** (MPR higher) | **0.20** |
| B | malignant CLDN4 mean log1p(CP10k), NMPR vs MPR | 3 vs 3 | NMPR 0.121 vs MPR 0.209; **Δ=−0.09** (MPR higher) | **0.40** |
| C | malignant TACSTD2 vs lineage T/NK fraction | 6 | **ρ=−0.54** | **0.27** |
| D | malignant CLDN4 vs lineage T/NK fraction | 6 | **ρ=−0.83** | **0.042** |

## Other honest numbers

| Test | n | ρ or Δ | p |
|---|---|---|---|
| TACSTD2 %pos, NMPR vs MPR | 3 vs 3 | Δ=−0.038 | 0.70 |
| CLDN4 %pos, NMPR vs MPR | 3 vs 3 | Δ=−0.012 | 0.70 |
| TACSTD2 mean vs CD3E/CD8A/NKG7 UMI fraction | 6 | ρ=−0.37 | 0.47 |
| TACSTD2 %pos vs lineage T/NK | 6 | ρ=−0.43 | 0.40 |
| CLDN4 mean vs UMI T/NK | 6 | ρ=−0.54 | 0.27 |
| CLDN4 %pos vs lineage T/NK | 6 | ρ=−0.89 | 0.019 |
| all-epithelial TACSTD2 vs MPR | 3 vs 3 | Δ=−0.19 (MPR higher) | 0.20 |

## Read this as extra, not as a copy of the given slide

- **MPR contrast:** malignant TACSTD2 and CLDN4 are **higher in MPR**, not NMPR. Both p≥0.20. At n=3 vs 3 this cannot be called significant.
- **vs T/NK:** both genes trend **negative**. TACSTD2 ρ=−0.54 is NS. CLDN4 ρ=−0.83 / %pos ρ=−0.89 have nominal p<0.05; **n=6 Spearman is one-or-two-point fragile** and is not a confirmatory cohort.
- Cell-level p-values are not reported (pseudoreplication).

## Files

- `fig_extra_malignant_tacstd2_cldn4.png` — the extra 4-panel figure
- `per_sample.tsv` / `stats.tsv` / `lineage_counts.tsv` / `summary.json`
- `cell_calls.tsv.gz` — QC barcodes with lineage flags
- `METHODS.md`
