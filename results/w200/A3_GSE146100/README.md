# W200-A3 · GSE146100 — malignant TACSTD2, R vs NR

**Task:** A3 analog on the only public multi-nodule pembrolizumab LUAD scRNA-seq
set: compare TACSTD2 in malignant cells between the responding nodule and the
two non-responding nodules.

## Verdict

**No patient-level R-vs-NR claim is possible.** This series is one 72-year-old
patient (Zhang et al., *JITC* 2021, PMID 33820821). The responding nodule is
W2 (KRAS G12C). The two non-responding nodules are W1 and W3 (both EGFR).
Response, driver genotype, and clone identity are the same contrast.

Pooled cluster-epithelial TACSTD2 looks **higher in NR**
(mean log-norm R=1.282 vs NR=1.624;
% expressing R=78.0 vs NR=80.2;
linear-scale log2FC(R/NR)=-0.601; cell-level MWU p=1.24e-07,
n=236 R / 847 NR cells). **Do not use that p-value.**
Cells are nested in 1 vs 2 nodules. Per-nodule means show the pooled NR
elevation is **W3-driven**: W1 (NR) ≈ W2 (R); W3 (NR) is higher and holds
most NR epithelial cells. A nodule-level test cannot be run (n=1 vs 2).

The marker-rule epithelial sensitivity set also looks higher in NR in the
pool (log2FC=-0.804)
and is likewise W3-driven (W1 is not higher than W2).

## What was used

- GEO [GSE146100](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146100),
  `GSE146100_NormData.txt.gz` (229 MB; author Seurat log-normalized matrix).
- 11612 cells × 19800 genes; barcodes prefixed
  W1_/W2_/W3_.
- No author cell-type labels. Malignant = Leiden clusters whose highest lineage
  score is epithelial (EPCAM/KRT/SFT* vs PTPRC/immune/stroma). Sensitivity:
  any of EPCAM/KRT8/KRT18/KRT19 > 0 and PTPRC = 0.
- Raw SRA FASTQ was not downloaded.

## Honest limits vs a real A3 / ICI claim

- n_patients = 1. This is a within-patient nodule comparison, not a cohort.
- R = KRAS; both NR = EGFR. Any TACSTD2 difference may be genotype, not response.
- All three nodules are post-pembrolizumab. There is no pre-treatment arm.
- "Malignant" is an epithelial proxy. Without CNV or a matched normal we cannot
  exclude residual normal epithelial cells.
- Cell-level statistics are pseudoreplicated.
- The pooled “NR higher” mean is not a shared NR property: W1 ≈ W2; W3 is the high nodule.

## Results (cluster epithelial, primary)

| Nodule | Response | Genotype | n malignant | % TACSTD2+ | mean log-norm |
|---|---|---|---|---|---|
| W1 | NR | EGFR L858R | 124 | 75.0 | 1.290 |
| W2 | R | KRAS G12C | 236 | 78.0 | 1.282 |
| W3 | NR | EGFR L858R/R77H | 723 | 81.1 | 1.681 |

Pooled R (W2) vs NR (W1+W3): log2FC=-0.601, cell-level AUC=0.388.

TACSTD2 is epithelial-restricted in this matrix (sanity): epithelial 79.7%
positive / mean 1.55 vs T/NK 2.2% / mean 0.03.

## Files

| File | Role |
|---|---|
| `sample_metadata.tsv` | Nodule / GSM / genotype / cell counts |
| `sample_composition.tsv` | Lineage mix per nodule |
| `cluster_lineage_summary.tsv` | Leiden → lineage + marker means |
| `cell_annotation.tsv` | Per-cell labels and TACSTD2 |
| `tacstd2_malignant_by_nodule.tsv` | Per-nodule and pooled R/NR stats |
| `tacstd2_R_vs_NR_celllevel.tsv` | Cell-level contrast (flagged) |
| `tacstd2_by_lineage.tsv` | Compartment sanity |
| `summary.json` / `sanity_checks.json` / `audit.json` | Verdict |
| `fig_umap_lineage_sample.png` / `fig_umap_markers.png` | UMAP |
| `fig_tacstd2_malignant_violin.png` / `fig_tacstd2_malignant_bar.png` | TACSTD2 |

## Reproduce

```bash
python3 scripts/w200/A3_GSE146100/download.py
python3 scripts/w200/A3_GSE146100/analyze.py
```
