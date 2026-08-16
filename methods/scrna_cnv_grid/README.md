# scRNA malignant-definition grid (`scrna_cnv_grid`)

Additive public-only grid of **how malignant cells are called**, then the same two sample-level scores the user uses for A3:

1. TACSTD2 in the called-malignant compartment vs **T/NK fraction** (Spearman ρ)
2. TACSTD2 in that compartment, **NMPR vs MPR** (exact Wilcoxon; pCR counted as MPR)

A3 (Hu et al. GSE207422; NMPR>MPR and ρ negative vs T/NK) is taken as given. This slice does not re-litigate that claim. It reports **which public definitions recover that direction**, with honest n / ρ / p.

## Series

| Accession | Design | Author malignant? |
|---|---|---|
| **GSE207422** | Neoadjuvant chemo-IO NSCLC, 12 post-treatment samples (MPR n=4 incl. pCR; NMPR n=8) | **No.** Author CopyKAT IDs are not on GEO. DRMref `Malignant cells` is a third-party label. |
| **GSE241934** | NEOTIDE IIT (EGFR-mut, 11 pts) + real-world EGFR-WT (34 pts) | **Yes.** Deposited `major.cell.type==Epi` (IIT n=1,699; paper CopyKAT malignant N=1,669). |

## Definitions (every row is reported)

- **Author** (if deposited): GSE241934 `Epi`. GSE207422: none.
- **Third-party** (GSE207422 only): DRMref malignant.
- **Marker epithelium:** argmax lineage on EPCAM/KRT8/18/19/7/CDH1; optional normal-lung filter (zero UMI or below epithelial p75 of SFTPA/AGER/SCGB/TPPP3/FOXJ1).
- **CNV aneuploid (CopyKAT/inferCNV-style):** window-smoothed expression CNV vs stromal reference, win=25; epithelial AND score > reference p95 (also p90 / 2-means).
- **Intersections:** author ∩ marker, author ∩ CNV-p95, marker-NNL ∩ CNV-p95, triple.

Unit of inference is the **patient/sample**. Cell-level p-values are not reported.

## Run

```bash
python3 methods/scrna_cnv_grid/download.py
python3 methods/scrna_cnv_grid/analyze_gse207422.py
python3 methods/scrna_cnv_grid/analyze_gse241934.py
python3 methods/scrna_cnv_grid/summarize.py
```

Outputs: `results/scrna_cnv_grid/`.
