# Results · Milo neighborhood DA (GSE207422)

Public-only run of `methods/scrna_milo/` on Hu et al. 2023 (GEO GSE207422).

**Verdict:** 0 neighborhoods at SpatialFDR < 0.1 or < 0.2 on both embeddings and both contrasts (TACSTD2 high vs low, n=5 vs 5; MPR vs NMPR, n=4 vs 8). Compositional T/NK-depleted / TACSTD2-high-epi interface nhoods exist (24 per embedding) and are listed; they are not DA-supported. Details in `FINDING.md`.

## Files

| File | Role |
|---|---|
| `FINDING.md` | Honest SpatialFDR / n |
| `summary.json` | Machine-readable n, lineage totals, verdict |
| `da_counts.tsv` | Per-embedding × contrast SpatialFDR counts |
| `sample_table.tsv` | Patient malignant TACSTD2 split + MPR |
| `lineage_counts_by_sample.tsv` | Marker lineages in the epi+immune graph |
| `nhoods_{pca,harmony}_{tacstd2_high_vs_low,mpr_vs_nmpr}.tsv` | All neighborhoods + DA |
| `kept_*.tsv` | Interface, TACSTD2-high epi, T/NK-depleted |
| `kept_index_cells_by_embedding.tsv` | Index-cell overlap across PCA/Harmony |
| `fig_volcano_*.png` | logFC vs −log10 SpatialFDR |
| `fig_interface_*.png` | epi TACSTD2 vs T/NK fraction |
| `fig_pvalue_histograms.png` | Nominal P (shows the test is not degenerate) |

## Reproduce

```bash
python3 methods/scrna_milo/00_download.py
python3 methods/scrna_milo/01_stream_matrix.py
python3 methods/scrna_milo/02_run_milo.py
```
