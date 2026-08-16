# Manifest — results/w200/C4_GSE334497

| Path | What |
|---|---|
| `WRITEUP.md` | Honest C4 analog call (zh+en) |
| `summary.json` | Machine-readable verdict + per-gene / per-set numbers |
| `tables/focal_genes.csv` | Tacstd2, Cldn4, Cldn7, Cxcl9, CORE6, T-cell genes |
| `tables/geneset_stats.csv` | Competitive MWU + exact permutation set tests |
| `tables/sample_set_scores.csv` | Per-sample mean z-scores |
| `tables/key_genes_log2.csv` | Per-sample log2(norm+1) for key genes |
| `tables/de_all_genes.csv` | Genome-wide Welch / MWU (no gene passes BH-FDR 0.05) |
| `figures/focal_boxplots.png` | Tacstd2, Cldn4, Cldn7, Cxcl9, Isg15, Cd8a |
| `figures/ifn_apm_heatmap.png` | Row-z heatmap, WT left / KO right |
| `figures/set_score_boxplots.png` | Pre-specified set scores |
| `raw/GSE334497_normalized_counts.csv.gz` | GEO supplementary (author-normalized counts) |
| `raw/ensembl_to_symbol.tsv` | NCBI GRCm38 Ensembl → official symbol (compact) |

Provenance: GSE334497, GPL21103, 10 samples (5 Trop2 KO + 5 WT 4T1 tumors, BALB/c, 3 weeks). Counts MD5 is in `summary.json`.
