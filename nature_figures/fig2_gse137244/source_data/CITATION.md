# Source data — Fig. 2 GSE137244

Bulk RNA-seq of mouse lung-cancer **cell lines**. Not single-cell RNA-seq. Not the private 8KL atlas. The normal-lung library (GSM4073826) is not plotted.

## Study

- GEO series: [GSE137244](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137244)
- Deng et al., PMID [34142094](https://pubmed.ncbi.nlm.nih.gov/34142094/)
- Matrix named in the upstream analysis: `GSE137244_counts.fpkm.csv.gz`
- Transform: log2(FPKM + 1)
- Contrast: mean of KL libraries minus mean of KP libraries
- Test plotted: two-sided exact Mann–Whitney U on the 10 tumor libraries

## Libraries (n = 5 vs 5)

| Arm | Genotype | Libraries | GSM |
|---|---|---|---|
| KP | KrasG12D;Trp53 | B6AL10-1, B6AL10-2, B6AL10-3, B6AL10-4, B6AL10-5 | GSM4073816–GSM4073820 |
| KL | KrasG12D;Lkb1 | KL155, KL47, KLC, KLD, KLE | GSM4073821–GSM4073825 |

The five KP columns are five libraries of the B6AL10 series, not five independent KP lines. The Mann–Whitney unit is the library, which is the unit that gives P = 2/252 = 0.00794 when the two groups do not overlap.

## Files in this folder

Copied verbatim from branch `cursor/gse137244-kl-kp-gsea-3033` (PR 606), paths under `methods/gse137244_kl_kp_deep/`:

| File | Upstream path | Role |
|---|---|---|
| `watch_genes_log2.tsv` | `tables/watch_genes_log2.tsv` | Per-library log2(FPKM+1) for Tacstd2, claudins, Ocln, Tjp1. These are the only gene-level points drawn. |
| `mean_signature_scores.tsv` | `tables/mean_signature_scores.tsv` | Row `EPITHELIAL_TJ_CORE`: per-library mean of the 18-gene set. Panel c uses this row and does not recompute it. |
| `samples.tsv` | `tables/samples.tsv` | Library, short name, GSM, arm, genotype |
| `tj_gene_callouts.tsv` | `tables/tj_gene_callouts.tsv` | Published gene-level Δ and Mann–Whitney results, used as a check |
| `locked_genes.tsv` | `tables/locked_genes.tsv` | Locked Tacstd2 and Cldn4 Δ and P |
| `epithelial_tj_core_genes.txt` | `gene_sets/epithelial_and_isg.gmt` (`EPITHELIAL_TJ_CORE`) | The 18 symbols in the panel c score |

Locked summaries also reported on branch `cursor/gse137244-kl-kp-genesets-385b` (PR 617): Tacstd2 Δ = +3.24, Cldn4 Δ = +5.57, two-sided Mann–Whitney P = 0.00794, n = 5 vs 5, log2(FPKM+1).

`fig2_plotted_values.tsv` and `fig2_summary.tsv` are written by `plot_fig2_gse137244.py` from the tables above. No sample value in the figure is typed in by hand.
