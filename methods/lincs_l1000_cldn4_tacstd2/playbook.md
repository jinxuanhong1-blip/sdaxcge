# Playbook — LINCS L1000 CLDN4 shRNA and TACSTD2 CRISPR (IFN / APM)

Additive public evidence. Does not use private knockdown coculture matrices. Does not revise locked CosMx, concordant-4, GSE137244, TCGA, or TISMO numbers.

```bash
python3 methods/lincs_l1000_cldn4_tacstd2/download.py
python3 methods/lincs_l1000_cldn4_tacstd2/analyze.py
```

`download.py` needs network. `analyze.py` reads the committed slice only.

| Locked choice | Value |
|---|---|
| Catalogs | LINCS 2020 `siginfo_beta.txt`; GSE92742 sig + pert info; GSE70138 sig info (2017-03-06) |
| CLDN4 | Any `cmap_name` / `pert_iname` hit. If zero, stop. No enrichment test. |
| TACSTD2 shRNA | `trt_sh`, `trt_sh.cgs`, `trt_sh.css`. If zero, say so. |
| TACSTD2 matrix | Level 5 moderated z, `trt_xpr` only, from the CRISPR gctx via HTTP range reads |
| Signatures | All 20 `qc_pass` exemplar signatures (2 guides × 10 cell lines, 96 h) |
| Primary rank | Within each cell line, mean z of its two guides; then median across the 10 lines |
| Gene dropped from the rank | TACSTD2 (the CRISPR target). CLDN4 stays in the background. |
| Duplicate symbol | MIA2 is both best inferred (Entrez 4253) and inferred (117153). Keep best inferred. |
| Engine | `gsea_core.py` weighted KS *p*=1, 1000 gene-set permutations, seed=42 |
| Headline sets | Hallmark IFN-γ, custom MHC-I / APM (21) |
| Secondary | Hallmark IFN-α, reported, outside the BH family |
| FDR | BH inside the two headline sets, separately on each consensus rank |
| Positive NES | Set is enriched at the high-z end (up after CRISPR) |
| Aliases | `WARS1→WARS`, `MARCHF1→MARCH1` only, and only when the current symbol is absent |
| Sensitivity | Median of 20 signatures; cell-median of signatures with `cc_q75 ≥ 0.2`; cell-median of `is_hiq==1`; landmark genes only on the primary z |
| Concordance | Per-signature NES and the mean of the two guides per cell. Wilcoxon vs 0 is descriptive. |
| On-target | TACSTD2 z and CLDN4 z are reported. They are not gene-set tests. |

Outputs: `FINDING.md`, `tables/gsea_headline.tsv`, `figures/fig_headline_nes.png`.
