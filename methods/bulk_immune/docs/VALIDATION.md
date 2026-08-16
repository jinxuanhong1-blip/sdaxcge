# Numerical validation

## Unit tests (no network)

```bash
cd methods/bulk_immune
PYTHONPATH=. python3 tests/test_ssgsea.py   # or: pytest tests/test_ssgsea.py
```

Covers rank-tie handling, within-sample monotone invariance of ssGSEA, CPM
scale detection, duplicate-gene collapse, Ensembl version stripping, and
family-wise BH-FDR.

## Optional: Python ssGSEA vs Bioconductor GSVA

Requires `GSVA` in R. Skips cleanly if it is not installed.

```bash
cd methods/bulk_immune/tests
Rscript test_ssgsea_vs_gsva.R
```

Agreement criterion: `max|python − GSVA| < 1e-4` on a 80-gene × 6-sample
Gaussian matrix, `ssgsea.norm=FALSE`, `alpha/tau=0.25`,
`tie_method=average_int` (GSVA truncates average ranks).

## Method-level notes

| Method | What is validated here | What is not |
|---|---|---|
| ssGSEA | statistic vs GSVA (optional); unit tests | GSVA's `ssgsea.norm=TRUE` global-range scaling (implemented, not the default) |
| ESTIMATE | same ssGSEA integral + official SI gene sets from estimate 1.0.13 | Affymetrix purity transform on RNA-seq (emitted only on request) |
| MCP-counter | official `genes.txt`; mean of log2 markers | nothing to fit |
| xCell | official 489 signatures, `fv`, `K` extracted from `xCell.data`; NNLS spillover | bit-identity with `GSVA::gsva` + `pracma::lsqlincon` on every cohort |
| TIDE | official `tidepy` 1.3.x weights | the web server's UI-only options |
| TIP | published annotation table; ssGSEA(pos)−ssGSEA(neg) | unpublished server-side weights, if any |
| CIBERSORT core | nu-SVR + clip + simplex, as in Newman 2015 | CIBERSORTx B/S-mode and absolute mode (use `scripts/R/run_cibersortx.sh`) |
| quanTIseq-style | NNLS + mRNA scaling + Other compartment | the full quanTIseq read-level pipeline |

## Demo cohorts actually fetched (2026-08-16)

| Cohort | File | n | Targets |
|---|---|---|---|
| GSE126044 | `GSE126044_counts.txt.gz` + series matrix | 16 (5 R / 11 NR; 11 fresh / 5 FFPE) | TACSTD2, CLDN4 present |
| GSE135222 | `GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz` + series matrix | 27 (21 PFS events) | TACSTD2, CLDN4 present |
| TCGA LUAD+LUSC | UCSC Xena HiSeqV2 + clinicalMatrix | 1017 primary tumours (515 LUAD / 502 LUSC) | TACSTD2, CLDN4 present |

Sample matching is by normalised ID, not column order. In GSE126044 the count
table and the series matrix disagree on the order of Dis_06 / Dis_07 / Dis_10.

Exact statistics cited in `playbook.md` are stored in
`results/demo/KEY_STATS.json` and the accompanying TSV files.
