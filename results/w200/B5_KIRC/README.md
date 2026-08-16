# B5_KIRC — open RCC ICI RNA, CLDN4 vs response

Honest analog of claim B5 (CLDN4-high → worse ICI response, user OR = 0.42)
restricted to **public** kidney-cancer ICI expression.

**Answer:** not supported. Braun nivolumab CR/PR vs PD: **OR = 0.84, n = 108,
p = 0.38**. Median-split OR = 1.27, p = 0.69. See `RESULTS.md`.

| File | What |
|---|---|
| `ANALYSIS_PLAN.md` | Pre-specified plan (committed before outcome tests) |
| `RESULTS.md` | Full write-up |
| `summary.json` | Machine-readable numbers |
| `braun_primary_patients.tsv` | n=108 primary table (enough to refit the OR) |
| `braun_nivo_patients.tsv` | All 181 nivo RNA cases |
| `javelin_patients.tsv` | 726 pts, PFS only |
| `gse67501_patients.tsv` | n=11 microarray |
| `reference_genes_primary.tsv` | Pre-named gene panel on the primary endpoint |
| `genomewide_top15.tsv` | Top 15 of 43,893 genes by logistic p (all FDR 0.93) |
| `fig_braun_primary_boxplot.png` | CLDN4 by CR/PR vs PD |
| `fig_or_forest.png` | Pre-specified logistic ORs |
| `provenance.md` | Download URLs and checksums |
