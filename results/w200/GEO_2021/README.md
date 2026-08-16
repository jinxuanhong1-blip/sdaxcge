# GEO 2021 lung ICI leftovers: TACSTD2 / CLDN4

Eight previously unresolved 2021 GEO candidates were manually checked against their deposited assays and sample metadata. Six cannot answer the question: they are miRNA, cell-line, T-cell, targeted-panel, or tumor datasets without a matched ICI outcome. `leftover_triage.csv` records every exclusion.

Two cohorts are usable. [GSE190265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190265) has tumor TPM plus PFS for 43 samples. Per doubling-like unit of log2(TPM+1), TACSTD2 had HR 1.08 (95% CI 0.90–1.31; p=0.392) and CLDN4 HR 1.08 (0.83–1.40; p=0.580). [GSE190266](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190266) has six-month PFS for 70 tumors but omits TACSTD2; CLDN4 had HR 0.92 (0.82–1.02; p=0.125). Median-split log-rank p-values were 0.826, 0.939, and 0.210. None survives even the three-test BH correction (smallest q=0.374).

Honest conclusion: these leftovers provide no evidence that either marker predicts ICI PFS. The analyses are univariate, small, retrospective, platform-heterogeneous, and not a pooled validation. GSE190265’s supplementary files contain 43 aligned records although GEO exposes only 34 GSM entries; this provenance discrepancy is retained, not hidden.

`GSE146100_note.md` separately reports a descriptive single-cell check; it is not patient-level biomarker evidence. Run `python3 analyze_leftovers.py` to reproduce the survival table and checksums.
