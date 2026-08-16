# GEO 2021 leftover lung ICI: TACSTD2 / CLDN4

A 2021-only GEO search found 38 real human series. After removing two already-analyzed cohorts, **only GSE190265 and GSE190266** still have tumor mRNA plus an ICI endpoint. Full exclusions: `leftover_triage.csv`. Interpretation: `WRITEUP.md`.

Leftover result: **no support** for TACSTD2 or CLDN4 as ICI PFS/DCB markers. France4 TACSTD2 is missing (16,383-gene Excel-style cap, last gene MTMR14). France4 CLDN4 DCB Wilcoxon p=0.052 does not survive BH correction and is absent in France3. The same files do detect CXCL10/CD274–DCB signal, so the target-gene null is not a dead assay.

```bash
python3 search_2021.py
python3 analyze_leftovers.py
```
