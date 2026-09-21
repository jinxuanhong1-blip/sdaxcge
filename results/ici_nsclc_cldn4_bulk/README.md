# Open NSCLC ICI bulk: CLDN4-high versus response

Median-split odds ratios and PFS hazard ratios for CLDN4 in every open GEO/ArrayExpress NSCLC tumor-bulk cohort with a clear responder label or PFS, plus the same contrast after residualizing CLDN4 on KRT8, KRT18, and KRT19.

The PDF figure OR = 0.42 (k = 11) is reported only if this meta actually returns those numbers. See `WRITEUP.md`.

```bash
python3 -m pip install -r results/ici_nsclc_cldn4_bulk/requirements.txt
python3 results/ici_nsclc_cldn4_bulk/analyze.py
```

Public matrices are read from `/tmp/ici_raw` (GEO FTP, EuropePMC PMC11669362 Table S1, Springer Supplementary Data 8). They are not committed.
