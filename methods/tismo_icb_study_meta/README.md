# TISMO ICB study-level REML

Paired pre/post meta-analysis of Tacstd2, Cldn4, and the tight-junction score on the locked 64 TISMO slices (49/64 Tacstd2). Studies are the GEO/ENA accessions. τ² is REML. Outputs are forest plots and leave-one-study-out fits.

LLC is labeled LLC. There is no KL or KP lung line in this ICB table.

```bash
python3 methods/tismo_icb_study_meta/analyze.py
python3 -m unittest tests/test_tismo_icb_study_meta.py
```

Write-up: `RESULTS.md`.
