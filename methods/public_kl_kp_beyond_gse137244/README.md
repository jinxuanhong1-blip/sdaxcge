# Public mouse KL vs KP beyond GSE137244

KL is Kras-mutant, Lkb1/Stk11-null, Trp53 intact. KP is Kras-mutant, Trp53-null, Lkb1 intact. KPL (both tumor suppressors lost) is scored only as a labeled secondary contrast.

```bash
python3 methods/public_kl_kp_beyond_gse137244/analyze.py
python3 methods/public_kl_kp_beyond_gse137244/sweep.py
```

Downloads GEO supplementary files into `cache/` (not committed). Tables land in `tables/`. The test is a two-sided exact Mann–Whitney, delta = mean(KL) − mean(KP), matching the locked GSE137244 cell-line scale. `sweep.py` adds histology filters, GPL8321 probes, epithelial normalization, tight-junction modules, and NHEJ/STING/IFN scores. Welch p is reported beside the rank test. Squamous-versus-adenocarcinoma contrasts are stored with `fair_genotype=False`.

Private 8KL matrices are not used.
