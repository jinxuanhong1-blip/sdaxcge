# Public mouse KL vs KP beyond GSE137244

KL is Kras-mutant, Lkb1/Stk11-null, Trp53 intact. KP is Kras-mutant, Trp53-null, Lkb1 intact. KPL (both tumor suppressors lost) is scored only as a labeled secondary contrast.

```bash
python3 methods/public_kl_kp_beyond_gse137244/analyze.py
```

Downloads GEO supplementary files into `cache/` (not committed). Tables land in `tables/`. The test is a two-sided exact Mann–Whitney on log2(abundance+1), delta = mean(KL) − mean(KP), matching the locked GSE137244 cell-line scale.

Private 8KL matrices are not used.
