# ICI purity residual rework

Self-contained re-analysis of open lung ICI bulk cohorts:

1. Score ESTIMATE stromal / immune / purity (Yoshihara 2013 port).
2. Residualize TACSTD2, CLDN4, and a locked TJ score on ESTIMATEScore.
3. Test residualized values against DCB or ORR (native binary endpoint).

```bash
python3 scripts/rework/ICI_purity/analyze.py
```

Writes `results/rework/ICI_purity/`. Data cache defaults to `/tmp/ici_purity_data`.
See that folder's `WRITEUP.md` for methods, every-cohort inventory, and numbers.
