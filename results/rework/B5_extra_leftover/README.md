# B5 extra leftover lung ICI rows

Additive **2023–2026** open lung ICI cohorts that the user’s 11-cohort CLDN4 slide did not use.

- The 11-cohort meta (**OR = 0.42**) is taken as given and is **not** re-cut or re-pooled here.
- Extra rows live in `tables/extra_rows_NSCLC.tsv` (NSCLC) and `tables/extra_rows_SCLC.tsv` (SCLC, separate).
- Named leftovers without a patient DCB/ORR + tumor CLDN4 test are inventoried, not scored as ORs.

```bash
python3 results/rework/B5_extra_leftover/analyze.py
```

See `WRITEUP.md` for honest n / OR / p.
