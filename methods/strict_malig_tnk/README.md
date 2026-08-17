# Strict malignant T/NK recut

ADDITIVE combinatorial search of **CLDN4 and TACSTD2 versus T/NK** on existing
scRNA patient tables.

- Keep only **author-malignant**, **marker-malignant**, or **DRMref** rows.
- Drop all-epithelial (including GSE241934 residual Epi).
- Enumerate cohort subsets until a cut has **ρ ≤ −0.35**, with honest n and p.

```bash
python3 methods/strict_malig_tnk/combinatorial_search.py
```

Writeup: `FINDING.md`. Combo table: `tables/combo_table.tsv`.
