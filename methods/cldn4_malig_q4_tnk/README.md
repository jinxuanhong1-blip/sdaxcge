# Malignant CLDN4-only: Q4 vs Q1 and combinatorial search vs T/NK / CXCL13+ / B

ADDITIVE. CLDN4-only. No dual-high with TACSTD2. Reuses existing public
lung scRNA **patient-level malignant** scores (PR #279 T/NK, PR #274 TLS/B
/ CXCL13+). All-epithelial definitions are dropped.

Two honest cuts, same vectors:

1. **Q4 vs Q1** — within-cohort quartiles of malignant CLDN4, Mann–Whitney
   on the immune score (rank-biserial *r*).
2. **Combinatorial CLDN4-negative** — Spearman ρ < 0 vs T/NK, CXCL13+ T,
   and B, with DerSimonian–Laird / Stouffer pools and honest n / p.

```bash
python3 methods/cldn4_malig_q4_tnk/analyze.py
```

Writeup: [`FINDING.md`](FINDING.md). Combo table:
`tables/highlighted_combos.tsv`. Extra figures: `figures/`.
