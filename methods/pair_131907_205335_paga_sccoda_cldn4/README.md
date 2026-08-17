# Pair GSE131907 + GSE205335 — PAGA + scCODA, CLDN4 only

ADDITIVE high-end on the combo that already differs (PR #320:
author-malignant CLDN4 %pos vs T/NK, Q4 vs Q1 n=23 r=−0.705).

- CLDN4 only. No dual-high TACSTD2×CLDN4.
- Do not add GSE148071. Do not run the 7-cohort pool.
- Patient is the unit (GSE131907 T/NK extract is sample-level; that is
  stated on those rows).
- Previous 53k-cell Harmony/UMAP OOM'd. This path **subsamples per unit**.

## Done criterion

This PR has:

1. PAGA tables (`malignant_sample_level_spearman.tsv`, connectivities)
2. scCODA tables (`sccoda_tests_primary.tsv`, `composition_table.tsv`)
3. `FINDING.md` with honest n

## Reproduce

```bash
python3 methods/pair_131907_205335_paga_sccoda_cldn4/scripts/run_all.py
```

scCODA-only (no GEO download; uses PR #320 given tables):

```bash
python3 methods/pair_131907_205335_paga_sccoda_cldn4/scripts/run_sccoda.py
```
