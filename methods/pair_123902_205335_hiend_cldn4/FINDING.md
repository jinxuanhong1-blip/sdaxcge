# FINDING — CLDN4-only high-end CellChat on GSE123902 + GSE205335

ADDITIVE. **CLDN4 only. No dual-high.** Patient (GSE123902 donor) is the unit.
Do **not** add GSE148071. Do **not** pile the seven-cohort pool.

The pair %pos Spearman is **taken as given** from PR #459 and is **not
re-audited**:

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE123902+GSE205335 | %pos | 35 | −0.522 (0.002, 0%) | −0.802 (0.005; 9 vs 9) |

Singles (given, same PR): GSE123902 n=13 ρ=−0.659; GSE205335 n=22 ρ=−0.435.

This folder adds CellChat-style **outgoing CLDN4-high malignant → same-patient
T/NK**. Jin et al. 2021 Hill probability on CellChatDB v2 protein pairs.
CellChat R is not run.

Primary high-end: within-patient malignant CLDN4 Q4 vs Q1. Sensitivity:
median split. Extra: between-patient Q4 vs Q1 using the given %pos scores.

The LR table (`results/lr_table.tsv`) is written by `scripts/analyze.py`
after the public UMIs are scored. Extra figures go in `figures/`.

## Reproduce

```bash
python3 methods/pair_123902_205335_hiend_cldn4/scripts/download.py
python3 methods/pair_123902_205335_hiend_cldn4/scripts/analyze.py
```
