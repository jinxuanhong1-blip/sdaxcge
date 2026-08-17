# Pair GSE189357 + GSE205335 — CLDN4-only high-end CellChat-style

ADDITIVE. **CLDN4 only.** No dual-high TACSTD2∩CLDN4.

The combo that **differs** is taken from PR #459 and is not re-audited:

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p) |
| --- | --- | ---: | --- | --- |
| GSE189357+GSE205335 | %pos | 31 | −0.478 (0.009, 0%) | −0.750 (0.010; 8/8) |

High-end work here: outgoing **CLDN4-high malignant → same-patient T/NK** (Jin et al. 2021 Hill *P*). Patient is the unit. Honest paired n.

```bash
python3 methods/pair_189357_205335_hiend_cldn4/scripts/download.py
python3 methods/pair_189357_205335_hiend_cldn4/scripts/analyze.py
```

Write-up: `FINDING.md`. Primary table: `results/ligand_table.tsv`.
