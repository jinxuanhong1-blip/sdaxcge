# Harmony n=73 — CLDN4-only high-end (CellChat + LIANA style)

Additive analysis on the given Harmony / multi-cohort patient table
(malignant-like CLDN4 vs T/NK n=73, ρ=−0.27; not re-audited).

```bash
python3 scripts/00_download.py   # public processed files only; no GSE253013 RDS
python3 scripts/01_run_lr.py
```

See `FINDING.md`.
