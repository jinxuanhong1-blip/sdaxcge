# GSE334497 — Trop2 KO vs a CLDN4-KD direction

Public 4T1 Trop2 knockout tumors (Wu *et al.*, JITC 2026). Question: does Trop2 KO also drop *Cldn4* and raise IFN/immune RNA?

The pre-specified call is in [FINDING.md](FINDING.md). The method / filter / gene-set sweep, including STING and NHEJ, is in [SWEEP.md](SWEEP.md).

```bash
python3 scripts/gse334497_cldn4_match/analyze.py
python3 scripts/gse334497_cldn4_match/sweep.py
```
