# FINDING — CLDN4-only high-end CellChat on the strict-malignant 6-unit combo

ADDITIVE. **CLDN4 only. No dual-high.** Patient is the unit.

The strict-malignant 6-unit mean-vs-T/NK Spearman is **taken as given** and is
**not re-audited**: n=95, ρ=−0.260, p=0.0195, I²=0% (PR #312).

Units: GSE207422 (12) + GSE205335 (22) + GSE291670 (6) + GSE253013 (9) +
GSE131907 (21) + GSE325414 (25).

GSE253013’s only public processed matrix is a 9.3 GB RDS and is **not
downloaded**. LR/meta tables are written by `scripts/analyze.py`.

See `METHODS.md`. Reproduce:

```bash
python3 methods/sixunit95_cldn4_hiend/scripts/download.py
python3 methods/sixunit95_cldn4_hiend/scripts/analyze.py
```
