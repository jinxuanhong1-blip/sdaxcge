# Scripts

- `download.py` — public processed GEO files only (GSE189357 10x MTX tar; GSE205335 UMI RDS + identity + SOFT).
- `analyze.py` — within-patient CLDN4-high vs low malignant → same-patient T/NK (CellChat-style Hill *P*).

```bash
python3 methods/pair_189357_205335_hiend_cldn4/scripts/download.py
python3 methods/pair_189357_205335_hiend_cldn4/scripts/analyze.py
```
