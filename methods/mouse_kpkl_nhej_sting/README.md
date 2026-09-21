# Public mouse lung KP/KL beyond GSE137244

Cldn4, an NHEJ module, and a cGAS–STING module on open public matrices. Each accession is tested alone. GSE137244 is not re-opened. Private 8 KL mice are not used and are not merged with any public series.

```bash
python3 methods/mouse_kpkl_nhej_sting/analyze.py
```

Processed GEO files are cached under `/tmp/kpkl` (override with `KPKL_DATA`). Tables land in `results/mouse_kpkl_nhej_sting/`.
