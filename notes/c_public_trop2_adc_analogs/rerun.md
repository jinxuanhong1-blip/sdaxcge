# Rerun

```bash
python3 -m pip install -r scripts/c_public_trop2_adc_analogs/requirements.txt
python3 scripts/c_public_trop2_adc_analogs/download.py
python3 scripts/c_public_trop2_adc_analogs/analyze.py
```

Outputs: `results/c_public_trop2_adc_analogs/{WRITEUP.md,key_stats.json,tables/,figures/}`.

Raw GEO files stay in `results/c_public_trop2_adc_analogs/raw/` (git-ignored).
