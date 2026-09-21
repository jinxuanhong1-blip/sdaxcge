# DepMap lung TROP2 vs CLDN1/4/7

```bash
pip install -r scripts/depmap_lung_trop2_cldn147/requirements.txt
python3 scripts/depmap_lung_trop2_cldn147/download.py
python3 scripts/depmap_lung_trop2_cldn147/analyze.py
```

`download.py` reads the public portal catalog
`https://depmap.org/portal/api/no-captcha/download/files` and saves extracts
under `results/depmap_lung_trop2_cldn147/`. The full RNA matrix stays in the
cache directory (default `/tmp/depmap_trop2_cldn147`).
