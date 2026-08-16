# How to rerun (this slice only)

```bash
python3 -m pip install pandas numpy scipy matplotlib statsmodels
python3 scripts/grok_cptac_luad/00_download.py --outdir data/grok_cptac_luad
python3 scripts/grok_cptac_luad/01_analyze.py --data data/grok_cptac_luad --outdir results/grok_cptac_luad
```

Raw matrices stay in `data/grok_cptac_luad/` (local cache, not committed).
Outputs: `results/grok_cptac_luad/tables/`, `results/grok_cptac_luad/figures/`, `results/grok_cptac_luad/WRITEUP.md`.
