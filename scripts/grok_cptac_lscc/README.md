# CPTAC LSCC protein + RNA + immune (parallel slice)

Outputs are restricted to `notes/grok_cptac_lscc/`, `scripts/grok_cptac_lscc/`, and `results/grok_cptac_lscc/`.

```bash
pip install -r scripts/grok_cptac_lscc/requirements.txt
python3 scripts/grok_cptac_lscc/00_download_freeze.py
python3 scripts/grok_cptac_lscc/01_analyze_lscc.py
```

Targets: TACSTD2 `ENSG00000184292`, CLDN4 `ENSG00000189143`.
Cohort is treatment-naive (no ICI labels). See `results/grok_cptac_lscc/WRITEUP.md`.
