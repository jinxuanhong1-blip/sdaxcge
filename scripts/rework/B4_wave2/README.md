# B4 extras (public lung ICI TJ/CLDN4 + TCGA)

GSE126044 TJ/NR is the index cohort (taken as given). This folder adds other
public lung ICI matrices and TCGA-LUAD/LUSC TJ vs CD8/GEP.

```bash
python3 -m pip install -r scripts/rework/B4_wave2/requirements.txt
python3 scripts/rework/B4_wave2/download.py
python3 scripts/rework/B4_wave2/extra_cohorts.py
```

Optional GSE126044 index recompute: `python3 scripts/rework/B4_wave2/analyze.py`

Outputs: `results/rework/B4_wave2/extra/`.
