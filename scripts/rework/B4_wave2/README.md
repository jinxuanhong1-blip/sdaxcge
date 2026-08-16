# B4 wave-2 (GSE126044 TJ p=0.019)

Public GSE126044 counts only. Tries author/paper TJ lists, one-sided MW,
fresh-only, ESTIMATE residualization, tertile splits, and CLDN4 vs
CLDN1/4/7/F11R/PARD3 vs keratin.

```bash
python3 -m pip install -r scripts/rework/B4_wave2/requirements.txt
python3 scripts/rework/B4_wave2/download.py
python3 scripts/rework/B4_wave2/analyze.py
```

Outputs: `results/rework/B4_wave2/`.

Additive paper extras (other public lung ICI cohorts + TCGA):

```bash
python3 scripts/rework/B4_wave2/download.py
python3 scripts/rework/B4_wave2/extra_cohorts.py
```

Outputs: `results/rework/B4_wave2/extra/`.
