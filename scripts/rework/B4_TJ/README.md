# B4_TJ recompute

Self-contained recompute of user claim **B4**: GSE126044 non-responders have a higher tight-junction score, p=0.019.

Also reports **CLDN4** alone and the same features in **GSE135222**.

```bash
python3 scripts/rework/B4_TJ/download.py
python3 scripts/rework/B4_TJ/analyze.py
```

Outputs land in `results/rework/B4_TJ/`. Nothing is tuned to p=0.019.
