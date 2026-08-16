# B2_n118 — download + recompute

```bash
python3 -m pip install -r scripts/rework/B2_n118/requirements.txt
python3 scripts/rework/B2_n118/download.py --outdir data/B2_n118
python3 scripts/rework/B2_n118/analyze.py --data data/B2_n118 --outdir results/rework/B2_n118
```

Raw matrices stay in `data/B2_n118/` (gitignored). Results are self-contained under `results/rework/B2_n118/`.
