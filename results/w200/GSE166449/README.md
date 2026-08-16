# w200 / GSE166449 — TACSTD2 and CLDN4 vs pembrolizumab response

Reproduce from the repo root:

```bash
python3 results/w200/GSE166449/download.py
python3 results/w200/GSE166449/analyze.py
```

or:

```bash
python3 scripts/w200/GSE166449.py
```

The honest result is in `WRITEUP.md` and `summary.json`.
Neither gene, the two-gene score, nor the CD3D/CD3E/CD8A T-cell proxy
is associated with response in this 7 vs 15 cohort.
