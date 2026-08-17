# Additive leftover: mouse lung ICI sc/spatial Cldn4

TISMO 49/64 and PR #285 (no n≥2 epithelial ICB vs true control) are **taken as given**. This folder does not re-score that pool.

Public processed GEO only:

| Series | What it is | Contrast in *these* files |
|---|---|---|
| GSE297632 | LLC ± anti-PD-1 scRNA (tolerant residual) | n=1 vs 1; subcutaneous |
| GSE261890 | NSCLC ICI spatial (BMKMANU S1000) | 2 public libraries: Sen.Res.Veh vs Sen.Res.ICI |
| GSE285606 | Kras/p53 344SQ parental vs PD1R1 (140P) scRNA | n=2 vs 2; **both IgG**, not ICB vs control |
| GSE303943 | CMT167R subcutaneous PKCi vs solvent | **skip** — not ICB |

Question: epithelial **Cldn4** versus T-cell / exclusion scores. Honest n. Extra figures for any Cldn4-high + immune-low cut.

```bash
python3 methods/mouse_ici_cldn4_leftover/01_download.py
python3 methods/mouse_ici_cldn4_leftover/02_score_scrna.py
python3 methods/mouse_ici_cldn4_leftover/03_score_spatial.py
python3 methods/mouse_ici_cldn4_leftover/05_summarize.py
python3 methods/mouse_ici_cldn4_leftover/04_figures.py
```

Primary write-up: `FINDING.md`. Per-series table: `results/mouse_ici_cldn4_leftover/per_series_table.tsv`.
