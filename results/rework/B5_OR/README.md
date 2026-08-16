# B5 rework: open ICI RNA, CLDN4 vs response, all cutoffs

User-reported claim: **CLDN4-high ICI OR = 0.42 [0.18–0.95], k=11**.
Prior note: open lung ICI was NS.

This slice downloads public ICI RNA matrices (lung GEO + Zenodo ICB
IMvigor210/Gide/Liu/Riaz/Braun), extracts CLDN4, and computes **median,
tertile, quartile, and continuous** associations. Every cutoff is written.
The script does not pick the cutoff that lands on 0.42.

## Run

```bash
python3 -m pip install -r results/rework/B5_OR/requirements.txt
python3 results/rework/B5_OR/download.py
python3 results/rework/B5_OR/analyze.py
```

Raw archives go to `results/rework/B5_OR/raw/` (gitignored).
Patient-level CLDN4 tables go to `processed/`.
Stats and forest plots go to `tables/` and `figures/`.

## Locked primary (declared before pooling)

- Cutoff: median
- Endpoint: curated R vs NR
- Model: DerSimonian–Laird random effects
- Studies: the independent public matrices listed in `analyze.py` (`PRIMARY_COHORTS`)

All other cutoff × endpoint × subset cells are still reported in
`tables/pooled_or.csv`.
