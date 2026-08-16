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

## Honest result (this run)

| Pool | Cutoff | Endpoint | k | RE OR [95% CI] | p |
|---|---|---|---:|---|---:|
| all independent (locked primary) | median | curated R vs NR | 9 | **0.70 [0.41–1.18]** | 0.18 |
| all independent | continuous | curated R vs NR | 9 | 0.92 [0.78–1.09] | 0.35 |
| all independent | tertile T3 vs T1 | curated R vs NR | 9 | 0.79 [0.51–1.22] | 0.28 |
| all independent | quartile Q4 vs Q1 | curated R vs NR | 9 | 0.80 [0.49–1.29] | 0.35 |
| lung only | median | curated R vs NR | 4 | 0.40 [0.15–1.04] | 0.059 |

The lung-only median OR (0.40) is the cell closest to the user 0.42. It is
**not adopted**. It is NS, k=4 not 11, and the locked primary is 0.70.
IMvigor210 median OR is 1.63 (opposite direction). Public k=9, not 11.
No pooled cell matches 0.42 at 2 decimal places.
