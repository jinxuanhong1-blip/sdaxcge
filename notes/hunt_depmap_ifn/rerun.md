# How to rerun (this slice only)

```bash
python3 -m pip install pandas numpy scipy matplotlib statsmodels
python3 scripts/hunt_depmap_ifn/00_download.py --outdir data/hunt_depmap_ifn
python3 scripts/hunt_depmap_ifn/01_analyze.py --data data/hunt_depmap_ifn --outdir results/hunt_depmap_ifn
```

Raw DepMap matrices stay in `data/hunt_depmap_ifn/` (local cache, not committed; expression is ~483 MB).
Outputs: `results/hunt_depmap_ifn/tables/`, `results/hunt_depmap_ifn/figures/`, `results/hunt_depmap_ifn/WRITEUP.md`.

Citation: DepMap, Broad (2024). DepMap 24Q4 Public. Figshare+. https://doi.org/10.25452/figshare.plus.27993248.v1
Hallmark gene sets: MSigDB 2024.1 Hs.
