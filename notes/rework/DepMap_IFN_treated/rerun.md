# How to rerun

```bash
python3 -m pip install pandas numpy scipy matplotlib statsmodels
python3 scripts/rework/DepMap_IFN_treated/00_download.py --outdir data/rework/DepMap_IFN_treated
python3 scripts/rework/DepMap_IFN_treated/01_analyze.py \
  --data data/rework/DepMap_IFN_treated \
  --outdir results/rework/DepMap_IFN_treated
```

Raw matrices stay in `data/rework/DepMap_IFN_treated/` (not committed).

Citations:

- DepMap, Broad (2024). DepMap 24Q4 Public. Figshare+. https://doi.org/10.25452/figshare.plus.27993248.v1
- DepMap, Broad; Kocak, M. (2023). Repurposing Public 23Q2. figshare. https://doi.org/10.6084/m9.figshare.23600310.v4
- Nusinow et al., *Cell* 2020. CCLE proteomics. https://gygi.hms.harvard.edu/publications/ccle.html
