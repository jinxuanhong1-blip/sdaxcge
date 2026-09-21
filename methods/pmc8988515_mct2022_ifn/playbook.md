# PMC8988515 IFN side

Yamamoto, Webb, Bitler et al., Mol Cancer Ther 2022 (PMC8988515, MCT-21-0827). The paper's own data-availability line is "upon request." This page hunts GEO, ArrayExpress, PRIDE, and NODE, then scores interferon on the expression tables that actually download.

```bash
python3 methods/pmc8988515_mct2022_ifn/download.py
python3 methods/pmc8988515_mct2022_ifn/analyze.py
```

Raw matrices stay in `data/` and are not committed. Result tables are under `tables/`.
