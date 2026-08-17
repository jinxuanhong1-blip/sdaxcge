# Scripts

- `00_download.py` — GET PRIDE SEARCH `proteinGroups_Cldn4.txt`; verify SHA-1. No RAW.
- `01_rank_partners.py` — filter, log2 LFQ, rank Cldn4 vs GFP, map IFN/MHC/TJ/trafficking/immune.

```bash
pip install -r requirements.txt
python3 00_download.py
python3 01_rank_partners.py
```
