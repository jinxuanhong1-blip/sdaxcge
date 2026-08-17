# Winning pair GSE131907+GSE205335 — propeller + CellPhoneDB (CLDN4 only)

ADDITIVE. CLDN4 only. No dual-high. No GSE207422.

Writeup: [`FINDING.md`](FINDING.md).

| Table | File |
|---|---|
| Propeller / speckle-style T/NK/B | `tables/propeller_composition.tsv` |
| CellPhoneDB-style MHC-I + T-recruit | `tables/cpdb_focus.tsv` |

```bash
python3 methods/winpair_131907_205335_propeller_cpdb_cldn4/scripts/download.py
python3 methods/winpair_131907_205335_propeller_cpdb_cldn4/scripts/analyze.py
```
