# sdaxcge

Public-data analyses around TACSTD2 (TROP2) and CLDN4.

## B1 analog — TCGA-COAD surface rank

Honest ranking of CLDN4 among all in-silico surfaceome genes by
co-expression with TACSTD2 in TCGA-COAD.

Outputs: [`results/w200/B1_COAD/`](results/w200/B1_COAD/).

```bash
pip install -r requirements.txt
python3 scripts/download_data.py
python3 scripts/coexpression.py
```
