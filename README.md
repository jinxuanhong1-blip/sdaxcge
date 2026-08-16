# sdaxcge

Public-data analyses of TACSTD2 (TROP2) and CLDN4.

This branch adds the **B1 analog** for TCGA-PRAD: an honest surfaceome
co-expression ranking of `CLDN4` against `TACSTD2`. See
`results/w200/B1_PRAD/`.

```bash
pip install -r requirements.txt
python scripts/download_data.py
python scripts/coexpression.py
```
