# Rerun — C4 analog GSE245459

```
pip install pandas numpy scipy matplotlib statsmodels
python3 scripts/w200/C4_GSE245459/download.py   # GEO suppl FPKM, ~14 MB
python3 scripts/w200/C4_GSE245459/analyze.py
```

Source file (not committed; re-download):

`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245459/suppl/GSE245459_fpkm.anno.txt.gz`

Slim gene × 12-sample FPKM is written to `results/w200/C4_GSE245459/tables/fpkm_gene_matrix.tsv.gz`.
