# scripts

| file | role |
| --- | --- |
| `download_gse207422.py` | GEO processed UMI + metadata |
| `download_priors.py` | refresh TRRUST / DoRothEA / CollecTRI snapshot |
| `analyze.py` | gate epithelium, Pearson±GRNBoost2, AUCell, LUAD vs all, figures |

`analyze.py` looks for the matrix in `/tmp/scrna_scenic/geo` or `data/scrna_scenic/GSE207422`.
