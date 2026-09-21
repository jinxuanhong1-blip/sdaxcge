# Tacstd2, Cldn4, and immune readouts at the mouse unit

Public data only. GSE137244 cell-line libraries, plus unsorted public mouse lung scRNA. Private 8 KL matrices are not inputs and are not merged with these series.

`analyze.py` downloads nothing by itself. Place the GEO files listed in `download.sh` under `/tmp/geo_dl`, then run:

```bash
python3 analyze.py
```

The epithelial gate is the one used for the public GEMM rescore. On GSE264739 and GSE295824 it matches the published mouse table exactly (22/22 mice) before Tacstd2 is interpreted. Details and the primary numbers are in `FINDING.md`.
