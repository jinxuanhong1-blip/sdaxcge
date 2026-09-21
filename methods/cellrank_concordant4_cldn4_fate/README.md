# CellRank 2 fate toward barrier CLDN4-high (concordant-4)

ADDITIVE. **CLDN4-only.** Tumor epithelium from **GSE123902 + GSE131907 + GSE205335 + GSE189357**.

This does not re-estimate the locked patient-level result (malignant CLDN4 %pos vs T/NK, n=65, ρ=−0.531). Cell counts are not that n.

CellRank 2 is the primary fate model. A kNN contact matrix and a multinomial logistic of neighbor state are secondary descriptions of the same graph. RNA velocity is off.

```bash
python3 methods/cellrank_concordant4_cldn4_fate/download.py
python3 methods/cellrank_concordant4_cldn4_fate/analyze.py
```

Numbers and the figure are written by `analyze.py` into `FINDING.md` and `results/`.
