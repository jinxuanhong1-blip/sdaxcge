# GSE180581 DNA-PK subunit GSEA

Matched HEK293T knockdowns of DNA-PKcs (PRKDC), Ku70 (XRCC6), and Ku80 (XRCC5). The reported IFN increase is Ku70 knockdown versus parental siControl, ranked by log2 fold-change. A second rank compares each siRNA with its matched siControl. Positive NES is up in the knockdown.

```bash
python3 scripts/gse180581_dnapk_subunit_gsea/analyze.py
```

Needs numpy, pandas, scipy, and matplotlib. Counts and gene sets are already under `data/`. The result write-up is `FINDING.md`.
