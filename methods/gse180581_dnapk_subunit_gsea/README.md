# GSE180581 DNA-PK subunit GSEA

Matched HEK293T knockdowns of DNA-PKcs (PRKDC), Ku70 (XRCC6), and Ku80 (XRCC5). Prerank GSEA of IFN, STING, antigen-processing, and cytosolic DNA-sensing sets. Positive NES is up in the knockdown versus the matched siControl.

```bash
python3 scripts/gse180581_dnapk_subunit_gsea/analyze.py
```

Needs numpy, pandas, scipy, and matplotlib. Counts and gene sets are already under `data/`. The result write-up is `FINDING.md`.
