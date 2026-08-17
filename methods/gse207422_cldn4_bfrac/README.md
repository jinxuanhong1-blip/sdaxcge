# GSE207422 malignant CLDN4 vs B-cell / TLS fraction

ADDITIVE public-UMI slice. Patient-level B-cell and TLS-like fractions versus
malignant CLDN4. Honest n is the 12 post-treatment patients (not 15 including
pre-biopsy; not cell-level). TACSTD2 is a companion gene, never a gate.

See `FINDING.md`.

```bash
python3 methods/gse207422_cldn4_bfrac/scripts/download.py
python3 methods/gse207422_cldn4_bfrac/scripts/extract.py
python3 methods/gse207422_cldn4_bfrac/scripts/analyze.py
```
