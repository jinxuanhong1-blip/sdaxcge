# GSE285029: TACSTD2 / CLDN4 vs PD-L1, IFN/MHC-I, IL-6/STAT3, CD8/GEP

Public bulk RNA from Koh et al., *JITC* 2025 (PMID 40050048; GEO
[GSE285029](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285029)).
Pre-ICI NSCLC tumor WTS, **n = 234** (author ICI–RNA-seq cohort).

This folder is **additive combination-rationale evidence only**. User A11
(galectin / nectin / TGF-β / CD47 with TROP2) and user C (SKB264 raises
PD-L1, so ADC+ICI is required) are taken as given. The question here is
whether **TROP2-high / CLDN4-high tumors already sit on a high PD-L1 or
high IFN program** in this public matrix.

GEO does not release RECIST, histology, PD-L1 IHC, TMB, or purity.
This analysis does not use paywalled supplemental patient-level response
labels.

```bash
python3 scripts/w200/A11_GSE285029/download.py
python3 scripts/w200/A11_GSE285029/analyze.py
```

Outputs land in `results/w200/A11_GSE285029/`.
