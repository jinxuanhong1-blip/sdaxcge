# Concordant-4 CLDN4-only (GSE123902 + GSE131907 + GSE205335 + GSE189357)

ADDITIVE. CLDN4 only. No TACSTD2∩CLDN4 dual-high. Not the 7-cohort pool.

Two primary analyses, patient/donor/sample as the unit:

1. Malignant CLDN4 %pos vs same-unit T/NK (Spearman + Q4 vs Q1 + LOO).
2. Malignant patient-pseudobulk DE, within-cohort Q4 vs Q1 (IFN / MHC-I/APM / chemokine / TJ).

```bash
python3 methods/concordant4_123902_131907_205335_189357_cldn4/analyze.py
```

Writeup: `FINDING.md`.
