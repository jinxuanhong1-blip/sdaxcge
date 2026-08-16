# A11: nectin genes vs TACSTD2-high public lung

Pre-specified public-lung test of whether **nectin-family genes** are higher in
**TACSTD2-high** tumors, or merely correlated with TACSTD2.

Primary cohort: TCGA LUAD + LUSC primary tumors (Xena GDC STAR TPM, log2(TPM+1)).
Sensitivity: DepMap Public 24Q4 lung cell lines (log2(TPM+1)), if the extract is present.

The four nectin genes are **NECTIN1, NECTIN2, NECTIN3, NECTIN4**.
**PVR** (CD155) is nectin-like and is reported separately.
**CLDN4** is a positive-control junction gene, not a nectin.

No filter is tuned to produce a positive class effect. Nulls and histology
splits are reported.

```bash
python3 scripts/w200/A11_nectin/download.py
python3 scripts/w200/A11_nectin/analyze.py
```

Outputs land in `results/w200/A11_nectin/`.
