# A11: TGF-β genes vs TACSTD2-high public lung

Pre-specified public-lung test of whether **TGF-β ligands / receptors / SMADs**
and a **HALLMARK TGF-β signaling score** are higher in **TACSTD2-high** tumors.

Primary cohort: TCGA LUAD + LUSC primary tumors (Xena GDC STAR TPM, log2(TPM+1)).
Sensitivity: DepMap Public 24Q4 lung cell lines, if the extract is present.

Primary TGF-β ligands: **TGFB1, TGFB2, TGFB3**.
Receptors and SMADs are reported separately.
**HALLMARK_TGF_BETA_SIGNALING** z-mean is the pathway-score test named in claim A11.
**CLDN4** is a junction positive control. **CD8A** is immune context only.

No filter is tuned to produce a positive class effect. Nulls, negatives, and
histology splits are reported.

```bash
python3 scripts/w200/A11_TGFB/download.py
python3 scripts/w200/A11_TGFB/analyze.py
```

Outputs land in `results/w200/A11_TGFB/`.
