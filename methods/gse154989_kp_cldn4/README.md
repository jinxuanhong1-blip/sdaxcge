# GSE154989 KP GEMM — Cldn4-only (additive public mouse)

Public processed Smart-seq2 of FACS `tdTomato+ / CD45− / CD11b− / TER119− / CD31−`
KP (and K / T) lung epithelium (Marjanovic et al., *Cancer Cell* 2020).

- Matrix: GEO `GSE154989_mmLungPlate_fQC_dSp_normTPM.h5` only.
- Unit: biological mouse (collapse `_T#` tumors). Honest n, not cell n.
- Cldn4 only. No dual-high. No T/NK fraction (immune cells were FACS-excluded).
- Testable arm: epithelial Cldn4 vs IFN / MHC-I/APM / TJ (Cldn4 held out).

```bash
python3 methods/gse154989_kp_cldn4/analyze.py
```
