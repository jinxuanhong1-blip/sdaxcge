# Scripts

- `download.py` — GSE123902 RAW tar + GSE131907 UMI / author annotation.
  Skips GSE148071, the 7-cohort pool, the 36.5 GB Laughney H5, and the
  2.9 GB GSE131907 log2TPM text.
- `analyze.py` — CellChat-style Hill *P* and LIANA-style mean-of-means
  outgoing CLDN4-high malignant → T/NK. Patient/donor (GSE131907 sample)
  is the unit.
- `lib.py` — Hill probability, CellPhoneDB-style mean-of-means, high-end
  splits, patient-ΔP meta.
