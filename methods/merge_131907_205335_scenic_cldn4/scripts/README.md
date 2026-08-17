# scripts

- `download.py` — GEO public UMI + author labels for GSE131907 and GSE205335
  (skip 2.86 GB log2TPM, EGA, GSE207422).
- `download_priors.py` — TRRUST / DoRothEA / CollecTRI for IFN/MHC-I/TJ/keratin
  TFs (not ChIP). Snapshot is committed.
- `analyze.py` — author-malignant gate, within-patient CLDN4 split, AUCell
  proxy, writes FINDING.md and the regulon table with patient n.
