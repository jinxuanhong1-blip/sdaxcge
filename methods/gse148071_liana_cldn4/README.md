# GSE148071 LIANA / LR — CLDN4-high malignant → T/NK

Additive ligand–receptor package on public TISCH2 `NSCLC_GSE148071`
(Wu et al. 2021, GEO GSE148071).

**Done criterion:** `FINDING.md` contains an LR table and honest patient/cell n.

```bash
python3 scripts/00_download.py
python3 scripts/01_run_liana.py
```

Primary method is a documented CellPhoneDB-style score. LIANA is secondary if
importable. CellChat is not faked.
