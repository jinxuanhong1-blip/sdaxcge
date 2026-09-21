# CPTAC LUAD/LSCC protein: TROP2, partial CLDN4, EPCAM control

Public TMT freeze v1.2 tumor protein. LUAD and LSCC stay separate. No imputation.

```bash
python methods/cptac_trop2_cldn4_mediation/download.py
python methods/cptac_trop2_cldn4_mediation/analyze.py
```

`download.py` checks the protein and phenotype SHA-256 values against the manifest from the earlier MHC/IFN protein page. Matrices are gitignored. Coefficients are written by `analyze.py` into `FINDING.md` and `tables/`.
