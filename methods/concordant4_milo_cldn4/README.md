# Concordant-4 Milo neighbourhoods vs patient-level CLDN4

Additive differential abundance on the locked concordant-4 cohort
(GSE123902 + GSE131907 + GSE205335 + GSE189357). The tested covariate is
the patient-level malignant CLDN4 group. The graph is the committed
Harmony embedding from the concordant-4 Seurat object.

Neighbourhoods are transcriptional states. The unit is the patient,
donor, or tumour-bearing sample already used for that cohort (n = 65).
Do not quote the cell count as n.

## Design

- High-rich vs low-rich: within-cohort quartiles on the locked table.
  High-rich = Q3+Q4 (31). Low-rich = Q1+Q2 (34).
- Primary model: `~ dataset + cldn4_high`, edgeR quasi-likelihood F test
  (the GLM inside `miloR::testNhoods`), TMM, k-distance SpatialFDR.
- A patient intercept, `~ patient_id + cldn4_high`, is rank 65/66.
  Each patient contributes one sample, so that term is not estimated.
- A within-patient malignant contrast is unidentified for the same reason
  (0 patients with two samples). It is not fit.

## Reproduce

```bash
pip install -r methods/concordant4_milo_cldn4/requirements.txt
# R packages edgeR, limma, jsonlite (user library is fine)
python3 methods/concordant4_milo_cldn4/scripts/run_milo.py
python3 methods/concordant4_milo_cldn4/scripts/deepen_milo.py
```

`scripts/nhood_gene_scores.py` reads the public GEO count matrices
from `/tmp/geo_c4` (GSE123902 raw tar, GSE131907 UMI matrix,
GSE189357 raw tar, GSE205335 UMI RDS and cell-identity table) and
writes the within-unit IFN / NHEJ / STING tables. The abundance
scripts do not need those matrices.

Write-up: `FINDING.md`.
