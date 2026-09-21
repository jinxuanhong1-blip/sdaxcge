# Concordant-4 Augur, Dirichlet, and ILR

Patient-level compositional models for how CLDN4-high tumors sit in the T / NK / myeloid simplex, plus an Augur-style ranking of which cell type's transcriptome tracks the CLDN4 quartile.

Cohorts are the locked four: GSE123902, GSE131907, GSE205335, GSE189357. n = 65 units. The malignant CLDN4 score is the locked one.

```bash
bash methods/concordant4_augur_dirichlet_ilr/scripts/download.sh /tmp/geo_c4
python3 methods/concordant4_augur_dirichlet_ilr/scripts/test_compositional.py
python3 methods/concordant4_augur_dirichlet_ilr/scripts/build_counts.py
python3 methods/concordant4_augur_dirichlet_ilr/scripts/fit_compositional.py
Rscript methods/concordant4_augur_dirichlet_ilr/scripts/extract_gse205335.R
python3 methods/concordant4_augur_dirichlet_ilr/scripts/extract_augur.py
python3 methods/concordant4_augur_dirichlet_ilr/scripts/run_augur.py
python3 methods/concordant4_augur_dirichlet_ilr/scripts/make_figures.py
```

See `METHODS.md` and `FINDING.md`.
