# Claim A10 analysis scripts

Run from the repo root, with the project venv:

```
.venv/bin/python scripts/01_human_bulk.py
.venv/bin/python scripts/02_perturbation.py
.venv/bin/python scripts/03_census_scrna.py
.venv/bin/python scripts/04_figures.py
.venv/bin/python scripts/05_write_report.py
```

`01` needs recount3 files in `data/raw/` (GTEx lung + TCGA-LUAD gene sums + G026 GTF).
`02` needs the GEO processed matrices in `data/raw/geo/`.
`03` needs network access to CELLxGENE Census (`cellxgene-census`).

Raw downloads are gitignored. Tables, figures, and `REPORT.md` are committed.
