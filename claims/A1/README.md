# Claim A1

TACSTD2 (TROP2) vs immune / cytotoxic / exhaustion signatures in TCGA + OncoSG
NSCLC, purity-adjusted partial Spearman correlation.

- `claim.md` — full claim statement, scope, and prediction
- `signatures.json` — target gene, signature gene sets, scoring, references
- `config.json` — cohorts, cBioPortal profiles, purity sources

Run the analysis from the repo root:

```
pip install -r requirements.txt
python code/run_claim_A1.py
python code/make_figures.py
```

Results are written to `results/claim_A1/` (see `results/claim_A1/report.md` for
the verdict). Public data only.
