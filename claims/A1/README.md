# Claim A1

TACSTD2 (TROP2) vs immune / cytotoxic / exhaustion signatures in TCGA
(LUAD+LUSC; PanCanAtlas and Firehose) + OncoSG NSCLC, purity-adjusted
partial Spearman under ESTIMATE, ABSOLUTE, CPE, and ESTIMATE-or-ABSOLUTE.

- `claim.md` — claim statement, scope, prediction
- `signatures.json` — target gene, signature gene sets, scoring, references
- `config.json` — cohorts, cBioPortal profiles, purity sources

```
pip install -r requirements.txt
python code/run_claim_A1.py
python code/make_figures.py
```

Results: `results/claim_A1/` (see `report.md` for n, ρ, p, and match/mismatch).
Public data only.
