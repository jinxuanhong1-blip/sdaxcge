# Concordant-4 CLDN4 vs T/NK maximum

Locked cohorts only: GSE123902, GSE131907, GSE205335, GSE189357.
The sweep searches expression cuts, percent-positive versus mean, keratin partials, and histology for the largest |Spearman ρ| that stays negative in all four cohorts, and reports the matching Q4 vs Q1 Cliff delta.

GEO matrices are not in this repo. `scripts/extract_scores.py` reads them from `/tmp/geo_dl`.

```bash
python3 methods/concordant4_cldn4_maxrho/scripts/extract_scores.py
python3 methods/concordant4_cldn4_maxrho/scripts/sweep.py
```

Numbers are written by the sweep to `FINDING.md`. Do not hand-edit them.
