# Concordant-4 scVI/scANVI: max patient-level CLDN4 vs T/NK effect

Same 65 units and the same patient-batch scVI CLDN4 as the concordant-4 scVI/scANVI integration. The script searches scVI/scANVI summaries, two T/NK denominators, and three codings. The primary row is the one that raises both the patient-level binomial likelihood-ratio chi-square and the absolute pooled Spearman relative to malignant CLDN4 % positive.

```bash
python3 methods/scanvi_c4_max_effect/maximize.py
```

Write-up: `FINDING.md`.
