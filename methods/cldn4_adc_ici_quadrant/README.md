# CLDN4-only ADC+ICI quadrant

Additive public-bulk count. **No dual-high.** CLDN4 Q4 ∩ (IFN-γ Q4 ∪ CD274 Q4) versus CLDN4 Q4 ∩ IFN Q1 on GSE285029, GSE218989, GSE126044, GSE166449.

```bash
python3 -m pip install -r methods/cldn4_adc_ici_quadrant/requirements.txt
python3 methods/cldn4_adc_ici_quadrant/analyze.py
```

Primary table: `tables/quadrant_table.tsv`. Write-up: `FINDING.md`.
