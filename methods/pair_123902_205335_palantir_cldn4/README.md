# Pair GSE123902+GSE205335 — CLDN4-only REAL Palantir

ADDITIVE. **CLDN4 only.** No dual-high. No GSE148071.

This is the pair that differs: malignant IFN family DE is already **−1.05**
(Q4 vs Q1). That DE is taken as given. This folder runs **real Palantir**
(Setty 2019) and writes destiny tables for CLDN4, barrier (CLDN4 held out),
and IFN. Root is not CLDN4-high. Donor/patient is the unit.

Done when `results/tables/destiny_*.tsv` exist.

```bash
pip install -r methods/pair_123902_205335_palantir_cldn4/requirements.txt
python3 methods/pair_123902_205335_palantir_cldn4/scripts/download.py
python3 methods/pair_123902_205335_palantir_cldn4/scripts/extract.py
python3 methods/pair_123902_205335_palantir_cldn4/scripts/analyze.py \
  --input /tmp/pair_123902_205335_palantir/epithelium.h5ad \
  --outdir methods/pair_123902_205335_palantir_cldn4/results \
  --finding methods/pair_123902_205335_palantir_cldn4/FINDING.md
```
