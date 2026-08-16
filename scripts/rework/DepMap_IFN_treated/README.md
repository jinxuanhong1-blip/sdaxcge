# DepMap IFN-treated conceptual rework

Closed mismatch from `results/hunt_depmap_ifn`: unstimulated NSCLC
TACSTD2 vs Hallmark IFN-α RNA was **ρ = +0.377**, not negative.

This slice does not recompute that RNA test. It asks four leftover
questions that the basal transcriptome cannot answer:

1. PRISM IFN-pathway viability (JAK inhibitors; IFN inducers / IFNR agonists).
   Recombinant IFN-α/γ proteins are **not** in PRISM.
2. CRISPR Chronos gene effects of IFN / APM genes vs TACSTD2 RNA.
3. CCLE MS protein ISG scores.
4. MHC-I **protein** (HLA-A/B/C, B2M), not RNA.

```bash
python3 -m pip install pandas numpy scipy matplotlib statsmodels
python3 scripts/rework/DepMap_IFN_treated/00_download.py
python3 scripts/rework/DepMap_IFN_treated/01_analyze.py
```

Outputs: `results/rework/DepMap_IFN_treated/`.
