# Concordant-4 signatures on TCGA-LUAD

Public-only check of bulk **CLDN4** against a **CD8 fraction** estimated from
the concordant-4 single-cell cohorts (GSE123902, GSE131907, GSE205335,
GSE189357), with the keratin-adjusted CLDN4–CD8A correlation on the same
TCGA-LUAD matrix as the comparator.

The fraction engine is Newman 2015 ν-SVR (the CIBERSORTx fractions core, not
the licensed S-mode container). The second engine is a simplex mixture MAP
with a flat Dirichlet prior (BayesPrism-style, not the BayesPrism Gibbs
sampler). Rules are in `PROTOCOL.md`. The run is in `FINDING.md`.

```bash
python3 methods/c4_bulk_deconv/build_reference.py
python3 methods/c4_bulk_deconv/deconvolve.py
```

Raw GEO/Xena matrices stay in `data/raw/` and are not committed.
