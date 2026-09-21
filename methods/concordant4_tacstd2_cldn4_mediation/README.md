# Concordant-4: TACSTD2 → T/NK through CLDN4

Patient-level mechanism test on the locked 65 units (GSE123902, GSE131907, GSE205335, GSE189357).

Malignant TACSTD2 and CLDN4 are scored on the same gates as the locked CLDN4 % positive versus T/NK analysis. The grid asks how much of the TACSTD2 association attenuates after CLDN4, and what fraction of the linear coefficient is the product-of-coefficients path TACSTD2 → CLDN4 → T/NK.

```bash
bash methods/concordant4_tacstd2_cldn4_mediation/scripts/download.sh /tmp/geo_dl
python3 methods/concordant4_tacstd2_cldn4_mediation/scripts/extract_scores.py
python3 methods/concordant4_tacstd2_cldn4_mediation/scripts/analyze.py
```

GEO matrices are not stored in this repo. `results/unit_scores.tsv` is enough to rerun `analyze.py`.

Numbers are in `FINDING.md`.
