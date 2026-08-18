# Concordant-4 UMAP split by unit CLDN4 quartile

ADDITIVE visualization only. **CLDN4-only.**

Splits the same Harmony UMAP from the four-set atlas
(GSE123902 + GSE131907 + GSE205335 + GSE189357) by the unit-level
malignant CLDN4 %pos quartile (PR #503 rule: within-cohort rank then Q1/Q4).

This figure visualizes the already-reported inverse association
(higher malignant CLDN4 ↔ fewer T/NK). It is not a new test.

```bash
python3 methods/concordant4_atlas_umap_annotate/download.py
python3 methods/concordant4_umap_cldn4_tnk_split/analyze.py
```
