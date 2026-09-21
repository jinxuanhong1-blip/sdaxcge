# Concordant-4: TACSTD2, CLDN4, and IFN/MHC

Malignant pseudobulk from the locked concordant four (GSE123902, GSE131907, GSE205335, GSE189357).

On these units TACSTD2-high is IFN-enriched. The continuous pseudobulk IFN-γ NES stays positive after CLDN4 is in the model (+2.32, patient interval +0.74 to +2.96; residualized +2.36, +1.04 to +3.18). The CLDN4-low half does not show a stable IFN decrease. Percent-positive mediation has a negative TACSTD2 → CLDN4 → IFN product and a positive direct slope. Writeup: `FINDING.md`.

Expression n = 64. P4001 is in the percent-positive vector and out of the UMI sum. The gene universe is the intersection across the four cohorts, then log2(TMM-CPM+1). That intersection is what reproduces the PR #503 IFN family score (logFC −0.584, 18 vs 16).

```bash
python3 methods/concordant4_tacstd2_cldn4_ifn_mediation/analyze.py
```

Writeup: `FINDING.md`. Tables: `results/tables/`. Figures: `results/figures/`.
