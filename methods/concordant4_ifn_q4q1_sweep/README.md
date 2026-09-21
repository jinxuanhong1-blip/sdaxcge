# Concordant-4 Q4 vs Q1 IFN / MHC / chemokine sweep

Malignant pseudobulks from the locked concordant four
(GSE123902, GSE131907, GSE205335, GSE189357). Quartiles are the locked
within-cohort CLDN4 %pos labels. The contrast is Q4 versus Q1.

The locked comparison is the cohort-adjusted OLS t on zero-filled
log2(TMM-CPM+1), the ranking behind Hallmark IFN-gamma NES about −3.85.
This script reproduces that enrichment score and then sweeps pre-specified
rank statistics and interferon, MHC, and chemokine gene sets. Dropping
GSE205335 is a separate scope: counts are re-filtered and TMM is re-fit.

```bash
python3 methods/concordant4_ifn_q4q1_sweep/analyze.py
```

Needs numpy, scipy, and matplotlib. No network calls at runtime.

Outputs:

- `FINDING.md` — reproduction, most negative IFN NES, MHC and chemokine
  maxima, and signature logFC, each with the GSE205335-excluded number.
- `results/tables/nes_sweep.tsv` — the full grid.
- `results/tables/signature_logfc.tsv` — mean-logCPM OLS beta per set.
- `results/tables/gene_logfc.tsv` — canonical and IFN-set genes.
- `results/figures/` — forests and the canonical-gene logFC plot.

Headline IFN sets have 15–500 genes in the ranked universe. The grid
maximum is descriptive. The OLS-t Hallmark IFN-gamma row is the comparison
to the prior fgsea result.
