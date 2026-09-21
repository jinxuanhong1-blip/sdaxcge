# Concordant-4 CLDN4-high malignant markers

Public only: GSE123902, GSE131907, GSE205335, GSE189357.
Not GSE148071, GSE127465, GSE154826, GSE200563, or E-MTAB-13526.

COMET-style XL-mHG (xlmhg 2.5.4) and within-unit Q4 versus Q1 deltas for genes
that mark CLDN4-high malignant cells, then Enrichr through gseapy. A marker has
to agree in all four cohorts. Gene lists are under `results/gene_lists/`.

```bash
bash methods/concordant4_cldn4_comet_modules/scripts/run_all.sh /tmp/concordant4_raw /tmp/concordant4_malig
```

See `METHODS.md` for the rules and `FINDING.md` for the numbers from the run.
