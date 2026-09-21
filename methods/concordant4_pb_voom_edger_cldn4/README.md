# Concordant-4 malignant CLDN4 Q4 vs Q1, patient pseudo-bulk

Additive methods pass on the locked four cohorts only:

GSE123902, GSE131907, GSE205335, GSE189357.

Not GSE148071, GSE127465, GSE207422, GSE154826, or GSE200563.
Not a TACSTD2×CLDN4 dual-high score. Not a T/NK re-audit.

The malignant count matrices are already unit-level UMI sums (the muscat `pbDS` collapse). Differential expression uses the two engines `pbDS` calls:

- limma-voom with robust empirical Bayes
- edgeR quasi-likelihood (`glmQLFit`, robust)

Quartile labels are the within-cohort malignant CLDN4 percent-positive split from the concordant-4 analysis. They are not re-cut on the DE subset. P4001 is Q1 on the 22-patient vector and is absent from the UMI-sum.

Pooled layers:

- stacked voom / edgeR with a cohort intercept
- DerSimonian–Laird random-effects meta of the cohort log fold-changes (Knapp–Hartung on the family scores)
- family-score linear mixed model with a random intercept for cohort

Families: IFN (Hallmark α ∪ γ), MHC-I/APM, chemokine, tight junction (CLDN4 held out), keratin.

```bash
bash methods/concordant4_pb_voom_edger_cldn4/run.sh
```

Writeup: `RESULTS.md` at the repo root.
