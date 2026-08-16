# A10: KLF4 vs TACSTD2 / CLDN4 in public lung

Tests the **direction** of KLF4 versus TACSTD2 (TROP2) and CLDN4
in open lung matrices. This is a **bulk / cell-line correlation** check,
not a ChIP or knockdown experiment. Filters are not tuned to force a sign.

```bash
python3 scripts/w200/A10_KLF4/download.py
python3 scripts/w200/A10_KLF4/analyze.py
```

Outputs land in `results/w200/A10_KLF4/`.

## Primary question (pre-specified)

In public lung, does KLF4 associate with TACSTD2 and with CLDN4 **in the
same direction**?

Primary statistic: Spearman ρ on continuous expression.
Primary cohort: TCGA-LUAD primary tumors (Xena GDC STAR TPM, log2(TPM+1)).
Pre-specified lung second histology: TCGA-LUSC (same statistic; does not
replace the LUAD verdict).
Replication: CPTAC LUAD and LSCC RNA + protein, DepMap 24Q4 lung cell lines.

BH-FDR is applied **within** the pre-specified KLF4 vs TACSTD2/CLDN4 rows only.
Context genes (KRT5, CLDN7, PECAM1) and the TACSTD2–CLDN4 pair are exploratory.

## What this cannot test

- Direct KLF4 binding or transcriptional control at TACSTD2/CLDN4.
- ICI response (TCGA/CPTAC/DepMap are not ICI-outcome cohorts).
- Protein from RNA, or IHC from RNA.

## Honest rule

Do not drop a cohort because the sign is inconvenient.
Do not switch to Pearson, residuals, or a subtype filter to manufacture
concordance. Report LUAD and LUSC even if they disagree.
