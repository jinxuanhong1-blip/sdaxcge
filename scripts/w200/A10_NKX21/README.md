# A10: NKX2-1 inverse with TACSTD2 / CLDN4 in public LUAD

Tests the **direction** of NKX2-1 (TTF-1) versus TACSTD2 (TROP2) and CLDN4
in open LUAD matrices. This is a **bulk correlation** check, not a ChIP or
knockdown experiment. Filters are not tuned to force an inverse.

```bash
python3 scripts/w200/A10_NKX21/download.py
python3 scripts/w200/A10_NKX21/analyze.py
```

Outputs land in `results/w200/A10_NKX21/`.

## Primary claim (pre-specified)

In public LUAD, NKX2-1 is **inversely** correlated with TACSTD2 and with CLDN4.

Primary statistic: Spearman ρ on continuous expression (not IHC).
Primary cohort: TCGA-LUAD primary tumors (Xena GDC STAR TPM, log2(TPM+1)).
Replication (same statistic, not used to pick the primary): CPTAC LUAD RNA,
CPTAC LUAD protein (if both proteins are quantified), DepMap 24Q4 LUAD cell lines.

BH-FDR is applied **within** the pre-specified primary pair list only.
Context genes (SFTPB, NAPSA, KRT5, CLDN7) are exploratory.

## What this cannot test

- Direct transcriptional repression (no public NKX2-1 ChIP-seq peak call at
  TACSTD2/CLDN4 is required or claimed here).
- ICI response (TCGA/CPTAC/DepMap are not ICI-outcome cohorts).
- Protein from RNA, or TTF-1 IHC from RNA.

## Honest rule

Do not drop a cohort because the sign is positive or null.
Do not switch to Pearson, residuals, or a subtype filter to manufacture an inverse.
