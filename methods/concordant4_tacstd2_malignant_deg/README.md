# Concordant-4 TACSTD2 Q4 vs Q1 → TJ / adhesion → T/NK

Paper-funnel public page. Split is **TACSTD2** (malignant pseudobulk), not CLDN4.

## Cohorts (locked)

GSE123902 + GSE131907 + GSE205335 + GSE189357 only.

## Pipeline

```bash
python3 methods/concordant4_tacstd2_malignant_deg/analyze.py
```

Writes `FINDING.md`, `tables/`, `figures/`.

## Design

1. Patient/donor/sample malignant UMI-sum → log2(TMM-CPM+1)
2. Within-cohort TACSTD2 Q4 vs Q1 DEG (OLS, cohort covariates)
3. Hypergeometric ORA + prerank GSEA (TACSTD2 held out of rank/query)
4. Enrichr on UP genes (GO BP / KEGG / Hallmark)
5. Broad TJ score and 7-gene TJ core vs same-unit T/NK (DL meta)

Honest DE n is the quartile tails in the count matrix, not n=65.
