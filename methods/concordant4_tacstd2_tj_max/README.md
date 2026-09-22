# Concordant-4 TACSTD2-high vs low, max-effect enrichment

Malignant pseudobulk from the locked concordant-4 cohorts only
(GSE123902, GSE131907, GSE205335, GSE189357).

The question is which pre-specified contrast and ranking statistic puts a
tight-junction or apical-junction set as far up the ORA and GSEA list as the
data support. The grid and the selection rule are in `analyze_max.py`.
Counts and the Hallmark / KEGG / GO collection are read from
`../concordant4_tacstd2_malignant_deg/data/`. Added MSigDB sets are in
`data/msigdb_junction_adhesion.json`.

```bash
python3 methods/concordant4_tacstd2_tj_max/analyze_max.py
```

Write-up: `FINDING.md`. Sweep tables: `tables/gsea_sweep.tsv`, `tables/ora_sweep.tsv`.
