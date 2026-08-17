# Playbook — GSE166449 CLDN4 Q4 vs Q1 prerank GSEA

Additive extra. **CLDN4 only.** Honest **n = 22** (Q4 vs Q1 is 6 vs 6). Does not audit or retract `results/w200/GSE166449`.

```bash
python3 methods/gse166449_cldn4_gsea/download.py
python3 methods/gse166449_cldn4_gsea/analyze.py
```

| Locked choice | Value |
|---|---|
| Cohort | GSE166449 pretreatment LUAD pembrolizumab bulk, n=22 |
| Anchor | CLDN4 (not TACSTD2) |
| Contrast | Q4 vs Q1 on deposited log2(TPM+1) CLDN4 |
| Rank | Welch *t*, Q4 minus Q1 |
| Comparator | CD274 Spearman vs CLDN4; CD274 Q4 vs Q1 GSEA on the same three sets |
| Engine | `gsea_core.py` weighted KS p=1, 1000 gene-set permutations, seed=42 |
| Headline | Hallmark IFN-γ, MHC-I/APM (21), KEGG tight junction |
| FDR | BH inside the three headline sets, per rank |
| Honesty | n=22 / 6 vs 6; CLDN4 is in KEGG TJ; also report NES after dropping CLDN4 |

Outputs: `FINDING.md`, `tables/gsea_headline.tsv`.
