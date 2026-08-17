# Playbook — GSE285029 CLDN4 Q4 vs Q1 prerank GSEA

Additive extra. **CLDN4 only.** Does not audit or retract any slide.

```bash
python3 methods/gse285029_cldn4_gsea/download.py
python3 methods/gse285029_cldn4_gsea/analyze.py
```

| Locked choice | Value |
|---|---|
| Cohort | GSE285029 pre-ICI NSCLC WTS, n=234 |
| Anchor | CLDN4 (not TACSTD2) |
| Contrast | Q4 vs Q1 on log2(clip0+1) CLDN4 |
| Rank | Welch *t*, Q4 minus Q1 |
| Engine | `gsea_core.py` weighted KS p=1, 1000 gene-set permutations, seed=42 |
| Headline | Hallmark IFN-γ, MHC-I/APM (21), KEGG tight junction |
| FDR | BH inside the three headline sets, per rank |
| Honesty | CLDN4 is in KEGG TJ; also report NES after dropping CLDN4 from the rank |

Outputs: `FINDING.md`, `tables/gsea_headline.tsv`.
