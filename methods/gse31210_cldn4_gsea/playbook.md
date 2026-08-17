# Playbook — GSE31210 CLDN4 Q4 vs Q1 prerank GSEA

Additive extra. **CLDN4 only.** CD8 ρ is already known (`methods/gse31210_cldn4_immune`). Does not audit or retract any slide.

```bash
python3 methods/gse31210_cldn4_gsea/download.py
python3 methods/gse31210_cldn4_gsea/analyze.py
```

| Locked choice | Value |
|---|---|
| Cohort | GSE31210 Japanese stage I–II LUAD, primary tumors, honest n=226 |
| Platform | GPL570 U133 Plus 2.0; unique-mapped max-mean collapse |
| Transform | `log2(MAS5+1)` |
| Anchor | CLDN4 (not TACSTD2; not CD8) |
| Contrast | Q4 vs Q1 on CLDN4 |
| Rank | Welch *t*, Q4 minus Q1 |
| Engine | `gsea_core.py` weighted KS p=1, 1000 gene-set permutations, seed=42 |
| Headline | Hallmark IFN-γ, MHC-I/APM (21), KEGG tight junction |
| FDR | BH inside the three headline sets, per rank |
| Honesty | CLDN4 is in KEGG TJ; also report NES after dropping CLDN4 from the rank |

Outputs: `FINDING.md`, `tables/gsea_headline.tsv`.
