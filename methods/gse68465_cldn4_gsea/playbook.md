# Playbook — GSE68465 CLDN4 Q4 vs Q1 prerank GSEA (IFN / MHC / TJ)

Additive extra. **CLDN4 only.** Does not audit or retract `methods/gse68465_cldn4` (CD8 / ESTIMATE purity).

```bash
python3 methods/gse68465_cldn4_gsea/download.py
python3 methods/gse68465_cldn4_gsea/analyze.py
```

| Locked choice | Value |
|---|---|
| Cohort | GSE68465 Director's Challenge LUAD microarray |
| Honest n | GEO 462 arrays; **443** LUAD; 19 Normal held out; GSEA contrast = Q4 vs Q1 of the 443 |
| Anchor | CLDN4 (not TACSTD2) |
| Matrix | Public series matrix MAS5. No CEL reprocess. |
| Transform | `log2(MAS5)` |
| Collapse | max-mean probe → first `///` HUGO from official GPL96 |
| Contrast | Q4 vs Q1 on log2 CLDN4 (ties stay in the tail) |
| Rank | Welch *t*, Q4 minus Q1 |
| Engine | `gsea_core.py` weighted KS p=1, 1000 gene-set permutations, seed=42 |
| Headline | Hallmark IFN-γ, MHC-I/APM (21), KEGG tight junction |
| FDR | BH inside the three headline sets, per rank |
| Honesty | Report Q4 vs Q1 n, not 443, as the GSEA n. CLDN4 is in KEGG TJ; also drop CLDN4 from the rank. U133A coverage is incomplete; n in rank is the intersection. |

Outputs: `FINDING.md`, `tables/gsea_headline.tsv`.
