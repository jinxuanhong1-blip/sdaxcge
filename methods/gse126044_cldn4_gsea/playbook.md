# Methods: GSE126044 CLDN4 Q4 vs Q1 prerank GSEA (IFN / MHC / TJ)

**Additive only.** Public Cho 2020 anti-PD-1 NSCLC counts. Single-gene CLDN4 vs ESTIMATE ImmuneScore is already known (n=16, Spearman r=−0.524, p=0.037) and is not re-tested. This folder adds prerank GSEA on CLDN4 **Q4 vs Q1**.

## Scope

| In | Out |
|---|---|
| `GSE126044_counts.txt.gz` (GEO suppl., HGNC symbols) | FASTQ / SRA / salmon |
| CLDN4 Q4 vs Q1 (honest n after the cut) | ImmuneScore / ESTIMATE re-test |
| Hallmark IFN-γ, MHC-I/APM, KEGG tight junction | TACSTD2 median A8 rerun |
| Sensitivity: Spearman vs continuous CLDN4 (n=16) | Sample-permutation GSEA |

## Engine

Same prerank as `gsea_core.py` in `methods/cldn4_ko_gsea` / `scrna_pseudobulk_gsea_meta`: Subramanian 2005 weighted KS *p*=1, 1000 gene-set permutations, seed=42, min set size 8. Positive NES = enriched in **CLDN4 Q4**.

Primary rank = Welch *t* on log2(CPM+1), Q4 minus Q1. Sensitivity ranks = mean log2(CPM+1) difference (Q4−Q1) and Spearman ρ vs continuous CLDN4 (uses all 16). BH-FDR inside the three headline sets.

Gene lists: `gene_sets.json` (same freeze as `methods/gse289287_trop2ko_gsea`).

## Reproduce

```bash
pip install -r methods/gse126044_cldn4_gsea/requirements.txt
python3 methods/gse126044_cldn4_gsea/analyze.py
```

Outputs: `FINDING.md`, `tables/gsea_headline.tsv`, `figures/`.
