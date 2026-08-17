# Methods: GSE289287 T-47D Trop-2 KO xenograft prerank GSEA

**Additive only.** Public author DESeq2. Not lung. Not SKB264. CLDN4 padj ~0.47 is known.

## Scope

| In | Out |
|---|---|
| `GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz` | FASTQ / SRA / salmon |
| T-47D Trop-2 KO vs WT xenografts (4 vs 3) | DSG2 KO arms |
| Hallmark IFN-γ, MHC-I/APM, KEGG tight junction | Lung series; SKB264 |

## Engine

Same prerank as `gsea_core.py` in `methods/cldn4_ko_gsea`: Subramanian 2005 weighted KS *p*=1, 1000 gene-set permutations, seed=42, min set size 8. Positive NES = enriched in Trop-2 KO.

Primary rank = author DESeq2 Wald `stat`. Sensitivity rank = `log2FoldChange`. BH-FDR inside the three headline sets.

Gene lists: `gene_sets.json` (subset of the `a8_sets.json` freeze used by the CLDN4-loss extra).

## Reproduce

```bash
pip install -r methods/gse289287_trop2ko_gsea/requirements.txt
python3 methods/gse289287_trop2ko_gsea/analyze.py
```

Outputs: `FINDING.md`, `tables/gsea_headline.tsv`, `figures/`.
