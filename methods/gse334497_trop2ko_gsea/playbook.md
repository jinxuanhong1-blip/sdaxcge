# Methods: GSE334497 4T1 Trop2 KO prerank GSEA

**Additive only.** Not lung. Not SKB264. Not a CLDN4 KO. Cldn4 gene-level *p* = 0.32 is given and is not re-tested. This folder reports NES/FDR for four headline sets and Cldn4 log2FC.

## Scope

| In | Out |
|---|---|
| GSE334497 author-normalized counts, 4T1 Trop2 KO vs WT tumors, 5 vs 5 | FASTQ / SRA / DESeq2 from reads |
| Prerank GSEA: IFN-γ, MHC-I/APM, KEGG TJ, keratin | Re-audit of the Cldn4 gene-level *p* |
| Cldn4 log2FC as effect size | Lung models, SKB264, CLDN4-loss matrices |

## Engine

Same prerank as `scripts/gse334497_trop2ko_gsea/gsea_core.py` (copy of the scRNA-pseudobulk / CLDN4-KO extra): Subramanian 2005 weighted KS *p*=1, 1000 gene-set permutations, seed=42, min set size 8. Rank = Welch *t* on log2(norm+1), KO minus WT. Positive NES = enriched in Trop2 KO.

Headline sets (BH-FDR within these four):

1. `HALLMARK_INTERFERON_GAMMA_RESPONSE` (MSigDB mouse Hallmark `mh.all.v2023.2.Mm`)
2. `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` (human freeze → mouse; HLA→H2, ERAP2→Lnpep)
3. `KEGG_TIGHT_JUNCTION`
4. `GOBP_KERATINIZATION`

Human lists: `data/genesets/a8_sets.json`.

## Reproduce

```bash
pip install -r requirements.txt
python3 scripts/gse334497_trop2ko_gsea/analyze.py
```

Outputs live under `methods/gse334497_trop2ko_gsea/` (`FINDING.md`, `tables/`, `figures/`).
