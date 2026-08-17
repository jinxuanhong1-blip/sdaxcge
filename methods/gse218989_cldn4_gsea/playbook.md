# Methods: GSE218989 CLDN4 Q4 vs Q1 prerank GSEA

**Additive only.** Public SMC–KAIST NSCLC ICI bulk TPM. Does not audit or retract `methods/gse218989_cldn4_ici`. Response vs CLDN4 is already known NS there (GEO MWU p=0.40, AUC=0.47) — do not headline an odds ratio.

## Scope

| In | Out |
|---|---|
| `GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz` | FASTQ / SRA |
| CLDN4 Q4 vs Q1 prerank (Welch *t*) | Invented histology (LUAD/LUSC not deposited) |
| Hallmark IFN-γ, MHC-I/APM, KEGG tight junction | Response OR as a headline |
| Honest NES / FDR | Spinning a null IFN/MHC result |

## Engine

Same prerank as `gsea_core.py` in `methods/cldn4_ko_gsea` / `gse289287_trop2ko_gsea`: Subramanian 2005 weighted KS *p*=1, 1000 gene-set permutations, seed=42, min set size 8. Positive NES = enriched in **CLDN4 Q4**. BH-FDR inside the three headline sets.

Primary rank = Welch *t* on `log2(TPM+1)`, Q4 minus Q1. Sensitivity ranks = Spearman ρ vs continuous CLDN4, and Welch *t* after dropping CLDN4 from the rank (CLDN4 is a KEGG TJ member; the split is on CLDN4).

## Reproduce

```bash
pip install -r methods/gse218989_cldn4_gsea/requirements.txt
python3 methods/gse218989_cldn4_gsea/analyze.py
```

Outputs: `FINDING.md`, `tables/gsea_headline.tsv`, `figures/`.
