# Methods: GSE253564 leftover CLDN4 Q4 vs Q1 GSEA

**Additive only.** CLDN4 only. Does not re-run the TACSTD2 keratin/TJ/EMT leftover in `a8_ici_gsea`. Does not re-audit the MPR/PFS leftover in `opus_geo_leftover`.

## Scope

| In | Out |
|---|---|
| GSE253564 pre-treatment FPKM | GSE248378 post-treatment residual tumours |
| CLDN4 quartile Q4 vs Q1 | TACSTD2 splits |
| Hallmark IFN-γ, MHC-I/APM, KEGG tight junction | Hallmark EMT / keratin (already in A8 leftover) |
| CLDN4 vs CD274 Spearman on all samples with both genes | FASTQ / SRA; invented MPR labels |

## Engine

Same prerank as `gsea_core.py` in `methods/cldn4_ko_gsea`: Subramanian 2005 weighted KS *p*=1, 1000 gene-set permutations, seed=42, min set size 8. Positive NES = enriched in **CLDN4 Q4** (high minus low). BH-FDR inside the three headline sets.

Primary rank = Welch *t* on log2(FPKM+1), Q4 vs Q1. Sensitivity rank = Spearman ρ of each gene vs continuous CLDN4 (uses every sample). CLDN4 is a KEGG TJ member; a TJ sensitivity drops CLDN4 from the rank.

## Honest n

GEO deposits **32** pre-treatment tumours. Quartile arms are counted after the split (ties at the 25th/75th percentile stay in the arm). Do not quote n=32 as the GSEA contrast n. Correlation n is the number of samples with finite CLDN4 and the paired gene.

## Reproduce

```bash
pip install -r methods/gse253564_cldn4_gsea/requirements.txt
python3 methods/gse253564_cldn4_gsea/analyze.py
```

Outputs: `FINDING.md`, `tables/gsea_headline.tsv`, `tables/correlations.tsv`, `figures/`.
