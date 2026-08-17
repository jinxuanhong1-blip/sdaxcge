# Methods: public CLDN4 / Cldn4-loss prerank GSEA

**Additive only.** User thesis (CLDN4 loss shares IFN / MHC-I; CLDN4-high is the TJ barrier) is taken as given. No existing slide is re-run or retracted. Not SKB264. GSE50927 is mouse whole lung, not lung cancer.

## Scope

| In | Out |
|---|---|
| Public processed matrices on GEO | FASTQ / SRA / salmon / DESeq2 from reads |
| GSE207704 T47D + MCF7 CLDN4-/- FPKM | Lung-cancer KO (none public) |
| GSE50927 author edgeR, naive whole lung | VILI-only contrasts as primary |
| GSE22493 series matrix if probes map to symbols | Invented probe maps |

## Engine

Same prerank as `scripts/scrna_pseudobulk_gsea_meta/gsea_core.py` (copied to `scripts/cldn4_ko_gsea/gsea_core.py`): Subramanian 2005 weighted KS *p*=1, 1000 gene-set permutations, seed=42, min set size 8. Positive NES = enriched after CLDN4 / Cldn4 loss.

Headline sets (BH-FDR within these five, per contrast):

1. `HALLMARK_INTERFERON_GAMMA_RESPONSE`
2. `HALLMARK_INTERFERON_ALPHA_RESPONSE`
3. `CUSTOM_MHC_I_ANTIGEN_PRESENTATION`
4. `KEGG_TIGHT_JUNCTION`
5. `GOBP_KERATINIZATION`

Human lists: `data/genesets/a8_sets.json` (same freeze as the scRNA pseudobulk extra). Mouse Hallmark IFN: MSigDB `mh.all.v2023.2.Mm`. Mouse MHC-I / TJ / keratinization: title-case + HLA→H2 / ERAP2→Lnpep.

## Reproduce

```bash
pip install -r requirements.txt
python3 scripts/cldn4_ko_gsea/analyze.py
```

Outputs live under `methods/cldn4_ko_gsea/` (`FINDING.md`, `tables/`, `figures/`). Processed GEO files are in `methods/cldn4_ko_gsea/data/`.
