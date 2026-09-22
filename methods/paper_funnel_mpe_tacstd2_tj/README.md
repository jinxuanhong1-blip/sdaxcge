# PAPER FUNNEL: Human MPE / GSE131907 — TACSTD2-high vs immune, TJ, CLDN4

Corroborate PPT middle-slide claims with public data only. **Never fabricate.**

HRA006761 (Clin Transl Med 2024 recurrent-MPE CLDN4 study) is GSA-Human **controlled access** (DAC HDAC002197). This run uses public **GSE131907** pleural effusion + primary tumor (Kim et al., Nat Commun 2020).

## Claims tested

1. TACSTD2 is high in malignant-like cells and low in immune cells  
2. Tight-junction (TJ) program is enriched in TACSTD2-high malignant-like cells  
3. CLDN4 is co-detected / coexpressed with TACSTD2  

## Run

```bash
bash methods/paper_funnel_mpe_tacstd2_tj/scripts/download.sh /tmp/gse131907
python3 methods/paper_funnel_mpe_tacstd2_tj/scripts/extract_genes.py \
  --datadir /tmp/gse131907 --outdir /tmp/gse131907/paper_funnel_mpe \
  --tj-list methods/paper_funnel_mpe_tacstd2_tj/data/tj_core_genes.txt
python3 methods/paper_funnel_mpe_tacstd2_tj/scripts/analyze.py \
  --datadir /tmp/gse131907 --subset /tmp/gse131907/paper_funnel_mpe \
  --outdir methods/paper_funnel_mpe_tacstd2_tj/results
```

Numbers and verdicts: [RESULTS.md](RESULTS.md). Raw tables under `results/tables/`.

## Notes

- Author `Cell_subtype` “Malignant cells” is **0** in PE and tLung. MPE “malignant-like” = PE epithelial EPCAM+ WT1− CALB2−.  
- Related: PR #635 (CLDN4-first MPE coexpression), PR #230 (atlas epithelial vs immune + TJ GSEA).  
- Does not touch locked CosMx / concordant-4 / private KL results.
