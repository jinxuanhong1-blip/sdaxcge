# GEO brute wave: NHEJ down, IFN up

Additive public screen. It does not replace the locked CosMx, concordant-4, GSE137244, TCGA, or TISMO calls, and it does not replace the GSE50927 no-VILI EdgeR table.

## Search

1. Record the strict GEO query: series text contains CLDN4 and one of DNA repair, NHEJ, PRKDC, STING, cGAS, interferon.
2. Record the brute token query (claudin / CLDN as well as CLDN4; STING1 and "cGAS" instead of the MeSH-expanded bare tokens).
3. Download queue, expression series only:
   - CLDN4 / claudin-4 text series
   - brute token hits
   - NHEJ, PRKDC, STING1 / TMEM173 / "stimulator of interferon genes" / "cGAS-STING"
   - cGAS or cyclic GMP-AMP in the title
   - "DNA repair" anywhere in the series record
4. `interferon` alone is counted and not downloaded. The interferon series that also carry a claudin token are already in the brute query.

## Open matrix

An open matrix is an FTP series matrix or a supplementary txt/csv/tsv/xlsx/zip table between 8 KB and 25 MB whose name looks like expression (counts, FPKM, TPM, normalized, RNA-seq). Raw CEL/tar, BAM, and 10x mtx bundles are skipped and listed.

## Score

Direction is CLDN4-high minus CLDN4-low.

- NHEJ core: PRKDC, XRCC5, XRCC6, XRCC4, LIG4, NHEJ1, DCLRE1C, PAXX. Need at least 5.
- IFN ISGs: ISG15, MX1, MX2, OAS1, OAS2, OAS3, IFIT1, IFIT2, IFIT3, IFI44, IFI44L, EIF2AK2, IRF7, STAT1, RSAD2, BST2, ISG20, CXCL10. Need at least 8.
- cGAS (MB21D1), STING1 (TMEM173), TBK1, and IRF3 are reported beside the score and are not inside the IFN mean.

Continuous series with at least 6 samples: Spearman of CLDN4 versus the mean gene-wise z-score. Positive composite is IFN z minus NHEJ z.

Paired titles (WT/control/overexpression versus KO/knockdown/silencing): mean log2 difference, high arm minus low arm, averaged across genes. A joint NHEJ-down IFN-up call is the sign pair, including when the NHEJ limb is near zero.

Two-channel log-ratio series: one contrast, oriented so the CLDN4 probe mean is the high-minus-low direction. No Spearman on replicate ratios.

A series enters the Spearman census only when coverage, sample size, CLDN4 variance, and a single-channel design all pass. BH-FDR is computed inside that census.
