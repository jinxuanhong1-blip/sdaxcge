# RESULTS — public TACSTD2 knockdown and sacituzumab RNA-seq

Question: in open RNA-seq, do CLDN4, classical NHEJ, cGAS–STING, and IFN/ISG genes move after TACSTD2/Trop2 knockdown or sacituzumab govitecan (IMMU-132)?

Cell-intrinsic answers come from three lines. SKOV3 shTACSTD2 lowers CLDN4 (and CLDN3/CLDN7), lowers STING1 and a pre-specified ISG set, and raises XRCC5/XRCC6/PRKDC/PAXX inside an otherwise split NHEJ list. Acute sacituzumab moves CLDN4 in both cell lines, down in CX-1 and up in KYSE30, and raises STING1 in both. The 8-gene NHEJ mean and the 14-gene ISG mean do not share one direction across those SG arms. Bulk 4T1 Trop2-KO tumors and day-29 CRC PDXs are reported separately because immune cells are in the RNA.

## 1. What was scored

Depositor-processed matrices. Expression entered as log2(x+1). Unpaired contrasts use a two-sided Welch t-test; Mann–Whitney is stored beside it. With n=3 vs 3 the two-sided Mann–Whitney p cannot fall below 0.1, so the p quoted below is Welch. Benjamini–Hochberg q is computed inside the 31-gene panel of that contrast, not across the genome. The CRC PDX contrast is paired (one-sample t and Wilcoxon on within-model deltas).

| contrast | accession | n | matrix | role |
|---|---|---|---|---|
| SKOV3 shTACSTD2 vs shNC | GSE245459 | 3 vs 3 | FPKM | primary knockdown |
| SKOV3 shTACSTD2+cisplatin vs shNC+cisplatin | GSE245459 | 3 vs 3 | FPKM | secondary |
| 4T1 Trop2 KO vs WT tumors | GSE334497 | 5 vs 5 | normalized counts | primary knockout, bulk tumor |
| CX-1 IMMU132 vs control, 2 days | GSE312098 | 3 vs 3 | FPKM | primary ADC, cell line |
| KYSE30 IMMU132 vs control, 1 day | GSE304294 | 2 vs 3 | FPKM | primary ADC, cell line, SG arm n=2 |
| CRC PDX IMMU132 vs matched control, day 29 | GSE311016 | 5 pairs | FPKM | primary ADC, bulk tumor |

Sets, fixed before scoring (`methods/trop2_kd_sg_rnaseq/gene_sets.py`):

- CLDN4 alone. CLDN3 and CLDN7 are context, not the set.
- NHEJ core: XRCC6, XRCC5, PRKDC, LIG4, XRCC4, NHEJ1, DCLRE1C, PAXX.
- STING axis: CGAS, STING1, TBK1, IRF3.
- IFN/ISG: STAT1, STAT2, IRF7, IRF9, ISG15, MX1, OAS1, OAS2, IFIT1, IFIT3, CXCL10, RSAD2, IFI44, BST2.
- Checks, not part of the asked sets: TACSTD2 (did the knockdown land) and CDKN1A (did the SN-38 payload register).

In the 4T1 matrix, OAS1 is the mouse gene Oas1a (ENSMUSG00000052776).

## 2. Knockdown checks

SKOV3 shTACSTD2 vs shNC: TACSTD2 log2FC = −2.63 (control mean log2(FPKM+1) = 2.68, knockdown = 0.05), Welch p = 5.3×10⁻⁵, panel q = 5.4×10⁻⁴.

4T1 Trop2 KO vs WT: Tacstd2 log2FC = −3.82 (7.07 → 3.25), p = 0.0011, q = 0.033. The knockout is visible in bulk tumor RNA.

Sacituzumab does not knock the transcript down. CX-1 TACSTD2 log2FC = +0.68 (p = 5.7×10⁻⁴). KYSE30 TACSTD2 log2FC = +1.16 (p = 0.0099).

CDKN1A, the acute payload check: CX-1 +1.17 (p = 0.0025). KYSE30 +1.54 (p = 0.066, n=2). Day-29 PDX CDKN1A log2FC = −0.51 (paired p = 0.19). The PDX profile is residual tumor after three weeks, not the day-1/day-2 damage snapshot.

## 3. CLDN4

| contrast | log2FC | Welch or paired p | panel q | other claudins |
|---|---:|---:|---:|---|
| SKOV3 shTACSTD2 | −1.92 | 3.3×10⁻⁵ | 5.1×10⁻⁴ | CLDN3 −1.48 (p=1.1×10⁻⁴); CLDN7 −1.39 (p=0.0010) |
| 4T1 Trop2 KO, bulk | −0.82 | 0.25 | 0.64 | CLDN3 −1.44 (p=0.082); CLDN7 −0.82 (p=0.091) |
| CX-1 IMMU132, 2 days | −0.86 | 1.7×10⁻⁵ | 5.2×10⁻⁴ | CLDN3 −0.70 (p=0.0059); CLDN7 −0.14 (p=0.064) |
| KYSE30 IMMU132, 1 day | +0.91 | 3.9×10⁻⁵ | 0.0010 | CLDN3 +1.46 (p=0.0026); CLDN7 +0.52 (p=8.8×10⁻⁴) |
| CRC PDX IMMU132, day 29 | −0.44 | 0.060 | 0.28 | 5/5 pairs negative (−0.94, −0.56, −0.04, −0.09, −0.55). CLDN7 −0.26 (p=0.033) |

SKOV3 control CLDN4 sits at log2(FPKM+1) = 2.05 and falls to 0.14, so the drop is not a floor artifact. The same knockdown under cisplatin only moves CLDN4 by −0.14 (p = 0.099), because shNC+cisplatin is already at 0.45.

CX-1 and KYSE30 are both acute cell-line SG contrasts with a CDKN1A rise, and CLDN4 moves in both, in opposite directions. The day-29 PDX mean is downward in every model, with a paired p of 0.060.

## 4. NHEJ

The 8-gene mean is not a single shift.

| contrast | mean log2FC | genes up / down | set Welch p |
|---|---:|---|---:|
| SKOV3 shTACSTD2 | +0.37 | 5 / 3 | 0.18 |
| 4T1 Trop2 KO | −0.04 | 3 / 5 | 0.58 |
| CX-1 IMMU132 | −0.09 | 3 / 5 | 0.29 |
| KYSE30 IMMU132 | −0.37 | 2 / 6 | 0.026 |
| CRC PDX IMMU132 | −0.16 | 4 / 4 | 0.24 |

Inside that average, SKOV3 knockdown raises the Ku/DNA-PKcs/PAXX group and leaves the ligase group flat: XRCC5 +0.98 (q=5.1×10⁻⁴), XRCC6 +0.80 (q=5.4×10⁻⁴), PAXX +0.88 (q=0.0010), PRKDC +0.67 (q=0.0026); LIG4 −0.23 (p=0.22), XRCC4 +0.05 (p=0.81).

KYSE30 SG (n=2) moves several of those same genes downward: PRKDC −0.80, XRCC4 −0.88, XRCC5 −0.71, XRCC6 −0.67 (each p<0.03), while LIG4 is +0.62. CX-1 SG is smaller and mixed (XRCC5 −0.35, PRKDC −0.31, LIG4 +0.42). 4T1 and the PDX stay within a few tenths on the set mean.

## 5. STING

| contrast | set mean log2FC | set p | STING1 log2FC | CGAS |
|---|---:|---:|---:|---|
| SKOV3 shTACSTD2 | −0.35 | 0.083 | −1.30 (q=0.0069) | both arms near 0 (log2 0.03 vs 0.01) |
| 4T1 Trop2 KO | +0.07 | 0.042 on the z-score, mean log2FC +0.07 | +0.31 (p=0.071) | −0.38 (p=0.10) |
| CX-1 IMMU132 | +0.35 | 0.017 | +0.50 (q=0.038) | +0.48 (q=0.024) |
| KYSE30 IMMU132 | +0.28 | 0.52 | +1.22 (q=0.0083) | −0.16 (p=0.27) |
| CRC PDX IMMU132 | +0.41 | 0.13 | +0.55 (p=0.32) | +0.17 (p=0.63) |

SKOV3 knockdown removes expressed STING1 (control log2 1.32 → 0.02). CGAS is already off, TBK1 is flat (+0.09), IRF3 moves −0.17. That is a STING1 drop, not a four-gene shutdown.

Both acute SG cell lines raise STING1. CX-1 also raises CGAS, so the four-gene mean reaches +0.35 (p=0.017). KYSE30’s STING1 rise does not carry CGAS, TBK1, or IRF3, and the set mean stays at +0.28 (p=0.52).

The 4T1 z-score p of 0.042 sits on a set mean of +0.07 log2. Gene shifts are a few tenths (IRF3 +0.32, p=0.018, q=0.28). Bulk RNA also mixes stroma and immune cells, so this is not read as tumor-cell STING activation.

## 6. IFN / ISG

| contrast | mean log2FC | genes up / down | set p |
|---|---:|---|---:|
| SKOV3 shTACSTD2 | −0.80 | 0 / 14 | 0.0013 |
| 4T1 Trop2 KO | +0.10 | 9 / 5 | 0.68 |
| CX-1 IMMU132 | +0.30 | 10 / 4 | 0.053 |
| KYSE30 IMMU132 | +0.13 | 8 / 6 | 0.16 |
| CRC PDX IMMU132 | +0.31 | 11 / 3 | 0.63 |

SKOV3 knockdown lowers expressed ISGs, not only genes sitting on the floor: IFIT1 −1.89 (control log2 4.06), OAS1 −1.46 (control 1.65), STAT2 −1.11, ISG15 −0.62, STAT1 −0.49, BST2 −0.50. All 14 pre-specified genes go down.

CX-1 SG is a partial ISG rise: ISG15 +0.63 (p=0.0033), IFIT1 +0.49 (p=0.0065), OAS1 +0.48 (p=0.025). The 14-gene mean is +0.30, p=0.053. KYSE30 splits (OAS1 +1.21, IFIT1 −0.61), set p=0.16. PDX pair deltas for ISG15 run from +4.93 to −1.81, so the mean of +0.31 (p=0.63) is not a shared induction. 4T1 IFN is flat in bulk RNA.

Under cisplatin, the SKOV3 knockdown ISG mean flips to +0.41 (p=0.016). That arm is secondary; the untreated knockdown is the ISG decrease.

## 7. Sets that were inventoried and not scored

Full table: `tables/inventory.tsv`.

- E-MTAB-16433, E-MTAB-16843, E-MTAB-16849 are real sacituzumab single-cell RNA-seq (CRC PDOX; organoid SG vs IgG1-SN-38 time course; liver-metastasis models). Processed counts are public (about 0.93 GB, 5.9 GB, and 5.0 GB). The EBI file host did not complete TLS from this run, so those matrices were not scored.
- GSE15212 is SW480 siTACSTD2 on GPL4133 microarray, not RNA-seq.
- GSE278664 discusses sacituzumab in the series text. The matrix is pre-treatment HGSOC biopsies (15 BRCAmut, 20 BRCAwt) from NCT02203513, not an SG-versus-control contrast.
- GSE309617 / GSE309616 are carboplatin-resistant TNBC PDXs. GSE303323 and related accessions are KRAS-MAPK inhibitor profiles. GSE292860 is Q901 ± topotecan. GSE302284 is osimertinib-versus-vehicle scRNA from a TROP2 CAR-T paper. GSE235812 is tumor-versus-PDX TACSTD2 correlation. None of these is a TACSTD2 knockdown or sacituzumab treatment arm.

## 8. Methods

GEO supplementary matrices were read from the NCBI FTP. Sample columns were matched to GEO library names: SKOV3 `sh1–3` vs `shNC1–3`; 4T1 KO libraries KO162, KO164, KO165, KO172, RESUB-KO163R versus WT RESUB-171R, RESUB-170R, RESUB-169R, RESUB-168R, control170; CX-1 X_4–X_6 vs X_1–X_3; KYSE30 OX2 vs OX1; PDX T_{36,82,83,114,196} vs C_{same}. GSE312098 and GSE311016 files are UTF-16. Where a symbol had more than one protein-coding row, the row with the higher mean was kept. Set scores in the tables are the mean of per-gene z-scores across the samples in that contrast; the mean log2FC above is the unweighted mean of member log2FCs and is the number used in the text. Combination arms (IMMU132+GSK2606414; IMMU132+IACS-010759) are in `tables/set_scores.tsv` and were not used as the sacituzumab monotherapy result.

## 9. Figures

- `figures/fig1_log2fc_heatmap.png` — gene-level log2FC for every contrast, including the two combination arms.
- `figures/fig2_set_mean_log2fc.png` — CLDN4, NHEJ, STING, and IFN mean log2FC.

## 10. What a paper can use

In SKOV3 RNA-seq (GSE245459, n=3 vs 3), shTACSTD2 drops TACSTD2 (log2FC −2.63, p=5.3×10⁻⁵) and CLDN4 (−1.92, p=3.3×10⁻⁵), with CLDN3 and CLDN7 alongside it, and it lowers a 14-gene ISG set (mean −0.80, 14/14 down, p=0.0013) and STING1 (−1.30, q=0.0069). XRCC5, XRCC6, PRKDC, and PAXX rise (log2FC +0.67 to +0.98); the full 8-gene NHEJ mean stays at +0.37 (p=0.18). In 4T1 bulk tumors (GSE334497, n=5 vs 5) Tacstd2 falls (log2FC −3.82, p=0.0011) and CLDN4’s mean is lower (−0.82, p=0.25). Acute sacituzumab (CDKN1A up) lowers CLDN4 in CX-1 at 2 days (−0.86, p=1.7×10⁻⁵; GSE312098) and raises it in KYSE30 at 1 day (+0.91, p=3.9×10⁻⁵; GSE304294, SG n=2). STING1 rises in both of those cell lines (+0.50 and +1.22). The ISG set is partial in CX-1 (mean +0.30, p=0.053) and mixed in KYSE30 (mean +0.13, p=0.16). Day-29 CRC PDXs (GSE311016, 5 pairs) have lower CLDN4 in 5/5 models (mean −0.44, paired p=0.060) without a shared NHEJ or IFN shift.
