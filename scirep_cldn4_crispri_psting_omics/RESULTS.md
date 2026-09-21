# Sci Rep 2025 CLDN4 CRISPRi / pSTING: TCGA supplement reanalysis

**One CLDN4-low row is a positive ISG association.** In Supplementary Table 2, IFI44-high co-occurs with CLDN4-low: log2 odds ratio 1.571, p=0.036, q=0.050. Tendency is printed “Co-ocurrence”. That is the only row in either open TCGA table that is CLDN4-low and sits in the paper’s own ISG or cGAS–STING gene list with a positive direction.

No q-value or odds ratio was recomputed. `<-3` was left censored. Six of the ten named ISGs are absent from the tables, and every named cGAS–STING gene is absent. Nothing was filled in for those genes.

Paper: Villagomez et al., *Scientific Reports* 2025;15:39257. DOI [10.1038/s41598-025-23137-1](https://doi.org/10.1038/s41598-025-23137-1). PMID 41214101. PMC12603150.

The ISG list is the ten genes named in the Results: ISG15, MX1, OAS1, IFIT1, EIF2AK2, IRF7, STAT1, RSAD2, BST2, IFI44. The cGAS–STING list is the set named in that same section: cGAS, STING, LC3, TBK1, IRF3, Beclin-1, and Rab7. Aliases searched are in `reanalysis_named_genes_absent.tsv`. Ranking is `analyze_tcga_tables.py` on the reprinted tables.

---

## CLDN4-low and the ISG table (MOESM2, n=617)

cBioPortal, TCGA ovarian serous cystadenocarcinoma, 617 samples. Five printed rows. The caption stops at “CLDN4 (Claudin-4)”.

| Rank by printed q | Genes as printed | Log2 odds | p | q | Direction | In the paper’s ISG list |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | IRF7 (high), CLDN4(high) | 1.991 | 0.004 | 0.007 | Co-ocurrence | yes, CLDN4-high |
| — | CLDN4 (high), CLDN4 (low) | <-3 | 0.010 | 0.016 | Mutual exclusivity | self row, not an ISG |
| 2 | BST2 (high), CLDN4(low) | <-3 | 0.029 | 0.042 | Mutual exclusivity | yes, CLDN4-low, negative direction |
| 3 | OAS1 (high), CLDN4 (high) | 1.194 | 0.035 | 0.049 | Co-occurrence | yes, CLDN4-high |
| 4 | IFI44(high), CLDN4(low) | 1.571 | 0.036 | 0.050 | Co-ocurrence | yes, CLDN4-low, positive direction |

Among the two CLDN4-low ISG rows, BST2 has the smaller q. Its log2 odds is censored at `<-3`, so the magnitude is larger than 3 and larger than IFI44’s 1.571, and the printed direction is mutual exclusivity: BST2-high does not co-occur with CLDN4-low. IFI44 is the CLDN4-low row whose printed direction is co-occurrence.

IRF7 with CLDN4-high is the strongest ISG row in the table (q=0.007). OAS1 is the other CLDN4-high co-occurrence (q=0.049). Those two are the associations the article text highlights. They are not CLDN4-low.

Absent from this table, with no result added: ISG15, MX1, IFIT1, EIF2AK2, STAT1, RSAD2.

Sign check: every co-occurrence has a positive log2 odds ratio, and every mutual-exclusivity row is negative or censored `<-3` (5/5).

---

## CLDN4-low proteins (MOESM5, 18 rows)

TCGA PanCancer Atlas ovarian serous cystadenocarcinoma. Proteins with a Benjamini–Hochberg q from a Student t-test, as the caption states. CLDN4 mRNA groups are z-score >1 and z-score <−1. Supplementary Figure 21 prints the group sizes as n=101/100 (high) and n=73/72 (low). Those slashes are quoted as printed. They were not used to recompute a test.

The Word header repeats `CLDN4: EXP>1` and `CLDN4: EXP<-1`. Those four columns stay unlabeled. The ranking uses the columns that have names: Log2 Ratio, p-Value, q-Value, and Higher expression in.

Seven proteins are higher when CLDN4 mRNA is low. None of them is in the ISG list or the cGAS–STING list. Within this table, order by q, then p, then absolute log2 ratio:

| Rank | Protein | Log2 ratio | p | q |
| --- | --- | --- | --- | --- |
| 1 | COL6A1 | -0.42 | 1.22E-04 | 8.77E-03 |
| 2 | YWHAE | -0.11 | 1.01E-03 | 0.0241 |
| 3 | YAP1 | -0.24 | 1.46E-03 | 0.0314 |
| 4 | CDH2 | -0.2 | 2.54E-03 | 0.0389 |
| 5 | MS4A1 | -0.1 | 3.14E-03 | 0.0389 |
| 6 | CHEK1 | -0.16 | 3.22E-03 | 0.0389 |
| 7 | YAP1_PS127 | -0.46 | 3.24E-03 | 0.0389 |

YAP1_PS127 has the largest absolute log2 ratio among CLDN4-low proteins (-0.46) and the weakest q in that set (0.0389, last by p among the four-way q tie). COL6A1 is first by q. CHEK1 is a DNA-damage checkpoint protein in this low set. It is not one of the cGAS–STING genes the paper names, so it is not counted as a STING hit.

The other eleven proteins are higher when CLDN4 mRNA is high, including MTOR (log2 ratio 0.31, q=0.0389) and MTOR_PS2448 (log2 ratio 0.2, q=0.0389). That matches the Figure 21 caption, which points the mTOR annotation at the CLDN4-high group. RAB25 is in the high group. Rab7 is not in the table. STAT3_PY705 is higher in the CLDN4-high group (log2 ratio 0.22, p=3.84E-04, q=0.0138). STAT3 is not one of the ten named ISGs.

STING, pSTING, cGAS, TBK1, IRF3, LC3, Beclin-1, and Rab7 do not appear as rows. Sign check: log2 ratio is negative for all 7 CLDN4-low proteins and positive for all 11 CLDN4-high proteins (18/18).

Protein q-values and ISG q-values are separate Benjamini–Hochberg families. They are not pooled into one ranking.

---

## Extra supplements scraped

Nature lists five supplementary materials. MOESM6 and a reporting-summary filename are not present (HTTP 403 on the Springer object URL). The Europe PMC package for PMC12603150 has 21 files. The five MOESM hashes match the Springer downloads. The other 16 files are HTML renders of main Figures 1–8, not tables.

MOESM1 (2,327,040 bytes) is nine pages of uncropped membranes. MOESM3 is the flow-antibody list (36 rows), not a proteome. MOESM4 is legends for Supplementary Figures 1–21 plus 21 PNGs. Those legends contain no second ISG or protein matrix.

---

## Raw RNA-seq and proteomics

Still **0 open raw RNA-seq files and 0 open raw proteomics files** for this DOI. GEO, SRA, PRIDE, MassIVE, BioStudies, OmicsDI, and Figshare have no deposit. The data-availability line is author-on-request. pSTING in the CRISPRi experiments is confocal microscopy and immunoblot. Queries and the negative repository screen are in `search_log.tsv`.

No new TCGA matrix was downloaded. The seven-cohort keratin screen is untouched. Private 8-KL, HIS FCS, and PDX proteome files were not used.
