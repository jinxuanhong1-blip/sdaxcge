# CLDN4-high malignant markers, concordant-4

ADDITIVE. CLDN4 only. Not a re-estimate of the locked patient-level T/NK correlation
(ρ = −0.531, n = 65). These lists are genes that move with the CLDN4-high versus
CLDN4-low split inside malignant cells, in every one of GSE123902, GSE131907,
GSE205335, and GSE189357. They are co-markers of that state. They are not a
surface-target ranking and they do not pin CLDN4 over EPCAM, CLDN7, or MUC1.

The split is within-unit malignant Q4 versus Q1 of log1p(CP10k CLDN4).
Units need ≥40 malignant cells and ≥10 CLDN4-positive malignant cells.
P4001-style thin units and units with essentially no CLDN4 range are out.
Patient / donor / sample remains the unit for the mean Δ. Cell counts below are
the cells that entered the COMET ranking, not n.

## Units used

61 of 65 locked units entered the contrast. Four were out on the pre-set gates, not on the result:

- GSE123902: 12 / 13. Out: LX699 (7 CLDN4-positive malignant cells).
- GSE131907: 19 / 21. Out: EBUS_13 (1 CLDN4-positive malignant cell), NS_16 (9).
- GSE205335: 21 / 22. Out: P4001 (27 malignant cells).
- GSE189357: 9 / 9.

11,890 symbols are present in all four cohorts after the GSE123902 intersection. Of those, 8,569 have a positive patient-mean Δ in every cohort and 153 have a negative one. Sign alone is a weak filter: the CLDN4-high arm is higher for most genes. The strict effect floor and the COMET decile are the lists to use.

## CLDN4 positive control

CLDN4 itself is excluded from the marker lists. In every used unit the Q4 mean
is above the Q1 mean, so the within-unit sign fraction is 1.

| cohort | units | patient Δ | cell Δ | AUC |
|---|---:|---:|---:|---:|
| GSE123902 | 12 | 1.882 | 1.799 | 0.984 |
| GSE131907 | 19 | 2.018 | 1.973 | 0.997 |
| GSE205335 | 21 | 2.439 | 2.323 | 0.999 |
| GSE189357 | 9 | 2.289 | 2.500 | 1.000 |

## 4/4 lists

- Strict patient-effect UP (Δ≥0.1, sign fraction≥0.6, AUC≥0.55, detection≥0.05 in all 4): **342** genes.
- Strict patient-effect DOWN: **0** genes.
- COMET top-decile UP (XL-mHG decile and patient Δ>0 in all 4): **196** genes.
- COMET top-decile DOWN: **12** genes.
- Same-sign patient Δ UP, no effect floor, ranked by the weakest cohort: **8569** genes.
- Same-sign patient Δ DOWN: **153** genes.

Strict UP (strongest minimum Δ first):

| gene | min patient Δ | max patient Δ | mt/ribo |
|---|---:|---:|---|
| ELF3 | 0.842 | 1.492 | False |
| JUN | 0.743 | 1.038 | False |
| SOX4 | 0.731 | 0.898 | False |
| TACSTD2 | 0.714 | 1.316 | False |
| FOS | 0.695 | 0.914 | False |
| NEAT1 | 0.576 | 1.120 | False |
| MUC1 | 0.566 | 1.319 | False |
| GPRC5A | 0.554 | 0.945 | False |
| CLDN3 | 0.535 | 0.811 | False |
| ZFP36L1 | 0.528 | 0.796 | False |
| KRT19 | 0.484 | 0.969 | False |
| RHOB | 0.476 | 0.923 | False |
| CD9 | 0.472 | 0.731 | False |
| EPCAM | 0.460 | 0.715 | False |
| SAT1 | 0.435 | 0.844 | False |
| … | 327 more in the TSV | | |


Strict DOWN:

_none_


Files: `results/gene_lists/`.

## What the lists are

The strict UP list is an epithelial / stress program that is higher in
CLDN4-high malignant cells in every cohort. The strongest minimum Δ values are
ELF3 (0.84 to 1.49), JUN, SOX4, TACSTD2, and FOS. Also on the strict list, and
on the COMET decile: MUC1, GPRC5A, CLDN3, KRT19, EPCAM, CLDN7, CDH1, F11R.
KRT8 and KRT18 are on it as well. 165 of the 342 strict UP genes are also in
the COMET top decile. One symbol on the strict list is mitochondrial/ribosomal.

This does not rank CLDN4 as a surface target. TACSTD2, MUC1, EPCAM, and CLDN7
move with it. NECTIN2 does not clear the strict rule.

No gene clears the strict DOWN floor (patient-mean Δ ≤ −0.10 in every cohort).
The 12 COMET down genes, and the top of the 153 same-sign down genes, are
ribosomal and elongation factors (RPL21, RPS12, EEF1A1, and the module
RPS12, EEF1A1, RPS15A, RPL13, RPS27A, RPL3, RPS3A, RPL26). That is the down
signal that actually agrees across cohorts. It is not an IFN or MHC-I marker
list. IFNG is same-sign down but far from the effect floor. STAT1, IRF1,
HLA-A, TAP1, IFI6, and IFIT1 are same-sign up. The locked IFN/MHC result was a
different contrast (patient pseudobulk of high-CLDN4 versus low-CLDN4
patients). This within-unit Q4 versus Q1 screen does not recover that program.

## Modules and enrichment

Co-expression modules use average linkage on the mean cross-cohort Spearman
(distance = 1 − r, cut at 0.5, modules smaller than 4 genes dropped). At that
cut almost every strict UP gene stays unclustered. The one UP module is the
immediate-early set JUN, FOS, EGR1, DUSP1, FOSB. The one DOWN module is the
ribosomal set above.

Enrichment was gseapy Enrichr (GO Biological Process 2023, KEGG 2021 Human,
MSigDB Hallmark 2020, Reactome 2022).

- Strict UP (n=342): 331 terms with adjusted P<0.05. Leading terms include
  Hallmark TNF-alpha signaling via NF-kB, p53, apoptosis, and unfolded protein
  response; KEGG protein processing in the endoplasmic reticulum, adherens
  junction, and tight junction.
- COMET UP (n=196): 242 terms with adjusted P<0.05.
- Same-sign DOWN top 100 (the strict DOWN list is empty, so this was the
  pre-set backup input): 124 terms, led by cytoplasmic translation and the
  KEGG ribosome. Chemokine and allograft terms in that table are carried by a
  minority of immune symbols in the top 100. Those symbols do not meet the
  strict Δ floor, so they are not called 4/4 down markers.
- COMET AND-pairs among the top 20 UP genes: 26 / 190 pass p<0.01, fold>1.5,
  and at least 10 CLDN4-high cells in the draw, in all 4 cohorts. TACSTD2
  pairs with CLDN3, CLDN7, and KRT19 are in that set. TM4SF1 pairs with most
  of the top of the list.

Marker-malignant gates in GSE123902 and GSE189357 require EPCAM or a KRT8/18/19
count above zero. A 4/4 gene still has to clear the two author-label cohorts
(GSE131907, GSE205335), so the gate alone cannot put a gene on the list.
