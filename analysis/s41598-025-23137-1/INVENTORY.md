# NHEJ–STING inventory: Sci Rep 2025 CLDN4 CRISPRi

Villagomez et al., Scientific Reports, 10 Nov 2025. [doi:10.1038/s41598-025-23137-1](https://doi.org/10.1038/s41598-025-23137-1). PMID 41214101, PMC12603150.

## What is deposited

The CRISPRi experiment that reports higher pSTING after CLDN4 loss has no RNA-seq, proteome, GEO series, SRA run, BioProject, GSA accession, or ProteomeXchange ID. Europe PMC lists supplementary files and no database cross-references. OmicsDI returns 0 datasets for the DOI. The data-availability statement says every raw file is available on request from the corresponding authors (Benjamin G. Bitler, benjamin.bitler@cuanschutz.edu; Fabian R. Villagomez, fromerovillag@utep.edu).

The pSTING result is confocal microscopy and immunoblot, not a count matrix. There is nothing to download that can score IFN, STING, or NHEJ transcripts inside those CRISPRi cells.

| Place searched | Query | Result |
|---|---|---|
| GEO / SRA / BioProject | elink from PMID 41214101 | no linked records |
| OmicsDI | `s41598-025-23137-1` | 0 datasets |
| GEO | `CLDN4` and (silencing or knockdown or CRISPRi) | only [GSE22493](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22493) (Gao, SKOV-3, 2010). Not this paper |
| GEO | `Bitler[Author]` and claudin | 0 series |
| GEO | `Villagomez[Author]` | 11 series; none are a CLDN4 CRISPRi. The ovarian series are EHMT inhibitor profiles from the Bitler lab |
| PRIDE | `claudin-4 OVCAR3`; `CLDN4 BioID` | empty |
| ProteomeXchange PXD032041 | opened the XML | prostate enzalutamide cistrome, unrelated |
| GSA / NGDC | DOI and CLDN4 search pages; OmicsDI | no CRA or HRA accession for this article |
| Figshare API | claudin-4 Villagomez / CRISPRi | no matching record |
| Zenodo API | same | HTTP 403 from this environment |
| Springer supplements | MOESM1–5 | downloaded. Tables are flow antibodies, a cBioPortal ISG export, and an RPPA export. MOESM4 is figure legends plus blot images, with no gene matrix |

## What the CRISPRi cells actually show

Perturbation, from the methods: dCas9-KRAB-MeCP2 (Addgene 110824) in OVCAR3 and OVCA429, guide `GCTGGCTTGCGCATCAGGAC` against human CLDN4. OVCAR8, which has no endogenous CLDN4, was used for overexpression.

Published directions, with no numeric fold-change in the text:

- CLDN4 overexpression in OVCAR8 lowered the pSTING confocal signal (Fig. 2a).
- CLDN4 knockdown in OVCAR3 and OVCA429 raised the pSTING confocal signal (Fig. 2b, c). The paper describes this as a noticeable increase. Total STING protein also changed on immunoblot (Supplementary Fig. 4d–f).
- The ISRE luciferase reporter moved the other way. Overexpression raised the basal type I interferon response in OVCAR8 (Fig. 3a). Knockdown produced a dramatic reduction of that response in OVCAR3 (Fig. 3c). cGAMP did not bring the knockdown cells back to the wild-type ISRE level (Fig. 3e).

The knockdown phenotype in this paper is higher pSTING staining with a lower ISRE reporter. The reporter is the functional interferon readout. It is not an RNA-seq gene score. Assay-by-assay wording, including the Fig. 2a versus Fig. 2e overexpression sentences and the pTBK1 blots, is in `RESULTS.md`. No CRISPRi count matrix was added there.

## Scores that are computable

### 1. The paper’s own ISG table

Supplementary Table 2 is a cBioPortal co-occurrence export on TCGA ovarian serous cystadenocarcinoma. The caption says 617 samples. Queried ISGs in the text: ISG15, MX1, OAS1, IFIT1, EIF2AK2, IRF7, STAT1, RSAD2, BST2, IFI44. The table prints five rows. The text names IRF7 and OAS1 as the high-CLDN4 associations.

| Gene 1 | Gene 2 | log2 odds ratio | p | q | Tendency |
|---|---|---:|---:|---:|---|
| IRF7 high | CLDN4 high | 1.991 | 0.004 | 0.007 | co-occurrence |
| CLDN4 high | CLDN4 low | < −3 | 0.010 | 0.016 | mutual exclusivity (same gene) |
| BST2 high | CLDN4 low | < −3 | 0.029 | 0.042 | mutual exclusivity |
| OAS1 high | CLDN4 high | 1.194 | 0.035 | 0.049 | co-occurrence |
| IFI44 high | CLDN4 low | 1.571 | 0.036 | 0.050 | co-occurrence |

Source file: `tables/paper_supp_table2_isg_cbioportal.tsv`.

Supplementary Table 3 is RPPA protein, stratified by CLDN4 mRNA z-score greater than 1 versus less than −1, Student t-test with Benjamini–Hochberg q. Eighteen proteins pass q < 0.05. The significant set contains no STING1, CGAS, TBK1, IRF3, PRKDC, XRCC5, XRCC6, LIG4, or NHEJ1. CHEK1 protein is higher when CLDN4 mRNA is low (log2 ratio −0.16, q = 0.0389). STAT3 pY705 is higher when CLDN4 mRNA is high (log2 ratio 0.22, q = 0.0138). The cBioPortal RPPA profile returns no row for CLDN4 itself. Full table: `tables/paper_supp_table3_rppa.tsv`. The first numeric pair differs by the printed log2 ratio, so those two columns are the group means. The docx repeats the group header on the next pair and does not name it.

### 2. TCGA PanCancer Atlas RNA-seq, continuous Spearman

Study `ov_tcga_pan_can_atlas_2018`, profile `rna_seq_v2_mrna` (RSEM). Pulled from the cBioPortal API on 21 Sep 2026. n = 300 tumors with RNA-seq values. Spearman correlation of log2(RSEM+1) against CLDN4. Benjamini–Hochberg q is across the 59 panel genes that returned values. This is the public cohort the paper used. It is a different summary from their 2020 high/low odds ratios, and the current RNA-seq profile has 300 tumors rather than the 617 samples in their table caption.

Panel mean (mean of per-gene z-scores, then Spearman versus CLDN4):

| Panel | Genes in the mean | rho | p | q |
|---|---:|---:|---:|---:|
| IFN / ISG | 30 | +0.037 | 0.53 | 0.53 |
| STING axis | 12 | +0.062 | 0.28 | 0.42 |
| NHEJ | 17 | −0.202 | 4.3×10⁻⁴ | 0.0013 |

Genes with q < 0.05:

| Panel | Gene | rho | p | q |
|---|---|---:|---:|---:|
| NHEJ | PRKDC | −0.258 | 5.9×10⁻⁶ | 3.5×10⁻⁴ |
| STING | TREX1 | +0.229 | 6.3×10⁻⁵ | 0.0019 |
| NHEJ | SHLD2 | −0.211 | 2.3×10⁻⁴ | 0.0045 |
| STING | STING1 | +0.192 | 8.5×10⁻⁴ | 0.012 |
| IFN | BST2 | +0.189 | 9.8×10⁻⁴ | 0.012 |
| STING | RAB7A | +0.169 | 0.0033 | 0.032 |
| NHEJ | TP53BP1 | −0.166 | 0.0040 | 0.032 |
| NHEJ | RIF1 | −0.164 | 0.0043 | 0.032 |

The paper’s two named ISGs on this continuous test: IRF7 rho = +0.081, q = 0.43; OAS1 rho = +0.068, q = 0.51. BST2, which their odds-ratio table places with CLDN4-high, is the ISG that stays positive here.

Core NHEJ transcripts that do not pass q < 0.05: XRCC5 (Ku80) rho = −0.118, q = 0.19; XRCC6 (Ku70) rho = −0.065, q = 0.51; LIG4 rho = +0.071, q = 0.50; NHEJ1 rho = −0.082, q = 0.43. CGAS rho = −0.093, q = 0.34. TBK1 rho = −0.133, q = 0.11. IRF3 rho = +0.128, q = 0.13.

Full table: `tables/tcga_ov_rnaseq_cldn4_spearman.tsv`.

### 3. CPTAC protein in the same cBioPortal study

Profile `protein_quantification`. CLDN4 protein is quantified in 104 tumors. Pairwise n is 82–104 depending on the gene. Spearman versus CLDN4 protein. Zero of 44 tested proteins have q < 0.05.

| Protein | n | rho | p | q |
|---|---:|---:|---:|---:|
| STING1 | 99 | −0.118 | 0.24 | 0.47 |
| CGAS | 86 | −0.088 | 0.42 | 0.64 |
| TBK1 | 104 | −0.268 | 0.0060 | 0.15 |
| PRKDC | 104 | +0.127 | 0.20 | 0.42 |
| TP53BP1 | 104 | −0.251 | 0.010 | 0.15 |
| XRCC5 | 104 | +0.129 | 0.19 | 0.42 |
| XRCC6 | 104 | +0.123 | 0.21 | 0.43 |
| RAB7A | 104 | −0.010 | 0.92 | 0.94 |
| STAT1 | 104 | +0.227 | 0.021 | 0.23 |
| RIGI (DDX58) | 104 | +0.214 | 0.029 | 0.25 |

Panel means: IFN rho = +0.164, q = 0.14; STING rho = −0.248, p = 0.011, q = 0.034; NHEJ rho = −0.003, q = 0.98. IRF7 protein has no values in this profile.

STING1 mRNA rises with CLDN4 mRNA. STING1 protein does not. PRKDC mRNA falls with CLDN4 mRNA. PRKDC protein does not. Full table: `tables/tcga_ov_cptac_protein_cldn4_spearman.tsv`.

### 4. Claudin-4 BioID protein list the 2025 paper cites

Neville et al., Heliyon 2022, [doi:10.1016/j.heliyon.2022.e10862](https://doi.org/10.1016/j.heliyon.2022.e10862), PMC9552118. This is the OVCAR3 BioID that the 2025 paper uses for the Rab7 interaction. It is an open PDF supplement (mmc1), not a PRIDE deposit. Supplementary Table 1 lists 54 enriched proteins.

RAB7A is in that list: control intensity 3.00, CLDN4-BioID intensity 40.25, average ratio 22.04, log2 ratio 9.03, 4 experiments.

The enriched list does not contain CGAS, STING1, TMEM173, TBK1, IRF3, TREX1, BECN1, IFNB1, or the core NHEJ proteins scored above. Parsed table: `tables/heliyon2022_bioid_supp_table1.tsv`. Presence calls: `tables/heliyon2022_bioid_pathway_presence.tsv`.

### 5. GSE22493, the only open CLDN4-silencing transcriptome

Gao, GEO series GSE22493, public 15 Dec 2010. SKOV-3-IP-Luc, CLDN4 siRNA lentivirus versus a CLDN4-overexpression channel, two-color, three arrays. Channel 1 / Cy3 is the overexpression control. Channel 2 / Cy5 is the knockdown. The series matrix stores a signed sample-to-control log-ratio, so a negative value means lower in the knockdown channel.

This is not the 2025 OVCAR3 / OVCA429 CRISPRi. The Operon v3 platform has no probe for STING1, TMEM173, CGAS, or TBK1. Replicate 1 is blank for 5,771 of 36,284 probes, including the single CLDN4 probe (ID 17169: missing, −1.737, −0.713; mean of the two present arrays −1.23). For 33 of 41 scored panel genes the replicate SD exceeds the absolute mean. No panel gene has q < 0.05. The array confirms that the CLDN4 channel is lower in the knockdown on the two arrays that have a spot. It does not score an IFN or STING transcriptional program. Table: `tables/gse22493_kd_over_oe_logratio.tsv`.

## Reading for the mechanism

The 2025 knockdown cannot be scored at gene level until the authors release the raw blots or a count matrix. On the assays they did publish, loss of CLDN4 raises the pSTING image and lowers the ISRE reporter.

In HGSOC tumors, CLDN4 mRNA does not track a broad IFN-gene mean (rho +0.037). It tracks higher STING1, TREX1, BST2, and RAB7A mRNA, and a lower NHEJ-gene mean driven by PRKDC, SHLD2, TP53BP1, and RIF1. The CPTAC protein matrix in the same study does not carry those mRNA correlations forward: the STING-protein mean is negative, and PRKDC protein is not anti-correlated with CLDN4 protein. The BioID neighborhood of CLDN4 contains RAB7A and does not contain the cGAS–STING or NHEJ enzymes.

## Reproduce

```bash
python3 analysis/s41598-025-23137-1/scripts/score_ifn_sting_nhej.py \
  --gse-matrix /path/GSE22493_series_matrix.txt.gz \
  --gse-platform /path/GPL10555_family.soft.gz
```

Requires Python 3 and SciPy. cBioPortal responses are cached under `--cache` (default `/tmp/nhej-sting/cbio_cache`). GEO files are the public series matrix and platform soft for GSE22493 / GPL10555.
