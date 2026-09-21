# Public omics for sacituzumab, sac-TMT / SKB264, and datopotamab

Query date: 2026-09-21. Archives: GEO, SRA, PRIDE, ArrayExpress.

## Lung

Open lung-cancer expression contrasts for these drugs: **0**.

`tables/lung_adc_expression_contrasts.tsv` has a header and no rows.

What the lung-token searches actually returned:

- GEO `(drug names) AND (lung OR NSCLC OR LUAD OR LUSC)` returned one Series, **GSE302284**. The deposited RNA is EGFR-TKI residual scRNA from a TROP2 CAR-T paper. It is not a sacituzumab, SKB264, or datopotamab arm.
- SRA `(drug names) AND (lung OR NSCLC)` returned one experiment inside **PRJNA754211**. That record is whole-exome DNA from the metastatic triple-negative breast cancer sacituzumab govitecan study (phs002555, runs in the dbGaP cluster). It has no expression matrix.

SKB264, sac-TMT, MK-2870, and tirumotecan: **0** records in GEO, SRA, PRIDE, and ArrayExpress.

Datopotamab, Dato-DXd, DS-1062, and Datroway: **0** GEO Series, **0** ArrayExpress experiments, **0** PRIDE projects. The only datopotamab SRA study is breast PDX **PRJNA1169861** (below).

No new GEO Series for these drug names was released between 2026-08-18 and 2026-09-21.

## Inventory

Full open/skip table: `tables/inventory.tsv`. Query counts: `tables/queries.tsv`.

Open expression matrices that really do treat with sacituzumab govitecan (IMMU132 / Trodelvy) are all colorectal or esophageal:

| Accession | Tissue | Contrast | Expression readout |
|---|---|---|---|
| GSE312098 | CRC CX-1 cell line | IMMU132 vs control, 48 h, 3 vs 3 | yes, scored below as a weak analogy |
| GSE311016 | CRC PDX | IMMU132 vs vehicle, day 29, 5 pairs | yes; scored earlier in PR 294 |
| GSE304294 | ESCC KYSE30 | IMMU132 vs vehicle, 1 day, 2 vs 3 | yes; scored earlier in PR 294 |
| E-MTAB-16433 | CRC PDOX | Trodelvy vs vehicle, 28 d, 4 vs 4 mice | yes; scored earlier in PR 294 |
| E-MTAB-16843 | CRC organoid | SG vs IgG1-SN-38 time course | yes, multi-GB matrix, not scored |
| E-MTAB-16849 | CRC liver-met PDOX | SG vs IgG1-SN-38 time course | yes, multi-GB matrix, not scored |

Other name hits are not drug-versus-control expression:

- **GSE278664** (PRJNA1168164), 35 ovarian biopsies. The Series design is BRCAmut versus BRCAwt on a prexasertib trial. There is no sacituzumab arm in the matrix.
- **PRJNA1169861**, datopotamab deruxtecan breast PDX study. Public data are 14 WXS, 8 genomic targeted-capture libraries, and 2 RNA-seq libraries (BCX.024, BCX.094). Those two RNA-seq biosamples have no treatment attribute and were collected in 2012 and 2014. They are baseline models, not a datopotamab-versus-control pair.
- **PRJNA754211 / phs002555**, clinical sacituzumab govitecan in metastatic TNBC. Public SRA experiments are WXS. **PRJNA754210** describes RNA and has 0 public SRA runs.
- **PRJNA1357832**, HCC1806 RNA-seq, 4 versus 4, is ASA (a TROP2-antibody BRD4 PROTAC) versus PBS. The abstract mentions sacituzumab govitecan and datopotamab as background.
- GEO text hits GSE309617/616 (carboplatin TNBC), the GSE303* KRAS/ERK pancreatic set, and GSE292860 (CDK7 inhibitor Q901) do not deposit an ADC arm.
- DS1062 GEO/SRA hits are yeast DNA polymerase delta. hRS7 SRA hits are 16S amplicons.
- PRIDE drug-name queries returned 0 projects. The TROP2 keyword returns **PXD065965** (prostate TROP2 extracellular-domain shedding), which is not an ADC treatment.

## Earlier PRs

| PR | What it already did |
|---|---|
| 294 | Scored GSE311016, GSE304294, and E-MTAB-16433. Left GSE312098 as a given colon result. |
| 341 | Second hunt. No new ADC-versus-control RNA series. |
| 393 | Third hunt. Leftover scored table empty. |
| 244 | GSE302284 is osimertinib residual RNA, not an ADC. |
| 233 | CD274 on the GSE312098 IMMU132 arm. |

This note does not replace those PRs. It re-queries the four archives on 2026-09-21 and computes GSE312098, which those hunts did not re-score.

## Weak analogy: GSE312098 colon CX-1

This is a colorectal cell line, 48 h, 3 µg/ml IMMU132 versus control, n=3 versus 3. It is the closest open ADC-versus-control matrix with TACSTD2, CLDN4, and interferon readouts. It is not a lung result.

Author FPKM, UTF-16 supplementary file. One protein-coding row per gene symbol (highest mean FPKM). Statistic: Welch t-test on log2(FPKM+1). Panel q is Benjamini-Hochberg inside the reported gene panel. IFNG is FPKM 0 in every sample, so its p-value is NA and it is left out of that BH correction.

IMMU132 versus control:

| Gene | Mean FPKM control | Mean FPKM IMMU132 | log2FC | Welch p | Panel q |
|---|---:|---:|---:|---:|---:|
| TACSTD2 | 41.5 | 67.2 | +0.683 | 5.65×10⁻⁴ | 0.0030 |
| CLDN4 | 91.8 | 50.1 | −0.861 | 1.69×10⁻⁵ | 5.42×10⁻⁴ |
| IFNG | 0 | 0 | 0 | NA | NA |
| STAT1 | 30.0 | 33.7 | +0.164 | 0.14 | 0.17 |
| ISG15 | 115 | 179 | +0.634 | 0.0033 | 0.012 |
| IFI6 | 97.3 | 187 | +0.948 | 0.0098 | 0.026 |
| OAS2 | 4.28 | 8.98 | +0.942 | 0.025 | 0.043 |
| IFIT1 | 5.70 | 8.42 | +0.492 | 0.0065 | 0.019 |
| CD274 | 0.23 | 0.54 | +0.324 | 0.019 | 0.042 |

Hallmark sets (MSigDB 2024.1.Hs), same contrast, tested against the other 19,740 protein-coding symbols:

| Set | Genes in matrix | Median log2FC | Fraction > 0 | Mann-Whitney p |
|---|---:|---:|---:|---:|
| Interferon gamma response | 196 | +0.131 | 0.755 | 2.0×10⁻¹⁷ |
| Interferon alpha response | 96 | +0.321 | 0.812 | 7.7×10⁻¹⁷ |

The Mann-Whitney p-values are small because the background is the rest of the protein-coding transcriptome. The interferon-gamma shift itself is a median log2FC of +0.13. IFNG RNA is absent. CXCL9 and CXCL10 stay below 0.4 FPKM. CD274 rises on the log scale and remains below 1 FPKM.

CLDN4 falls. CLDN3 (−0.70) and TJP1 (−0.48) fall with it. CLDN1 rises (+0.66). The matrix does not show one direction for every tight-junction gene.

The same file has an IMMU132 + GSK2606414 (PERK inhibitor) arm, also 3 versus 3. That arm is combination treatment: TACSTD2 log2FC +0.799, CLDN4 −0.614, interferon-alpha median log2FC +0.456. It is reported in `tables/gse312098_key_genes.tsv` and is not an ADC-only contrast.

Figure: `figures/gse312098_immu132_keygenes.png`.

Reproduce with `python3 methods/trop2_adc_public_omics_20260921/analyze_gse312098.py`. The script downloads the GEO FPKM file and the MSigDB hallmark GMT into `data/`, which is gitignored.
