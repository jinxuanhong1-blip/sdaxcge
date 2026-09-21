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

A second pass the same day, aimed at lung PDX and cell-line ADC RNA, is still empty. GEO queries for PDX/xenograft plus these drugs plus lung/NSCLC/LUAD, for Calu-3/HCC827/H1975 plus the drugs, and for TROP2/TACSTD2 plus antibody-drug conjugate plus lung, each returned either 0 Series or only **GSE302284**. SRA returned 0 for PDX plus drug plus lung, 0 for Calu-3 plus the drugs, and 0 for TROP2 plus the drugs plus lung restricted to RNA-Seq. BioProject drug-name AND lung/NSCLC returned **PRJNA1289808**, which is GSE302284. **GSE235812** is breast tumor versus matching PDX TACSTD2 correlation, with no ADC arm. `tables/lung_adc_expression_contrasts.tsv` stays header only.

## Inventory

Full open/skip table: `tables/inventory.tsv`. Query counts: `tables/queries.tsv`.

Open expression matrices that really do treat with sacituzumab govitecan (IMMU132 / Trodelvy) are all colorectal or esophageal:

| Accession | Tissue | Contrast | Expression readout |
|---|---|---|---|
| GSE312098 | CRC CX-1 cell line | IMMU132 vs control, 48 h, 3 vs 3 | yes, scored below as a weak analogy |
| GSE311016 | CRC PDX | IMMU132 vs vehicle, day 29, 5 pairs | yes; in the sweep below |
| GSE304294 | ESCC KYSE30 | IMMU132 vs vehicle, 1 day, 2 vs 3 | yes; in the sweep below |
| E-MTAB-16433 | CRC PDOX | Trodelvy vs vehicle, 28 d, 4 vs 4 mice | yes; sweep uses the PR 294 mouse pseudobulk |
| E-MTAB-16843 | CRC organoid | SG vs IgG1-SN-38 time course | not scored this round: ftp.ebi.ac.uk TLS failed |
| E-MTAB-16849 | CRC liver-met PDOX | SG vs IgG1-SN-38 time course | not scored this round: ftp.ebi.ac.uk TLS failed |

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

Reproduce the single-series table with `python3 methods/trop2_adc_public_omics_20260921/analyze_gse312098.py`. The script downloads the GEO FPKM file and the MSigDB hallmark GMT into `data/`, which is gitignored.

## Method sweep: CLDN4 down, IFN, APM, NHEJ

Grid: 78 rows. Contrasts are the ADC-versus-control arms above, plus the two combination arms (IMMU132+GSK2606414; IMMU132+IACS010759). Transforms are log2(x+1) and log2(x+0.1); for the PDOX counts also log2(CPM+1), log2(CPM+0.1), and log2(count+1). Set summaries are the median and the mean. Three set bundles are defined in `sweep_cldn4_ifn_apm_nhej.py`.

Pre-specified primary row, fixed before ranking: log2(x+1) or log2(CPM+1), median log2FC, Hallmark interferon-gamma, Reactome class I MHC antigen presentation (29 genes), Reactome NHEJ (34 of 68 symbols present; the absent symbols are histone genes under current names such as H2BC and H4C).

signed_score = (−CLDN4 log2FC) + IFN summary + APM summary + NHEJ summary.

A sign hit requires CLDN4 log2FC < 0 and all three summaries > 0. Moderate threshold: CLDN4 < −0.25 and the same sign rule. Strict threshold: CLDN4 < −0.5 and each summary > 0.25.

Primary ADC-monotherapy rows:

| Contrast | n | CLDN4 log2FC (Welch p) | IFN-γ median | APM median | NHEJ median | Sign hit |
|---|---:|---:|---:|---:|---:|---|
| GSE312098 IMMU132 vs control | 3 vs 3 | −0.861 (1.69×10⁻⁵) | +0.131 | +0.089 | +0.041 | yes |
| GSE311016 IMMU132 vs control | 5 pairs | −0.437 (Welch 0.15; paired 0.060) | +0.001 | −0.307 | −0.119 | no |
| GSE304294 IMMU132 vs control | 2 vs 3 | +0.911 (3.90×10⁻⁵) | +0.035 | +0.109 | −0.145 | no |
| E-MTAB-16433 SG vs vehicle | 4 vs 4 mice | −0.345 (0.063) | −0.034 | −0.145 | +0.057 | no |

On that primary row, GSE312098 is the only sign hit. Its NHEJ median is +0.041, and the sample-level Welch test on the mean NHEJ score is p = 0.91. The 8-gene classical NHEJ core (XRCC5, XRCC6, PRKDC, LIG4, XRCC4, NHEJ1, DCLRE1C, PAXX) has median log2FC **−0.233** on the same transform, so the NHEJ leg changes sign with the gene-set definition.

Strict threshold: **0 / 78** rows. Moderate threshold: **6 / 78**, all of them GSE312098 IMMU132 versus control.

The highest ADC-only sign-hit score in the grid is still GSE312098, using log2(x+0.1), the primary sets, and the median: CLDN4 −0.873, IFN +0.328, APM +0.159, NHEJ +0.042, score 1.40. The larger interferon median is the pseudocount. NHEJ stays near +0.04.

E-MTAB-16433 produces 6 sign hits only on log2(count+1), where CLDN4 log2FC is −0.149. Those rows miss the moderate threshold. The library-normalized CPM rows do not sign-hit. Mouse total counts in the pseudobulk are 11.3–21.3 million (SG) and 11.0–15.6 million (vehicle). The cell-level metadata was not re-opened: `ftp.ebi.ac.uk` and `ftp.sra.ebi.ac.uk` closed the TLS handshake from this environment.

Combination arms on the primary sets and median log2(x+1) are not sign hits. GSE312098 combination: CLDN4 −0.614, IFN +0.241, APM +0.154, NHEJ −0.070. GSE304294 combination: CLDN4 +1.36, NHEJ −0.180.

Datopotamab still has no public treatment-versus-control RNA matrix, so it is not in the sweep.

Figure: `figures/sweep_primary_adc_only.png`. Full grid: `tables/sweep_all_rows.tsv`.

Reproduce with `python3 methods/trop2_adc_public_omics_20260921/sweep_cldn4_ifn_apm_nhej.py`.
