# ImmPort / TCIA — open NSCLC immunotherapy transcriptomes for CLDN4 vs NHEJ, IFN, APM

**No open human NSCLC immunotherapy expression matrix is in ImmPort or The Cancer Imaging Archive.** CLDN4 vs NHEJ / IFN / APM was therefore not computed on an ICI cohort. Response was not tested.

The only open TCIA-linked NSCLC transcriptome that contains CLDN4 is surgical RNA-seq, not immunotherapy. On that matrix CLDN4 tracks keratins. The raw NHEJ and APM associations do not remain after a keratin residual. The IFN mean is flat. This is bulk co-expression. It is not a spatial-exclusion result and it does not revise the locked CosMx exclusion page.

Reproduce: `python3 methods/immport_tcia_cldn4_nhej/analyze.py`

---

## Immunotherapy gate

| Source | What was opened | Expression? | ICI? | Scored? |
|---|---|---|---|---|
| ImmPort shared-study search | 86 studies from Oncology plus NSCLC / PD-1 / checkpoint terms | File paths return HTTP 401 without a token | No human NSCLC tumor RNA-seq on ICI in the public titles | No |
| TCIA `Anti-PD-1_Lung` | 46 lung cases, anti-PD-1, 2016 | PT / CT / SC only. The only spreadsheet on the page is an NBIA image manifest, not expression | Yes, imaging | No |
| TCIA `S0819` | n=1,299 imaging | No RNA file | Cetuximab ± chemotherapy, not PD-1/PD-L1 | No |
| TCIA `CMB-LCA` | Moonshot Biobank lung imaging | No RNA matrix on the collection page | Not an ICI transcriptome | No |
| TCIA CPTAC-LUAD / CPTAC-LSCC, TCGA-LUAD / TCGA-LUSC | Imaging collections; genomics are external | Resection / TCGA, not ICI | No | No |
| TCIA NSCLC-Radiogenomics | [GSE103584](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103584) RNA-seq, 130 tumors | Yes | **No.** Surgically excised tumors, 2008–2014. Clinical file has 0 ICI drug tokens | Yes, labeled **NOT_ICI** |
| TCIA NSCLC-Radiomics-Genomics | [GSE58661](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58661) Lung3 microarray, 89 biopsies | Probe matrix yes; **CLDN4 symbol not in the matrix** | Surgery, not ICI | No |
| Other lung TCIA collections (Radiomics, PET, 4D, NLST, LIDC, RIDER, QIN) | Imaging | No | No | No |

ImmPort lung or checkpoint titles are mouse models (SDY901, SDY1064, SDY1079, SDY1298), CyTOF or flow (SDY3403, SDY2310), COVID or SARS-CoV-2 studies, or non-lung checkpoint studies (HCC, MC38, diabetes). SDY901 is the mouse EGFR PD-1-pathway paper, not patients on ICI. OAK / POPLAR / IMpower150 remain controlled-access and were not used.

`GSE28827` is the 26-subject microarray from the same radiogenomics surgery study. It was not scored again.

---

## One row that was actually scored

Cohort: **GSE103584**, TCIA NSCLC-Radiogenomics, Bakr et al. *Sci Data* 2018. Scale: `log2(x+1)` of the deposited linear matrix. NA stayed missing (not filled as 0). A gene enters a score only if it is non-missing in ≥80% of 130 samples. Score = mean of gene-wise z. Patient = Case ID. All 130 RNA columns match the clinical file (211 clinical rows total).

Histology of the 130: adenocarcinoma 96, squamous 31, NSCLC NOS 3. Days between CT and surgery are recorded for 130/130. Adjuvant treatment Yes 37, chemotherapy Yes 37, radiation Yes 14. ICI drug strings in the clinical file: 0. CLDN4 finite in **128/130** (median log2 4.49). Those 128 are the correlation n.

Primary family = 3 program Spearman tests. BH q is across those 3 only. Keratin residual is sensitivity (Pearson of rank residuals on the KRT8/KRT18/KRT19 mean-z; df = n − 3).

| Program | Genes in score | ρ (p) | BH q | ρ \| KRT (p) |
|---|---|---|---|---|
| NHEJ | 7/8 (PAXX absent): XRCC6, XRCC5, PRKDC, LIG4, XRCC4, NHEJ1, DCLRE1C | **+0.327 (1.61×10⁻⁴)** | **4.84×10⁻⁴** | +0.053 (0.553) |
| IFN | 9/12: STAT1, IRF1, CXCL9, CXCL10, IDO1, GBP4, GBP5, CXCR6, IL2RG | **−0.002 (0.982)** | 0.982 | −0.077 (0.387) |
| APM | 7/12: HLA-A, HLA-B, HLA-C, B2M, PSMB10, NLRC5, ERAP1 | **+0.279 (1.41×10⁻³)** | **2.11×10⁻³** | +0.161 (0.071) |

Numeric rows: `tables/one_row.tsv`, `tables/spearman.tsv`.

### How to read the row

**Not immunotherapy.** Do not use these ρ values as ICI response, resistance, or “CLDN4 rises after ICI.”

**Keratins.** CLDN4 vs KRT8 / KRT18 / KRT19 is ρ = +0.629 / +0.490 / +0.531 (n = 128 / 128 / 121). The epithelial axis is the main CLDN4 companion in this bulk matrix.

**NHEJ.** The raw positive score is real and is carried by XRCC6 (ρ=+0.509, n=128), XRCC5 (+0.340), and NHEJ1 (+0.297). After the keratin residual the program association is gone (+0.053, p=0.553). This is not a keratin-independent NHEJ program.

**IFN.** The pre-specified mean is flat. It is not an IFNG test: IFNG is finite in 51/130 and was left out of the score. Among those 51, IFNG vs CLDN4 is ρ=−0.023 (p=0.88). Inside the score, STAT1 is positive (ρ=+0.319, n=128) and GBP5 is negative (ρ=−0.265, n=122). CXCL9 is −0.096 (n=125) and CXCL10 is −0.008 (n=120). The mean cancels. Gene-level BH q values are in `tables/gene_spearman.tsv` and are not the primary family.

**APM.** TAP1 is finite in 9/130, TAP2 in 48, TAPBP in 42, PSMB8 in 42, PSMB9 in 48. They were excluded by the 80% rule. The score that was tested is HLA-A/B/C, B2M, PSMB10, NLRC5, and ERAP1, not a full TAP-containing APM. The raw positive ρ shrinks to +0.161 (p=0.071) after keratin.

**Not spatial exclusion.** Same-sample bulk correlation is not a neighbor-count or Ripley result. The locked CosMx exclusion statement is untouched.

---

## GSE58661

Lung3 (Aerts et al.), TCIA NSCLC-Radiomics-Genomics, n=89, surgery, not ICI. The series matrix has 60,607 probes named `AFFX-*` or `merck-<accession>_at`. GPL10379 (the symbol table for this HuRSTA array) lists 2 exact CLDN4 probes (`100123349_TGI_at`, `100303030_TGI_at`). Neither ID is in the matrix, and CLDN4’s RefSeq accession `NM_001305` is not a `merck-` probe name. Programs were not scored. No second claudin was substituted.

---

## What was not done

- No ICI response, DCB, or “CLDN4 up after resistance” test.
- No merge with GSE126044, GSE135222, GSE248249, GSE248378, or any other ICI bulk series.
- No Fisher combination of the surgical ρ with an ICI cohort.
- No CPTAC or TCGA re-analysis.
- No mouse ImmPort matrix (SDY901 and the other oncology models).
- No TACSTD2 column on the claim table.
- NA was not imputed as zero. Multi-gene microarray probes were not used.

---

## Files

- `FINDING.md` — this page
- `analyze.py` — live ImmPort and TCIA queries, download-once matrices, scores, tables, figures
- `tables/one_row.tsv`, `spearman.tsv`, `gene_spearman.tsv`, `gene_coverage.tsv`
- `tables/immport_studies.tsv`, `tcia_collections.tsv`, `summary.json`
- `figures/gse103584_cldn4_vs_programs.png`, `figures/forest_not_ici.png`
