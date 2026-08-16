# B5_SKCM — Pre-specified analysis plan

**Question.** Is tumor `CLDN4` mRNA expression in pre-treatment biopsies associated with
response to immune checkpoint inhibition (ICI) in metastatic cutaneous melanoma (SKCM)?

**Status.** This plan is written and committed *before* any test that involves the response
outcome is run. Only outcome-blind steps (file format inspection, column/sample-name mapping,
gene annotation) preceded it. The commit history of this repository is the audit trail.

**Design.** Retrospective secondary analysis of publicly downloadable bulk RNA-seq from ICI-treated
melanoma cohorts. This is an observational, hypothesis-testing exercise, not a biomarker
qualification study.

---

## 1. Cohorts

Inclusion rule: a public dataset is eligible if (a) it profiles bulk RNA from human metastatic
melanoma tumors, (b) patients received an approved ICI (anti-PD-1, anti-CTLA-4, or combination),
(c) a per-sample response label is distributed with the data, and (d) a processed expression matrix
is downloadable without controlled access (no dbGaP/EGA application).

Eligible and included:

| Cohort | Accession | Therapy | Matrix used |
|---|---|---|---|
| Hugo 2016 (UCLA/VIC) | GSE78220 | pembrolizumab | author FPKM (gene symbols) |
| Riaz 2017 (BMS-038) | GSE91061 | nivolumab | author FPKM (Entrez IDs, hg19 knownGene) |
| Auslander 2018 (MGH) | GSE115821 | anti-PD-1 / anti-CTLA-4 | author raw counts + gene length -> TPM |

Screened and excluded, with reason (recorded so the cohort set is not silently curated):

- Gide 2019 (PRJEB23709), Liu 2019, Van Allen 2015, Nathanson 2017: no processed expression
  matrix in an open repository; raw data behind controlled access. Excluded.
- GSE282471, GSE300793: melanoma ICI series with no per-sample response label in GEO metadata
  (GSE282471 also ships only raw FASTQ/tar). Excluded.
- GSE93157 and other NanoString/panel series: targeted panels that do not assay `CLDN4`. Excluded.
- Mouse/cell-line/mechanistic series returned by the GEO query: not human patient cohorts. Excluded.

Non-ICI reference cohort (context only, no response analysis): TCGA-SKCM (UCSC Xena GDC hub
`TCGA-SKCM.star_tpm.tsv.gz`), used solely to characterize the `CLDN4` expression distribution in
melanoma and its correlates.

## 2. Sample selection

- **Pre-treatment (baseline) biopsies only.** On-treatment biopsies are excluded from the primary
  analysis (Hugo `Pt16.OnTx`; Riaz `visit == On`; MGH `treatment state == On`).
- **One observation per patient.** Where a patient has replicate baseline biopsies (Hugo
  `Pt27A`/`Pt27B`), the log2 expression values are averaged. Where a patient has baseline biopsies
  from more than one distinct treatment episode (possible in MGH, which contains sequential
  anti-CTLA-4 then anti-PD-1 lines), the earliest episode is used for the primary analysis.
- No sample is excluded on the basis of its expression values or its response label.
- Library-size / QC exclusion: none pre-specified. Total counts and detected-gene counts are
  reported for transparency.

## 3. Exposure variable

- `CLDN4` (Entrez 1364; HGNC symbol CLDN4). Resolved via NCBI `Homo_sapiens.gene_info`, not by
  hard-coded row indices.
- Transformed as `log2(x + 0.1)` where `x` is FPKM (Hugo, Riaz) or TPM (MGH, TCGA).
- For cross-cohort pooling, standardized within cohort to mean 0 / SD 1 (`z(CLDN4)`).
- **Detectability is reported before interpretation**: fraction of baseline samples with
  `CLDN4 > 0`, and the cohort median. If `CLDN4` is undetectable in a large fraction of samples,
  the association test is still run and reported, but the result is flagged as
  low-information rather than quietly dropped or quietly interpreted.

## 4. Outcome variable

Binary clinical response at the primary analysis:

- Hugo: `Complete Response` / `Partial Response` -> R; `Progressive Disease` -> NR.
- Riaz: `PRCR` -> R; `PD` -> NR; **`SD` excluded from the primary analysis**; `UNK` always excluded.
- MGH: `R` -> R; `NR` -> NR (as labeled by the authors).

Pre-specified sensitivity definition S1: Riaz `SD` recoded as NR (i.e. R = PRCR, NR = PD+SD).

## 5. Statistical analysis

**Primary endpoint.** Two-sided Mann–Whitney U test of `log2 CLDN4` in R vs NR, computed
independently in each of the three cohorts, and a fixed-effect inverse-variance meta-analysis of
the per-cohort logistic-regression log-odds ratio per 1 SD of `log2 CLDN4`.

- Significance threshold: alpha = 0.05 two-sided for the meta-analytic estimate. The three
  per-cohort tests are treated as descriptive components of the meta-analysis, not as three
  independent chances to declare a hit.
- Effect sizes reported with 95% CIs: odds ratio per SD, rank-biserial correlation, AUC (=
  probability that a random responder has higher `CLDN4` than a random non-responder), and
  difference in medians.
- Heterogeneity: Cochran's Q, I-squared; DerSimonian–Laird random-effects estimate reported
  alongside the fixed-effect estimate.

**Genome-wide contextualization (the main guard against a false positive).** The identical
per-cohort test and meta-analysis is run for every gene measured in all three cohorts. `CLDN4` is
then reported as a rank and percentile within that empirical null, together with a permutation-free
empirical p-value (fraction of genes with a smaller meta-analytic p-value) and the
Benjamini–Hochberg FDR of `CLDN4` in the genome-wide sweep. A gene that is nominally significant
but sits in the bulk of the genome-wide distribution is reported as such.

**Positive controls (assay sensitivity).** The same pipeline is run for genes with a
well-replicated association with ICI benefit in melanoma: `CD8A`, `CD274`, `IFNG`, `GZMB`, `PRF1`,
`CXCL9`, `CXCL10`, `STAT1`, `HLA-DRA`, plus the 6-gene Ayers IFN-gamma score
(`IFNG`, `STAT1`, `CXCL9`, `CXCL10`, `IDO1`, `HLA-DRA`). If these controls are not detectable here,
the analysis is underpowered and a null result for `CLDN4` is uninformative rather than negative.

**Negative / specificity controls.** `CLDN3`, `CLDN7`, `EPCAM`, `CDH1`, `KRT5`, `KRT14`, `KRT1`,
`KRT10`, `SFN` — epithelial and keratinocyte genes. Because `CLDN4` is an epithelial tight-junction
gene and melanoma is neural-crest derived, a `CLDN4` signal in a skin/subcutaneous biopsy may
reflect contaminating normal epidermis or adnexal structures rather than tumor biology. Reported:
Spearman correlation of `CLDN4` with each of these genes and with a keratinocyte score, and whether
the `CLDN4`-response association survives adjustment for that score.

**Confounder-adjusted model.** Logistic regression `R ~ z(CLDN4) + covariates`, run per cohort where
covariates exist: keratinocyte score (all cohorts), biopsy site category (Hugo), library
strandedness batch (Hugo), sequencing batch (MGH), antibody class (MGH). Reported as a secondary,
clearly-labeled analysis; the primary estimate remains unadjusted.

**Exploratory time-to-event.** Overall survival (Hugo only; days + vital status in GEO metadata)
by continuous `z(CLDN4)` and by median split, via Cox proportional hazards and log-rank. Labeled
exploratory; OS is not available for the other two cohorts.

**Additional pre-specified sensitivity analyses.**

- S1: Riaz `SD` as NR (above).
- S2: median-dichotomized `CLDN4` (high vs low) instead of continuous.
- S3: Riaz analyzed from raw counts (upper-quartile normalized log-CPM) instead of author FPKM.
- S4: MGH analyzed with all baseline treatment episodes rather than the first per patient.
- S5: pooled individual-participant logistic regression with cohort as a fixed effect, replacing the
  two-stage meta-analysis.

## 6. What would count as a positive result

All of the following, stated in advance:

1. Meta-analytic p < 0.05 for `z(CLDN4)`, **and**
2. consistent direction of effect in at least 2 of 3 cohorts, **and**
3. `CLDN4` in the top 5% of the genome-wide meta-analytic ranking, **and**
4. the association not abolished by adjustment for the keratinocyte/epithelial-contamination score.

Anything less is reported as a null or an inconclusive result. No alternative outcome definition,
subgroup, or transformation will be substituted to reach significance; every sensitivity analysis
listed above is reported whatever it shows.

## 7. Known limitations (stated in advance, not after seeing the answer)

- Total baseline sample size across all three open cohorts is on the order of 100 patients, with
  roughly a third responders. This gives limited power: a single-gene continuous predictor needs a
  large effect (roughly OR per SD >= 2, AUC >= 0.68) to be detectable at 80% power at this n.
- The three cohorts differ in therapy (pembrolizumab / nivolumab / mixed anti-CTLA-4 and anti-PD-1),
  response criteria, quantification pipeline (two genome builds and three annotation sets), and
  library chemistry. Batch and pipeline effects are not removable.
- Bulk RNA cannot attribute `CLDN4` to tumor cells; skin biopsies contain epidermis.
- Riaz `SD` handling materially changes the responder definition and is a known lever for
  p-hacking, which is why both variants are reported.
- Prior ipilimumab exposure (a real prognostic stratifier in Riaz) is not in the GEO metadata and is
  therefore not adjusted for.
- No independent validation cohort is held out; all open cohorts are used in the meta-analysis.
