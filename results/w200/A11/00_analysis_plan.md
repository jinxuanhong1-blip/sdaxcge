# A11 - Pre-registered analysis plan

**Question.** In human lung cancer, do the CD47, Galectin, Nectin and TGF-beta
immune-evasion axes differ between TACSTD2 (TROP2)-high and TACSTD2-low tumours?

**Why it matters.** TACSTD2/TROP2 is the target of sacituzumab govitecan and
datopotamab deruxtecan in NSCLC. If TROP2-high tumours carried a distinct
innate-checkpoint profile, that would motivate specific ADC + immunotherapy
combinations. The claim is only useful if it survives the obvious confounder
(below), so the confounder analysis is the primary analysis, not a footnote.

This plan was written and committed **before** any result was inspected.

---

## The confounder that drives the whole design

TACSTD2 is an epithelial gene. In bulk RNA-seq, a "TACSTD2-high" tumour is,
to a first approximation, a tumour with **more epithelial/tumour content**.
Therefore, in bulk data and with no adjustment:

- genes expressed mainly by tumour cells (e.g. NECTIN4) will correlate
  *positively* with TACSTD2, and
- genes expressed mainly by stroma/immune cells (e.g. TGFB1, LGALS1) will
  correlate *negatively* with TACSTD2,

**purely as a composition artefact**, with no biological relationship at all.
Any analysis that reports only the unadjusted comparison will produce a
confident, publishable-looking, and largely meaningless answer.

So A11 is designed as a contrast between:

1. **Naive analysis** - what the standard "high vs low" comparison reports; and
2. **Confounder-adjusted analysis** - the same comparison after removing
   tumour-content variation, plus single-cell confirmation *within tumour cells*.

The reported finding is the difference between (1) and (2).

## Datasets (all public, all downloaded, none simulated)

| Role | Dataset | Scale |
|---|---|---|
| Bulk discovery | TCGA-LUAD (UCSC Xena GDC hub, STAR TPM) | ~590 samples |
| Bulk replication | TCGA-LUSC (same pipeline) | ~553 samples |
| Purity | TCGA ABSOLUTE calls (GDC PanCanAtlas) | 10,787 samples |
| Single cell discovery | GSE131907 (Kim et al. 2020, LUAD) | ~208k cells, 44 samples |
| Single cell replication | GSE127465 (Zilionis et al. 2019, NSCLC) | ~55k cells, 7 patients |

## Pre-specified analyses

### Bulk

- Primary tumours only. Expression = log2(TPM+1). LUAD and LUSC analysed
  separately; LUSC is treated as replication, not pooled.
- **Exposure:** TACSTD2, as (a) continuous and (b) top vs bottom tertile.
- **Outcomes:** the 30-gene primary panel (CD47, Galectin, Nectin, TGFB axes),
  each axis's ligand-only module score, and the F-TBRS TGF-beta *activity*
  signature.
- **A1 (naive):** Spearman rho and tertile Wilcoxon + Cliff's delta, BH-FDR
  across the primary panel.
- **A2 (adjusted):** partial Spearman controlling ABSOLUTE purity; and
  OLS `gene ~ TACSTD2 + purity`. Pre-specified summary = the **attenuation**
  of each effect from A1 to A2.
- **A3:** repeat within purity strata (tertiles of ABSOLUTE purity).
- **A4 (calibration):** the identical pipeline is run on 10 housekeeping
  negative-control genes. Whatever signal appears there is the analysis's
  own composition-driven noise floor, and is subtracted from interpretation.

### Single cell

- **S1:** per-cell-type mean expression and percent-expressing for every panel
  gene. This directly measures how much of each axis is even tumour-cell-derived.
- **S2:** the primary test. Restrict to **malignant/epithelial cells from tumour
  tissue**, aggregate to **patient-level pseudobulk**, and correlate axis genes
  with TACSTD2 across patients. Patient is the unit of analysis; per-cell tests
  are reported only as a descriptive and explicitly flagged as pseudoreplicated.
- **S3:** does TME composition (cell-type fractions) differ by patient TACSTD2?

## Interpretation rules, fixed in advance

- An axis counts as **genuinely associated** with TACSTD2 only if it is
  (i) significant after purity adjustment in bulk LUAD, (ii) directionally
  consistent in LUSC, and (iii) directionally consistent in within-malignant-cell
  single-cell pseudobulk. Anything less is reported as **not established**.
- Effect sizes and confidence intervals are reported alongside every p-value.
  With n~590 the study is powered to detect rho~0.11 at p<0.05, so statistical
  significance alone carries almost no information; |rho| >= 0.3 is the
  pre-specified threshold for "worth acting on".
- Negative and null results are reported with the same prominence as positive
  ones. If the honest answer is "this is mostly a purity artefact", that is
  the headline.
