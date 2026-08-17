# Methods — winning-pair DoRothEA TF activity (CLDN4-only)

ADDITIVE. **CLDN4 only.** TACSTD2 does not define groups (no dual-high).
Patient is the unit. **GSE131907 + GSE205335 only.** GSE207422 and
GSE148071 are not loaded, not scored, and not pooled.

decoupleR / dorothea R were not installed. TF activity is the documented
DoRothEA **weighted mean (wmean)** from Badia-i-Mompel et al. 2022.

## Datasets (public processed UMI only)

| Cohort | Citation | Matrix | Not used |
| --- | --- | --- | --- |
| GSE131907 | Kim et al., *Nat Commun* 2020, PMID 32385277 | `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz` + author annotation + series matrix | 2.86 GB log2TPM text; EGA FASTQ |
| GSE205335 | Hu et al., palliative ICI biopsy/effusion scRNA | `GSE205335_Lung_IO_UMI_matrix.rds.gz` (dgCMatrix) + author identity + GEO SOFT | EGA raw |

No FASTQ. TPM text is skipped because UMI exists.

## Compartments (author labels)

| Cohort | Malignant cells used |
| --- | --- |
| GSE131907 | `Cell_type == Epithelial cells` and `Cell_subtype ∈ {Malignant cells, tS1, tS2, tS3}` on tumor origins `tLung, tL/B, mLN, mBrain` |
| GSE205335 | `lineage.sub == Malignant cells` |

Honest label limits:

- Primary tLung epithelium is annotated **tS1 / tS2 / tS3**, not `Malignant cells`. Those tS* barcodes are included.
- Author `Malignant cells` in GSE131907 sit in tL/B, mLN, mBrain.
- PE unlabeled epithelium, nLung, and nLN are excluded.
- CopyKAT was not re-run.
- GEO `patient_id` (GSE131907) or GEO `patient` (GSE205335) is the unit. Same-patient tumor-origin samples are pooled.

T/NK is **not** required for this TF test (malignant cells only).

## CLDN4 splits (malignant cells, within patient)

Score = `log1p(CP10k)` CLDN4 on that patient's malignant cells.

| Split | High | Low |
| --- | --- | --- |
| **Primary: Q4 vs Q1** | equal-count top quartile after a stable sort on CLDN4 (`n//4` cells) | bottom quartile (`n//4`) |
| Companion: median | ≥ median (if median is 0: CLDN4>0 vs =0) | < median (or =0) |

`pd.qcut` on ranks is **not** used. Zero-inflated CLDN4 (many exact zeros)
collapses rank bins and would drop patients with hundreds of malignant
cells. Equal-count tails keep those patients. The high tail must have a
strictly higher mean CLDN4 than the low tail.

**Honest paired n.** A patient is scored only if:

- malignant cells ≥ 20
- each compared bin ≥ 8 cells (Q4/Q1) or ≥ 10 cells (median)

Patients missing a tail are out. Cells are not n. Thin n (pooled n<8) is flagged, not hidden. p-values are descriptive.

## DoRothEA network

Human A+B+C edges from OmniPath (`resources/dorothea_hs_ABC.tsv`;
Garcia-Alonso et al. 2019). `mor` = +1 activation / −1 repression.
Unsigned edges default to +1 (decoupleR convention). A TF is scored
only if ≥5 of its ABC targets are present in **both** UMI matrices.

## TF activity (documented wmean; decoupleR not run)

Per cell, for TF \(t\) with signed weights \(w_g\) on targets \(g\):

\[
\mathrm{wmean}_{t,c} = \frac{\sum_g w_g \, x_{g,c}}{\sum_g |w_g|}
\]

where \(x_{g,c} = \log1p(\mathrm{CP10k})\) of target \(g\) in cell \(c\).
This is `decoupleR::run_wmean` / `decoupler.run_wmean` without the
permutation NES. Patient-bin activity = mean of per-cell wmean in that
bin. Within-patient high vs low does not need a cross-patient z-score.

## Focused programs (primary table)

TFs are pre-specified, not data-mined. A TF is listed only if it is in
DoRothEA A+B+C. Missing from A+B+C (not imputed): CIITA, NLRC5, GRHL1,
GRHL3, OVOL1/2, RFXANK, RFXAP.

| Program | TFs |
| --- | --- |
| IFN | STAT1, STAT2, IRF1, IRF2, IRF3, IRF7, IRF8, IRF9, STAT3 |
| MHC | RFX5, IRF1 (shared with IFN) |
| TJ | GRHL2, KLF4, ELF3, TFAP2A |
| keratin | TP63, KLF5, SOX2 |

CLDN4 is the splitter and is never a scored target-set member for the
companion gene-set modules.

## Companion gene-set modules (not TF activity)

Mean `log1p(CP10k)` of public lists (CLDN4 excluded from TJ):

- IFN ISG core (40 genes)
- MHC-I APM
- TJ (no CLDN4)
- keratin (simple + basal; genes absent from a matrix are dropped)

These are extra, not the TF table.

## Tests

**Primary.** Wilcoxon signed-rank (two-sided) on per-patient
\(\mathrm{wmean}_{\mathrm{Q4}} - \mathrm{wmean}_{\mathrm{Q1}}\) for each
focused TF. Report GSE131907, GSE205335, and the additive pool
(patients concatenated; not a cell merge). Minimum paired n after the
floor = 6. Exact p is used when n ≤ 25.

**Companion.** Same test on the median split; patient-mean CLDN4 vs
patient-mean TF Spearman (all floor-eligible patients, not tails).

**Not done.** VIPER NES, permutation NES, pySCENIC, dual-high
TACSTD2∩CLDN4, GSE207422, GSE148071.

## Software

Python (numpy / pandas / scipy / matplotlib / rdata). No R CellChat,
no decoupleR, no dorothea Bioconductor package.
