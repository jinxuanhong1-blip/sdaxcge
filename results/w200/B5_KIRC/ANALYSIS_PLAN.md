# B5 analog — KIRC/RCC: CLDN4 expression vs ICI response

**Pre-specified analysis plan. Written and committed BEFORE any association statistics were run.**
(Only data inventory — cohort sizes, response-category counts, CLDN4 expression range — was inspected first,
to confirm the analysis was feasible at all. No test of CLDN4 against outcome was run before this commit.)

## Question

In advanced clear cell renal cell carcinoma (ccRCC / TCGA "KIRC" lineage), is pre-treatment tumor
`CLDN4` mRNA expression associated with response to immune checkpoint inhibition (ICI)?

This is the kidney-cancer analog of project claim **B5** (CLDN4 vs ICI response / odds ratio).

## Datasets (all fully public, no controlled access)

| Cohort | Treatment | Assay | n (RNA) | Outcome available |
|---|---|---|---|---|
| Braun et al. 2020, *Nat Med* — CheckMate-009/010/025 | nivolumab (anti-PD-1); everolimus control arm | bulk RNA-seq | 311 (181 nivo, 130 evero) | RECIST best response (ORR), clinical benefit, PFS, OS |
| Motzer et al. 2020, *Nat Med* — JAVELIN Renal 101 | avelumab+axitinib; sunitinib control arm | bulk RNA-seq TPM | 726 | PFS only (**no response data in the public supplement**) |
| Ascierto et al. 2016 — GSE67501 | nivolumab | microarray | 11 | responder / non-responder |

Braun is the **primary** cohort because it is the only public RCC ICI cohort with both transcriptome and
RECIST response. JAVELIN is a PFS-only replication that additionally allows a randomized
treatment-by-biomarker interaction test. GSE67501 is a severely underpowered exploratory check.

## Primary analysis (one pre-specified test)

- **Population:** Braun nivolumab-treated patients with pre-treatment RNA (n = 181).
- **Endpoint:** binary objective response — responder = best response CR/PR (recorded as `CR`, `PR`, `CRPR`)
  vs non-responder = `PD`. Patients with `SD` or `NE` are excluded from the primary endpoint.
  Expected n ≈ 39 responders vs 69 PD.
- **Test:** logistic regression `response ~ z(CLDN4)`, where CLDN4 is the log2 expression value supplied by
  the authors, z-scored within the analysis population. Report odds ratio per +1 SD, Wald 95% CI, two-sided p.
- **Significance:** α = 0.05, two-sided, single primary test.
- **Direction:** the project hypothesis is that CLDN4-high tumors respond *worse* (OR < 1). The test is
  nonetheless two-sided and the observed direction will be reported as-is.

## Secondary / sensitivity analyses (all clearly labeled as secondary; no alpha protection claimed)

1. Mann-Whitney U on CLDN4 between responders and non-responders; AUC (= rank-biserial equivalent).
2. Alternative response definitions: (a) CR/PR vs SD+PD (all evaluable, n ≈ 172); (b) Braun's
   clinical-benefit variable, CB vs NCB (excluding ICB).
3. Median-split and tertile-split 2×2 odds ratios (Fisher exact, conditional MLE OR + CI).
4. Multivariable logistic model adjusting for trial cohort (CM-009/010/025) and tumor purity
   (CLDN4 is epithelial, so purity is a plausible confounder), and for MSKCC risk group.
5. **Control arm (specificity):** identical analysis in the everolimus arm. If CLDN4 associates with
   outcome in both arms, it is prognostic rather than ICI-predictive.
6. Time-to-event: Cox PH for PFS and OS on z(CLDN4), within each arm.
7. **JAVELIN replication:** Cox PFS on z(CLDN4) in the avelumab+axitinib arm; same in the sunitinib arm;
   and a randomized arm × CLDN4 interaction test on PFS.
8. **GSE67501:** Mann-Whitney CLDN4 responder vs non-responder (n = 11; reported with explicit power caveat).
9. Correlation of CLDN4 with tumor purity, and with CD8A/immune markers, for interpretation.

## Honesty controls (pre-specified)

- **Transcriptome-wide null calibration.** The primary logistic test will be run for *every* gene in the
  Braun matrix. CLDN4's percentile rank among all genes and its Benjamini-Hochberg FDR in that
  transcriptome-wide context will be reported. A gene that looks "significant" at p<0.05 but sits in the
  middle of the genome-wide p-value distribution is not evidence of anything.
- **Minimum detectable effect.** Power for the primary test will be computed so that a null result can be
  interpreted (i.e. what effect sizes this cohort could and could not have detected).
- **Positive/negative reference genes.** The identical pipeline will be run on a small pre-named panel
  (CD8A, GZMB, PDCD1, CD274, PBRM1, EPCAM, TACSTD2, CLDN3, CLDN7, ACTB) purely to show what this
  pipeline produces on genes with and without expected signal in this cohort. These are context, not claims.
- **All results reported regardless of direction or significance.** If the answer is null, the summary will
  say so plainly. No endpoint will be swapped for another after seeing results; any post-hoc analysis added
  after unblinding will be explicitly labeled POST-HOC in the summary.

## Known limitations (stated up front)

- Only one public RCC ICI cohort has RECIST response with transcriptome; there is no independent public
  RCC ICI *response* replication cohort of meaningful size. JAVELIN can only replicate on PFS, and its
  ICI arm is a combination (avelumab + axitinib), so a PFS effect there is not attributable to ICI alone.
- Braun expression values are ComBat batch-corrected, upper-quartile-normalized log2 TPM with an additive
  constant, as distributed by the authors; they are not raw counts and cannot be re-normalized.
- Bulk expression confounds tumor content with tumor-cell-intrinsic expression.
- CheckMate-009/010/025 are heavily pre-treated, second-line-plus populations; JAVELIN is first-line.
