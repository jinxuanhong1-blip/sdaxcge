# Result

Neither public urothelial atezolizumab RNA cohort supports **TACSTD2** or
**CLDN4** as a pretreatment response marker. The well-powered cohort
(IMvigor210) is a **null**, and its CLDN4 median-split odds ratio goes
**opposite** the B5 claim that CLDN4-high patients respond worse (claimed
OR ≈ 0.42). Snyder 2017 is too small to decide the question.

## Primary endpoint

Pretreatment `log2(CPM + 1)` versus RECIST v1.1 **CR/PR vs SD/PD**.
Not-evaluable cases were excluded. AUC is for responders having *higher*
expression. OR is response odds in the expression-high half versus the
expression-low half (cohort-internal median). Bootstrap 95% CIs use 2000
resamples.

| Cohort | Gene | R / NR | Δ median | AUC (R higher) | 95% CI | MW p | OR (high vs low) | 95% CI | Fisher p |
|---|---|---:|---:|---:|---|---:|---:|---|---:|
| IMvigor210 | TACSTD2 | 68 / 230 | +0.14 | 0.556 | 0.471–0.634 | 0.161 | 1.47 | 0.85–2.53 | 0.214 |
| IMvigor210 | CLDN4 | 68 / 230 | +0.35 | 0.553 | 0.467–0.625 | 0.188 | 1.59 | 0.92–2.74 | 0.129 |
| Snyder 2017 | TACSTD2 | 7 / 14 | +0.49 | 0.531 | 0.250–0.812 | 0.856 | 1.33 | 0.21–8.29 | 1.00 |
| Snyder 2017 | CLDN4 | 7 / 14 | −0.45 | 0.449 | 0.176–0.735 | 0.743 | 0.56 | 0.09–3.52 | 0.659 |

IMvigor210 point estimates are slightly *higher* expression in responders.
The IMvigor210 CLDN4 OR confidence interval **excludes 0.42**. A large
adverse CLDN4-high effect of the size stated in the B5 claim is not present
in this public trial RNA set.

Snyder’s CLDN4 point OR (0.56) sits nearer the claimed direction, but it is
3 high-responders / 8 high-nonresponders versus 4 / 6 in the low half. The
interval runs from 0.09 to 3.5. That is noise, not confirmation.

## The assay is not broken

The same IMvigor210 matrix and the same response labels recover the expected
immune-gene direction:

| Gene | AUC | MW p | Notes |
|---|---:|---:|---|
| CXCL9 | 0.646 | 0.00026 | strongest of the prespecified controls |
| CD8A | 0.590 | 0.024 | |
| GZMA/PRF1 mean (CYT) | 0.579 | 0.048 | |

If TACSTD2 or CLDN4 had a CXCL9-sized effect here, the test would have seen
it. They do not.

## Overall survival

IMvigor210 (n=348, 232 deaths), Cox per +1 log2(CPM+1):

- TACSTD2 HR 0.98 (0.90–1.07), p=0.70
- CLDN4 HR 0.97 (0.88–1.08), p=0.61
- CD8A HR 0.88 (0.81–0.97), p=0.006 (control works)

Median-split log-rank: TACSTD2 p=0.14, CLDN4 p=0.37. Snyder OS (n=25, 17
deaths) is null for both genes (HR ≈ 0.90–0.92, p>0.5) and for CD8A.

## Context, not a rescue analysis

TACSTD2 and CLDN4 track each other (Spearman ρ=0.59 in IMvigor210, 0.65 in
Snyder) and are higher in TCGA luminal subtypes I/II than in III/IV. That is
biology, not a response effect. Within-subtype and within-immune-phenotype
tests were exploratory. Almost all were null. The single raw p<0.05
(TACSTD2 inside subtype I, MW p=0.028) is one cell in a post-hoc stratum
grid and is not a result.

Switching IMvigor210 to DESeq size-factor-normalized `log2(count/sf + 1)`,
or dropping SD and comparing CR/PR versus PD only, does not change the
primary conclusion.

Snyder’s paper-defined clinical-benefit label (9 vs 16, including two
durable SD and counting four non-evaluable patients as no-benefit) gives
TACSTD2 AUC 0.67, p=0.18. That is still not significant, uses a different
endpoint, and should not be quoted as a positive finding.

## What this does *not* show

- It does not test protein, IHC, or TROP2-ADC activity.
- It does not test whether either gene is prognostic in untreated
  urothelial cancer.
- It does not reproduce an 11-cohort CLDN4 meta-analysis. It only asks
  whether the two public urothelial anti–PD-L1 RNA sets that were named
  for this extra look support the claim. They do not.
- Snyder n=21 evaluable cannot rule in or rule out a modest effect.

Bottom line: public IMvigor210 RNA is the honest answer for this pair of
genes in urothelial ICI. The answer is no association. Snyder does not
override that.
