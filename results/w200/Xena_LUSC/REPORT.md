# W200 Xena LUSC — TACSTD2 / CLDN4 vs immune after ESTIMATE purity

**Self-contained. Public data only. LUSC only. Written to be read without the rest of the repo.**

**Verdict: PARTIAL — TACSTD2–CD8 remains; CLDN4 is null; do not cite ImmuneScore**

In TCGA-LUSC, ESTIMATE-purity partial Spearman leaves a negative TACSTD2–CD8 association (and TACSTD2–GEP18 is negative; TACSTD2–CYT is null). CLDN4 is null on CD8, CYT and GEP18 — raw and adjusted. The circular ImmuneScore | ESTIMATE-purity row is positive for both genes and is not an immune-cold claim. ABSOLUTE (DNA) keeps TACSTD2–CD8/CYT/GEP18 negative and stronger; CLDN4 stays null. Methylation leukocyte fraction is not negative.

## Why this analysis exists

Prior A1 work on pooled TCGA NSCLC reported a purity-adjusted negative
Spearman between TACSTD2 and immune / cytotoxic signatures, and the
histology split showed that GEP18 and ESTIMATE ImmuneScore were
LUSC-only. This file is the LUSC-only follow-up with the covariate
the user asked for: **ESTIMATE purity**, and with **CLDN4** tested
the same way as TACSTD2.

It is not a LUAD analysis. It is not an ICI-response analysis. It is
not a claim that TROP2 or claudin-4 *causes* immune exclusion.

## Analysis set

| Filter | n |
|---|---:|
| TCGA-LUSC HiSeqV2 primary tumors (`-01`, 15-char barcode) | 502 |
| With official ESTIMATE RNAseqV2 scores | 501 |
| With Yoshihara ESTIMATE purity (unclipped) | 501 |
| ESTIMATE purity clipped to [0, 1] | 0 |
| With PanCanAtlas ABSOLUTE purity | 493 |
| With methylation leukocyte fraction | 501 |
| Complete for primary ESTIMATE-purity tests (gene + feature + ESTIMATE purity) | 501 |

Primary tumors only. One row per 15-character barcode. Replicate
aliquots that collapse to the same `-01` barcode are averaged.

## Pre-specified features

| Name | Definition | Honest limitation |
|---|---|---|
| **TACSTD2** | Xena HiSeqV2 log2(RSEM+1) | RNA, not protein / IHC. |
| **CLDN4** | same matrix | RNA, not tight-junction function. |
| **CD8** | `CD8A` | Single gene, not a deconvolution fraction. |
| **CYT** | mean(log2 GZMA, log2 PRF1) = Rooney 2015 | Tracks cytotoxic mRNA, not killing. |
| **GEP18** | unweighted mean of within-LUSC z-scores of the 18 Ayers 2017 genes | Merck NanoString weights are **not public**. This is the open surrogate. |
| **ESTIMATE** | official **ImmuneScore** (MD Anderson RNAseqV2 table) | Built from immune genes. **Algebraically circular** with ESTIMATE purity. Not a primary claim under this covariate. |
| **LEUK** | DNA-methylation leukocyte fraction (PanImmune) | Orthogonal to RNA. Different missingness. |

**Primary covariate:** Yoshihara ESTIMATE tumor purity,
`cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)`, from the official
MD Anderson LUSC RNAseqV2 table. We do not recompute ESTIMATE.

ESTIMATE purity is systematically higher than ABSOLUTE in this set
(median ~0.80 vs ~0.51). That is a known scale difference between the
two algorithms, not a merge error. Ranks, not absolute purity values,
enter the partial Spearman.

**Rank identity (important):** in the observed ESTIMATEScore range the
cosine transform is strictly decreasing, so ESTIMATE purity and
ESTIMATEScore have identical reversed ranks. Partial Spearman of X,Y
| ESTIMATE_purity **equals** partial Spearman of X,Y | ESTIMATEScore.
Calling the covariate 'purity' does not add information beyond the
ESTIMATE score itself.

**Sensitivity covariate:** PanCanAtlas ABSOLUTE (DNA). Method: first-order
partial Spearman (algebraic formula, same as rework A1). Residual-rank
partial Spearman is a sensitivity column. 95% CIs are Fisher-z with
variance `1/(n-4)` for partial correlations. BH-FDR is across the 8
primary tests (2 genes × CD8/CYT/GEP18/ESTIMATE, ESTIMATE-purity partial).

## Circularity / overlap (locked)

ESTIMATE immune signature size: 141 genes (from `ESTIMATE_SI_geneset.gmt`).

| Feature | Genes overlapping the ESTIMATE immune signature |
|---|---|
| CD8 (`CD8A`) | none |
| CYT | PRF1 |
| GEP18 | CCL5, CD27, HLA-E, NKG7 (4/18) |
| ESTIMATE ImmuneScore | the score **is** that signature |
| LEUK | none (methylation, not a gene set) |

If CD8A / GZMA / PRF1 / GEP18 genes sit inside the ESTIMATE immune
signature, then 'after ESTIMATE purity' is partly 'after a score
that already contains the endpoint'. ABSOLUTE does not have that
problem. LEUK does not have that problem.

## Direct answer

In TCGA-LUSC, ESTIMATE-purity partial Spearman leaves a negative TACSTD2–CD8 association (and TACSTD2–GEP18 is negative; TACSTD2–CYT is null). CLDN4 is null on CD8, CYT and GEP18 — raw and adjusted. The circular ImmuneScore | ESTIMATE-purity row is positive for both genes and is not an immune-cold claim. ABSOLUTE (DNA) keeps TACSTD2–CD8/CYT/GEP18 negative and stronger; CLDN4 stays null. Methylation leukocyte fraction is not negative.

- **TACSTD2** is the gene with a residual CD8 signal. It is not a
  general cytotoxic / GEP / leukocyte-fraction finding under ESTIMATE
  purity. ABSOLUTE (DNA) makes TACSTD2 look more immune-cold than
  ESTIMATE purity does, because ESTIMATE purity already contains the
  immune score.
- **CLDN4** does not show an immune-cold pattern in LUSC on these
  readouts, despite ρ≈0.41 with TACSTD2. Do not treat the two genes
  as interchangeable immune correlates.
- **ImmuneScore | ESTIMATE purity** is circular and here it is
  *positive*. Citing it as 'immune-cold after purity' would be wrong.
- **LEUK** (methylation) is not negative. The TACSTD2–CD8 result is
  an RNA-signature result, not a DNA-methylation leukocyte result.

## Primary result — ESTIMATE-purity partial Spearman

| Gene | CD8 | CYT | GEP18 | ESTIMATE ImmuneScore (circular) |
|---|---|---|---|---|
| TACSTD2 | -0.166 (p=2.01e-04, n=501) | -0.062 (p=0.163, n=501) | -0.114 (p=0.011, n=501) | 0.183 (p=3.81e-05, n=501) |
| CLDN4 | 0.041 (p=0.365, n=501) | -0.015 (p=0.741, n=501) | 0.052 (p=0.242, n=501) | 0.181 (p=4.68e-05, n=501) |

Unadjusted (no purity) for attenuation:

| Gene | CD8 | CYT | GEP18 | ESTIMATE ImmuneScore |
|---|---|---|---|---|
| TACSTD2 | -0.199 (p=7.02e-06, n=502) | -0.126 (p=0.005, n=502) | -0.158 (p=3.75e-04, n=502) | -0.055 (p=0.216, n=501) |
| CLDN4 | -0.017 (p=0.699, n=502) | -0.055 (p=0.220, n=502) | -0.024 (p=0.598, n=502) | -0.004 (p=0.925, n=501) |

ABSOLUTE-purity partials are in the next section and in
`correlations.tsv`. Per-row FDR for the 8 primary tests is
`partial_fdr_8tests`.

## Orthogonal checks

### Methylation leukocyte fraction (not RNA)

| Gene | Unadjusted vs LEUK | ESTIMATE-purity partial vs LEUK | ABSOLUTE partial vs LEUK |
|---|---|---|---|
| TACSTD2 | 0.065 (p=0.146, n=501) | 0.158 (p=4.06e-04, n=500) | 0.012 (p=0.795, n=492) |
| CLDN4 | 0.021 (p=0.634, n=501) | 0.069 (p=0.126, n=500) | 0.021 (p=0.646, n=492) |

### ABSOLUTE (DNA) purity partial — same immune features

| Gene | CD8 | CYT | GEP18 | ESTIMATE ImmuneScore |
|---|---|---|---|---|
| TACSTD2 | -0.244 (p=4.30e-08, n=493) | -0.172 (p=1.27e-04, n=493) | -0.226 (p=3.95e-07, n=493) | -0.131 (p=0.003, n=493) |
| CLDN4 | -0.031 (p=0.488, n=493) | -0.074 (p=0.099, n=493) | -0.046 (p=0.306, n=493) | -0.019 (p=0.675, n=493) |

### Do ESTIMATE-purity and ABSOLUTE-purity partials agree?

ESTIMATE purity vs ESTIMATEScore Spearman ρ=-1.000 (p=<1e-300, n=501); expected ≈ −1 if the cosine transform is strictly monotone in this cohort. 6/10 gene×feature calls (sign + p<0.05) agree between ESTIMATE-purity partial and ABSOLUTE partial.

- TACSTD2 vs CD8: ESTIMATE-adj negative ρ=-0.166; ABSOLUTE-adj negative ρ=-0.244
- TACSTD2 vs CYT: ESTIMATE-adj null ρ=-0.062; ABSOLUTE-adj negative ρ=-0.172  **call differs**
- TACSTD2 vs GEP18: ESTIMATE-adj negative ρ=-0.114; ABSOLUTE-adj negative ρ=-0.226
- TACSTD2 vs ESTIMATE: ESTIMATE-adj positive ρ=0.183; ABSOLUTE-adj negative ρ=-0.131  **call differs**
- TACSTD2 vs LEUK: ESTIMATE-adj positive ρ=0.158; ABSOLUTE-adj null ρ=0.012  **call differs**
- CLDN4 vs CD8: ESTIMATE-adj null ρ=0.041; ABSOLUTE-adj null ρ=-0.031
- CLDN4 vs CYT: ESTIMATE-adj null ρ=-0.015; ABSOLUTE-adj null ρ=-0.074
- CLDN4 vs GEP18: ESTIMATE-adj null ρ=0.052; ABSOLUTE-adj null ρ=-0.046
- CLDN4 vs ESTIMATE: ESTIMATE-adj positive ρ=0.181; ABSOLUTE-adj null ρ=-0.019  **call differs**
- CLDN4 vs LEUK: ESTIMATE-adj null ρ=0.069; ABSOLUTE-adj null ρ=0.021

## TACSTD2 vs CLDN4 in this LUSC set

Raw Spearman TACSTD2 vs CLDN4: ρ=0.410 (p=8.93e-22, n=502). ESTIMATE-purity partial: ρ=0.402 (p=6.91e-21, n=501). ABSOLUTE partial: ρ=0.403 (p=1.32e-20, n=493). They are correlated in LUSC; they are not interchangeable. Every immune test is run on each gene separately.

## What this is not

- Not LUAD, not pooled NSCLC, not OncoSG.
- Not protein, IHC, spatial, or single-cell.
- Not ICI response, DCB, ORR, or survival (survival file was not used).
- Not a claim that ESTIMATE 'purity adjustment' is independent of the
  immune score. It is not. See circularity section.
- Not a re-implementation of ESTIMATE; official MD Anderson tables only.
- CD8B-mean and StromalScore rows are secondary and in the TSV only.

## How to rerun

```
pip install -r requirements.txt
python scripts/w200_xena_lusc.py
```

Downloads land in `data/` (gitignored). This report is rewritten from
the tables on every run.

## Provenance

See `provenance.json` for URL, bytes and sha256 of every input.

