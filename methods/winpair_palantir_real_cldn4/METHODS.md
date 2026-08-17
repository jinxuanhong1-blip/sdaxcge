# Methods — winning-pair REAL Palantir, CLDN4 only

ADDITIVE. Does **not** rewrite `methods/scrna_paga_cldn4/` (PR #325) or
`methods/winpair_131907_205335_slingshot_cldn4/` (PR #449 DPT fallback).
GSE148071 is not added. GSE207422 is not added. TACSTD2 does not define groups.
No dual-high gate.

**Question.** On the winning public pair GSE131907+GSE205335 **epithelium only**,
do Palantir destinies place **CLDN4-high** at a CLDN4-excluded barrier end and
**CLDN4-low** at an IFN-higher end of the same malignant trajectory?

## Datasets (public processed UMI only)

| Cohort | Citation | Matrix | Epithelium kept | Dropped |
| --- | --- | --- | --- | --- |
| GSE131907 | Kim et al., *Nat Commun* 2020, PMID 32385277 | raw UMI text + author annotation | `Cell_type == Epithelial cells` and `Sample_Origin` ∈ {tLung, nLung, tL/B, mLN, mBrain} | PE unlabeled epithelium; 2.86 GB log2TPM; EGA FASTQ |
| GSE205335 | Ahn / Lee, *eLife* 2024 | UMI dgCMatrix + author identity + SOFT | `lineage.total == Epithelial cells` on non-normal tissues | Normal Lung / LN / Brain; EGA FASTQ |

Winning-pair membership is taken as given from PR #290. This folder asks a
**Palantir destiny** question, not a T/NK re-rank.

## Palantir (required)

Palantir (Setty et al., 2019) is installed and run. DPT is a companion only.
Empty or failed DPT is recorded and **is not a stop**.

| Item | Choice |
| --- | --- |
| Graph | Seurat-v3 HVG 3000 → PCA → Harmony (`batch = dataset`) → k-NN 30 |
| Clusters | Leiden 0.6 (display / early-cell fallback only) |
| Early cell | GSE131907 nLung author AT2, **never CLDN4-high** |
| Terminals | Palantir auto. Not defined by CLDN4, barrier, or IFN |
| Scores | CLDN4 continuous + tertiles; barrier/keratin **without CLDN4**; IFN core **without CLDN4** |
| Inferential n | **Patient** (GSE131907 patient number + GSE205335 patient). Cell-level ρ is descriptive |
| Cap | ≤350 cells / sample after protecting nLung AT2 |
| Extra figures | Within-patient CLDN4-high vs low; PT-decile trends; fate vs programs |
| Unused | GSE148071; GSE207422; dual-high TACSTD2∩CLDN4; ICI / MPR / RECIST as a trajectory label |

## Tests

Primary: patient-level Spearman of mean CLDN4 / barrier / IFN / Palantir PT / fate
probabilities. BH inside that list only.

Sensitivity (not BH): per-cohort; malignant-cell–restricted patient means.

Thesis is tested, not assumed. Destinies are phenotypic terminals, not lineages.
