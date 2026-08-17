# Methods — GSE148071 REAL Slingshot/PAGA, CLDN4 only

ADDITIVE. Does not rewrite `methods/gse148071_paga_cldn4/` (PR #394, PAGA/DPT).
PR #394 barrier/keratin patient Spearman **ρ=+0.590 is given and is not re-audited**.
TACSTD2 does not define groups. No dual-high gate.

**Question.** On public advanced-NSCLC malignant-like epithelium (GSE148071, Wu 2021), where do **CLDN4**, a CLDN4-excluded barrier/keratin program, and an IFN program sit on a **real Slingshot** lineage (PAGA geometry) relative to an Alveolar/AT2 root that is **not** CLDN4-high?

## Data

| Item | Choice |
| --- | --- |
| Cohort | GSE148071 (Wu et al. 2021); 42 advanced NSCLC biopsies |
| Counts | Public GEO raw UMI (`GSE148071_RAW.tar`) |
| Labels | TISCH2 major-lineage ∈ {Malignant, Alveolar, Basal, Epithelial} |
| Unused | TISCH expression.h5; ICI / RECIST / MPR; GEO histology; dual-high TACSTD2∩CLDN4 |

## Trajectory

| Item | Choice |
| --- | --- |
| Graph | log1p(CP10k) → Seurat HVG 3000 → PCA 30 → k-NN 30 |
| Clusters | Leiden 0.6; PAGA on Leiden (geometry) |
| Clock | **Bioconductor slingshot** 2.10 on PCA 1–10, `start.clus` = Leiden of the root cell |
| Comparator | scanpy DPT (same root); not the claim clock |
| Root | TISCH Alveolar, median AT2 score, **CLDN4-high tertile excluded** |
| Cap | max 500 cells/patient if n>20k (seed 0) |
| Score | CLDN4 continuous + tertiles |
| Barrier/keratin | KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**no CLDN4**) |
| IFN | Hallmark IFNα ∩ IFNγ intersection (73 genes; WARS/WARS1 alias) |
| Inferential n | **Patient** (one sample each). Cell-level ρ is descriptive. |
| Extra figures | Along-PT CLDN4/barrier/IFN; within-patient CLDN4-high vs low (min 8 cells/arm) |

## Tests

Primary (BH inside this list only): patient-mean Spearman of CLDN4 / barrier / IFN vs principal Slingshot PT, plus CLDN4 vs IFN, AT2/SFTPC vs PT (root controls), DPT comparator.

Sensitivity (not BH): malignant-only, drop leftover-dominant, per-Slingshot-lineage, mean-across-lineages PT.

PR #394 CLDN4 vs barrier (n=42, ρ=+0.590) is cited, not re-fit as a discovery contrast.

## Done criterion

`tables/lineage_level.tsv` and `tables/patient_level.tsv` exist.
