# Claim A1

## Statement

In non-small-cell lung cancer, tumor-intrinsic **TACSTD2** (TROP2) expression is
**negatively correlated** with the abundance of anti-tumor immune infiltration —
specifically with **immune (T-cell effector)**, **cytotoxic**, and **T-cell
exhaustion / immune-checkpoint** transcriptional signatures — and this negative
association **persists after adjusting for tumor purity** (partial Spearman
correlation controlling for purity).

## As reported by the user

> "still negative after purity, n≈1031"

i.e. the purity-adjusted partial Spearman correlations between TACSTD2 and the
immune/cytotoxic/exhaustion signatures remain **negative** (and statistically
significant), on a combined sample of roughly **n ≈ 1031** tumors drawn from
TCGA and OncoSG lung-cancer cohorts.

## Scope of this replication (A1 only)

- **Public data only.** All expression, purity, and clinical data are pulled
  from public resources (cBioPortal REST API and the TCGA PanCanAtlas ABSOLUTE
  supplemental table on the NCI GDC). No controlled-access data is used.
- **Cohorts**
  - TCGA lung adenocarcinoma — `luad_tcga_pan_can_atlas_2018`
  - TCGA lung squamous cell carcinoma — `lusc_tcga_pan_can_atlas_2018`
  - OncoSG lung adenocarcinoma (Nat Genet 2020) — `luad_oncosg_2020`

  "TCGA" in the claim is interpreted as TCGA NSCLC (LUAD + LUSC), because that is
  the combination whose expression+purity sample count lands near the reported
  n ≈ 1031; per-histology results are also reported so the reader can see every
  slice. This is stated explicitly rather than hidden.
- **Signatures** (see `signatures.json` for exact gene lists and citations):
  immune / T-cell effector, cytotoxic, and T-cell exhaustion / checkpoint.
- **Method:** first-order **partial Spearman correlation** between TACSTD2 and
  each signature score, controlling for tumor purity. Four purity constructions
  are reported: **ESTIMATE** (Aran 2015), **ABSOLUTE** (GDC PanCanAtlas),
  **CPE** (Aran consensus of ESTIMATE+ABSOLUTE+LUMP+IHC), and
  **ESTIMATE-or-ABSOLUTE** (ESTIMATE preferred, ABSOLUTE fallback). OncoSG has
  no public ESTIMATE table; its cBioPortal clinical `PURITY` is used when
  OncoSG is included. Firehose LUAD/LUSC is run as a sensitivity because that
  is the ESTIMATE-complete published set (n = 1014), closest to the claimed
  n ≈ 1031. Reported per cohort, pooled, and via a fixed-effect meta-analysis.

## Prediction (what "replicates" means here)

The claim replicates if, after ESTIMATE and/or ABSOLUTE purity adjustment,
TACSTD2 vs each of the three signatures is **negative**. n ≈ 1031 is treated
as a separate, checkable number — not forced. A sign match with an n mismatch
is reported as such.
