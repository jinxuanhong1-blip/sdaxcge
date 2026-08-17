# FINDING — pairwise combo: malignant CLDN4 vs T/NK (GSE123902 + GSE131907)

**Additive. CLDN4 only. Not a pile.** Two public processed LUAD scRNA atlases, patient as the unit. No TACSTD2∩CLDN4 dual-high gate. No other cohorts.

GSE123904 is the SuperSeries. Human matrices are GSE123902 (Laughney et al., *Nat Med* 2020, PMID 32042191). Mouse GSE123903 is not used. GSE131907 is Kim et al., *Nat Commun* 2020, PMID 32385277.

## Public-processed gate (<2 GB)

| File | Public? | Size | Used? |
|---|---|---:|---|
| GSE123902_RAW.tar (17 dense UMI CSVs) | yes | **90.4 MB** | **yes** |
| Author `PATIENT_LUNG_ADENOCARCINOMA_ANNOTATED.h5` (dpeerlab S3) | yes | **36.5 GB** | **no** — over the 2 GB gate |
| TISCH `LUAD_GSE123902` / `NSCLC_GSE123902` | **no** | — | 404 |
| GSE131907 raw UMI `txt.gz` + author annotation | yes | **389.8 MB + 1.8 MB** | **yes** |
| GSE131907 log2TPM `txt.gz` | yes | **2.9 GB** | **no** — UMI exists |
| GSE123904 SuperSeries RAW.tar | yes | 211.6 MB | not needed (human is in GSE123902) |

Nothing required for this pair is missing. Analysis did not stop.

## Combo table — patient-level malignant CLDN4 vs same-patient T/NK

Spearman of malignant-cell mean `log1p(CP10k)` CLDN4 vs T/NK fraction. Eligible if ≥20 malignant and ≥20 T/NK cells. Combo is Fisher-z DerSimonian–Laird (k=2). Machine table: `results/combo_cldn4_tnk.tsv`. Forest: `figures/forest_CLDN4_tnk.png`.

| cohort | assay | n | CLDN4 vs T/NK ρ [95% CI] | p | I² | note |
|---|---|---:|---|---:|---:|---|
| GSE123902 | GEO dense UMI; marker epithelial in tumor | **13** | −0.132 [−0.637, +0.452] | 0.668 | — | 13/13 tumor patients eligible |
| GSE131907 | GEO UMI; author malignant + tS1/tS2/tS3 | **31** | −0.327 [−0.610, +0.031] | 0.073 | — | 31/32 tumor-origin patients; P0009 has 5 malignant cells |
| **RE combo** | Fisher-z DerSimonian–Laird | **44** | **−0.277 [−0.539, +0.033]** | **0.079** | **0%** | not a pile |

Stouffer *z* = −1.76 (p = 0.078). Fisher combined p = 0.196 (direction-agnostic). Fixed-effect ρ = RE ρ (τ² = 0, Q = 0.314).

Both arms are negative. They do **not** differ by the pre-specified rule (same sign; |Δρ| = 0.195 < 0.25; I² = 0%; neither arm p < 0.05). **CellChat-style was not run.**

## Honest n

### GSE123902 (Laughney human)

- GEO CSVs: **42,847** barcodes / 17 samples / 14 LX IDs (paper QC atlas is 41,384; we do not substitute that n).
- Tumor samples: 8 primary + 5 metastasis = **13 patients**, **30,126** cells. LX685 is normal-only and is out of the Spearman.
- Matched normals (LX675/682/684/685) are not scored as malignant.
- Marker epithelial in tumor (EPCAM/KRT8/KRT18/KRT19/KRT7 vs T/NK/myeloid/B; CLDN4 not used to assign): **4,433** cells.
- T/NK in those 13 tumors: **13,305** cells.
- Eligible Spearman n = **13** (all 13 have ≥20 malignant and ≥20 T/NK).
- Author cell-type H5 was not used (36.5 GB). This is **not** a CNV-malignant call.

Primary-only n=8 ρ=−0.238 (p=0.570). Mets-only n=5 ρ=−0.500 (p=0.391). Both are underpowered; they are not a second combo.

### GSE131907 (Kim)

- UMI: **208,506 × 29,634**. Series: 58 samples / 44 patients. Treatment-naive. **No ICI / MPR / RECIST labels.**
- Tumor-origin only: tLung + tL/B + mLN + mBrain. PE, nLung, and nLN are out (PE epithelium is not author-malignant).
- Malignant = author `Cell_subtype` ∈ {Malignant cells, tS1, tS2, tS3}. tLung has almost no `Malignant cells` label; tS* is the primary-tumor epithelial call used in the GSE131907 CellChat slice.
- Catalog tumor-origin patients = **32**. Eligible = **31**. Dropped: P0009 (tLung, 5 malignant / 1,981 T/NK).
- Eligible cells: **31,131** malignant + **32,760** T/NK. Site mix among the 31: tLung 10, mBrain 10, mLN 7, tL/B 4.
- Do not write n=44 for this test.

## Restriction (not the infiltration test)

Malignant CLDN4 is higher than same-patient T/NK CLDN4 in **13/13** GSE123902 patients (median 0.75 vs 0.011) and **30/31** eligible GSE131907 patients (median 1.35 vs 0.046). P1013 is the exception (CLDN4 almost empty in both compartments). That is epithelial restriction. It is not an immune-cold claim.

Figures: `figures/scatter_gse123902.png`, `figures/scatter_gse131907.png`, `figures/paired_restriction_gse123902.png`, `figures/paired_restriction_gse131907.png`.

## What is not claimed

- The pair does **not** support a significant anti-correlation. Pooled p = 0.079; the 95% CI includes 0.
- GSE123902 n=13 is small; the CI is wide. Marker epithelial ≠ author / inferCNV malignant.
- GSE131907 is mets-heavy once tS* + `Malignant cells` are taken together. Cell-pool dominance (EBUS_28, mBrain) is a problem for cell-pooled tests; this Spearman is patient-level and does not pool cells across patients.
- No dual-high. No CellChat / LIANA / NicheNet on this pair (the two rhos do not differ).
- No ICI response, MPR, or survival test. Both atlases are treatment-naive (GSE123902 has one post-neoadjuvant primary in the paper; GEO CSVs do not carry a usable RECIST column).
- Not a multi-cohort pile. Leftover lung scRNA from other agents is not added here.

## Reproduce

```bash
python3 methods/pair_123902_131907_cldn4/download.py --out /tmp/pair_123902_131907
tar -xf /tmp/pair_123902_131907/gse123902/GSE123902_RAW.tar -C /tmp/pair_123902_131907/gse123902
python3 methods/pair_123902_131907_cldn4/analyze.py --skip-cellchat
```

`--skip-cellchat` is the default outcome for this pair. CellChat-style code is in `analyze.py` and runs only if the two rhos differ.
