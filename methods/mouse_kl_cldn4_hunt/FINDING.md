# Mouse KL / Cldn4 public-data hunt (CLDN4-only)

Additive catalog of **public** mouse lung scRNA/bulk with **Cldn4** plus immune or GEMM labels (KL, KP, Kras, Stk11/Lkb1, LLC, TISMO). Private 8 KL mice were not used. FASTQ was not downloaded.

## What is done

- Catalog: `methods/mouse_kl_cldn4_hunt/catalog.tsv` (23 rows including exclusions).
- Scored mouse-level tables (Cldn4 vs T/NK and epithelial IFN/MHC):
  - `analysis/GSE154989/mouse_level_scores.tsv` (K/KP plate scRNA; **n=28** mice with ≥10 cells)
  - `analysis/GSE6135/mouse_level_scores.tsv` (KL vs K/KP/Ink4a microarray; **n=20** mice / **n=25** tumors)
  - `analysis/E-MTAB-5311/mouse_level_scores.tsv` (LLC bulk, TISMO-class; **n=16** mice)

Assigned-elsewhere accessions (GSE76628, GSE179502, GSE267321, GSE274477) are cataloged only.

## Hunt verdict (usable vs not)

| Accession | Why keep / drop |
|---|---|
| **GSE154989** | Best unused GEMM scRNA: Cldn4 present, processed h5, **30 mice / 28 scored**. K and KP, not KL. |
| **GSE6135** | Best unused **KL** bulk: Cldn4 probe present, series matrix, **13 KL vs 12 other** tumors. |
| **E-MTAB-5311** | Best unused **LLC** bulk with processed Atlas counts (**n=16**). TISMO website dump was JS-gated; this is the public LLC matrix that was actually scored. |
| GSE127465 mouse | Processed MTX + T/NK labels, **n=4** mice, but **CD45+ only** so Cldn4 is the wrong compartment. Not scored. |
| GSE165641 | KL mixed-lineage scRNA, processed Rdata, **n=2**. Too small for high-vs-low. |
| GSE180963 | K vs KL scRNA, **n=2**. Too small. |
| GSE244452 | KL vs KP bulk **n=6**, DEG-only (no counts). No FASTQ. |
| GSE137669 | **Heart**, not lung. |
| GSE131907 | **Human only**. |
| GSE50927 | Cldn4 KO lung injury; DE tables only; Cldn4 is deleted, not a continuous tumor trait. |
| TISMO portal | Immune labels exist; full expression dump not retrieved without FASTQ reprocessing. LLC scored via E-MTAB-5311 instead. |

## Methods (what was scored)

**Cldn4 only** (no other claudins).

**T/NK score:** mean z-score of present genes among Cd3d/e/g, Cd2, Cd8a/b1, Cd4, Nkg7, Gzma/b, Prf1, Klrb1c, Ncr1, Klrd1, Klrc1, Ifng.

**IFN/MHC score:** mean z-score of present Stat1/2, Irf1/7/9, Isg15, Ifit1/2/3, Mx1, Oasl2, Rsad2, Ifih1, Ddx58, Ifnb1, B2m, H2-K1/D1/Q4/Q6/Q7, H2-Aa/Ab1/Eb1, Tap1/2, Psmb8/9, Nlrc5, Ciita.

**Unit:** one row per biological mouse (tumors from the same mouse collapsed). GSE6135 also reports tumor-level n because some mice contribute 2–3 tumors. High vs low is a **median split** of Cldn4 (or detected vs undetected when Cldn4 is zero-inflated). Tests: Spearman and two-sided Mann–Whitney / Welch. No multiple-testing theater.

Processed inputs only: GEO h5 / series matrix / Expression Atlas `E-MTAB-5311-raw-counts.tsv`.

## Results (honest n)

### 1. GSE154989 — K/KP epithelium (n=28 mice)

Plate-seq of GEMM lung tumor cells (Marjanovic et al.). Cldn4 is present. T/NK genes are expressed at low level in this cancer-cell matrix; this is **not** a T/NK fraction.

| Endpoint | n | Spearman r | p | High-vs-low MWU p |
|---|---:|---:|---:|---:|
| T/NK gene score | 28 | +0.12 | 0.53 | 0.77 |
| IFN/MHC gene score | 28 | **−0.59** | **0.00099** | **0.00087** |

The IFN/MHC inverse association **does not hold within genotype** (K n=9, r=−0.25, p=0.52; KP n=15, r=+0.01, p=0.96). Mean Cldn4 is higher in KP than K than T_early, and IFN/MHC is higher in T_early. The all-mouse signal is largely **genotype/stage composition**, not a within-model Cldn4 effect. No KL arm.

### 2. GSE6135 — KL vs other bulk tumors (n=25 tumors / 20 mice)

Ji et al. Affymetrix lung tumors. Cldn4 probe `1418283_at` present.

| Contrast | n | Result |
|---|---|---|
| Cldn4 vs T/NK (tumors) | 25 | r=−0.36, p=0.081 (trend, not significant) |
| Cldn4 vs IFN/MHC (tumors) | 25 | r=−0.30, p=0.15 |
| Cldn4 vs T/NK (mice) | 20 | r=−0.36, p=0.11 |
| Cldn4 vs IFN/MHC (mice) | 20 | r=−0.16, p=0.51 |
| KL-only Cldn4 vs T/NK | 13 tumors | r=−0.35, p=0.24; high-vs-low MWU p=0.10 |
| KL vs other mean Cldn4 | 13 vs 12 tumors | 6.18 vs 5.64 (array log2), MWU p=0.29 |

Direction is consistently **higher Cldn4, lower T/NK**, including inside KL, but **n is small and no test is significant**. This is the only scored public KL bulk with a processed matrix.

### 3. E-MTAB-5311 — LLC tumors (n=16 mice)

Expression Atlas raw counts. Cldn4 is **zero in 9/16** tumors (max 0.51 CPM). Split is therefore detected vs undetected, not a real high/low continuum.

| Endpoint | n | Spearman r | p | Detected (7) vs zero (9) MWU p |
|---|---:|---:|---:|---:|
| T/NK | 16 | +0.39 | 0.13 | 0.30 |
| IFN/MHC | 16 | +0.34 | 0.19 | 0.30 |

No support for a Cldn4–immune association in LLC at this n. Cldn4 abundance is too low for a strong high-vs-low claim.

## Bottom line

Public mouse data with **Cldn4 + immune/GEMM labels + processed matrices + honest n>2** are scarce once assigned GSEs are set aside.

- **KL:** GSE6135 is usable (n=13 KL tumors). Cldn4 vs T/NK is a non-significant negative trend. GSE165641/GSE180963 are KL scRNA but n=2.
- **K/KP scRNA:** GSE154989 gives n=28 and a strong Cldn4–IFN/MHC inverse **across** genotypes that **collapses inside** K or KP.
- **LLC / TISMO-class:** E-MTAB-5311 n=16, Cldn4 mostly undetected; no significant association.

No result here should be over-read as a KL Cldn4–T/NK law. The honest KL n in this hunt is **13 tumors / ~10 mice** (GSE6135), not 8 private KL.

## Re-run

Scripts in `methods/mouse_kl_cldn4_hunt/scripts/` expect processed files in `/tmp/geo_dl` (GEO/Atlas URLs listed in the catalog). Do not fetch SRA FASTQ.
