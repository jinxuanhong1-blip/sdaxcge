# Results — public lung ICI scRNA+TCR

Computed from public GEO files. Unit = sample/patient. No fabricated statistics.
E2 (tumor-epithelial TACSTD2/CLDN4 vs expanded CXCL13+) is **not estimable**
on any cohort here.

---

## 1. GSE243013 (primary)

434,458 TCR-annotated T cells; **231** samples with TCR (paper n=234 patients;
GEO has more sampleIDs). All 231 have ≥50 TCR cells (min 59, median 1,704).
MPR-any = GEO `pCR`+`MPR`: **n=126**. `non-MPR`: **n=104**. One sample
(`P433`, label `unknowm`) excluded from tests.

CXCL13+ = author clusters **`CD4T_Tfh_CXCL13` + `CD4T_Th1-like_CXCL13`**
(21,373 + 13,338 cells). There is **no `CD8T_Tex_CXCL13`** in the public
TCR table. CD8 Tex (`HAVCR2` + `LAYN`) is secondary.

### 1.1 E1 vs MPR (Wilcoxon, two-sided)

| Endpoint | Median MPR-any (n=126) | Median non-MPR (n=104) | p | Cliff’s δ |
| --- | ---: | ---: | ---: | ---: |
| Expanded CXCL13+ clones / 1k CD8 | 16.0 | 28.3 | **7.64×10⁻⁶** | −0.343 |
| CXCL13+ fraction of all TCR cells | 0.0657 | 0.0828 | 4.63×10⁻⁴ | −0.268 |
| CXCL13+ among **expanded** cells | 0.0291 | 0.0560 | **2.85×10⁻⁶** | −0.359 |
| CXCL13+ among **non-expanded** cells | 0.0905 | 0.0981 | 0.0735 | −0.137 |
| Expanded-cell fraction (any clone size ≥2) | 0.469 | 0.477 | **0.977** | +0.002 |
| Expanded CD8 Tex clones / 1k CD8 | 25.7 | 34.7 | 1.01×10⁻⁵ | −0.338 |

Same direction in LUAD (n=21 vs 40) and LUSC (n=105 vs 64) for the primary
endpoint (p=0.00442 and 0.00161). pCR vs non-pCR is weaker
(16.6 vs 21.1 / 1k CD8, p=0.0413).

**Overall clone expansion does not differ by MPR.** What differs is the
**CXCL13+ (CD4 Tfh/Th1-like) share of expanded cells**, which is **lower**
in MPR-any at resection. This is a post-treatment immune atlas: MPR leaves
less residual viable tumor by definition. It is **not** evidence that CXCL13+
clones cause non-response, and it is **not** a CD8 Tex-CXCL13 result.

### 1.2 Combinatorial: expanded vs non-expanded (paired, n=231)

| Comparison | Median among expanded | Median among non-expanded | paired Wilcoxon p |
| --- | ---: | ---: | ---: |
| CXCL13+ cell fraction | 0.0398 | 0.0946 | **1.93×10⁻²⁴** |
| CD8 Tex cell fraction | 0.128 | 0.0205 | **5.64×10⁻³⁸** |

CXCL13+ (CD4 clusters) is **enriched in non-expanded** cells. CD8 Tex is
**enriched in expanded** cells. These are opposite combinatorial patterns
and must not be pooled.

### 1.3 Residual TACSTD2 / CLDN4 (not E2)

Immune-compartment detection from the public MTX. Spearman, n=231 depth-pass.

| Residual metric | TCR metric | ρ | p |
| --- | --- | ---: | ---: |
| TACSTD2 frac pos | expanded CXCL13+ / 1k CD8 | −0.040 | **0.541** |
| CLDN4 frac pos | expanded CXCL13+ / 1k CD8 | +0.044 | **0.502** |
| TACSTD2 frac pos | expanded-cell fraction | −0.003 | 0.970 |
| CLDN4 frac pos | expanded CD8 Tex / 1k CD8 | +0.212 | 0.00120 |
| TACSTD2 frac pos | expanded CD8 Tex / 1k CD8 | +0.140 | 0.0331 |

**Residual TACSTD2/CLDN4 does not track expanded CXCL13+.** The weak
positive residual-CLDN4 vs CD8 Tex correlation is compatible with shared
non-MPR / residual-tumor confounding (both residual CLDN4 and Tex are higher
in non-MPR). It is **not** a tumor-epithelial TACSTD2 program.

E2 status: `not estimable`.

Figures: `fig_gse243013_mpr_boxplots.png`,
`fig_gse243013_residual_tacstd2_vs_cxcl13.png`.

---

## 2. Caushi — GSE176022 processed exists; EGA skipped

**GSE176022 GEO processed exists** (`GSE176022_RAW.tar`, 10.9 MB): MiXCR
tables from MANAFEST / influenza / HIV / CEF **cultures**. Not scRNA, not
CXCL13, not TACSTD2. Used only as an open TCR leftover
(`gse176022_bulk_tcr_summary.tsv`).

**GSE176021** public processed used: 110 per-sample `.vdj.tar.gz` (0 failed)
+ CD3 annotation RDS. **493,408** paired TRA+TRB cells. EGA FASTQ and the
4.1 GB GEX `RAW.tar` were **not** downloaded. CXCL13 RNA is **not** in the
RDS (columns: `nCount_RNA`, `nFeature_RNA`, `barcode`, `imid`, `CellType`,
UMAP).

Tumor captures collapsed to **15 patients** (6 MPR, 9 non-MPR).
NY016-016 is normal-only in the series map and is not in the tumor table.

| Endpoint | Median MPR (n=6) | Median non-MPR (n=9) | p | Cliff’s δ |
| --- | ---: | ---: | ---: | ---: |
| Mean expanded-cell fraction (tumor) | 0.605 | 0.592 | **0.272** | +0.370 |
| Tfh-proxy fraction (RDS CellType) | 0.109 | 0.126 | 0.388 | −0.296 |

n=6 vs 9 is underpowered. Do not interpret the point estimate.

Combinatorial Tfh-proxy (paired, n=15): median 0.091 among expanded vs
0.164 among non-expanded, p=1.83×10⁻⁴ — same direction as GSE243013 CD4
CXCL13 clusters (Tfh-like cells are less expanded).

E2: not estimable (CD3-sorted T cells).

---

## 3. Leftover open TCR

### GSE179994

86,402 scTCR cells; 47 samples; 36 patients; barcode+sample join to T-cell
metadata **80.0%**. Cluster `Tex` used as a CXCL13-associated author label
(paper text; CXCL13 RNA not in the public T-cell table).

| Test | n | p | note |
| --- | ---: | ---: | --- |
| Tex fraction: expanded vs non-expanded (paired) | 47 | **1.30×10⁻⁸** | Tex enriched in expanded |
| Expanded-cell fraction: post vs pre | 14 vs 33 | 0.569 | null |
| Tex fraction: post vs pre (joined) | 14 vs 33 | 0.436 | null |

**Response / MPR is not on the GEO series matrix** (patient id + pre/on-treatment
only). Vs-response is **not estimable** from the public tables. Paper
Supplementary Table 1 was not redistributed.

### GSE185204

n=3 patients, public MTX+TCR ~900 MB. **Listed only; not downloaded.**
Descriptive n is too small for E1.

### Not lung ICI

GSE236581 (CRC ICB VDJ) and GSE162025 (NPC) excluded.

---

## 4. What is supported / not supported

**Supported (GSE243013, honest n/p):** at post-chemo-IO resection, MPR-any
samples have **fewer expanded CXCL13+ (CD4 Tfh/Th1-like) clones per CD8**
than non-MPR (n=126 vs 104, p=7.64×10⁻⁶). Overall TCR expansion fraction
is unchanged (p=0.977). Residual immune-compartment TACSTD2/CLDN4 does
**not** correlate with that CXCL13+ endpoint (p=0.54 / 0.50).

**Not supported / not estimable:**

- Tumor-epithelial TACSTD2/CLDN4 vs expanded CXCL13+ (E2) — no paired tumor score.
- A public `CD8T_Tex_CXCL13` cluster on GSE243013 — absent.
- Caushi CXCL13 RNA — not in the open RDS; EGA skipped.
- GSE176022 as scRNA+TCR — it is bulk culture TCR.
- GSE179994 vs MPR/RECIST — labels not on GEO.
- “Expanded clones are depleted in MPR” — false; only the CXCL13+ subset moves.

Tables: `GSE243013/`, `Caushi/`, `leftover/`.
