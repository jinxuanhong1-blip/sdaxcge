# PAPER FUNNEL — Human MPE/scRNA vs PPT middle slides

**Never fabricate.** HRA006761 was not downloadable. Numbers below are from public **GSE131907** only (Kim et al., Nat Commun 2020).

## PPT middle-slide claims → public verdict

| Middle-slide claim | Public verdict | Where |
|---|---|---|
| TACSTD2-high in malignant vs immune | **Supported** (compartment detection) | MPE carcinoma-like vs PE T/NK/B/myeloid; also primary tS and author malignant |
| TJ enrichment in TACSTD2-high | **Supported** (cell-level module / Q4 vs Q1) | MPE carcinoma-like, primary tS, author malignant; holds after dropping CLDN4 from the module |
| CLDN4 coexpression with TACSTD2 | **Supported** (co-detection; primary sample levels) | MPE co-detection OR≈5; primary tS sample ρ=0.77. MPE double-positive *levels* do not rise together |
| MPE TACSTD2-high → fewer T/NK (sample) | **Not supported / underpowered** | Only 3 PE samples with ≥20 carcinoma-like cells |

These three middle-slide molecular claims are what this funnel corroborates. Sample-level “cold” MPE neighborhood is **not** claimed here.

## HRA006761 (not analyzed)

| Field | Value |
|---|---|
| Accession | HRA006761 (GSA-Human) |
| BioProject | PRJCA023797 |
| DAC | HDAC002197 |
| Page status | Controlled access / Request Data |
| Analyzed | **No** — no public matrix |

The publication (Clin Transl Med 2024;14:e1649) states CLDN4 correlates with ELF3, EpCAM, and TACSTD2 in recurrent MPE. That is the **authors’ claim**. It was not recomputed here.

## GSE131907 setup

- Matrix: author raw UMI; normalization `log1p(UMI / full-library UMI × 10000)`.
- Positive = raw UMI > 0.
- Five PE samples, 20,304 cells. Author `Cell_subtype` “Malignant cells” in PE: **0**.
- MPE malignant-like gate: PE epithelial **and** EPCAM>0 **and** WT1=0 **and** CALB2=0 → **259** carcinoma-like cells (3 samples ≥20 cells: EFFUSION_06/11/12).
- Mesothelial-like (EPCAM=0 and WT1+ or CALB2+): 76 cells; TACSTD2 detected in 6.6%.
- Primary tumor states tS1/tS2/tS3: 6,352 cells (10 samples with ≥20 after dropping LUNG_T09).
- TJ module = mean of CUSTOM_TJ_CORE genes present in the matrix (25 genes; PATJ absent, INADL kept). Secondary score drops CLDN4.

## 1. TACSTD2-high malignant vs immune

| Compartment | n | TACSTD2 %pos | mean log1p |
|---|---:|---:|---:|
| MPE carcinoma-like | 259 | **83.8** | 1.65 |
| MPE T | 11,974 | 1.4 | 0.020 |
| MPE NK | 1,297 | 3.1 | 0.047 |
| MPE B | 3,285 | 1.3 | 0.015 |
| MPE myeloid | 3,063 | 5.6 | 0.053 |
| Primary tS | 6,352 | 88.2 | 1.67 |
| Primary tLung T | 18,587 | 3.2 | 0.051 |
| Author malignant (all sites) | 24,784 | 72.1 | 1.06 |
| All T | 79,676 | 1.6 | 0.025 |

Mann–Whitney on log1p expression: MPE carcinoma-like vs MPE T, p≈0 (n=259 vs 11,974). The informative contrast for the slide is the **%pos gap**, not the asymptotic p-value.

**Not claimed:** sample-level anti-correlation of carcinoma-like TACSTD2 with T/NK fraction inside MPE. With ≥20 carcinoma-like cells, n=3 (Spearman undefined). With any carcinoma-like cells, n=4, TACSTD2 %pos vs T/NK ρ=+0.40, exact p=0.75.

Primary tS TACSTD2 mean vs T/NK: ρ=+0.10, p=0.78, n=10. Author-malignant TACSTD2 %pos vs T/NK: ρ=−0.15, p=0.51, n=21. (Author-malignant **CLDN4** %pos vs T/NK is ρ=−0.52, p=0.015, n=21 — that is the locked GSE131907 direction for CLDN4, not for TACSTD2.)

## 2. TJ enrichment in TACSTD2-high

Within-compartment TACSTD2 Q4 vs Q1, mean TJ-core module:

| Compartment | n Q4 / Q1 | mean TJ high | mean TJ low | Δ | MW p |
|---|---:|---:|---:|---:|---:|
| MPE carcinoma-like | 65 / 65 | 0.442 | 0.279 | +0.162 | 9.5×10⁻⁸ |
| Primary tS | 1,588 / 1,588 | 0.365 | 0.243 | +0.123 | 1.6×10⁻¹⁶¹ |
| Author malignant | 6,196 / 6,917 | 0.317 | 0.193 | +0.124 | ≈0 |

Holding CLDN4 out of the module does not remove the MPE contrast (Δ=+0.160, p=4.7×10⁻⁸).

Cell-level Spearman TACSTD2 vs TJ-core: MPE ρ=0.418 (p=2.4×10⁻¹²); primary tS ρ=0.389; author malignant ρ=0.345. Without CLDN4: MPE ρ=0.438.

Top MPE Q4-vs-Q1 detection odds ratios among TJ genes include MAGI1 (OR 9.7), MAGI3 (8.8), TJP3 (7.6), CLDN7 (5.0), CLDN4 (3.5).

**Caveat:** primary tS *sample-mean* TACSTD2 vs TJ-core is ρ=0.31, p=0.38, n=10 — underpowered between patients. Author-malignant sample means: TACSTD2 vs TJ-core ρ=0.46, p=0.037, n=21; without CLDN4 ρ=0.34, p=0.13.

This is a **module / cell-level** TJ enrichment, not a full prerank GSEA on the 2.9 GB log2TPM matrix (that atlas-level GSEA was already reported in PR #230: CUSTOM_TJ_CORE NES=2.21).

## 3. CLDN4 coexpression with TACSTD2

Co-detection (Fisher on UMI>0), TACSTD2 as anchor:

| Compartment | CLDN4 % in TACSTD2+ | in TACSTD2− | OR | Fisher p | level ρ (all) | level ρ (double+) |
|---|---:|---:|---:|---:|---:|---:|
| MPE carcinoma-like | 88.0 | 59.5 | 5.00 | 3.7×10⁻⁵ | 0.047 (p=0.45) | −0.025 (n=191, p=0.74) |
| Primary tS | 91.7 | 59.8 | 7.45 | 1.9×10⁻¹⁰² | 0.383 | 0.402 |
| Author malignant | 93.9 | 62.9 | 9.11 | ≈0 | 0.282 | 0.289 |

CLDN4-as-anchor TACSTD2 co-detection in MPE matches PR #635 (88.4% vs 60.5%, OR 5.00).

Primary tS sample-mean TACSTD2 vs CLDN4: **ρ=0.770, p=0.0092, n=10**.

ELF3 co-detection with TACSTD2 in MPE carcinoma-like: 87.6% vs 45.2%, OR 8.52, p=1.1×10⁻⁸; level ρ=0.248.

## What this does / does not do for the PPT

**Usable on middle slides (public, recomputed):**

1. TACSTD2 is a malignant/epithelial marker vs immune in MPE and primary compartments of GSE131907.  
2. TACSTD2-high malignant-like cells score higher on a TJ-core module, including after removing CLDN4.  
3. CLDN4 is co-detected with TACSTD2 in MPE carcinoma-like cells; primary tumor-state sample means track each other.

**Do not put on the slide as a Public MPE result:**

- HRA006761 recurrent-MPE coefficients (controlled; authors’ paper only).  
- MPE sample-level TACSTD2 vs T/NK anti-correlation (n≤4).  
- MPE double-positive TACSTD2–CLDN4 *level* correlation (null).  
- Replacing the locked CosMx exclusion or concordant-4 CLDN4–T/NK results.

## Files

- `results/tables/compartment_detection.tsv`, `malignant_vs_immune_mw.tsv`
- `results/tables/tj_enrichment.tsv`, `cldn4_coexpression.tsv`
- `results/tables/sample_compartment.tsv`, `sample_spearman.tsv`, `mpe_carcinoma_like_detail.tsv`
- `results/figures/fig_tacstd2_malig_vs_immune.png`, `fig_tj_q4q1.png`, `fig_cldn4_codetect.png`, `fig_mpe_tacstd2_vs_tnk.png`
- `results/summary.json`
