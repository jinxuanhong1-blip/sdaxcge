# GSE271689 — GeoMx WTA tumour-segment TACSTD2/CLDN4 vs ICI OS

## TL;DR (honest)

**High tumour-segment TACSTD2 is nominally associated with worse OS in the Yale first-line ICI discovery set. It does not replicate in the Greek WTA tumour table. UQ is the same direction and not significant. This is not a confirmed ICI-OS biomarker.**

| Cohort | n (OS events) | TACSTD2 OS HR / +1 SD | p | vs prior “high TROP2 worse” |
|---|---:|---:|---:|---|
| Yale (GSE271689, 1L ICI) | 37 (17) | **1.71 [1.02–2.88]** | **0.043** | same direction, nominal |
| Yale LUAD only | 25 (12) | 1.67 [0.95–2.91] | 0.073 | same direction, NS |
| UQ CTA (mixed line) | 30 (15) | 1.37 [0.87–2.17] | 0.18 | same direction, NS |
| Greek WTA tumour ROI | 61 (33) | 0.92 [0.69–1.23] | 0.59 | **opposite / null** |

- Yale **median-split** TACSTD2 OS is weaker (HR_hi/lo = 2.16, log-rank p = 0.12). Do not quote the continuous p as a KM result.
- **CLDN4** tracks TACSTD2 in Yale (ρ = 0.57) but is NS for OS (HR = 1.38, p = 0.15).
- **CD8-ish signatures are not associated with TACSTD2 in Yale** (ρ ≈ 0). CLDN4 vs CD8A is inverse (ρ = −0.39, p = 0.016). UQ CD8A itself predicts better OS (HR = 0.52, p = 0.012); TACSTD2 does not.
- Greek TACSTD2 vs cytotoxic signature is modestly inverse (ρ = −0.30, p = 0.018) **without** an OS effect.

Do not cite GSE271689 as independent confirmation that TROP2-high NSCLC is ICI-resistant. Yale is small (17 deaths). Greece is the better-powered WTA OS set and is null.

## Dataset

Aung, Monkman, Warrell et al., *Nature Genetics* 2025 (PMID 41073787): spatial proteomics + GeoMx WTA in advanced NSCLC on PD-1–based immunotherapy.

- **GSE271689** = Yale discovery + Greek validation WTA DCC files (586 AOIs). GEO has **no gene matrix and no survival table**.
- Tumour segments are PanCK / CK AOIs. Immune segments (CD45, CD68) were built but are not the primary TROP2 readout.
- Outcomes for the 37 Yale patients used in the paper’s Fig. 6b KM are in Source Data Fig. 6b (public XLSX). All 37 are first-line ICI, pre-treatment biopsy (32 pembrolizumab).
- UQ is GeoMx **CTA** (GSE221733; ~1.8k genes). **CLDN4 is not on the panel.**
- Greek tumour WTA used here is Source Data Fig. 6d (61 ROIs; almost 1:1 with unique OS/PFS pairs). Values are already centred (can be negative).

## What was used

| Input | Why |
|---|---|
| `GSE271689_RAW.tar` DCC files | Public raw counts for Yale (+ Greek wells; survival is not in GEO) |
| `Hs_R_NGS_WTA_v1.0.pkc` (Zenodo 12752405) | RTS_ID → gene (WTA v1.0) |
| GEO series matrix | AOI → spot id / CK vs CD45 vs CD68 |
| Source Data Fig. 6b | Yale OS/PFS + clinical |
| Source Data Fig. 6c | UQ CTA tumour + 2-year-censored OS |
| Source Data Fig. 6d | Greek WTA tumour + OS |

Reprocessing check vs the paper’s `CK_av_3B.csv` (37 shared patients): TACSTD2 Spearman ρ = 0.96, CLDN4 ρ = 0.90.

## What was skipped (honest)

- SRA/FASTQ re-alignment. DCC files are the deposited count level.
- Training a new LASSO on ~18k genes. That would overfit n = 37.
- Treating UQ as WTA. It is CTA; CLDN4 cannot be tested.
- Calling Greek `File_number` a patient ID (one value has 22 ROIs and mixed histotypes). Rows are used as the paper used them (ROI ≈ patient).
- Multiple-testing correction across every gene. TACSTD2, CLDN4, and the CD8-ish scores were pre-specified.

## Methods

1. Map DCC RTS IDs with WTA v1.0 PKC; drop NTC / failed AOIs (raw reads, alignment, saturation, ≥10k counts, CK/CD45/CD68).
2. Q3-normalise, log2(x+1). Patient score = mean of tumour (CK) AOIs (2–4 cores; 37 patients with paper outcomes).
3. Pre-specified features: TACSTD2, CLDN4, CD8A, mean z-score of {CD8A, CD8B, GZMA, GZMB, PRF1, NKG7}, Ayers-like IFNG set, CXCL9.
4. Cox OS/PFS on z-scored features (primary). Median-split KM is secondary. Yale multivariable: TACSTD2 + CD8A + adenocarcinoma indicator.
5. Same Cox/correlation on UQ and Greek public tumour tables.

Reproduction: `python3 scripts/hunt_gse271689/03_analyze_trop2_os.py` after the fetch/build steps in `00_fetch_data.sh`.

## Results

### 1. Yale OS — TACSTD2 worse, barely

Continuous Cox, n = 37 / 17 deaths:

- TACSTD2 HR = 1.71 (1.02–2.88), p = 0.043, C-index = 0.63
- CLDN4 HR = 1.38 (0.89–2.13), p = 0.15
- CD8A HR = 0.92 (0.59–1.44), p = 0.72
- cytotoxic signature HR = 1.04, p = 0.87
- Ayers-like HR = 0.82, p = 0.42

Median-split TACSTD2 OS log-rank p = 0.12. PFS is the same direction and NS (HR = 1.48, p = 0.078).

Multivariable OS ~ TACSTD2_z + CD8A_z + adeno: TACSTD2 HR = 1.84 (1.04–3.26), p = 0.036. CD8A drops out. This is still 17 events / 3 coefficients.

### 2. Immune — not a cold-TROP2 story in Yale tumour AOIs

| pair | Spearman ρ | p |
|---|---:|---:|
| TACSTD2–CLDN4 | 0.57 | 2.1e-4 |
| TACSTD2–CD8A | −0.19 | 0.27 |
| TACSTD2–cytotoxic | −0.08 | 0.65 |
| TACSTD2–Ayers | 0.01 | 0.96 |
| CLDN4–CD8A | −0.39 | 0.016 |

CD8B / GZMA / PRF1 / IFNG sit near the GeoMx LOQ in CK AOIs. CD8A and CXCL9 are detected. The cytotoxic mean is therefore mostly CD8A plus noise. Do not over-interpret “no CD8 association” as biology vs assay floor.

### 3. UQ — CD8 works, TACSTD2 does not

CTA tumour, 2-year-censored OS, n = 30 / 15 deaths (10 first-line; first-line-only Cox not estimable, 3 deaths).

- TACSTD2 HR = 1.37, p = 0.18
- CD8A HR = 0.52 (0.31–0.87), p = 0.012
- TACSTD2 vs CD8A ρ = 0.08

CLDN4 is absent. This cohort can sanity-check a CD8-ish signal; it cannot test claudin-4.

### 4. Greek — null OS, weak inverse with cytotoxic score

WTA tumour, n = 61 / 33 deaths.

- TACSTD2 HR = 0.92, p = 0.59
- CLDN4 HR = 1.10, p = 0.61
- TACSTD2 vs CLDN4 ρ = 0.35
- TACSTD2 vs cytotoxic ρ = −0.30, p = 0.018

If the claim is “TROP2-high tumours are CD8-low,” Greece is the only WTA set that even hints at it, and it still has no OS effect.

## Files

| File | Contents |
|---|---|
| `stats.txt` | Numeric dump |
| `cox_os_pfs.csv` | All Cox rows (Yale / UQ / Greek, OS and PFS) |
| `cox_multivariable_yale.csv` | Yale OS ~ TACSTD2 + CD8A + adeno |
| `correlations.csv` | Pre-specified Spearman pairs |
| `yale_patient_compact.csv` | 37 patients, focus genes + OS/PFS |
| `uq_patient_compact.csv` | 30 UQ samples |
| `greek_roi_compact.csv` | 61 Greek tumour ROIs |
| `km_*.png` / `forest_os_continuous.png` / `scatter_yale_*.png` | Figures |
| `provenance.json` | URLs and what was not committed |

Large DCC / Q3 matrices / paper XLSX are **not** in git. Rebuild with `scripts/hunt_gse271689/00_fetch_data.sh`.
