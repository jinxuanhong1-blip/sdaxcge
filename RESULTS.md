# Additive Visium LUAD histology — CLDN4-only

CLDN4 (single gene, not a signature) versus CD8A on open histology-stratified LUAD Visium. No GSE307534, GSE277206, or Zenodo 13337961. No private 8-KL.

## Data search

### 1. Lepidic vs acinar Visium (Clin Transl Med / PMC10844893)

Paper: Wang et al., *Clin Transl Med* 2024 (PMC10844893). Three FFPE Visium sections of stage IA LUAD with lepidic and acinar regions.

**Accession is not GEO or Zenodo.** Data availability points only to GSA-Human **HRA005794** (NGDC), listed as **controlled access** (DAC: Yang Dawei, Zhongshan Hospital). No processed matrices, GEO series, or Zenodo record was found for these Visium slides.

Open substitute with **pathologist lepidic and acinar spot labels**: **GSE273378** (below).

### 2. Open Visium of STAS lung cancer

| Source | Platform | STAS labels | Access |
|---|---|---|---|
| **GSE273378** (Taube et al., vascular invasion / stage I LUAD) | Visium FFPE, 16 sections | **Yes** — 57 STAS spots on 5 sections; also Lepidic / Acinar / Papillary / Solid / Micropapillary | Open GEO processed matrices + `pathology.csv` |
| J Transl Med 2024 S100P+TFF1+ STAS Visium (PMC11462816) | Visium FFPE | STAS vs non-STAS | **Request-only** (“available upon request”) |
| Frontiers 2025 STAS DSP (HRA012209) | GeoMx DSP, not Visium | STAS vs NSTAS | GSA-Human controlled |
| 10x Visium HD Human Lung Cancer FFPE (post-Xenium Prime 5K Exp. 2) | Visium HD 16 µm bins | STAS focus described in Long et al. *J Exp Clin Cancer Res* 2025; no public spot-level STAS mask | Open (CC BY 4.0) |

**Runnable open STAS Visium with labels: GSE273378.** HD LUAD was run as the public STAS-documented slide without a STAS mask.

### 3. Other histology-stratified LUAD Visium (not the excluded accessions)

| Dataset | Histology | n sections | Notes |
|---|---|---|---|
| **GSE273378** | Lepidic, acinar, papillary, solid, micropapillary, STAS | 16 | Spot-level pathology |
| **GSE300676** | Mixed lepidic / filigree mPAP / overt mPAP in each section (4 cases × 2) | 8 | No spot-level histology file; computational lepidic vs invasive scores |
| **GSE189487** | AIS, MIA, IAC (section-level) | 6 | AIS ≈ lepidic growth; IAC typically invasive patterns |
| **KERO Ad-SpatialAnalysis 2024** (Takano et al. *Nat Commun*) | AIS/MIA (TSU) vs IA (LUAD_No) | 4 downloaded (TSU-23, TSU-36, LUAD_No_17, FFPE_LUAD_No3_A) | Open processed Space Ranger; full set at https://kero.hgc.jp/Ad-SpatialAnalysis_2024.html |

HRA001238 (Cell Discovery ST_LUAD, 5 patients with lepidic/acinar/papillary/MP/solid annotations) is GSA-Human controlled.

## Methods (CLDN4-only)

Per section, after QC (UMI ≥100 and genes ≥50):

1. Log-normalize (log1p of 10⁴-scaled counts).
2. **Spearman(CLDN4, CD8A)** on all QC spots and on epithelial spots.
3. Epithelial spots: pathologist tumor labels when present (Lepidic, Acinar, Papillary, Solid, Micropapillary, Cribriform, STAS); otherwise top 40% of (KRT8+EPCAM+KRT19)/3.
4. Among epithelial spots, CLDN4 **Q4 vs Q1**:
   - **Nearest CD8A-high distance** (µm). CD8A-high = spots at/above the 75th percentile of CD8A. Array pitch 100 µm (standard Visium) or 16 µm (HD).
   - **Neighbor CD8A**: mean CD8A in the first hexagonal ring (~1.65× pitch). CD8A is sparse; first-ring means are often zero, so nearest-distance is the primary spatial metric.
   - **KRT8 residual**: residual of neighbor CD8A after OLS on KRT8 (epithelial-content control), plus Spearman of CLDN4 residual after KRT8 vs CD8A.
5. Stratify the same Q4/Q1 tests inside lepidic vs acinar labels, and STAS vs non-STAS tumor, when those labels exist.

Script: `scripts/run_cldn4_visium.py`. Tables: `results/section_metrics.csv`, `results/pooled_contrasts.csv`. Maps: `maps/`.

## Results

### Section-wide CLDN4 vs CD8A Spearman

| Dataset | n sections | median ρ | sections with ρ < 0 |
|---|---:|---:|---:|
| GSE273378 | 16 | +0.013 | 7 / 16 |
| GSE300676 | 8 | −0.014 | 5 / 8 |
| GSE189487 | 6 | −0.018 | 4 / 6 |
| KERO | 4 | +0.005 | 2 / 4 |
| Visium HD 16 µm | 1 | −0.008 | 1 / 1 |

No dataset shows a consistent section-wide anti-correlation. Sign tests against 50% negative are all p > 0.3.

### Q4 vs Q1 CLDN4 among epithelial spots — nearest CD8A-high

**GSE300676** (mixed lepidic–filigree–mPAP; strongest spatial signal):

| Section | ρ (all spots) | Q4−Q1 nn CD8A-high (µm) | p |
|---|---:|---:|---:|
| mPAP3_A | −0.034 | −72 | 0.026 |
| mPAP3_B | −0.014 | +132 | 1.3×10⁻¹⁶ |
| mPAP4_A | +0.010 | +27 | 0.031 |
| mPAP4_B | −0.014 | 0 | 0.35 |
| mPAP1_A | +0.021 | +65 | 0.071 |
| mPAP1_B | +0.005 | +91 | 6.6×10⁻⁹ |
| mPAP2_A | −0.043 | +229 | 3.3×10⁻⁹ |
| mPAP2_B | −0.052 | +98 | 4.0×10⁻⁴ |

In 6/8 sections, CLDN4-high (Q4) epithelial spots sit farther from CD8A-high spots than Q1 (positive delta). Four of those are p < 0.05. This is a local exclusion geometry, not a strong gene–gene Spearman.

**Visium HD LUAD (STAS-documented slide, 16 µm bins):** ρ = −0.008; Q4−Q1 nearest CD8A-high = +0.77 µm (p = 0.001). Direction matches exclusion; the absolute shift is small at 16 µm resolution because n is large (~10⁵ bins).

### Lepidic vs acinar (pathologist labels, GSE273378)

Five sections contain both lepidic and acinar spots: GSM8427428, GSM8427432, GSM8427436, GSM8427439, GSM8427442.

| Metric (epithelial, within pattern) | Lepidic median | Acinar median | Wilcoxon p (paired n=5) |
|---|---:|---:|---:|
| CLDN4–CD8A Spearman | −0.026 | −0.074 | 1.00 |
| Q4−Q1 nearest CD8A-high (µm) | ~0 | +41 | 0.81 |
| Mean CLDN4 | 1.48 | 1.31 | 0.81 |
| Mean CD8A | 0.150 | 0.112 | 0.31 |

Acinar point estimates are more CD8-cold (more negative Spearman, Q4 farther from CD8A-high, lower mean CD8A), but **n = 5 paired sections does not support a claim of difference**. Across all labeled spots (not only paired sections): 10 lepidic subsets (5/10 Spearman < 0) and 7 acinar subsets (5/7 Spearman < 0).

Lepidic-predominant LMP section GSM8427436 is an outlier in the opposite direction (epithelial Spearman +0.13).

### STAS vs non-STAS (GSE273378)

STAS spots are rare (57 total). Only GSM8427443 had enough STAS spots (n=35) for a within-STAS Q4/Q1 test (Spearman −0.024; Q4−Q1 nn +91 µm, p = 0.47).

Section-level (epithelial metrics; 5 STAS-containing vs 11 without):

| Metric | STAS+ median | STAS− median | MWU p |
|---|---:|---:|---:|
| CLDN4–CD8A Spearman | −0.013 | −0.005 | 0.58 |
| Q4−Q1 nearest CD8A-high | ~0 | ~0 | 0.83 |
| Mean CD8A | 0.047 | 0.100 | 0.58 |

No STAS vs non-STAS difference at this sample size. STAS+ sections trend toward lower mean CD8A.

### AIS vs IAC (GSE189487, section-level histology)

AIS: TD5 ρ = −0.097 (p = 6.6×10⁻⁵), Q4−Q1 nn +85 µm (p = 0.051); TD8 ρ = +0.077 (p = 0.001).  
IAC: TD1 ρ ≈ 0; TD2 ρ = −0.010.  
MIA in between and non-significant. Two AIS sections disagree in sign.

### KERO AIS/MIA vs IA

TSU-23 (AIS/MIA) section-wide ρ = +0.092; LUAD_No_17 (IA) epithelial Spearman −0.076 with Q4 farther from CD8A-high (p = 0.038). FFPE_LUAD_No3_A ρ = −0.048. Four sections only; treated as additive, not a histology test with labels.

### KRT8 residual

Q4 CLDN4 epithelial spots are also KRT8-higher in almost every section (epithelial-content confounder). After residualizing neighbor CD8A on KRT8, Q4−Q1 deltas collapse toward zero. Spearman of KRT8-residualized CLDN4 vs CD8A tracks the unadjusted epithelial Spearman and does not create a new exclusion signal.

## Maps

- Per-section four-panels (CLDN4, CD8A, epithelial CLDN4 quartiles, pathology or lepidic−invasive score): `maps/GSE273378__*.png`, `maps/GSE300676__*.png`, `maps/GSE189487__*.png`, `maps/KERO_AdSpatial2024__*.png`, `maps/VisiumHD_10x_LUAD_STAS__*.png`
- Pooled: `maps/summary_lepidic_acinar_stas.png`, `maps/summary_spearman_all_sections.png`

Key example with lepidic + acinar + STAS on one slide: `maps/GSE273378__GSM8427428_LM_SD_1216_1.png`. Highest STAS count: `maps/GSE273378__GSM8427443_LM_SD_15.png`.

## What this run adds

Open, labeled lepidic vs acinar and STAS vs non-STAS Visium (GSE273378) plus three other histology-stratified LUAD Visium resources. CLDN4-only vs CD8A is **not a uniform section-wide anti-correlation**. The clearest spatial pattern is **greater distance from CLDN4-high epithelial spots to CD8A-high spots** in GSE300676 mixed lepidic/mPAP tumors. Lepidic vs acinar and STAS vs non-STAS contrasts are directionally consistent with a colder acinar/STAS niche and are **underpowered** with the labeled n available here.
