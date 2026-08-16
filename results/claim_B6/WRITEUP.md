# Claim B6 — CLDN4/TACSTD2-high spots and T/B-cell spatial exclusion

**Claim (as given):** In public Visium / GeoMx (and CosMx if available) lung data, CLDN4/TACSTD2-high spots anti-colocalize with T and B cells.

**This folder:** `results/claim_B6/` only. Scripts: `scripts/run_*.py`, `scripts/common.py`, `scripts/neighborhood.py`.

**PR #110 is a separate radius rework** (`results/rework/B6_radius/`, branch `cursor/b6-radius-rework-bea8`). It used a different immune definition (broad immune / T-effector), dropped immune-rich index spots, and reported tumor-only hex-ring **|median partial ρ| ≤ 0.038** (claim not supported). This analysis does **not** replace or edit PR #110.

No statistic below was invented. Every number is in `tables/` and `summary.json`.

---

## Honest verdict

**PARTIAL — same-spot / compartment anti-colocalization is real in tumor Visium and GeoMx; single-cell CosMx neighborhoods are not.**

| Test | What it actually measures | Result | Supports “avoids immune niches”? |
|---|---|---|---|
| Visium same-spot Spearman (CLDN4+TACSTD2 vs T/B) | Do the two programs occupy the **same 55 µm spot**? | Tumor: 26/26 sections ρ < 0. Median ρ(epi, T+B) = **−0.255** (E-MTAB, n=20) and **−0.217** (GSE189487, n=6). Adjacent/healthy ≈ 0. | **Compositional yes.** Expected if CLDN4 is epithelial and CD3/MS4A1 are lymphoid. |
| Visium hex rings 1/2/3 (CLDN4 vs neighbor T+B) | Are neighbors of CLDN4-high spots lymphocyte-low? | Tumor ring 1 median ρ = **−0.333** (E-MTAB) / **−0.213** (GSE189487), all tumor sections negative. Shrinks at rings 2–3. Adjacent near 0. | **Yes, but confounded** by tumor-core vs stroma autocorrelation. |
| Same rings, partial ρ controlling for KRT/EPCAM | Extra CLDN4 signal after epithelial content | E-MTAB tumor ring 1 median partial ρ = **−0.185** (20/20 negative). | **Weaker, still negative** with a T/B neighborhood. |
| PR#110-like index (epi-high ∩ not T+B-rich) + partial | Closest comparison to PR #110, **T/B not broad-immune** | E-MTAB ring 1 median partial ρ = **−0.194** (19/20 negative). GSE189487 ring 1 median **−0.063**. Effect decays with ring. | **Supports T/B exclusion more than PR #110’s broad-immune test.** Different immune definition. |
| GeoMx GSE271689 | AOIs pre-cut into CK / CD45 / CD68 | All AOIs ρ(epi, T+B) = **−0.429**. Within CK still **−0.368**. CK mean epi +0.46 vs CD45 −0.35; paired ROI epi CK>CD45 and T+B CK<CD45. | **Compartment yes.** Antibody segmentation builds in tumor vs leukocyte separation. |
| CosMx GSE287472 (1 LUAD, 38 574 QC cells) | Same-cell and 25/50/100 µm balls | Neighborhood ρ ≈ **0.00** at every radius. Same-cell ρ(CLDN4, T+B) = **+0.046** (opposite of exclusion). CLDN4 detected in **3.6%** of cells. | **Does not support.** n=1 patient; sparse 1k-plex counts. |

**Bottom line.** The claim is true as **spot-level / compartment-level anti-colocalization of an epithelial program with lymphocyte programs in tumor Visium and GeoMx**. That is not the same as a CLDN4-specific immune-exclusion niche after architecture is removed. PR #110’s broad-immune rings stay near zero. This analysis’s T/B rings stay negative. CosMx (the only public CosMx LUAD matrix with CLDN4 + coordinates we found) shows no neighborhood exclusion.

---

## Data (processed files only)

| Accession | Platform | Tissue | What was used | Decision |
|---|---|---|---|---|
| [E-MTAB-13530](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13530) | 10x Visium | 8 NSCLC tumors (20 sections), 8 adjacent (16), 2 healthy (4) | `filtered_feature_bc_matrix.h5` + `spatial.tar` | Ran |
| [GSE189487](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189487) | 10x Visium | 6 early LUAD (AIS/MIA/IAC) | MTX + features + barcodes + `tissue_positions_list` | Ran |
| [GSE271689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271689) | GeoMx DSP WTA | Advanced NSCLC, CK / CD45 / CD68 AOIs | GEO DCC + [Hs_R_NGS_WTA_v1.0.pkc](https://zenodo.org/records/12752405) (probe→gene) + family.soft annotations | Ran |
| [GSE287472](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287472) | CosMx SMI (~1.2k genes) | 1 LUAD (GSM8746298) | `exprMat` + `metadata` (not the 144 MB transcript file) | Ran; n=1 |

FASTQ/SRA were not downloaded. Raw caches live under `data/` and are gitignored.

No other public CosMx **human lung cancer** matrix with CLDN4 + cell coordinates was usable in this run. Other GEO CosMx “lung” hits were mouse metastasis, IPF, infection, organoids, or non-lung.

---

## Methods (pre-specified)

**Signatures** (`scripts/common.py`)

- Epithelial (claim genes): `CLDN4`, `TACSTD2`
- Epithelial broad (control only): `EPCAM`, `KRT8`, `KRT18`, `KRT19`, `CDH1`, `KRT7`
- T: `CD3D`, `CD3E`, `CD3G`, `CD2`, `CD8A`, `CD8B`, `TRAC`, `CD247`, `IL7R`
- B: `MS4A1`, `CD79A`, `CD79B`, `CD19`, `BANK1`, `CD22`

Visium: in-tissue spots, ≥200 genes, CP10K + log1p, `scanpy.tl.score_genes`.

**Same-spot:** Spearman of epi score vs T / B / T+B across spots of one section. Also rank-biserial of T+B in the top epi tertile vs the rest. Cross-section inference: Wilcoxon signed-rank and two-sided sign test on the per-section ρ values.

**Neighborhood (Visium):** honeycomb rings. Neighbors of `(row, col)` are `(row, col±2)` and `(row±1, col±1)`. Ring *k* = graph distance *k*, self excluded, ≥3 neighbors. Primary: Spearman of spot **CLDN4** vs mean ring T+B. Sensitivity: partial Spearman controlling for the broad epithelial score.

**PR#110-like sensitivity (this repo only):** index = broad-epi ≥ section median **and** T+B < section Q3. Partial Spearman of index CLDN4 vs ring T+B, controlling for broad-epi. Immune-rich spots stay in the **neighborhood** average. This is **not** PR #110 (different immune list, no T-effector score, no published-matrix GeoMx).

**GeoMx:** sum probes to genes via WTA PKC; drop NTC / `cell type=NA`; total counts ≥ 1000; Q3 normalize; log1p; same signatures. AOIs are antibody compartments, not a hex grid.

**CosMx:** QC `nCount_RNA ≥ 20`, `Area.um2 ∈ [20, 2000]`; µm/px = median `sqrt(Area.um2 / Area)` = 0.120; radii 25 / 50 / 100 µm; ≥5 neighbors; log1p(CP10K); z-mean T/B scores. One patient.

Thresholds were not tuned after seeing results.

---

## Results

### 1. Visium same-spot

Figure: `figures/visium_samespot_rho.png`. Table: `tables/visium_samespot_all.csv`.

**E-MTAB-13530 tumor (20/20 sections, 8 patients)**

| Contrast | median ρ | n negative | Wilcoxon p |
|---|---:|---:|---:|
| epi vs T+B | −0.255 | 20/20 | 1.9×10⁻⁶ |
| epi vs T | −0.243 | 20/20 | 1.9×10⁻⁶ |
| epi vs B | −0.160 | 20/20 | 1.9×10⁻⁶ |
| rank-biserial T+B in epi-high tertile | −0.259 | 20/20 | 1.9×10⁻⁶ |

**GSE189487 (6/6; AIS TD5/TD8, MIA TD3/TD6, IAC TD1/TD2)**

| Contrast | median ρ | n negative | Wilcoxon p |
|---|---:|---:|---:|
| epi vs T+B | −0.217 | 6/6 | 0.031 |
| epi vs T | −0.201 | 6/6 | 0.031 |
| epi vs B | −0.106 | 5/6 | 0.062 |

E-MTAB adjacent (n=16) median ρ(epi, T+B) = **−0.016** (Wilcoxon p = 1.0). Healthy (n=4) = **−0.005**. The anti-colocalization is a **tumor-section** phenomenon, not a Visium artefact.

### 2. Visium hex-ring neighborhood

Figure: `figures/visium_neighborhood_rings.png`. Table: `tables/visium_neighborhood_rings.csv`.

**E-MTAB tumor, CLDN4 vs neighbor T+B (all QC spots)**

| Ring | median ρ | n negative | Wilcoxon p | median partial ρ (epi-controlled) |
|---|---:|---:|---:|---:|
| 1 | −0.333 | 20/20 | 1.9×10⁻⁶ | −0.185 |
| 2 | −0.269 | 20/20 | 1.9×10⁻⁶ | −0.144 |
| 3 | −0.217 | 20/20 | 1.9×10⁻⁶ | −0.097 |

**GSE189487 ring 1:** median ρ = −0.213 (6/6 negative); median partial = −0.108.

Adjacent-normal ring 1 is mixed and near zero for several patients (P10, P16, P17, P24); a few (P11, P15) are negative. Healthy ring 1 is **positive** in D1_1 (ρ = +0.15). The tumor-restricted pattern argues against a generic “high-expression spots anti-correlate with everything” artefact.

**Caveat.** Neighbors of a tumor-core spot are usually more tumor-core spots. Raw neighborhood ρ therefore repeats the same-spot composition test with spatial smoothing. Partialling the KRT/EPCAM score shrinks |ρ| but does not zero it when the neighborhood is T/B.

### 3. PR#110-like filter (comparison only)

Figure: `figures/visium_pr110like_partial.png`. Table: `tables/visium_pr110like_rings.csv`.

Index = epithelial-high ∩ not T+B-rich; partial Spearman of CLDN4 vs ring T+B.

| Dataset | Ring | median partial ρ | n negative | Wilcoxon p |
|---|---:|---:|---:|---:|
| E-MTAB tumor | 1 | −0.194 | 19/20 | 3.8×10⁻⁶ |
| E-MTAB tumor | 2 | −0.117 | 19/20 | 9.5×10⁻⁶ |
| E-MTAB tumor | 3 | −0.088 | 17/20 | 2.1×10⁻⁴ |
| GSE189487 | 1 | −0.063 | 6/6 | 0.031 |
| GSE189487 | 2 | −0.064 | 6/6 | 0.031 |
| GSE189487 | 3 | −0.082 | 6/6 | 0.031 |

P15_T2 is the one E-MTAB tumor section with ring-1 partial ρ ≈ 0 (+0.005). GSE189487 effects are small (|median| ≈ 0.06–0.08), closer to PR #110’s band.

**Why this is not a contradiction of PR #110.** PR #110’s primary neighborhood was **broad immune**, not T+B, and they also reported T-effector rings that were weakly **positive**. Swapping the immune definition changes the sign and the size. Both results can be true at once.

### 4. GeoMx GSE271689

Figures: `figures/geomx_compartment_scores.png`. Tables: `gse271689_geomx_*.csv`.

577 AOIs after QC (CK 211, CD45 174, CD68 192).

| Scope | n | ρ(epi, T+B) | ρ(epi, T) | ρ(epi, B) |
|---|---:|---:|---:|---:|
| all AOIs | 577 | −0.429 | −0.439 | −0.270 |
| within CK | 211 | −0.368 | −0.376 | −0.221 |
| within CD45 | 174 | −0.212 | −0.260 | −0.062 (n.s.) |
| within CD68 | 192 | +0.031 | +0.038 | −0.042 |

Compartment means (score_genes, log1p Q3):

| Compartment | epi | T | B | T+B |
|---|---:|---:|---:|---:|
| CK | +0.462 | −0.173 | −0.080 | −0.132 |
| CD68 | −0.294 | −0.007 | −0.008 | −0.003 |
| CD45 | −0.348 | +0.074 | +0.095 | +0.082 |

Paired ROIs with both CK and CD45 (n=59): epi CK−CD45 median = **+0.814** (Wilcoxon p = 8.4×10⁻¹⁰); T+B CK−CD45 median = **−0.179** (p = 1.9×10⁻⁸).

**Caveat.** Segmentation used PanCK / CD45 / CD68 **protein**, not CLDN4. Transcript CLDN4/TACSTD2 vs CD3/MS4A1 is semi-independent, but the AOI map is still tumor vs leukocyte by design. This is not a hex-ring test. PR #110 additionally analyzed published tumor-only Yale/Greek matrices and found Yale CLDN4–CD8A ρ = −0.39 that did **not** replicate in Greek (ρ = −0.05). We did not re-download those supplementary sheets.

### 5. CosMx GSE287472

Figure: `figures/cosmx_neighborhood_rho.png`. Tables: `gse287472_cosmx_*.csv`.

One LUAD, 64 420 cells → 38 574 after QC. Panel includes CLDN4, TACSTD2, CD3D/E/G, CD8A/B, MS4A1, CD79A. CLDN4 is nonzero in **3.6%** of cells (TACSTD2 1.4%). A top-tertile “CLDN4-high” contrast is undefined (2/3 quantile = 0); that is why `rb_TBnb_in_CLDN4high` is NaN.

| Radius | scope | n | ρ(CLDN4, neighbor T+B) | p |
|---|---|---:|---:|---:|
| same cell | all QC | 38 574 | +0.046 | 3.4×10⁻¹⁹ |
| same cell | PanCK-high | 19 300 | +0.041 | 1.0×10⁻⁸ |
| 25 µm | all QC | 26 450 with ≥5 neighbors | +0.003 | 0.66 |
| 50 µm | all QC | 32 541 | +0.007 | 0.21 |
| 100 µm | all QC | 37 303 | +0.003 | 0.51 |
| 25–100 µm | PanCK-high | — | −0.004 to +0.002 | all p > 0.5 |

Same-cell ρ is small and **positive** (CLDN4 and lymphocyte transcripts are not mutually exclusive in this 1k-plex count matrix). Neighborhood ρ is indistinguishable from zero. This does **not** support spatial avoidance. It is one patient and a sparse panel; it also does not strongly refute the Visium tumor result.

---

## What would over-claim

- Treating same-spot ρ ≈ −0.25 as proof that CLDN4 **causes** immune exclusion. It is a cell-type composition correlation.
- Ignoring PR #110. Their broad-immune rings are near zero on the same E-MTAB tumor slides.
- Treating GeoMx CK vs CD45 as an independent spatial neighborhood.
- Treating CosMx n=1, 3.6% CLDN4 detection as a definitive cohort test.
- Inventing a CosMx “CLDN4-high tertile” effect (the tertile is empty).

---

## Files

```
results/claim_B6/
  WRITEUP.md
  summary.json
  figures/
    visium_samespot_rho.{png,pdf}
    visium_neighborhood_rings.{png,pdf}
    visium_pr110like_partial.{png,pdf}
    geomx_compartment_scores.{png,pdf}
    cosmx_neighborhood_rho.{png,pdf}
  tables/
    emtab13530_per_section.csv
    gse189487_per_section.csv
    visium_samespot_all.csv
    visium_neighborhood_rings.csv
    visium_pr110like_rings.csv
    gse271689_geomx_*.csv
    gse287472_cosmx_*.csv
```
