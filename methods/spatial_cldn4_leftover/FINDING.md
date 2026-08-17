# FINDING — leftover public spatial CLDN4 vs T/B neighborhood

Additive slice only. **Do not re-score** already-reported [GSE265899](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE265899), [GSE289483](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289483), [GSE273378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273378), or prior B6 series (GSE221322, E-MTAB-13530, GSE189487, GSE271689, GSE287472).

Question: in **other** public Visium / GeoMx / CosMx / Xenium lung series that actually measure **CLDN4**, how does CLDN4 sit relative to a T/B neighborhood?

Every ρ below is computed. n and p are reported. Series that lacked CLDN4 after a real matrix check are **empty for CLDN4** (not imputed).

Hunt log: `results/spatial_cldn4_leftover/tables/hunt_catalog.tsv`.

---

## Verdict

Leftover CLDN4 spatial is **not empty**. Three public leftover matrices have CLDN4 plus a T/B score:

| Series | Platform | Unit | CLDN4 vs T+B | n | p |
|---|---|---|---|---:|---:|
| [GSE292299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292299) NSCLC Visium | same-spot | section median ρ | **−0.027** (9/16 neg) | 16 | **0.50** |
| GSE292299 | hex ring 1 | section median ρ | **−0.101** (11/16 neg) | 16 | **0.083** |
| GSE292299 | ring 1, partial vs broad epi | section median ρ | **−0.042** (12/16 neg) | 16 | **0.039** |
| [GSE318867](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318867) SCLC Visium | same-spot | section median ρ | **−0.049** (3/3 neg) | 3 | **0.25** |
| GSE318867 | hex ring 1 | section median ρ | **−0.064** (3/3 neg) | 3 | **0.25** |
| [GSE299786](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299786) LUAD CosMx TMA1 | FOV mean | FOV ρ | **+0.422** | 36 | **0.010** |
| GSE299786 TMA2 | FOV mean | FOV ρ | **+0.552** | 38 | **3.3×10⁻⁴** |

**What holds as leftover:** NSCLC Visium GSE292299 is a real n=16 leftover. Same-spot CLDN4 vs T+B is **null** at the section level. The hex-ring neighborhood is a **small negative** that stays weakly negative after partialling a broad epithelial score (ring-1 partial p=0.039). CosMx LUAD FOVs go the **other way** (positive). SCLC Visium n=3 cannot support a signed-rank claim (two-sided Wilcoxon floor at n=3 is p=0.25).

**What is empty:** the best leftover GeoMx ICI series ([GSE329813](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329813), neoadjuvant chemo-IO, 22 patients, MPR labels) has TACSTD2 and a full T/B panel but **no CLDN4** in the public normalized matrix (1,812 genes). Xenium leftovers GSE319755 (568-plex) and GSE311609 (5k) also lack CLDN4.

Do not pool Visium ρ with CosMx FOV ρ.

---

## Methods (locked before leftover ρ)

| Score | Genes used if present |
|---|---|
| Target | `CLDN4` (required); `TACSTD2` reported as extra |
| Broad epithelial (control) | `EPCAM`, `KRT8`, `KRT18`, `KRT19`, `CDH1`, `KRT7` |
| T | `CD3D`, `CD3E`, `CD3G`, `CD2`, `CD8A`, `CD8B`, `TRAC`, `CD247`, `IL7R` |
| B | `MS4A1`, `CD79A`, `CD79B`, `CD19`, `BANK1`, `CD22` |
| T+B | mean of available T and B members |

- **Visium:** in-tissue spots, ≥200 genes, log1p(CP10K). Same-spot Spearman. Neighborhood = hex rings 1/2/3 from `(array_row, array_col)` neighbors `(row, col±2)` and `(row±1, col±1)`; require ≥3 neighbors. Partial Spearman = rank residual controlling for the broad epithelial score. **Section is the unit.** Report median ρ, n negative, two-sided Wilcoxon signed-rank p on the per-section ρ values.
- **CosMx:** public exprMat + `CenterX/Y_global_px`. QC: `cell_ID>0`, nCount≥20. Same-cell Spearman is exploratory (pseudoreplication). Neighborhood units: FOV mean, and within-FOV kNN-20 in pixel space (no invented micron scale).
- **GeoMx / Xenium:** gene check only unless CLDN4 is present.

Reproduce: `python3 methods/spatial_cldn4_leftover/analyze.py`

---

## GSE292299 — leftover NSCLC Visium (immunotherapy spatial biology)

Public filtered H5 + `tissue_positions` for **16 NSCLC** sections (DLBCL members not used). CLDN4 and TACSTD2 present. 121,087 QC spots.

### Section-level (primary)

| Test | n sections | median ρ | n neg / n pos | Wilcoxon p |
|---|---:|---:|---:|---:|
| same-spot CLDN4 vs T+B | 16 | −0.027 | 9 / 7 | 0.50 |
| same-spot CLDN4 vs T | 16 | −0.014 | 9 / 7 | 0.60 |
| same-spot CLDN4 vs B | 16 | +0.029 | 6 / 10 | 0.27 |
| same-spot TACSTD2 vs T+B | 16 | −0.024 | 9 / 7 | 0.74 |
| same-spot CLDN4 vs T+B \| broad epi | 16 | −0.028 | 10 / 6 | 0.38 |
| ring 1 CLDN4 vs neighbor T+B | 16 | −0.101 | 11 / 5 | 0.083 |
| ring 1 partial \| broad epi | 16 | −0.042 | 12 / 4 | **0.039** |
| ring 2 CLDN4 vs neighbor T+B | 16 | −0.092 | 11 / 5 | 0.18 |
| ring 3 CLDN4 vs neighbor T+B | 16 | −0.084 | 11 / 5 | 0.30 |

Same-spot sign is mixed. Strong negatives: P13 ρ=−0.308 (n=5,242, p=8.1×10⁻¹¹⁶); P11 ρ=−0.260 (n=13,202); P9 ρ=−0.181 (n=12,242). Strong positives: P12 ρ=+0.159 (n=2,042); P7 ρ=+0.121 (n=5,701). Ring-1 follows the same split (P13 ring-1 ρ=−0.449; P12 ring-1 ρ=+0.166). This leftover does **not** recover a uniform same-spot exclusion.

---

## GSE318867 — leftover SCLC Visium (STING paper spatial member)

Three sections with matched positions (VA↔VA2, VB↔VB1, VISA↔VA1). CLDN4 and TACSTD2 present. 9,610 QC spots.

| Section | n spots | same-spot CLDN4 vs T+B | ring-1 CLDN4 vs nei T+B |
|---|---:|---|---|
| VA | 2,736 | ρ=−0.049, p=0.011 | ρ=−0.065, p=7.2×10⁻⁴ |
| VB | 2,586 | ρ=−0.040, p=0.042 | ρ=−0.023, p=0.25 |
| VISA | 4,288 | ρ=−0.052, p=7.1×10⁻⁴ | ρ=−0.064, p=3.0×10⁻⁵ |

Section-level Wilcoxon n=3: p=0.25 for every all-negative contrast. **Honest n: too small to call.** TACSTD2 same-spot is weakly **positive** (3/3; median ρ=+0.054), the opposite of CLDN4.

---

## GSE299786 — leftover CosMx LUAD TMA

LUAD TMA1/TMA2 only (mesothelioma TMAs not used). 1,000-plex; CLDN4 and TACSTD2 present. Public metadata has `CenterX/Y_global_px`.

| TMA | cells (QC) | FOVs | CLDN4 %pos | same-cell ρ(CLDN4, T+B) | FOV ρ(CLDN4, T+B) | FOV partial \| epi | kNN-20 ρ |
|---|---:|---:|---:|---|---|---|---|
| TMA1 | 37,419 | 36 | 20.6% | +0.080 (p=6.0×10⁻⁵⁴) | **+0.422 (p=0.010)** | +0.202 (p=0.24) | +0.050 |
| TMA2 | 39,123 | 38 | 34.8% | +0.156 (p=6.2×10⁻²¹³) | **+0.552 (p=3.3×10⁻⁴)** | +0.395 (p=0.014) | +0.129 |

FOV is the confirmatory unit. Both TMAs are **positive**. TMA1 FOV partial vs epithelium is NS. Cell-level p-values are not confirmatory.

---

## Empty leftover (CLDN4 absent after a real check)

| Accession | Why it looked leftover | CLDN4 |
|---|---|---|
| GSE329813 | GeoMx NSCLC neoadjuvant pembrolizumab + chemo; 126 ROIs / 22 patients; MPR vs NMPR; tumor bed + LN | **No** (TACSTD2 yes; T/B yes) |
| GSE319755 | Xenium lung tumor / invasive front / adjacent | **No** (568 genes; TACSTD2 yes) |
| GSE311609 | Xenium 5k archival NSCLC | **No** (CLDN1/5/7/18 yes; CLDN4 no) |
| GSE246011 LUAD-14C | Visium counts | genes yes; **no coordinates** |

GSE329813 is the leftover GeoMx ICI series that would have been the additive GeoMx n. It is empty for CLDN4.

---

## What this does not claim

- It does not re-analyze GSE265899 / GSE289483 / GSE273378.
- It does not claim a pooled leftover meta-ρ.
- It does not claim CosMx micron balls (pixel kNN only).
- It does not claim GSE329813 CLDN4 (gene missing).
- Cell-level CosMx p-values are not confirmatory.

---

## Files

- `methods/spatial_cldn4_leftover/analyze.py`
- `methods/spatial_cldn4_leftover/hunt_catalog.tsv`
- `results/spatial_cldn4_leftover/tables/summary.json`
- `results/spatial_cldn4_leftover/tables/gse292299_samespot.csv`, `gse292299_rings.csv`, `gse292299_summary.json`
- `results/spatial_cldn4_leftover/tables/gse318867_samespot.csv`, `gse318867_rings.csv`, `gse318867_summary.json`
- `results/spatial_cldn4_leftover/tables/gse299786_summary.json`
- `results/spatial_cldn4_leftover/tables/gse329813_gene_check.json`
