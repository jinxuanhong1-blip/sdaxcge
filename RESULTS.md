# RESULTS — Additive public Visium, CLDN4-only (Molecular Cancer 2025)

Post-neoadjuvant chemoimmunotherapy NSCLC Visium from Cui / Liu *Mol Cancer* 2025 (doi:10.1186/s12943-025-02287-w). **CLDN4-only.** No private 8-KL. No claim-failed.

## Accessions

Paper Data Availability names **SRA**, not GEO or Zenodo. Spatial raw is public FASTQ. There is no processed Visium matrix and no H&E.

| Role | Accession | Status |
|---|---|---|
| This paper, spatial raw | [PRJNA1139087](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1139087) / [SRP521546](https://www.ncbi.nlm.nih.gov/sra?term=SRP521546) | **Public** FASTQ (ENA). Downloaded. |
| This paper, scRNA raw | PRJNA1068179 | Public (not used) |
| Paper validation scRNA | GSE207422 | Public scRNA (not Visium) |
| Open pre-ICI NSCLC Visium | [GSE292299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292299) | Public filtered H5 + `tissue_positions` (16 pre-treatment NSCLC). Raw withheld. |

Sample map (BioSample isolate):

| Patient | Alias | SRA | Pathology (paper) | Counted read pairs | Target-panel UMIs |
|---|---|---|---|---:|---:|
| PA08 | Tumor-1 | SRR29925401 | NMPR, “cold” | 369,145,160 | 647,939 |
| PA09 | Tumor-4 | SRR29925398 | NMPR, “cold” | 366,823,687 | 249,125 |
| PA10 | Tumor-2 | SRR29925400 | pCR | 275,367,577 | 35,372 |
| PA12 | Tumor-3 | SRR29925399 | NMPR, therapy-responsive histology | 315,790,119 / 398,010,594 SRA | 260,989 |

Paper ST: 12,741 spots after their QC. This count: **12,979** tissue-filtered barcodes (≈3,245 / section). All four sections are **post** neoadjuvant PD-L1 blockade + cisplatin. There is **no matched pre Visium** in this deposit.

PA12 FASTQ bytes matched ENA exactly; the pair iterator stopped at 79% of the SRA spot count. Target-gene UMI discovery had already plateaued (259,932 UMIs at 300M reads → 260,989 at 315.8M), so the locked CLDN4 / CD8A / KRT8 panel is not read-depth limited.

## Other ICI / neoadjuvant NSCLC Visium (open vs controlled)

| Dataset | Pre / post | Open? |
|---|---|---|
| PRJNA1139087 (this paper) | Post only (n=4) | Yes, raw FASTQ |
| GSE292299 Behera / Raman | Pre-ICI biopsies (16 NSCLC) | Yes, processed matrices |
| Yan *Nat Genet* 2024 | Pre GeoMx + post Visium | **Controlled** Zenodo 10.5281/zenodo.8227624 / 13901289; NGDC PRJCA009862 |
| JITC 2025 B-cell case (Zenodo 15053161) | Pre brain met + post lung | Record public, **files restricted** |
| NADIM Visium (SITC 2025 abstract) | Pre + post | No public accession found |
| E-MTAB-13530 De Zuani | Treatment-naive NSCLC | Open, **not ICI** |
| GSE189487 Zhu | AIS→IAC | Open, **not ICI** |
| GSE329813 | Neoadjuvant chemo-IO | GeoMx, **no CLDN4** in public matrix |

If a later deposit of Yan or NADIM is released without a click-through, rerun the same locked metrics. Do not treat restricted Zenodo as open.

## Question

Are CLDN4-high **tumor** spots farther from CD8 / lower CD8 neighbor? After ICI, does that exclusion still hold?

## Methods (locked)

- **PRJNA1139087:** stream official Visium Human Transcriptome Probe Set v1.0 sequences for CLDN4, CD8A/B, KRT8/18/19, EPCAM, CDH1, KRT7, CD3D/E/G, NKG7, PTPRC. R1 16-mer → Visium v1 whitelist (1-edit). R2 insert before polyA → probe. UMI collapse. Library size = whitelist reads / barcode (probe-targeted count is not whole-transcriptome). Tissue = barcodes above a reads quantile. Coordinates from `visium-v1_coordinates.txt` (no H&E).
- **GSE292299:** public filtered H5, in-tissue, ≥200 genes, log1p CP10k from full library.
- Tumor spots: epithelial score (`EPCAM,KRT8,KRT18,KRT19,CDH1,KRT7`) ≥ section median.
- CLDN4-high / low: top vs bottom quartile of CLDN4 among tumor spots.
- CD8+ spots: CD8A count > 0. Nearest distance in µm from hex array (100 µm center-to-center).
- Radius neighbors: hex rings 1/2/3 as `(row, col±2)` and `(row±1, col±1)`.
- KRT8 residual: rank residual / partial Spearman controlling KRT8.
- **Section is the unit.** Report n, median, sign count, Wilcoxon p. n=4 post cannot reject a two-sided Wilcoxon at 0.05 even if all signs agree (floor p=0.125).

Reproduce: `python3 methods/molcancer_visium_cldn4/count_fastq.py --aria2 && python3 methods/molcancer_visium_cldn4/analyze.py`

## Results

### GSE292299 — open pre-ICI NSCLC Visium (n=16 sections, 121,087 QC spots)

Not matched to PA08–12. CD8A is the Gene Expression feature (the H5 also has an empty antibody `CD8A` column; first ENSG hit is used).

| Locked test | n sections | median | n neg / n pos | Wilcoxon p |
|---|---:|---:|---:|---:|
| Tumor same-spot CLDN4 vs CD8A | 16 | ρ = **−0.0035** | 9 / 7 | 0.56 |
| Tumor same-spot CLDN4 vs CD8A \| KRT8 | 16 | ρ = **+0.0017** | 7 / 9 | 0.67 |
| Tumor CLDN4 vs nearest CD8A+ (µm) | 16 | ρ = **+0.0063** | 6 / 10 | 0.94 |
| CLDN4-high − low nearest CD8A (µm) | 16 | Δ = **0** | 4 / 5 (rest 0) | 0.95 |
| Tumor CLDN4 vs hex ring-1 CD8A | 16 | ρ = **−0.033** | 13 / 3 | **0.011** |
| Ring-1 partial \| KRT8 | 16 | ρ = **−0.042** | 12 / 4 | **0.021** |
| CLDN4-high − low ring-1 CD8A | 16 | Δ = **0** | 7 / 2 | 0.37 |

**Pre-ICI, computed:** same-spot and nearest-CD8 distance do **not** show CLDN4-high tumor spots farther from CD8. The only signed section-level signal is a **small negative** hex ring-1 neighbor CD8A that stays negative after KRT8 residual. Quartile high-vs-low deltas sit at 0. This is not a same-spot exclusion claim.

Maps: `results/molcancer_visium_cldn4/maps/GSE292299_NSCLC_P1–P4_maps.png`.

### PRJNA1139087 — post-chemoIO (this paper; n=4 sections, 12,979 tissue spots)

| Locked test | n sections | median | n neg / n pos | Wilcoxon p |
|---|---:|---:|---:|---:|
| Tumor same-spot CLDN4 vs CD8A | 4 | ρ = **−0.0027** | 2 / 2 | 0.625 |
| Tumor same-spot CLDN4 vs CD8A \| KRT8 | 4 | ρ = **−0.0090** | 3 / 1 | 0.875 |
| Tumor CLDN4 vs nearest CD8A+ (µm) | 4 | ρ = **+0.011** | 2 / 2 | 1.0 |
| CLDN4-high − low nearest CD8A (µm) | 4 | Δ = **0** | 1 / 2 (one 0) | 1.0 |
| Tumor CLDN4 vs hex ring-1 CD8A | 4 | ρ = **−0.017** | 2 / 2 | 0.875 |
| Ring-1 partial \| KRT8 | 4 | ρ = **−0.022** | 2 / 2 | 0.875 |
| CLDN4-high − low ring-1 CD8A | 4 | Δ = **−0.0006** | 3 / 0 | 0.25 |

Per-section (tumor spots unless noted):

| Section | Pathology | Tumor same-spot ρ | \| KRT8 | nn CD8A ρ | high−low nn µm | ring-1 ρ | \| KRT8 |
|---|---|---:|---:|---:|---:|---:|---:|
| PA08 | NMPR cold | −0.027 (p=0.27) | −0.018 | +0.028 | ≈0 | −0.041 | −0.046 |
| PA09 | NMPR cold | +0.024 (p=0.34) | −0.000 | **−0.075** (closer) | **−65** (p=0.015) | +0.066 | +0.047 |
| PA10 | pCR (CLDN4 almost gone) | +0.022 (p=0.38) | +0.295 | −0.005 | 0 | +0.008 | +0.003 |
| PA12 | NMPR, responsive-like | **−0.166** (p=1.7e−11) | **−0.135** | +0.064 | ≈0 | **−0.186** | **−0.175** |

PA08 all-spot CLDN4–CD8A is ρ=−0.183 (p=6e−26) and dies on the tumor mask / KRT8 residual. PA10 pCR has 708 CLDN4 UMIs total; high vs low is poorly defined. PA09 nearest-CD8 is the **opposite** of exclusion (CLDN4-high closer). PA12 is the only section with a tumor-restricted same-spot and ring-1 negative that survives KRT8; nearest-distance high−low is still ≈0.

**Post-ICI, computed:** CLDN4-high tumor spots are **not** farther from CD8 at the section-level unit. Signs are mixed. n=4 cannot reject a two-sided Wilcoxon at 0.05 (floor p=0.125). This is not a post-ICI exclusion claim, and it is not a matched pre/post claim.

Maps: `results/molcancer_visium_cldn4/maps/PA08_maps.png`, `PA09_maps.png`, `PA10_maps.png`, `PA12_maps.png`.

Tables: `results/molcancer_visium_cldn4/tables/prjna1139087_sections.csv`, `gse292299_sections.csv`.

## Answer

1. **Pre-ICI (open, unmatched GSE292299, n=16):** no same-spot or nearest-distance exclusion. Small ring-1 neighbor CD8A deficit after KRT8 (median ρ=−0.042, 12/16 negative, p=0.021). High−low deltas are 0.
2. **Post-ICI (this paper, n=4):** no section-level exclusion. Median tumor same-spot ρ≈0; nearest-CD8 high−low Δ=0. One of four sections (PA12) shows tumor same-spot and ring-1 negatives that survive KRT8; PA09 shows the opposite on distance.
3. There is **no open matched pre/post** NSCLC Visium for these patients. Yan 2024 post Visium and the JITC 2025 case remain controlled / restricted.

## What this does not claim

- It does not claim a matched pre/post pair.
- It does not claim Space Ranger image-based tissue detection.
- It does not re-score leftover T+B ρ from GSE292299 as the primary endpoint (CD8A distance / KRT8 residual are the locked tests here).
- n=4 post sections cannot reject a two-sided Wilcoxon at 0.05 even if all four signs agree.
- It does not treat PA12’s within-section Spearman as a cohort result.
