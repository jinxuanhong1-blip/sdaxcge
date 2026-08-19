# RESULTS — Additive public Visium, CLDN4-only (Molecular Cancer 2025)

Post-neoadjuvant chemoimmunotherapy NSCLC Visium from Cui / Liu *Mol Cancer* 2025 (doi:10.1186/s12943-025-02287-w). **CLDN4-only.** No private 8-KL. No claim-failed.

## Accessions

| Role | Accession | Status |
|---|---|---|
| This paper, spatial raw | [PRJNA1139087](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1139087) / [SRP521546](https://www.ncbi.nlm.nih.gov/sra?term=SRP521546) | **Public** FASTQ (ENA). No GEO processed matrix. No images. |
| This paper, scRNA raw | PRJNA1068179 | Public (not used) |
| Paper validation scRNA | GSE207422 | Public scRNA (not Visium) |
| Open pre-ICI NSCLC Visium | [GSE292299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292299) | Public filtered H5 + `tissue_positions` (16 pre-treatment NSCLC). Raw withheld. |

Sample map (BioSample isolate):

| Patient | Alias | SRA | Pathology (paper) |
|---|---|---|---|
| PA08 | Tumor-1 | SRR29925401 | NMPR, “cold” |
| PA09 | Tumor-4 | SRR29925398 | NMPR, “cold” |
| PA10 | Tumor-2 | SRR29925400 | pCR |
| PA12 | Tumor-3 | SRR29925399 | NMPR but therapy-responsive histology |

Paper ST: 12,741 spots after QC. All four sections are **post** neoadjuvant PD-L1 blockade + cisplatin chemotherapy. There is **no matched pre Visium** in this deposit.

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
- **Section is the unit.** Report n, median, sign count, Wilcoxon p. n=4 post cannot support a strong signed-rank claim (Wilcoxon two-sided floor at n=4 is p=0.125).

Reproduce: `python3 methods/molcancer_visium_cldn4/count_fastq.py && python3 methods/molcancer_visium_cldn4/analyze.py`

## Results

*Filled after the public FASTQ count and GSE292299 matrices finish.*

## What this does not claim

- It does not claim a matched pre/post pair.
- It does not claim Space Ranger image-based tissue detection.
- It does not re-score leftover T+B ρ from GSE292299 as the primary endpoint (CD8A distance / KRT8 residual are the locked tests here).
- n=4 post sections cannot reject a two-sided Wilcoxon at 0.05 even if all four signs agree.
