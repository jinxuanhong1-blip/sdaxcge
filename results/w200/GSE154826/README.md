# GSE154826 — LCAM vs TACSTD2 (processed/open only)

## TL;DR (honest)

**TACSTD2 (TROP2) is not associated with LCAM-hi vs LCAM-lo in the author immune-pseudobulk DE.**

- Author file `DE_LCAMhi_vs_LCAMlo_pseudobulk.rd`: TACSTD2 log2FC = **+0.327**, p = **0.251**, padj = **0.359**.
- TACSTD2 is vanishingly rare in this matrix (UMI frequency ~1e-5; 342 cells above the author threshold). That is expected: **GSE154826 is CD45+ immune-enriched CITE-seq**, not a tumor-epithelium atlas.
- TACSTD2 is an **epithelial** gene here (immune vs epithelial-gated DE: l2fc = **−5.53**), together with EPCAM (−7.83) and CLDN4 (−7.42). It is not an LCAM module gene.
- LCAM sanity checks pass: CXCL13, SPP1, IGHG1, PDCD1 are strongly LCAM-hi.
- **Tumor-cell TACSTD2 vs LCAM cannot be measured in the clustered Mount Sinai object.** Epithelial/endothelial/fibroblast barcodes were gated as `epi_endo_fibro_doublet` and excluded from lineage analysis. Per-cell tumor TROP2 would require the huge unannotated 10x MTX dump or the full `lung_ldm.rd` UMI object; both were **honestly skipped**.

Do not cite this series as evidence that TROP2-high tumors are LCAM-hi or LCAM-lo.

## Dataset

- **GSE154826** / Leader, Grout et al., *Cancer Cell* 2021 (PMID 34767762): CITE-seq of 361,929 cells from 35 early-stage NSCLC lesions (Mount Sinai).
- Design: tumor and non-involved lung, **CD45+ bead or FACS enrichment** (plus some dead-cell depletion / CD2+ TCR libraries). Series matrices contain **metadata only** (8.6K + 6.3K).
- LCAM-hi = co-occurring **T_activated (PDCD1+CXCL13+)**, **IgG plasma**, and **MoMac-II (SPP1+)**; LCAM-lo = B, AM, cDC2, AZU1_mac, Tcm/naive_II, cDC1 (author `figure_5abcd_s5a.R`).

## What was used (processed, open, compact)

| File | Why |
|---|---|
| `DE_LCAMhi_vs_LCAMlo_pseudobulk.rd` | Author LCAM-hi vs LCAM-lo immune pseudobulk DE (20,235 genes) |
| `immune_vs_ep_de.csv` | Author immune vs epithelial-gated DE (confirms TACSTD2 is epithelial) |
| `cell_metadata.csv` + `annots_list.csv` + `table_s1_sample_table.csv` | Cell→cluster→sample; composition LCAM scores |
| GEO `GSE154826_sample_annots.csv.gz` + series matrices | Public metadata / accession check |

URLs and MD5s: `provenance.json`.

## What was skipped (honest)

- **SRA/FASTQ** SRP251372 / PRJNA609924 (huge raw).
- **GEO `GSE154826_amp_batch_ID_*.tar.gz`** (~3 GB raw feature-barcode MTX, no cluster labels). Processed in the GEO sense, but they do not add a labeled epithelial TACSTD2 measurement without re-clustering 362k cells.
- **Dropbox `lung_ldm.rd`** (full `umitab` for 361,929 cells). Not required once the author LCAM DE table is in hand.

## Methods

1. Download the compact public tables listed above (script caches under `/tmp/gse154826_open`).
2. Read `DE_total` from the author LCAM-hi vs LCAM-lo RDS (`fg` = LCAM-hi, `bg` = LCAM-lo).
3. Extract TACSTD2 and control genes; rank by |log2FC|.
4. Confirm epithelial restriction in `immune_vs_ep_de.csv` (negative l2fc = higher in epithelial-gated cells).
5. Recompute composition LCAM scores from cell metadata using the author within-lineage normalization and cluster sets. A median split of `lcam_difference` on V2-beads tumor samples is only a descriptive label; it is **not** paired to TACSTD2 expression (no compact expression matrix).

Reproduction: `python3 scripts/w200_gse154826/analyze_lcam_tacstd2.py`

## Results

### 1. TACSTD2 is null in LCAM-hi vs LCAM-lo immune DE

| gene | log2FC | p | padj | freq LCAM-lo | freq LCAM-hi |
|---|---:|---:|---:|---:|---:|
| **TACSTD2** | **0.327** | **0.251** | **0.359** | 9.6e-6 | 1.2e-5 |
| CXCL13 | 2.70 | ~0 | ~0 | 3.4e-5 | 2.3e-4 |
| SPP1 | 1.99 | ~0 | ~0 | 4.6e-4 | 1.8e-3 |
| IGHG1 | 4.74 | ~0 | ~0 | 3.5e-4 | 9.4e-3 |
| PDCD1 | 0.86 | ~0 | ~0 | 1.0e-5 | 2.0e-5 |
| EPCAM | 0.52 | 0.0037 | 0.0095 | 1.9e-6 | 3.1e-6 |
| CLDN4 | 0.88 | ~0 | ~0 | 2.2e-6 | 5.0e-6 |

EPCAM/CLDN4 “hits” are residual transcripts at frequencies even lower than TACSTD2. They are not a tumor-cell TROP2 measurement.

### 2. TACSTD2 is epithelial, not LCAM

Immune vs epithelial-gated DE: TACSTD2 l2fc = −5.53 (immune 1.1e-5 vs epi 5.7e-4). PTPRC is the opposite (+3.68). TACSTD2 appears in the paper’s Figure 2A DC gene list next to CLEC9A; that is a heatmap annotation list, not evidence that TROP2 is a cDC1 or LCAM marker.

### 3. Composition LCAM axis (no TACSTD2)

V2-beads tumor samples only (author primary LCAM axis): **30 samples / 26 patients**. LCAM-hi and LCAM-lo composition scores are anti-correlated (Spearman ρ = −0.47, p = 0.0086). A median split of the difference gives 15 / 15 samples. Scores are in `v2_beads_tumor_lcam_scores.csv`. This restates the paper’s immune-module split; it does not test TROP2.

## Files

| File | Contents |
|---|---|
| `stats.txt` | Numeric summary |
| `lcam_hi_vs_lo_focus_genes.csv` | TACSTD2 + controls from author DE |
| `lcam_hi_vs_lo_pseudobulk_de.csv` | Full 20,235-gene author DE |
| `immune_vs_epithelial_focus_genes.csv` | TACSTD2/EPCAM/CLDN4 vs immune genes |
| `sample_lcam_composition_scores.csv` | Per-sample LCAM composition |
| `v2_beads_tumor_lcam_scores.csv` | Author-matched V2-beads tumor subset |
| `lcam_hi_vs_lo_volcano.png` | Volcano; TACSTD2 highlighted |
| `focus_genes_lcam_log2fc.png` | Focus-gene log2FC bars |
| `tacstd2_immune_vs_epithelial.png` | Epithelial restriction |
| `v2_beads_tumor_lcam_scores.png` | Composition LCAM-hi vs LCAM-lo |
| `provenance.json` | URLs, MD5s, skip list |

## Caveats

- **Wrong compartment for tumor TROP2.** CD45 enrichment plus doublet gating removes the cells that actually express TACSTD2.
- Author DE p-values of 0 are R underflow, not exact zeros.
- Composition LCAM calls are a median split of `lcam_hi − lcam_lo` on V2-beads tumor samples after author normalization. They are a convenience label, not the paper’s published patient table (Table S1 on GitHub has no LCAM column).
- No ICI outcome in this resection series; LCAM–response claims in the paper use other cohorts (including non-public POPLAR).
- CITE-seq ADT panel is immune-focused; there is no TROP2 protein antibody in the public record used here.
