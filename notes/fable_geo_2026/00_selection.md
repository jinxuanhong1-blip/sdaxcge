# Dataset selection — GEO 2024–2026 human lung ICI series

Parallel slice. Outputs restricted to `notes/fable_geo_2026/`,
`scripts/fable_geo_2026/`, `results/w200/GEO_2026/`.

## Search
- NCBI E-utilities (`db=gds`), query: lung-cancer terms AND ICI terms AND
  `Homo sapiens` AND `gse[Entry Type]` AND PDAT 2024/01/01–2026/12/31.
- 156 hits → 155 GSE candidates after dropping the pre-specified exclusion set
  (GSE126044, GSE135222, GSE136961, GSE166449, GSE93157, GSE207422, GSE205335).
- Raw ESummary + candidate table: `results/w200/GEO_2026/geo_search_*.{json,tsv}`.
- Full sample-characteristic label scan over all 155:
  `results/w200/GEO_2026/label_scan_all.tsv` (17 series carry response/outcome
  fields).

## Target genes
TACSTD2 (TROP2) and CLDN4 are epithelial/tumor-cell genes → best measured in
tumor-tissue expression (bulk RNA-seq, spatial WTA/CTA, or tumor scRNA), not in
blood/PBMC assays.

## Datasets carrying BOTH tumor-tissue expression AND ICI response labels

| GSE | Cohort | Assay | Genes present | Response label | Used |
|-----|--------|-------|---------------|----------------|------|
| GSE261345 | ES-SCLC, chemo-immunotherapy (CANTABRICO) | GeoMx DSP, Cancer Transcriptome Atlas | TACSTD2 ✓, CLDN4 ✗ (not on CTA panel) | best RECIST response + progression/death + survival dates | yes |
| GSE261348 | ES-SCLC, chemo-immunotherapy (IMfirst) | GeoMx DSP, CTA | TACSTD2 ✓, CLDN4 ✗ | best RECIST response + progression/death + survival dates | yes |
| GSE233203 | NSCLC (EGFR-TKI resistant) ABCP atezolizumab combo | 10x scRNA-seq (7 patients) | TACSTD2 ✓, CLDN4 ✓ | therapeutic response (Response / Non-response) | yes |
| GSE292421 | pan-cancer incl. NSCLC, ICB resistance | bulk RNA-seq FPKM | TACSTD2 ✓, CLDN4 ✓ | TME immunophenotype (inflamed/excluded/desert) — NOT direct response | supportive only |
| GSE329813 | NSCLC neoadjuvant pembro+chemo (2026 leftover) | GeoMx DSP | TACSTD2 ✓, CLDN4 ✗ | MPR/NMPR in sample titles | yes (post-tx residual) |
| GSE292299 | NSCLC pretreatment Visium (2026 leftover) | Visium H5 (per-GSM, <2 GB) | TACSTD2 ✓, CLDN4 ✓ | Tx_Response R/NR | yes (n_NR=4, exploratory) |

## Checked but not usable for TACSTD2/CLDN4 vs response
- GSE309652 (NSCLC, anti-PD-(L)1, clean R/NR, n=72): NanoString **metabolism**
  panel (768 genes) — TACSTD2/CLDN4/EPCAM absent. Cannot measure target genes.
- GSE329813 (NSCLC neoadjuvant chemo-immuno spatial, n=127): first-pass
  characteristics scan only saw `batch`/`tissue`. **Correction (2026 leftover
  pass):** MPR/NMPR is in `!Sample_title`. Analyzed as leftover 2026
  (`leftover_2026_WRITEUP.md`). TACSTD2 ✓, CLDN4 ✗. Post-treatment residual
  tissue — not a pretreatment predictor.
- GSE253564 (durvalumab±RT, bulk lung tumor FPKM, full transcriptome): both
  genes present, but no per-sample response label in GEO metadata (treatment
  Arm1/Arm2 only). Excluded from response association.
- Blood/PBMC response-labeled series (GSE266219, GSE295969, GSE306542,
  GSE295601, GSE310370): epithelial genes not meaningfully expressed → excluded.
- Cell-line / mechanistic (GSE252437, GSE255144, GSE271377, GSE253718) and the
  pan-cancer scRNA atlas GSE218989: no per-patient ICI response usable here.

## Constraints honored
- Only open, processed supplementary files; each downloaded file < 2 GB.
- Excluded prior-slice accessions.
- No invented accessions — every GSE/GSM/URL comes from live E-utilities / GEO
  FTP responses saved under `results/w200/GEO_2026/`.
