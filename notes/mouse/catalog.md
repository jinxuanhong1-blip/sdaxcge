# Mouse Lung ICI Dataset Catalog / 小鼠肺癌免疫检查点抑制 (ICI) 数据集目录

Focus genes: **Tacstd2** (TROP2, ENSMUSG00000051397) and **Cldn4** (ENSMUSG00000047501).
Scope: mouse (Mus musculus) lung tumor datasets with ICI / immune context.
Rule: processed matrices only; skip FASTQ / raw MS; skip files > 2 GB; no invented accessions.

| Accession | Source | Model / Tissue | Assay | ICI context | Groups | Processed file(s) | Decision |
|---|---|---|---|---|---|---|---|
| GSE239485 | GEO | Lewis Lung Carcinoma (LLC), subcutaneous; lung tumor | bulk RNA-seq | anti-PD-1 (+ poly I:C +/- anti-C5aR1) | Control vehicle (n=8); Poly I:C+anti-PD-1 (n=8); Poly I:C+anti-PD-1+anti-C5aR1 (n=8) | GSE239485_Processed_data.xlsx (7.1 MB, normalized log2, gene_name keyed) | ANALYZED (core, replicated) |
| GSE133604 | GEO | KrasG12D;p53-/- (KP) lung tumor; lung tumor | scRNA-seq (10x) | anti-PD-1 (+ Asf1a KO) | Ctrl; Ctrl+anti-PD-1; Asf1a-KO; Asf1a-KO+anti-PD-1 (1 sample/arm, ~6k cells each) | GSE133604_RAW.tar -> per-sample matrix.mtx.gz + shared genes/barcodes (78 MB) | ANALYZED (scRNA pseudobulk + single-cell) |
| GSE222158 | GEO | Oncogene-driven NSCLC (murine); lung tumor, FACS-sorted CD45+/CD3+ immune cells | scRNA-seq (10x) | anti-PD-1 (+ DC-CCL21 in situ vaccine) | CD45: Ctrl/PD1/DCvax/Combo; CD3: Ctrl/Combo (1 sample/arm) | GSE222158_RAW.tar -> per-sample barcodes/features/matrix (162 MB) | ANALYZED (immune-compartment control; tumor markers ~absent) |
| GSE129297 | GEO | SCLC (RPM) lung tumor; lung tumor | scRNA-seq (10x) | anti-PD-1 (+ CDK7i YKL-5-124) | Ctrl; anti-PD-1; CDK7i(YKL); anti-PD-1+CDK7i (1 sample/arm, unfiltered droplet matrices) | GSE129297_RAW.tar -> per-sample matrix.mtx.gz + shared features/barcodes (156 MB) | ANALYZED (scRNA pseudobulk + single-cell) |
| GSE241978 | GEO | CMT167 LUAD (C57BL/6 syngeneic) cell lines; lung adenocarcinoma cells | bulk RNA-seq (DE summary table) | immune-checkpoint pathway (AhR KO -> PD-L1/IDO1); no ICI drug | AhR-KO vs Cas9 control (DE statistics table, no per-sample matrix) | GSE241978_...CMT_KO_vs_Cas9Ctrl.xlsx (14 MB, DE table) | SUPPORTING (genetic checkpoint-pathway perturbation, not ICI treatment) |
| GSE330658 | GEO | Egfr-mutant lung cancer (C57BL/6J subcutaneous); lung tumor | bulk RNA-seq (per-sample TPM tables) | anti-PD-L1 context; deposited arms = PTX / anti-VEGF | Control; PTX; anti-VEGF; PTX+anti-VEGF (n=2 each). anti-PD-L1 arm NOT in processed deposit | GSE330658_RAW.tar -> 8 per-sample xlsx with TPM (41 MB) | ANALYZED (descriptive, n=2/arm) |
| GSE197260 | GEO | Egfr-mutant NSCLC (syngeneic C57BL/6J); subcutaneous lung tumor | bulk RNA-seq (TPM) | anti-PD-1 (4H2) + anti-VEGFR2 (DC101) after EGFR-TKI (gefitinib) | vehicle_d3, gef_d3, gef_d14, gef_vehicle_d21, gef_dc101_d21, gef_4h2_d21, gef_comb_d21 (n=1/arm) | GSE197260_RNAseqTPM_MM_EGFR-TKI-CD8.txt.gz (0.3 MB, symbol-keyed TPM) | ANALYZED (descriptive trajectory, n=1/arm) |
| E-MTAB-13704 | ArrayExpress/BioStudies | Lung GEMM (LUGEMM14); lung tumor in situ | bulk RNA-seq | anti-PD-L1 mono + combos (ATRi, VEGFRi, Cisplatin/aPD-L1/aCTLA4) | vehicle(5); aPD-L1(5); ATRi(4); ATRi/aPD-L1(5); VEGFRi/aPD-L1(4); Cisplatin/aPD-L1/aCTLA4(4) | GEMMS_raw_counts.csv (14 MB, Ensembl-keyed counts) + sdrf | ANALYZED (core, replicated, richest ICI design) |
| GSE297630 | GEO | LLC subcutaneous; lung tumor cells (recovered post-treatment) | microarray (Clariom S Mouse) | anti-PD-1 (tolerant/surviving cells vs control) | Control (C, n=3); anti-PD-1 tolerant (P, n=3) | GSE297630_processed_data.xlsx (5.3 MB, per-sample log2 RMA + author stats) | ANALYZED (core, replicated) |
| GSE297632 | GEO | LLC subcutaneous; lung tumor (whole) | scRNA-seq (10x) | anti-PD-1 (tolerant vs control) | Control (~8.4k cells); anti-PD-1 (~9.4k cells) | GSE297632_RAW.tar -> per-sample barcodes/features/matrix (232 MB) | ANALYZED (scRNA pseudobulk + single-cell; companion of GSE297630) |
| PXD059688 | PRIDE | NSCLC syngeneic (C57BL/6J); lung tumor | LC-MS/MS proteomics | anti-PD-1 +/- high-dose ascorbic acid (Con/P/AA/AP) | Con, P(aPD1), AA(ascorbic acid), AP(AA+aPD1) x ~3 (raw MS only) | Only raw .raw/.mgf/.msf (>2GB each) + one single-condition 93-protein mzTab (0.3 MB, no Tacstd2/Cldn4) | CATALOGED, NOT ANALYZED (no <2GB cross-condition protein-quant table; raw MS skipped per rules) |
| GSE309199 | GEO | RPM SCLC (Rb1/Trp53/MycT58A) lung tumor in situ | bulk RNA-seq | anti-PD-1 +/- HDAC inhibitor entinostat | Ctrl (n=3); aPD-1 (n=3); entinostat (n=3); aPD-1+ent (n=3) | GSE309199_Mouse_Azam.TPMcalculator.raw_counts.tsv.gz (0.5 MB) | ANALYZED (treatment, n=3; NOT per-mouse ICB R/NR) |
| GSE330941 | GEO | LLC subcutaneous (ICI-refractory WT vs Ago2KO ICI-sensitized) | bulk RNA-seq | model-level ICI sensitivity; RNA at d12 **without ICI on these samples** | WT (n=4); Ago2KO (n=4) | GSE330941_filtered_tablecounts_tpm.csv.gz (0.7 MB) | ANALYZED (model-level ICI-sensitivity proxy, NOT R vs NR) |
| GSE261890 | GEO | mouse NSCLC, immunotherapy-sensitive vs resistant | spatial RNA | ICI-sensitive vs resistant + veh/aPD-1 | 2 spatial slides | RAW.tar 1.8 GB; rds 730 MB | CATALOGED, NOT ANALYZED (spatial, n=2; not a bulk R/NR matrix) |
| GSE76628 | GEO | **EXCLUDED** — Ad-VEGF-A164 flank angiogenesis / gastric-cancer stromal signatures | microarray | anti-VEGFR (DC101/G6), not PD-1/PD-L1 | flank angiogenic sites (n=78) | CEL only | **EXCLUDED: not lung, not ICI, not ICB response** |

## Notes
- **ANALYZED (core, replicated)**: E-MTAB-13704, GSE239485, GSE297630 — replicated treated-vs-control designs used for formal statistics.
- **ANALYZED (scRNA)**: GSE129297, GSE133604, GSE297632, GSE222158 — pseudobulk + single-cell correlations.
- **ANALYZED (descriptive)**: GSE330658 (n=2/arm), GSE197260 (n=1/arm) — reported as fold changes/trends, underpowered for p-values.
- **ANALYZED (added for ICB-response search)**: GSE309199 (treatment), GSE330941 (model-level ICI sensitivity).
- **SUPPORTING**: GSE241978 — AhR-KO (genetic checkpoint-pathway perturbation), DE table only.
- **NOT ANALYZED**: PXD059688 (raw MS >2 GB); GSE261890 (spatial only).
- **EXCLUDED**: GSE76628 — gastric/flank VEGF stroma, not mouse lung ICI RNA.
- **Honest ICB-response gap**: no public mouse *lung* ICI RNA with per-mouse responder vs non-responder labels was found. Search log: `notes/mouse/icb_response_search.md`.

Provenance: metadata in `notes/mouse/raw_meta/`; download URLs + sizes + checksums in `notes/mouse/data/MANIFEST.tsv`; gene ID map in `notes/mouse/gene_map.json`.
