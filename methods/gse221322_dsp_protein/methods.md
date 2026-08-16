# Methods — GSE221322 GeoMx protein DSP (additive spatial-protein readout)

Paper-ready paragraph. Numbers are from `results/gse221322_dsp/stats.txt` after `scripts/gse221322_dsp/analyze.py`. This analysis is extra spatial-protein evidence for the B6 spatial pattern (CLDN4-high epithelium sitting away from immune niches). It does not re-analyse Visium or CosMx.

---

## Dataset

Public GeoMx Digital Spatial Profiler (DSP) protein counts were obtained from GEO series **GSE221322** (Monkman, Kim, Mayer et al., *Immunology* 2023; PMID 37022147; platform GPL29263, NanoString GeoMx Human Protein for nCounter 2020). The series abstract describes Nanostring protein DSP on 42 TMA cores from immunotherapy-treated NSCLC; the deposited tables contain **48 paired tumour/stroma areas of illumination (AOIs)** from **41 patients** (16 ICI responders, 24 non-responders, 1 response not available). Tumour and stroma segments were the authors’ cytokeratin-positive and cytokeratin-negative masks on each core. No restricted-access files were used.

Two author-deposited matrices were downloaded: ERCC-corrected counts (`GSE221322_4301_protein_QC.csv`) and the analysis-ready normalized matrix (`GSE221322_4301_protein_norm.csv`). Patient identifier, ICI response, follow-up, vital status, AOI area and nuclei count were taken from the GEO series matrix. The primary matrix was the deposited normalized table, declared before immune correlations were inspected. Area- and nucleus-scaled log2 QC counts were a prespecified normalization sensitivity analysis only.

## Panel and target mapping

The deposited panel has **68 antibodies** (Immune Cell Profiling, IO Drug Target, Immune Cell Typing, Immune Activation Status, PI3K/AKT, MAPK and Cell Death modules, including three isotype controls). **CLDN4 and TROP2/TACSTD2 protein are not on the panel.** No claudin antibody is present. The closest barrier / epithelial proteins on the array are **EpCAM (TACSTD1)** and **PanCk**. Prespecified immune proteins were CD8, CD3 and PD-1; the neighborhood set was CD45, CD4, GZMB, CD20, FOXP3, CD68, HLA-DR, PD-L1 and CD56. Isotype controls (Ms IgG1, Rb IgG, Ms IgG2a) and housekeepers (Histone H3, GAPDH, S6) were used for QC, not as biological endpoints.

## Quality control

One tumour AOI (`4301 Protein | 026 | Tumour`) had 0 nuclei and an area of 54 µm² and was dropped (nuclei floor = 20). The paired stroma AOI from that core was retained for stroma-only tests and excluded from paired-core tests. After QC, **47 tumour AOIs (40 patients)** and **48 stroma AOIs (41 patients)** remained. Signal-to-noise was the target QC count divided by the mean of the three isotype controls. PD-1 sat near isotype background (median S/N ≈ 1.1; 0% of AOIs > 2× mean IgG) and was reported with that caveat.

## Statistics

The spatial unit was the AOI. Patient-mean AOIs (mean of cores from the same `patient_id` within a compartment) were the prespecified sensitivity analysis to limit pseudo-replication. Associations were two-sided Spearman rank correlations. Approximate 95% intervals used the Fisher *z* transform. Benjamini–Hochberg FDR was applied within each barrier protein to the 12 neighborhood tests (column `q_bh_within_barrier`). Tumour AOIs were also split at the median of EpCAM (or PanCk) and compared with two-sided Mann–Whitney tests. Partial Spearman correlations residualized ranks of EpCAM and each immune protein on PanCk, nuclei count, or Histone H3. Cross-compartment tests correlated tumour-segment EpCAM or PanCk with the same-core stroma immune protein (47 paired cores; 40 patients after averaging). ICI response contrasts were a published-cohort sanity check (Monkman et al. reported higher stromal EpCAM in responders) and were not used to test B6.

Software: Python 3, pandas, SciPy, matplotlib. Reproduction: `bash scripts/gse221322_dsp/00_fetch.sh && python3 scripts/gse221322_dsp/analyze.py`.
