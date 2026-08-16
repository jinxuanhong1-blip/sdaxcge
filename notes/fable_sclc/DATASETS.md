# SCLC immunotherapy datasets — verified accessions

All accessions below were manually verified on **2026-08-16** by resolving the
repository landing page (GEO/EGA/cellxgene/cBioPortal) or the paper's data
availability statement. "Used" marks datasets analyzed in this slice
(`scripts/fable_sclc/`, total processed download **1.49 GB < 2 GB**).

## Public (open-access processed data)

| Dataset | Modality | Accession / ID | Content | ICI context | Used | Size used |
|---|---|---|---|---|---|---|
| George et al. 2015, *Nature* | bulk RNA-seq (FPKM) | cBioPortal `sclc_ucologne_2015` (processed; raw is controlled, see below) | 81 SCLC tumors with expression | none (pre-ICI cohort); used for SCLC-I subtype biology | yes | 9 MB (full FPKM matrix via REST API) |
| Jiang et al. 2016, *PLoS Genet* | bulk RNA-seq (log2-normalized) | GEO [GSE60052](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60052) (SRA PRJNA257389) | 79 SCLC tumors + 7 normal lungs | none; subtype/immune context | yes | 9.9 MB |
| CANTABRICO trial DSP, 2024 | spatial (GeoMx CTA ~1.8k genes) | GEO [GSE261345](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261345) | 121 "Full ROI" segments, 26 ES-SCLC patients, pre-treatment FFPE | first-line **durvalumab** + platinum/etoposide (EudraCT 2020-002328-35); per-patient RECIST + PFS/OS dates in GEO characteristics | yes | 4.1 MB counts + series matrix |
| IMfirst trial DSP, 2024 | spatial (GeoMx CTA) | GEO [GSE261348](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261348) | 175 ROIs, 32 ES-SCLC patients | first-line **atezolizumab** + platinum/etoposide; clinical in GEO characteristics | yes | 5.4 MB counts + series matrix |
| Chan et al. 2021, *Cancer Cell* (HTAN MSK) | scRNA-seq | cellxgene collection [62e8f058-9c37-48bc-9200-e767f318a8ec](https://cellxgene.cziscience.com/collections/62e8f058-9c37-48bc-9200-e767f318a8ec); HTAN portal | 147,137 cells / 42 donors ("Combined samples"); 77,143 SCLC cells / 20 donors | 9/20 SCLC donors ICI-exposed (treatment metadata) | yes | 1.46 GB h5ad |
| CheckMate-032 WES somatic calls | WES (VCF) | EVA PRJEB25807 / PRJEB25808 | somatic variant calls (TMB analysis) | nivolumab ± ipilimumab, 2L+ | no (no expression) | — |
| George et al. 2015 murine arrays | microarray | GEO GSE69091 | mouse SCLC cell lines | none | no | — |

Notes on verification: GSE261345/GSE261348 supplementary listings and the
per-GSM clinical characteristics were confirmed directly against the GEO FTP
(`ftp.ncbi.nlm.nih.gov/geo/series/GSE261nnn/...`); the cellxgene dataset list
and file sizes were confirmed via the cellxgene curation API; cBioPortal
profiles were confirmed via `https://www.cbioportal.org/api`.

## Controlled-access (application required)

| Dataset | Modality | Accession | ICI context | Access route |
|---|---|---|---|---|
| **IMpower133** (Gay et al. 2021 *Cancer Cell*; Nabet et al. 2024 *Cancer Cell*) | bulk RNA-seq + trial clinical | EGA [EGAS50000000138](https://ega-archive.org/studies/EGAS50000000138) — EGAD50000000195 (clinical: arm, OS, PFS, BOR, PD-L1 IHC, subtype), EGAD50000000196 (log2(TPM+1), n=271), EGAD50000000312 (raw FASTQ) | first-line atezolizumab + CE vs CE (randomized) — the key public-trial dataset for TACSTD2/CLDN4 × ICI-response testing | DAC: devsci-dac-d@gene.com (Roche/Genentech); additional clinical via vivli.org |
| **George et al. 2015 raw data** | WGS / SNP6 / transcriptome FASTQ | EGA [EGAS00001000925](https://ega-archive.org/studies/EGAS00001000925) | none | EGA DAC (U Cologne) |
| **CheckMate-032 RNA-seq** (Rimm/Schalper JTO 2023, antigen presentation) | bulk RNA-seq, n=286 pre-treatment | no public accession; not deposited | nivolumab ± ipilimumab (2L+) | request via BMS / study authors |
| **Chan et al. 2021 raw scRNA** | scRNA FASTQ/BAM (Level 1–2) | dbGaP [phs002371](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002371) / HTAN MSK | mixed treatment | dbGaP application (processed counts are open, see above) |
| Cai et al. 2024 *Cell Discovery* | GeoMx DSP RNA+protein, 44 SCLC | GSA-Human [HRA004312](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA004312) | none (resected, treatment-naive); NE-immune heterogeneity | GSA-Human DAC (restricted) |
| Sun et al. 2025 *npj Precis Oncol* | CosMx SMI single-cell spatial, 18 ES-SCLC (31 biopsies) | GSA [PRJCA036290](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA036290) | chemo-immunotherapy, PR vs PD | GSA (raw restricted under Chinese regulations) |
| 2025 Visium FFPE SCLC study | Visium spatial | GSA-Human [HRA015375](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA015375) | none | GSA-Human DAC |
| **IMfirst / CANTABRICO patient-level trial data** beyond GEO fields | clinical | not deposited | atezolizumab / durvalumab first-line | Roche data sharing (go.roche.com/data_sharing) |

## Gaps

- No fully open bulk-RNA + ICI-response SCLC cohort exists; IMpower133
  (EGAS50000000138) is the highest-value controlled target, and CheckMate-032
  RNA-seq has no accession at all.
- The GeoMx cohorts used here are the only open expression datasets with
  per-patient ICI outcome we could verify; their CTA panel lacks CLDN4
  (TACSTD2 is present).
