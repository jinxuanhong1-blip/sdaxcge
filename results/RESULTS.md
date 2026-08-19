# Open GeoMx lung WTA/CTA hunt: tumor-ROI CLDN4 vs immune-ROI CD8

**CLDN4 is on WTA, not CTA.** Tumor-ROI CLDN4 vs matched immune-ROI CD8A was runnable in five open WTA series besides GSE221733. The only significant CD8A correlation is GSE265899 LUAD (ρ=0.438, p=0.002, n=47 pairs). Yale/Greek NSCLC ICI WTA (GSE271689 / GSE292098) and LUSC-IPF / pleomorphic WTA are compatible with no linear association for CD8A. CTA lung series (GSE174743, GSE174749, GSE249568, GSE261345, GSE261348) were downloaded and lack CLDN4. STAS GeoMx (HRA012209), lepidic-vs-acinar DSP (HRA005794), and Yan et al. GeoMx+Visium (controlled Zenodo) were skipped as controlled. GSE221733 was not analyzed here.

Plots: `results/GSE265899_CLDN4_vs_CD8A.png`, `results/GSE271689_CLDN4_vs_CD8A.png`, `results/GSE292098_CLDN4_vs_CD8A.png`, `results/GSE289483_CLDN4_vs_CD8A.png`, `results/GSE305762_CLDN4_vs_CD8A.png`.

CLDN4-only. GSE221733 left to another agent. Controlled-access series skipped. Private 8-KL not used.

## Accessions tried

| Accession | Panel | Disease | Open | CLDN4 | n AOI (tumor / immune) | n pairs | Tumor CLDN4 vs immune CD8A Spearman | vs CD8 T signature | Notes |
|---|---|---|---|---|---|---|---|---|---|
| GSE265899 | WTA | LUAD (n=8 patients, PD-L1 spatial) | yes | yes (95/95 AOI>0) | 95 (48 / 47) | 47 | ρ=0.438 p=0.0021 n=47 | ρ=0.308 p=0.0352 n=47 | tumor-ROI vs immune-ROI |
| GSE292098 | WTA | NSCLC ICI (Greek validation cohort) | yes | yes (136/136 AOI>0) | 136 (68 / 68) | 68 | ρ=0.080 p=0.5170 n=68 | ρ=0.151 p=0.2197 n=68 | Companion WTA of GSE271689 paper; UQ CTA is GSE221733 (skipped). |
| GSE271689 | WTA | NSCLC ICI (Yale + Greek DCC; Greek also in GSE292098 xlsx) | yes | yes (566/579 AOI>0) | 579 (212 / 175) | 166 | ρ=-0.070 p=0.3699 n=166 | ρ=0.107 p=0.1711 n=166 | tumor-ROI vs immune-ROI |
| GSE174749 | CTA | NSCLC (single tumor, 191-ROI grid, PanCK vs TME) | yes | no | 376 (185 / 191) | 175 | ρ=NA p=NA n=0 | ρ=NA p=NA n=0 | tumor-ROI vs immune-ROI |
| GSE174743 | CTA | NSCLC (5 tumors, PanCK vs TME) | yes | no | 223 (105 / 103) | 97 | ρ=NA p=NA n=0 | ρ=NA p=NA n=0 | tumor-ROI vs immune-ROI |
| GSE289483 | WTA | Pulmonary pleomorphic carcinoma (n=9) | yes | yes (114/114 AOI>0) | 114 (37 / 24) | 20 | ρ=-0.163 p=0.4926 n=20 | ρ=-0.585 p=0.0067 n=20 | tumor-ROI vs immune-ROI |
| GSE305762 | WTA | LUSC associated with IPF (n=6) | yes | yes (52/52 AOI>0) | 52 (7 / 7) | 7 | ρ=-0.214 p=0.6445 n=7 | ρ=-0.393 p=0.3833 n=7 | tumor-ROI vs immune-ROI |
| GSE309894 | WTA | ALK+ NSCLC alectinib early resistance (n=6 cases, 34 ROIs) | yes | yes (34/34 AOI>0) | 34 (0 / 0) | 0 | ρ=0.037 p=0.8342 n=34 | ρ=-0.322 p=0.0636 n=34 | Within-ROI mixed compartments, not tumor-vs-immune segmentation. |
| GSE261348 | CTA | ES-SCLC IMfirst chemo-ICI | yes | no | 175 (0 / 0) | 0 | ρ=NA p=NA n=0 | — | CTA panel does not include CLDN4. |
| GSE261345 | CTA | ES-SCLC CANTABRICO chemo-ICI | yes | no | 121 (0 / 0) | 0 | ρ=NA p=NA n=0 | — | CTA panel does not include CLDN4. |
| GSE249568 | CTA (+ custom CLDN18/MET probes) | METex14 NSCLC, 1 patient, pre/post tepotinib, 94 ROIs | yes | no | 94 (0 / 0) | 0 | ρ=NA p=NA n=0 | — | CLDN4 not on this CTA custom panel (CLDN18 variants present; CD8A present). |
| GSE250509 | GeoMx protein (immune panels) | METex14 NSCLC paired with GSE249568 RNA | yes | no | 94 (0 / 0) | 0 | ρ=NA p=NA n=0 | — | CD8 protein present; CLDN4 protein absent and paired CTA GSE249568 lacks CLDN4 RNA, so tumor-CLDN4 vs immune-CD8 protein is not runnable. |
| GSE221733 | CTA | NSCLC ICI (UQ) | yes | not tested | — | — | — | — | Assigned to another agent (UQ GeoMx CTA NSCLC ICI). |
| GSE221322 | GeoMx protein | NSCLC ICI (UQ) | yes | not tested | — | — | — | — | Protein companion of GSE221733 (same UQ study); skipped with that series. |
| HRA012209 | CTA | NSCLC STAS vs non-STAS | no | not tested | — | — | — | — | STAS GeoMx CTA (Frontiers Pharmacol 2025). GSA-Human controlled; skip controlled. Figshare Data Sheet 1 is supplementary tables, not a public count matrix. |
| HRA005794 | DSP (protein + limited RNA) | Stage IA LUAD lepidic vs acinar | no | not tested | — | — | — | — | PMC10844893 lepidic vs acinar DSP. GSA-Human controlled; skip controlled. No GEO series for the DSP counts. |
| Zenodo 13901289 / 14728962 | GeoMx + Visium | NSCLC neoadjuvant chemo-ICI | no | not tested | — | — | — | — | Yan et al. Nat Genet 2024 (often cited from Nat Commun spatial papers) GeoMx+Visium NSCLC ICI. Application-controlled Zenodo; skip controlled. No GEO WTA/CTA release. |
| private 8-KL |  |  | no | not tested | — | — | — | — | Out of scope (no private 8-KL). |

## Runnable paired analyses (tumor CLDN4 vs immune CD8)

### GSE265899

- Panel: WTA; LUAD (n=8 patients, PD-L1 spatial)
- Pairs: 47 tumor-ROI / immune-ROI
- CLDN4 detected in 95/95 AOIs; CD8A in 95/95
- Spearman log2(Q3+1) tumor CLDN4 vs immune CD8A: ρ=0.438, p=0.0021, n=47
- Same vs immune CD8 T signature (14 genes): ρ=0.308, p=0.0352, n=47

### GSE292098

- Panel: WTA; NSCLC ICI (Greek validation cohort)
- Pairs: 68 tumor-ROI / immune-ROI
- CLDN4 detected in 136/136 AOIs; CD8A in 136/136
- Spearman log2(Q3+1) tumor CLDN4 vs immune CD8A: ρ=0.080, p=0.5170, n=68
- Same vs immune CD8 T signature (14 genes): ρ=0.151, p=0.2197, n=68

### GSE271689

- Panel: WTA; NSCLC ICI (Yale + Greek DCC; Greek also in GSE292098 xlsx)
- Pairs: 166 tumor-ROI / immune-ROI
- CLDN4 detected in 566/579 AOIs; CD8A in 570/579
- Spearman log2(Q3+1) tumor CLDN4 vs immune CD8A: ρ=-0.070, p=0.3699, n=166
- Same vs immune CD8 T signature (14 genes): ρ=0.107, p=0.1711, n=166
- Patient-level means (59 patients): ρ=-0.080, p=0.5484, n=59

### GSE289483

- Panel: WTA; Pulmonary pleomorphic carcinoma (n=9)
- Pairs: 20 tumor-ROI / immune-ROI
- CLDN4 detected in 114/114 AOIs; CD8A in 114/114
- Spearman log2(Q3+1) tumor CLDN4 vs immune CD8A: ρ=-0.163, p=0.4926, n=20
- Same vs immune CD8 T signature (13 genes): ρ=-0.585, p=0.0067, n=20

### GSE305762

- Panel: WTA; LUSC associated with IPF (n=6)
- Pairs: 7 tumor-ROI / immune-ROI
- CLDN4 detected in 52/52 AOIs; CD8A in 52/52
- Spearman log2(Q3+1) tumor CLDN4 vs immune CD8A: ρ=-0.214, p=0.6445, n=7
- Same vs immune CD8 T signature (14 genes): ρ=-0.393, p=0.3833, n=7

## Methods

- Counts used as deposited Q3/Neg-normalized matrices when provided; DCC files were mapped through `Hs_R_NGS_WTA_v1.0.pkc` and Q3-scaled to 1000.
- Tumor vs immune pairing uses PanCK+/CK vs CD45+/TME segments when both exist for the same ROI/patient.
- CD8 T signature is the mean z-score of log2(count+1) for CD8A/B, CD3D/E/G, GZMA/B/H, PRF1, NKG7, CCL5, CST7, CD2, CD7 (genes present on the panel).
- CTA Human Cancer Transcriptome Atlas PKC (`GeoMx_Hs_CTA_v1.0`) has 1812 targets and **does not include CLDN4**.
- Correlations use Spearman on log2(x+1). No claim is made for series that were not downloaded.

