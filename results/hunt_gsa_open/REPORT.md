# Hunt: Chinese GSA/NGDC open processed lung ICI — TACSTD2 / CLDN4 vs MPR

**Verdict: no.** There is no truly open, processed NGDC dataset in which human lung ICI samples have both **TACSTD2** or **CLDN4** expression and **MPR** (or any pathologic-response) labels.

That is a catalog fact, not a statistics failure. The Chinese neoadjuvant lung ICI RNA-seq studies that *would* support this test are almost all **HRA Controlled**. The few Open records are the wrong data type, the wrong species, the wrong treatment setting, or raw FASTQ.

Date of catalogs: 2026-08-16. Sources: [OMIX release list](https://ngdc.cncb.ac.cn/omix/releaseList) (10,216 records) and GSA-Human `finished.json` (7,032 records).

## What “open processed” meant here

| Allowed | Not allowed / not counted as a hit |
|---|---|
| OMIX **Open-access** tables with a public HTTPS file | HRA **Controlled** / DAC (listed only as near-misses) |
| GSA-Human **Open** *if* a processed matrix is actually there | Open HRA that is **FASTQ only** |
| Human lung + ICI + gene + MPR on the same samples | Mouse, cell line, metabolome, PK/ADA, mutation-only |
| | Title says Open but no downloadable file (HGRMP hold) |

HRA Controlled data were **not** requested or downloaded.

## Catalog counts

| Set | n |
|---|---|
| OMIX total / Open / Controlled | 10216 / 6059 / 4157 |
| OMIX title contains lung | 440 (138 Open) |
| OMIX title contains lung **and** ICI | 72 (**7 Open**) |
| HRA total / Open / Controlled | 7032 / 1385 / 5647 |
| HRA title contains lung **and** ICI | 50 (**2 Open**, both raw FASTQ) |
| HRA lung+ICI **Controlled** (near-miss list) | 48 (17 titles also say neoadjuvant/resectable) |

Title search is leaky (false PCR / “Open-Label trial” hits). Descriptions were read for every Open OMIX that matched lung, ICI, TACSTD2/TROP2/CLDN4, or MPR/pCR.

## Closest Open records (inspected, not imagined)

### Human lung ICI — processed, but cannot test the genes vs MPR

**OMIX001548 / OMIX001549** (PRJCA011037, Zhejiang Cancer Hospital; *Cell Death & Disease* 2022).  
Open. Cohort A: Nanostring TIME on 57 BRAF-mutant/WT NSCLC tumors. Cohort B: 417-patient ICI efficacy table (PFS/OS/`effect` 272 vs 145).  
- Nanostring panel is immune genes (BLK, CD19, GZMB, …). **TACSTD2 and CLDN4 are not on the panel.**  
- Cohort A has no ICI response labels. Cohort B has response and **no expression**.  
- Setting is advanced ICI, not neoadjuvant **MPR**.

**OMIX004639** (PRJCA018645). Open Phase II CRF: sintilimab + paclitaxel/platinum, first-line **advanced LUSC**. Clinical workbook, no expression, not resectable/MPR.

**OMIX002797 / OMIX002798** (PRJCA014383). Open HLX10 (serplulimab) PK/ADA/TMB/MSI, including ES-SCLC and an NSCLC303 nested zip. Pharmacology, not transcriptome.

**OMIX010537** (PRJCA041569, Guangdong Provincial People’s Hospital; *BMC Medicine* 2025). Open.  
- `OMIX010537-01.xlsx`: 14,144 mutation rows, 652 genes, ~922 patients. **No TACSTD2, no CLDN4.** No MPR column.  
- `OMIX010537-03` is a mislabeled TSV of **7 untreated NSCLC cell-line** FPKMs. Both genes are present (H441 TACSTD2 FPKM ≈ 914; H460 CLDN4 FPKM ≈ 14). That is not a patient ICI/MPR test.  
- 140 MB cell WES (`-02`) was not downloaded.

**OMIX013717**. Labeled Open “processed counts” after LINC00467 KD in LUAD cells. **No file URL on the page** (typical MOST HGRMP hold). Not actually downloadable.

**OMIX014557**. Open metabolome `.tar` for gut LDRT + PD-1 in mNSCLC. Not a gene matrix.

### Human lung ICI — Open HRA, but raw FASTQ only

**HRA003748** (PRJCA009346, Nanfang Hospital). This is the only Open **human NSCLC anti-PD-1 RNA-seq** on GSA-Human.  
- 65 tumor samples, 92 runs, Illumina PE150. HTTPS listing is `HRR*/HRR*_f1.fq.gz` (~5 GB/sample).  
- Metadata export has empty Treatment / Stage / Survival. Description says pretreatment specimens before anti-PD-1.  
- **No processed count matrix on NGDC.** BioProject page lists no companion OMIX.  
- Advanced ICI biomarker study, **not MPR**.  
Using it for TACSTD2/CLDN4 vs MPR would require aligning FASTQ *and* getting response labels that are not in the open metadata. That is outside “open processed.”

**HRA006427**. Open cell-line RNA-seq (CDK4/6 + LDRT + PD-1, Rb-deficient SCLC). 6 FASTQ runs. Not patients, not MPR.

### TACSTD2 exists on NGDC Open — wrong question

**HRA011882** (PRJCA041484). Open. Title is Trop2 in EGFR-TKI drug-tolerant persister cells (PC9, H1975; 16 FASTQ runs). TKI, not ICI. Raw only.

**OMIX019394**. Open mouse PDAC RNA-seq after a TROP2 aptamer-drug conjugate. `Tacstd2` and `Cldn4` are in the count table. Pancreas, not lung ICI/MPR.

**OMIX001948**. Trop-2 trial data — **Controlled**. Not used.

### Mouse lung + PD-1 with both genes (not MPR)

**OMIX014319** (PRJCA055741). Open mouse lung-cancer brain-metastasis RNA-seq after CFRT / HFRT / SF-HDRT (paper is RT-induced TLS + PD-1). Both genes are in the matrix at very low counts:

| gene | Ctrl 1–3 | CFRT 1–3 | HFRT 1–3 | SF-HDRT 1–3 |
|---|---|---|---|---|
| Tacstd2 | 0, 0.16, 0.04 | 0.10, 0.14, 0.07 | 0.06, 0.15, 0 | 0, 0.08, 0.16 |
| Cldn4 | 0, 0, 0.03 | 0.03, 0.07, 0 | 0.11, 0.17, 0 | 0, 0, 0.04 |

n=3/group, mouse, brain mets, RT groups. **Not human MPR.**

## What would have been the right studies (Controlled — not used)

48 HRA lung+ICI titles are Controlled. At least 17 are neoadjuvant/resectable on the title, including:

- HRA001033 — sc/bulk NSCLC during neoadjuvant immunotherapy  
- HRA002071, HRA002904, HRA003360 — neoadjuvant PD-1 ± chemo  
- HRA004391, HRA004588, HRA005191, HRA006493, HRA007118 — neoadjuvant / stereo-seq  
- HRA008054, HRA019979 — neoadjuvant camrelizumab / sintilimab  

Those are the cohorts where MPR is even defined. They require a GSA-Human DAC. They are listed in `hra_controlled_lung_ici_near_misses.tsv` and were **not** downloaded.

OMIX titles for Chinese neoadjuvant lung ICI (CheckMate-style, camrelizumab, sintilimab, tislelizumab, …) are likewise **Controlled**.

## CLDN4 specifically

No OMIX or HRA **title** mentions CLDN4 / claudin-4.  
CLDN4 appears in Open processed files only as:

1. mouse genes (`OMIX014319`, `OMIX019394`)  
2. 7 cell-line FPKMs (`OMIX010537-03`)

No Open human lung ICI expression table contains CLDN4.

## Limits (stated)

- GSA **CRA** (non-human) was not fully paged. Human clinical ICI lives in HRA/OMIX, not CRA.  
- Description crawl covered Open OMIX matching lung / ICI / gene / MPR keywords (198 pages), not all 6,059 Open OMIX records. A silent lung ICI dataset with none of those words in title or description would be missed.  
- “MPR” title hits were almost all PCR-assay false positives.  
- Paper supplements / GEO companions were not treated as GSA/NGDC hits. HRA003748’s BioProject does not link an OMIX matrix.

## Files in this folder

| file | what |
|---|---|
| `stats.json` | catalog counts + verdict |
| `verdicts.tsv` | one row per inspected candidate |
| `omix_open_shortlist.tsv` / `omix_open_details.tsv` | Open OMIX crawl |
| `hra_all_lung_ici_gene.tsv` / `hra_details.tsv` | HRA lung/ICI/gene |
| `hra_controlled_lung_ici_near_misses.tsv` | DAC-only list |
| `hra003748_samples.tsv` | 65 Open FASTQ samples (no response labels) |
| `omix014319_tacstd2_cldn4.tsv` | mouse gene rows |
| `omix010537_cell_tacstd2_cldn4.tsv` | 7 cell-line FPKMs |
| `gene_phenotype_scans.tsv` | file-level scan log |

Reproduce: `python3 src/parse_catalogs.py && python3 src/hunt_gsa_open.py` (needs network to ngdc.cncb.ac.cn; caches under `data/cache/`, ignored by git).
