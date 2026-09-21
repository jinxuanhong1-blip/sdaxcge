# Open HRA/CNGB Chinese human lung scRNA: malignant CLDN4 vs STING / IFN / APM

Date of GSA-Human `finished.json`: pulled 2026-09-21 (7,312 records). OMIX release list pulled the same day (10,740 rows). Controlled HRA files were not requested.

## Access

Lung single-cell / single-nucleus / spatial **titles** in GSA-Human: **74**. Open: **2**. Controlled: **72**.

Title search misses a study whose title never says lung and never says single-cell. Every Open lung *disease* title in the same dump was read; the only Open single-cell/single-nucleus lung titles are the two below.

### Open HRA (the whole open list)

| Accession | What it is | Per-cell CLDN4 contrast |
|---|---|---|
| **HRA009335** (PRJCA031820) | Primary pulmonary lymphoepithelioma-like carcinoma, 4 snRNA runs (3 tumor, 1 adjacent). Sun Yat-sen. Open FASTQ only. | **Not scored.** PE150. In 4,000 read pairs from HRR2054880, R1 has a poly(T) stretch in the first 60 bp in 3,780/4,000 reads; R2 carries the Smart-seq TSO `AAGCAGTGGTATCAACGCAGAGTAC` in 215/4,000. Neither read is a 10x cell-barcode + UMI. One sample index (`TAGGACGT`) for the whole run. Processed matrix is not in the open deposit (paper: corresponding author on request). |
| **HRA011368** (PRJCA038849) | Single-nucleus RNA-seq of human **fetal** lung, twin-twin transfusion. Open FASTQ (HRR2360386 R1 ~10 GB, R2 ~25 GB). | Not a malignant-cell contrast. Not downloaded. |

HRA009335 is the only Open human **lung tumor** single-nucleus accession. It cannot support a malignant CLDN4-high vs CLDN4-low score.

### Controlled HRA lung sc/sn/spatial (not used)

72 titles. Raw files were not downloaded. Full table: `tables/hra_lung_scrna.tsv`. Examples that would have been the LUAD test if they were open: HRA001130 (lineage LUAD), HRA000154 / HRA000156 (subsolid / early LUAD), HRA002376 (preinvasive LUAD), HRA004391 (neoadjuvant immunotherapy LUAD), HRA012291 is not in the title-sc list under that accession; HRA008879 is the early-stage LUAD scRNA on PRJCA030984.

### Open processed tables (OMIX), same Chinese archives

These are Open OMIX releases, not Open HRA. The HRA raw run for the SCLC matrix is **Controlled**.

| OMIX | BioProject | HRA of that BioProject | Used? |
|---|---|---|---|
| **OMIX002441** | PRJCA006026 | **HRA001232 Controlled** | **Yes.** 5,025 cells, author counts + cell info. SCLC primary tumors, matched adjacent, one relapse. |
| OMIX011746 | PRJCA045688 | no HRA row in finished.json | No. 1.67 GB Seurat RDS of pleural effusion and blood, TB vs LUAD **immunity**. Not a malignant-cell matrix in the file title. |
| OMIX003147 | PRJCA015245 | not a tumor HRA | No. Embryonic lung. |
| OMIX008248 | PRJCA033575 | COVID white lung | No. Neutrophil / COVID, not tumor CLDN4. |
| OMIX004145 | PRJCA017221 | tuberculous pleural effusion | No. Pleural immune comparison, not lung-tumor malignant cells. |
| OMIX007208 | PRJCA029532 | pediatric Mycoplasma BAL | No. Bronchoalveolar immune atlas, not tumor. |

### CNGB / CNSA

CNP0005129 (Tianjin Chest Hospital; LUAD scRNA + spatial, STAS; 14 samples, 535 GB) is listed on CNGBdb as **Apply for data**, with metadata download only. That is controlled. It was not downloaded.

CNGBdb project search `lung single-cell` (235 hits) and `lung adenocarcinoma single cell` (25 hits) returned NCBI BioProject mirrors, not a public CNP tumor matrix. No other CNP lung-tumor scRNA page was confirmed Open in this pass. Absence from that search is not proof that no public CNP exists under a title that omits those words.

## Analysis that was possible: OMIX002441 SCLC malignant cells

Source: `OMIX002441-01.csv` (gene × cell counts) and `OMIX002441-02.csv` (cell info). Paper: Signal Transduction and Targeted Therapy, DOI 10.1038/s41392-022-01150-4 (Fuchou Tang / Peking University).

This is **small-cell lung cancer**, not LUAD. It is not added to the concordant-4 LUAD T/NK correlation. A result here does not rewrite that LUAD result.

**Malignant** = author `cell_type == tumor` (same cells have `NT == tumor`). n = 2,104 / 5,025. Immune and stromal cells are labeled `NT == normal` in this table even when the sample code is a tumor region, so they are not used as the low-CLDN4 malignant group.

Matrix QC: 24153 genes, header cells 5025, header matches metadata: True. First 30 genes look like integers: True. Sum of all genes in matrix column 0 = 34516; metadata nCount_RNA for that cell = 34516.

Score = **mean of log1p(count / nCount_RNA × 10,000)** over genes present in the matrix. Library size is the author `nCount_RNA`, checked against the column sum above.

**Primary contrast:** within each patient, Spearman correlation of malignant-cell CLDN4 log1p CPM with the score. Patients with ≥ 40 malignant cells are included. This keeps patients in whom almost every malignant cell is CLDN4-positive.

**Secondary contrast:** CLDN4 count > 0 vs = 0, only when both groups have ≥ 15 cells. Quartiles are not used. In this matrix most malignant cells are CLDN4-positive, so the zero class is a small tail and the patients with the highest CLDN4 (too few zeros) drop out of the secondary contrast.

A positive Spearman means higher CLDN4 goes with a higher score inside that patient's malignant cells. A positive pos-minus-zero delta means the CLDN4-positive malignant cells score higher than the CLDN4-zero malignant cells from the same patient. p is a two-sided Wilcoxon signed-rank test across patients (SciPy `zero_method='wilcox'`). It is descriptive. For 5 non-zero pairs that all share a sign, the smallest two-sided p this test can return is 0.0625. That floor is not a significance claim.

Patients with any malignant cells: 11. Spearman set (n≥40): **9** (P10, P11, P12, P13, P2, P3, P4, P5, P7). Pos-vs-zero set: **5** (P10, P11, P12, P13, P4).

### CLDN4 detection

| compartment | n cells | % CLDN4 > 0 | mean log1p CPM |
|---|---:|---:|---:|
| malignant_tumor | 2104 | 0.901 | 1.255 |
| author_normal_epithelial | 371 | 0.685 | 0.931 |
| Tcell | 1039 | 0.254 | 0.138 |
| myeloid | 840 | 0.396 | 0.230 |
| Bcell | 443 | 0.440 | 0.244 |

### Patient-level result

Delta = CLDN4-positive minus CLDN4-zero. For Spearman rows, the number is the within-patient Spearman of CLDN4 log1p CPM vs the score (not a high-minus-low delta).

| feature | kind | n patients | median | mean | n down / up | Wilcoxon p |
|---|---|---:|---:|---:|---|---:|
| spearman_CLDN4_STING | within_patient_spearman | 9 | -0.000373382 | 0.0292919 | 5/4 | 0.570312 |
| spearman_CLDN4_IFN | within_patient_spearman | 9 | 0.0705791 | 0.077596 | 2/7 | 0.128906 |
| spearman_CLDN4_APM | within_patient_spearman | 9 | 0.096975 | 0.081874 | 3/6 | 0.128906 |
| spearman_CLDN4_TJ | within_patient_spearman | 9 | 0.304251 | 0.286331 | 0/9 | 0.00390625 |
| spearman_CLDN4_log_nCount | within_patient_spearman | 9 | 0.0529089 | 0.095859 | 1/8 | 0.0195312 |
| spearman_CLDN4_nFeature | within_patient_spearman | 9 | 0.10928 | 0.130275 | 1/8 | 0.0117188 |
| STING | pos_minus_zero_log1p_cpm | 5 | 0.0154844 | 0.00839498 | 1/4 | 0.1875 |
| IFN | pos_minus_zero_log1p_cpm | 5 | 0.0273548 | 0.0317003 | 0/5 | 0.0625 |
| APM | pos_minus_zero_log1p_cpm | 5 | 0.066262 | 0.0712974 | 0/5 | 0.0625 |
| TJ | pos_minus_zero_log1p_cpm | 5 | 0.132694 | 0.150498 | 0/5 | 0.0625 |
| CLDN4 | gene_pos_minus_zero | 5 | 1.16967 | 1.28051 | 0/5 | 0.0625 |
| STING1 | gene_pos_minus_zero | 5 | 0.0257041 | 0.0240541 | 2/3 | 0.8125 |
| TBK1 | gene_pos_minus_zero | 5 | -0.0445202 | -0.0513628 | 3/2 | 0.3125 |
| STAT1 | gene_pos_minus_zero | 5 | 0.0843301 | 0.0782366 | 1/4 | 0.125 |
| ISG15 | gene_pos_minus_zero | 5 | 0.0612248 | 0.14985 | 0/5 | 0.0625 |
| MX1 | gene_pos_minus_zero | 5 | 0.059943 | 0.0465915 | 0/5 | 0.0625 |
| HLA-A | gene_pos_minus_zero | 5 | 0.02037 | 0.081286 | 2/3 | 0.3125 |
| B2M | gene_pos_minus_zero | 5 | 0.0339145 | 0.0892559 | 2/3 | 0.4375 |
| TAP1 | gene_pos_minus_zero | 5 | 0.011243 | 0.00867835 | 1/4 | 0.4375 |
| EPCAM | gene_pos_minus_zero | 5 | 0.230322 | 0.397811 | 0/5 | 0.0625 |
| PTPRC | gene_pos_minus_zero | 5 | -0.00841324 | -0.00984112 | 3/2 | 0.625 |
| CD3D | gene_pos_minus_zero | 5 | 0.00338655 | 0.00363874 | 2/3 | 0.8125 |
| CD68 | gene_pos_minus_zero | 5 | -0.0182194 | -0.0302414 | 4/1 | 0.3125 |

Genes present in each score are in `tables/sclc_family_summary.tsv`.

### How to read it

Author-malignant SCLC cells are mostly CLDN4-positive (about 90% of 2,104 cells; mean log1p CPM about 1.25). Adjacent epithelial cells labeled `normal` are lower (about 68%). Immune compartments are lower still (T cells about 25%). CLDN4 is measured in this matrix. It is not a LUAD cohort, and these patients are not added to concordant-4.

Primary result, within-patient Spearman, 9 patients with ≥ 40 malignant cells (P2, P3, P4, P5, P7, P10, P11, P12, P13):

- **STING** median ρ = −0.0004 (5/9 negative, 4/9 positive), Wilcoxon p = 0.57. Flat.
- **IFN** median ρ = +0.071 (7/9 positive), p = 0.13.
- **APM** median ρ = +0.097 (6/9 positive), p = 0.13.
- **TJ** (CLDN4 held out) median ρ = +0.30 (9/9 positive), p = 0.0039.
- CLDN4 vs log nCount median ρ = +0.053 (8/9 positive), p = 0.020. CLDN4 vs nFeature median ρ = +0.11 (8/9 positive), p = 0.012. The APM ρ is the same size as the nFeature ρ.

The two patients with the highest CLDN4 detection (P2, 98.5% positive, median count 48; P7, 95% positive) have negative IFN Spearman (−0.15 and −0.01). The positive median is not those high-CLDN4 tumors.

Secondary pos-vs-zero contrast is 5 patients (P4, P10, P11, P12, P13). CLDN4-zero malignant cells are the low-UMI tail: median nCount in the zero class is about half the positive class (P10 11,386 vs 32,438; P4 13,400 vs 23,516; P11 23,016 vs 42,187). IFN and APM means are slightly higher in the positive class (median deltas +0.027 and +0.066; 5/5 up; p = 0.0625, which is the floor of this test at n = 5). That contrast is library-size confounded. STING score delta median is +0.015 (4/5 up, p = 0.19). STING1, TBK1, and IRF3 do not move as one induced program (TBK1 median delta −0.045, 3/5 down).

Immune-marker check on the same 5 patients, pos minus zero: PTPRC median delta −0.008 (3/5 down), CD3D +0.003, CD68 −0.018 (4/5 down). The CLDN4-positive class is not an immune doublet.

TJ is a control family with **CLDN4 held out**. It is not the claim.

## What was not done

- No DAC application, no controlled FASTQ, no Cell Ranger on HRA009335 (no cell barcode in the reads).
- OMIX011746 RDS was not loaded.
- CNP0005129 was not downloaded.
- Concordant-4 GEO cohorts were not re-run and were not mixed with this SCLC table.

## Files

| file | |
|---|---|
| `tables/hra_lung_scrna.tsv` | every lung sc/sn/spatial HRA title, Open or Controlled |
| `tables/omix_open_human_lung_scrna.tsv` | Open human OMIX lung single-cell titles |
| `tables/sclc_patient_deltas.tsv` | one row per SCLC patient |
| `tables/sclc_family_summary.tsv` | paired deltas and Spearman |
| `tables/sclc_detection.tsv` | CLDN4 detection by author compartment |
| `tables/qc.json` | matrix QC |

Cells in the metadata table: 5025.
