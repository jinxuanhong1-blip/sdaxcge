# E-MTAB-13526 De Zuani NSCLC atlas — CLDN4 join test

**Decision: NO-JOIN** the concordant pool (GSE123902+GSE131907+GSE205335+GSE189357).

ADDITIVE. **CLDN4-only.** Public ArrayExpress / BioStudies accession
[E-MTAB-13526](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526)
(De Zuani, Xue, Park, et al., *Nat Commun* 2024, PMID 38782901). Treatment-naive
NSCLC resections plus matched non-involved lung. **Not an ICI cohort.** No public
MPR / RECIST / R-vs-NR labels. No dual-high TACSTD2×CLDN4.

The 45.51 GB / 58.74 GB author-annotated h5ads were **not downloaded**. This page
uses only processed files that are actually public and <~5 GB.

---

## Usable file list (downloaded, public, <~5 GB)

Source: BioStudies FTP
`https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/526/E-MTAB-13526/Files/`.
Machine table: `tables/used_files.tsv`. Full accession catalog with skip reasons:
`tables/public_file_catalog.tsv` (247 deposited files).

| File | Role | Size |
|---|---|---:|
| `E-MTAB-13526.sdrf.txt` | metadata (patient, FACS, disease, site) | 1.46 MB |
| `E-MTAB-13526.idf.txt` | study metadata | 6.2 KB |
| `P4_T2-features.tsv.gz` | shared Cell Ranger GRCh38 3.1.0 features (33,538 genes) | 0.30 MB |
| `P4_T2-matrix.mtx.gz` | CD235a− tumor mtx | 37.8 MB |
| `P4_T3-matrix.mtx.gz` | CD235a− tumor mtx | 47.3 MB |
| `P8_T2-matrix.mtx.gz` | CD235a− tumor mtx | 248.2 MB |
| `P15_T2-matrix.mtx.gz` | CD235a− tumor mtx | 103.9 MB |
| `P16_T2-matrix.mtx.gz` | CD235a− tumor mtx | 70.1 MB |
| `P17_T2-matrix.mtx.gz` | CD235a− tumor mtx | 249.9 MB |
| `P17_T3-matrix.mtx.gz` | CD235a− tumor mtx | 315.2 MB |
| `P18_T2-matrix.mtx.gz` | CD235a− tumor mtx | 200.9 MB |
| `P19_T2-matrix.mtx.gz` | CD235a− tumor mtx | 158.8 MB |
| `P20_T2-matrix.mtx.gz` | CD235a− tumor mtx | 170.1 MB |
| `P21_T1-matrix.mtx.gz` | CD235a− tumor mtx | 65.8 MB |
| `P21_T2-matrix.mtx.gz` | CD235a− tumor mtx | 11.4 MB |
| `P22_T1-matrix.mtx.gz` | CD235a− tumor mtx | 41.1 MB |
| `P23_T1-matrix.mtx.gz` | CD235a− tumor mtx | 305.6 MB |
| `P24_T1-matrix.mtx.gz` | CD235a− tumor mtx | 349.6 MB |
| **Usable total** | 15 CD235a− tumor matrices + features + SDRF/IDF | **2.38 GB** |

These 15 lanes are the only deposited **tumor** matrices with FACS = `CD235a-`
(RBC-depleted, not CD45-sorted). They are the only small processed files that can
contain both epithelium and T/NK. CLDN4 is in the features table (index 12425).

## Skipped (not used; not invented)

| File / class | Size / n | Why skipped |
|---|---|---|
| `10X_Lung_Tumour_Annotated_v2.h5ad` | **58.74 GB** | over the ~5 GB cap; user-locked skip |
| `10X_Lung_Healthy_Background_Annotated_v2.h5ad` | **45.51 GB** | over the ~5 GB cap; user-locked skip |
| CD45+ / MDSC / CD45++MDSC tumor matrices | 23 tumor lanes | FACS-sorted immune; no malignant epithelium; T/NK inflated |
| Background (`*_B*`) and donor (`D*`) matrices | 43 lanes | not a tumor malignant+T/NK contrast |
| All `*-barcodes.tsv.gz` | 81 × 20.9 MB | 6,794,880-barcode Cell Ranger whitelist; not required |
| Duplicate `*-features.tsv.gz` | 80 copies | identical 304,728-byte features table |
| ENA `ERP153272` FASTQ | raw | not a processed matrix |

No GEO processed subset of this atlas was used. No matrix was invented.

---

## Is there a usable malignant/epithelial + T/NK + CLDN4 compartment?

**Yes**, in the CD235a− tumor Cell Ranger matrices. There is **no** author
malignant label in the small files (that lives only in the skipped h5ads).
Malignant-like = marker epithelial (`EPCAM`, `KRT8/18/19`) minus normal-lung
markers (`SFTPA2`, `AGER`, `SCGB1A1`, `SCGB3A1`, `TPPP3` score ≤ 0.05). T/NK =
argmax T (`CD3D/E`, `CD2`) or NK (`NKG7`, `GNLY`, `FGFBP2`) among the same
lineage panel. Author-like QC while streaming: UMI 400–100,000, genes 180–6,000,
mito ≤ 20%.

| Item | n |
|---|---:|
| QC cells (15 CD235a− tumor lanes) | **276,018** |
| Patients with a CD235a− tumor lane | **12** |
| Eligible (≥10 malignant-like and ≥20 T/NK) | **10** |
| Excluded: P16 (2 malignant-like / 2 epithelial) | 1 |
| Excluded: P22 (1 malignant-like / 969 epithelial) | 1 |
| CLDN4 present | yes |
| IFN genes present (Hallmark IFNα ∪ IFNγ) | **222 / 224** (missing `MARCHF1`, `WARS1`) |
| MHC-I/APM genes present | **21 / 21** |
| Public ICI / MPR / R labels | **0** |

Eligible cell counts (patient-pooled lanes):

| Patient | histology | n_cells | n_malig-like | n_T/NK | T/NK frac | CLDN4 mean log1p |
|---|---|---:|---:|---:|---:|---:|
| P4 | LUSC | 5,456 | 228 | 3,976 | 0.729 | 1.097 |
| P8 | LUSC | 46,202 | 164 | 14,502 | 0.314 | 0.395 |
| P15 | LUAD | 8,911 | 66 | 3,900 | 0.438 | 1.405 |
| P17 | NSCLC_other | 50,938 | 10,368 | 19,110 | 0.375 | 0.724 |
| P18 | LUSC | 30,122 | 116 | 10,924 | 0.363 | 0.240 |
| P19 | LUSC | 22,112 | 94 | 8,354 | 0.378 | 0.228 |
| P20 | LUSC | 28,219 | 2,541 | 8,266 | 0.293 | 0.076 |
| P21 | LUAD | 4,709 | 31 | 377 | 0.080 | 0.192 |
| P23 | NSCLC_other | 33,771 | 1,778 | 3,121 | 0.092 | 0.493 |
| P24 | LUAD | 32,325 | 9,719 | 4,881 | 0.151 | 1.259 |

P16 and P22 are in the all-12 row only. They are not in the eligible n=10.

---

## CLDN4 vs T/NK (patient-level, CLDN4-only)

Unit = **patient**. Score = malignant-like mean log1p UMI. T/NK fraction among
all QC cells in CD235a− tumor lanes.

| Contrast | n | ρ | p |
|---|---:|---:|---:|
| **Eligible malig-like CLDN4 mean log1p vs T/NK** | **10** | **+0.418** | **0.229** |
| Eligible malig-like CLDN4 %pos vs T/NK | 10 | +0.309 | 0.385 |
| All 12 CD235a− patients (includes P16/P22) | 12 | −0.046 | 0.888 |

Sign vs concordant pool: pool is **inverse** (ρ < 0). This atlas is **weakly
positive and non-significant**. **Does not match.**

---

## Malignant IFN/MHC, Q4 vs Q1 (n allows a thin tail only)

IFN = mean log1p of 222 present Hallmark IFNα ∪ IFNγ genes in malignant-like
cells. MHC = mean log1p of the locked 21-gene MHC-I/APM panel. IFN/MHC = union.
Quartiles of malignant-like CLDN4 mean log1p on the same n=10 eligible patients.
Tails are **n_Q1=3 / n_Q4=3**. That is the honest quartile n; do not quote it as
a well-powered Q4-vs-Q1 test. Mann–Whitney two-sided p cannot fall below 0.1 at
3 vs 3.

| Contrast | n | n_Q1 / n_Q4 | mean Q1 | mean Q4 | Δ (Q4−Q1) | p |
|---|---:|---|---:|---:|---:|---:|
| CLDN4 Q4 vs Q1 · IFN | 10 | 3 / 3 | 0.117 | 0.216 | **+0.099** | 0.10 |
| CLDN4 Q4 vs Q1 · MHC-I/APM | 10 | 3 / 3 | 0.292 | 0.597 | **+0.305** | 0.10 |
| CLDN4 Q4 vs Q1 · IFN/MHC | 10 | 3 / 3 | 0.120 | 0.225 | **+0.105** | 0.10 |
| CLDN4 Q4 vs Q1 · T/NK fraction | 10 | 3 / 3 | 0.250 | 0.439 | +0.189 | 0.40 |

Continuous (same n=10; better powered than the 3-vs-3 tails):

| Contrast | n | ρ | p |
|---|---:|---:|---:|
| Malig-like CLDN4 vs IFN | 10 | **+0.806** | **0.00486** |
| Malig-like CLDN4 vs MHC-I/APM | 10 | **+0.842** | **0.00222** |
| Malig-like CLDN4 vs IFN/MHC | 10 | **+0.806** | **0.00486** |

Sign vs concordant pool: pool is **IFN/MHC down in CLDN4-high**. Here IFN/MHC is
**up** in CLDN4-high (continuous ρ > 0; Q4−Q1 Δ > 0). **Does not match.**

---

## Join rule

Join only if sign matches: T/NK inverse **and/or** IFN/MHC down in CLDN4-high.

| Gate | Result | Match? |
|---|---|---|
| T/NK inverse (eligible n=10) | ρ = **+0.418**, p = 0.229 | no |
| IFN/MHC down in high (continuous n=10) | ρ = **+0.806**, p = 0.00486 | no |
| IFN/MHC down in high (Q4 vs Q1, 3 vs 3) | Δ = **+0.105**, p = 0.10 | no |

**NO-JOIN.** The small processed files are usable and were used. The signs go
the other way from the concordant pool. This treatment-naive atlas is extra n
that does **not** enter GSE123902+GSE131907+GSE205335+GSE189357.

Reproduce: `python3 methods/emtab13526_cldn4_join/download.py && python3 methods/emtab13526_cldn4_join/extract.py && python3 methods/emtab13526_cldn4_join/analyze.py`.
Primary numbers: `tables/summary.json`, `tables/association_statistics.tsv`,
`tables/q4q1.tsv`, `tables/per_patient_metrics.tsv`.
