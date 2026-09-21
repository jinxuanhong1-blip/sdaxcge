# Open pieces of the 11 ICI cohorts

## What is citable

The cohort list called “11 ICI cohorts” is Supplementary Data S1 of Lee et al., Science Advances 2024 (PMID [38295179](https://pubmed.ncbi.nlm.nih.gov/38295179), DOI 10.1126/sciadv.adj0785). That paper builds cell–cell communication models. It does not test CLDN4.

No published CLDN4 immune-checkpoint meta-analysis of these 11 cohorts was found. The citable Bessede ICI paper is about TROP2, not CLDN4: Bessede et al., Clinical Cancer Research 2024 (PMID [38048058](https://pubmed.ncbi.nlm.nih.gov/38048058)). It uses POPLAR and OAK RNA-seq (EGA EGAS00001005013, DAC EGAC00001002120; 891 tumors, atezolizumab arm n=405 in the paper) plus the institutional BIP study (NCT02534649). Those matrices are not open, so Bessede was not recomputed here. The closest published claudin-and-checkpoint paper is the 2026 TROP2/claudin-7 breast study (PMID 41932810); it is CLDN7, two breast cohorts, and not this list.

## Cohort access

Labels for the nine bulk cohorts are Lee Data S9 (CR/PR = responder, SD/PD = non-responder), saved in `data/ici_11_cohorts/lee2024_data_s9_labels.tsv`. Rank tests use a Mann–Whitney AUC for non-responder expression above responder expression. AUC 0.5 is no separation. CD8A is a direction check: responders should sit higher, which is AUC below 0.5.

| Cohort | Cancer | Accession | Open processed matrix | CLDN4 recomputed |
|---|---|---|---|---|
| VanAllen | melanoma | phs000452; cBio `skcm_dfci_2015` | yes | yes, 40/42 Data S9 patients have RNA |
| Liu | melanoma | cBio `mel_dfci_2019` | yes | yes, 119/119 |
| Gide | melanoma | PRJEB23709 FASTQ; cBio `mel_iatlas_gide_2019` | yes, iAtlas reprocess | yes |
| Hugo | melanoma | GSE78220 | yes, FPKM | yes, 27 baseline |
| Jung | NSCLC | GSE135222 | yes | yes, 27 |
| Cho | NSCLC | GSE126044 | yes, counts | yes, 16 |
| Kim | gastric | PRJEB25780 | no (78 RNA-seq FASTQ, no counts) | no |
| Mariathasan | urothelial | IMvigor210; cBio `blca_iatlas_imvigor210_2017` | yes | yes, 298 with RECIST |
| Prat | melanoma panel | GSE93157 | yes, NanoString | no, CLDN4 and TACSTD2 absent (CD8A present) |
| Sade-Feldman | melanoma scRNA | GSE120575 | yes, CD45+ TPM | gene row only |
| Jerby-Arnon | melanoma scRNA | GSE115978 | yes | no R/NR test (label is resistant vs untreated) |
| POPLAR+OAK | NSCLC | EGAS00001005013 | no | no |
| BIP | NSCLC | NCT02534649 | no | no |

Machine-readable copy: `results/ici_11_cohorts/cohort_table.tsv`.

## Recomputed CLDN4 and TACSTD2

Full precision is in `results/ici_11_cohorts/open_gene_tests.tsv`. Per-sample values are in `results/ici_11_cohorts/per_sample_open.tsv`.

| Cohort | n (R / NR) | CLDN4 AUC | CLDN4 p | TACSTD2 AUC | TACSTD2 p | CD8A AUC |
|---|---|---|---|---|---|---|
| Hugo GSE78220 | 27 (15/12) | 0.48 | 0.86 | 0.64 | 0.23 | 0.49 |
| Jung GSE135222 | 27 (8/19) | 0.48 | 0.90 | 0.50 | 1.00 | 0.22 |
| Cho GSE126044 | 16 (4/12) | 0.88 | 0.030 | 0.73 | 0.21 | 0.10 |
| Liu | 119 (47/72) | 0.57 | 0.18 | 0.52 | 0.67 | 0.48 |
| VanAllen | 40 (13/27) | 0.55 | 0.63 | 0.33 | 0.094 | 0.37 |
| Gide, pre+on (Data S9 n) | 91 (49/42) | 0.62 | 0.037 | 0.58 | 0.17 | 0.19 |
| Gide, pretreatment | 73 (40/33) | 0.62 | 0.070 | 0.57 | 0.29 | 0.20 |
| Mariathasan IMvigor210 | 298 (68/230) | 0.45 | 0.22 | 0.45 | 0.21 | 0.41 |

CD8A is higher in responders for Jung, Cho, Gide, and Mariathasan, so those label orientations are pointing the right way. Hugo and Liu CD8A are flat.

Open NSCLC bulk: Jung is null for CLDN4. Cho is the only open NSCLC cohort with higher CLDN4 in non-responders (4 vs 12). GEO’s own `patient response` string disagrees on GSM3589680 (Dis_17): GEO says responder, Data S9 says non-responder, and that sample is high for both CLDN4 and CD8A. Calling it a responder moves the Cho CLDN4 AUC to 0.73 (p=0.18). The open atezolizumab cohort (IMvigor210, urothelial) does not show higher CLDN4 or TACSTD2 in non-responders. That does not re-test Bessede on OAK.

Sade-Feldman GSE120575 is CD45-sorted. CLDN4 is nonzero in 48/16,291 cells (median 0). That is not a tumor CLDN4 measurement.

## Reproduce

```bash
python3 scripts/ici_11_cohorts/recompute_open.py
```

Needs pandas, scipy, and openpyxl. Downloads land in `/tmp/ici11` (override with `ICI11_CACHE`). cBioPortal is queried live for VanAllen, Liu, Gide, and IMvigor210.
