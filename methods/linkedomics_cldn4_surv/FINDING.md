# Finding — LinkedOmics CPTAC CLDN4 RNA/protein vs OS/PFS (LUAD, LSCC)

Additive slice. **CLDN4 only.** Open LinkedOmics CPTAC pan-cancer freeze v1.2.
ImmuneScore / ESTIMATE / xCell / CIBERSORT correlations are **already reported** in the LUAD and LSCC CPTAC slices and were **not re-audited** here.

## Verdict

Clinical tables are **open** (`LUAD_survival.txt`, `LSCC_survival.txt`, plus meta).
These OS/PFS times are **treatment-naive surgical follow-up**, not immunotherapy outcomes. No ICI labels exist in this freeze.

**No CLDN4 RNA or protein association with OS or PFS reaches p<0.05** on the pre-specified continuous Cox (per 1 SD) in either histology. Median-split log-rank is likewise null at α=0.05. Closest: LUAD CLDN4 protein OS log-rank p=0.087 (15 events / 75). Event counts are small; a null is **inconclusive**, not proof of no prognostic value.

Primary model = continuous Cox **per 1 SD** of CLDN4 (pairwise-complete). Median split is sensitivity only (ties at the median coded high). HR >1 means higher CLDN4, higher hazard.

## Survival table

| Cohort | Layer | Endpoint | n | Events | Cox HR per 1 SD (95% CI) | Cox p | Median-split log-rank p | High / low n (events) |
|---|---|---|---:|---:|---|---:|---:|---|
| LUAD | CLDN4 RNA | OS | 105 | 23 | 1.13 (0.76–1.70) | 0.543 | 0.455 | 53/52 (10/13) |
| LUAD | CLDN4 RNA | PFS | 104 | 23 | 1.00 (0.66–1.52) | 0.993 | 0.974 | 51/53 (11/12) |
| LUAD | CLDN4 protein | OS | 75 | 15 | 1.15 (0.71–1.86) | 0.563 | 0.0875 | 37/38 (9/6) |
| LUAD | CLDN4 protein | PFS | 74 | 16 | 1.30 (0.77–2.18) | 0.325 | 0.514 | 37/37 (8/8) |
| LSCC | CLDN4 RNA | OS | 94 | 23 | 0.99 (0.64–1.51) | 0.945 | 0.885 | 49/45 (12/11) |
| LSCC | CLDN4 RNA | PFS | 92 | 20 | 0.80 (0.50–1.27) | 0.337 | 0.361 | 48/44 (9/11) |
| LSCC | CLDN4 protein | OS | 66 | 14 | 0.87 (0.47–1.59) | 0.643 | 0.581 | 33/33 (8/6) |
| LSCC | CLDN4 protein | PFS | 66 | 14 | 0.74 (0.42–1.32) | 0.313 | 0.924 | 33/33 (7/7) |

Machine-readable copy: `results/survival.tsv`.

## Data (open, HEAD HTTP 200)

Filenames from the LinkedOmics CPTAC-pancan-LUAD / CPTAC-pancan-LSCC download tables, served from `cptac-pancancer-data` S3 freeze `data_freeze_v1.2_reorganized`. Phenotype / ImmuneScore files were **not** downloaded.

| Cohort | File | HTTP | Bytes | Role |
|---|---|---:|---:|---|
| LUAD | `LUAD_survival.txt` | 200 | 2485 | downloaded |
| LUAD | `LUAD_meta.txt` | 200 | 12232 | downloaded |
| LUAD | `LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt` | 200 | 35757367 | downloaded |
| LUAD | `LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | 200 | 20132337 | downloaded |
| LSCC | `LSCC_survival.txt` | 200 | 2471 | downloaded |
| LSCC | `LSCC_meta.txt` | 200 | 12051 | downloaded |
| LSCC | `LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt` | 200 | 35179517 | downloaded |
| LSCC | `LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | 200 | 20631271 | downloaded |

LUAD: Gillette *Cell* 2020 (PMID 32649874), treatment-naive resected adenocarcinoma. RNA CLDN4 `ENSG00000189143.9` n=110; protein `ENSG00000189143.9` n=79 (NA=31). Survival table OS events=24, PFS events=23.

LSCC: Satpathy *Cell* 2021 (PMID 34358469), newly diagnosed resected squamous carcinoma, no prior chemo/RT. RNA CLDN4 `ENSG00000189143.9` n=108; protein `ENSG00000189143.9` n=78 (NA=30). Survival table OS events=25, PFS events=20.

CLDN4 protein missingness is TMT dropout, not a join error. All protein tests are pairwise-complete.

Analytic n is smaller than the matrix n because follow-up time is missing (or ≤0) for some cases, including deaths deposited without `OS_days` (LUAD 1, LSCC 2). Those rows cannot enter a Cox or log-rank model. LUAD RNA OS uses 105/110 cases (23 events; table lists 24). LSCC RNA OS uses 94/108 (23 events; table lists 25).

## Methods

- Gene row matched by Ensembl prefix `ENSG00000189143` (versioned IDs in the freeze).
- Join key = case ID (matrix columns × `case_id` in survival). Time unit = days as deposited.
- Cases with time ≤0 or missing time/event/CLDN4 are dropped for that test.
- **Primary:** Cox PH on CLDN4 z-scored within the pairwise-complete set (HR per 1 SD). 95% CI = exp(β ± 1.96 SE). Two-sided Wald p from `statsmodels.PHReg`.
- **Sensitivity:** median split (ties at median → high); two-sample log-rank.
- Per-unit Cox is in `results/survival.tsv` for audit; RNA and protein live on different log2 scales, so per-SD is the comparable number.
- No ImmuneScore residualization, no stage-adjusted Cox, no optimal cutpoint search.
- Not ICI survival. Do not meta-analyse these HRs with atezolizumab / pembrolizumab series.

## What this does not claim

- It does not re-test CLDN4 vs ImmuneScore / CD8 / IFN signatures.
- It does not claim a prognostic or predictive ICI biomarker.
- It does not treat a null with ~15–23 events as evidence of no effect.
- It does not use LinkedOmics LinkFinder web p-values; all numbers are recomputed from the freeze files.

## Outputs

- `results/survival.tsv` — n / events / HR / CI / p for every test
- `results/coverage.json` — row IDs, missingness, event totals
- `results/sample_level.tsv` — case-level CLDN4 + OS/PFS (no immune columns)
- `results/figures/km_cldn4_os_pfs.png` — eight median-split KM panels
- `data/manifest.json` — URLs, HTTP status, SHA256 (local cache, not committed)

