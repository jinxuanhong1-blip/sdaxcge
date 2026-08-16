# REWORK A3 — GSE205335 malignant TACSTD2: MPR vs NMPR, 0-cell drop, %positive cutoffs

**Self-contained. Public processed GEO files only. Written to be read without the rest of the repo.**

**Verdict: NO — MPR not labeled; after dropping 0-malignant patients, malignant TACSTD2 is still null at every pre-specified %positive cutoff**

GSE205335 is advanced / palliative ICI, labeled with **RECIST 1.1 (PR / SD / PD / NE) only**.
There is **no MPR or NMPR field** in GEO sample characteristics or in the eLife paper.
The requested MPR vs NMPR contrast **cannot be run**. After dropping the 4 patients
with 0 author-annotated malignant cells, malignant TACSTD2 is still null for
responders vs non-responders: mean log1p(CP10K) p = 0.958
(this is the published hunt p ≈ 0.96), %positive (UMI ≥ 1) p = 0.713.
Every pre-specified cell-level %positive cutoff and every patient-level high/low
cutoff is also null on the primary set. An ADC+SQ-only sensitivity at UMI ≥ 1
is nominally p = 0.048 (n = 4 vs 8) and
is **not** used as the verdict. This is an **underpowered primary null
(n = 6 R vs 10 NR), not proof of no association**.

## Why this rework exists

A prior hunt (`results/hunt_gse205335/`) reported malignant-cell TACSTD2
R vs NR p = 0.96 (mean log1p CP10K; 6 vs 10 patients). Claim A3 in this project
is actually about **GSE207422** (neoadjuvant NSCLC, MPR vs NMPR). This folder
asks whether the GSE205335 null changes if we:

1. use **MPR vs NMPR if labeled**,
2. **drop patients with 0 malignant cells**,
3. score TACSTD2 as **%positive at pre-specified cutoffs** instead of only the mean.

## MPR vs NMPR — not labeled

| Source | Pathologic response (MPR / NMPR)? | What is labeled |
|---|---|---|
| GEO `!Sample_characteristics_ch1` (33 GSM, fetched for this run) | **No** | `recist`: PR / SD / PD / NE |
| eLife preprint 98366 (Ahn / Lee) | **No** | RECIST 1.1; R = PR, NR = SD+PD; no CR |
| Author cell table `GSE205335_Lung_IO_CellIdentity.txt` | **No** | cell lineage only |

This is a **stage III–IV / extensive-disease biopsy cohort on palliative ICI**,
not a neoadjuvant resection cohort. MPR requires a resected primary and a
percent residual viable tumor. That endpoint does not exist here. Substituting
RECIST for MPR would be a silent endpoint swap; it is not done.

Paper coding, reused: **R = PR** (no CR), **NR = SD + PD**, **NE excluded**
from response tests. Sensitivity: PR vs PD (drop SD).

## Analysis set and the 0-malignant drop

Author malignant call: `lineage.sub == "Malignant cells"` (n = 28,512 cells).
CNV was **not** re-inferred. Patient is the unit; samples from the same patient
are pooled.

| | n |
|---|---:|
| Patients in GEO | 26 |
| Dropped (0 malignant cells) | 4 |
| Patients with 1–19 malignant cells | 0 |
| Remaining with malignant cells | 22 |
| Of those, R / NR / NE | 6 / 10 / 6 |
| Primary test (drop 0, R vs NR) | 6 vs 10 |

Dropped zero-malignant patients: P2001 (PR, ADC, Normal LN), P2009 (PR, ADC, Normal LN), P2016 (PD, SQ, Normal LN), P3032 (PR, ADC, Normal Brain).

Dropout is **asymmetric**: 3 responder vs 1 non-responder
(0 NE). Three of the four are **normal tissue** (normal LN / normal brain)
from early-stage or metastatic patients — they have no tumor cells by sampling,
not by a failed annotation. Absence of captured tumor in some PR patients can
also be treatment or sampling; it shrinks the R arm and is not imputed.

**Drop-0 and the original hunt's ≥20-cell floor are the same set** in this
cohort: nobody has 1–19 malignant cells (P4001 is the smallest remaining, 27 cells).
Rework item 2 therefore **does not change n or the p = 0.96 mean-expression test**.
It only makes the exclusion rule explicit.

## Pre-specified %positive cutoffs

Locked before looking at cutoff p-values. A malignant cell is TACSTD2+ if it
meets the threshold. Patient score = 100 × (positive malignant cells / malignant cells).
Primary cutoff = **UMI ≥ 1** (detection). The other four are robustness, not a
search for a significant threshold.

Cell-level thresholds → patient %positive → two-sided Mann–Whitney U, R vs NR:

| Cell cutoff | median %pos R | median %pos NR | p | BH-FDR (5 cutoffs) | Cliff's δ (R−NR) |
|---|---:|---:|---:|---:|---:|
| UMI ≥ 1 (detection) | 70.9 | 52.0 | 0.713 | 1.000 | 0.133 |
| UMI ≥ 2 | 48.8 | 41.2 | 1.000 | 1.000 | 0.000 |
| UMI ≥ 3 | 34.9 | 37.2 | 0.792 | 1.000 | -0.100 |
| log1p(CP10K) ≥ 0.5 | 63.7 | 49.4 | 0.713 | 1.000 | 0.133 |
| log1p(CP10K) ≥ 1.0 | 45.4 | 42.2 | 0.958 | 1.000 | 0.033 |

Patient-level TROP2-high (using % UMI ≥ 1) → two-sided Fisher exact:

| High if | R high / R | NR high / NR | OR (high in R vs NR) | p |
|---|---:|---:|---:|---:|
| ≥ 10% UMI≥1 | 5/6 | 10/10 | 0.000 | 0.375 |
| ≥ 25% UMI≥1 | 4/6 | 9/10 | 0.222 | 0.518 |
| ≥ 50% UMI≥1 | 4/6 | 5/10 | 2.000 | 0.633 |
| ≥ 75% UMI≥1 | 2/6 | 1/10 | 4.500 | 0.518 |

Continuous metrics on the same primary patients:

| Metric | median R | median NR | p | Cliff's δ |
|---|---:|---:|---:|---:|
| mean log1p(CP10K) | 1.098 | 0.910 | 0.958 | -0.033 |
| %pos UMI ≥ 1 | 70.9 | 52.0 | 0.713 | 0.133 |
| pseudobulk CPM | 481.3 | 574.3 | 0.635 | -0.167 |

No primary cell-cutoff p < 0.05. No primary Fisher p < 0.05. Cliff's δ values
are small. Do not quote a "best" cutoff.

## Sensitivity (not used for the verdict)

| Set | n R vs NR | median %pos R | median %pos NR | p |
|---|---|---:|---:|---:|
| All histologies, drop 0 malignant, R vs NR | 6 vs 10 | 70.9 | 52.0 | 0.713 |
| ≥20 malignant cells, R vs NR | 6 vs 10 | 70.9 | 52.0 | 0.713 |
| ADC+SQ only, drop 0, R vs NR | 4 vs 8 | 74.1 | 52.0 | 0.048 |
| Author Core malignant cells, R vs NR | 4 vs 7 | 74.1 | 56.6 | 0.073 |
| All histologies, drop 0, PR vs PD | 6 vs 7 | 70.9 | 47.4 | 0.731 |
| ADC+SQ only, drop 0, PR vs PD | 4 vs 6 | 74.1 | 52.0 | 0.114 |

NSCLC-only (ADC+SQ) removes two SCLC PR patients with very low malignant
TACSTD2 (P1016, P1115) plus SCLC PD (P1025) and NUT SD (P1056). That
**raises** the R-arm median %positive. At UMI ≥ 1 this sensitivity is
**nominally p = 0.048**
(n = 4 vs 8;
Cliff's δ = 0.750).
It is **not the verdict**: n = 4 vs 8, uncorrected across six analysis sets
and five cutoffs; UMI ≥ 2 in the same NSCLC set is p = 0.214; NSCLC PR vs PD
is p = 0.114;
mean log1p in NSCLC is p = 0.368.
A single n=4 Mann–Whitney that kisses 0.05 after dropping neuroendocrine
tumors is a hypothesis, not a finding.

Author `core.patient == Core` is the paper's 14-sample / 11-patient clinical
core. It is a sensitivity, not the primary, because the rework asked for all
patients after a 0-malignant drop.

## Exploratory: T/NK (not the rework question)

Paired malignant vs T/NK TACSTD2 remains large in every patient with both
compartments (n = 22): median %pos
58.4% vs
1.2%, Wilcoxon
p = 4.77e-07. TROP2 is tumor-restricted
here regardless of ICI response.

Claim A3's GSE207422 number was ρ ≈ −0.40 to −0.50 between malignant TACSTD2
and T/NK presence. In GSE205335, Spearman of malignant %pos (UMI ≥ 1) vs
T/NK fraction among (malignant + T/NK):

| Set | n | ρ | p |
|---|---:|---:|---:|
| drop0_any_response | 22 | 0.185 | 0.411 |
| primary_drop0_RvsNR | 16 | 0.306 | 0.249 |
| nsclc_drop0 | 17 | -0.265 | 0.305 |

That is **not** the −0.40 to −0.50 claim, and this cohort is not GSE207422.

## Marker sanity

| Marker | % positive malignant | % positive T/NK |
|---|---:|---:|
| EPCAM (UMI > 0) | 69.3 | 2.5 |
| PTPRC / CD45 (UMI > 0) | 4.2 | 80.9 |

Consistent with the authors' malignant vs lymphocyte labels.

## Honest interpretation

1. **MPR vs NMPR is impossible here.** Do not write "MPR" on GSE205335 figures.
2. **Dropping 0-malignant patients does not move p = 0.96.** Those four patients
   already contributed no tumor-cell score; no one sits between 1 and 19 cells.
3. **On the primary set, %positive cutoffs do not rescue a response
   association.** Detection, UMI ≥ 2/3, and two log1p thresholds are all
   null; so are 10/25/50/75% patient high/low splits. This was pre-specified,
   not p-hacked.
4. **n = 6 vs 10 (4 vs 8 if NSCLC-only) is underpowered.** Only a large effect
   could have been seen. A null here is inconclusive, not "TROP2 is unrelated
   to ICI".
5. **SCLC in the R arm is a real confounder.** Two PR SCLC tumors are
   TACSTD2-low, as expected for neuroendocrine histology. Restricting to
   ADC+SQ makes UMI ≥ 1 nominally p = 0.048 (4 vs 8). That is a sensitivity,
   not a confirmed NSCLC effect: other cutoffs and PR-vs-PD in the same
   subset are not significant, and n = 4 cannot carry a claim.
6. **Asymmetric tumor-cell dropout.** 3 PR vs 1 PD had zero malignant cells.
   Two PR dropouts are normal LN. Do not treat the remaining R arm as a
   random sample of responders.
7. **mRNA ≠ protein ≠ ADC target occupancy.** UMI ≥ 1 is a transcript detection
   call, not an IHC H-score.
8. **Mixed sites and 3′/5′ chemistry** are confounded with patient; n is too
   small to adjust.
9. **Malignant labels are the authors'.** No independent inferCNV in this rework.

## Reproduce

```
pip install -r requirements.txt
# R with Matrix is required for the RDS extract
Rscript scripts/extract_gse205335_tacstd2.R
python scripts/rework_A3_GSE205335.py
```

The Python script will download the two GEO processed supplements into
`data/GSE205335/` (~500 MB UMI matrix, gitignored) if they are missing, gunzip
the RDS once (GEO file is double-gzipped), and call the R extract.

## Files

- `sample_table.tsv` — one row per patient: RECIST, histology, n malignant, all %pos cutoffs
- `patient_coverage.tsv` — the 0-cell drop and who is in the primary set
- `cutoff_tests.tsv` — Mann–Whitney for every set × metric
- `patient_highlow_fisher.tsv` — Fisher exact at 10/25/50/75%
- `tnk_fraction_spearman.tsv` — exploratory ρ
- `gsm_sample_metadata.csv` — GEO characteristics as fetched
- `summary.json` / `provenance.json`
- `figures/primary_r_vs_nr_and_cutoffs.png`
- `figures/per_patient_pctpos_bars.png`
- `figures/patient_highlow_50pct.png`

## Data

- GEO **GSE205335** processed: `GSE205335_Lung_IO_UMI_matrix.rds.gz` (dgCMatrix,
  33,714 genes × 96,505 cells, raw UMI) and `GSE205335_Lung_IO_CellIdentity.txt.gz`.
- Sample characteristics: NCBI GEO GSM6210624–GSM6210656 (`recist`, subtype, tissue).
- Raw FASTQ is EGA **EGAD00001008703** (controlled). Not used.
