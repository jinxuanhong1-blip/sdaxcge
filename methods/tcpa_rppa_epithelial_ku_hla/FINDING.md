# Finding — TCPA RPPA LUAD/LUSC: CLDN4 or epithelial antibodies vs DNA-PKcs / Ku / STAT1 / HLA

**Protein vs protein.** TCPA RPPA500, release 5.0, **disease-specific level 4** (replicates-based normalization). **LUAD** and **LUSC** kept separate. Public patient-tumor matrices from `https://tcpa.drbioright.org/rppa500/` (`TCGA-LUAD-L4`, `TCGA-LUSC-L4`).

These are TCGA resection RPPA profiles. **No ICI labels.** Do not read the correlations as immunotherapy outcomes, spatial exclusion, or a CLDN4 result.

## Verdict

CLDN4 is **not on the panel**. The only claudin antibody is **Claudin-7**, and it does **not** track Ku80.

| Pair | LUAD | LUSC |
|---|---|---|
| Claudin-7 vs Ku80 | ρ=**−0.030**; p=**0.57**; q=**0.78**; n=**365** | ρ=**+0.035**; p=**0.52**; q=**0.52**; n=**328** |
| E-cadherin vs Ku80 | ρ=**+0.236**; p=**5.2×10⁻⁶**; q=**1.5×10⁻⁵**; n=**365** | ρ=**+0.528**; p=**5.4×10⁻²⁵**; q=**1.6×10⁻²⁴**; n=**328** |
| EMA (MUC1) vs Ku80 | ρ=**−0.015**; p=**0.78**; q=**0.78**; n=**354** | ρ=**−0.160**; p=**0.0042**; q=**0.0062**; n=**318** |
| Any of the three vs HLA-DQA1 | **not tested** (0/365 quantified) | **not tested** (0/328 quantified) |

**What is present and tested:** Ku80 (`KU80`, XRCC5) against the three epithelial antibodies that are actually quantified: Claudin-7 (`CLAUDIN7`, CLDN7), E-cadherin (`ECADHERIN`, CDH1), and EMA (`EMA`, epithelial membrane antigen / MUC1).

**What holds:** E-cadherin protein is positively associated with Ku80 protein in both histologies. The LUSC correlation is large (ρ=+0.528, n=328). The LUAD correlation is smaller (ρ=+0.236, n=365) and in the same direction. Both survive BH within the three Ku80 tests that could be run.

**What does not hold:** Claudin-7 is null vs Ku80 in both histologies. EMA is null in LUAD and **inverse** in LUSC (q=0.006). Do not write “epithelial markers track Ku80,” and do not borrow the E-cadherin result for CLDN4. The three antibodies are not one axis in LUSC (E-cadherin vs EMA ρ=−0.064, p=0.25, n=318).

**What cannot be tested on this slide:**

| Requested protein | LUAD / LUSC |
|---|---|
| CLDN4 | no antibody id |
| DNA-PKcs (PRKDC) | no antibody id. DNA-named ids on the list are DNA ligase IV and DNA polymerase gamma |
| Ku70 (XRCC6) | no antibody id. Ku80 is the only Ku subunit |
| STAT1 | no antibody id. STAT3, phospho-STAT3 Y705, and STAT5-alpha are different antibodies and were not substituted |
| HLA class I (HLA-A/B/C, B2M) | no antibody id |
| HLA-DQA1 | id exists, but **every lung tumor is NA** (LUAD 0/365, LUSC 0/328). The same case IDs are empty in the pan-cancer RPPA500 file. Other tumor types in that file do have HLA-DQA1 values (5,433 observed), so this is a lung-row gap, not a failed join |

Do not pool LUAD and LUSC. Do not fill DNA-PKcs, STAT1, or HLA with XRCC1, PARP1, ATM, ATR, IRF1, STING, cGAS, CIITA, or PD-L1.

## Honest n

TCPA documents LUAD as 365 patient tumors / 416 antibodies, and LUSC as 328 / 416. The CouchDB protein-id list has 433 documents in each histology; the extra documents include all-NA placeholders (HLA-DQA1, alpha-catenin).

| | LUAD | LUSC |
|---|---:|---:|
| Tumors in disease L4 | **365** | **328** |
| Claudin-7 quantified | 365 | 328 |
| E-cadherin quantified | 365 | 328 |
| EMA quantified | **354** (11 NA) | **318** (10 NA) |
| Ku80 quantified | 365 | 328 |
| HLA-DQA1 quantified | **0** | **0** |
| Claudin-7 ∩ Ku80 (tested) | **365** | **328** |
| E-cadherin ∩ Ku80 (tested) | **365** | **328** |
| EMA ∩ Ku80 (tested) | **354** | **318** |
| Any ∩ HLA-DQA1 | **0** | **0** |

Do not write n=365 for the LUAD EMA test, or n=328 for the LUSC EMA test. Sample ids are 12-character TCGA case barcodes. The release labels the matrix patient tumor. No duplicate case/protein rows.

BH q is within histology across the **3** primary tests that had data (the Ku80 column). The three HLA-DQA1 tests were not entered, because there was no endpoint variance.

Recomputed Spearman matches the ρ stored on the TCPA protein documents. Maximum absolute difference on the tested pairs: **4.4×10⁻⁸**.

## Epithelial antibodies are not interchangeable

| Pair | LUAD | LUSC |
|---|---|---|
| Claudin-7 vs E-cadherin | ρ=+0.463; p=8.4×10⁻²¹; n=365 | ρ=+0.287; p=1.2×10⁻⁷; n=328 |
| Claudin-7 vs EMA | ρ=+0.281; p=7.5×10⁻⁸; n=354 | ρ=+0.234; p=2.6×10⁻⁵; n=318 |
| E-cadherin vs EMA | ρ=+0.302; p=7.1×10⁻⁹; n=354 | ρ=**−0.064**; p=**0.25**; n=318 |

N-cadherin was not used as an epithelial marker. P-cadherin and p63 were left out of the primary set. Alpha-catenin’s id is present and entirely NA in both lung histologies, so it was not scored.

Sensitivity, not in the BH family: within-histology z-mean of Claudin-7, E-cadherin, and EMA (complete trio only).

| | LUAD | LUSC |
|---|---|---|
| z-mean vs Ku80 | ρ=+0.067; p=0.21; n=**354** | ρ=+0.217; p=9.7×10⁻⁵; n=**318** |
| z-mean vs HLA-DQA1 | n=0 | n=0 |

The LUSC z-mean correlation is the E-cadherin signal sitting inside a score whose components disagree. It is not a second epithelial result.

## Methods (this slice)

- Primary matrix: TCPA disease-specific level 4 abundances (`TCGA-LUAD-L4-<protein>`, `TCGA-LUSC-L4-<protein>`), release 5.0. Single-histology L4 is the TCPA recommendation when a disease was profiled on more than one RPPA batch.
- The pan-cancer file `TCPA_TCGA_RPPA500.tsv` (7,828 samples) was used only to confirm missingness for the same case ids. It was not the correlation matrix.
- Predictor if CLDN4 is absent: Claudin-7, E-cadherin, EMA. No RNA stand-in.
- Endpoints kept only when that antigen is the antibody: Ku80 for Ku; HLA-DQA1 for the only HLA id. DNA-PKcs, Ku70, STAT1, and HLA class I had no id, so those tests were not run.
- Test: two-sided Spearman, pairwise complete. BH within histology on the tests that ran.
- Cohorts kept separate. No meta-analysis. No ICI labels.

## What this does not claim

- It does not test CLDN4 protein. Claudin-7 is not a CLDN4 result.
- It does not test DNA-PKcs, Ku70, STAT1, or HLA class I.
- It does not treat HLA-DQA1 as measured in lung. The id is an all-NA document here.
- It does not treat E-cadherin–Ku80 as an epithelial-program result. EMA and Claudin-7 do not follow it.
- It does not use N-cadherin, STAT3, IRF1, STING, XRCC1, PARP1, ATM, or ATR as substitutes.
- It does not test ICI response, PFS, or OS.
- It does not pool LUAD and LUSC.
- It does not describe spatial exclusion or effector-cell state. This is bulk tumor RPPA.

## Outputs

- `results/presence.tsv` — requested ids vs the disease L4 catalog
- `results/panel_name_audit.tsv` — ids matching claudin / Ku / DNA-P / STAT / HLA / cadherin substrings
- `results/not_standins.tsv` — neighboring antibodies that were not substituted
- `results/n_table.tsv` — quantified n
- `results/pancan_lung_missingness.tsv` — same case ids in the pan-cancer file
- `results/spearman.tsv` — primary, sensitivity, and absent-endpoint rows
- `results/epithelial_concordance.tsv`
- `results/sample_scores.tsv`
- `results/summary.json`
- `results/fig_epithelial_vs_ku80_hladqa1.png`
- `results/protein_ids_LUAD.txt`, `results/protein_ids_LUSC.txt`

```bash
python3 methods/tcpa_rppa_epithelial_ku_hla/download.py --outdir data/tcpa_rppa_epithelial_ku_hla
python3 methods/tcpa_rppa_epithelial_ku_hla/analyze.py --data data/tcpa_rppa_epithelial_ku_hla --outdir methods/tcpa_rppa_epithelial_ku_hla
```
