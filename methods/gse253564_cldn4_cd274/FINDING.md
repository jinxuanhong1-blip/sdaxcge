# FINDING — GSE253564 leftover CLDN4 vs CD274 / HLA continuous

**Additive only. CLDN4 only.** Public pre-treatment FPKM from the leftover neoadjuvant **durvalumab ± SBRT** series ([GSE253564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564); Altorki et al., *Cell Reports Medicine* 2024, PMID 38401548). This folder does **not** run GSEA (that leftover is `methods/gse253564_cldn4_gsea`) and does **not** re-audit MPR/PFS (`opus_geo_leftover`). TACSTD2 is a companion column only.

**Question.** Does continuous CLDN4 track PD-L1 (`CD274`) or classical MHC-I (`HLA-A`, `HLA-B`, `HLA-C`) on this leftover ICI-lung matrix? Partial Spearman residualizes ranks on `CD8A`.

**Scale.** Author FPKM → `log2(FPKM+1)`. MHC-I mean-z = gene-wise z of HLA-A/B/C only (**3/3**), not B2M/TAP. Patient/tumour is the unit (one pre-treatment column per tumour).

Reproduce: `python3 methods/gse253564_cldn4_cd274/analyze.py` (downloads the public FPKM if missing).

---

## Honest n

GEO deposits **32** pre-treatment tumours. Continuous Spearman uses every complete pair. Do **not** write the GSEA quartile n (8 vs 8) as the n for these ρ values.

| Item | n |
|---|---|
| GEO pre-treatment FPKM columns | **32** |
| Unique tumours (1 column = 1 tumour) | **32** |
| CLDN4 + CD274 + HLA-A/B/C non-NA | **32 / 32** |
| CD8A non-NA (partial covariate) | **32 / 32** |
| HLA-I genes in the mean-z | HLA-A,HLA-B,HLA-C (3/3) |
| GSEA Q4 vs Q1 (other leftover; not used here) | 8 vs 8 |
| RECIST / DCB on this page | **0** (not an ICI-response test) |

---

## One row (CLDN4, n=32)

| Cohort | n | scale | CD274 ρ (p) | HLA-A ρ (p) | HLA-B ρ (p) | HLA-C ρ (p) | MHC-I mean-z ρ (p) | CD274 \| CD8A ρ (p) | MHC-I \| CD8A ρ (p) |
|---|---:|---|---|---|---|---|---|---|---|
| **GSE253564** | **32** | log2(FPKM+1) | -0.228 (0.209) | -0.113 (0.538) | -0.181 (0.322) | +0.064 (0.729) | -0.091 (0.622) | -0.095 (0.612) | +0.274 (0.135) |

Full numeric row: `tables/one_row.tsv`. Axis-level tests: `tables/spearman.tsv`.

---

## How to read the row

Continuous CLDN4 vs **CD274** is **null** at n=32 (ρ=-0.228, p=0.209). Classical HLA-A/B/C and the HLA-A/B/C mean-z are also **null** (|ρ|≤0.181; MHC-I ρ=-0.091, p=0.622). Residualising ranks on CD8A does not create a CD274 or MHC-I hit (partial CD274 ρ=-0.095, p=0.612; partial MHC-I ρ=+0.274, p=0.135).

HLA-C after the CD8A residual is nominally positive (partial ρ=+0.397, p=0.027). Raw HLA-C is null (ρ=+0.064, p=0.729). That partial is one of five primary axes, unadjusted, and the MHC-I mean-z residual stays null. Do not upgrade it to an HLA-C claim.

The sibling GSEA leftover (`gse253564_cldn4_gsea`) already printed CLDN4 vs CD274 ρ=−0.228, p=0.209, n=32 on the same matrix. This page keeps that CD274 number and adds the **gene-level HLA-A/B/C** row plus the **CD8A partial**. Do not upgrade the GSEA MHC-I/APM set (which mixed B2M/TAP) into an HLA-gene claim. HLA-A/B/C mean-z here is -0.091.

Companion (not a CD274/HLA claim): CLDN4 vs CD8A ρ=-0.438, p=0.012; vs B2M ρ=-0.250, p=0.168; vs TACSTD2 ρ=+0.687, p=1.39e-05. CLDN4 is CD8A-low on this leftover matrix; that is a companion immune-low direction, not the CD274/HLA question.

This is pre-treatment leftover durvalumab ± SBRT bulk, not an MPR test. CLDN4 vs MPR was already null in `opus_geo_leftover`.

---

## Locked design (before ρ)

| Item | Choice |
|---|---|
| Matrix | Public author FPKM only. No FASTQ / SRA. |
| File | `GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz` |
| Scale | log2(FPKM+1); duplicate symbols collapsed to the highest-mean locus |
| Unit | one GEO pre-treatment column = one tumour |
| Honest n | complete cases with CLDN4 + CD274 + HLA-A/B/C |
| Predictor | CLDN4. TACSTD2 is a companion only. |
| HLA-I score | unweighted mean of per-gene z-scores for HLA-A, HLA-B, HLA-C (3/3). Not B2M/TAP. |
| Continuous test | two-sided Spearman |
| Partial | first-order rank residual on CD8A |
| Out of scope | GSEA; Q4 vs Q1 NES; GSE248378 post-treatment; MPR/PFS re-test; TACSTD2 splits |

---

## What this does not claim

- Not GSEA. The IFN/MHC/TJ NES leftover is a different agent and a different folder.
- Not a TACSTD2 analysis and not a dual-high gate.
- Not an MPR / PFS / recurrence test.
- Not CLDN4 vs PD-L1 protein or HLA protein. These are RNA rows.
- Not a claim that CLDN4 *induces* PD-L1. Association only. Bulk FFPE mixes epithelium and infiltrate.
- Do not pool this ρ with the four ICI-bulk CD274 rows in `methods/cldn4_cd274_ici`.
- Do not treat the HLA-C | CD8A partial (nominal p=0.027) as a pre-specified HLA-C hit. Raw HLA-C and MHC-I mean-z stay null.

---

## Files

- `FINDING.md` — this page
- `analyze.py` — download-once, score, write tables/figures
- `tables/one_row.tsv` — the one-row table
- `tables/spearman.tsv` / `n_table.tsv` / `samples.tsv` / `summary.json`
- `figures/cldn4_vs_cd274.png`, `cldn4_vs_hla_i.png`
