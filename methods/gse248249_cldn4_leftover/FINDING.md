# GSE248249 leftover — CLDN4 vs CD274 / IFN / MHC-I / CXCL9/10

Additive **CLDN4-only** leftover on the public acquired-IO-resistance NSCLC **array** [GSE248249](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248249) (Memon et al., *Cancer Cell* 2024, [PMID 38215748](https://pubmed.ncbi.nlm.nih.gov/38215748/)). Platform [GPL23126](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL23126) Affymetrix Clariom D Human (transcript/gene version). Values = RMA from the deposited series matrix (Affymetrix Expression Console). Not NanoString.

**CLDN4 is on the panel.** Official-symbol cluster `TC0700007993.hg.1` (`NM_001305 // CLDN4`). Matrix feature count = **138,745**.

This page does **not** re-audit the given paired post-AR CLDN4 contrast (p=0.017 on n=13 pairs). It does **not** re-score TACSTD2. New cuts only: CLDN4 vs CD274 / IFN / MHC-I / CXCL9/10.

Reproduce: `python3 methods/gse248249_cldn4_leftover/analyze.py`.

---

## Honest n

GEO annotates timepoint (Pre-treatment / Post-treatment), patient, tumor site, sex. It does **not** annotate RECIST or a sensitive vs resistant label. Every patient in this molecular set has acquired resistance. The leftover unit is **one post-AR sample per patient**.

| Item | n |
|---|---:|
| Samples | **42** |
| Patients | **29** |
| Pre-treatment | 13 |
| Post-treatment (acquired resistance) | **29** |
| Patients with both timepoints | 13 |
| Sensitive vs resistant in GEO | **0** (no such field) |
| **Primary leftover n (post, 1/patient)** | **29** |
| CLDN4 / CD274 / HLA-A/B/C / CXCL9 / CXCL10 finite | 29 / 29 |

Do not write the given paired p=0.017 as the leftover n. Do not treat 13 pre vs 29 post as the CD274 n.

---

## One-row leftover table (CLDN4, post n=29)

| Cohort | n | scale | CD274 ρ (p) | IFN mean-z ρ (p) | MHC-I mean-z ρ (p) | CXCL9 ρ (p) | CXCL10 ρ (p) | CXCL9/10 mean-z ρ (p) |
|---|---:|---|---|---|---|---|---|---|
| **GSE248249 post-AR** | **29** | RMA | -0.059 (0.761) | +0.082 (0.673) | +0.018 (0.925) | -0.013 (0.945) | +0.181 (0.347) | +0.155 (0.422) |

Sensitivity (same genes; not the primary n):

| Stratum | n | CD274 ρ (p) | IFN ρ (p) | MHC-I ρ (p) | CXCL9 ρ (p) | CXCL10 ρ (p) | CXCL9/10 ρ (p) |
|---|---:|---|---|---|---|---|---|
| All samples | 42 | +0.044 (0.780) | +0.000 (0.999) | +0.040 (0.800) | -0.058 (0.716) | +0.101 (0.524) | +0.064 (0.688) |
| Pre | 13 | +0.324 (0.280) | +0.176 (0.566) | +0.165 (0.590) | +0.082 (0.789) | +0.148 (0.629) | +0.055 (0.859) |
| Post \| CD8A | 29 | +0.093 (0.639) | +0.317 (0.100) | +0.170 (0.388) | +0.223 (0.254) | +0.367 (0.055) | +0.418 (0.027) |

Full numbers: `tables/one_row.tsv`, `tables/leftover_tests.tsv`. Sample table: `tables/sample_level.tsv`.

---

## What is present

| Item | Present? | Probe / score |
|---|---|---|
| CLDN4 | yes | `TC0700007993.hg.1` |
| CD274 | yes | `TC0900006559.hg.1` |
| CXCL9 | yes | `TC0400011052.hg.1` |
| CXCL10 | yes | `TC0400011053.hg.1` |
| HLA-A / HLA-B / HLA-C | yes / yes / yes | 3/3 official-symbol clusters |
| IFN (Ayers IFNG 6-gene) | yes | IDO1, CXCL9, CXCL10, STAT1, HLA-DRA, IFNG (6/6) |
| IFNG single gene | yes | `TC1200011177.hg.1` |
| CD8A (partial covariate) | yes | `TC0200013323.hg.1` |
| TACSTD2 | ignored | leftover is CLDN4 only |
| ICI response label | no | GEO has timepoint only |

---

## Finding

On the leftover **post-AR** slice (n=29, one sample per patient):

- **CD274:** CLDN4 vs CD274 Spearman -0.059, p=0.761.
- **IFN:** CLDN4 vs Ayers IFNG 6-gene mean-z +0.082, p=0.673. Single-gene IFNG -0.005, p=0.980.
- **MHC-I:** CLDN4 vs HLA-A/B/C mean-z +0.018, p=0.925.
- **CXCL9/10:** CLDN4 vs CXCL9 -0.013, p=0.945; vs CXCL10 +0.181, p=0.347; vs CXCL9/10 mean-z +0.155, p=0.422.

The leftover is a **table, not an empty**: CD274, IFN, MHC-I, and CXCL9/10 are on the panel; the post-AR associations are **null**. Residualising ranks on CD8A is a sensitivity, not a new axis. The CXCL9/10 | CD8A partial is nominally positive (ρ=+0.418, p=0.027). Raw CXCL9/10 is null (ρ=+0.155, p=0.422). That partial is one of six leftover axes, unadjusted, and is not upgraded to a CXCL claim.

The all-sample (n=42) and pre (n=13) rows mix or shrink timepoints; they are not the leftover n.

The given paired post−pre CLDN4 contrast (p=0.017, 11/13 post lower) is **not re-tested** here. This page asks whether leftover CLDN4 tracks PD-L1 / IFN / MHC-I / CXCL9/10 on the public AR matrix, not whether CLDN4 itself drops at AR.

n=29 is small. This page does not claim a general CLDN4–immune rule in acquired-resistance NSCLC. It claims **no detectable leftover association in this GEO Clariom D slice**.

---

## Locked design (before ρ)

| Item | Choice |
|---|---|
| Matrix | Public GEO series matrix only. No CEL re-RMA. |
| Scale | Deposited RMA. No extra transform. |
| Unit | one post-AR GSM per patient (every patient has a post sample) |
| Honest n | complete cases with CLDN4 + leftover axes |
| Predictor | CLDN4. TACSTD2 ignored. |
| IFN | unweighted mean of per-gene z-scores for Ayers IFNG 6-gene (6/6) |
| MHC-I | unweighted mean of per-gene z-scores for HLA-A, HLA-B, HLA-C (3/3). Not B2M/TAP. |
| CXCL9/10 | single genes plus mean-z of CXCL9+CXCL10 |
| Continuous test | two-sided Spearman |
| Partial | first-order rank residual on CD8A (sensitivity) |
| Out of scope | re-audit of paired p=0.017; TACSTD2; sensitive vs resistant; GSEA |

---

## What this does not claim

- Not a re-audit of the given paired post-AR CLDN4 p=0.017.
- Not a TACSTD2 analysis and not a dual-high gate.
- Not sensitive vs resistant. GEO has no such labels.
- Not CLDN4 vs PD-L1 protein or HLA protein. These are array RMA rows.
- Not a claim that CLDN4 induces PD-L1 or IFN. Association only. Bulk FFPE mixes epithelium and infiltrate. Site is mixed (only 4/13 pairs are the same organ on the given page).
- Do not pool this ρ with leftover ICI-bulk CD274 rows from other series.
- Do not treat the CXCL9/10 | CD8A partial (nominal p=0.027) as a pre-specified chemokine hit. Raw CXCL9, CXCL10, and CXCL9/10 mean-z stay null.

---

## Files

- `FINDING.md` — this page
- `analyze.py` — download-once series matrix, score leftover axes, write tables
- `tables/one_row.tsv` — primary leftover row
- `tables/leftover_tests.tsv` — all strata / axes
- `tables/sample_level.tsv` / `n_table.tsv` / `panel_presence.tsv` / `summary.json`
