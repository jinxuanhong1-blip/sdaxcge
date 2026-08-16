# Hunt: the durvalumab / anti–PD-L1 NSCLC RNA behind "TACSTD2–immune ρ≈−0.65 (purity-partial ≈−0.46)"

**Goal.** Find the *public* durvalumab / PACIFIC / anti–PD-L1 NSCLC RNA dataset behind a
reported TACSTD2 (TROP2) vs immune Spearman ρ≈−0.65 (purity-partial ≈−0.46); also OncoSG.
Download only open, processed data (< 2 GB), and **repeat the purity-partial Spearman**.

**Bottom line (honest).** The reported numbers are reproducible.
- In the open **durvalumab (anti–PD-L1)** dataset **GSE248378**, TACSTD2 vs ESTIMATE
  ImmuneScore gives Spearman **ρ = −0.649** and the ESTIMATE-purity-partial gives
  **ρ = −0.464** (n = 29). This matches the reported ≈−0.65 / ≈−0.46 almost exactly.
- In **OncoSG** LUAD (n = 169), using the cohort's **own measured tumour purity** and IMSIG
  immune signatures, TACSTD2 vs immune is **ρ = −0.504** and the purity-partial is
  **ρ = −0.467** — an independent corroboration that the association is not merely a purity artefact.

No accession IDs were invented; every ID below was verified against its live source.

---

## What is actually public vs. not (the "hunt")

The specific anti–PD-L1 NSCLC cohort in which TROP2/TACSTD2 was linked to immune exclusion
and ICI resistance is **POPLAR + OAK (atezolizumab)** — Bessede et al., *Clin Cancer Res*
2024 (PMID via DOI 10.1158/1078-0432.CCR-23-2566). That RNA-seq is **controlled access**, so
it is **not** an open <2 GB download:

| Candidate | Drug / class | RNA data status | Verdict |
|---|---|---|---|
| **POPLAR + OAK** (Bessede 2024) | atezolizumab (anti–PD-L1) | EGA **controlled**: `EGAS00001005013` (reprocessed), source FASTQ `EGAC00001002120` | Not openly downloadable |
| **PACIFIC** (NCT02125461) | durvalumab (anti–PD-L1) | Tissue optional; only PD-L1 IHC. Access via AstraZeneca request portal only | No public transcriptome |
| **SUBMARINE / WJOG11518L** (JTO 2023) | durvalumab (PACIFIC regimen) | Transcriptome generated; **no public GEO/EGA accession found** | Not openly downloadable |
| **GSE248378** | **durvalumab (anti–PD-L1)** | **Open** GEO bulk RNA-seq FPKM (resected tumours) | ✅ Used |
| **luad_oncosg_2020 (OncoSG)** | none (LUAD resection cohort) | **Open** on cBioPortal; ships purity + IMSIG | ✅ Used |
| GSE207422 | anti–**PD-1** + chemo (neoadjuvant) | Open, but **not** durvalumab/anti–PD-L1 | Not used (misattributed as "durvalumab" by a web summary; verified against GEO record) |

So the openly downloadable durvalumab/anti–PD-L1 NSCLC RNA that reproduces the claim is
**GSE248378**, and **OncoSG** was analysed as requested.

---

## Datasets used (all open, processed, total ≈ 4 MB)

1. **GSE248378** — "Neoadjuvant durvalumab ± radiation in stage I–III NSCLC" (anti–PD-L1),
   bulk RNA-seq **FPKM** of resected (post-treatment) tumours, 15,166 genes × 29 samples.
   File `GSE248378_Durva_Post_FPKMs.txt.gz` (~1 MB) from GEO.
2. **luad_oncosg_2020** — OncoSG LUAD (Chen et al., *Nat Genet* 2020), cBioPortal.
   Per-sample **PURITY** and **IMSIG** immune-cell signatures (`data_clinical_sample.txt`),
   plus TACSTD2 (EntrezID 4070) mRNA z-scores via the cBioPortal API. n = 169 with
   TACSTD2 + IMSIG + PURITY all present.

Exact URLs, sizes, and SHA-256 are in [`data/DOWNLOADS_manifest.tsv`](data/DOWNLOADS_manifest.tsv).

## Methods

- **Immune signature (GSE248378):** ESTIMATE ImmuneScore via ssGSEA, re-implemented in Python
  directly from the ESTIMATE R package (`estimate_1.0.13`, R/estimateScore.R): filter to the
  10,412 common genes, per-sample average-rank normalisation ×10000/Ng, ssGSEA walk with weight
  exponent 0.25, ES = Σ(Fn−F0). Gene sets = ESTIMATE `SI_geneset.gmt` (141 immune, 141 stromal).
- **Purity (GSE248378):** ESTIMATEScore → Affymetrix purity formula
  `cos(0.6049872018 + 0.0001467884·ESTIMATEScore)` used as a **proxy** (see caveats).
- **Immune signature + purity (OncoSG):** the cohort's own published per-sample IMSIG immune
  signatures (composite = mean of z-scored B/T/NK/macrophage/monocyte/neutrophil/plasma/interferon)
  and its own measured tumour **PURITY** — a purity estimate **independent** of the immune score.
- **Purity-partial Spearman:** rank-based first-order partial correlation,
  `r_xy·z = (r_xy − r_xz·r_yz)/√((1−r_xz²)(1−r_yz²))`, with two-sided t-test on df = n−3.

## Results

**GSE248378 (durvalumab, anti–PD-L1; n = 29)**

| Comparison | Spearman ρ | p |
|---|---|---|
| TACSTD2 vs ImmuneScore | **−0.649** | 1.4e-04 |
| TACSTD2 vs ESTIMATEScore | −0.532 | 3.0e-03 |
| TACSTD2 vs ESTIMATE purity (proxy) | +0.532 | 3.0e-03 |
| ImmuneScore vs ESTIMATE purity (proxy) | −0.903 | 2.1e-11 |
| **TACSTD2 vs ImmuneScore, purity-partial** | **−0.464** | 1.3e-02 |

**OncoSG luad_oncosg_2020 (LUAD; n = 169; real purity)**

| Comparison | Spearman ρ | p |
|---|---|---|
| TACSTD2 vs IMSIG immune composite | −0.504 | 2.9e-12 |
| TACSTD2 vs PURITY | +0.238 | 1.8e-03 |
| immune composite vs PURITY | −0.619 | — |
| **TACSTD2 vs immune, purity-partial** | **−0.467** | 1.7e-10 |

Per-IMSIG-signature (OncoSG), Spearman / purity-partial: T-cells −0.455/−0.404,
NK −0.471/−0.421, macrophages −0.462/−0.408, neutrophils −0.503/−0.456,
monocytes −0.440/−0.381, B-cells −0.322/−0.240, plasma −0.329/−0.285, interferon −0.130/−0.063 (ns).

## Interpretation & honest caveats

- The GSE248378 result reproduces the reported ρ≈−0.65 / partial≈−0.46 to two decimals,
  strongly suggesting the original claim came from **ESTIMATE ImmuneScore + ESTIMATE-purity-partial
  on this (or a very similar) durvalumab bulk RNA-seq cohort**.
- **Caveat (GSE248378 purity):** ESTIMATE purity is Affymetrix-calibrated and, by construction,
  derived from the immune+stromal scores (here ImmuneScore vs purity-proxy ρ = −0.90). So its
  "purity-partial" partly conditions on a function of the immune score itself, and n = 29 is small
  (post-treatment resected tumours, not the reported discovery cohort). Treat as a reproduction of
  the *method*, not an independent biological purity control.
- **OncoSG is the clean control:** it uses an **independent, measured** purity. TACSTD2 correlates
  *positively* with purity (+0.24; TACSTD2 is tumour-cell expressed) yet the TACSTD2–immune
  association stays strongly negative after removing purity (−0.504 → −0.467). The negative
  TACSTD2–immune relationship is therefore **not** an artefact of tumour purity.

## Reproduce

```bash
pip install pandas scipy numpy
cd results/hunt_durva
python3 analyze.py        # reads data/, writes *_scores.csv, *_merged_table.csv, results_summary.json
```

## Files
- `analyze.py`, `estimate_ssgsea.py` — analysis code (ESTIMATE ssGSEA + partial Spearman).
- `results_summary.json` — all ρ / p values.
- `GSE248378_estimate_scores.csv`, `GSE248378_merged_table.csv`, `OncoSG_merged_table.csv`.
- `data/` — downloaded inputs + `DOWNLOADS_manifest.tsv` (URLs, sizes, SHA-256).
