# GSE260770 — CLDN4-only plasma-exosome FPKM + counts

Additive **CLDN4-only** page on public [GSE260770](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE260770) (Cheng / Feng / Liang, *Signal Transduct Target Ther* 2024, PMID 38637495; NCT04026841 sintilimab on high-risk GGO in MPLC). GEO deposits **FPKM and counts**. No dual-high. TACSTD2 is not scored and no prior TACSTD2 page on this series is re-audited.

Tissue is **peripheral-blood plasma exosome RNA**, not tumor. Patient is the unit. GEO has **no patient ID and no timepoint**. Tests below are therefore **library-level** and are **not** 50 independent patients.

Reproduce: `python3 methods/gse260770_cldn4/analyze.py`.

---

## One-row table

| dataset | tissue | n (honest) | CLDN4 vs CD8A | CLDN4 vs T/NK | CLDN4 vs IFN | CLDN4 vs CD274 | CLDN4 vs response | epithelial residual | verdict |
|---|---|---|---|---|---|---|---|---|---|
| **GSE260770** | plasma exosome (not tumor) | **50 libraries**; patient ID **absent** (paper exosome subset **n=10**; trial ITT **n=36**) | ρ=**−0.10** p=**0.48** | ρ=**+0.16** p=**0.28** | ρ=**+0.18** p=**0.20** | ρ=**+0.05** p=**0.75** | Cliff δ (R−NR)=**+0.18** p=**0.074** (26 / 24 libraries) | computed as sensitivity (19/50 any EPI>0); not a tumor-purity correction | **FLOOR / NOT_TUMOR / PATIENT_ID_ABSENT** |

Primary matrix is author **log2(FPKM+1)**. Counts **log2(CPM+1)** give the same CLDN4–CD8A, CLDN4–CD274, and response numbers (CLDN4 FPKM vs counts ρ=0.999). Full numbers: `tables/primary.tsv`, `tables/all_tests.tsv`, `tables/label_inventory.tsv`. Extra figure: `figures/extra_CLDN4_panel.png`.

---

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| GEO libraries (GSM) | yes | **50** | GSM8124173–GSM8124222; titles `E-lncA411-##h` all unique |
| Patient / donor / timepoint on GEO or BioSample | **no** | **0** | characteristics are tissue, disease, agent, group only |
| Trial ITT patients (paper) | text | **36** | NCT04026841; not a GEO field |
| Paper exosome-RNA patients | text | **10** | 5 R + 5 NR matched; blood T1–T4 |
| Computable patient n from GEO | **no** | **0** | collapsing 50 libraries would invent IDs |
| GEO response label | yes | **50** | 26 “Responsed” / 24 “Non-responsed” **libraries** |
| CLDN4 / CD8A / CD274 present | yes | 50 | official symbols on both FPKM and counts |
| CLDN4 FPKM > 0 | yes | **7 / 50** | median 0; max 0.64 FPKM |
| CD8A FPKM > 0 | yes | **3 / 50** | median 0 |
| CD274 FPKM > 0 | yes | **6 / 50** | median 0 |
| Any epithelial gene FPKM > 0 | yes | **19 / 50** | EPCAM 2, KRT8 7, KRT18 9, KRT19 2, CDH1 8, KRT7 0 |

Do not write n=50 patients. Do not write n=10 as if GEO mapped the paper subset. The only computable public unit is the **library**.

---

## CLDN4 vs CD8A / T/NK / IFN / CD274 (library-level)

T/NK = mean-z of CD3D, CD3E, CD8A, NKG7, GNLY, PRF1, GZMA, GZMB (8/8 present). IFN = Ayers-like mean-z of IFNG, STAT1, CXCL9, CXCL10, IDO1, HLA-DRA (IDO1 is all-zero; the score still uses the six-gene list). Spearman on log2(FPKM+1). Bootstrap 95% CI, 2,000 resamples, seed `20260817` (floor genes make the CD8A interval unhelpful; the p-value is the claim).

| pair | n | ρ | p | counts ρ (p) | residual / partial \| EPI | p_adj |
|---|---:|---:|---:|---|---:|---:|
| CLDN4 vs **CD8A** | 50 | **−0.10** | **0.48** | −0.10 (0.48) | +0.09 | 0.54 |
| CLDN4 vs **T/NK** | 50 | **+0.16** | **0.28** | +0.15 (0.30) | +0.12 | 0.39 |
| CLDN4 vs **IFN** | 50 | **+0.18** | **0.20** | +0.10 (0.50) | +0.12 | 0.42 |
| CLDN4 vs **CD274** | 50 | **+0.05** | **0.75** | +0.05 (0.75) | +0.08 | 0.56 |

NKG7 (42/50) and HLA-DRA (45/50) are detected in exosome RNA; CLDN4, CD8A, and CD274 are not. This is not a tumor CD8 / PD-L1 test.

---

## CLDN4 vs response (if labeled)

GEO `group` is “Responsed to sintilimab” vs “Non-responsed to sintilimab” on every library. Paper ORR is 2/36 ITT lesions (total response 5/36). The GEO label is **not** remapped from the paper table.

| contrast | n | Metric | Effect | p | Note |
|---|---:|---|---:|---:|---|
| CLDN4, R vs NR | **50 (26 / 24) libraries** | Cliff δ (R−NR) | **+0.18** | **0.074** | both medians 0; 7 nonzero libraries |
| same, log2(CPM+1) | 50 (26 / 24) | Cliff δ (R−NR) | +0.18 | 0.074 | identical ranks |

A p=0.074 on a floor gene in **repeated blood libraries** is not a patient-level ICI association. It is reported so the leftover is a table, not an empty.

---

## Epithelial residual

Asked for if possible. An epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7) can be formed; 19/50 libraries have any of those genes >0. Partial Spearman after that score does **not** change the nulls (table above). This is **not** tumor-purity residual: the series is plasma exosome, EPCAM is 2/50, KRT7 is 0/50.

---

## Finding

CLDN4 is at the floor in public GSE260770 plasma-exosome FPKM and counts (7/50 libraries nonzero; median 0). Versus CD8A, T/NK, IFN, and CD274 the library-level Spearmans are null. Versus the GEO sintilimab label, Cliff δ = +0.18, p = 0.074, both group medians 0. Patient n is **not computable** from GEO. This is not a tumor-epithelium or tumor-purity page.

n=50 libraries from an unknown number of patients (paper subset 10; trial 36) is the honest public n. This page does not claim “no CLDN4–immune or CLDN4–ICI association in GGO/MPLC.” It claims **no usable patient-level tumor CLDN4 test in this GEO deposit**.

---

## Files

- `FINDING.md` — this page
- `analyze.py` — GEO download-once (FPKM + counts + series matrix), score, write tables/figures
- `tables/primary.tsv` — the one-row table
- `tables/all_tests.tsv` / `label_inventory.tsv` / `gene_detection.tsv` / `summary.json`
- `processed/GSE260770_library.tsv`
- `figures/extra_CLDN4_panel.png` (and `.pdf`) plus single-contrast PNGs
