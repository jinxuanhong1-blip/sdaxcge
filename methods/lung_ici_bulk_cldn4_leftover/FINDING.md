# Additive CLDN4 rows — leftover OPEN lung ICI bulk

Public GEO processed matrices only. Patient is the unit. This page adds **CLDN4** (TACSTD2 as companion) to leftover OPEN lung ICI bulk series that were previously scored mainly for TACSTD2 or not scored for CLDN4 vs response.

**Not re-audited (taken as given):** GSE135222 TJ-7 ρ=−0.48; GSE190266 CLDN4 DCB p=0.055; GSE285029; GSE253564 / GSE248378.

Two axes are reported **separately** (not pooled):

1. **Barrier / ICI-poor / immune-low** — CLDN4 higher where response is worse, or CLDN4 anti-correlates with CD8A / IFN / ESTIMATE ImmuneScore.
2. **ADC+ICI candidate** — CLDN4-high with IFN or MHC-I (raw or after ImmuneScore). High target and an IFN/MHC-I program can coexist; that is a combination-therapy axis, not a contradiction of the barrier axis.

ESTIMATE = Yoshihara 2013 stromal + immune gene-list **mean z** on the author matrix (not the Affymetrix-calibrated purity transform). Residual / partial Spearman uses ImmuneScore as the covariate.

Reproduce: `python3 methods/lung_ici_bulk_cldn4_leftover/analyze.py` (GEO files cached under `/tmp/lung_ici_leftover/`).

---

## One row per cohort (CLDN4)

| Cohort | Tissue / assay | n | Endpoint | Metric | Effect | p | Note |
|---|---|---:|---|---|---:|---:|---|
| **GSE126044** | NSCLC tumor, RNA-seq counts, anti-PD-1 | 16 (5 R / 11 NR) | GEO responder vs non-responder | Cliff δ (R−NR) | **−0.53** | **0.115** | CLDN4 higher in NR (med 2.54 vs 3.77 log2CPM). TACSTD2 δ=−0.27, p=0.44. |
| **GSE136961** | NSCLC tumor, Oncomine Immune 395 | 21 (9 DCB / 12 NDB) | title D vs N (DCB); PFS days in GEO | CLDN4 not on panel | — | — | CD8A / IFNG / HLA-A/B/C present. No CLDN4–DCB test. |
| **GSE166449** | LUAD tumor TPM, pembrolizumab | 22 (7 R / 15 NR) | GEO title Responder vs nonResponder | Cliff δ (R−NR) | **+0.029** | **0.945** | med R 1.21 vs NR 1.33 log2(TPM+1). TACSTD2 δ=+0.24, p=0.41. |
| **GSE182328** | NSCLC tumor counts, ICI-treated | 44 (22 / 22) | GEO Akkermansia detectable vs not (**not** DCB/ORR/PFS/MPR) | Cliff δ (Akk+ − Akk−) | **−0.36** | **0.045** | CLDN4 higher in Akk−. No RECIST/PFS/MPR field on GEO. TACSTD2 δ=−0.14, p=0.43. |
| **GSE161537** | advanced NSCLC, HTG OBP ~2560, 2L PD-1/PD-L1 | 76 RECIST-evaluable (20 CR/PR / 56 SD/PD); 82 with PFS | GEO RECIST + PFS | CLDN4 not on panel | — | — | Panel has CLDN3 / EPCAM / CD8A / IFNG. 6 RECIST=NA. |
| **GSE162520** | early NSCLC, HTG OBP, **surgery** | 92 | no ICI DCB/MPR/PFS/ORR | not on panel; no ICI label | — | — | GEO design: surgically treated NSCLC. Surgical PFS/OS are not ICI endpoints. |
| **GSE93157 NSCLC** | nCounter PanCancer Immune 730 | 35 (9 CR/PR / 26 SD/PD) | GEO best.resp + PFS | CLDN4 not on panel | — | — | 22 non-squamous + 13 squamous. CD8A / IFNG / HLA present. |
| **GSE280232** | resectable KRAS-mutant lung, neoadjuvant ICI | — | no MPR/ORR/PFS/DCB in GEO | no processed bulk | — | — | FTP is 10x MTX/TCR RAW.tar only. Characteristics: regimen + KRAS/STK11. |
| **GSE283829** | NSCLC tumor counts | 17 CR vs PD (7 / 10); 27 with SD | GEO disease stage CR vs PD | Cliff δ (CR−PD) | **−0.029** | **0.962** | DCR (CR+SD vs PD) δ=−0.059, p=0.82. TACSTD2 CR vs PD δ=−0.086, p=0.81. |
| **GSE260770** | **peripheral blood** FPKM, sintilimab | 50 (26 R / 24 NR) | GEO Responsed vs Non-responsed | Cliff δ (R−NR) | **+0.18** | **0.074** | CLDN4 nonzero 7/50; median 0. Not tumor epithelium. |
| **GSE293591 lung** | LUAD+LUSC TPM (20+14) | 34 | no ICI label; CLDN4 vs IFN | Spearman ρ | **−0.22** | **0.206** | vs MHC-I ρ=+0.21, p=0.23. After ImmuneScore, vs MHC-I partial ρ=**+0.46**, p=**0.007**. |

Full tests: `tables/cohort_primary.tsv`, `tables/all_tests.tsv`, `tables/inventory.tsv`. Patient tables: `processed/`.

---

## Extra figures (direction-matched)

`figures/extra_CLDN4_two_axes.png` — four panels, two axes, not a meta-analysis.

| Panel | Contrast | n | Result | p |
|---|---|---:|---|---:|
| A | GSE126044 CLDN4, NR vs R | 11 vs 5 | δ=−0.53 (higher in NR) | 0.115 |
| B | GSE182328 CLDN4 vs CD8A | 44 | ρ=**−0.40** | **0.0067** |
| C | GSE283829 CLDN4 vs ESTIMATE ImmuneScore | 27 | ρ=**−0.57** | **0.0020** |
| D | GSE293591 lung CLDN4 vs MHC-I **after ImmuneScore** | 34 | partial ρ=**+0.46** | **0.0075** |

Single-cohort plots used to build the extra page: `GSE126044_CLDN4_ORR.png`, `GSE182328_CLDN4_vs_CD8A.png`, `GSE283829_CLDN4_vs_IFN.png` / `_vs_MHC.png`, `GSE293591_lung_CLDN4_vs_MHC.png`.

---

## Axis 1 — barrier / immune-low / ICI-poor (tumor, CLDN4 present)

| Cohort | Test | n | ρ or δ | p |
|---|---|---:|---:|---:|
| GSE126044 | CLDN4 vs ORR | 5 vs 11 | δ=−0.53 | 0.115 |
| GSE126044 | CLDN4 vs ImmuneScore | 16 | ρ=−0.57 | 0.022 |
| GSE126044 | CLDN4 vs CD8A / IFN / MHC-I | 16 | ρ=−0.46 / −0.45 / −0.44 | 0.076 / 0.080 / 0.087 |
| GSE126044 | same three, partial \| ImmuneScore | 16 | ρ=+0.17 / −0.03 / +0.19 | all p>0.48 |
| GSE182328 | CLDN4 vs Akkermansia (GEO surrogate) | 22 vs 22 | δ=−0.36 | 0.045 |
| GSE182328 | CLDN4 vs CD8A / IFN / ImmuneScore | 44 | ρ=−0.40 / −0.30 / −0.32 | 0.0067 / 0.045 / 0.036 |
| GSE182328 | CLDN4 vs CD8A, partial \| ImmuneScore | 44 | ρ=−0.27 | 0.081 |
| GSE283829 | CLDN4 vs ImmuneScore | 27 | ρ=−0.57 | 0.0020 |
| GSE283829 | CLDN4 vs IFN / CD8A | 27 | ρ=−0.36 / −0.35 | 0.064 / 0.078 |
| GSE283829 | CLDN4 CR vs PD | 7 vs 10 | δ=−0.029 | 0.96 |
| GSE166449 | CLDN4 vs ORR | 7 vs 15 | δ=+0.029 | 0.95 |
| GSE166449 | CLDN4 vs CD8A / IFN | 22 | ρ=−0.25 / −0.21 | 0.27 / 0.35 |

Raw CLDN4 tracks **immune-low** in the three whole-transcriptome tumor ICI series with an ESTIMATE ImmuneScore (GSE126044, GSE182328, GSE283829). After ImmuneScore, the CD8A/IFN/MHC-I partials in those ICI series sit near zero (purity/infiltrate accounts for most of the raw anti-correlation). GSE182328’s Akkermansia contrast is the only leftover **GEO-deposited grouping** with a nominal CLDN4 p<0.05; it is not RECIST/DCB/PFS/MPR.

---

## Axis 2 — CLDN4-high IFN / MHC-I (ADC+ICI candidates)

| Cohort | Test | n | ρ | p |
|---|---|---:|---:|---:|
| GSE293591 lung | CLDN4 vs MHC-I, partial \| ImmuneScore | 34 | **+0.46** | **0.0075** |
| GSE293591 lung | CLDN4 vs MHC-I (raw) | 34 | +0.21 | 0.23 |
| GSE293591 lung | CLDN4 vs IFN (raw) | 34 | −0.22 | 0.21 |
| GSE166449 | TACSTD2 vs MHC-I (companion) | 22 | +0.37 | 0.092 |
| GSE126044 / GSE182328 / GSE283829 | CLDN4 vs MHC-I, partial \| ImmuneScore | 16 / 44 / 27 | +0.19 / −0.14 / +0.04 | 0.49 / 0.36 / 0.84 |

The leftover ICI tumor matrices do not show a raw CLDN4-high IFN program. The **positive** leftover MHC-I residual is in GSE293591 lung (no ICI labels). That is an antigen-presentation / ADC-target coexistence row, not an ICI-response row.

---

## Companion TACSTD2 (same patients)

| Cohort | Endpoint | n | δ or ρ | p |
|---|---|---:|---:|---:|
| GSE126044 | ORR | 5 vs 11 | δ=−0.27 | 0.44 |
| GSE166449 | ORR | 7 vs 15 | δ=+0.24 | 0.41 |
| GSE182328 | Akkermansia | 22 vs 22 | δ=−0.14 | 0.43 |
| GSE283829 | CR vs PD | 7 vs 10 | δ=−0.086 | 0.81 |
| GSE283829 | vs ImmuneScore | 27 | ρ=−0.41 | 0.032 |

---

## Honest exclusions (still one table row)

- **GSE136961, GSE161537, GSE93157 NSCLC:** ICI labels exist; CLDN4/TACSTD2 are absent from the targeted immune panels.
- **GSE162520:** HTG panel without CLDN4; GEO is a **surgical** early-NSCLC series, not ICI-treated.
- **GSE280232:** no public processed bulk; no MPR/ORR/PFS/DCB in GEO characteristics.
- **GSE260770:** blood, CLDN4 at floor. Reported so the leftover is closed; not a tumor-barrier test.
- **GSE293591:** no ICI field. Lung subset used only for the immune / MHC-I axis.

No statistics were taken from the named “given” series. No labels were read out of PDFs when GEO lacked the field.

---

## Files

- `FINDING.md` — this page
- `analyze.py` — download-once, score, write tables/figures
- `estimate_gene_sets.csv` — Yoshihara stromal/immune lists (tidyestimate)
- `hgnc_symbol_ensembl.tsv` — symbol ↔ Ensembl for GSE283829
- `tables/cohort_primary.tsv` — the one-row-per-cohort table
- `tables/all_tests.tsv` / `inventory.tsv` / `summary.json`
- `processed/*_patient.tsv`
- `figures/extra_CLDN4_two_axes.png` and per-cohort PNGs
