# GSE165252 leftover — CLDN4 vs immune / CD274

**EMPTY. Not public lung.**

**Slice:** `methods/gse165252_cldn4/`
**Date opened:** 2026-08-17
**Scope:** additive leftover only **if** [GSE165252](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165252) is a public **lung** RNA series. Then score CLDN4 vs immune / CD274. Honest n or empty.
**This slice is a gate, not a re-score.** No CLDN4–CD274 ρ was computed. No esophageal numbers are reported as lung.

---

## Verdict

**GSE165252 is public, but it is not lung.** Honest lung n = **0**. The leftover table is empty.

The deposited series is the phase II **PERFECT** trial: atezolizumab plus neoadjuvant chemoradiotherapy in **resectable esophageal adenocarcinoma** (van den Ende et al., *Clin Cancer Res* 2021, PMID [33504550](https://pubmed.ncbi.nlm.nih.gov/33504550/)). GEO overall design: “RNA sequencing of **esophageal** tumor biopsies at baseline … on-treatment (3rd week) and at resection.”

The gate “if public lung” fails. CLDN4 vs CD274 / ImmuneScore / CD8A is **not run**.

---

## Gate

| item | value |
|---|---|
| accession | [GSE165252](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165252) |
| public | yes (released 2021-02-04) |
| processed matrix on GEO | yes (`GSE165252_norm.cnt_PERFECT.txt.gz`, `GSE165252_vst_PERFECT.txt.gz`) |
| organism | Homo sapiens |
| platform | GPL20301 Illumina HiSeq 4000 |
| tissue / disease | **esophageal adenocarcinoma**, not NSCLC / LUAD / LUSC |
| lung samples on the series matrix | **0** |
| honest n for this leftover | **0** |
| CLDN4 vs CD274 / immune | **not run** |

`tables/one_row.tsv` is header plus one empty lung row. `tables/inventory.tsv` records the GEO identity used to fail the gate.

---

## What the series actually is (not the analysis n)

Opened from `GSE165252_series_matrix.txt.gz` (GEO FTP, 2026-08-17). These counts describe the esophageal record. They are **not** a lung analysis n.

| item | n | rule |
|---|---:|---|
| arrays / GSM on the series matrix | 77 | `!Sample_geo_accession` |
| unique patients (title `sample_<id>_`) | 40 | patients 1–4, 6–12, 15, 17–26, 28–45 |
| baseline / on-treatment / resection | 35 / 31 / 11 | sample title token |
| endoscopic vs resection biopsy | 66 / 11 | `!Sample_source_name_ch1` |
| sample-level response | 23 R / 48 NR / 6 NA | `response:` characteristic |
| patients with a baseline biopsy | 35 | 20 NR / 12 R / 3 NA; 5 patients have no baseline GSM |
| ICI labels | yes | atezolizumab + nCRT; not a lung ICI arm |
| lung / LUAD / LUSC / NSCLC token in title, source, or characteristics | **0** | none |

Patient is not interchangeable with GSM. Do not write n=77 as a patient n, and do not write n=40 as a lung n.

---

## Why some LUAD papers list this accession

A few LUAD immunotherapy papers treat GSE165252 as a second bulk ICI cohort next to GSE126044 (example: ferroptosis / EMT deconvolution, *Front. Cell Dev. Biol.* 2022). That is a **literature mis-use of an esophageal series**, not a GEO lung deposit.

This leftover does **not** follow that mis-label. It does not score CLDN4 on esophagus and then report the ρ as public lung. It does not pool GSE165252 with GSE126044, GSE218989, or other OPEN lung ICI leftovers (PR #296).

---

## Empty leftover table

| dataset | public lung? | n_lung | CLDN4–CD274 ρ (p) | CLDN4–CD8A ρ (p) | CLDN4–ImmuneScore ρ (p) | note |
|---|---|---:|---|---|---|---|
| GSE165252 | **no** | **0** | — | — | — | esophageal EAC PERFECT trial; gate fail |

No Spearman, no ESTIMATE residual, no Q4 vs Q1, no ICI ORR test.

---

## What is not done

- No download of the PERFECT count / VST matrix for gene-level tests.
- No CLDN4 / TACSTD2 / CD274 / HLA / CD8A correlation on esophagus.
- No re-label of PERFECT responders as lung ICI.
- No recut of already-scored lung leftovers (GSE126044, GSE218989, GSE166449, GSE182328, GSE253564, GSE31210).

---

## Files

- `FINDING.md` — this write-up
- `tables/inventory.tsv` — GEO identity + gate
- `tables/one_row.tsv` — empty lung leftover row
