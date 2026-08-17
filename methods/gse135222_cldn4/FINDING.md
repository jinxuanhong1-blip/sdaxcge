# GSE135222 leftover — CLDN4 vs CD8A / CD274 / response

Additive **CLDN4-only** leftover on the public NSCLC ICI RNA-seq slice [GSE135222](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135222) (Jung et al., *Nat Commun* 2019, PMID 31537801; RNA-seq subset of the SMC cohort also in Kim et al., *Clin Epigenetics* 2020, PMID 32762727).

**Honest n = 27 patients** (one pre-ICI tumor RNA-seq sample each). CLDN4, CD8A, and CD274 are all on the GEO RSEM TPM matrix. ICI response is present as author-assigned durable clinical benefit (DCB), not as a GEO characteristic.

This page does **not** re-score TACSTD2. Patient is the unit. Reproduce: `python3 methods/gse135222_cldn4/analyze.py`.

---

## One-row leftover table

| Contrast | n | Endpoint / genes | Metric | Effect | p | Note |
|---|---:|---|---|---:|---:|---|
| **CLDN4 vs CD8A** | **27** | both present (log2 TPM+1) | Spearman ρ | **+0.12** | **0.54** | no association |
| **CLDN4 vs CD274** | **27** | both present (log2 TPM+1) | Spearman ρ | **−0.21** | **0.30** | no association |
| **CLDN4 vs response** | **27 (8 DCB / 19 NDB)** | Kim 2020 Table 1 benefit Y/N | Cliff δ (DCB−NDB) | **+0.039** | **0.90** | AUC 0.52; med 7.45 vs 7.23 |

Full numbers: `tables/primary.tsv`, `tables/leftover_tests.tsv`, `tables/inventory.tsv`. Patient table: `processed/GSE135222_patient.tsv`. Figure: `figures/CLDN4_vs_CD8A_CD274_DCB.png`.

---

## What is present

| Item | Present? | Honest n / source |
|---|---|---|
| CLDN4 | yes (`ENSG00000189143`) | 27 / 27 |
| CD8A | yes (`ENSG00000153563`) | 27 / 27 |
| CD274 | yes (`ENSG00000120217`) | 27 / 27 |
| ICI response | yes | **8 DCB / 19 NDB** author benefit (Kim 2020 Table 1). GEO itself has PFS event + `pfs.time` only. |
| TACSTD2 | ignored | leftover is CLDN4 only |

Jung/Kim define DCB as RECIST v1.1 PR or SD lasting >6 months. A GEO-only PFS ≥ 183 d rule is **7 / 20** and disagrees on one patient (NSCLC1708: PFS 168 d, event, author DCB=Y). Sensitivity on that PFS split: Cliff δ = +0.11, p = 0.69. Same null.

---

## Finding

CLDN4 does not track CD8A, CD274, or author DCB in this public n=27 slice. The leftover is a **table, not an empty**: the genes and the response label are present; the associations are null.

n=27 (8 vs 19) is small. This page does not claim “no CLDN4–immune or CLDN4–ICI association in NSCLC.” It claims **no detectable association in this GEO RNA-seq slice**.

---

## Files

- `FINDING.md` — this page
- `analyze.py` — GEO download-once, score, write tables/figure
- `tables/primary.tsv` — the three leftover rows
- `tables/leftover_tests.tsv` — primary + PFS≥183 sensitivity
- `tables/inventory.tsv` / `summary.json`
- `processed/GSE135222_patient.tsv`
- `figures/CLDN4_vs_CD8A_CD274_DCB.png`
