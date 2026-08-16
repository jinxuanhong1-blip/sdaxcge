# Hunt log — GEO 2025–2026 lung IO scRNA (TACSTD2 / CLDN4)

All accessions below were returned by live NCBI E-utilities (`esearch`/`esummary` on `db=gds`) on 2026-08-16. No GSE/GSM IDs were invented.

## Queries (verbatim)

1. **IO + scRNA + lung + 2025–2026** — `results/hunt_geo2026_sc/geo_esearch_raw.json`  
   72 GSE hits.
2. **Broad 2026 lung scRNA** (IO keywords not required) — `results/hunt_geo2026_sc/geo_esearch_2026_broad.json`  
   193 GSE hits.

## Already-listed exclusion set

`already_listed.txt` (n=172) = every accession in PR #35 `results/fable_geo_2026/geo_search_candidates.tsv` plus the prior-slice classics named in that PR (GSE126044, GSE135222, GSE136961, GSE166449, GSE93157, GSE207422, GSE205335, GSE99254, GSE176021, GSE179994, GSE131907, …).

Of the 72 IO+scRNA hits, **49 were already listed** and **23 are NEW**. Of those 23, **14 have PDAT in 2026**.

## NEW 2026 IO-keyword series — why not a lung-ICI malignant×response test

| Accession | PDAT | Why not used as lung-ICI response cohort |
|-----------|------|------------------------------------------|
| GSE295599 | 2026/04/07 | Bulk RNA-seq companion; mixed solid tumors, not lung-primary. |
| GSE295600 | 2026/04/07 | scRNA of nivo+ipi+LDRT (RACIN). Cancer types on GEO samples: ovarian, prostate, colon, pancreatic, breast. **No lung-primary patient.** Patient 10 biopsy site = lung, cancer type = colon. |
| GSE316195 | 2026/07/13 | Metastatic PDAC (cemiplimab), not lung. |
| GSE325414 | 2026/06/10 | **NSCLC scRNA, NEW.** PEF ablation treat-and-resect, not ICI. No response labels. **Analyzed for TACSTD2/CLDN4 vs T/NK only.** |
| GSE276139 | 2026/04/03 | LUAD leptomeningeal-met CSF. TROP2 mentioned as CSF diagnostic marker. Not ICI-treated / no response labels. 745 MB RAW; not required for the ICI question. |
| GSE318408 / GSE318080 / GSE318867 | 2026/02/15 | SCLC STING / NK biology (human+mouse / ChIP / Visium). Not ICI-treated patients. |
| GSE306609 | 2026/02/26 | RCC CAR-T / CD70, not lung. |
| GSE311521 | 2026/03/02 | Pan-cancer cell lines. |
| GSE327497 | 2026/07/21 | PD-1 CAR-T in neuroinflammation, not lung cancer. |
| GSE317884 | 2026/07/09 | A549/DDP cell-line GNA experiment. |
| GSE255516 / GSE255517 | 2026/01 | CRISPR screens / mouse TME, not patient ICI. |
| GSE295518 | 2026/03/01 | synNotch CAR-macrophage, n=2. |
| GSE326791 | 2026/08/07 | CAR-NKT, n=2. |
| GSE301248 | 2026/07/15 | Engineered HSPC / T cells. |
| GSE335372 | 2026/06/26 | CAR-T spatial, n=2. |
| GSE320507 | 2026/03/05 | Long COVID NK cells. |
| GSE336437 | 2026/08/13 | BK virus endothelial cells. |
| GSE294051 | 2026/04/06 | SWI/SNF chromatin, not ICI. |
| GSE317816 | 2026/01/28 | Cell-confinement / chromosome loss. |

**Hunt conclusion:** as of this live search, **no NEW (not-already-listed) 2026 GEO series is a human lung-primary ICI-treated scRNA cohort with per-sample response labels and open processed matrices.**

## Listed-but-never-analyzed lung ICI scRNA used here

| Accession | PDAT | Why analyzed |
|-----------|------|----------------|
| GSE291670 | 2025/03/31 | In the PR #35 candidate TSV but **no TACSTD2/CLDN4 vs MPR / T-NK analysis exists**. 6 NSCLC, neoadjuvant anlotinib + camrelizumab, sample titles = MPR-1/2/3 and Non-MPR-1/2/3. Open 10x MTX. Malignant + immune cells present. |

## Already-listed 2026 lung ICI scRNA *not* re-analyzed here

| Accession | Reason skipped |
|-----------|----------------|
| GSE233203 | Already analyzed in PR #35 (pseudobulk TACSTD2/CLDN4 vs Response). |
| GSE243013 | Immune-only atlas; 6.6 GB MTX; no malignant cells (PR #62). |
| GSE317309 | CCL21-DC + pembrolizumab NSCLC scRNA. Listed. RAW.tar = 8.1 GB (over budget). |
| GSE337519 | Listed. n=1 patient paired pre/post; cannot test response. |
| GSE304741 | Listed. Sorted NK/T only, treatment-naive; no ICI response; no malignant cells in the public RDS. |

## Analyzed in this slice

1. **GSE291670** — malignant TACSTD2/CLDN4 vs MPR and vs T/NK fraction.
2. **GSE325414** — NEW 2026 NSCLC; malignant TACSTD2/CLDN4 vs T/NK (and vs mTLS). No ICI labels.
