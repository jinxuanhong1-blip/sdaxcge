# Hunt: public TROP2+CLDN4 lung IHC tables

**Verdict: NONE**

No public Figshare, Zenodo, or journal-supplement table gives **paired TROP2 + CLDN4** H-scores or a **co-localization** metric in human lung.

Do not treat the near-misses below as substitutes. They fail at least one inclusion rule.

## Inclusion rules

A hit must satisfy all of:

1. Public file on Figshare, Zenodo, or a journal supplement (CSV/XLSX/TSV or equivalent table).
2. Human lung tissue (normal, LUAD, LUSC, NSCLC, or other primary lung tumor).
3. Protein assay: IHC, multiplex IF, or dual stain — not RNA-seq / scRNA / proteomics MS alone.
4. Both markers on the **same cases**: TROP2 (TACSTD2 / Trop-2) **and** CLDN4 (claudin-4).
5. Numeric table: per-case or per-core H-score (or equivalent intensity×extent score), or a co-localization / dual-positive fraction.

## Explicitly out of scope (private)

| Label | Why excluded |
| --- | --- |
| user 52 WSI | Private whole-slide set; not a public deposit. |
| Zhejiang 25-pair | Private paired cohort; not a public Figshare/Zenodo/supplement table. |

A published **TROP2-only** 25-pair (carcinomatous vs sarcomatous components) from Nanchang University PSC (PMID 40739577) is **not** the Zhejiang 25-pair and is **not** TROP2+CLDN4.

## What was searched (2026-08-16)

### Repository APIs

- **Zenodo** `GET /api/records`: `TROP2 AND CLDN4`, `TACSTD2 AND CLDN4`, `TROP2 AND lung AND "H-score"`, `CLDN4 AND lung AND "H-score"`, `TROP2 AND CLDN4 AND immunohistochemistry`, `TACSTD2 AND claudin-4 AND IHC`, plus broader `TROP2+CLDN4`, `TACSTD2`, `TROP2 lung`, `CLDN4 lung immunohistochemistry`.
- **Figshare** `POST /v2/articles/search`: `TROP2 CLDN4`, `TACSTD2 CLDN4`, `TROP2 CLDN4 lung`, `TROP2 H-score lung`, `CLDN4 H-score lung`, `claudin-4 TROP2 immunohistochemistry`, `Trop-2 claudin-4 NSCLC`, `CLDN4 lung immunohistochemistry`, `claudin-4 NSCLC H-score`, `TROP2 NSCLC supplementary H-score`, `Moldvay claudin lung`.

### Literature / supplement index

- Europe PMC: title intersection `TROP2` ∩ `CLDN4` = **0**. Abstract intersection of both proteins + lung + IHC = **0**.
- Europe PMC abstract `(TROP2 OR Trop-2 OR TACSTD2) AND (CLDN4 OR "claudin-4")` returned 13 records; none are lung dual-IHC H-score tables.
- Web / publisher pages for Figshare collections, Zenodo records, PLOS / Frontiers / AACR / Springer supplements.

### Other stores checked only to rule out mis-attribution

Dryad, OSF, Mendeley Data, Human Protein Atlas, CPTAC-LUAD/TCIA. These are **outside** the requested Figshare/Zenodo/supplement lane and still do not provide paired TROP2+CLDN4 H-score tables.

## Near-misses (do not use as the requested table)

See `near_misses.tsv`. Summary:

| Item | What it actually is | Why it fails |
| --- | --- | --- |
| Zenodo 10.5281/zenodo.18543127 (and 18494664) | NSCLC TMA; AI H-scores for **TROP-2** (membrane/cytoplasm) and **cMET** | No CLDN4. Files **restricted**. |
| Figshare 10.6084/m9.figshare.19665135 | Dum et al. TROP2 TMA on 18,563 tumors; lung **category %** in the paper | TROP2 only; supplement is IHC **figures**, not a patient-level H-score table; no CLDN4. |
| Ahmed et al. PLOS ONE 10.1371/journal.pone.0321555 | Trop-2 membrane H-scores in NSCLC (summaries / plots) | TROP2 only; no CLDN4; patient-level table not deposited (data on request to Gilead). |
| Liu et al. J Transl Med PMID 40739577 | Trop-2 H-score in 35 PSC; 25 CaC/SaC pairs | TROP2 only; Nanchang, not Zhejiang; no CLDN4 table on Figshare/Zenodo. |
| Báthori et al. Pathol Oncol Res 10.3389/pore.2023.1611328 | CLDN1–5,7,18 H-scores in 35 rare lung ACC/MEC | CLDN4 present; **no TROP2**. |
| Moldvay / Jung-type CLDN NSCLC IHC | CLDN family scores in lung subtypes | No TROP2; no public paired table found on Figshare/Zenodo. |
| Nakashima et al. Am J Pathol 2010 | TACSTD2 binds CLDN1/7; **CLDN4 did not co-IP** | Cell-line biochemistry, not lung H-score table. |
| Zhao et al. JITC 2026 / bioRxiv 2024 | TROP2 / **claudin-7** immune-exclusion program | Breast, CLDN7 not CLDN4, not a lung H-score table. |
| Sun et al. npj Precis Oncol 10.1038/s41698-026-01523-w | TACSTD2 **RNA** co-expression with CLDN4/CLDN7 | Transcriptomic, not IHC H-score / co-localization. |
| Chen et al. Clin Transl Med 10.1002/ctm2.1649 | scRNA: CLDN4+ MPE cells; TACSTD2 correlation | RNA, not IHC H-score table. |
| Human Protein Atlas TACSTD2 / CLDN4 | Separate lung IHC images and ordinal stains | Not Figshare/Zenodo/supplement; not paired H-scores; not co-localization. |
| CPTAC-LUAD | MS protein abundance + H&E WSI | Not TROP2+CLDN4 IHC H-scores. |

## Files in this folder

| File | Contents |
| --- | --- |
| `VERDICT.txt` | One-line result: NONE |
| `README.md` | This note |
| `search_log.md` | Query-by-query log |
| `near_misses.tsv` | Structured near-misses |
| `excluded_private.tsv` | Private sets named in the request |
