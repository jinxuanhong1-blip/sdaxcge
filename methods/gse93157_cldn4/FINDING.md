# GSE93157 leftover — CLDN4 vs response / CD8 / CD274

**Public lung ICI: yes. CLDN4: absent. Tests: empty.**

Prat et al., *Cancer Res* 2017 ([PMID 28487385](https://pubmed.ncbi.nlm.nih.gov/28487385/); GEO [GSE93157](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93157)). FFPE tumor RNA on the nCounter PanCancer Immune Profiling Panel ([GPL19965](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL19965)). Patients received nivolumab or pembrolizumab. This slice is the leftover CLDN4 question on that named series: **CLDN4 vs ICI response, vs CD8, vs CD274**.

The series is public lung ICI. The gene is not on the deposited panel. No Mann–Whitney, Spearman, or AUC was run. A “non-significant” p-value would be invented.

## Honest n

Counted from the GEO series matrix (`!Sample_source_name_ch1`, `best.resp`, `drug`, `pfs`).

| Item | n | Source |
|---|---:|---|
| GEO arrays (all histologies) | **65** | GSM2445676–GSM2445740 |
| Melanoma | 25 | `source_name` |
| Head and neck | 5 | `source_name` |
| **Lung (this leftover)** | **35** | 22 non-squamous + 13 squamous |
| Lung with `best.resp` | **35** | CR 1 / PR 8 / SD 12 / PD 14 |
| Lung ORR (CR/PR vs SD/PD) | **9 vs 26** | `best.resp` |
| Lung DCB (CR/PR/SD vs PD) | **21 vs 14** | `best.resp` |
| Lung with PFS | **35** | `pfs` / `pfse` |
| Lung drug | 18 nivo / 17 pembro | `drug` |
| CLDN4-measured lung samples | **0** | gene not on panel |

Do not use the GEO characteristic `response`. Every one of the 35 lung samples is labeled `RC_RP_SD`, including the 14 with `best.resp=PD`. The usable clinical field is `best.resp` (and `pfs` / `pfse`).

## Panel check

Platform table: 784 probes (730 endogenous, 40 housekeeping, 8 negative, 6 positive). Series matrix: 765 quantified IDs. Supplementary `GSE93157_raw_data_values.txt.gz` searched the same way.

| Query | On GPL19965 | In series matrix | In raw suppl | Aliases / accessions searched |
|---|---|---|---|---|
| **CLDN4** | **no** | **no** | **no** | CLDN4, CLAUDIN4, CLAUDIN-4, CPE-R, CPER, CPETR, CPETR1, WBSCR8, NM_001305 |
| TACSTD2 (companion) | no | no | no | TACSTD2, TROP2, TROP-2, EGP1, EGP-1, GA733-1, M1S1, NM_002353 |
| Any CLDN* | no | no | no | claudin family grep |

Immune / epithelial genes that *are* on the panel (not substitutes): **CD8A**, **CD8B**, **CD274**, PDCD1, IFNG, EPCAM, CDH1.

## Assigned tests (empty)

| Test | n | Status | Why |
|---|---:|---|---|
| CLDN4 vs ICI response (`best.resp` ORR / DCB / PFS) | **0** | empty | CLDN4 not measured |
| CLDN4 vs CD8 (CD8A) | **0** | empty | CLDN4 not measured |
| CLDN4 vs CD274 | **0** | empty | CLDN4 not measured |

CD8A and CD274 can be scored against response on these 35 lung tumors. That is a panel-gene control, not the leftover question, and is not reported here as a CLDN4 result.

## Why stop

The leftover is closed as **missing measurement**, not as a negative association. Nearby panel genes were not used as stand-ins for CLDN4.

## Reproduce

```bash
python3 methods/gse93157_cldn4/check_panel.py
```

Writes `methods/gse93157_cldn4/panel_check.json`. GEO files cache under `methods/gse93157_cldn4/cache/` (git-ignored).
