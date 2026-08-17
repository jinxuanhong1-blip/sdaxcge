# FINDING — GSE161537 leftover: CLDN4-only on NIVOBIO HTG EdgeSeq

**Verdict: panel-missing.** GSE161537 is public 2L PD-1/PD-L1 NSCLC (NIVOBIO; n=82 patients). **CLDN4 is not on the deposited HTG EdgeSeq Oncology Biomarker Panel.** Therefore CLDN4 vs DCB / ORR / PFS and CLDN4 vs CD8A / IFN / CD274 are **n = 0**. No values were invented. CLDN3 is on the panel and was **not** used as a proxy.

This page is additive CLDN4-only. Prior TACSTD2 leftover pages were not opened or re-scored.

---

## Honest n

Patient is the unit. One GSM = one `patient id` = one matrix column. Do not write n>82.

| Item | n | Note |
| --- | ---: | --- |
| Public GEO series | **1** | [GSE161537](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161537), public 2022-07-07; PMID [36111282](https://pubmed.ncbi.nlm.nih.gov/36111282/), [36038492](https://pubmed.ncbi.nlm.nih.gov/36038492/) |
| GSM / titles | **82 / 82** | GSM4909684–GSM4909765; title `Targeted RNA of NSCLC patient number *` |
| Unique `patient id` | **82** | no duplicate patients |
| Matrix columns (`Patient_number`) | **82** | `GSE161537_nivobio_log2cpm.csv.gz`; IDs match GEO 1:1 |
| Deposited genes | **2559** | HTG EdgeSeq OBP; unique symbols = 2559 |
| Histology (deposited) | 55 AD / 17 SC / 3 CNEGC / 3 CSARC / 3 NOS / 1 ADSC | `histology:` |
| Sex | 53 M / 29 F | `Sex:` |
| Stage | 54 IV / 10 IIIB / 9 IIIA / 4 IIB / 2 IA / 2 IB / 1 IIA | `Stage:` |
| Immunotherapy line (characteristic) | 47 `1` / 17 `2` / 13 `3` / 2 `0` / 1 each of `4`,`5`,`6` | series design is 2L PD-1/PD-L1; this field is the deposited IO-line token, not re-coded |
| HOT phenotype | 44 HOT / 38 COLD | author score; not a CLDN4 test |
| RECIST (`best response on immunotherapy (recist)`) | **76 labeled / 6 NA** | CR 1 / PR 19 / SD 22 / PD 34 / NA 6 |
| ORR if labeled (CR/PR vs SD/PD) | **20 vs 56** | among 76 with RECIST; 6 NA excluded |
| PFS months | **82** | `pfs (month)`; **no** `pfs (event)` field |
| OS months / events | **82 / 61** | `os (month)` / `os (event)` 61 events, 21 censored |
| DCB / NDB field | **0** | not deposited; not recoded from PFS |
| **CLDN4 measured** | **0** | also absent: TACSTD2, CLDN7, CLDN18, CD8B |
| CLDN4 vs DCB | **0** | empty (no CLDN4; no DCB label) |
| CLDN4 vs ORR | **0** | empty |
| CLDN4 vs PFS | **0** | empty |
| CLDN4 vs CD8A / IFN / CD274 | **0** | empty |
| CD8A / IFNG / CD274 on panel | 82 / 82 / 82 | present; **not substitutes** |
| CLDN3 on panel | 82 | present; **not a CLDN4 proxy** |

Do not treat immunotherapy-line `1` as a contradiction of the series 2L design. Do not invent DCB as PFS ≥ 6 months. Do not write n=2559 as a patient n.

---

## Panel check

Assay (GEO overall design and extract protocol): HTG EdgeSeq technology, **Oncology Biomarker Panel (OBP)**, FFPE pre-IO biopsy or archival resection, Illumina NextSeq 500 ([GPL18573](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL18573) is the sequencer, not a gene table). The gene list that matters is the deposited matrix.

Source: GEO FTP `GSE161537_nivobio_log2cpm.csv.gz` (semicolon CSV; 2559 rows × 82 patients). Aliases searched as exact symbol tokens.

| Query | On deposited OBP matrix | IDs | Aliases / tokens searched |
| --- | --- | --- | --- |
| **CLDN4** | **no** | — | CLDN4, CLAUDIN4, CLAUDIN-4, CPE-R, CPER, CPETR, CPETR1, WBSCR8, HCPER |
| TACSTD2 / TROP2 | no | — | TACSTD2, TROP2, TROP-2, EGP1, EGP-1, GA733-1, M1S1 |
| Any CLDN* | **CLDN3 only** | `CLDN3` | family grep `^CLDN` |
| CD8A | yes | `CD8A` | — |
| CD8B | no | — | — |
| CD274 | yes | `CD274` | — |
| IFNG | yes | `IFNG` | — |
| IFNGR1 | yes | `IFNGR1` | — |
| PDCD1 | yes | `PDCD1` | — |
| PDCD1LG2 | yes | `PDCD1LG2` | — |
| CXCL9 / CXCL10 | yes / yes | `CXCL9`, `CXCL10` | — |

Vendor OBP gene-list PDF (Omixys HTG EdgeSeq OBP list) also lists **CLDN3** and does not list **CLDN4**. The stop rule uses the **deposited NIVOBIO matrix**, not the brochure.

Re-probe:

```bash
python3 methods/gse161537_cldn4/check_panel.py
```

Writes `panel_check.json`, `presence_table.tsv`, `panel_symbols.tsv`, `tables/inventory.tsv`, `tables/one_row.tsv`. GEO files cache under `methods/gse161537_cldn4/cache/` (git-ignored).

---

## Assigned tests (empty)

| Test | n | Status | Why |
| --- | ---: | --- | --- |
| CLDN4 vs DCB | **0** | empty | CLDN4 not measured; DCB not a GEO field |
| CLDN4 vs ORR (RECIST CR/PR vs SD/PD) | **0** | empty | CLDN4 not measured |
| CLDN4 vs PFS (`pfs (month)`) | **0** | empty | CLDN4 not measured |
| CLDN4 vs CD8A | **0** | empty | CLDN4 not measured |
| CLDN4 vs IFN (`IFNG`) | **0** | empty | CLDN4 not measured |
| CLDN4 vs CD274 | **0** | empty | CLDN4 not measured |

CD8A, IFNG, and CD274 can be scored against RECIST / PFS on these 82 patients. That is a panel-gene control, not the leftover question, and is **not** reported here as a CLDN4 result.

**No extra figures.** A CLDN4 boxplot, forest, or Spearman scatter would require a CLDN4 vector. None was drawn. CLDN3 was not plotted as a stand-in.

---

## Why stop

The leftover is closed as **panel-missing**, not as a negative association. Nearby panel genes (CLDN3, CD8A, IFNG, CD274) were not used as stand-ins for CLDN4. No dual-high cut. No TACSTD2 audit.

---

## What this is not

- Not a private or missing series. FTP log2CPM and series matrix downloaded and opened.
- Not a CLDN4 vs DCB / ORR / PFS / CD8 / IFN / PD-L1 result. Those tests have no CLDN4 vector.
- Not a CLDN3 leftover. CLDN3 is on the OBP; it is a different gene.
- Not a CD8A / IFNG / CD274 vs response substitute.
- Not a re-open of any prior TACSTD2 leftover page on this or another accession.

---

## Empty leftover table

| dataset | public 2L PD-1/PD-L1 NSCLC | n_patient | CLDN4 n | CLDN4–DCB | CLDN4–ORR | CLDN4–PFS | CLDN4–CD8A | CLDN4–IFN | CLDN4–CD274 | note |
| --- | --- | ---: | ---: | --- | --- | --- | --- | --- | --- | --- |
| GSE161537 NIVOBIO | **yes** | **82** | **0** | — | — | — | — | — | — | HTG OBP **panel-missing** CLDN4; no proxy |

Numeric row: `tables/one_row.tsv`.

---

## 结论

GSE161537（NIVOBIO，HTG EdgeSeq OBP，二线 PD-1/PD-L1 NSCLC）公开、患者单位 n=82，且有 RECIST / PFS。但沉积 panel **不含 CLDN4**。CLDN4 对 DCB / ORR / PFS / CD8A / IFN / CD274 的检验为 **空，n=0**。未用 CLDN3 或其他 panel 基因代替。未编造表达值。
