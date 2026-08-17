# Hunt2: additive public TROP2-ADC RNA for CLDN4 / IFN / TJ

**Slice:** `methods/trop2_adc_hunt2_cldn4/`
**Date opened:** 2026-08-17
**Scope:** GEO and ArrayExpress only. Additive to hunt1 (`methods/trop2_adc_cldn4_extra/`).
**Skip (not re-opened for scoring):** GSE312098, GSE311016, E-MTAB-16433, GSE302284, EGA, dbGaP.
**This slice is a hunt, not a re-score.** No accessions, sample counts, or fold-changes were invented.

Honest n is the number of biological units in the deposited record (GEO `n_samples`, GEO sample IDs, or ArrayExpress ENA/SDRF libraries), not hashed cell counts.

---

## Verdict

**No new public GEO or ArrayExpress TROP2-ADC vs control RNA series.**

`tables/accession_table.tsv` is honestly empty (header only). Hunt2 did not add a series that can test CLDN4 / IFN / MHC-I / APM / tight-junction after TROP2-ADC treatment beyond what hunt1 already catalogued.

That is not “no public TROP2-ADC RNA exists.” Hunt1 already scored three SG matrices (GSE311016, GSE304294, E-MTAB-16433) and took GSE312098 as given. Hunt2 asked whether anything else is sitting in GEO/ArrayExpress. After opening the leftover hits, the answer is no.

---

## What was searched

NCBI GEO `esearch` (`db=gds`) and EBI BioStudies `collection=arrayexpress`, 2026-08-17. Query log: `tables/queries.tsv`. Every hit listed below was opened (esummary and/or SOFT brief / ArrayExpress study + file list) before it was called empty or leftover.

**GEO, TROP2-ADC names**

| Query | Hits | After opening |
|---|---:|---|
| `(sacituzumab OR IMMU132 OR IMMU-132 OR Trodelvy OR govitecan) AND gse[ETYP]` | **13** | Same 13 as hunt1. No 14th series. |
| `datopotamab` (any GEO type) | **0** | No Dato-DXd RNA. |
| `"DS-1062"`, `DS-1062a`, `Datroway` | **0** | — |
| `Dato-DXd OR "DS-1062" OR DS1062 OR Datroway` | 2 | GSE173065, GSE141884 = yeast Pol-δ / Okazaki. Token collision, not datopotamab. |
| `SKB264`, `"MK-2870"`, `tirumotecan` | **0** | No sac-TMT RNA. |
| `"sac-TMT"` | 2 | GSE329257, GSE8614 = tandem-mass-tag / hashtag, not sac-TMT. |
| `PF06664178`, `RN927C`, `BAT8003`, `ESG401`, `hRS7` | **0** | Hyphenated codes were not trusted (NCBI splits `PF-06664178`). |
| `"TROP2 ADC"` / `"TROP2-ADC"` | 1 | GSE302284 only (skipped). |
| `"anti-TROP2" AND (ADC OR antibody-drug)` | 2 | GSE312098, GSE311016 only. |
| `(ADC OR "antibody-drug conjugate") AND (TROP2 OR TACSTD2)` | 14 | 13-set plus FAT1 / Lkb1-Pten false hits. |
| `(SN-38 OR SN38) AND (TROP2 OR TACSTD2)` | 4 | 13-set members plus **GSE213103** (irinotecan, not ADC). |

**ArrayExpress collection**

| Query | Hits |
|---|---|
| sacituzumab / govitecan | **E-MTAB-16433, E-MTAB-16843, E-MTAB-16849** only |
| Trodelvy | E-MTAB-16433 only |
| IMMU-132 / IMMU132 / datopotamab / Dato-DXd / DS-1062 / SKB264 / tirumotecan | **0** |
| TROP2 | 16 studies: the ORCA-HD SG trio plus Trop2+ cell-sorting and GEO mirrors. No extra ADC-vs-control RNA. |

Full open/skip log: `tables/hunt_catalog.tsv`.

---

## Already known (not additive)

These are real TROP2-ADC RNA records. Hunt2 confirms they still exist and does not add them.

| Accession | Tissue | Honest n | Deposited contrast | Hunt2 action |
|---|---|---|---|---|
| GSE312098 | CRC CX-1 | 12 samples; IMMU vs vehicle **3 vs 3** | IMMU132 ± PERK i, 2 d | Skip (requested) |
| GSE311016 | CRC PDX | **5 pairs / 10 samples** | IMMU132 vs control, day 29 | Skip (requested) |
| GSE304294 | ESCC KYSE30 | 11 samples; IMMU vs control **2 vs 3** | IMMU132 ± IACS, 1 d | Already scored in hunt1; not new |
| E-MTAB-16433 | CRC PDOX HD42466 | **4 vs 4 mice** (2 hashed libraries; 7330 cells) | Trodelvy vs vehicle, 28 d | Skip (requested) |
| E-MTAB-16843 | CRC organoid | **4 ENA samples** (GEX+HTO multiplex) | SG vs IgG1-SN-38 time course | Known; full counts 4.35 GB / log1p 5.49 GB |
| E-MTAB-16849 | CRC liver-met PDOX | **6 ENA samples** (GEX+HTO multiplex) | SG vs IgG1-SN-38 time course | Known; full counts 4.62 GB. Time-split counts 1.17 GB and 1.44 GB are trajectory splits, not a new contrast |

No organoid RNA from the GSE312098 / GSE311016 CRC paper. GSE311016 is still a regular Series, not a SuperSeries.

---

## Opened leftovers that are not TROP2-ADC vs control RNA

Every GEO series in the 13-hit sacituzumab/IMMU132 union that is not in the skip list was opened.

| Accession | Honest n | What the deposited RNA actually is |
|---|---|---|
| **GSE278664** | **35** (15 BRCAmut + 20 BRCAwt) | HGSOC biopsies, NCT02203513. Sample titles are `Tumor biopsy, P*mut` / `S*wt`. Paper is about SG + berzosertib; matrix is genotype, not SG vs vehicle. |
| **GSE309617** | **22** | TNBC PDX carboplatin-resistant vs sensitive (4 models). |
| **GSE309616** | **2** | scRNA of one of those carboplatin PDX pairs. |
| **GSE303051** | **30** | PDAC lines, ulixertinib vs DMSO. |
| **GSE302981** | **8** | Same PDAC paper, ulixertinib vs DMSO. |
| **GSE303323** | **6** | Capan-1 TRIM22 OE vs vector. |
| **GSE303106** | **9** | CUT&RUN / ChIP, not RNA ADC. |
| **GSE303105** | **2** | Mouse PDAC scRNA, RMC-6236. |
| **GSE303104** | **2** | Mouse PDAC scRNA, ulixertinib. |
| **GSE235812** | **36** | Untreated breast P0 vs matching PDX. A later Dato-DXd BCX paper cites this accession; the matrix is still not ADC vs vehicle. |
| **GSE213103** | **8** | Colon PDX FACS CSC vs non-CSC after **irinotecan (CPT-11)** vs saline. SN-38 prodrug chemotherapy, not TROP2 ADC. |

False-token hits (opened, then dropped): GSE173065, GSE141884, GSE158506, GSE54352, GSE54351, GSE329257, GSE8614.

TROP2-loss analogs (GSE289287, GSE334497, GSE245459) are not ADC treatment and were not scored.

---

## What this hunt cannot claim

- It cannot claim hunt1 was wrong. Hunt1’s three scored SG matrices still stand.
- It cannot turn E-MTAB-16843 / E-MTAB-16849 into a new CLDN4/IFN/TJ score. Those accessions were already listed; the full processed matrices are still >2 GB.
- It cannot use EGA / dbGaP (including ICARUS-LUNG01 Dato-DXd spatial and phs002555 clinical SG RNA).
- It cannot treat Trop2 KO / shTACSTD2 as TROP2-ADC on-treatment RNA.
- It cannot invent a public datopotamab or SKB264 expression series. GEO and ArrayExpress searches for those names are empty after false hits are removed.

---

## Reproduce the hunt

```bash
# GEO (example; 2026-08-17 counts are in tables/queries.tsv)
curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=%28sacituzumab%20OR%20IMMU132%20OR%20IMMU-132%20OR%20Trodelvy%20OR%20govitecan%29%20AND%20gse%5BETYP%5D&retmode=json"

# ArrayExpress collection
curl -sS "https://www.ebi.ac.uk/biostudies/api/v1/search?query=sacituzumab&collection=arrayexpress&pageSize=50"
```

Open each returned accession before adding a row to `tables/accession_table.tsv`. If a later search adds a real TROP2-ADC vs control matrix, that row belongs there with honest n. Today that table has zero rows.
