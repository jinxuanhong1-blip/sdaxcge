# Hunt 3: leftover public TROP2-ADC RNA (CLDN4 / IFN / TJ)

**Slice:** `methods/trop2_adc_hunt3_cldn4/`
**This slice is additive.** It does not re-score hunt 1 or hunt 2 matrices.
**Search date:** 2026-08-17. Accessions and file sizes were opened on that date. No fold-changes were invented.

**Hard skips (requested):** GSE312098, GSE311016, E-MTAB-16433, GSE302284.
**Hunt 2 hits (also skipped, not re-scored):** GSE311016, GSE304294, E-MTAB-16433.

IMMU132 = sacituzumab govitecan (SG / Trodelvy). **Not SKB264 / sac-TMT.** Honest n is the number of biological units in a deposited ADC-vs-control matrix, not hashed cell counts.

---

## Accession table (scored leftover series)

Master file: `tables/scored_series.tsv`.

| Series | Tissue | n | CLDN4 vs vehicle | IFN / MHC-I / APM | KEGG TJ |
|---|---|---|---|---|---|
| *(empty)* | — | **0** | — | — | — |

**Honest n scored in hunt 3 = 0.** No leftover public processed TROP2-ADC vs vehicle / vs untargeted-ADC matrix was both (i) outside the skip + hunt 2 hit list and (ii) under 2 GB.

---

## Hunt (what was opened)

NCBI GEO `esearch` (gds, 2026-08-17):

| Query | GSE count | What it actually was |
|---|---:|---|
| `sacituzumab[All Fields] AND gse[ETYP]` | 10 | Same 10 as hunt 2 minus the IMMU132-only pair (GSE312098 / GSE311016 / GSE304294 appear only when IMMU132 is in the query) |
| `IMMU132 OR IMMU-132 AND gse[ETYP]` | 13 | Hunt 2 set: GSE312098, GSE311016, GSE304294 + the 10 above |
| `Trodelvy` / `govitecan` AND `gse[ETYP]` | 10 | Same 10; no new GSE |
| `datopotamab` / `Dato-DXd` AND `gse[ETYP]` | **0** | No GSE |
| `DS1062[All Fields] AND gse[ETYP]` | 2 | **False hits:** GSE173065, GSE141884 (S. cerevisiae Pol δ / Okazaki, not DS-1062a) |
| `SKB264 OR sac-TMT OR tirumotecan` | **0** | No GSE |
| `TROP2 ADC OR Trop-2 ADC` | 1 | GSE302284 only (skipped) |
| `hRS7` / `hRS7-SN-38` / `ESG401` / `BAT8003` / `OBI-992` / `PF-06664178` / `SHR-A1921` / `JS108` | **0** | No GSE |
| `deruxtecan AND gse[ETYP]` | 7 | HER2/HER3 ADCs (T-DXd, HER3-DXd), not TROP2 |
| `ADC AND TACSTD2 AND gse[ETYP]` | 3 | GSE158506, GSE54352, GSE54351 — FAT1 / Lkb1-Pten, not ADC treatment |

SRA / BioProject leftovers that are **not** a new open GSE:

| Accession | What it is | Scored? |
|---|---|---|
| PRJNA1372170 | GSE312098 BioProject | No — skip |
| PRJNA1368476 | GSE311016 BioProject | No — skip |
| PRJNA1300296 | GSE304294 BioProject | No — hunt 2 hit |
| PRJEB105687 / ERP186838 | E-MTAB-16433 | No — skip |
| PRJEB111028 / ERP191666 | E-MTAB-16843 | No — >2 GB |
| PRJEB111088 / ERP191723 | E-MTAB-16849 | No — >2 GB full contrast |
| PRJNA1168164 | GSE278664 | No — prexasertib biopsies |
| PRJNA1169861 / SRP536975 | GSE235812 Dato-DXd paper | No — untreated P0/P1 |
| PRJNA754211 / 754210 | SG resistance WES+RNA | No — dbGaP/EGA phs002555 |
| PRJNA1357832 | TROP2 antibody–PROTAC (ASA), 8 SRA runs | No — not an ADC; raw only |
| PRJEB102638 / ERP184030 | TROP2.Saci/YW597.10 **bispecific** (anti-Frizzled), HPAF-II 48 h | No — not an ADC; FASTQ + Zenodo scripts, no count matrix |

ArrayExpress / BioStudies (`type=study`, 2026-08-17): the only E-MTAB hits for Trodelvy / sacituzumab / govitecan are the ORCA-HD `TROP2_in_CRC` family already opened in hunt 2.

PubMed title/abstract hits that mention GSE / GEO / ArrayExpress / RNA-seq with sacituzumab (PMID 42218654 SURE-01; 41686301 UTUC IHC; 38986529 / 37913680 pan-cancer TACSTD2; 40509933 CCA PDOs) do **not** deposit a new public ADC-vs-vehicle matrix.

Full open/skip log: `tables/hunt_catalog.tsv`.

---

## Leftover matrices that are real SG RNA but still not scored

These are the only leftover **processed** SG treatment objects. Hunt 2 already catalogued them as >2 GB. Hunt 3 re-opened the file lists and metadata (not the 5 GB count matrices).

### E-MTAB-16843 — CRC organoid SG vs IgG1-SN-38

Metadata (`HD4246_SG_org_metadata.txt`, 46,906 cells): `drug` is SG (23,788) vs UNADC (23,118). Timepoints 0–12 h plus washout. **Both arms are present.**

Processed files:

| File | Size | Usable as SG vs control? |
|---|---:|---|
| `HD4246_SG_org_counts.txt` | **4.67 GB** | Yes in principle; over the 2 GB processed-matrix cutoff |
| `HD4246_SG_org_counts_log1p.txt` | **5.90 GB** | Same |

No smaller split exists.

### E-MTAB-16849 — CRC liver-met SG vs IgG1-SN-38

Full mets metadata (`HD4246_SG_mets_metadata.txt`, 44,779 cells): `treatment` is ADC (23,560) vs SG (21,219). **Both arms are present in the full object.**

| File | Size | Treatment in metadata |
|---|---:|---|
| `HD4246_SG_mets_counts.txt` | **4.96 GB** | ADC + SG |
| `HD4246_SG_mets_magic.csv` | 23.4 GB | — |
| `HD4246_SG_0h-48h_counts.txt` | 1.26 GB | **SG only** (11,385 cells; 0 / 6 / 24 / 48 h) |
| `HD4246_SG_48h-120h_counts.txt` | 1.55 GB | **SG only** (13,984 cells; 48 / 72 / 120 h) |

The <2 GB pieces cannot test SG vs untargeted ADC. Scoring them as a time course inside SG would not be an ADC-vs-control contrast and would not be additive CLDN4/IFN/TJ evidence.

---

## What hunt 3 did not pretend to find

- No new GEO Series with sacituzumab / IMMU132 / Trodelvy / govitecan beyond the hunt 2 list.
- No open datopotamab / Dato-DXd / DS-1062a expression Series (the two `DS1062` GSE records are yeast).
- No public SKB264 / sac-TMT / MK-2870 RNA.
- No public processed matrix for clinical SG resistance (phs002555) or TROPION Dato-DXd lung RNA (EGA).
- GSE278664 was re-read from SOFT: every sample `treatment:` field is Prexasertib; overall design is 15 BRCAmut + 20 BRCAwt HGSOC biopsies. Honest n for SG vs vehicle = 0.
- PRJNA1357832 is a TROP2-targeted **PROTAC** conjugate, raw SRA only.
- PRJEB102638 sample attributes say `TROP2.Saci/YW597.10-bispecific-48 h`, not sacituzumab govitecan.

---

## Verdict

Hunt 3 is empty on purpose, not because the search was not run.

The public TROP2-ADC RNA that can actually be scored for CLDN4 / IFN / TJ remains the four series already claimed:

| Already claimed | By | Why hunt 3 left it |
|---|---|---|
| GSE312098 CX-1 2 d | hunt 1 (given) | requested skip |
| GSE311016 CRC PDX day 29 | hunt 2 | requested skip + hunt 2 hit |
| GSE304294 ESCC 1 d (n=2 vs 3) | hunt 2 | hunt 2 hit |
| E-MTAB-16433 CRC PDOX 28 d | hunt 2 | requested skip + hunt 2 hit |

The only other real SG-vs-control processed RNA (E-MTAB-16843 organoid; E-MTAB-16849 full liver-met) is still **>2 GB**. The <2 GB E-MTAB-16849 splits are SG-only. That is not a hidden fifth series.

Do not fill CLDN4 / IFN / TJ numbers from an empty table.
