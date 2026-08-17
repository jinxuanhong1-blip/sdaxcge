# Extra public TROP2-ADC RNA (beyond GSE312098 CX-1 2d)

**Slice:** `methods/trop2_adc_cldn4_extra/`
**Given and not re-audited:** GSE312098 CX-1 IMMU132 2 day, CLDN4 log2FC −0.86 / IFN-APM up.
**This slice is additive.** Numbers below were computed from deposited processed matrices on 2026-08-17. No accessions or fold-changes were invented.

IMMU132 = sacituzumab govitecan (SG). **Not SKB264 / sac-TMT.** No lung ADC-vs-vehicle matrix was scored.

---

## Hunt (what was opened)

NCBI GEO `esearch` (gds, `gse[ETYP]`, 2026-08-17) for `sacituzumab OR IMMU132 OR IMMU-132 OR Trodelvy OR govitecan` returned **13 series**. Each was opened (esummary and/or accession page) before any scoring.

| Accession | What the deposited RNA actually is | Scored? |
|---|---|---|
| **GSE312098** | CRC CX-1 ± IMMU132 ± GSK, 2 d, FPKM | **No — given** |
| **GSE311016** | CRC PDX ± IMMU132, day 29, FPKM 5.1 Mb | **Yes** |
| **GSE304294** | ESCC KYSE30 ± IMMU132 ± IACS, 1 d, FPKM 4.4 Mb | **Yes** |
| GSE302284 | osimertinib DTP NSCLC | **No — skipped by request** |
| phs002555 / phs2555 | dbGaP clinical SG resistance WES+RNA | **No — skipped by request** |
| GSE278664 | HGSOC biopsies, BRCAmut vs wt (SG paper, wrong treatment) | No |
| GSE309617 / GSE309616 | TNBC PDX parental vs carboplatin-resistant | No |
| GSE303323 / 303051 / 303106 / 303105 / 303104 / 302981 | PDAC KRAS/ERK-inhibitor / TRIM22 OE (paper *mentions* SG combo in PDX; deposited RNA is not SG vs vehicle) | No |
| GSE235812 | untreated breast P0 vs matching PDX (SG used in other experiments) | No |

`datopotamab[All Fields] AND gse[ETYP]` = **0 hits**. SKB264 / sac-TMT were not requested and have no public GSE.

**GSE311016 is not a SuperSeries.** GEO type is a regular Expression Series (10 samples, GPL24676). `esearch GSE311016[ACCN] AND gse[ETYP]` returns one GSE. Series relations list only BioProject PRJNA1368476 — no SubSeries. The CRC paper key-resources table is: PDX RNA = GSE311016; cell-line RNA = GSE312098. Organoids in that paper were IncuCyte / IHC, not a deposited RNA matrix.

**ArrayExpress / ORCA-HD `TROP2_in_CRC` (different paper):**

| Accession | Model | Processed size | Scored? |
|---|---|---|---|
| **E-MTAB-16433** | CRC PDOX SG vs vehicle, 28 d, scRNA | counts 714 Mb; log1p 885 Mb | **Yes** (both <2 GB) |
| E-MTAB-16843 | CRC organoid SG vs IgG1-SN-38 time course | counts 4.67 GB; log1p 5.90 GB | **No** (>2 GB) |
| E-MTAB-16849 | CRC liver-met SG vs untargeted ADC | counts 4.96 GB; MAGIC 23.4 GB | **No** (full matrix >2 GB; <2 GB files are SG-only trajectory splits) |
| E-MTAB-16583 / 16585 / 16836 / 16835 | untreated cohort / untreated PDOX / FOLFIRI / mouse MDO | — | **No** (not TROP2-ADC vs vehicle) |

Full open/skip log: `tables/hunt_catalog.tsv`.

---

## Methods (scored matrices only)

- **GSE311016:** author FPKM, UTF-16, symbol `gene_name`. Columns `C_{model}` / `T_{model}` for PDX-36, 82, 83, 114, 196. Paired *t* on `log2(FPKM+1)`, n=5 pairs. One tumor per model per arm.
- **GSE304294:** author FPKM, UTF-8. Columns `OX1_*` (n=3) / `OX2_*` (n=2) / `OX3_*` (n=3) / `OX4_*` (n=3) match series-matrix sample order Control / IMMU132 / IACS010759 / Combination. Welch *t* on `log2(FPKM+1)` for IMMU132 vs Control only. Combination arm not used.
- **E-MTAB-16433:** author processed UMI counts (genes × 7,330 cells) + metadata. Mouse-level pseudobulk (sum UMIs) → CPM → `log2(CPM+1)`. Welch *t*, n=4 Trodelvy mice vs 4 vehicle mice. **Library is aliased with treatment** (all vehicle in S01, all SG in S02).
- Gene sets frozen in `gene_sets.py` before looking at extra fold-changes: CLDN4, TACSTD2; compact IFN/MHC-I panel; MHC-I; APM; broader IFN; KEGG tight-junction core. Set test = Mann–Whitney of gene-level log2FC vs the rest of the matrix (not preranked NES; FPKM/pseudobulk is not DESeq2).
- Genome-wide BH *q* is exploratory. Honest n is the number of biological units in the deposited matrix, not hashed cell counts.

Reproduce:

```bash
python3 methods/trop2_adc_cldn4_extra/download.py
python3 methods/trop2_adc_cldn4_extra/analyze.py
```

---

## Extra series actually scored

Master table: `tables/scored_series.tsv`. Gene-level: `tables/key_genes_extra_series.tsv`. Set-level: `tables/geneset_stats_extra_series.tsv`.

### 1. GSE311016 — CRC PDX, IMMU132 vs vehicle, day 29

**Tissue: CRC. n=5 pairs. Not lung. Not SKB264. Not GSE312098.**

| Gene | mean FPKM vehicle | IMMU132 | log2FC | paired p | q |
|---|---:|---:|---:|---:|---:|
| TACSTD2 | 29.2 | 15.9 | −0.42 | 0.27 | 0.60 |
| **CLDN4** | 129 | 97.8 | **−0.44** | **0.060** | 0.60 |
| CLDN7 | — | — | −0.26 | 0.033 | 0.60 |
| OCLN | — | — | −0.49 | 0.070 | 0.60 |
| ISG15 | — | — | +1.26 | 0.32 | 0.60 |
| HLA-A | — | — | +0.005 | 0.98 | 0.99 |
| B2M | — | — | −0.36 | 0.43 | 0.65 |
| TAP1 | — | — | −0.70 | 0.038 | 0.60 |

Per-model CLDN4 FPKM (vehicle → IMMU132): PDX-36 110→57; 82 179→121; 83 120→117; 114 132→124; 196 103→70. All five models go down; the paired mean log2FC does not reach p<0.05.

| Set | n tested | up/down | median log2FC | MW p |
|---|---:|---:|---:|---:|
| C4 IFN/MHC-I panel | 6 | 5/1 | +0.18 | 0.056 |
| MHC-I | 8 | 6/2 | +0.11 | 0.16 |
| APM | 19 | 10/9 | +0.005 | 0.56 |
| IFN (broad) | 55 | 35/19 | +0.082 | **0.0078** |
| KEGG TJ | 26 | 9/16 | −0.11 | **0.0085** |

**Honest read:** same paper as the given CX-1 2d series, later (day 29) CRC PDX. CLDN4 trends down and the TJ set is down; APM/MHC-I are null; broad IFN is weakly up and is not a clean on-treatment opening. Do not pool with GSE312098. ISG15 is up in some models and down in others (PDX-83 235→67 FPKM).

Figures: `figures/fig_extra1_gse311016_pdx_pairs.png`, `figures/fig_extra2_gse311016_keygenes.png`.

### 2. GSE304294 — ESCC KYSE30, IMMU132 vs vehicle, 1 day

**Tissue: ESCC (esophageal squamous), not lung, not breast. n=2 vs 3. Underpowered. Not SKB264.**

Confirmed real GEO Series (GPL24676) with processed `GSE304294_gene_fpkm.txt.gz` (4.4 Mb). IMMU132 has only two deposited replicates.

| Gene | mean FPKM vehicle | IMMU132 | log2FC | Welch p | q |
|---|---:|---:|---:|---:|---:|
| TACSTD2 | 199 | 446 | **+1.16** | 0.0099 | 0.14 |
| **CLDN4** | 73.7 | 139 | **+0.91** | **3.9×10⁻⁵** | **0.045** |
| CLDN7 | — | — | +0.52 | 8.8×10⁻⁴ | 0.066 |
| OCLN | — | — | +1.17 | 0.0019 | 0.082 |
| IFIT1 | — | — | −0.61 | 2.9×10⁻⁴ | 0.054 |
| ISG15 | — | — | +0.21 | 0.036 | 0.24 |
| HLA-A | — | — | +0.62 | 0.011 | 0.15 |
| B2M | — | — | +0.51 | 0.0084 | 0.13 |
| TAP1 | — | — | +0.50 | 0.0061 | 0.12 |

| Set | n tested | up/down | median log2FC | MW p |
|---|---:|---:|---:|---:|
| C4 IFN/MHC-I panel | 6 | 5/1 | +0.27 | 0.022 |
| MHC-I | 8 | 8/0 | +0.50 | 4.4×10⁻⁶ |
| APM | 19 | 17/2 | +0.50 | 3.4×10⁻⁸ |
| IFN (broad) | 55 | 33/20 | +0.012 | 0.026 |
| KEGG TJ | 26 | 19/6 | +0.15 | 2.1×10⁻⁴ |

**Honest read:** this 1-day ESCC SG analog does **not** reproduce CLDN4/TJ down. CLDN4 and the TJ set go **up**. APM and MHC-I go up. n=2 in the IMMU arm. Do not average with CX-1 or with the CRC PDX.

Figures: `figures/fig_extra3_gse304294_keygenes.png`, `figures/fig_extra3b_gse304294_points.png`.

### 3. E-MTAB-16433 — CRC PDOX scRNA, Trodelvy vs vehicle, 28 days

**Tissue: CRC PDOX (HD42466). Different paper from GSE311016. Processed and <2 GB.**

Metadata (7,330 cells): 4 vehicle mice (3,042 cells, library S01) vs 4 Trodelvy mice (4,288 cells, library S02). Each sample_id is one mouse (`pooled_by_multiplexing_N=4` in the SDRF). Scoring is mouse pseudobulk, not a single-cell test. **Treatment and 10x library are completely aliased.**

Mean CPM: CLDN4 475 → 370; TACSTD2 256 → 189.

| Gene | log2FC | Welch p | q |
|---|---:|---:|---:|
| TACSTD2 | −0.42 | 0.034 | 0.32 |
| **CLDN4** | **−0.35** | **0.063** | 0.38 |
| CLDN1 | −0.54 | 0.0058 | 0.19 |
| OCLN | −0.70 | 0.0041 | 0.18 |
| TJP1 | −0.58 | 0.0045 | 0.18 |
| F11R | −0.45 | 0.0011 | 0.14 |
| ISG15 | +0.029 | 0.92 | 0.97 |
| HLA-A | +0.11 | 0.68 | 0.86 |
| B2M | +0.025 | 0.90 | 0.96 |

Per-mouse CLDN4 CPM: vehicle 427, 607, 426, 441 vs SG 379, 396, 347, 358.

| Set | n tested | up/down | median log2FC | MW p |
|---|---:|---:|---:|---:|
| C4 IFN/MHC-I panel | 6 | 4/2 | +0.034 | 0.99 |
| MHC-I | 8 | 4/4 | −0.060 | 0.24 |
| APM | 19 | 9/10 | −0.066 | 0.093 |
| IFN (broad) | 52 | 28/24 | +0.034 | 0.84 |
| KEGG TJ | 24 | 6/18 | −0.31 | **5.2×10⁻⁵** |

**Honest read:** another CRC SG-vs-vehicle matrix, not from the GSE312098 paper. CLDN4 trends down (p=0.063) and the TJ set is down. IFN / MHC-I / APM are null. Library batch is aliased with treatment, so this is not an independent confirmation of the CX-1 2d IFN-APM opening.

Figures: `figures/fig_extra4_emtab16433_keygenes.png`, `figures/fig_extra4b_emtab16433_mouse_points.png`.

---

## What was not scored (verified, not empty-handed)

- No additional GSE311016 SubSeries (it is not a SuperSeries).
- No organoid RNA from the IMMU132 CRC Cell Reports Medicine paper.
- E-MTAB-16843 organoid SG time course is real processed RNA but **>2 GB**.
- E-MTAB-16849 liver-met SG time course is real processed RNA but **>2 GB**.
- No open datopotamab RNA series.
- No public SKB264 RNA (out of scope).
- GSE302284 and phs2555 left untouched as requested.

---

## Verdict

| Extra series | Tissue | n | CLDN4 vs vehicle | IFN / MHC-I / APM | KEGG TJ |
|---|---|---|---|---|---|
| GSE311016 | CRC PDX | 5 pairs | −0.44, p=0.060 | IFN set up (MW p=0.0078); APM/MHC-I null | down (MW p=0.0085) |
| GSE304294 | ESCC KYSE30 | 2 vs 3 | **+0.91**, p=3.9×10⁻⁵ | APM/MHC-I up; IFN set weak | **up** (MW p=2.1×10⁻⁴) |
| E-MTAB-16433 | CRC PDOX | 4 vs 4 mice | −0.35, p=0.063 | all three null | down (MW p=5.2×10⁻⁵) |

The given CX-1 2d CLDN4-down / IFN-APM-up pattern is **not** a general public TROP2-ADC RNA rule. The same paper’s day-29 CRC PDX is a CLDN4/TJ trend without a clean APM opening. The only other open bulk IMMU132-vs-vehicle matrix (ESCC, n=2) goes the other way on CLDN4/TJ.

These are extra public SG matrices. They are not SKB264, not internalization, and not a license to pool.

Set-level forest: `figures/fig_extra5_geneset_median_log2fc.png`.
