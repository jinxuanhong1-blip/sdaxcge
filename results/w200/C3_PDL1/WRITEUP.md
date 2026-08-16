# C3 — PD-L1/CD274 after CLDN4 loss, or after TROP2 ADC

**Verdict: not supported.**

Claim C3 has two arms. Public literature and the open matrices that actually perturb CLDN4 or a TROP2 ADC do not establish that PD-L1/CD274 rises in either case.

中文：文献里没有“CLDN4 缺失后 PD-L1 升高”的实验论文。公开转录组里两条乳腺癌 CLDN4-KO 细胞系方向相反；小鼠肺 Cldn4-KO 的 Cd274 不变。TROP2 ADC 只有一套可用的 bulk：CX-1 细胞 2 天 IMMU132 上 CD274 有统计升高但 FPKM<1，同一工作的 CRC PDX 29 天则是平的。不能把 SG+ICI 联用试验或 claudin-low / CLDN18.2 论文算进 C3。

## What would have counted

- Arm A: CLDN4/Cldn4 genetic loss (KD/KO/si/sh/CRISPR) with a PD-L1/CD274 readout that goes up.
- Arm B: a TROP2 ADC (sacituzumab govitecan / datopotamab deruxtecan / SKB264) with an on-treatment PD-L1/CD274 readout that goes up.

Association, claudin-low subtype, CLDN18.2, naked TROP2 antibody, TROP2 KO, SN-38 alone, and ICI+ADC trials without an induction assay do not count.

## Literature (Europe PMC, 2026-08-16)

`scripts/01_harvest_literature.py` pulled 25 queries. Union: **2,235** unique records (`catalog/records.tsv`). Title+abstract screen: `catalog/screened.tsv`.

| Query class | What came back |
|---|---|
| TITLE_ABS CLDN4 **and** PD-L1 | 2 (A01) + 6 (A02). All incidental co-mentions or IHC panels. |
| Any paper that states CLDN4 loss raises PD-L1 | **0** after reading the tight set and the auto-PRIORITY_A hit. |
| TROP2 ADC **and** PD-L1 | Almost all combo-trial / review / baseline-CPS papers. |
| TROP2 ADC with a measured PD-L1 *induction* endpoint | **0**. |

The one auto-flagged “PRIORITY_A” record (PMID 40531557) is a gastric-cancer IHC classification that stains Claudin-4 and PD-L1 on the same TMA. That is not loss-of-function.

Closest real CLDN4–immunity paper is PMID 41214101 (HGSC; claudin-4 OE/KD; Rab7 / type I IFN / TCR-ζ). It never reports PD-L1. PMID 41016339 (H1688 CLDN4 KO) and PMID 40892111 (pancreatitis CLDN4 shRNA) report RNA-seq and deposit no matrix; neither is a PD-L1 paper.

Arm B lookalikes that must stay out of the claim:

- PMID 41932810 / GSE334497 — TROP2/claudin-7 **barrier**. TROP2 *loss* opens T cells. Naked hRS7 ± anti-PD-1. Not “ADC raises PD-L1”.
- PMID 42526440 — TROP2 antibody-PROTAC **lowers** PD-L1 via BRD4/c-Myc.
- PMID 38048058 — high TROP2 associates with ICI resistance. Association.
- SN-38 payload papers disagree: PMID 41443824 and 39418821 say DNA damage / STING **raise** PD-L1; PMID 36641765 says SN-38 **lowers** PD-L1 via FoxO3a. That contradiction is a payload prior, not a TROP2-ADC result.

Human-read table: `catalog/human_reviewed.tsv`.

## Public omics that can actually test the claim

Live GEO / BioStudies / OmicsDI / SRA search: `scripts/03_search_public_omics.py`, log in `omics/search_log.md`.

GEO `CLDN4 AND (PD-L1 OR CD274)`: **0** series. No extra CLDN4 KD/KO transcriptome beyond the C4 set (GSE50927, GSE207704, GSE22493). That matches the C4_more hunt in this repo.

### Arm A — CD274 after CLDN4 loss

| Accession | What it is | CD274 |
|---|---|---|
| **GSE207704** | MCF7 / T47D CLDN4-KO FPKM (replicates already averaged) | MCF7 **down** (0.93 → 0.11 FPKM). T47D **up** (0.27 → 0.80). CLDN4 RNA only ~2-fold down. |
| **GSE50927** | Mouse lung Cldn4 KO (not cancer; VILI mixed in) | Cd274 logFC **+0.14**, p=0.46, FDR=1. Cldn4 logFC −6.06 (KO is real). |
| **GSE22493** | SKOV-3 siCLDN4 vs CLDN4 overexpression | **CD274 not on GPL10555**. Cannot test. C4 already showed IFN/APM does not open. |

GSE274940 (complete Cldn-null EpH4) is not CLDN4-specific; Cldn4 RNA is not even down in the KO columns. GSE22421 is C-CPE, not genetic loss.

**Arm A does not hold.** The only cancer-line KO that reports CD274 goes in opposite directions in two lines, on a collapsed FPKM table, with leftover CLDN4 RNA.

### Arm B — CD274 after TROP2 ADC

Usable bulk:

| Accession | What it is | CD274 |
|---|---|---|
| **GSE312098** | CX-1 CRC, 2-day IMMU132 (sacituzumab govitecan) vs control, n=3 | log2FC **+0.32**, Welch p=0.019. Raw FPKM 0.22 / 0.26 / 0.21 → 0.56 / 0.44 / 0.64. |
| **GSE311016** | Five CRC PDX, 29-day IMMU132 vs control (same paper) | Paired mean log2FC **−0.08**, p=0.61. **Flat.** |

So the only statistically “up” ADC contrast is a 2-day cell-line bump at FPKM < 1. The matched in-vivo 29-day PDX series does not reproduce it. That is not a general PD-L1 rise after TROP2 ADC.

GSE278664’s title mentions SG; the samples are baseline BRCA-mut vs wt biopsies. Excluded.

Open and **not quantified here** (honest gap, not a hidden negative):

- **E-MTAB-16849** — CRC liver-met scRNA-seq, SG vs untargeted IgG1-SN-38, 0–120 h. Processed files 5–23 GB.
- **E-MTAB-16433** — CRC PDOX scRNA-seq, SG vs vehicle, 28 days. log1p matrix ~0.9 GB.

Those are real TROP2-ADC series. They were not downloaded. Anyone who wants a single-cell CD274 test should start there, not invent a new accession.

Payload-only GSE222450 (SN-38, not ADC): Cd274 RNA vs vehicle is flat (log2FC −0.11, p=0.85) after mapping A=Vehicle / E=SN-38 from GEO `!Sample_description`. The paper’s claim is protein PD-L1 *down*.

### Association is not C3

cBioPortal PanCancer Atlas, no purity correction:

| Cohort | n | Spearman CLDN4–CD274 |
|---|---|---|
| LUAD | 510 | +0.05 (p=0.26) |
| LUSC | 484 | +0.02 (p=0.61) |
| BRCA | 1082 | −0.12 (p=6.9e-5) |
| OV | 300 | −0.05 (p=0.43) |
| STAD | 412 | −0.07 (p=0.14) |
| COADREAD | 592 | −0.30 (p=1.0e-13) |

CCLE cell lines (no infiltrate): **+0.20** (n=1156); lung **+0.22** (n=119, p=0.018). Tumor-intrinsic CLDN4 and CD274 move together, weakly. That is the opposite of “CLDN4-low ⇒ PD-L1-high” as a cell-autonomous rule.

## What this does not say

- It does not say PD-L1 cannot rise after DNA damage in some models. Two SN-38 papers claim that; one claims the reverse. That is not TROP2-ADC-specific and is not CLDN4-loss.
- It does not re-analyze E-MTAB-16849 / 16433.
- It does not measure PD-L1 **protein**. Every public number here is RNA except the papers that already disagree on protein.
- GSE207704 and GSE312098 are small-n. The p-value on GSE312098 is real for that table and still does not carry the claim.

## Reproduce

```bash
python3 results/w200/C3_PDL1/scripts/01_harvest_literature.py
python3 results/w200/C3_PDL1/scripts/02_screen_literature.py
python3 results/w200/C3_PDL1/scripts/03_search_public_omics.py
python3 results/w200/C3_PDL1/scripts/04_fetch_geo_metadata.py
python3 results/w200/C3_PDL1/scripts/05_analyze_cd274.py
python3 results/w200/C3_PDL1/scripts/06_tcga_cbioportal.py
python3 results/w200/C3_PDL1/scripts/07_cd274_contrasts_clean.py
```

Raw GEO downloads stay in `data/` (gitignored). Numbers used in the verdict are in `omics/cd274_contrasts/` and `omics/tcga_cbioportal/`.
