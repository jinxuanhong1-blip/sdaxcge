# Public CLDN4 / TACSTD2 loss does not cleanly open IFN / MHC-I / APM

**Claim under test.** Private CLDN4 knockdown opens an IFN / MHC-I / APM program (`IFI27`, `OAS2`, `IFIT1`, `MX1`, `ISG15`, `HLA-A`) and correlates with SKB264 (r = 0.52, private). This write-up asks two public questions, honestly:

1. In public **CLDN4 KD/KO** sets (GSE207704, GSE22493, GSE50927, GSE274940), does that same IFN/APM set go **up**?
2. In public **TACSTD2 KD/KO** sets (GSE334497, GSE289287, GSE245459), which way does IFN go?

**Short answers.**

- **CLDN4, cancer-cell data: no.** The two human cancer-cell accessions that actually lose CLDN4 (GSE207704 T47D/MCF7 CRISPR KO; GSE22493 SKOV-3 siRNA vs CLDN4-overexpression) do **not** open IFN/APM. Measured ISGs are down or mixed; the six-gene panel is mostly missing or split 3/3 with no gene at p < 0.05.
- **CLDN4, the one “up” hit is not a cancer-cell test.** GSE50927 uninjured *Cldn4*-KO whole lung is IFN/APM-up, but it is **n = 1 per group**, bulk lung (immune-cell content), and not a tumour cell. After ventilator injury the signal fades.
- **GSE274940 is not a CLDN4 test.** The deposited “Cldn-null” EpH4 line still has *Cldn4* RNA (slightly up). IFN is null.
- **TACSTD2: direction is not consistent.** T-47D Trop-2 KO xenografts (GSE289287, NRG mice) are the only clean **up**. 4T1 Trop2-KO tumours (GSE334497) have a weak Hallmark IFNα competitive shift that fails sample-level permutation; CORE6 is flat. SKOV3 shTACSTD2 without drug (GSE245459) **closes** IFN/APM, including a collapse of HLA-A.

The private r = 0.52 with SKB264 is not testable here (no public SKB264 molecular response table was used).

---

## What was tested

Pre-specified readout (not mined after looking):

| Set | Contents |
|---|---|
| USER_CORE6 | IFI27, OAS2, IFIT1, MX1, ISG15, HLA-A (mouse: Ifi27/Ifi27l2a/b, Oas2, Ifit1, Mx1, Isg15, H2-K1/H2-D1) |
| ISG_CORE | 40 canonical type-I ISGs |
| MHC1_APM | HLA-A/B/C/E/F, B2M, TAP1/2, TAPBP, PSMB8/9/10, NLRC5, ERAP1/2, … |
| HALLMARK_IFN_ALPHA / GAMMA | MSigDB Hallmark |
| REACTOME IFN-α/β and MHC-I peptide loading | MSigDB C2 |
| CTRL MYC-V1, OXPHOS | negative controls for a global “everything moves” artefact |

**Statistics.** Where per-sample matrices exist: log2 size-factor / FPKM+1 values, Welch t on genes, competitive Mann–Whitney on log2FC (set vs background), and an exact label-permutation test on a sample-level mean-z set score. Competitive p-values treat genes as independent and are **optimistic** inside co-regulated IFN modules; the permutation p is the one that respects sample n. GSE207704 and GSE50927 have no usable replicates in the GEO deposit, so those p-values are descriptive only.

SRA re-quantification of GSE207704 was attempted (AWS Open Data SRA objects download; `fasterq-dump` 3.4.1 segfaults on this image). The analysis therefore uses the author-collapsed cuffdiff FPKM.

Nothing outside `results/align_cldn4_ifn/` was written.

---

## 1. CLDN4 KD/KO

### GSE207704 — T47D and MCF7 CLDN4−/− (human breast; best cancer-cell test)

PMID 37059993. CRISPR KO vs WT, 2 reps × 2 lines. GEO only deposited **one FPKM column per genotype** (replicates already pooled). CLDN4 is down but residual mRNA remains (T47D 42.5 → 20.4, log2FC **−1.05**; MCF7 85.9 → 50.8, log2FC **−0.75**).

CORE6 coverage is poor: **IFI27, OAS2, MX1, HLA-A are absent from the file**, as are B2M, TAP1/2, NLRC5, PSMB8/9. Of the two CORE6 genes that are present:

| gene | T47D log2FC | MCF7 log2FC |
|---|---:|---:|
| IFIT1 | −0.04 | **−1.60** |
| ISG15 | **−1.94** | +0.11 |

Other measured ISGs go **down** in both lines (OAS1 −0.77 / −1.99; OAS3 −0.89 / −0.59; IFI6 −1.36 / −0.47; IFI44 −3.39 / −1.87; BST2 −4.68 / −1.38). The only both-line IFN-related **up** is HERC5 (+1.19 / +2.02). HLA-C is discordant (T47D −0.52, MCF7 +1.39). Hallmark IFNα is shifted **down** in T47D (Δmedian −0.23, competitive p_up = 1.0) and null in MCF7 (Δmedian +0.02).

**Call: does not open IFN/APM.** Residual CLDN4 mRNA, collapsed replicates, and a hole in MHC-I annotation all weaken a formal test — but the genes that *are* measured move the wrong way.

### GSE22493 — SKOV-3 CLDN4 siRNA vs CLDN4-overexpression (human ovarian)

Three two-colour Operon v3 arrays (2006). Primary numbers are GEO-deposited VALUE = log2(KD / control). The control is **CLDN4-overexpressing** cells, not parental/scramble. CLDN4 probe 17169 is missing on GSM558700; the other two arrays give deposited log2FC −1.74 and −0.71 (mean **−1.23**). Raw ScanArray intensities disagree in sign on GSM558701, so knockdown confirmation is shaky.

CORE6 on deposited VALUE: **3 up / 3 down, none at nominal p < 0.05**.

| gene | mean log2FC | t-test p |
|---|---:|---:|
| IFI27 | −0.71 | 0.73 |
| OAS2 | +0.88 | 0.46 |
| IFIT1 | **−2.03** | 0.055 |
| MX1 | +0.36 | 0.56 |
| ISG15 (G1P2) | +0.48 | 0.62 |
| HLA-A (median of probes) | −0.85 | 0.41 |

Hallmark IFNα Δmedian +0.01 (p_up = 0.61). MHC1_APM Δmedian **−0.42**. IFIT1 is the most consistent CORE6 gene and it is down.

**Call: does not open IFN/APM.** Ovarian, old array, wrong control, dirty CLDN4 QC.

### GSE50927 — Cldn4 KO whole mouse lung ± VILI

PMID 25106430. Five GSMs, **n = 1 per group**. Author edgeR tables; sign checked against *Cldn4* itself (logFC −6.1 / −10.2 / −12.7). Uninjured KO vs WT is the only CLDN4-specific public contrast where IFN/APM clearly goes **up** (CORE6 median +0.60; Isg15 +1.08, Oas2 +1.49, Ifit1 +0.58; Hallmark IFNα Δmedian +0.19). After high-injury VILI the shift is smaller; after low-injury VILI it is gone.

This is **whole lung**, so an IFN increase can be extra leukocytes, not epithelial MHC-I. Author p-values come from a no-replicate edgeR fit and should not be treated as real false-discovery control.

**Call: IFN up in uninjured KO lung, but this cannot carry a cancer-cell-intrinsic claim.**

### GSE274940 — EpH4 “Cldn-null” vs WT (mouse mammary epithelium)

PMID 41171911. The paper describes a complete claudin-family null, **not a CLDN4-only KO**. In the deposited raw counts *Cldn4* CPM is **not down** (WT ~500–700, KO ~500–960; log2FC **+0.32**). *Cldn3* and *Cldn7* are down; *Cldnd1* is almost gone. IFN/APM is null (Hallmark IFNα Δmedian −0.02; permutation p_up = 0.75).

**Call: exclude from the CLDN4 claim.** Either the KO leaves *Cldn4* RNA, or the deposit does not match the paper’s genotype.

---

## 2. TACSTD2 KD/KO — IFN direction

### GSE334497 — 4T1 Trop2 KO tumours in BALB/c (n = 5 vs 5)

The depositing paper’s own TROP2 / claudin immune-exclusion dataset (not an independent test of that paper). *Tacstd2* log2FC **−3.82** (p = 0.001). CORE6 is **flat** (Ifit1 −0.02, Isg15 −0.09, Oas2 −0.10, Mx1 +0.16, H2-K1 +0.13). Hallmark IFNα has a small competitive up-shift (Δmedian +0.10, MW p = 1.8×10⁻⁴) that **fails** the exact 252-split sample permutation (p_up = 0.21). Any IFN/T-cell signal in this bulk tumour can be infiltrate.

**Call: CORE6 does not open. Weak set-level IFN is not sample-significant.**

### GSE289287 — T-47D Trop-2 KO xenografts in NRG mice (4 KO vs 3 WT)

Author DESeq2 table plus per-sample normCounts. *TACSTD2* log2FC **−3.26**. NRG hosts have no T/B/NK cells, so this is the closest public read of a **tumour-cell-intrinsic** IFN program after Trop-2 loss. CORE6 is **6/6 up**:

| gene | log2FC | p |
|---|---:|---:|
| IFI27 | +0.66 | 0.050 |
| OAS2 | +0.81 | 0.007 |
| IFIT1 | +0.91 | 0.015 |
| MX1 | +0.93 | 0.008 |
| ISG15 | +1.44 | 1×10⁻⁵ |
| HLA-A | +0.63 | 0.071 |

Hallmark IFNα Δmedian +0.24 (competitive p_up ~ 0). Sample-level permutation on n = 4 vs 3 is 0.20 — underpowered, not a contradiction. MYC control moves **down**, so this is not a global up-shift. GEO has no Trop-2 KO *in-vitro* arm (in-vitro samples are WT and DSG2-KO only).

**Call: IFN/ISG up. This is the only public TACSTD2 set that supports the direction of the private CLDN4 claim, and it is Trop-2, not CLDN4.**

### GSE245459 — SKOV3 shTACSTD2 ± cisplatin (n = 3 vs 3)

Author FPKM. Knockdown is real (untreated TACSTD2 mean 5.40 → 0.034; log2(FPKM+1) FC **−2.63**). CLDN4 falls with it (3.15 → 0.10). Untreated shTACSTD2 **closes** IFN/APM: CORE6 5/6 down (MX1 −0.25), HLA-A **−4.17**, Hallmark IFNα Δmedian −0.41, sample permutation p_up = 1.0. The cisplatin-treated arm is mixed (IFIT1 +2.29, ISG15 +0.50, HLA-A still −2.11) and is confounded by drug × genotype; it is not a clean C4 analog.

**Call: untreated TACSTD2 KD goes the opposite way. Ovarian, FPKM, n = 3.**

---

## Verdict table

| contrast | KD/KO real? | IFN direction | Honest call |
|---|---|---|---|
| GSE207704 T47D CLDN4 KO | yes (−1.05) | down / null | does **not** open IFN/APM |
| GSE207704 MCF7 CLDN4 KO | yes (−0.75) | down / null | does **not** open IFN/APM |
| GSE22493 SKOV3 CLDN4 KD | shaky (−1.23 on 2/3 arrays) | 3/3 CORE6, set null | does **not** open IFN/APM |
| GSE50927 lung Cldn4 KO, no VILI | yes (−6.1) | **up** | up, but n=1 bulk lung, not cancer |
| GSE50927 VILI-high | yes | smaller up | same limit |
| GSE50927 VILI-low | yes | null | same limit |
| GSE274940 EpH4 Cldn-null | **no** (*Cldn4* +0.32) | null | not a CLDN4 test |
| GSE334497 4T1 Trop2 KO | yes (−3.82) | CORE6 flat; set weak | does **not** open CORE6 |
| GSE289287 T47D Trop2 KO xeno | yes (−3.26) | **up** (6/6 CORE6) | only clean TACSTD2-up |
| GSE245459 SKOV3 shTACSTD2 | yes (−2.63) | **down** | opposite |
| GSE245459 + cisplatin | yes | mixed | drug-confounded |

---

## What this does and does not say

**Does not support** “public CLDN4 KD/KO opens the same IFN/MHC-I/APM set.” The human cancer-cell CLDN4 losses go the other way or nowhere. The mouse-lung up is a different biology.

**Does not support** a single TACSTD2 → IFN direction. One xenograft is up, one in-vitro ovarian KD is down, one immunocompetent tumour is flat on CORE6.

**Does not refute** a private CLDN4-KD experiment in a different line (e.g. lung) with a different reagent. Public data simply do not reproduce that direction in the accessions named in the request.

**Cannot test** the private r = 0.52 with SKB264.

Limits that stay in force: no lung CLDN4 CRISPR cancer-cell RNA-seq is in this list; GSE207704 MHC-I genes are mostly missing from the deposit; GSE22493 is a 2006 two-colour array vs overexpression; GSE50927 is n = 1; GSE274940 is the wrong genotype; TACSTD2 sets are analogs, not CLDN4.

---

## Reproduce

```bash
python3 results/align_cldn4_ifn/scripts/00_fetch_geo_meta.py
bash   results/align_cldn4_ifn/scripts/02_download.sh
python3 results/align_cldn4_ifn/scripts/01_parse_samples.py
python3 results/align_cldn4_ifn/scripts/05_analyze.py
python3 results/align_cldn4_ifn/scripts/06_report.py
```

Outputs: `tables/verdict.tsv`, `core6_per_gene.tsv`, `geneset_results.tsv`, `ifn_apm_per_gene.tsv`, `perturbation_qc.tsv`; `figures/01–04_*.png`.

---

## 中文摘要

公开 CLDN4 KD/KO **不能**复现“IFN/MHC-I/APM 上调”。人癌细胞两套（GSE207704、GSE22493）是下调或无；唯一上调是小鼠全肺 n=1（GSE50927），不能当肿瘤细胞证据。GSE274940 的 Cldn4 RNA 没掉，不是 CLDN4 实验。TACSTD2 方向不一致：T47D 移植瘤上调（GSE289287），SKOV3 shRNA 下调（GSE245459），4T1 瘤 CORE6 持平（GSE334497）。私有 SKB264 r=0.52 这里测不了。
