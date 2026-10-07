# TROP2-ADC → tumor-cell immune-signature screen

**Goal.** Find PUBLIC proteomics or spatial-transcriptomics datasets in which a
**TROP2-ADC** (sacituzumab govitecan / IMMU-132, SKB264 / sac-TMT, datopotamab
deruxtecan) or a **bare anti-TROP2 antibody** raises tumor-cell immune signatures
— antigen-presentation machinery (APM) up, type-I/II interferon ISGs up, STING
axis up, with CLDN4 / tight-junction down — as mechanistic support for an
**SKB264 (TROP2-ADC) + ICI** lung-cancer paper.

**Date:** 2026-10-07. **Analyst:** Cloud agent (opus). All numbers below are
computed from author-provided processed tables; none are fabricated. Null and
opposite-direction results are reported as-is.

---

## 1. Headline verdict

- **No qualifying PROTEOMICS dataset exists** in PRIDE / ProteomeXchange /
  MassIVE / jPOST / iProX for any TROP2-ADC or anti-TROP2-antibody perturbation.
  Every TROP2 proteomics hit is about TROP2 *biology* (shedding, interactome,
  xenograft abundance, LYTAC/degrader), not an ADC/antibody treatment arm.
- **No qualifying SPATIAL-transcriptomics dataset exists** (GeoMx/CosMx/Xenium/
  Visium). The only TROP2 spatial data found (GSE345112/GSE345045) is descriptive
  GeoMx profiling of rare salivary/odontogenic carcinomas with **no** TROP2-ADC
  perturbation.
- **Three qualifying TROP2-ADC perturbation datasets DO exist, as bulk RNA-seq**,
  all using **IMMU-132 (sacituzumab govitecan)** with a clean treated-vs-control
  arm. These are the closest available public evidence and directly measure every
  target immune gene. They were scored.

**What the RNA-seq shows (TROP2-ADC vs vehicle):**

| Dataset | Model | APM | IFN/ISG | STING | CLDN4 (want ↓) | Verdict |
|---|---|---|---|---|---|---|
| **GSE312098** | CX-1 CRC cells, in vitro | **↑ +0.42, p=6.2e-5** | **↑ +0.36, p=2.1e-3** | **↑ +0.23, p=4.3e-3** | **↓ −0.86, p=1.7e-5** | **Supports** (all 4 in hypothesized direction, all significant) |
| **GSE304294** | KYSE30 ESCC cells, in vitro | **↑ +0.54, p=3.9e-4** | +0.06, p=0.36 (null) | **↑ +0.20, p=9.0e-3** | **+0.91, p=3.9e-5 (UP, opposite)** | **Partial** (APM+STING up; IFN flat; CLDN4 wrong way) |
| **GSE311016** | CRC PDX, in vivo (paired) | +0.09, p=0.80 (null) | +0.18, p=0.60 (null) | +0.20, p=0.35 (null) | −0.44, p=0.060 (trend ↓) | **Directionally consistent, not significant** |

Δ = mean difference of per-sample signature score in log2(FPKM+1) units
(treated − control); p = per-signature score test (Welch two-sample; paired for
GSE311016).

**Bottom line for the paper.** Public data provide **one clean positive**
(GSE312098: a TROP2-ADC induces the full APM↑ / IFN↑ / STING↑ / CLDN4↓ program in
colorectal tumor cells), **one partial** (GSE304294: APM and STING up, but IFN
unchanged and CLDN4 up), and **one null-but-consistent** in-vivo PDX (GSE311016).
The mechanism the SKB264+ICI paper proposes is **supported in vitro and
directionally echoed in vivo, but not reproduced in a lung model and not in
proteomic or spatial media** — those datasets do not exist publicly yet.

---

## 2. Search coverage (what was queried, and the outcome)

| Repository | Query terms | Result |
|---|---|---|
| **GEO** (NCBI gds) | sacituzumab; IMMU-132/IMMU132; datopotamab; Dato-DXd; DS-1062; SKB264; sac-TMT; tirumotecan; "TROP2 ADC"; TACSTD2+ADC; +spatial/Visium/CosMx/Xenium/GeoMx; +proteomics | 3 qualifying (GSE304294/312098/311016); rest disqualified (below). datopotamab/SKB264/tirumotecan: 0 hits. |
| **PRIDE / ProteomeXchange** | sacituzumab; TROP2; TACSTD2; datopotamab; SKB264; "TROP2 ADC" | 0 TROP2-ADC perturbation. TROP2 hits are shedding/interactome/xenograft/LYTAC only. |
| **MassIVE** (full catalog grep) | sacituzumab; TROP2; TACSTD2; datopotamab; SKB264; IMMU-132; Dato-DXd | 0 (1 TROP2 EV biomarker set, no ADC arm). |
| **jPOST** | TROP2 | 0 relevant. |
| **iProX** (PROXI) | SKB264; TROP2; sacituzumab; datopotamab; TACSTD2 | 0 returned. |
| **Zenodo** | "SKB264"; "sacituzumab govitecan"; "datopotamab deruxtecan"; "TROP2 antibody-drug conjugate" | Only clinical/review PDFs; no omics tables. |
| **figshare** | SKB264; sacituzumab; datopotamab; "TROP2 ADC" | No qualifying omics tables. |
| **BioStudies / ArrayExpress** | sacituzumab; datopotamab; SKB264; "TROP2 ADC" | Token-inflated hit counts; no TROP2-ADC perturbation omics beyond the GEO series. |
| **Europe PMC** | 13 drug terms × 11 assay terms, full-text + text-mined accessions (`scripts/search_epmc.py`) | Literature cross-check; the drug×omics papers point back to the same GEO series. Dataset *existence* is established authoritatively by the direct repository searches above, which found no qualifying proteomics/spatial TROP2-ADC data. |

---

## 3. Datasets examined and DISQUALIFIED (honest audit)

| Accession | Why it looked relevant | Disqualifying reason (Constraint 1 or no perturbation) |
|---|---|---|
| **GSE278664** | Title: "…sacituzumab govitecan and berzosertib…DNA damage response in ovarian cancer" | All 35 deposited samples are **Prexasertib-treated patient biopsies**; the processed matrices contain **no sacituzumab (nor control) arm**. No usable TROP2-ADC-vs-control contrast. |
| **GSE309617 / GSE309616** | "sacituzumab" keyword; TNBC therapeutic synergy | Deposited arms are **Untreated vs Carboplatin** (TNBC PDX). No TROP2-ADC arm. |
| **GSE302284** | "Targeting TROP2" in EGFR-mutant NSCLC (lung!) | Perturbation is **anti-TROP2 CAR-T cellular therapy**, not an ADC or bare antibody → outside Constraint 1. |
| **GSE345112 / GSE345045** | TROP2 + spatial (GeoMx) | Descriptive comparative profiling of CCOC/HCCC; TROP2 is only a marker; **no TROP2-ADC perturbation**. |
| **GSE69160** | "TACSTD2 + ADC" keyword | **EpCAM**-aptamer-toxin conjugate, not TROP2. |
| PRIDE PXD065965, PXD039272, PXD028335, PXD058514/23/39, MassIVE EV set | TROP2 proteomics | TROP2 biology (shedding/interactome/xenograft abundance/degrader), **no ADC or antibody treatment arm**. |

---

## 4. Qualifying datasets (pass Constraint 1)

All three use **IMMU-132 = sacituzumab govitecan**, an explicitly allowed TROP2-ADC
(humanized anti-TROP2 hRS7 IgG conjugated to SN-38). Comparator is vehicle/untreated
("Control"). Non-ADC arms present in the same experiments are mechanistic
**specificity controls**, not analogs/off-target-payload ADCs, and are **never
scored as TROP2-ADC evidence** (see §7).

| Accession | System | Tissue | TROP2-ADC arm | Control | Same-experiment non-ADC arms | Design |
|---|---|---|---|---|---|---|
| GSE304294 | KYSE30 cells | Esophageal SCC | IMMU-132, 1 day (n=2) | Vehicle (n=3) | IACS-010759 (OXPHOS inhibitor); IMMU-132+IACS combo | in vitro |
| GSE312098 | CX-1 cells | Colorectal | IMMU-132, 2 day (n=3) | Vehicle (n=3) | GSK2606414 (PERK inhibitor); IMMU-132+GSK combo | in vitro |
| GSE311016 | Patient-derived xenografts | Colorectal | IMMU-132 10 mg/kg ×3 wk (n=5) | Vehicle (n=5) | — | in vivo, paired by PDX model |

Column→arm mapping was confirmed from each sample's GEO `Library name` field
(e.g. GSE304294 OX1=Control, OX2=IMMU-132; GSE312098 X_1-3=Control, X_4-6=IMMU-132;
GSE311016 C_*=Control, T_*=IMMU-132, paired by PDX id).

---

## 5. Methods

- Input: author-provided **processed `*_gene_fpkm.txt.gz`** matrices only
  (4–6 MB each; re-downloadable via `scripts/fetch_data.sh`). **No raw/multi-GB
  matrices were downloaded.**
- Transform: **log2(FPKM + 1)**. Duplicate gene symbols collapsed by max.
- Per gene: **log2FC = mean_log2(treated) − mean_log2(control)**; two-sided
  **Wilcoxon rank-sum (Mann-Whitney U)** treated vs control.
- Per signature: per-sample **score = mean log2(FPKM+1) across detected genes**,
  then **Welch two-sample t-test** (plus **paired t-test** for GSE311016).
- Signatures (genes detected / requested):
  - **APM**: B2M, HLA-A, HLA-B, HLA-C, HLA-E, HLA-F, TAP1, TAP2, PSMB8, PSMB9, NLRC5, TAPBP.
  - **IFN/ISG**: STAT1, IRF1, CXCL9, CXCL10, CXCL11, GBP1, GBP2, GBP4, IRF7, ISG15, MX1, OAS1, IFIT1, IFIT3, IDO1.
  - **STING axis**: CGAS, STING1, TBK1, IKBKE, IRF3, IFNB1 (aliases MB21D1, TMEM173 handled).
  - **Sanity**: CLDN4 (expect ↓), TACSTD2/TROP2 (target sanity).

Outputs: `signature_scores.csv`, `per_gene_log2fc.csv`,
`specificity_context_arms.csv`.

---

## 6. Per-dataset results and verdict

### GSE312098 — CX-1 colorectal cells, IMMU-132 vs vehicle (n=3 vs 3) — **SUPPORTS**
- APM score **+0.419 log2, p = 6.2e-5** (↑ as hypothesized).
- IFN/ISG score **+0.356 log2, p = 2.1e-3** (↑).
- STING score **+0.234 log2, p = 4.3e-3** (↑).
- CLDN4 **−0.861 log2, p = 1.7e-5** (↓, as hypothesized).
- TACSTD2 (target) +0.683, p = 5.7e-4 (mRNA not reduced; ADC acts on protein).
- Per-gene fingerprint (log2FC): GBP1 +1.08, GBP2 +0.95, B2M +0.83, HLA-F +0.80,
  HLA-C +0.74, HLA-B +0.70, ISG15 +0.63, GBP4 +0.63, PSMB8 +0.60, IFIT1 +0.49,
  OAS1 +0.48, HLA-E +0.47, IFIT3 +0.45, TAP2 +0.44, TAP1 +0.43, IRF1 +0.35,
  PSMB9 +0.32, STAT1 +0.16; flat/down: CXCL9 +0.04, IDO1 0, IRF7 −0.04, MX1 −0.04,
  CXCL10 −0.07, TAPBP −0.46. (Per-gene Wilcoxon p floors at 0.1 at n=3-vs-3, so no
  single gene reaches p<0.05; the coherent direction across ~20 genes and the
  aggregate score tests carry the signal.)
- **Verdict:** clean, internally coherent positive. A TROP2-ADC reproduces the
  entire proposed program (APM↑, IFN/ISG↑, STING↑, CLDN4↓) in tumor cells.

### GSE304294 — KYSE30 esophageal-SCC cells, IMMU-132 vs vehicle (n=2 vs 3) — **PARTIAL**
- APM **+0.539 log2, p = 3.9e-4** (↑).
- IFN/ISG +0.063 log2, **p = 0.36 — null** (no interferon induction).
- STING **+0.196 log2, p = 9.0e-3** (↑).
- CLDN4 **+0.911 log2, p = 3.9e-5 — UP, OPPOSITE** to the tight-junction-down hypothesis.
- TACSTD2 +1.159, p = 9.9e-3.
- **Verdict:** APM and STING rise, but the IFN/ISG module does not move and CLDN4
  goes the *wrong* way. Partial / mixed. Note only 2 IMMU-132 replicates (reduced
  power), 1-day timepoint.

### GSE311016 — colorectal PDX, IMMU-132 vs vehicle, in vivo, paired (n=5 vs 5) — **DIRECTIONALLY CONSISTENT, NOT SIGNIFICANT**
- APM +0.086 (Welch p=0.87, paired p=0.80) — null.
- IFN/ISG +0.184 (Welch p=0.58, paired p=0.60) — null, direction up.
- STING +0.203 (Welch p=0.29, paired p=0.35) — null, direction up.
- CLDN4 −0.437 (Welch p=0.15, **paired p=0.060**) — trend down, as hypothesized.
- TACSTD2 −0.421 (paired p=0.27).
- **Verdict:** every point estimate points the hypothesized way (APM/IFN/STING up,
  CLDN4 down) but none is significant; large inter-PDX variability and a 29-day
  endpoint dilute the signal. This is tumor-cell-intrinsic signal only (human reads;
  mouse stroma excluded), i.e. not immune-infiltrate-driven.

---

## 7. Specificity controls (same-experiment non-ADC arms — context only, NOT scored as TROP2 evidence)

These address: is the immune program **specific to the TROP2-ADC**, or a generic
stress response? (`specificity_context_arms.csv`.)

- **GSE312098**, PERK inhibitor **GSK2606414 alone** vs control: APM +0.385
  (p=1.5e-3) — i.e. APM rises with the non-ADC drug too — but IFN +0.10 (ns),
  STING **−0.12** (down), CLDN4 −0.06 (ns). So the **IFN↑ / STING↑ / CLDN4↓** part
  of the program is **IMMU-132-specific**, not reproduced by the PERK inhibitor;
  APM↑ is partly shared (generic). The IMMU-132+GSK **combination** amplifies
  everything (APM +0.81, IFN +0.39, STING +0.46, CLDN4 −0.61) but is confounded by
  the partner.
- **GSE304294**, OXPHOS inhibitor **IACS-010759 alone**: APM +0.22 (p=0.015, weaker
  than IMMU-132's +0.54), IFN +0.02 (ns), STING **−0.25** (down), CLDN4 +0.50 (up).
  IMMU-132 raises APM/STING more than the OXPHOS inhibitor.

**Interpretation:** where the TROP2-ADC works (GSE312098), the interferon/STING/
tight-junction arm of the signature is attributable to the ADC rather than generic
cytotoxic stress — the strongest single piece of mechanistic support for the
SKB264+ICI rationale.

---

## 8. Caveats / limits (read before citing)

1. **Medium mismatch:** the preferred media (proteomics, spatial) yielded **zero**
   qualifying TROP2-ADC datasets. All positive evidence here is **bulk RNA-seq**.
2. **No lung model:** qualifying data are CRC and ESCC; none are NSCLC. Support is
   cross-tumor and mechanistic, not lung-specific.
3. **Payload identity:** all three use SN-38-payload **sacituzumab govitecan**, not
   SKB264 (belotecan-analog/Kthiol-linker) or Dato-DXd (DXd). Mechanism (TROP2-ADC
   → tumor-cell DNA-damage/immune program) is shared-class, but payload/linker
   differences are a real extrapolation gap.
4. **Small n / underpowered per-gene tests:** n=2–5 per arm; per-gene Wilcoxon
   p-values floor at 0.1 (n=3) or 0.2 (n=2 vs 3). The per-signature score tests are
   the primary readout; treat single-gene p-values as descriptive.
5. **In-vivo null:** GSE311016 is directionally consistent but non-significant.

---

## 9. Reproducibility

```
methods/trop2adc_proteomics_spatial_opus_20261007b/
├── RESULTS.md                      # this file
├── signature_scores.csv           # per-signature Δlog2 + Welch/paired p, all 3 datasets
├── per_gene_log2fc.csv            # per-gene log2FC + Wilcoxon p, all genes/signatures
├── specificity_context_arms.csv   # non-ADC specificity arms (context only)
├── scripts/
│   ├── fetch_data.sh              # re-download the 3 processed FPKM tables (small)
│   ├── score_trop2adc.py          # scoring pipeline (log2, log2FC, Wilcoxon, score t-test)
│   └── search_epmc.py             # Europe PMC drug×assay accession miner
└── data/                          # FPKM tables (gitignored; fetch via fetch_data.sh)
```

Reproduce: `bash scripts/fetch_data.sh && python3 scripts/score_trop2adc.py`.

**Data sources (processed tables, author-provided):**
- GSE304294 — Esophageal SCC, IMMU-132 ± IACS-010759: `GSE304294_gene_fpkm.txt.gz`
- GSE312098 — Colorectal CX-1, IMMU-132 ± GSK2606414: `GSE312098_gene_fpkm.txt.gz`
- GSE311016 — Colorectal PDX, IMMU-132 vs vehicle: `GSE311016_gene_fpkm.txt.gz`
