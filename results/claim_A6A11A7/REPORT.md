# Claims A6, A11, A7 — public analogs (honest)

Date: 2026-08-16  
Code: `scripts/claim_A6A11A7/`  
Machine-readable verdicts: `VERDICTS.json`

These three claims were tested on **public** data only. The Zhejiang paired ICI TROP2 cohort is treated as private and was not used.

Expression in scRNA is **log1p(UMI)** on extracted gene rows (not library-size CPM). Tests are Spearman and two-sided Mann–Whitney. No multiple-testing theater.

---

## A6. TROP2-high samples have higher tumor fraction and lower CD8/NK

**Verdict: only the composition-confounded form is supported. The biological form is not.**

TACSTD2 is an epithelial gene. If “TROP2-high sample” means high *all-cell* TACSTD2, the sample is almost guaranteed to have more epithelium and therefore a lower immune (and often CD8) fraction. That is not evidence that TROP2-high *tumor cells* exclude CD8/NK.

### What the public atlases actually label

- **GSE131907** (Kim et al., 208,506 cells, 58 samples): author `Cell_subtype == Malignant cells` is **absent from all 11 primary tLung samples**. Malignant calls sit in metastases / EBUS / some tL/B. For primary tumors the honest “tumor fraction” proxy is **epithelial fraction**.
- **GSE207422** (Hu et al., neoadjuvant PD-1 + chemo): no author cell-type table. Lineages were assigned by marker-score argmax (epithelial = EPCAM/KRT; CD8 = T lineage + CD8A; NK = NKG7/KLRD1/GNLY). These are not CNV-malignant calls.

### GSE131907 primary tLung (n=11)

| TROP2 definition | vs epithelial fraction | vs CD8 fraction | vs NK fraction |
|---|---|---|---|
| All-cell TACSTD2 | **rho=0.945, p=1.1e-5** | rho=−0.21, p=0.54 | rho=+0.04, p=0.92 |
| Epithelial TACSTD2 | rho=0.15, p=0.67 | rho=+0.38, p=0.25 | rho=−0.65, p=0.032 |

All-cell TACSTD2 vs immune fraction: rho=−0.89, p=2.3e-4. That is the confound.  
CD8 *among immune cells* is not lower (rho=+0.18, p=0.59).  
Partial Spearman of all-cell TACSTD2 vs CD8 given epithelial fraction: rho=+0.19, p=0.57.

The NK hit for epithelial TACSTD2 (p=0.032) is n=11 and does **not** replicate in the mixed-origin set.

### GSE131907 all tumor-like origins (tLung, tL/B, mLN, mBrain, PE; n=37)

| TROP2 definition | vs tumor (malignant) fraction | vs CD8 | vs NK |
|---|---|---|---|
| All-cell TACSTD2 | rho=0.56, p=2.9e-4 (epithelial frac 0.89) | **rho=−0.52, p=0.001** | rho=−0.11, p=0.50 |
| Malignant-cell TACSTD2 (n=21 with ≥10 malignant cells) | rho=0.06, p=0.80 | rho=+0.19, p=0.40 | rho=+0.39, p=0.08 |
| Epithelial TACSTD2 | rho=−0.27, p=0.11 | rho=+0.26, p=0.13 | rho=+0.27, p=0.10 |

Partial all-cell TACSTD2 vs CD8 given malignant fraction: rho=−0.36, p=0.029. A residual CD8 association remains in the mixed set after a crude composition adjustment, but it is small, NK is gone, and it is not seen in primary tLung or when TROP2 is measured inside malignant/epithelial cells.

### GSE207422 (n=15; ICI-treated, unpaired)

All-cell TACSTD2 vs epithelial fraction rho=0.80, p=3e-4; vs CD8 rho=−0.61, p=0.016; vs NK rho=−0.02, p=0.95.  
Partial given epithelial fraction: CD8 rho=+0.06, p=0.83.  
Epithelial TACSTD2 vs CD8 rho=−0.39, p=0.15.

Same pattern: the CD8 drop is the epithelial-fraction drop.

### A6 bottom line

Do not write “TROP2-high tumors are CD8/NK-cold” from public scRNA unless TROP2 is computed **inside epithelium/malignant cells** (or you residualize composition). Public data do not support that stricter claim. NK is not generally lower. Primary-tumor n is too small to rescue the claim.

Tables: `A6_GSE131907_tests.csv`, `A6_GSE131907_sample_composition.csv`, `A6_GSE207422_tests.csv`.  
Figures: `A6_GSE131907_tLung_scatter.png`, `A6_GSE131907_tLung_boxplots.png`.

---

## A11. Galectin / Nectin / TGFb / CD47 with TACSTD2-high

**Verdict: CD47 and NECTIN2 travel with TACSTD2. Galectins are mixed. TGFb is not a consistent partner. The cassette is weaker than the epithelial/junctional program (EPCAM/CLDN4).**

### Cell-level (zero-inflated; overstates concordance)

GSE131907 malignant cells (n=24,784) Spearman vs TACSTD2:

| gene | rho |
|---|---|
| CLDN4 | 0.63 |
| EPCAM | 0.56 |
| **CD47** | **0.55** |
| **NECTIN2** | **0.45** |
| LGALS3 | 0.41 |
| LGALS9 | 0.34 |
| LGALS1 | 0.28 |
| PVR | 0.23 |
| TGFB1 | 0.22 |
| CD274 | 0.14 |

GSE207422 epithelial cells (n=18,205) repeat CD47 0.42, NECTIN2 0.55, NECTIN4 0.58, LGALS3 0.72, TGFB1 only 0.12, CLDN4 0.83.

### Sample-level (more honest)

GSE131907 malignant-cell means, tumor origins with ≥10 malignant cells (n=21; **no primary tLung**):

| gene | rho | p |
|---|---|---|
| CD47 | 0.86 | 5.9e-7 |
| EPCAM | 0.69 | 4.9e-4 |
| CLDN4 | 0.59 | 0.005 |
| NECTIN2 | 0.56 | 0.009 |
| LGALS9 | 0.51 | 0.019 |
| PVR | 0.44 | 0.046 |
| TGFB1 | −0.004 | 0.99 |

TCGA bulk tumors (Xena HiSeqV2, `*-01` only):

| gene | LUAD n=515 | LUSC n=502 |
|---|---|---|
| CLDN4 | 0.46 | 0.41 |
| LGALS3 | 0.38 | 0.29 |
| CD47 | 0.26 | 0.07 (NS) |
| TGFB2 | 0.26 | −0.11 |
| TGFB1 | 0.15 | 0.20 |
| LGALS9 | 0.19 | −0.06 (NS) |
| PVR | **−0.18** | 0.12 |
| CD8A | −0.10 | −0.20 |

LUAD CD47 and TGFB1 still look positive after residualizing EPCAM. LUSC CD47 does not. PVR sign-flips.

GSE207422 pretreatment bulk (n=24) is a warning, not a confirmation: TACSTD2 vs TGFB1 **rho=−0.45**, vs LGALS1 −0.48, vs EPCAM +0.70, vs CLDN4 +0.86.

### A11 bottom line

Safe public sentence: **TACSTD2-high epithelial/malignant cells co-express CD47 and NECTIN2, and often LGALS3/9; this sits inside a broader epithelial/CLDN4 program.**  
Unsafe sentence: “TACSTD2-high tumors upregulate a Galectin/Nectin/TGFb/CD47 immune-evasion module.” TGFb fails the sample-level test. NECTIN1/3/4 were not in the GSE131907 matrix (only `PVRL2`). Cell-level p-values are not biology.

Tables: `A11_GSE131907_correlations.csv`, `A11_GSE207422_epithelial_correlations.csv`, `A11_TCGA_correlations.csv`.

---

## A7. Public paired pre/post ICI TROP2 IHC or RNA

**Verdict: none found. Zhejiang stays private.**

Full inventory: `A7_paired_ICI_TROP2_hunt.csv`.

### What would count

Same patient, pre-ICI and post-ICI tumor, TROP2 protein (IHC/mIF H-score) or TACSTD2 RNA, with a downloadable matrix or scoring table.

### What exists

1. **No public paired ICI TROP2 IHC.**  
   - Bessede et al. *Clin Cancer Res* 2024 scored TROP2 mIF on pretreatment NSCLC (BIP), not paired, no patient-level table.  
   - Inoue et al. scored 160 paired lung TROP2 IHC across mixed treatments; **only 5 patients (3%) got ICI**; 76% of pairs unchanged.  
   - JCO 2025 abstract 8591: 51 paired NSCLC biopsies, TROP2 H-score changed in 29.4%; mixed systemic therapy, not ICI-stratified, data not public.  
   - Sci Rep 2025 TROP2 IHC in 110 Nivo+Ipi patients is baseline-only.

2. **No public same-patient paired ICI TACSTD2 RNA matrix.**  
   Closest: **GSE207422** (Hu et al., *Genome Med* 2023).  
   - scRNA: 3 pre-treatment biopsies + 12 post-treatment resections, **15 different patients, 0 pairs**.  
   - Unpaired TACSTD2_all median 0.33 pre vs 0.10 post, MW p=0.23; epithelial TACSTD2 0.95 vs 0.60, p=0.54.  
   - Deposited bulk log2TPM is **pretreatment only** (24 samples). MPR vs NMPR TACSTD2 p=0.34.

3. **Public ICI RNA that is pretreatment-only** (cannot replace paired Zhejiang IHC):  
   - GSE126044 n=16: TACSTD2 responder median 1595 vs non-responder 1577, p=0.58.  
   - GSE135222 n=27: TACSTD2 present (`ENSG00000184292`); no post RNA.  
   - POPLAR/OAK TACSTD2 vs atezolizumab (Bessede 2024) is pretreatment and EGA-controlled.

4. **Paired ICI RNA that is the wrong cancer or locked:** BELLINI / TONIC (TNBC, EGA DAC); Bassez 2021 (breast); Caushi 2021 / CheckMate (no TROP2 IHC table). Tempus CRC-24-0605 has some paired ICI RNA but is commercial and does not report TROP2.

### A7 bottom line

Do not imply a public paired ICI TROP2 validation set. The honest substitute is unpaired GSE207422 RNA (n_pre=3) plus pretreatment-only ICI RNA (GSE126044/GSE135222/POPLAR-OAK). Protein-level pairing is either private (Zhejiang), ICI-rare (Inoue), or abstract-only and not ICI-specific (JCO 8591).

---

## Methods (short)

- GSE131907: public annotation + raw UMI gene rows from GEO. CD8 subtypes = Exhausted/Naive/Cytotoxic/CD8-low T. NK = author `NK`.  
- GSE207422: public UMI matrix + sample metadata; marker-score lineages.  
- TCGA LUAD/LUSC: UCSC Xena `HiSeqV2`, tumor `*-01` only.  
- Scripts: `scripts/claim_A6A11A7/extract_genes.py`, `analyze_claims.py`.  
- Not done: inferCNV malignant calls on tLung; library-size CPM; Zhejiang data.

## Do not overclaim

A6 is a purity/composition story unless TROP2 is lineage-restricted.  
A11 is CD47/NECTIN2 ± galectins inside an epithelial program, not a clean four-family checkpoint module.  
A7 has no public paired ICI TROP2 IHC/RNA analog for the Zhejiang cohort.
