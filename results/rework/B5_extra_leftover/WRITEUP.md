# Extra leftover 2023–2026 lung ICI rows (additive)

The user’s 11-cohort CLDN4-high ICI meta (**OR = 0.42**) is **taken as given**. These rows are **additive only** for a future expanded meta. They are **not** pooled into that OR and the original 11 cohorts are **not** re-cut.

## Extra NSCLC rows (do not mix SCLC)

Locked primary = within-cohort **median** split (ties at the median → low). Continuous = rank-biserial (responders minus non-responders). OR is high vs low × responder vs non-responder; zero cells get Haldane–Anscombe +0.5; p is Fisher’s exact; Woolf 95% CI.

DCB = PFS ≥ 6 months. Patients censored before 6 months are **excluded** from the DCB 2×2 (not scored as non-DCB).

| Series | Gene | Endpoint | n | Effect | Value | 95% CI | p |
|---|---|---|---:|---|---:|---|---:|
| GSE274975 NSCLC | CLDN4 | ORR (CR/PR vs SD/PD), median | 56 | OR high vs low | **1.00** | 0.30–3.35 | 1.00 |
| GSE274975 NSCLC | TACSTD2 | ORR (CR/PR vs SD/PD), median | 56 | OR high vs low | **1.00** | 0.30–3.35 | 1.00 |
| GSE274975 NSCLC | CLDN4 | ORR, continuous | 56 (14 R / 42 NR) | rank-biserial | −0.034 | — | 0.86 |
| GSE274975 NSCLC | TACSTD2 | ORR, continuous | 56 (14 R / 42 NR) | rank-biserial | +0.150 | — | 0.41 |
| GSE274975 NSCLC | CLDN4 | DCB (PFS ≥ 6 mo), median | 54 | OR high vs low | **1.00** | 0.33–3.02 | 1.00 |
| GSE274975 NSCLC | TACSTD2 | DCB (PFS ≥ 6 mo), median | 54 | OR high vs low | 0.73 | 0.24–2.20 | 0.78 |
| GSE274975 NSCLC | CLDN4 | PFS, Cox per 1 SD | 60 | HR | 1.02 | — | 0.89 |
| GSE274975 NSCLC | TACSTD2 | PFS, Cox per 1 SD | 60 | HR | 0.98 | — | 0.90 |
| GSE283829 NSCLC | CLDN4 | CR vs SD+PD, median | 27 | OR high vs low | **0.75** | 0.13–4.25 | 1.00 |
| GSE283829 NSCLC | TACSTD2 | CR vs SD+PD, median | 27 | OR high vs low | **0.75** | 0.13–4.25 | 1.00 |
| GSE283829 NSCLC | CLDN4 | CR vs PD, median | 17 | OR high vs low | **0.75** | 0.11–5.24 | 1.00 |
| GSE283829 NSCLC | TACSTD2 | CR vs PD, median | 17 | OR high vs low | **0.75** | 0.11–5.24 | 1.00 |

Tertile T3-vs-T1 and LUAD/LUSC splits are in `tables/all_tests.tsv` and `tables/extra_rows_NSCLC_histology.tsv`.

### What these series are

- **GSE274975** (Poddubskaya et al., *Front Immunol* 2024, PMID 39723204). Pretreatment FFPE bulk RNA-seq, STAR counts → `log2(CPM+1)`. RECIST and PFS from EuropePMC Table S1 (join on `OB_pat_LuC_*`). NSCLC with both RECIST and genes: 14 responders (3 CR + 11 PR) vs 42 SD/PD. Six NSCLC patients are DCB-ineligible (PFS < 6 mo and censored). One GEO-annotated SCLC (`Small cell lung cancer` source name) is **excluded** from this NSCLC table. Regimens are mixed PD-1/PD-L1 ± chemo ± CTLA-4 ± PARPi.
- **GSE283829** (Lindberg et al., *J Thorac Oncol* 2025). GEO `disease stage` is RECIST-like best response (CR 7 / SD 10 / PD 10). No PR labels, no PFS. Not in the original 11-cohort slide.

## Named leftover series the slide did not use

| Accession | What GEO actually is | Extra CLDN4 vs DCB/ORR row? |
|---|---|---|
| **GSE328294** | H23 KRAS-mutant NSCLC cell line, 48 h DMSO vs IBI351, n=6 FPKM | **No.** Cell-line treatment contrast only (below). |
| **GSE308745** | PBMC scRNA / scTCR / ADT, CD3/CD4 T cells, before neoadjuvant PD-1 (PB cohort2; PMID 41904165). RAW.tar is 10x MTX (~1.6 GB). | **No.** Not tumor bulk. GEO SOFT has no DCB/ORR field. Titles use RE/I/G prefixes; that is **not** treated as a response codebook. |
| **GSE302284** | scRNA of residual EGFR-mutant tumor/LN plus DFCI282/PC9 osimertinib vs vehicle (TROP2 CAR-T / DTP paper). | **No.** No ICI or ADC response labels. Treatment = TKI vs vehicle. |
| **GSE244944** | SCLC H3K27ac ChIP-seq, COR-L88 / NCI-H82 | **No.** Separate SCLC table. No patient ICB. |
| **GSE244945** | SCLC cell-line bulk RNA-seq | **No patient OR.** Separate SCLC table. IMpower133 RNA is **not** on GEO. |
| **GSE244946** | COR-L88 scRNA, n=3 | **No.** Separate SCLC table. |

### GSE328294 (not an ICI OR)

| Gene | Contrast | n | rank-biserial (IBI351 − DMSO) | exact MWU p |
|---|---|---:|---:|---:|
| CLDN4 | IBI351 vs DMSO | 3 vs 3 | −1.00 (down) | 0.10 |
| TACSTD2 | IBI351 vs DMSO | 3 vs 3 | +1.00 (up) | 0.10 |

Direction is a drug-on-cell-line effect. n=3 vs 3 cannot reject a null at α=0.05 even when ranks do not overlap.

## Extra SCLC table (separate; never mixed into the NSCLC OR)

Must-try GSE244944 / GSE244945 / GSE244946 **do not contain patient ICB response labels**. IMpower133 RNA (the clinical claim in that paper) is not on GEO.

| Series | Gene | Contrast | n | rank-biserial | p | Use |
|---|---|---|---:|---:|---:|---|
| GSE244945 H82 | CLDN4 | NOTCH1-dox vs parental | 3 vs 3 | +1.00 | 0.10 | mechanistic only |
| GSE244945 H82 | TACSTD2 | NOTCH1-dox vs parental | 3 vs 3 | +1.00 | 0.064 | mechanistic only |
| GSE244945 H69 | CLDN4 | NOTCH1-dox vs parental | 3 vs 3 | −1.00 | 0.10 | mechanistic only |
| GSE244945 H69 | TACSTD2 | NOTCH1-dox vs parental | 3 vs 3 | +1.00 | 0.077 | mechanistic only |
| GSE244945 H524 | CLDN4 | NOTCH1-dox vs parental | 3 vs 3 | −1.00 | 0.10 | mechanistic only |
| GSE244945 H524 | TACSTD2 | NOTCH1-dox vs parental | 3 vs 3 | +1.00 | 0.064 | mechanistic only |
| GSE244944 | — | ChIP-seq | 18 | — | — | no expression-vs-ICB test |
| GSE244946 | — | COR-L88 scRNA | 3 | — | — | not patient ICB |

## Other leftover GEO that was opened and is not an extra OR row

| Accession | Why it is not an extra CLDN4 vs DCB/ORR row |
|---|---|
| GSE309652 | Author R/NR is on GEO (24 R / 48 NR), but **CLDN4 and TACSTD2 are absent** from the 768-gene NanoString panel. |
| GSE309446 | TCR-seq of blood/FFPE after sitravatinib + tislelizumab + docetaxel. No gene-expression matrix and no per-patient DCB/ORR in GEO. |
| GSE202417 | Clariom D of **PBMC CD8 T cells** (nivo ± bezafibrate), not tumor bulk. |
| GSE207422 | Already in the original 11; not re-cut. |

## How to reproduce

```bash
python3 -m pip install -r results/rework/B5_extra_leftover/requirements.txt
python3 results/rework/B5_extra_leftover/analyze.py
```

Public GEO FTP + EuropePMC supplementary zip only. Raw downloads stay in `/tmp/b5_extra_raw` (not committed).
