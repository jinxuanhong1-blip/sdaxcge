# FINDING — pairwise CLDN4-only combo GSE127465 + GSE148071

**The two cohorts do not differ. CellChat-style outgoing was not run.**

Patient-level malignant CLDN4 vs T/NK is weakly **positive** in both public matrices and is not a CLDN4-high = T/NK-low story.

| cohort | n eligible / deposited | ρ [95% CI] | p |
|---|---:|---|---:|
| GSE127465 | **7 / 7** | **+0.214 [−0.642, +0.833]** | 0.645 |
| GSE148071 | **25 / 42** | **+0.189 [−0.223, +0.544]** | 0.365 |

Fisher-z difference *z* = 0.048, **p = 0.962**. Median-split T/NK is also not lower in CLDN4-high (GSE127465 rank-biserial *r* = 0, p = 1.0, 4 vs 3; GSE148071 *r* = +0.256, p = 0.289, 13 vs 12).

Additive. **CLDN4 only.** Not a triple/quad. No dual-high. Not a cell-level merge. Honest test n is **7** and **25**, not 40,362 or 89,887 cells.

---

## Gate

Both public processed matrices are **<2 GB**. Both have **CLDN4** and both have malignant/epithelial plus T/NK lineages. The stop-if-missing rule does **not** apply.

| Series | Processed expression | Bytes | <2 GB | CLDN4 | Malignant / epithelium | T/NK |
|---|---|---:|:---:|:---:|:---:|:---:|
| [GSE127465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465) | `GSE127465_human_counts_normalized_54773x41861.mtx.gz` | 528,303,938 | yes | present | author `PatientN-specific` in tumor | author `tT cells` + `tNK cells` |
| [GSE148071](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148071) | `GSE148071_RAW.tar` (42 `*_exp.txt.gz`) | 180,193,280 | yes | present | no GEO label; marker-argmax epithelium (putative) | marker-argmax T + NK |

GSE127465 mouse MTX (166.2 Mb) was not used. Blood (14,411 cells) was not used for malignant CLDN4 or T/NK. GSE127465 `RAW.tar` was not needed once the series MTX + metadata were present.

---

## Honest n

Patient is the unit. Do not write n = 54,773 or n = 89,887.

| Item | n | Note |
|---|---:|---|
| GSE127465 human patients deposited | **7** | p1–p7; Zilionis et al., *Immunity* 2019, PMID 30979687 |
| GSE127465 tumor / blood cells | 40,362 / 14,411 | blood excluded |
| GSE127465 own-patient `PatientN-specific` in tumor | **3,574** | primary malignant (p1 1,289; p2 98; p3 527; p4 172; p5 757; p6 214; p7 517) |
| GSE127465 tumor T / NK / T+NK | 12,776 / 1,116 / **13,892** | all 7 have ≥20 malignant and ≥20 T/NK |
| GSE127465 eligible patients | **7** | this is the GSE127465 test n |
| GSE148071 patients deposited | **42** | Wu et al., *Nat Commun* 2021, PMID 33953163; one biopsy each |
| GSE148071 cells in tar | 89,887 | not the test n |
| GSE148071 putative epithelial / T/NK cells | 51,213 / 4,161 | marker-argmax |
| GSE148071 eligible (≥20 epi **and** ≥20 T/NK) | **25** | this is the GSE148071 test n |
| GSE148071 excluded | **17** | almost no T/NK (P17 7,322 epi / 0 TNK; P3 7,107 / 17; P41 5,471 / 3; P16/P20/P33 0 TNK) |
| Pairwise combo (sum of eligible) | **32** | independent patients; not one object |
| Dual-high TACSTD2×CLDN4 | **0** | not defined |
| Extra series | **0** | pairwise only |
| CellChat-style outgoing | **0** | cohorts do not differ |

GSE127465 p2 is thin on NK (6 NK) but T/NK = 110 and malignant = 98, so it stays in the n=7 set. GSE148071 T/NK is thin in the eligible pool: P7 is 1,060 / 4,161 of all T/NK.

---

## Patient-level malignant CLDN4 vs T/NK

Computed from the public matrices. Not a taken-as-given ρ.

| cohort | assay | n deposited | n eligible | ρ [95% CI] | p |
|---|---|---:|---:|---|---:|
| GSE127465 | inDrops, author-normalized MTX; author malignant | 7 | 7 | +0.214 [−0.642, +0.833] | 0.645 |
| GSE148071 | 10x UMI; putative epithelium | 42 | 25 | +0.189 [−0.223, +0.544] | 0.365 |

Differ? **No.** Opposite-sign rule does not fire. Fisher-z two-sample p = 0.962. Neither |ρ| ≥ 0.30.

### Median split (descriptive)

| cohort | n high / low | median T/NK high | median T/NK low | rank-biserial *r* | p |
|---|---:|---:|---:|---:|---:|
| GSE127465 | 4 / 3 | 0.296 | 0.286 | 0.00 | 1.00 |
| GSE148071 | 13 / 12 | 0.106 | 0.047 | +0.256 | 0.289 |

Positive *r* = CLDN4-high has **higher** T/NK. GSE127465 n=7 makes the CI and the 4-vs-3 split uninformative except as a sign check.

Forest: `results/forest_CLDN4_tnk.png`. Per-patient: `results/gse127465_patients.tsv`, `results/gse148071_patients.tsv`. Machine combo: `results/combo_cldn4_tnk.tsv`.

---

## CellChat-style outgoing

**Not run.** Assignment: *If it differs, CellChat-style outgoing.* The two patient-level associations do not differ (same sign, Fisher-z p = 0.96, both compatible with zero). No ligand table was invented.

This is not a claim that ligand–receptor structure is identical. It is a claim that the assigned trigger for outgoing CellChat-style on this pair was not met.

---

## Methods

- Predictor: malignant **CLDN4 only** (mean log1p). TACSTD2 is not a gate.
- GSE127465 malignant = `Major cell type == Patient{N}-specific` inside that patient's **tumor**. Type I/II, club, ciliated, fibroblasts, and endothelium are not malignant. Matrix values are the author total-count-normalized MTX; CLDN4 score = mean log1p of that value in malignant cells. T/NK fraction = (tT + tNK) / tumor cells.
- GSE148071 epithelium / T / NK = marker-argmax on log1p(CP10k) modules (same markers as `methods/gse148071_cellchat_cldn4`). Epithelium is **putative** malignant (no GEO / CopyKAT label). T/NK fraction = (T + NK) / all cells in that biopsy.
- Eligible: ≥20 scored malignant/epithelial **and** ≥20 T/NK.
- Spearman on eligible patients. Fisher-z 95% CI. Difference = opposite sign (with |ρ| floor) **or** Fisher-z p < 0.05 **or** one |ρ| ≥ 0.30 and the other |ρ| < 0.15.
- CellChat-style outgoing (Hill / truncated mean, CellChatDB v2, 100 high/low permutations) only if they differ. CellChat R was not called.
- p-values on n=7 are descriptive.

---

## What this is not

- Not a CLDN4-cold-TME confirmation on this pair. Both point estimates go the other way and neither is significant.
- Not n=42 for GSE148071. Seventeen biopsies cannot be scored.
- Not n=7 as a precise effect size. The GSE127465 CI covers almost the full [−1, 1] range.
- Not a triple/quad (no GSE131907 / GSE205335 / GSE207422).
- Not dual-high TACSTD2×CLDN4.
- Not a merged Seurat/AnnData object (inDrops normalized MTX vs 10x UMI).
- Not CellChat R / LIANA / NicheNet.
- Not ICI response (no usable RECIST/MPR table in the files used here).

---

## 结论

GSE127465 与 GSE148071 的公开 processed 矩阵都 **<2 GB**，都有 CLDN4 和恶性/上皮 + T/NK，因此按任务做了配对，没有 stop。患者层面恶性 CLDN4 vs T/NK：**7 例** ρ=+0.21（p=0.65）和 **25/42 例** ρ=+0.19（p=0.37），方向相同，Fisher-z p=0.96，**不构成差异**。按任务门控：**不跑** CellChat-style outgoing。不能写成 n=42 或 n=54,773。不是 dual-high，不是三/四组合并。

Reproduce: `python3 methods/pair_127465_148071_cldn4/scripts/analyze.py`
