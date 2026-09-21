# FINDING — LINCS L1000 CLDN4 shRNA hunt, TACSTD2 CRISPR IFN / APM

Additive public evidence only. This folder does not touch private CLDN4 knockdown cocultures, and it does not revise locked CosMx, concordant-4, GSE137244, TCGA keratin, or TISMO numbers.

Level 5 values are moderated z-scores from the LINCS 2020 beta release. Positive z is up versus the plate control. Positive NES is enrichment at that up end. Engine is `gsea_core.py`: weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. BH-FDR is inside the two headline sets (Hallmark IFN-γ and custom MHC-I / APM), on each consensus rank separately. Hallmark IFN-α is secondary and is outside that BH family. Numbers below are written from `tables/`.

---

## 一句话 / TL;DR

CLDN4 在 LINCS 2020、phase I、phase II 里没有 shRNA / CRISPR / 过表达签名，做不了敲低富集。TACSTD2 没有 shRNA；有 20 个 CRISPR（10 个细胞系 × 2 个 guide，96 h）。细胞系中位数秩上，IFN-γ NES **+1.818**（FDR 0.002），IFN-α NES **+1.804**（nom p 0.001），MHC-I / APM NES **+0.885**（FDR 0.393）。靶基因 TACSTD2 的中位 z 是 **0**，同一细胞系的两个 guide 不相关（中位 Spearman **−0.058**）。Broad `is_hiq` 的 3 个签名把 IFN-γ NES 翻成 **−0.885**。

| Question | Honest n | Result |
|---|---|---|
| CLDN4 shRNA, CRISPR, OE, or siRNA | 0 signatures | No enrichment test |
| TACSTD2 shRNA (`trt_sh` / `.cgs` / `.css`) | 0 signatures | No shRNA enrichment test |
| TACSTD2 CRISPR cell-median rank | 10 lines, 20 signatures, 12,326 genes | IFN-γ NES +1.818 FDR 0.002; IFN-α NES +1.804 nom p 0.001; MHC-I / APM NES +0.885 FDR 0.393 |

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Catalog | LINCS 2020 `siginfo_beta.txt` (1,201,944 rows). Phase I GSE92742 sig info (473,647) and pert info (51,383). Phase II GSE70138 sig info 2017-03-06 (118,050). |
| CLDN4 match | Exact `cmap_name` (2020) or `pert_iname` (GEO). Substring `CLDN4` in the 2020 file is also 0. |
| TACSTD2 matrix | CRISPR Level 5 gctx, HTTP range slice. 175,024,881 bytes read from the 6,516,812,529-byte file. |
| Included signatures | All 20 `qc_pass=1`, `is_exemplar_sig=1`, `is_ncs_sig=1`, `is_null_sig=0`. Time is 96 h. Guides `BRDN0001498584` (well L07) and `BRDN0001498879` (well D23). |
| Cell lines | A375, A549, AGS, BICR6, ES2, HT29, MCF7, PC3, U251MG, YAPC. A549 is the only lung line (NSCLC). |
| Primary rank | Mean of the two guides inside each cell line, then median across 10 lines. TACSTD2 removed from the rank. |
| Headline | Hallmark IFN-γ (181/200 in the L1000 universe), MHC-I / APM (20/21). |
| Secondary | Hallmark IFN-α (84/97). |
| Aliases | `WARS1→WARS`, `MARCHF1→MARCH1`. |
| Duplicate symbol | MIA2 Entrez 117153 (inferred) dropped; Entrez 4253 (best inferred) kept. Not in either set. |

Series / build: LINCS 2020 beta, clue.io S3 `builds/LINCS2020`. Assay paper: Subramanian et al., *Cell* 2017. Data use: https://clue.io/connectopedia/data_use_policy .

---

## CLDN4 is absent, so there is no knockdown enrichment

| Source | Rows | CLDN4 | TACSTD2 |
|---|---:|---:|---:|
| LINCS 2020 `siginfo_beta` | 1,201,944 | 0 | 20 |
| GSE92742 signature info | 473,647 | 0 | 0 |
| GSE92742 perturbagen info | 51,383 | 0 | 0 |
| GSE70138 signature info (2017-03-06) | 118,050 | 0 | 0 |

LINCS 2020 class sizes, exact `cmap_name` (from `tables/inventory_pert_type.tsv`):

| Class | Signatures | Unique names | CLDN4 | TACSTD2 |
|---|---:|---:|---:|---:|
| `trt_sh` (shRNA) | 177,263 | 4,917 | 0 | 0 |
| `trt_sh.cgs` | 36,720 | 4,345 | 0 | 0 |
| `trt_sh.css` | 24,368 | 3,807 | 0 | 0 |
| `trt_xpr` (CRISPR) | 140,945 | 5,158 | 0 | 20 |
| `trt_oe` | 34,171 | 4,040 | 0 | 0 |
| `trt_si` | 162 | 23 | 0 | 0 |

Nearby genes that do exist, so the symbol field is not blank for claudins: CRISPR for CLDN1, CLDN2, CLDN11, CLDND1; overexpression for CLDN18; shRNA and CRISPR for EPCAM and MUC1. CLDN4 itself was not a LINCS perturbagen. The shRNA gctx was not downloaded.

---

## TACSTD2 CRISPR headline

Primary contrast: cell-median z, **10** lines. Rank length **12,326** after dropping TACSTD2 and the extra MIA2 row.

| Set | Role | NES | nom p | BH FDR | n in rank | mean z | leading-edge n |
|---|---|---:|---:|---:|---:|---:|---:|
| Hallmark IFN-γ | headline | +1.818 | 0.001 | 0.002 | 181 | +0.149 | 93 |
| Hallmark IFN-α | secondary | +1.804 | 0.001 | — | 84 | +0.170 | 35 |
| MHC-I / APM | headline | +0.885 | 0.393 | 0.393 | 20 | +0.055 | 6 |

Nominal p = 0.001 is the permutation floor: 0 of 1000 nulls reached the observed ES, so p = 1/1001. BH on the two headline p-values (0.001 and 0.393) gives FDR 0.002 and 0.393. IFN-α was not in that family.

IFN-γ leading edge starts IRF7, MX1, RIPK1, IL15, PTGS2, PLA2G4A, ISG20, UPP1, ARID5B, RSAD2 (first 25 in `tables/gsea_headline.tsv`). MHC-I / APM leading edge is 6 genes (TAPBPL, IRF1, PSMB10, ERAP1, TAP1, PSMB9) and the set NES is null.

### Sensitivities (same engine, own BH)

| Rank | n sig / n cell | IFN-γ NES (FDR) | IFN-α NES (nom p) | MHC-I / APM NES (FDR) |
|---|---|---|---|---|
| Median of 20 signatures | 20 / 10 | +2.106 (0.002) | +2.001 (0.001) | +0.938 (0.331) |
| `cc_q75 ≥ 0.2` cell-median | 10 / 8 | +1.471 (0.006) | +1.545 (0.011) | +0.687 (0.670) |
| Broad `is_hiq=1` cell-median | 3 / 3 | −0.885 (0.485) | −0.793 (0.506) | +0.954 (0.462) |
| Primary z, landmark genes only | 20 / 10 | +1.603 (0.010) | +1.019 (0.207) | skipped (2 genes) |

`cc_q75 ≥ 0.2` drops MCF7 and PC3 (neither guide clears 0.2). The three `is_hiq` signatures are ES2 L07, BICR6 L07, and A549 L07. Their joint rank does not keep the positive IFN NES.

Per-signature NES is scattered (figure, right). IFN-γ is positive in 13/20 signatures (median NES +1.163; Wilcoxon p = 0.044 on signatures). Averaged to the cell line, 8/10 lines have mean NES > 0, median NES falls to +0.240, Wilcoxon p = 0.105. MHC-I / APM is 11/20 positive, cell-line split 5/5, Wilcoxon p = 0.432. Those Wilcoxon tests describe sign counts. The primary number is the single cell-median GSEA.

---

## On-target z and guide agreement

TACSTD2 is best inferred, not a landmark. Its Level 5 z is the on-target readout available in this matrix.

| | Median z | Signatures with z < 0 |
|---|---:|---:|
| TACSTD2 | 0.000 | 9 / 20 |
| CLDN4 (not the perturbed gene) | +0.042 | 10 / 20 |

A549, the only lung line, has TACSTD2 z **+0.462** (L07, `is_hiq`) and **+0.503** (D23). Both are up. BICR6 is 0.000 on both guides. Within a cell line, Spearman correlation of the two full z vectors has median **−0.058** (range −0.274 to +0.194). The two reagents do not share a signature.

A549 guide-level NES, from `tables/gsea_per_signature.tsv`:

| Guide | hiq | cc_q75 | TAS | IFN-γ NES (nom p) | MHC-I / APM NES (nom p) |
|---|---|---:|---:|---|---|
| L07 `BRDN0001498584` | 1 | 0.256 | 0.140 | −1.350 (0.017) | −1.735 (0.003) |
| D23 `BRDN0001498879` | 0 | 0.300 | 0.195 | +1.622 (0.002) | −0.926 (0.211) |

Median TAS across the 20 signatures is 0.147. Median `cc_q75` is 0.210. BICR6 was run at n = 2 replicates; the other lines at n = 3.

---

## L1000 coverage

| Set | In universe | Landmark | Best inferred | Inferred | Absent |
|---|---:|---:|---:|---:|---|
| Hallmark IFN-γ (200) | 181 | 33 | 138 | 10 | 19, including `CD274` and `NLRC5` |
| Hallmark IFN-α (97) | 84 | 8 | 74 | 2 | 13 |
| MHC-I / APM (21) | 20 | 2 (`PSMB8`, `PSMB10`) | 17 | 1 (`HLA-G`) | `NLRC5` |

Most of the IFN-γ set is inferred. The landmark-only sensitivity still gives IFN-γ NES +1.603 (FDR 0.010) on 33 genes. MHC-I / APM has two landmarks, under the size floor of 8, so that sensitivity is not a GSEA. `CD274` and `NLRC5` are absent from the 12,328-gene space, same hole as on U133A for those two symbols.

---

## What this page can carry

The public L1000 catalog has no CLDN4 knockdown, knockout, or overexpression signature, so it supplies no CLDN4 IFN / APM enrichment number.

TACSTD2 has a CRISPR series. On the pre-specified cell-median rank, Hallmark IFN-γ and IFN-α sit at the up end and MHC-I / APM does not. The same matrix does not show TACSTD2 itself going down, the two guides in a cell line do not correlate, and the three Broad high-quality signatures do not keep the positive IFN NES. A549's high-quality guide is negative for both IFN-γ and MHC-I / APM. This is a cell-line L1000 result (one NSCLC line among nine other lineages), not a tumour-cell IFN call and not an ICI series.
