# FINDING — public mouse KL vs KP/K Cldn4 contrasts

**ADDITIVE public MOUSE only. Cldn4-only. No private 8-KL matrices. No Tacstd2∩Cldn4 dual-high.** Thesis is already taken as given: Cldn4-high epithelium is barrier / immune-cold; KL (Stk11/Lkb1 loss) is ICI-resistant relative to KP/K. This folder does **not** audit that thesis as claim-failed. It scores **separate** combinatorial contrasts (not one mega-merge) on four public processed matrices.

Primary endpoints, when the same biological unit can support them:

1. epithelial / malignant **Cldn4** (mean or %pos)
2. **T/NK** fraction or bulk T/NK score
3. epithelial **IFN / MHC** if the same cells exist

Honest n = **mice or tumors**, never cells.

---

## Contrast summary (KL minus KP / KL minus K)

Full machine table: `tables/contrast_summary.tsv`. Headline rows below use the honest unit.

| accession | contrast | unit | n | Cldn4 Δ | T/NK Δ | IFN Δ | MHC Δ | verdict |
|---|---|---|---|---|---|---|---|---|
| **GSE6135** | KL − KP | **7 vs 5 mice** | 7 / 5 | **+0.224** (p_mwu=0.88) | **−0.658** (p=0.030) | +0.590 (p=0.030) | −0.278 (p=0.64) | **only public KL vs KP with Cldn4 + T/NK**. Cldn4 up (ns). T/NK down vs KP. |
| **GSE6135** | KL − K | **7 vs 5 mice** | 7 / 5 | **+0.505** (p=0.64) | −0.052 (p=0.88) | +0.734 (p=0.030) | −0.725 (p_w=0.047 / p_mwu=0.11) | Cldn4 up (ns). T/NK ≈ K. MHC down vs K (trend). |
| **GSE154989** | KL − KP / KL − K | — | 0 / 15 / 9 | — | — | — | — | **no-go**: K/KP epithelium only; no KL arm. T/NK design no-go (CD45− FACS). |
| **GSE154989** | KP − K (allowed) | **15 vs 9 mice** | — | **+1.73** (p=0.020) | no-go | −0.012 (p=0.81) | **−0.374** (p=0.0005) | Not a KL contrast. KP epithelium is Cldn4-higher and MHC-lower than K. |
| **GSE179502** | KL − KP | — | 3 / 0 | — | — | — | — | **no-go**: KT;Lkb1 XTR only; no p53 / KP arm. |
| **GSE179502** | NonRestored − Restored ≈ KL − K | **3 vs 3 mice** | 3 / 3 | **+0.157** (p_w=0.036; p_mwu=0.10) | no-go | +0.004 (p=0.40) | −0.085 (p=0.40) | Cldn4 complete rank separation (KL-like > Restored). T/NK unscorable. IFN/MHC flat. |
| **GSE267321** | KL − KP | — | 2 / 0 | — | — | — | — | **no-go**: K / KK / KLK only; no KP. |
| **GSE267321** | KLK − K | **2 vs 2 tumors** | 2 / 2 | **empty** (6/7956 cells >0) | **−0.299** (0.097 vs 0.397; p_mwu=0.33 floor) | leftover epi; underpowered | leftover epi; underpowered | Cldn4 no-go (floor). T/NK both KLK tumors < both K tumors. KLK is STK11+KEAP1, not KL-only. |

Δ is always **resistant-like minus sensitive-like**. Positive Cldn4 = higher on the KL / KL-like arm. Negative T/NK = colder on that arm. p values are two-sided Mann–Whitney unless marked Welch.

Do not stack these four series into one n. Do not quote cell counts as n.

---

## GSE6135 — only public KL bulk with Cldn4

Ji et al., *Nature* 2007 ([PMID 17676035](https://pubmed.ncbi.nlm.nih.gov/17676035/); [GSE6135](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6135)). Affymetrix Mouse 430A 2.0 (GPL8321) series matrix. Human A549/H2126 arrays on GPL5373 were **not** used.

**KL** = Lkb1 **L/L or L/−** primary lung tumors. **Excluded** from primary KL: L/+ heterozygotes (mouse 392) and the 592 metastasis. Ink4a is catalogued, not a primary contrast.

| arm | mice | tumors (primary) | histology |
|---|---:|---:|---|
| KL | **7** | 10 | Ad (4 mice), Sq (2), Ad-sq (1) |
| KP | **5** | 5 | Ad |
| K | **5** | 5 | Ad |
| KL_het (excluded) | 1 | 2 | Ad |
| K_Ink4a (not primary) | 2 | 2 | Ad |

Cldn4 = probe-collapsed log2 intensity. T/NK and IFN/MHC = mean z-score of present genes across the 25 mouse arrays (gene lists in `gene_sets.py` / `tables/genes_used.tsv`).

### Mouse means (primary tumors only)

| mouse | arm | n_primary | hist | Cldn4 | T/NK | IFN | MHC |
|---|---|---:|---|---:|---:|---:|---:|
| KL_113 | KL | 1 | Sq | 7.29 | −0.090 | 0.383 | −1.063 |
| KL_452 | KL | 2 | Ad-sq | 7.04 | −0.425 | 1.229 | −1.085 |
| KL_540 | KL | 2 | Sq | 7.29 | −0.607 | 1.054 | −0.309 |
| KL_459 | KL | 1 | Ad | 5.41 | 0.446 | 0.078 | 0.277 |
| KL_547 | KL | 1 | Ad | 4.97 | −0.266 | 0.059 | 0.238 |
| KL_592 | KL | 2 | Ad | 6.19 | −0.548 | −0.216 | −0.295 |
| KL_861 | KL | 1 | Ad | 4.65 | −0.007 | −0.223 | 1.021 |
| KP_186 | KP | 1 | Ad | 5.76 | 0.806 | −0.131 | 0.182 |
| KP_196 | KP | 1 | Ad | 5.71 | 0.627 | −0.199 | 0.005 |
| KP_197 | KP | 1 | Ad | 5.78 | 0.775 | −0.234 | −0.090 |
| KP_498 | KP | 1 | Ad | 6.30 | −0.219 | −0.296 | 1.029 |
| KP_500 | KP | 1 | Ad | 5.92 | 0.230 | −0.399 | −0.606 |
| K_268 | K | 1 | Ad | 5.15 | −0.539 | −0.721 | 0.453 |
| K_287 | K | 1 | Ad | 5.81 | 0.210 | −0.477 | 0.757 |
| K_405 | K | 1 | Ad | 5.17 | 0.883 | −0.076 | 0.502 |
| K_484 | K | 1 | Ad | 5.49 | −0.773 | −0.155 | 0.813 |
| K_498 | K | 1 | Ad | 6.45 | −0.591 | −0.555 | 0.231 |

### What holds on this series

- **KL vs KP T/NK is the clean combinatorial immune contrast.** All 5 KP mice except KP_498 sit above the KL cloud. Δ = −0.66, MWU p=0.030, rank-biserial −0.77. This is the ICI-relevant direction (KL colder than KP).
- **Cldn4 is directionally higher in KL** vs KP and vs K. It is **not** significant at 7 vs 5. Sq / Ad-sq KL mice carry the high tail (Cldn4 ~7.0–7.3). Ad-only KL mice (592, 861, 459, 547) sit at 4.65–6.19, overlapping KP/K. Report the direction; do not over-call a Cldn4 genotype law.
- Combined IFN/MHC is flat. Split scores: bulk IFN is **higher** in KL (Sq/Ad-sq driven); MHC is **lower** vs K. This is bulk tumor, not purified epithelium. Do not read bulk IFN as epithelial IFN-cold.

Table: `tables/GSE6135.tsv` (mice) and `tables/GSE6135_tumors.tsv`.

---

## GSE154989 — K/KP Smart-seq2 epithelium (KL no-go)

Marjanovic et al., *Cancer Cell* 2020 ([PMID 32707077](https://pubmed.ncbi.nlm.nih.gov/32707077/); [GSE154989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154989)). Public TPM COO `GSE154989_mmLungPlate_fQC_dSp_normTPM.h5`. FACS `tdTomato+ / CD45− / CD11b− / TER119− / CD31−`.

**KL vs KP / KL vs K: no-go** (no Stk11-loss arm).

**T/NK: design no-go.** Residual Cd3d / Nkg7 / Ptprc is leak, not a compartment. Do not invent a T/NK fraction.

Allowed, because the user said we may:

| contrast | n | result |
|---|---|---|
| K vs KP Cldn4 mean | 9 K vs 15 KP mice (≥20 cells) | KP higher: 3.73 vs 2.01 TPM; Δ=+1.73; p=0.020 |
| K vs KP Cldn4 %pos | same | 0.586 vs 0.315; p=0.0035 |
| K vs KP IFN | same | flat (Δ=−0.012, p=0.81) |
| K vs KP MHC | same | KP lower: 0.923 vs 1.298; p=0.0005 |
| KP-only Cldn4 vs IFN | 15 KP | Spearman ρ=+0.49, p=0.062 (**not down**) |
| KP-only Cldn4 vs MHC | 15 KP | Spearman ρ=−0.50, p=0.058 (negative trend on the compact MHC set) |

Biological mouse = deposited `mouseID` with trailing `_T#` stripped. Dropped `KP_2w_ND_m1` (4 cells). T (normal AT2) excluded from K vs KP. Do not mix T AT2 into a Cldn4–MHC claim (they are Cldn4-low / MHC-high).

Table: `tables/GSE154989.tsv`.

---

## GSE179502 — KT;Lkb1 XTR neoplastic epithelium

Murray / Winslow, *Nat Commun* 2022 ([PMID 35228570](https://pubmed.ncbi.nlm.nih.gov/35228570/); [GSE179502](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179502)). Public Cell Ranger aggr mtx. Author cohort: **NonRestored** (Lkb1 off ≈ KL) vs **Restored** (Lkb1 on ≈ K). Score = mean log1p(CP10k) on QC cells (UMI≥500, genes≥200, mito≤20%). Honest n = **3 vs 3 mice**. MWU cannot beat p=0.1.

**KL vs KP: no-go** (no Trp53 / KP arm).

**T/NK: design no-go** (FACS-sorted neoplastic epithelium).

| mouse | cohort | arm | n_QC | Cldn4 | Cldn4 %pos | Stk11 | IFN | MHC |
|---|---|---|---:|---:|---:|---:|---:|---:|
| CM0875 | NonRestored | KL | 2522 | 0.261 | 0.249 | 0.005 | 0.074 | 0.631 |
| ZR1932 | NonRestored | KL | 972 | 0.340 | 0.302 | 0.013 | 0.075 | 0.835 |
| ZR1966 | NonRestored | KL | 2889 | 0.300 | 0.298 | 0.004 | 0.064 | 0.657 |
| CM0879 | Restored | K-like | 3322 | 0.134 | 0.114 | 0.294 | 0.074 | 0.711 |
| CM0884 | Restored | K-like | 836 | 0.215 | 0.225 | 0.311 | 0.063 | 0.864 |
| ZR1969 | Restored | K-like | 2102 | 0.080 | 0.078 | 0.309 | 0.064 | 0.802 |

Cldn4: every NonRestored mouse > every Restored mouse (Δ=+0.157, Welch p=0.036, MWU p=0.10, rank-biserial +1). Stk11 restore check: 0.007 vs 0.305. IFN/MHC do not separate. ZR1966 is restorable but vehicle-treated — author-labeled NonRestored; Stk11 stays on the floor.

Table: `tables/GSE179502.tsv`.

---

## GSE267321 — LKR13 non-malignant K / KK / KLK

Qian / Skoulidis / Heymach ([GSE267321](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267321)). Public normalized CSV. Title is **non-malignant** cells. Genotypes: **K**, **KK** (KEAP1), **KLK** (STK11+KEAP1). There is **no LKR13-KL** and **no KP**.

**KL vs KP: no-go.**

**Cldn4: no-go.** Row exists; **6 / 7956 cells > 0** (0.075%). Floor, not a genotype test.

**T/NK: usable, n=2 vs 2 tumors.** T/NK = Cd3d/e/g or Cd8a or Nkg7/Ncr1/Klrb1c.

| tumor | geno | n cells | n T/NK | frac T/NK | n Cldn4+ | n leftover epi |
|---|---|---:|---:|---:|---:|---:|
| LKR13-K-1 | K | 1767 | 800 | **0.453** | 2 | 129 |
| LKR13-K-2 | K | 1323 | 451 | **0.341** | 1 | 71 |
| LKR13-KK-1 | KK | 2172 | 392 | 0.180 | 1 | 70 |
| LKR13-KK-2 | KK | 681 | 48 | 0.070 | 0 | 13 |
| LKR13-KLK-1 | KLK | 1232 | 177 | **0.144** | 0 | 58 |
| LKR13-KLK-2 | KLK | 781 | 40 | **0.051** | 2 | 5 |

KLK minus K T/NK Δ = **−0.299**. Both KLK tumors sit below both K tumors. MWU p floor at 2 vs 2 is **1/3**. IFN/MHC on leftover marker epithelium is underpowered (KLK-2 has 5 epi cells) and is not a claim.

Table: `tables/GSE267321.tsv`.

---

## How to read this (locked rules)

- **Combinatorial, not merged.** Four contrasts, four n’s. No cross-series z-score mega-merge.
- **GSE6135 is the only series that can score KL vs KP for both Cldn4 and T/NK.** Cldn4 direction matches barrier-high; T/NK vs KP matches immune-cold. n=7 vs 5 is honest and thin. Histology is mixed in KL.
- **GSE179502** is the only epithelial Cldn4 KL-like vs Lkb1-on contrast. Cldn4 goes down when Lkb1 is restored. No T/NK. No KP.
- **GSE154989** is K vs KP epithelium only. KP is Cldn4-higher than K. T/NK cannot be scored. Within-KP IFN is not down.
- **GSE267321** is a non-malignant digest. Cldn4 is empty. T/NK is lower in KLK than K at 2 vs 2. KLK ≠ KL.
- Cldn4-only. Tacstd2 is not a gate. Private 8-KL matrices were not opened. GSE76628, GSE137669, GSE127465 mouse, and GSE274477 were not reopened.
- Reproduce: `python3 methods/public_kl_vs_kp_cldn4/analyze.py` (downloads GEO processed files to `/tmp/geo`).

## Files

| path | what |
|---|---|
| `tables/contrast_summary.tsv` | KL−KP / KL−K (and allowed proxies) across series |
| `tables/honest_n.tsv` | which series can contrast what |
| `tables/GSE6135.tsv` | mouse units |
| `tables/GSE154989.tsv` | mouse units |
| `tables/GSE179502.tsv` | mouse units |
| `tables/GSE267321.tsv` | tumor units |
| `tables/genes_used.tsv` | genes present per series |
| `analyze.py` / `gene_sets.py` | locked methods |
