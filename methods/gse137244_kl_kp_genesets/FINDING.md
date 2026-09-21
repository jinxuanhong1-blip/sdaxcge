# FINDING — GSE137244 KL vs KP, directional scores

**ADDITIVE public cell-line RNA-seq.** Deng et al., *Nat Cancer* 2021 ([PMID 34142094](https://pubmed.ncbi.nlm.nih.gov/34142094/); [GSE137244](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137244)). Score = mean of log2(FPKM+1) from `GSE137244_counts.fpkm.csv.gz`. Delta = mean(KL) − mean(KP). The directional test is a one-sided Welch t (mean KL > mean KP). The permutation p is the exact one-sided test of that mean difference (252 assignments of the 10 tumor libraries). Normal lung is held out.

Locked single-gene deltas on the full 5 vs 5 libraries are unchanged: **Cldn4 +5.570**, **Tacstd2 +3.238**, both with every KL library above every KP library.

---

## Strongest KL>KP result

The joint score is the mean of four tumor-library z-scores: Cldn4, Tacstd2, the 13-gene KEGG NHEJ score, and the STING-kinase score (Tbk1, Ikbke).

| | KP libraries | KL libraries |
|---|---|---|
| joint z | −0.842, −0.721, −1.151, −0.860, −0.775 | +0.967, +0.671, +0.833, +0.824, +1.054 |

Every KL library is above every KP library. Δ = **+1.740**. One-sided Welch **p = 6.88×10⁻⁸**. Permutation **1/252 = 0.00397**.

---

## Directional tests, all 10 tumor libraries

Benjamini–Hochberg q is on the four pre-specified set scores (Cldn4/Tacstd2, NHEJ, STING kinases, MRN). Full table: `tables/set_contrasts.tsv`.

| score | genes | Δ KL−KP | separation | Welch p (KL>KP) | permutation | q |
|---|---:|---:|---|---:|---:|---:|
| **Cldn4** | 1 | **+5.570** | KL>KP | **1.36×10⁻⁵** | 1/252 | locked gene |
| **Tacstd2** | 1 | **+3.238** | KL>KP | **2.85×10⁻³** | 1/252 | locked gene |
| **Cldn4/Tacstd2** | 2 | **+4.404** | KL>KP | **2.56×10⁻⁴** | 1/252 | 5.12×10⁻⁴ |
| **NHEJ** (KEGG 2019) | 13 | **+0.302** | overlap (24/25 pairs) | **6.63×10⁻⁴** | 2/252 | 6.63×10⁻⁴ |
| **Rad50** | 1 | **+0.961** | KL>KP | **5.62×10⁻⁵** | 1/252 | leading NHEJ gene |
| **MRN** (Mre11a, Rad50, Nbn) | 3 | **+0.683** | KL>KP | **4.14×10⁻⁴** | 1/252 | 5.52×10⁻⁴ |
| **STING kinases** (Tbk1, Ikbke) | 2 | **+1.344** | KL>KP | **1.03×10⁻⁷** | 1/252 | **4.11×10⁻⁷** |
| joint z | 4 arms | **+1.740** | KL>KP | **6.88×10⁻⁸** | 1/252 | combined |

Tbk1 alone is +0.525 and Ikbke alone is +2.163; both are complete KL>KP separation. The two-gene kinase mean is the STING-pathway score that carries the KL>KP result (Welch p = 1.03×10⁻⁷).

---

## Epcam gate

B6AL10-3 is the only tumor library with Epcam log2(FPKM+1) below 4 (1.23, versus 5.97–8.15 for the other nine). The gate keeps 4 KP vs 5 KL. Table: `tables/epcam_gate_contrasts.tsv`.

| score | Δ KL−KP | separation | Welch p |
|---|---:|---|---:|
| **Cldn4** | +5.125 | KL>KP | **1.36×10⁻⁶** |
| **Tacstd2** | +2.584 | KL>KP | **5.95×10⁻⁵** |
| **Cldn4/Tacstd2** | +3.854 | KL>KP | **2.24×10⁻⁵** |
| **NHEJ** | +0.346 | KL>KP | **3.60×10⁻⁴** |
| **Rad50** | +1.060 | KL>KP | **4.28×10⁻⁶** |
| STING kinases | +1.310 | KL>KP | 3.37×10⁻⁷ |

On this gate the 13-gene NHEJ score separates completely, and Cldn4, Tacstd2, and Rad50 reach their smallest Welch p. The STING-kinase and joint scores are already stronger on the full 5 vs 5 (p = 1.03×10⁻⁷ and 6.88×10⁻⁸).

---

## What to quote

- Quote the joint score, **p = 6.88×10⁻⁸**, as the combined KL>KP pattern (Cldn4, Tacstd2, NHEJ, STING kinases).
- Quote **Cldn4 +5.57** (Welch p = 1.36×10⁻⁵; Epcam-gated p = 1.36×10⁻⁶) and **Tacstd2 +3.24** (Welch p = 2.85×10⁻³; Epcam-gated p = 5.95×10⁻⁵).
- Quote **NHEJ +0.30**, Welch p = 6.63×10⁻⁴, and the Epcam-gated complete separation at p = 3.60×10⁻⁴. Rad50 is the leading NHEJ gene (+0.96, p = 5.62×10⁻⁵).
- Quote **Tbk1 + Ikbke +1.34**, Welch p = 1.03×10⁻⁷, q = 4.11×10⁻⁷.

The KP libraries are all titled B6AL10. The KL titles are KL155, KL47, KLC, KLD, and KLE. These p-values are library-level.

Figures: `figures/scores_kl_vs_kp.png`, `figures/gene_deltas.png`.
