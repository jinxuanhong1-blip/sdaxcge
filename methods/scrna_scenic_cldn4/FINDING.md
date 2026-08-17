# FINDING — CLDN4-high malignant regulons (GSE207422)

**Verdict: descriptive table only. No new ELF3–CLDN4 claim.** Lightweight GRN proxy (Pearson + public priors + AUCell) on public **GSE207422** malignant-like cells. **A10 ELF3–CLDN4 bulk RNA is taken as given** and is not re-tested. **CLDN4-high regulons only** — CLDN4-low edges were used solely to mark high-specific targets.

Full pySCENIC cisTarget was **not** run. No ChIP peaks were invented. Author CopyKAT IDs are **not** on GEO.

Primary table: [`results/scrna_scenic_cldn4/regulons.tsv`](../../results/scrna_scenic_cldn4/regulons.tsv).

---

## Honest n

| Item | n | Note |
| --- | ---: | --- |
| Cells (public UMI) | 92,330 | Hu et al. *Genome Med* 2023, PMID 36869384 |
| Epithelial (marker argmax) | 13,043 | same gate as `methods/scrna_scenic` |
| Malignant-like | **9,782** | epithelial AND normal-lung ≤ p75; **not CopyKAT** |
| Samples with ≥20 malignant-like | **12** | 6 LUAD + 6 LUSC; primary GRN samples |
| CLDN4-high (within-sample median) | **5,599** | cells |
| CLDN4-low (same samples) | 4,169 | cells |
| Samples with ≥20 high **and** ≥20 low | **5** | unit of prior-AUCell contrast; **UNDERPOWERED** (<6) |
| LUAD CLDN4-high | 1,633 / 6 samples | sensitivity |
| TFs screened in CLDN4-high | 861 | human TF list ∩ detection ≥5% |

**Split caveat (do not hide).** Within-sample median on a zero-inflated gene: **7/12 samples have n_low = 0** because the sample CLDN4 median is 0, so every malignant-like cell is labeled high. Those samples (BD_immune04/06/08/10/12/13/15) contribute no low arm. The five paired samples are BD_immune01, 03, 05, 07, 09.

**Cell-count imbalance.** BD_immune07 (LUSC, post) is 2,242 / 5,599 high cells (40%). Primary histology mix among the 12 samples: 7,414 LUSC vs 2,354 LUAD malignant-like cells. **Sample n is the claim.** Cell-level ρ is not reported as a patient-level result.

---

## Regulon table (CLDN4-high only)

CLDN4 is **held out** of every regulon (it is the split gene). Public priors list CLDN4 as a target of **GRHL2** and **SP1** (CollecTRI), **not ELF3**. `elf3_cldn4_in_prior = false`.

### Focus TFs — Pearson neighborhoods inferred in CLDN4-high cells

| TF | given? | prior n | Pearson n | mean r (top 50) | n high-specific (all genes) | high-specific in top 50 | prior ∩ Pearson top |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| **ELF3** | **yes (A10)** | 26 | 50 | 0.366 | 177 | **0** | CLDN7 |
| GRHL1 | no | 1 | 50 | 0.309 | 311 | 4 | ∅ |
| GRHL2 | no | 10 | 50 | 0.363 | 91 | **0** | ∅ |
| KLF4 | no | 112 | 50 | **0.551** | 86 | **0** | ATF3; HSPA1A; KRT19 |
| KLF5 | no | 75 | 50 | **0.518** | 114 | **0** | ANP32B; JUP; KLF4; MYH14 |
| OVOL1 | no | 9 | 50 | 0.275 | 313 | **42** | GADD45A |
| OVOL2 | no | 1 | 50 | 0.092 | 0 | 0 | ∅ |
| TFAP2A | no | 272 | 50 | 0.299 | 387 | 1 | VEGFA |
| SP1 | no | 1172 | 50 | 0.121 | 0 | 0 | 11 genes (weak r) |
| TP63 | no | 185 | 50 | **0.527** | 55 | **0** | CDKN2A; CLDN1; IGFBP2; JAG1; KRT5; PERP; SERPINB5; SFN; TCF4 |

High-specific = Pearson r_high ≥ 0.15 **and** r_low < 0.10. CLDN4-low regulons are **not** published.

### Pearson top targets (CLDN4-high; first 15)

| TF | top 15 targets |
| --- | --- |
| ELF3 *(given)* | TACSTD2; CLDN1; KLF5; HSPB1; ATF3; ANXA1; PERP; KLF4; KRT19; HES1; RAP2B; TIPARP; NET1; DSTN; RND3 |
| GRHL1 | JUP; CSTB; LRRC8A; CSTA; EMP1; PERP; SDC1; NECTIN4; RHCG; ANXA1; KRT16; RAP2B; DSP; DSC2; SLC20A2 |
| GRHL2 | ALDH3A1; NMRAL2P; TFRC; SDC1; CBR1; EIF4A2; NUCKS1; ABCC5; CALB2; CDKN2A; SOX2; ATP2B1; HSPA1B; SET; GPX2 |
| KLF4 | RAP2B; SERPINB5; ANXA1; CLDN1; RND3; MAL2; JUN; ATF3; PERP; HSPA1B; SFN; KLF5; TBL1XR1; LRRC8A; PMAIP1 |
| KLF5 | PERP; SERPINB5; RAP2B; ANXA1; KLF4; SFN; CLDN1; HSPB1; JUP; KRT19; MAL2; HSPA1B; SDC1; RND3; DSTN |
| OVOL1 | EMP1; KRT16; TACSTD2; KLK10; DSC2; LRRC8A; JUP; CSTB; LYPD3; MAL2; ANXA1; GRHL3; RHCG; RAB11FIP1; CD24 |
| OVOL2 | COL5A3; GFI1B; SLC25A39; SNRPD3; SHQ1; RPL15; SOD1; NKX2-8; WDR43; RHEB; SNU13; SNRPD2; EEF2; CUTA; SUMO2 |
| TFAP2A | EZR; ANXA2; PFKP; LDHA; SERINC2; SPINT1; PKM; PPP1R14B; F3; SLC2A1; SLC7A5; ITGB1; KRT81; GNG4; DSP |
| SP1 | RAB11FIP1; ACTN4; RAP2B; DVL3; CMIP; EIF4G1; WASL; MCL1; CLDN1; TBL1XR1; IRF2BP2; BPTF; MAPK6; DSC2; RARG |
| TP63 | SOX2; NTRK2; TFRC; GPX2; ABCC5; NMRAL2P; ALDH3A1; CDKN2A; TCF4; EIF4A2; PKP1; AKR1C2; KRT5; ATP1B3; ALDH1A1 |

Full membership: `regulons.tsv` / `regulons_focus.tsv`. Edge-level r: `pearson_top_edges_cldn4_high.tsv`.

---

## What this supports / does not

**Supports (descriptive, co-expression only):**

- In CLDN4-high malignant-like cells, **KLF4 / KLF5 / TP63 / ELF3 / GRHL2** have the tightest focus Pearson neighborhoods (mean r 0.36–0.55). The genes are an epithelial / tight-junction / keratin set (CLDN1, PERP, KRT19, JUP, SDC1, ANXA1, TACSTD2). That is the same program `methods/scrna_scenic` already saw on *all* malignant-like cells — now restricted to the CLDN4-high arm.
- **OVOL1** is the only focus TF whose top-50 is mostly high-specific (42/50): a squamous/keratin-leaning module (KRT16, KLK10, IVL, GRHL3, NECTIN4, TACSTD2).
- ELF3 Pearson top includes TACSTD2, CLDN1, KLF4, KLF5, KRT19. Public ELF3 prior ∩ Pearson top = **CLDN7 only**. Priors do **not** list ELF3→CLDN4.

**Does not support:**

- Binding of any TF at CLDN4 (no ChIP, no cisTarget).
- A **new** ELF3–CLDN4 discovery. ELF3 is labeled `elf3_given`. The bulk couple is A10, taken as given.
- “CLDN4-high has a *private* GRN.” For ELF3 / KLF4 / KLF5 / GRHL2 / TP63, **0 / 50** top targets meet the high-specific rule. Those neighborhoods also exist in CLDN4-low. Shared epithelial program, not a high-exclusive regulon.
- Public-prior **AUCell** being higher in CLDN4-high. Paired sample n = **5** (UNDERPOWERED). ELF3_prior Δ = +0.0004, p = 0.63. GRHL2 / KLF5 prior p = 0.0625 (unadjusted). ELF3 **RNA** Δ = +0.25, p = 0.0625 — direction matches A10, **not a discovery**.
- Data-driven screen hits ranked by n high-specific (BTF3, NPM1, ENO1, AES, …). Those are highly detected general factors, not a CLDN4 GRN.

**LUAD-only sensitivity (1,633 high cells / 6 samples):** ELF3 Pearson mean r = 0.41, but the top list is ribosomal + KRT7 / CLDN7 / TACSTD2. Small-n dropout. Do not treat LUAD neighborhoods as a second GRN.

---

## Method (locked)

- Dataset: GSE207422 public processed UMI (175 MB gzip). GSE253013 / GSE131907 not used.
- Malignant-like: marker epithelial AND normal-lung score ≤ epithelial 75th percentile.
- CLDN4-high: ≥ **within-sample** median of CLDN4 `log1p(CP10k)` among malignant-like; sample kept if ≥20 such cells.
- GRN: Pearson of 861 TFs vs 6,827 genes **in CLDN4-high cells only**; regulon = top 50 with r ≥ 0.10 (or top positive). Public priors = TRRUST v2 + DoRothEA + CollecTRI for the focus/control TFs.
- AUCell: Aibar-style recovery on **public priors** (external; not circular). Primary unit for high vs low = sample-paired means.
- pySCENIC: not importable in this run; cisTarget motif DBs were not downloaded either way.

---

## Files

| file | content |
| --- | --- |
| `results/scrna_scenic_cldn4/regulons.tsv` | **primary** CLDN4-high regulon table |
| `regulons_focus.tsv` | focus-TF slice |
| `n_table.tsv` | honest n |
| `tf_screen_cldn4_high.tsv` | 861-TF screen |
| `paired_high_vs_low.tsv` | sample-paired prior AUCell / TF RNA (n=5) |
| `sample_means_high_vs_low.tsv` | per-sample high/low means |
| `regulons_LUAD.tsv` | LUAD-only sensitivity |
| `fig_high_specific_counts.png` | high-specific counts |
| `fig_prior_aucell_paired.png` | underpowered paired AUCell |

---

## 中文

**结论：只交 CLDN4-high 调控子表，不宣称新的 ELF3–CLDN4。** GSE207422 恶性样细胞 9,782；样本水平 n=**12**（≥20 个恶性样）。组内中位数划分后 CLDN4-high **5,599** 个细胞。7/12 例 CLDN4 中位数为 0，低组为空；配对 AUCell 只有 **5** 例，**不够权**。BD_immune07（LUSC）占 high 细胞 40%，细胞数不能当患者重复。

KLF4 / KLF5 / TP63 / ELF3 / GRHL2 的 Pearson 邻域最紧（mean r 0.36–0.55），是上皮/紧密连接/角蛋白共表达，**不是结合**。这些 top50 在 CLDN4-low 里同样相关（high-specific 0/50）。OVOL1 是唯一 top50 多为 high-specific 的焦点 TF（42/50）。公共先验没有 ELF3→CLDN4；ELF3 先验 ∩ Pearson = CLDN7。ELF3–CLDN4 体 RNA 当作已知，不记新发现。
