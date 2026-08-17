# FINDING — pairwise merge GSE207422 + GSE148071 (CLDN4 only)

**Verdict: include 207422; patient-level T/NK is not a CLDN4-cold story; CellChat outgoing is mixed barrier-up + recruit-up.** Additive CLDN4-only merge of public **GSE207422** (Hu et al. 2023; 12 post-treatment) and **GSE148071** (Wu et al. 2021; 25 eligible of 42). Honest combo n = **37**, not 42+15. Patient-level malignant/epithelial CLDN4 vs T/NK fraction is weakly **positive** and not significant: **combo DL ρ = +0.145 (p = 0.416, I² = 0%, 95% CI −0.20 to +0.46)**. GSE207422 alone is ρ = +0.035 (n=12, p=0.914). The given author-DRMref / TISCH sensitivity is also null (combo ρ = +0.070, p=0.696).

CellChat-style kept splits (taken as given from PR #324 and PR #348) share **28** significant pairs with the same ΔP sign (**23 outgoing** Mal → T/NK). Consensus outgoing **higher from CLDN4-high**: **NECTIN2–TIGIT**, **MDK**, **HLA-E/F–CD8**, **JAM1/F11R–ITGAL/ITGB2**, **ICAM1**, **CXCL16–CXCR6**. Consensus outgoing **higher from CLDN4-low**: **SPP1–CD44**, laminins–CD44, **MIF**. **CD274–PDCD1 is not a consensus pair** (undetected in 148071; higher from low in 207422). **LGALS9–PTPRC is 148071-only.** TACSTD2 is not a gate. No dual-high.

Cell-pooled Hill probabilities are not a patient mixed model. Combo rho is the patient-level test.

---

## English

### Honest n

| Item | n | Note |
| --- | ---: | --- |
| GSE148071 deposited | **42** | Wu et al., *Nat Commun* 2021, PMID 33953163; not the test n |
| GSE148071 CellChat-eligible (≥25 epi **and** ≥25 T/NK) | **25** | marker-argmax epithelium = putative malignant |
| GSE148071 excluded | **17** | mostly epithelium with almost no T/NK |
| GSE207422 deposited samples | **15** | Hu et al., *Genome Med* 2023, PMID 36869384; 3 pre + 12 post |
| GSE207422 used (post, all pass 25/25) | **12** | pre-treatment excluded in the given CellChat folder |
| **Combo patient n (primary)** | **37** | **25 + 12; do not write 57** |
| Consensus LR (both sig, same sign) | **28** | 23 outgoing, 5 incoming |
| Dual-high / TACSTD2 gate | **0** | not used |

Per-sample counts: `results/patient_scores.tsv`. Dominance in the given CellChat pools is unchanged: GSE148071 eligible T/NK is thin (P7 = 26% of 4,025); GSE207422 post epithelium is dominated by BD_immune07 (60.9% of 8,407). Those facts belong to the cell-pooled LR, not to the patient rho.

### Combo rho — patient-level malignant CLDN4 vs T/NK

Full table: [`results/combo_rho.tsv`](results/combo_rho.tsv). Extra figures: [`figures/fig_combo_rho_scatter.png`](figures/fig_combo_rho_scatter.png), [`figures/fig_combo_rho_forest.png`](figures/fig_combo_rho_forest.png).

**Primary (same lineage as the CellChat folders): mean epithelial CLDN4 vs T/NK ÷ all cells, eligible ≥25/25.**

| Cohort | n | ρ | p | 95% CI | Definition |
| --- | ---: | ---: | ---: | --- | --- |
| GSE148071 | 25 | **+0.189** | 0.365 | −0.22 to +0.54 | marker-argmax epi mean vs TNK/all |
| GSE207422 | 12 | **+0.035** | 0.914 | −0.55 to +0.60 | same; 12 post |
| **GSE148071+GSE207422** | **37** | **+0.145** | **0.416** | **−0.20 to +0.46** | DL RE on Fisher-z; I² = 0% |

Stouffer (signed, weight √(n−3)): z = +0.82, p = 0.411. Fixed-effect ρ = +0.145 (τ² = 0).

Sensitivity (T/NK among epi+TNK only): combo ρ = +0.025, p = 0.887, n = 37.

**Given extracts (not re-audited):** TISCH locked 148071 n=25 ρ=+0.135 p=0.519; author DRMref 207422 n=12 ρ=−0.091 p=0.779; combo DL ρ = **+0.070** p=0.696 I²=0%. Same honest n=37. Mixed malignant definition — do not prefer this over the CellChat-matched primary.

**What this is not:** a significant negative CLDN4–T/NK correlation. Adding 207422 does not flip 148071 from weakly positive to cold. Do not write n=42 or n=15 as the rho n.

### Ligand table — outgoing CLDN4-high → T/NK (consensus)

Given kept splits: GSE148071 **median** (9,916 / 9,916 / 4,025 cells; 25 patients); GSE207422 **median_post** (4,203 / 4,204 / 33,760; 12 patients). CellChat R was not run. A pair is consensus if it is significant in **both** kept ligand tables and ΔP has the same sign.

Full consensus table: [`results/lr_table.tsv`](results/lr_table.tsv). Outgoing-only: [`results/lr_table_outgoing.tsv`](results/lr_table_outgoing.tsv). Extra figures: [`figures/fig_extra_ligand_table.png`](figures/fig_extra_ligand_table.png), [`figures/fig_extra_consensus_outgoing.png`](figures/fig_extra_consensus_outgoing.png), [`figures/fig_extra_key_pairs.png`](figures/fig_extra_key_pairs.png).

**Higher from CLDN4-high in both (outgoing):**

| Pair | Pathway | Class | ΔP 148071 | ΔP 207422 | mean ΔP |
| --- | --- | --- | ---: | ---: | ---: |
| NECTIN2–TIGIT | NECTIN | barrier\|inhibitory | +0.020 | +0.096 | +0.058 |
| MDK–NCL | MK | other | +0.064 | +0.023 | +0.043 |
| MDK–ITGA4/ITGB1 | MK | other | +0.058 | +0.027 | +0.042 |
| ICAM1–ITGAL/ITGB2 | ICAM | other | +0.009 | +0.054 | +0.032 |
| CD55–ADGRE5 | ADGRE | other | +0.019 | +0.042 | +0.031 |
| HLA-E–CD8A | MHC-I | inhibitory | +0.031 | +0.008 | +0.019 |
| JAM1/F11R–ITGAL/ITGB2 | JAM | barrier | +0.012 | +0.011 | +0.012 |
| **CXCL16–CXCR6** | CXCL | recruit | +0.002 | +0.019 | +0.011 |
| HLA-E–CD8B | MHC-I | inhibitory | +0.006 | +0.013 | +0.010 |

**Higher from CLDN4-low in both (outgoing):**

| Pair | Pathway | ΔP 148071 | ΔP 207422 | mean ΔP |
| --- | --- | ---: | ---: | ---: |
| SPP1–CD44 | SPP1 | −0.141 | −0.036 | −0.089 |
| SPP1–ITGA4/ITGB1 | SPP1 | −0.109 | −0.026 | −0.067 |
| LAMC1–CD44 | LAMININ | −0.022 | −0.085 | −0.053 |
| LAMB3–CD44 | LAMININ | −0.074 | −0.025 | −0.049 |
| LAMA3–CD44 | LAMININ | −0.026 | −0.060 | −0.043 |
| MIF–CD74/CD44 | MIF | −0.073 | −0.006 | −0.040 |
| MIF–CD74/CXCR4 | MIF | −0.064 | −0.009 | −0.036 |

**Incoming consensus (T/NK → Mal):** CD8A–CEACAM5 (+0.092 / +0.252), GZMA–PARD3, SEMA4D–PLXNB2, CD6–ALCAM higher into high; TNF–TNFRSF1A slightly lower into high in both.

### What is not supported

- **n = 42, n = 15, or n = 57** as the merge n. Eligible combo n = 37.
- A patient-level “CLDN4-high = T/NK-low” claim. Primary combo ρ is **positive** and non-significant.
- **CD274–PDCD1** as a shared outgoing hit. Not detected in 148071; higher from CLDN4-low in 207422.
- **LGALS9–PTPRC** as a pair-merge hit. Significant in 148071 only.
- Dual-high / TACSTD2-defined groups. This folder never uses TACSTD2 as a gate.
- CellChat R networks, LIANA, ICI response on 148071 (GEO has age/sex only), or a patient mixed model of ligand probability.
- Author CopyKAT IDs. Epithelium in both given CellChat folders is marker-argmax (putative malignant). 207422 author DRMref is sensitivity only.

### Files

`results/combo_rho.tsv`, `results/lr_table.tsv`, `results/lr_table_outgoing.tsv`, `results/honest_n.tsv`, `results/patient_scores.tsv`, `results/summary.json`. Extra figures under `figures/` and copied into `results/`.

---

## 中文

**结论：必须纳入 207422；患者水平 T/NK 不是 CLDN4 冷肿瘤故事；外向配体是屏障升高 + 招募也升高的混合结果。** GSE148071 合格 25/42，GSE207422 术后 12/15，合并诚实 n=**37**，不是 57。患者水平恶性/上皮 CLDN4 对 T/NK 比例的 combo DL ρ=**+0.145**（p=0.416，I²=0%）。207422 单独 ρ=+0.035（n=12）。给定 TISCH/DRMref 敏感性同样不显著（ρ=+0.070）。

两边都显著且同号的配体–受体 **28** 对（外向 23）。CLDN4 高上皮对 T/NK 共同升高：**NECTIN2–TIGIT、MDK、HLA-E、JAM1/F11R、CXCL16–CXCR6**。共同降低：SPP1、层粘连蛋白、MIF。**CD274–PDCD1 不是共识对。** LGALS9–PTPRC 只在 148071。不用 TACSTD2，无 dual-high。细胞池化概率不是患者混合模型。
