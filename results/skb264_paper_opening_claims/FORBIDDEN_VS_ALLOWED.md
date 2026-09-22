# SKB264 paper opening — forbidden vs allowed claims (one page)

**Use this page before drafting the opening.** Sources: PRs #712, #718, #726, #735, #741, #744. Public only. No private 8KL / SKB264 matrices.

Three lies this page exists to stop: **mediation**, **scRNA-merge**, **TJ-GSEA#1**.

---

## FORBIDDEN (do not write)

| Lie class | Forbidden sentence shape | Why locked out |
|---|---|---|
| **Mediation** | “TACSTD2/TROP2–immune cold is mediated by / depends on / acts through CLDN4.” | **#718 call = null** (0/12 primary rows support). Concordant-4 %pos ACME CI crosses 0; mean-log1p meta is **opposite**. |
| **Mediation** | “Short-range CosMx TACSTD2 exclusion shrinks once you account for CLDN4.” | **#726**: inside CLDN4-low, 10/20 µm immune-fraction still **8/8 & 5/5**; median rank attenuation only **0.124 / 0.103**. Not causal mediation. |
| **Mediation** | “CLDN4 is the TROP2 immune axis” / replace locked CLDN4 CosMx **0.36 / 0.52** with TACSTD2 numbers. | Axes are **parallel**, not nested. Locked CLDN4 ratios stay. TACSTD2 CosMx lock is **10/20 µm**, not 50/100. |
| **scRNA-merge** | Inflate concordant-4 with GSE148071 / GSE127465 / GSE154826 / GSE200563 / E-MTAB-13526 / GSE207422. | **#712 / #741**: locked set is only **GSE123902 + GSE131907 + GSE205335 + GSE189357** (n=65 for %pos). Extra sets dilute or flip. |
| **scRNA-merge** | Merge private 8KL with public mouse Harmony / call TISMO LLC a KL/KP lung line. | Public handoff + **#735/#744**: no true KL/KP in TISMO; private matrices are not in this repo. |
| **scRNA-merge** | Quote **n=65** as the TACSTD2 Q4-vs-Q1 DE *n*, or quote cell counts as *n*. | **#741**: DE contrast = **19 vs 15**; expression units = **64**. Unit = patient/donor/sample. |
| **TJ-GSEA#1** | “Tacstd2-high GSEA ranks tight junction / adhesion / Claudin **#1–3**.” | **#744**: `any_cohort_highlight_in_top3 = false` on GSE137244 / GSE164758 / GSE137396. Best TJ ranks **21 / not-positive / 23**; Claudin **123–296**. |
| **TJ-GSEA#1** | “KEGG Tight junction is the top NES hit” (human or mouse). | **#741**: human GSEA #1 = TNFA; KEGG TJ NES +2.04 but **NES-rank 22**. **#744**: mouse KEGG TJ ranks **31 / negative / 190**. |
| **TJ-GSEA#1** | “Broad TJ-high tumors exclude T/NK” / use TJ score in place of CLDN4 %pos. | **#741**: broad TJ vs T/NK **ρ=+0.137**, p=0.60; core **ρ=−0.191**, p=0.41. Immune exclusion stays on **CLDN4 %pos**. |

---

## ALLOWED (safe opening limbs)

| Limb | Allowed claim | Number to print | PR |
|---|---|---|---|
| **CLDN4 exclusion (primary)** | Malignant CLDN4 %pos (UMI>0) vs T/NK on locked concordant-4. | **ρ=−0.531**; p=1.65×10⁻⁵; I²=0%; N=**65**; Cliff δ=−0.724. Joint max of 6,435-panel sweep. | **#712** (locks #539) |
| **CosMx CLDN4 geography** | CLDN4-high tumor cells have fewer cytotoxic neighbors; exclusion, not muzzling. | Cytotoxic ratio **0.36 / 0.52** at 50/100 µm; 8/8 & 5/5; sign P=0.031. Effector genes hi/lo **1.11–1.22**, 0/8 down. | locked CosMx / handoff |
| **TACSTD2 CosMx (separate)** | TACSTD2-high tumor cells have fewer immune/CD8 neighbors at **short range**. | 10/20 µm CD8+NK ratios **0.455 / 0.670** (or immune-fraction **0.519 / 0.643**); 8/8 & 5/5; P=0.031. | **#726**, **#735** |
| **TACSTD2 bulk corroboration** | Same-sign fewer T/NK-class scores in named public bulks — with strength labels. | OncoSG CD8A ‖ purity **ρ=−0.309**; TCGA keratin-adj CD8 pooled **ρ=−0.069** (WEAK); CPTAC LUAD xCell-CD8 **PARTIAL**; TROP2–CD8A **protein NO**; TISMO = Tacstd2↑ **49/64**, **not** T/NK exclusion. | **#735** |
| **TACSTD2 → barrier DEG (human only)** | On concordant-4, TACSTD2 Q4 vs Q1 malignant DEG lifts keratin / adhesion / claudin **family scores**; ORA apical junction is top. | Family FDR: KERATIN 0.0055, adhesion 0.0042, CLAUDIN 0.0062, TJ 0.0042. ORA **HALLMARK_APICAL_JUNCTION rank 1**. GSEA keratin **rank 2** (not “TJ #1”). | **#741** |
| **Parallel, not nested** | State CLDN4 and TACSTD2 immune associations as **co-occurring public limbs**. | After CLDN4, CLDN4–T/NK partial still **ρ=−0.543**. TACSTD2–T/NK raw on %pos is **null** (ρ=−0.112). | **#718** |

---

## Opening recipe (copy-safe)

1. Lead with **CLDN4 %pos → fewer T/NK** (ρ=−0.531) and CosMx CLDN4 **0.36/0.52**.
2. Add **TACSTD2 short-range CosMx** (10–20 µm) and OncoSG as a **second** limb — never as “because CLDN4.”
3. If mentioning barrier programs after TACSTD2-high: say **human ORA apical junction / keratin near-top (#741)**; never mouse **TJ-GSEA #1–3 (#744)**.
4. Keep cohorts named. Never merge. Never mediate. Never swap TJ score for CLDN4 %pos.

**One-line veto:** If a draft sentence contains *mediated*, *through CLDN4*, *merged scRNA*, or *TJ ranks #1*, delete it.
