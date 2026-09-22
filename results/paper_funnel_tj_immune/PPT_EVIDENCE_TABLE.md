# PAPER FUNNEL — TJ / CLDN4 / junction → fewer T/NK / CD8 / IFN (public corroboration)

**One slide table.** Numbers copied from existing public-analysis PRs only. Nothing was re-fit or invented for this file.

**Funnel position:** After Tacstd2-high → TJ / keratin-barrier enrichment (A8; LUAD KEGG TJ NES **2.09** FDR 0.005; LUSC TJ assembly NES **2.18** FDR 0.004; PR #200 / #46), this slide asks whether the **TJ / CLDN4 / junction module** itself is inversely correlated with T/NK, CD8, and IFN at the patient / mouse unit.

**Headline (honest):** Human lung scRNA **strongly** supports CLDN4 / barrier-high → fewer T/NK and lower malignant IFN/MHC pathway scores (concordant-4 ρ=**−0.531**, I²=0%, N=65; fgsea IFNγ NES=**−3.85**). CosMx CLDN4 exclusion is locked (cytotoxic ratio **0.36 / 0.52**). Public mouse GEMM forest does **not** lock the same inverse (pooled ρ=**+0.48** vs T/NK). Small-n mouse max-|ρ| can hit −1 on 4 mice — exploratory ceiling, not confirmatory.

---

## PPT paste table

| # | Dataset | Unit | Exposure → endpoint | Number to show | Supports inverse? | One-line caveat |
|---|---|---|---|---|---|---|
| 1 | **Concordant-4** scRNA (GSE123902+131907+205335+189357) | patient / donor / sample | malignant CLDN4 %pos → T/NK fraction | **ρ=−0.531**; p=1.65×10⁻⁵; I²=**0%**; N=**65**; Q4/Q1 Cliff δ=**−0.724** (19/16) | **YES — STRONG** | Locked panel; joint max of 6,435-panel sweep (#712) |
| 2 | Concordant-4 (same) | same | largest \|ρ\| with I²=0 | **ρ=−0.533** (partial KRT18); Cliff δ=−0.697 | YES (I²=0 max) | Third-decimal gain; |Cliff| smaller than locked |
| 3 | Concordant-4 (same) | same | largest \|ρ\| sign-only (I² ignored) | **ρ=−0.552**; I²=67%; N=60 | YES but heterogeneous | GSE123902 ρ=−0.107; GSE189357 ρ=−0.940 — do not promote |
| 4 | **GSE131907** alone | patient (n=21) | malignant CLDN4 %pos → T/NK | **ρ=−0.478**; p=0.028 (mean CLDN4 ρ=−0.532) | **YES** | Sample-level comparison ρ=−0.522; IFN family near-null here |
| 5 | Concordant-4 malignant PB | patient Q4 vs Q1 (18/16) | pathway NES (fgsea; CLDN4 held out) | IFNγ **NES=−3.85**; IFNα **−3.42**; MHC-I **−2.75**; EMT **−3.42** | **YES — IFN/MHC DOWN** | Bootstrap CIs stay negative. KEGG TJ NES +1.13 **crosses zero** |
| 6 | Concordant-4 malignant PB | patient continuous (N=64) | GSVA Spearman vs CLDN4 %pos | IFNγ **ρ=−0.40**; IFNα **ρ=−0.32**; MHC **ρ=−0.25** (DL meta) | YES (pathway scores) | GSE131907 alone flat for IFN |
| 7 | Concordant-4 malignant family DE | patient Q4 vs Q1 (18/16) | median logFC | IFN **−0.51**; MHC-I/APM **−0.69**; TJ (CLDN4 out) **+0.11** | YES (IFN/MHC); TJ mild UP | Matches thesis extras; not KEGG NES |
| 8 | **CosMx He2022** NSCLC | section / patient | CLDN4-hi vs lo cytotoxic neighbors | **50 µm ratio 0.36; 100 µm 0.52**; 8/8 & 5/5; sign P=0.031 | **YES — spatial** | Exclusion, not muzzling (effector genes hi/lo 1.11–1.22) |
| 9 | TCGA-LUAD bulk (orthogonal) | tumor | structural 15-gene TJ → CD8 | **ρ=−0.29**; p=1.6×10⁻¹¹; purity-partial **−0.26** | YES (bulk TJ↔CD8) | KEGG TJ does **not** support (ρ≈+0.07). Not scRNA |
| 10 | Public mouse **GEMM forest** | mouse | epithelial Cldn4 %pos → T/NK | pooled **ρ=+0.477**; p=0.30; I²=72%; k=5; 34 mice | **NO** | Opposite sign from human; interval covers 0 |
| 11 | Public mouse GEMM forest | mouse | Cldn4 %pos → epithelial IFN | pooled **ρ=+0.092**; p=0.80; 56 mice | **NO** | Near null |
| 12 | Public mouse max-\|ρ\| (GSE154977/165641/180963) | mouse | Cldn4% → T fraction | **ρ=−1.000**; n=4; exact p=0.083 | EXPLORATORY only | Ceiling on 4 mice; not confirmatory. IFN/APM max is **+1.0** |

---

## Concordant-4 cohort ρ (row 1; locked %pos)

| Cohort | Unit | n | Spearman ρ | p |
|---|---|---:|---:|---:|
| GSE123902 | donor | 13 | −0.659 | 0.014 |
| GSE131907 | sample | 21 | −0.522 | 0.015 |
| GSE205335 | patient | 22 | −0.435 | 0.043 |
| GSE189357 | patient | 9 | −0.600 | 0.088 |
| **DL pool** | | **65** | **−0.531** | **1.65×10⁻⁵** (I²=0%) |

Sources: PR #503 / #539 / #712.

---

## Pathway-level scores to put next to the ρ (rows 5–7)

| Score | Statistic | Value | Source |
|---|---|---|---|
| Hallmark IFNγ fgsea NES (Q4 vs Q1) | NES (boot 95%) | **−3.85** (−4.01 to −2.70) | #691 |
| Hallmark IFNα fgsea NES | NES (boot 95%) | **−3.42** (−3.74 to −1.99) | #691 |
| MHC-I/APM fgsea NES | NES (boot 95%) | **−2.75** (−3.05 to −1.76) | #691 |
| KEGG tight junction fgsea NES | NES (boot 95%) | +1.13 (−1.12 to +1.42) — **crosses 0** | #691 |
| GSVA IFNγ vs %pos (DL) | Spearman ρ | **−0.399** | #691 |
| GSVA IFNα vs %pos (DL) | Spearman ρ | **−0.319** | #691 |
| Family DE IFN / MHC / TJ | median logFC Q4−Q1 | **−0.51 / −0.69 / +0.11** | #503 |

---

## Largest \|ρ\| callouts (for the slide footnote)

| Scope | Rule | Largest \|ρ\| | Keep as primary? |
|---|---|---|---|
| Human concordant-4 | locked %pos (joint max \|ρ\|+‖Cliff‖) | **−0.531** | **YES — primary** |
| Human concordant-4 | max \|ρ\| among I²=0 panels | **−0.533** | sensitivity |
| Human concordant-4 | max \|ρ\| sign-only | **−0.552** | **NO** (I²=67%) |
| Human GSE131907 | patient %pos | **−0.478** | cohort row |
| Mouse max-\|ρ\| grid | T fraction, n≥4 | **−1.000** (n=4, p_exact=0.083) | **NO** as confirmation |
| Mouse GEMM REML forest | pre-specified primary | **+0.477** | honest null/opposite |

---

## Verdict line for the slide footer

> After Tacstd2-high TJ enrichment (A8), malignant CLDN4 / barrier tracks fewer T/NK in human lung scRNA (concordant-4 ρ=−0.53, I²=0%, N=65) and CosMx exclusion (0.36/0.52). Pathway scores: IFNγ NES −3.85, MHC −2.75. Public mouse GEMM does not lock the inverse (ρ=+0.48). Do not sell mouse max-|ρ|=−1 on n=4 as confirmation.

---

## Do not say

- Do not invent a public-mouse ρ≈−0.5 that matches human.
- Do not merge private 8KL with public mouse Harmony.
- Do not quote KEGG TJ NES as “up” when the bootstrap crosses zero (#691).
- Do not quote cell counts (22,653 / 208,506 / 765,771) as n — unit is patient / mouse / section.
- Do not replace locked CosMx CLDN4 0.36/0.52 with TACSTD2 short-range ratios.
- Do not claim CD8/NK **state** (effector/exhaustion) is down — fraction exclusion stays; state is null (#649).
