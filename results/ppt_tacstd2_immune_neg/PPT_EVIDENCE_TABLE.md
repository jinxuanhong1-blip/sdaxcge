# PPT OPENING (corrected) — TROP2/TACSTD2-high ↔ immune NEGATIVE

**Claim on this slide:** TACSTD2 / TROP2-high associates with **fewer T/NK neighbors or lower CD8 / ImmuneScore**.  
**Not on this slide:** drug resistance, ICB non-response, Tacstd2↑ after ICB (TISMO 49/64), Visium same-spot co-localization as exclusion.

Numbers are **copied from existing public PRs only**. Nothing was re-fit or invented.

---

## Six-cohort paste table

| # | Cohort | Endpoint | Honest number | Call |
|---|---|---|---|---|
| 1 | **CosMx He2022** NSCLC | TACSTD2-hi vs lo: CD8+NK neighbors | **10 µm ratio 0.455; 20 µm 0.670; 8/8 & 5/5; sign P=0.031** | **YES** (short-range) |
| 1b | CosMx (same) | Immune-neighbor fraction | 10 µm **0.519**; 20 µm **0.643**; 8/8 & 5/5 | YES (immune) |
| 2 | **concordant-4** scRNA n=65 | Malignant TACSTD2 %pos vs T/NK | **ρ=−0.112; p=0.46; I²=15%** | **NULL** |
| 2b | concordant-4 (pipeline check) | Malignant **CLDN4** %pos vs T/NK | ρ=**−0.531**; p=1.6×10⁻⁵; I²=0% | CLDN4 yes — **not TACSTD2** |
| 3 | **TCGA** LUAD+LUSC | TACSTD2 vs CD8A ‖ ABSOLUTE | pooled **ρ=−0.191** (n=995); LUAD −0.104; LUSC −0.244 | **YES** (CD8; LUSC stronger) |
| 3b | TCGA (same) | TACSTD2 vs ESTIMATE **ImmuneScore** ‖ ABSOLUTE | pooled **ρ=−0.107**; LUAD **+0.072** (ns); LUSC **−0.131** | **PARTIAL** (LUSC-only) |
| 3c | TCGA 8 cohorts | TACSTD2 vs CD8 ‖ **KRT8/18/19** | pooled **ρ=−0.069**; 6/8 neg; p=0.043; I²=76% | **WEAK** |
| 4 | **OncoSG** LUAD n=169 | TACSTD2 vs CD8A ‖ purity | **ρ=−0.309**; p=4.7×10⁻⁵ | **YES** |
| 4b | OncoSG (same) | Immune / IMSIG T / NK ‖ purity | immune **−0.318**; T **−0.404**; NK **−0.421** | YES |
| 5 | **CPTAC** LUAD protein | TROP2 vs xCell CD8 / immune | CD8 **ρ=−0.289** (WES **−0.261**); immune **−0.309** (WES **−0.272**) | **YES** (LUAD protein) |
| 5b | CPTAC LSCC protein | same | CD8 **−0.080** (p=0.41); immune **−0.133** (p=0.17) | **NULL** (LSCC) |
| 6 | **Mouse public** (5 studies, 34 mice) | Tacstd2 %epi vs T fraction | pooled **c=−0.277**; perm p=**0.137** | **WEAK / ns** |
| 6b | Mouse (largest study) | GSE295824 Tacstd2 % vs T | **ρ=−0.612**; n=16; p=0.012 | Study-level yes; not the pool |

---

## Slide footer (one line)

> CosMx 10–20 µm and OncoSG CD8/NK support TACSTD2-high → fewer T/NK / lower immune scores. TCGA CD8 is negative (ImmuneScore LUSC-only). CPTAC TROP2 protein: LUAD yes, LSCC no. Concordant-4 TACSTD2 vs T/NK is **null** (CLDN4 −0.53). Public mouse pool is **ns**. **Not drug resistance.**

---

## Do not say on this slide

- Do not write Tacstd2↑ after ICB (TISMO 49/64) as immune-negative evidence.
- Do not write Visium co-localization as spatial exclusion.
- Do not merge private 8KL with public mouse.
- Do not replace locked CosMx **CLDN4** cytotoxic ratios 0.36 / 0.52 with TACSTD2 numbers.
- Do not quote concordant-4 CLDN4 ρ=−0.53 as if it were TACSTD2.
- Do not quote CPTAC LSCC or concordant-4 TACSTD2 as corroborating exclusion.
