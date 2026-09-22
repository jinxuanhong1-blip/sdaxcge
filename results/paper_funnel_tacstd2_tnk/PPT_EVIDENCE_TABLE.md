# PAPER FUNNEL — Tacstd2-high → fewer T/NK (public corroboration)

**One slide table.** Numbers are copied from existing public-analysis PRs only. Nothing was re-fit or invented for this file.

**Headline (honest):** Spatial CosMx and OncoSG support Tacstd2-high → fewer T/NK-class neighbors or bulk CD8/NK scores. TCGA keratin-adjusted CD8 is weak but same sign. CPTAC TROP2 protein exists: LUAD xCell-CD8 agrees; CD8A protein does not. TISMO locks Tacstd2 up after ICB — it does **not** support fewer T/NK.

---

## PPT paste table

| # | Dataset | Endpoint | Number to show | Supports fewer T/NK? | One-line caveat |
|---|---|---|---|---|---|
| 1 | **CosMx He2022** NSCLC | TACSTD2-hi vs lo: CD8+NK neighbors | **10 µm ratio 0.455; 20 µm 0.670; 8/8 & 5/5; sign P=0.031** | **YES** | Short-range; 50/100 µm not 8/8 |
| 2 | CosMx He2022 (same) | Immune-neighbor fraction | 10 µm **0.519**; 20 µm **0.643**; 8/8 & 5/5 | YES (immune) | Companion to row 1 |
| 3 | **TCGA** 8 cohorts | TACSTD2 vs CD8 after **KRT8/18/19** | pooled **ρ=−0.069**; 6/8 neg; p=0.043; I²=76% | **WEAK** | CESC +; PAAD flips after keratin |
| 4 | TCGA (sensitivity) | same vs CD3/CD8 after **KRT5/6/14+TP63** | pooled **ρ=−0.119**; 7/8; p=5e-5 | post-hoc stronger | Sweep max (#705); prefer row 3 |
| 5 | **OncoSG** LUAD | TACSTD2 vs CD8A ‖ purity | **ρ=−0.309**; n=169; p=4.7e-5 | **YES** | Surgical EA LUAD; not ICI |
| 5b | OncoSG (same) | IMSIG T / NK ‖ purity | T **−0.404**; NK **−0.421** | YES | Not the CD8A primary |
| 6 | **CPTAC** LUAD protein | TROP2 vs xCell CD8 | **ρ=−0.289**; WES partial **−0.261** | **PARTIAL** | LSCC ns (−0.08); treatment-naive |
| 7 | CPTAC protein–protein | TROP2 vs **CD8A protein** | LUAD −0.038; LSCC −0.103; **CIs cross 0** | **NO** | Protein present; exclusion not in CD8A protein |
| 8 | **TISMO** ICB | Tacstd2 treated − baseline | **49/64 up**; Wilcoxon **p=5.8e-5** | **NO (different claim)** | ICB induction only; baseline CD8 link is **positive** |

---

## Per-cohort TCGA keratin-adj TACSTD2–CD8 (row 3; partial ρ)

| Cohort | n | partial ρ | sign |
|---|---:|---:|:---:|
| LUAD | 516 | −0.104 | − |
| LUSC | 501 | −0.221 | − |
| BRCA | 1095 | −0.039 | − |
| CESC | 304 | **+0.118** | **+** |
| KIRC | 533 | −0.047 | − |
| STAD | 412 | −0.079 | − |
| BLCA | 406 | −0.158 | − |
| PAAD | 178 | **+0.029** | **+** (marginal was −0.139) |

6/8 negative after KRT8/18/19. Source: PR #593 `immune_exclusion.tsv`.

---

## Verdict line for the slide footer

> CosMx 10–20 µm and OncoSG CD8/NK support Tacstd2-high → fewer T/NK. TCGA keratin-adj CD8 is small (ρ≈−0.07). CPTAC TROP2 protein: LUAD xCell-CD8 yes, CD8A protein no. TISMO = Tacstd2↑ after ICB (49/64), not T/NK exclusion.

---

## Do not say

- Do not write Visium co-localization as spatial exclusion.
- Do not merge private 8KL with public mouse Harmony.
- Do not call TISMO LLC a KL/KP lung model.
- Do not quote CPTAC TROP2 vs CD8A protein as corroborating exclusion.
- Do not replace locked CosMx CLDN4 cytotoxic ratios 0.36 / 0.52 with TACSTD2 numbers.
