# Reviewer-prebuttal FAQ (four story gaps)

**Compile-only.** Numbers copied from cited PRs. No new matrices, no invented p-values, no private 8KL / KD / HIS folds in the answers.

Use this page when a reviewer (or co-author) asks why the paper does **not** claim X. Sibling gates: `FORBIDDEN_VS_ALLOWED.md` (#747), `FIGURE_PLAN.md` (#750), elevator pitch (#749).

---

## Q1. Why not say “TROP2 drives cold”?

**Short answer.** TROP2 / TACSTD2 marks a cold niche in some public limbs, but it does **not** carry a stable patient-level T/NK association on the locked concordant-4, and public mediation does **not** put the cold under CLDN4. Say **observation / parallel limb**, not driver / cause / “through CLDN4.”

| Reviewer push | Reply with number | PR |
|---|---|---|
| “TACSTD2-high tumors are cold on your main scRNA panel.” | Concordant-4 malignant TACSTD2 %pos vs T/NK **ρ=−0.112**, p=0.46 (null). CLDN4 %pos on the same units is **ρ=−0.531**. | **#718**, **#733**, **#739** |
| “TROP2→immune acts through CLDN4.” | Mediation matrix **null** (0/12 primary rows support). After CLDN4, CLDN4–T/NK partial still **ρ=−0.543**; TACSTD2 direct coefficient **changes sign**. | **#718**, **#733**, **#747** |
| “CosMx proves TROP2 drives exclusion.” | CosMx TACSTD2 cold is real at **10/20 µm** (CD8+NK **0.455 / 0.670**; 8/8 & 5/5) but is a **parallel** limb. Inside CLDN4-low it still holds 8/8 — not “accounted for by CLDN4.” Do **not** swap locked CLDN4 **0.36 / 0.52** (50/100 µm) for TACSTD2 numbers. | **#726**, **#735**, **#738** |
| “Bulk RNA says TROP2 drives cold.” | OncoSG yes (ρ=**−0.309**); TCGA keratin-adj CD8 **weak** (−0.069); CPTAC protein–CD8A **no**; TISMO locks Tacstd2↑ **49/64** after ICB — **not** a T/NK-exclusion lock. | **#735**, **#739**, **#542** |

**Allowed sentence.** “TROP2-high niches co-occur with fewer short-range cytotoxic neighbors (CosMx) and lower CD8 scores in named bulks; the reproducible patient-level malignant→T/NK inverse on concordant-4 is **CLDN4 %pos**, stated as a parallel limb.”

**Forbidden.** “TROP2 drives / causes / mediates immune cold”; “TROP2-immune through CLDN4.”

---

## Q2. Why CLDN4, not the whole tight-junction module?

**Short answer.** After Tacstd2-high lifts junction / keratin programs, the **broad TJ score does not track fewer T/NK**. The immune-inverse that survives is **gene-level CLDN4** (and private KD pins CLDN4 among claudins). Whole-TJ is the fork that fails.

| Reviewer push | Reply with number | PR |
|---|---|---|
| “Just score the TJ module.” | Concordant-4 broad TJ (206 genes) vs T/NK **ρ=+0.137**, p=0.60 (null). TJ core **ρ=−0.191**, p=0.41. Same panel: CLDN4 %pos **ρ≈−0.52**. | **#741** |
| “Protein TJ should be cold.” | CPTAC-LSCC: CLDN4 protein vs ImmuneScore **ρ=−0.432** (n=78); structural **TJ-15 ρ=−0.082**, p=0.40 (n=108). CLDN4 **within** TJ, not mean TJ. | **#289**, **#737** |
| “Why not CLDN7 / EPCAM / whole claudin panel?” | Public surface ranking does **not** uniquely nail CLDN4 over CLDN7 (difference permutation p=0.17). What still pins CLDN4 is **private KD co-culture**, not the bulk surface table (LUAD CLDN4 was rank 7). | **#644**; PUBLIC_HANDOFF |
| “F11R / PARD3 as TJ load-bearing partners?” | F11R is a weak unreplicated TCGA call; PARD3 is not histology-independent. CLDN4 is the stronger TACSTD2 correlate in the same bulks. | **#167** |

**Allowed sentence.** “Tacstd2-high DEG enriches junction/keratin programs (#741 ORA apical junction #1); the member screen that associates with fewer T/NK is **CLDN4**, not the broad TJ mean.”

**Forbidden.** “TJ-high tumors exclude T/NK”; replace CLDN4 %pos with a TJ score.

---

## Q3. Why GSE137244 (cell-line bulk) as the mouse entry?

**Short answer.** It is the only **open** KL vs KP expression contrast that completely separates on both Tacstd2 and Cldn4 at a rank-test floor. It is **library bulk RNA-seq of cell lines**, used as the **genotype entry DEG**, not as immune exclusion. No public scRNA has both KL and KP with usable mouse *n* on both arms.

| Reviewer push | Reply with number | PR |
|---|---|---|
| “Why this accession?” | Deng et al. GSE137244: KL n=5 vs KP n=5, log2(FPKM+1). Tacstd2 **Δ=+3.238**, Cldn4 **Δ=+5.570**, exact MW **p=0.00794** (complete separation). TJ7 mean **+3.269**. | **#606**, **#685**, **#694** |
| “Why not in vivo GEMM bulk / scRNA?” | Fair in vivo series that separate on **Tacstd2** (GSE164758, GSE6135) do **not** separate on raw **Cldn4**. **No open scRNA** contains both KL and KP with usable mouse *n* on both arms. | **#612**, **#746** |
| “Cell-line bulk cannot prove immune cold.” | Agreed — and we do not ask it to. GSE137244 is the **KL>KP Tacstd2/Cldn4/TJ entry**. Immune-cold limbs are CosMx / concordant-4 / named bulks. Cell-line ISG ≠ infiltrate exclusion. | **#606**, PUBLIC_HANDOFF |
| “Isn’t Tacstd2-high on this matrix genotype-confounded?” | Yes: median-split Tacstd2 is fully KL vs KP. Catalog GSEA still does **not** put TJ/adhesion/Claudin in ranks #1–3 (`any_cohort_highlight_in_top3=false`). | **#744** |
| “TISMO LLC as KL?” | **No.** TISMO has no true KL/KP lung lines; LLC stays LLC. | **#542**, **#543** |

**Allowed sentence.** “Public KL>KP entry is GSE137244 cell-line RNA-seq (Tacstd2 +3.24, Cldn4 +5.57, MW p=0.00794); in vivo Tacstd2 separation exists without a second fair Cldn4 hit.”

**Forbidden.** “GSE137244 shows immune exclusion”; “TISMO LLC = KL”; merge private 8KL with public mouse Harmony.

---

## Q4. Why keep a private TJ#1 (and refuse a public one)?

**Short answer.** **TJ#1** = claiming Tacstd2-high GSEA puts tight-junction / adhesion / Claudin at rank **#1 (or #1–3)**. Public catalog GSEA does **not** support that headline. A private 8KL TJ#1 is allowed **only if labeled PRIVATE**; the public paper must not smuggle it in.

| Reviewer push | Reply with number | PR |
|---|---|---|
| “Your mouse GSEA says TJ is #1.” | Public honesty lock: **#744** `any_cohort_highlight_in_top3 = false` (GSE137244 / GSE164758 / GSE137396). Best catalog TJ ranks **~21 / not-positive / ~23**; Claudin family **123–296**. | **#744**, **#747** |
| “Human GSEA then?” | Concordant-4 TACSTD2 Q4 vs Q1: ORA **apical junction rank 1** and keratin GSEA **rank 2** are allowed (**#741**). KEGG Tight junction NES +2.04 is enriched but **NES-rank 22** — not “TJ #1.” | **#741**, **#747** |
| “Why mention private TJ#1 at all?” | Private 8KL can carry a genotype/DEG TJ-top claim the public catalog cannot. **#750** Box 3: *Private 8KL TJ#1 only if labeled PRIVATE.* Without that label it is the same lie class as public TJ-GSEA#1. | **#750**, **#749** |
| “Custom TJ core ranked #1 somewhere public.” | Custom epithelial-TJ panels can top a **pre-specified** gene-set list on GSE137244 (#746). That is **not** a catalog Hallmark/KEGG/GOBP #1–3 claim and must not be rewritten as one. | **#746** vs **#744** |

**Allowed sentence.** “Public Tacstd2-high catalog GSEA does not place TJ/adhesion/Claudin in the top three ranks; human ORA apical junction / keratin near-top is the public barrier enrichment we print. Any 8KL TJ-rank-#1 claim stays **PRIVATE**.”

**Forbidden.** “Tacstd2-high GSEA ranks TJ #1”; paste private 8KL NES into a public figure without a PRIVATE tag.

---

## One-page cheat sheet

| Question | Do say | Do not say | Lead PRs |
|---|---|---|---|
| TROP2 drives cold? | Parallel observation; CosMx 10/20 µm; OncoSG | Drives / mediates / through CLDN4 | **#718**, **#726**, **#735**, **#739** |
| CLDN4 not whole TJ? | CLDN4 %pos −0.531; CPTAC CLDN4 vs TJ-15 null | TJ-high excludes T/NK | **#741**, **#289**, **#712**, **#644** |
| Why GSE137244 bulk? | Only open KL/KP with Tacstd2+Cldn4 complete separation | Immune exclusion; LLC=KL | **#606**, **#612**, **#685**, **#746** |
| Why private TJ#1? | Public catalog TJ not top-3; PRIVATE label only | Public TJ-GSEA#1 | **#744**, **#747**, **#749**, **#750** |

---

## Provenance

See `provenance.json`. All four answers are story / honesty gates over already-locked public PRs. No private matrices were read for this page.
