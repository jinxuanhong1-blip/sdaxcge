# Hunt: TACSTD2-high tumors vs tight junction / keratin-barrier / EMT

**Honest discovery report.** Pre-specified hypothesis. Every tested cohort is shown, including nulls and contradictions. Do not cite the supportive subset alone.

## Hypothesis (locked)

In public bulk and scRNA tumor data, **TACSTD2-high** samples GSEA-enrich **tight junction** and **keratinization / skin-barrier** programs and **deplete EMT**, with the leading-edge / co-expression intersection near **CLDN1, CLDN4, CLDN7, F11R, PARD3**.

## Short answer

**TJ / keratin-barrier up is common. EMT-down is not. The five-gene intersection is real in bulk for CLDN1/4 (±CLDN7, F11R), not for PARD3 as a stable member.**

| Claim | Holds? | Where |
|---|---|---|
| TJ GSEA up in TACSTD2-high | **Mostly yes** in bulk | TCGA most cancers; GSE31210; GSE207422 bulk (small n) |
| Keratinization / skin-barrier GSEA up | **Mostly yes** in bulk | Same; strongest in squamous / barrier histologies |
| EMT GSEA down | **No as a general rule** | Yes: LUSC, HNSC, BLCA, STAD, GSE31210. **Opposite** (Hallmark EMT up): TCGA-LUAD, BRCA, COAD, PAAD, SKCM, KIRC, GBM, … |
| Intersection CLDN1/4/7 F11R PARD3 | **Partial** | CLDN1 and CLDN4 are the stable pair. F11R and CLDN7 usually yes in bulk. **PARD3 is the weak link.** |
| GSE131907 malignant (sample pseudobulk) | **Null** | n=21, quartile 6 vs 6; no GSEA FDR<0.25; no focal q<0.05 |
| GSE207422 bulk | **Underpowered** | n=24; NES looks supportive; Hallmark EMT is *not* down; focal rhos are very high (treat as exploratory) |

LUAD — the disease of GSE207422 / GSE131907 — does **not** give a clean triad in TCGA RNA-seq: TJ/keratin up, Hallmark EMT **up**, CLDN7 and PARD3 **not** correlated with TACSTD2.

## What this is not

- Not evidence that TACSTD2 *causes* tight junctions. TACSTD2 (Trop-2) is itself an epithelial surface protein; GDLD literature already links it to claudins / F11R / ZO-1 in cornea.
- EMT-down in TACSTD2-high tumors is **partly tautological** when it happens: TACSTD2 tracks epithelial identity. Hallmark EMT is also a messy stromal/ECM set, so “EMT up” in bulk can mean stroma, not a transcription-factor EMT program.
- Within-cancer quartile splits remove LUAD-vs-LUSC confounding. They do **not** remove intra-cancer differentiation (basal vs luminal, squamous vs adeno components).
- GSEA FDR 0.25 is the classical exploratory cutoff. FDR 0.05 is also reported. This is a hunt.

## Methods (locked before looking at NES)

- **Split (primary):** TACSTD2 top vs bottom quartile. **Sensitivity:** median split (same direction unless noted).
- **Ranking:** Welch t, high minus low, on the supplied log-like matrix.
- **GSEA:** gseapy prerank, 1000 permutations, min_size=5, seed=42.
- **Primary sets (frozen in `data/genesets/primary_sets.json`):**
  - TJ: `KEGG_TIGHT_JUNCTION`, `GOBP_TIGHT_JUNCTION_ORGANIZATION`, `HALLMARK_APICAL_JUNCTION`, `CUSTOM_TJ_CORE`, `CUSTOM_CLAUDIN_PAR_FOCAL`
  - Keratin / barrier: `GOBP_KERATINIZATION`, `GOBP_CORNIFICATION`, `GOBP_ESTABLISHMENT_OF_SKIN_BARRIER`, `GOBP_KERATINOCYTE_DIFFERENTIATION`, `GOBP_EPIDERMAL_CELL_DIFFERENTIATION`, `CUSTOM_BARRIER_KERATIN`
  - EMT: `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`, `GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION`, `CUSTOM_EMT_CORE`
- TACSTD2 is **not** a member of any primary set (checked). NES is not self-inflated by the ranking gene.
- **Focal genes:** Spearman vs TACSTD2, BH-FDR within cohort (also context genes CDH1, VIM, ZEB1, …).
- **Complementary:** z-mean signature scores vs TACSTD2 (more honest than GSEA when n is small).
- **scRNA:** sample-level malignant/epithelial **pseudobulk** (log2 CPM+1), not mixed-cell bulk.
  - GSE131907: author `Cell_subtype == Malignant cells`, tumor/met sites only (not nLung/nLN); samples with <30 malignant cells dropped → **21 samples**.
  - GSE207422 scRNA: **no public cell-type table**. Epithelial gate: `(EPCAM+KRT8+KRT18+KRT19) ≥ 1 UMI` and `PTPRC=PECAM1=COL1A1/COL3A1=0` → 9,021 cells / 14 samples. Marker-based, not author annotation.
- **Verdict:** `supportive` = TJ-up AND keratin/barrier-up AND EMT-down (any primary set in that bucket, GSEA FDR cut) AND ≥3/5 focal genes rho>0 q<0.05. `n_high` or `n_low` < 8 → **underpowered** regardless of NES.
- Sources: UCSC Xena TCGA HiSeqV2 (primary tumors, sample type 01; SKCM 01+06); GEO GSE207422 (Hu et al. *Genome Med* 2023); GEO GSE131907 (Kim et al. *Nat Commun* 2020); GEO GSE31210 (Okayama stage I–II LUAD, GPL570). Gene sets: MSigDB v2023.2.Hs + KEGG v7.5.1 + custom modules.

## Verdicts (quartile split)

| Cohort | n | Q high/low | FDR&lt;0.25 | FDR&lt;0.05 | TJ↑ | K/barrier↑ | EMT↓ | focal+/5 |
|---|---:|---|---|---|:---:|:---:|:---:|---:|
| TCGA-LUAD | 515 | 129/129 | contradicts_EMT | contradicts_EMT | Y | Y | n | 3 |
| TCGA-LUSC | 502 | 126/126 | supportive | mixed | Y | Y | Y | 5 |
| TCGA-HNSC | 520 | 130/130 | supportive | supportive | Y | Y | Y | 5 |
| TCGA-BLCA | 407 | 102/102 | supportive | supportive | Y | Y | Y | 5 |
| TCGA-CESC | 303 | 76/76 | mixed | mixed | Y | Y | n | 5 |
| TCGA-BRCA | 1097 | 275/275 | contradicts_EMT | contradicts_EMT | Y | Y | n | 5 |
| TCGA-COAD | 286 | 72/72 | contradicts_EMT | contradicts_EMT | Y | Y | n | 0 |
| TCGA-ESCA | 184 | 46/46 | mixed | mixed | Y | Y | n | 3 |
| TCGA-STAD | 415 | 104/104 | supportive | supportive | Y | Y | Y | 5 |
| TCGA-PAAD | 178 | 45/45 | contradicts_EMT | contradicts_EMT | Y | Y | n | 5 |
| TCGA-OV | 304 | 76/76 | contradicts_EMT | contradicts_EMT | Y | Y | n | 4 |
| TCGA-UCEC | 176 | 44/44 | contradicts_EMT | contradicts_EMT | Y | Y | n | 3 |
| TCGA-SKCM | 472 | 118/118 | contradicts_EMT | contradicts_EMT | Y | Y | n | 4 |
| TCGA-KIRC | 533 | 134/134 | contradicts_EMT | contradicts_EMT | Y | Y | n | 3 |
| TCGA-GBM | 154 | 39/39 | contradicts_EMT | contradicts_EMT | Y | Y | n | 3 |
| GSE207422_bulk | 24 | 6/6 | underpowered | underpowered | Y | Y | Y* | 5 |
| GSE31210 | 226 | 57/57 | supportive | supportive | Y | Y | Y | 5 |
| GSE131907_malignant_pseudobulk | 21 | 6/6 | underpowered | underpowered | n | n | n | 0 |
| GSE207422_sc_epithelial_pseudobulk | 14 | 4/4 | underpowered | underpowered | Y | Y | n | 0 |

\*GSE207422 bulk: `CUSTOM_EMT_CORE` NES −1.63 FDR 0.023 (down), but **Hallmark EMT NES +1.16 FDR 0.239** (not down). The machine label “EMT↓” is the custom core only.

Adequately sized **supportive** (FDR&lt;0.25): **TCGA-LUSC, TCGA-HNSC, TCGA-BLCA, TCGA-STAD, GSE31210**.

## Key NES (quartile; selected sets)

Positive NES = enriched in TACSTD2-high. Full table: `gsea_key_sets_quartile.tsv`.

| Cohort | KEGG TJ | CUSTOM TJ | CLAUDIN-PAR | Keratinization | Skin barrier | Hallmark EMT | Custom EMT |
|---|---:|---:|---:|---:|---:|---:|---:|
| TCGA-LUAD | 1.55** | 2.01*** | 1.74** | 2.77*** | 2.58*** | **+1.55\*** | −0.77 ns |
| TCGA-LUSC | 1.31* | 1.48* | 1.87** | 3.22*** | 2.90*** | −1.36* | −1.36* |
| TCGA-HNSC | 1.61** | 2.47*** | 1.93** | 2.70*** | 2.59*** | −3.05*** | −2.53*** |
| TCGA-BLCA | 1.73** | 2.34*** | 1.74** | 2.05** | 2.00** | −2.78*** | −2.23*** |
| TCGA-STAD | 1.92** | 2.41*** | 1.84** | 2.95*** | 2.69*** | −1.20* | −1.66** |
| GSE31210 | 1.63* | 2.45*** | 1.83** | 2.00*** | 2.17*** | −2.16*** | −1.85** |
| GSE207422 bulk | 1.96*** | 2.60*** | 1.91*** | 2.46*** | 2.47*** | +1.16 ns | −1.63* |
| GSE131907 mal. pb | 1.24 ns | 1.53 ns | 1.27 ns | 0.99 ns | 1.27 ns | +0.94 ns | +0.82 ns |
| GSE207422 sc epi pb | 1.63* | −0.93 ns | 1.75** | 2.61*** | 2.50*** | **+1.73\*** | **+1.69\*** |

Stars: \* FDR&lt;0.05, \*\* FDR&lt;0.01, \*\*\* FDR&lt;0.001, ns FDR≥0.05. Exact FDR in the TSV.

## Focal genes vs TACSTD2 (Spearman, quartile table = all samples)

| Cohort | CLDN1 | CLDN4 | CLDN7 | F11R | PARD3 |
|---|---:|---:|---:|---:|---:|
| TCGA-LUAD | 0.37* | 0.46* | 0.05 | 0.18* | −0.02 |
| TCGA-LUSC | 0.42* | 0.41* | 0.40* | 0.25* | 0.15* |
| TCGA-HNSC | 0.28* | 0.52* | 0.50* | 0.35* | 0.35* |
| TCGA-BLCA | 0.19* | 0.53* | 0.36* | 0.28* | 0.13* |
| GSE31210 | 0.38* | 0.47* | 0.32* | 0.41* | 0.20* |
| GSE207422 bulk | 0.88* | 0.86* | 0.83* | 0.87* | 0.66* |
| GSE131907 mal. pb | 0.25 | 0.24 | 0.01 | 0.14 | 0.22 |
| GSE207422 sc epi pb | 0.52 | 0.49 | 0.41 | 0.22 | 0.38 |
| TCGA-COAD | −0.07 | −0.06 | 0.04 | 0.05 | 0.02 |

\* BH-FDR &lt; 0.05. GSE207422 bulk rhos are large because n=24; they are not a confirmation cohort.

## Leading-edge intersection (TJ sets with NES&gt;0 and FDR&lt;0.25)

Gene in the leading edge of **any** such TJ set:

| Gene | # cohorts (of 18 TJ-up) |
|---|---:|
| CLDN1 | 17 |
| CLDN7 | 16 |
| CLDN4 | 15 |
| F11R | 15 |
| PARD3 | 11 |

GSE131907 is absent: no TJ set reached FDR&lt;0.25. TCGA-COAD is TJ-up on some sets but **none** of the five genes are in those leading edges.

**Stable intersection: CLDN1–CLDN4 (usually +CLDN7 +F11R). PARD3 is optional, not core.**

## scRNA cell-level supplement (powered, but clustered)

Sample-pseudobulk GSEA is the fair comparison to bulk and was **null / underpowered**. Among individual cells, Spearman of log1p(UMI) vs TACSTD2 (no library-size correction — **depth-confounded**; do not over-read EMT):

| | CLDN1 | CLDN4 | CLDN7 | F11R | PARD3 | CUSTOM_CLAUDIN_PAR | CUSTOM_EMT | Hallmark EMT |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| GSE131907 malignant cells (n=24,784) | 0.33 | 0.63 | 0.48 | 0.42 | 0.29 | 0.63 | **+0.42** | **+0.61** |
| GSE207422 epithelial cells (n=9,021) | 0.61 | 0.70 | 0.62 | 0.54 | 0.47 | 0.75 | +0.10 | **+0.56** |

Directionally, TACSTD2 tracks the claudin–F11R–PARD3 module **inside** malignant/epithelial cells (CLDN4 strongest). EMT is **not** down at cell level in either scRNA cohort. Because these correlations are on raw UMI, part of the positive EMT signal is likely detection/depth. The library-normalized sample-pseudobulk remains the honest GSEA test and does not support the full hypothesis in GSE131907.

## Honest take

1. **The hunt finds a TJ / barrier program next to TACSTD2, not a universal EMT-low state.** Squamous and barrier epithelia (LUSC, HNSC, BLCA) plus STAD and one LUAD array (GSE31210) are the only adequately sized cohorts that match the *full* pre-specified triad.
2. **TCGA-LUAD disagrees with GSE31210 on EMT.** Same histology, different platform / population / stage mix. Do not collapse them into “LUAD supports EMT-down.”
3. **GSE131907 (the named LUAD atlas) does not support sample-level GSEA of the triad.** That is a real negative, not a missing analysis.
4. **GSE207422 bulk is too small for GSEA.** The huge focal rhos are compatible with the claudin module but are not independent confirmation.
5. **PARD3 should not be billed as part of a required intersection.** CLDN1/CLDN4 are the genes that actually recur.
6. Circularity: TACSTD2 is epithelial. Keratin/barrier enrichment in LUSC/HNSC/CESC can be **differentiation**, not a TACSTD2-specific mechanism. The more informative tests are LUAD and the malignant-cell atlas — and those are the ones that fail the EMT-down arm.

## Files

- `verdicts.tsv` — machine-readable labels
- `gsea_prerank_quartile.tsv` / `gsea_prerank_all.tsv` — all primary sets
- `gsea_key_sets_quartile.tsv` / `gsea_key_sets_wide.tsv` — the sets that matter
- `focal_gene_correlations.tsv`
- `signature_scores_vs_tacstd2.tsv`
- `leading_edge_focal.tsv`
- `celllevel_focal.tsv` / `celllevel_signatures.tsv`
- `cohort_inventory.tsv`
- `figures/nes_heatmap.png`, `figures/focal_rho_heatmap.png`
- `analysis.log`

Reproduce: `python3 scripts/hunt_tj_gsea.py` then `python3 scripts/hunt_tj_gsea_celllevel.py` (requires the GEO/Xena files under `data/raw/`, not committed).
