# Paper funnel: public mouse Tacstd2 (no private 8KL)

Public data only. Private 8-KL single-cell matrices were not read and were not merged with these series. Inferential unit is the mouse (or cell line / tumor as labeled). Tacstd2 is never used to call epithelium.

Three arms, all thesis-aligned where the numbers allow, and honest where they do not:

1. **KL > KP Tacstd2** (bulk / array / cell line)
2. **Tacstd2% vs T fraction** on the public integrate scRNA mice (GSE154977 + GSE180963 + GSE165641)
3. **Tacstd2-high epithelial DEG → tight junction** on those same mice

Detail files: `FINDING_KL_KP.md` (tables), `FINDING_TACSTD2_MAXRHO.md`, `FINDING_DEG_TJ.md`.

## 1. KL > KP Tacstd2

Locked GSE137244 cell-line reference reproduces: **Tacstd2 Δ = +3.238**, exact Mann–Whitney p = **0.00794** (5 vs 5; complete separation; the n=5 vs 5 floor). Cldn4 on the same matrix is +5.570 (locked). TJ7 mean is +3.269 (not substituted for the locked TJ +3.03).

Public series with n large enough for two-sided exact MW p < 0.05 also separate on Tacstd2 (raw Cldn4 does not):

| accession | setting | n KL vs KP | Tacstd2 Δ | MW p | role |
|---|---|---|---:|---:|---|
| GSE164758 | primary_bulk_tumor | 9 vs 8 | +0.858 | 8.227e-05 | new |
| GSE6135 | mouse_all_histology | 7 vs 5 | +1.795 | 0.002525 | new_tacstd2_same_cohort_as_prior_cldn4 |
| GSE137244 | cultured_cells | 5 vs 5 | +3.238 | 0.007937 | locked_reference |
| GSE6135 | mouse_adeno_only | 4 vs 5 | +0.989 | 0.01587 | sensitivity_adeno_only |
| GSE137396 | in_vivo_nodule | 5 vs 5 | +1.144 | 0.09524 | rescore |
| GSE244452 | syngeneic_bulk_tumor | 3 vs 3 | +7.228 | 0.1 | new |
| GSE274352 | cultured_cells_empty_vector | 3 vs 3 | -1.543 | 0.1 | new |
| GSE274351 | lcm_early_adenoma | 5 vs 5 | -0.860 | 0.5476 | rescore |

Headline public in vivo Tacstd2 KL>KP calls: **GSE164758** primary tumors (+0.858, p=8.2×10⁻⁵, 9 vs 8) and **GSE6135** mouse-level array (+1.795, p=0.00253, 7 vs 5; adenocarcinoma-only still +0.989, p=0.0159). Series that go the other way or sit on the n=3 floor are in `tables/kl_kp/focus_contrasts.tsv` and are not dropped from the table.

## 2. Tacstd2% vs T fraction (integrate mice)

Exposure = percent of locked-gate epithelial cells with Tacstd2 detection. GSE154977 is an AT2-lineage FACS sort and is a T-fraction design no-go. Eligible mixed-digest mice for T fraction: GSE165641 KL1/KL2 + GSE180963 K/KL (n=4).

**Prespecified locked gate, Tacstd2 count > 0, T/NK fraction:** ρ = **−1.000**, exact permutation p = **0.0833** (2/24), the floor for n=4. Dataset-residual ρ = −0.600.

**Maximum |ρ| on the fixed grid** (same rules as the Cldn4 max-|ρ| PR): also ρ = **−1.000**, exact p = **0.0833**, at Tacstd2 count ≥ 5 on the locked gate; dataset-residual ρ = **−1.000**. Every one of the 1800 eligible T-fraction tests that hit the |ρ|=1 ceiling is **negative** (0 positive).

Within each study the order agrees: higher Tacstd2% → lower T/NK (KL1 > KL2 on Tacstd2% and lower T; KL > K on Tacstd2% and lower T). Each study alone has only 2 mice.

IFN/APM is **not** sold as Tacstd2-cold here. Locked eight-mouse Tacstd2% vs epithelial IFN is ρ = **+0.643** (exact p = 0.096). The grid maximum |ρ| for IFN/APM is positive. That is reported in `FINDING_TACSTD2_MAXRHO.md` and is not inverted.

## 3. Tacstd2-high malignant/epithelial DEG → TJ

Locked-gate epithelium; Tacstd2 count > 0 vs = 0 within mouse. GSE180963_K falls below the high-cell floor (5 Tacstd2+ of 28 epi) and is inventoried, not scored. Seven mice score.

Frozen TJ_TISMO module (Cldn3/4/6/7, Cdh1, F11r, Ocln): **7/7** mice have Tacstd2-high > Tacstd2-low (mean Δ = +0.2608); **7/7** also have within-mouse one-sided MW p < 0.05; one-sided binomial on direction p = **0.007812**. TJ_EPITHELIAL (18) and CLDN4_TJ_EDGE (11) are the same 7/7 pattern.

Across-mouse mean pseudobulk gene delta ranks place TJ_TISMO above background (rank MW p = 3.6×10⁻⁹; hypergeometric top-500 k=5/7, enrichment 28.6, p = 1.9×10⁻⁷).

## Honest limits

- No private 8KL. No mega-merge of 8KL with these public mice.
- T-fraction n=4 cannot beat exact p = 0.0833.
- GSE154977 cannot enter T-fraction tests (FACS).
- Max-|ρ| is a grid extreme, not a single prespecified test; the locked pct>0 row is the prespecified T-fraction readout.
- IFN is not Tacstd2-cold on these public mice.
- Cell-level MW p inside a mouse is descriptive; mouse is the unit for binomial / ranking.

Private 8 KL mice used: **0**.
