# MAX EFFECT Part2: TJ gene screen — maximize |ρ| / margin so CLDN4 wins

**Placement:** Part2-ready specificity screen (which TJ gene tracks immune-cold).
**Do not force into Part1** (Tacstd2 → TJ → immune entry). Never fabricate.

Pre-specified panel: **CLDN4 vs CLDN3, CLDN7, CLDN1, OCLN, F11R, CDH1**.
Layers: concordant-4 (n=65) + TCGA (8 cohorts) + CosMx He2022.

## Slide verdict (honest)

| Layer | Does CLDN4 clearly win? | Max-effect call |
|---|---|---|
| **Concordant-4** | **Point-estimate yes; not perm-clear vs CLDN3/7** | Max \|ρ\| among CLDN4-lead: **ρ=−0.531** (unadj %pos, I²=0%). Max margin: **0.144** (KRT18/19 partial, ρ=−0.482 vs CLDN3 −0.337). |
| **TCGA** | **No — 0/60 specs** | Anchor CD8∥KRT8/18/19: lead **F11R** (−0.131); CLDN4 **#5/6** (−0.081). |
| **CosMx** | **Only vs CDH1; sister claudins OFF** | Anchor q4/q1@50 µm: **CDH1** colder (0.786 vs 0.800). Max-margin: detected≥1 @25 µm CLDN4 **0.623** vs CDH1 0.695 (7/8, not 8/8). CLDN3/7/OCLN/F11R **OFF panel**. |

Locked CosMx CLDN4 cytotoxic ratios **0.36 / 0.52** (50/100 µm, 8/8, 5/5) are **not** recomputed or replaced. Locked concordant-4 CLDN4 %pos ρ=−0.531 is reproduced (max \|Δ %\| vs PR #644 ≈ 5×10⁻⁵).

## 1. Concordant-4 (maximize |ρ| and margin)

Grid: 54 specs (filters × panels × pct/mean × none/KRT/KRT1819). Eligible = CLDN4 most negative and negative in every cohort in the filter. **33 eligible.**

| Objective | Spec | CLDN4 ρ | Runner | Margin |
|---|---|---:|---|---:|
| **Max \|ρ\|** (eligible) | all4 / full7 / %pos / none | **−0.531** | CLDN7 −0.409 | 0.122 |
| **Max margin** (eligible) | all4 / tj6 / %pos / KRT18+19 | **−0.482** | CLDN3 −0.337 | **0.144** |

Max-|ρ| ranks (unadj %pos): CLDN4 −0.531 → CLDN7 −0.409 → CLDN3 −0.408 → CDH1 −0.401 → OCLN −0.208 → F11R −0.186 → CLDN1 +0.162.

Label-swap (KRT18/19 max-margin, 2000 shuffles): CLDN4 separable from **OCLN / F11R / CDH1** (perm p < 0.05), **not** from CLDN3 (p=0.087) or CLDN7 (p=0.179). Same pattern as PR #644 / #748: public scRNA puts CLDN4 first on the point estimate; it does not pin CLDN4 over sister claudins by permutation.

## 2. TCGA (0 eligible CLDN4 leads)

Grid: 60 specs (5 immune scores × 4 keratin sets × 3 cohort sets). Eligibility: CLDN4 most negative and negative in ≥7/8 (or ≥6/7 funnel / 2/2 lung).

**Eligible CLDN4-lead count = 0.** Lead-gene census: F11R 26, CLDN7 14, CDH1 9, OCLN 8, CLDN3 3, CLDN4 **0**.

Anchor (PR #748-like): CD8 score, KRT8+18+19, lung+funnel:

| rank | gene | meta ρ |
|---:|---|---:|
| 1 | F11R | −0.131 |
| 2 | CDH1 | −0.126 |
| 3 | CLDN7 | −0.109 |
| 4 | OCLN | −0.101 |
| **5** | **CLDN4** | **−0.081** |
| 6 | CLDN3 | −0.060 |

Closest CLDN4 gets is rank 3 (CYT, no keratin, lung+funnel) — still not #1. **Do not promote TCGA as a CLDN4-wins Part1 claim.**

## 3. CosMx (panel-limited max-effect)

On 960 / h5ad: **CLDN4, CDH1** (plus ESAM/CDH11 context). **OFF (not fabricated):** CLDN3, CLDN7, OCLN, F11R.

Grid: 16 (mode × radius) with both on-panel TJ genes usable. CLDN4 coldest vs CDH1: **5/16**. None of those five are 8/8 & 5/5 for CLDN4.

| Spec | CLDN4 ratio | CDH1 ratio | CLDN4 sections Δ&lt;0 |
|---|---:|---:|---|
| q4/q1 @50 µm (anchor) | 0.800 | **0.786** | 4/8 |
| **detected≥1 @25 µm (max margin)** | **0.623** | 0.695 | 7/8 |
| detected≥1 @50 µm | 0.776 | 0.801 | 5/8 |

Short-radius detection cuts can make CLDN4 colder than CDH1 on the mean ratio; they do not recover the locked 8/8 exclusion call, and they cannot speak to CLDN3/7.

## Part1 vs Part2

- **Part1 (entry):** Tacstd2 → TJ up (GSE137244) and locked CLDN4↔T/NK / CosMx exclusion. Unchanged.
- **Part2 (this PR):** gene-level TJ screen after that entry. Concordant-4 supports a **CLDN4 point-estimate lead**. TCGA and CosMx **do not** uniquely pin CLDN4 over the full panel. Searched maxima inflate margins — report as sensitivity, not a new lock.

Private KD co-culture remains what functionally pins CLDN4 (handoff); it is out of this public repo.

## Reproduce

```bash
# Concordant-4 (committed units; no GEO re-download)
python3 methods/tj_gene_screen_max_effect/scripts/sweep_concordant4.py

# TCGA (Xena GDC STAR TPM → $TCGA_CACHE)
TCGA_CACHE=/tmp/tcga python3 methods/tj_gene_screen_max_effect/scripts/sweep_tcga.py

# CosMx h5ad (figshare 25976224; gitignored)
python3 methods/tj_gene_screen_max_effect/scripts/download_cosmx.py
COSMX_H5AD=methods/tj_gene_screen_max_effect/data/cosmx_nsclc/cosmx_human_nsclc_clustered.h5ad \
  python3 methods/tj_gene_screen_max_effect/scripts/sweep_cosmx.py

python3 methods/tj_gene_screen_max_effect/scripts/make_verdict.py
```

PPT paste: `ppt/PPT_SLIDE.md`.
