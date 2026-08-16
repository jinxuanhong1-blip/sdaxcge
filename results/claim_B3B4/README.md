# Claims B3 + B4 — independent reproduction (honest)

Two user-stated claims (PPT 2026-08-17; `claims/B3.md`, `claims/B4.md`):

| Claim | Statement | User number | This reproduction | Verdict |
|-------|-----------|-------------|-------------------|---------|
| **B3** | TCGA-LUAD: TJ-high tumors have **low CD8 / low GEP** | `p < 1e-6` | Official 15-gene TJ vs CD8 Spearman **ρ = −0.29, p = 1.6e-11**; vs GEP **ρ = −0.28, p = 9.0e-11**. Survives ESTIMATE purity adjustment (**partial ρ = −0.26, p = 1.7e-9**). | **Supported** (negative association, well below 1e-6). Sign and significance **depend on the TJ gene set**: the broad KEGG pathway does **not** support it. |
| **B4** | GSE126044: non-responders have **higher TJ** than responders | `p = 0.019` | The **7-gene module named on the B3 claim page** (`CLDN1, CLDN4, CLDN7, F11R, TJP1, TJP2, OCLN`), scored as z-mean, gives Mann-Whitney **two-sided p = 0.019**. Broader TJ sets (15-gene, KEGG, GOCC) and ssGSEA of the same 7 genes do **not** reach 0.019. | **Reproduced only for that 7-gene z-mean + MWU.** Direction is often NR>R; significance is **gene-set- and test-dependent**. Small n (5 R / 11 NR). |

We did **not** invent a gene list to hit the target p-values. The 7-gene module is the list written on the official B3 claim page; the 15-gene set is the structural TJ signature already used in this repo (PR #69).

---

## Methods (pre-specified, not tuned)

**TJ definitions**

| Name | Genes | Role |
|------|-------|------|
| `TJ_claim_module` | `CLDN1, CLDN4, CLDN7, F11R, TJP1, TJP2, OCLN` | Named on `claims/B3.md` |
| `TJ_15` (primary structural) | `CLDN1/3/4/7, OCLN, TJP1/2/3, F11R, JAM2/3, MARVELD2/3, CGN, CGNL1` | Repo standard (PR #69) |
| `KEGG_TIGHT_JUNCTION` | 132 genes, MSigDB C2:CP:KEGG (hsa04530) | Named in the B3 test plan; actin/myosin-heavy |
| `GOCC_TIGHT_JUNCTION` | 138 genes, MSigDB GOCC | Sensitivity |

**Scoring** (both reported): (1) **z-mean** = mean per-gene z-score; (2) **ssGSEA** = GSVA-style Barbie 2009 area-under-walk, τ=0.25, unnormalized. Rank tests (Spearman, Mann-Whitney) are invariant to a later min-max rescale.

**Immune scores:** CD8 = `CD8A`/`CD8B`; CYT = `GZMA`/`PRF1` (Rooney); GEP = Ayers 2017 18-gene T-cell-inflamed profile (unweighted z-mean).

**Data**
- B3: UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` (log2 norm_count+1), **515** primary tumors (`-01`). Purity = Aran et al. *Nat Commun* 2015 ESTIMATE column, **511** matched samples.
- B4: GEO `GSE126044` raw counts → log2(CPM+1). Response labels from the GEO series matrix: **5 responders, 11 non-responders**.

**Tests:** Spearman + 2000× bootstrap 95% CI; median-split Wilcoxon (TJ-high vs TJ-low); rank-residual partial Spearman controlling for ESTIMATE purity; Mann-Whitney U (two-sided and one-sided NR>R) plus Welch t.

---

## B3 — TCGA-LUAD

### Spearman (n = 515)

| TJ | vs CD8 ρ (p) | vs GEP ρ (p) | vs CYT ρ (p) |
|----|-------------:|-------------:|-------------:|
| TJ_15 z-mean | **−0.291 (1.6e-11)** | **−0.281 (9.0e-11)** | −0.325 (3.6e-14) |
| TJ_15 ssGSEA | −0.265 (1.1e-9) | −0.220 (4.3e-7) | — |
| TJ_claim 7-gene z-mean | **−0.314 (3.3e-13)** | **−0.304 (2.0e-12)** | — |
| KEGG_TIGHT_JUNCTION z-mean | **+0.078 (0.078)** | +0.195 (8.1e-6) | — |
| KEGG_TIGHT_JUNCTION ssGSEA | **+0.072 (0.10)** | +0.162 (2.3e-4) | — |

`p < 1e-6` is **comfortably met** for the structural sets (TJ_15 and the 7-gene claim module), both as continuous Spearman and as TJ-high vs TJ-low Wilcoxon (CD8 MWU p = 2.0e-9 for TJ_15 z-mean). The association is **negative**: TJ-high = immunologically colder.

The **broad KEGG pathway does not support the claim vs CD8** (near-zero, wrong-sign). That is the main gene-set caveat.

### Purity (Aran ESTIMATE, n = 511)

CD8 tracks purity strongly (ρ = −0.56, p = 1.2e-43). TJ_15 tracks purity only weakly (ρ = +0.14, p = 0.0015). After rank-residual adjustment:

| Partial Spearman \| ESTIMATE purity | ρ | 95% CI | p |
|-------------------------------------|--:|--------|--:|
| TJ_15 z-mean vs CD8 | **−0.263** | (−0.338, −0.181) | **1.7e-9** |
| TJ_15 z-mean vs GEP | **−0.277** | (−0.355, −0.200) | **1.8e-10** |

Purity confounding is real but **does not explain** the TJ–CD8/GEP inverse association.

See `B3_continue_scatter.png`, `B3_continue_purity_partial.png`, `B3_correlation_table.csv`.

---

## B4 — GSE126044 (5 R / 11 NR)

| TJ / gene | direction | MWU two-sided p | MWU one-sided NR>R | Welch two-sided p |
|-----------|-----------|----------------:|-------------------:|------------------:|
| **TJ_claim 7-gene z-mean** | NR>R | **0.019** | 0.010 | 0.103 |
| TJ_claim 7-gene ssGSEA | NR>R | 0.115 | 0.057 | 0.129 |
| TJ_15 z-mean / ssGSEA | NR>R | 0.221 | 0.111 | — |
| GOCC_TIGHT_JUNCTION z-mean | NR>R | 0.377 | 0.189 | — |
| KEGG_TIGHT_JUNCTION z-mean / ssGSEA | NR≤R | 0.913 | 0.457 | — |
| **CLDN4 alone** (logCPM) | NR>R | 0.115 | 0.057 | 0.166 |

**What this means, honestly**

1. The user-reported **p = 0.019 is recovered** when — and only when — we use the 7-gene module written on the B3 claim page and a z-mean + two-sided Mann-Whitney U. All 7 genes are present in GSE126044. This was a **pre-specified** list, not a search over gene subsets.
2. The same 7 genes scored by **ssGSEA** give p = 0.115. The repo's 15-gene structural set gives p = 0.221. KEGG is null. **CLDN4 alone** is a trend (p = 0.115 two-sided / 0.057 one-sided), not 0.019.
3. Welch's t-test on the 7-gene z-mean is **not** significant (p = 0.10). The claimed p-value is the **non-parametric** one.
4. n = 5 vs 11 is tiny. A first-pass 36-config sweep over broader TJ sets never reached p < 0.05 (best one-sided p ≈ 0.06). So B4 is **fragile**: one reasonable definition matches the slide; most others do not.
5. The cohort is not simply underpowered for immune signal. Responders have **higher** CD8 (MWU p = 4.6e-4), CYT (p = 4.6e-4) and GEP (p = 0.0055). That is consistent with B3 (TJ-high ↔ low CD8/GEP ↔ non-response).
6. **Confound:** all 5 responders are fresh-frozen; all 5 FFPE samples are non-responders. Sample-prep batch is partly aliased with response.

See `B4_continue_boxplots.png`, `B4_group_table.csv`, `B4_sensitivity_grid.csv`.

---

## Reproducing

```bash
pip install -r requirements.txt
bash scripts/download_data.sh
python3 scripts/analyze_B3B4.py              # first-pass z-mean / KEGG REST
python3 scripts/sensitivity_and_figures.py   # 36-config B4 sweep
python3 scripts/analyze_B3B4_continue.py     # official 15-gene + 7-gene + ssGSEA + purity
```

## Files

| File | Contents |
|------|----------|
| `results_B3B4_continue.json` | Full continuation statistics |
| `B3_correlation_table.csv` / `B4_group_table.csv` | Flat tables |
| `TCGA_LUAD_scores_continue.csv` / `GSE126044_scores_continue.csv` | Per-sample scores |
| `B3_continue_scatter.png` / `B3_continue_purity_partial.png` / `B4_continue_boxplots.png` | Figures |
| `results_B3B4.json`, `B4_sensitivity_*.csv/json` | First-pass outputs (kept) |

## Bottom line

- **B3 holds.** Structural TJ (15-gene or the claim-page 7-gene module) is inversely associated with CD8/GEP/CYT in TCGA-LUAD at p ≪ 1e-6, and the inverse association survives ESTIMATE purity adjustment. The broad KEGG tight-junction pathway is the wrong gene set for this claim.
- **B4's p = 0.019 is real for one pre-specified definition** (7-gene z-mean + Mann-Whitney) and **not a general property** of "the TJ signature" in GSE126044. Report it as a small-n, definition-dependent result, not as a robust cohort-level finding.
