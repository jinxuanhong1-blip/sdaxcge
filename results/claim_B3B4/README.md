# Claims B3 + B4 — independent reproduction (honest assessment)

This directory contains an **independent** attempt to reproduce two claims about a
tight-junction (TJ) gene-expression signature. The stated goal was honesty: the
numbers below are what the data actually give, whether or not they match the
values that were claimed.

| Claim | Statement | User's stated value | This reproduction | Verdict |
|-------|-----------|---------------------|-------------------|---------|
| **B3** | In TCGA-LUAD, the TJ signature is associated with CD8 and with the GEP signature | `p < 1e-6` | Spearman **p = 7.1e-9 (CD8)**, **p = 1.7e-9 (GEP)** | **Supported** (significance), but the association is **negative** and sign is definition-dependent |
| **B4** | In GSE126044, non-responders (NR) have a **higher** TJ signature than responders (R) | `p = 0.019` | Direction NR>R reproduces; **Mann-Whitney p = 0.22 two-sided / 0.11 one-sided** | **Direction supported, significance NOT reproduced** |

## What "TJ signature" means here (important caveat)

The exact gene membership of the "TJ signature" was **not provided** with the
task, and the signature choice materially changes the answer. To stay honest and
reproducible we did not invent a bespoke list to hit the target p-values; instead
we anchored to citable sources and reported sensitivity:

- **`TJ_core` (primary):** curated core *structural* tight-junction components —
  claudins (`CLDN1/2/3/4/5/7/8/10/11/15/18/23`), occludin (`OCLN`), MARVEL
  proteins (`MARVELD2/3`), zonula occludens (`TJP1/2/3`), junctional adhesion
  molecules (`F11R/JAM2/JAM3`), and polarity/scaffold proteins
  (`CGN, CGNL1, CRB3, PARD3, PARD6A, MPDZ, PATJ, SYMPK`).
- **`TJ_kegg` (sensitivity):** the full KEGG "Tight junction" pathway
  (`hsa04530`), fetched live from the KEGG REST API. This set is broad and
  contains many actin/myosin/tubulin genes, so it is *less* TJ-specific.
- **`TJ_claudin` (sensitivity):** claudin-family subset only.

**CD8** = mean of `CD8A`, `CD8B`. **GEP** = the 18-gene T-cell-inflamed Gene
Expression Profile (Ayers et al., *J Clin Invest* 2017). Signature score =
mean per-gene z-score (each gene standardized across samples) — a transparent,
widely used approach.

## B3 — TCGA-LUAD (n = 515 primary tumors)

Data: UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` (log2 norm_count+1, RSEM),
restricted to primary tumors (barcode sample-type `-01`).

| TJ definition | vs | Spearman ρ | Spearman p | Pearson r | Pearson p |
|---------------|----|-----------:|-----------:|----------:|----------:|
| TJ_core | CD8 | **−0.25** | **7.1e-9** | −0.22 | 7.3e-7 |
| TJ_core | GEP | **−0.26** | **1.7e-9** | −0.21 | 1.3e-6 |
| TJ_kegg | CD8 | +0.07 | 0.095 (n.s.) | +0.08 | 0.055 |
| TJ_kegg | GEP | +0.23 | 8.5e-8 | +0.26 | 1.9e-9 |

**Assessment:** The claim `p < 1e-6` is **comfortably met** for the primary
`TJ_core` signature (Spearman) against both CD8 and GEP. Two honest caveats:
1. **The association is negative** — higher tight-junction expression tracks with
   *lower* CD8/GEP (i.e., a more immunologically "cold" tumor). If the original
   claim implied a positive association, the sign disagrees.
2. **It is signature-dependent.** The broad KEGG set is essentially null against
   CD8 and flips positive against GEP, so the exact gene list matters.

See `B3_TCGA_LUAD_scatter.png`.

## B4 — GSE126044 (5 responders, 11 non-responders)

Data: GEO `GSE126044` raw counts (pre-treatment NSCLC, anti-PD-1; Cho et al.),
per-sample response labels taken directly from the GEO series matrix.
Normalization log2(CPM+1); test Mann-Whitney U.

| TJ definition | median NR | median R | direction | MWU p (two-sided) | MWU p (one-sided NR>R) |
|---------------|----------:|---------:|-----------|------------------:|-----------------------:|
| TJ_core | 0.100 | −0.039 | NR>R | 0.221 | 0.111 |
| TJ_kegg | 0.095 | −0.001 | NR>R | 0.441 | 0.220 |

**Assessment:** The **direction reproduces** (non-responders have a higher TJ
score), but the effect is **not statistically significant**, and we could **not
reproduce `p = 0.019`**.

A 36-configuration robustness sweep (3 TJ definitions × {all 16 / fresh-only} ×
{logCPM / logcounts} × {z-mean / mean-expr / rank-mean scoring}) found:
- NR>R direction in **36/36** configurations (robust direction), but
- **best two-sided p = 0.126**, **best one-sided p = 0.063**;
- **no** configuration reached p < 0.05, and **none** landed near 0.019.

See `B4_sensitivity_grid.csv`, `B4_sensitivity_summary.json`,
`B4_GSE126044_boxplot.png`.

### Why B4 likely differs from the claimed value
The cohort is *not* simply too small to detect immune signal. As a sanity check,
canonical immune signatures separate the groups strongly and in the expected
direction (responders higher):

| Signature | median R | median NR | MWU p (two-sided) |
|-----------|---------:|----------:|------------------:|
| CD8 | 1.25 | −0.49 | **0.00046** |
| GEP | 0.62 | −0.19 | **0.0055** |

So GSE126044 clearly detects that responders are more CD8/GEP-inflamed. The TJ
signature simply does not reach significance here under any standard choice we
tried. Reproducing `p = 0.019` would require the original TJ gene list and exact
scoring/test (e.g., a specific ssGSEA implementation, a one-sided test, or a
particular gene subset). Also note a **confound**: all 5 responders are
fresh-frozen samples while all 5 FFPE samples are non-responders, so sample-prep
batch is partly aliased with response.

## Reproducing

```bash
python3 scripts/analyze_B3B4.py            # B3 + B4 core results -> results_B3B4.json
python3 scripts/sensitivity_and_figures.py # B4 sweep, sanity check, figures
```

Inputs are downloaded to `data/` (UCSC Xena TCGA-LUAD, GEO GSE126044 counts +
series matrix, KEGG hsa04530).

## Files
- `results_B3B4.json` — full B3/B4 statistics and matched/missing signature genes.
- `TCGA_LUAD_signature_scores.csv`, `GSE126044_signature_scores.csv` — per-sample scores.
- `B4_sensitivity_grid.csv`, `B4_sensitivity_summary.json` — robustness sweep + immune sanity check.
- `B3_TCGA_LUAD_scatter.png`, `B4_GSE126044_boxplot.png` — figures.

## Bottom line
- **B3 is supported** (highly significant, `p < 1e-6`), with the honest caveats
  that the correlation is **negative** and depends on the TJ gene set.
- **B4 is not reproduced as stated**: the NR>R direction is robust, but the
  effect is not significant (best one-sided `p ≈ 0.06`) and we could not recover
  `p = 0.019` with any standard signature/scoring/test. The two findings are
  biologically consistent (TJ-high ↔ low CD8/GEP ↔ non-response), but B4's
  statistical claim overstates what this dataset supports.
