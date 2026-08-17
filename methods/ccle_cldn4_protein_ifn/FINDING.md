# CCLE/DepMap lung protein — CLDN4 vs CD274 / IFN proteins

**Additive only.** Lung **RNA** Spearman ρ for CLDN4 vs Hallmark IFN-γ / MHC-I / CD274 is already reported in [`methods/depmap_cldn4_ifn`](https://github.com/jinxuanhong1-blip/sdaxcge/blob/cursor/depmap-cldn4-ifn-c8ad/methods/depmap_cldn4_ifn/FINDING.md) ([PR #304](https://github.com/jinxuanhong1-blip/sdaxcge/pull/304)): all-lung **n = 214**, ρ = **+0.284 / +0.217 / +0.330**. That RNA table is **not** re-estimated here. This page is the **protein extra**.

**Honest n.** Protein n is **not** 214. Public CCLE/DepMap protein matrices were inventoried; only Gygi/Nusinow TMT-MS can pair CLDN4 with CD274 or IFN proteins. Pairwise-complete n is reported on every row.

## Which protein table can do this pair?

| Matrix | CLDN4 protein | CD274 / PD-L1 protein | IFN / MHC-I proteins | Can compute CLDN4 vs CD274/IFN protein ρ? |
|---|---|---|---|---|
| **Gygi/Nusinow CCLE TMT-MS** (Cell 2020; `protein_quant_current_normalized.csv.gz`) | **Yes** (`sp\|O14493\|CLD4_HUMAN`) | **Yes** (`sp\|Q9NZQ7\|PD1L1_HUMAN`) | Yes (HLA-A/B/C, B2M, STAT1/2, MX1, ISG15, …; Hallmark IFN-γ **154/200** proteins with ≥8 lung values) | **Yes** |
| CCLE RPPA 20180123 | **No** (Claudin-7 only) | **No** (no PD-L1 / CD274 antibody) | STAT3, ADAR1 only among IFN-ish antibodies | **No** |
| ProCan-DepMapSanger 8498 (Gonçalves 2022) | **No** (CLDN1/3/7, not CLDN4) | **No** | Not scored (no CLDN4 to pair) | **No** |

Gygi lung cell lines mapped to DepMap 24Q4 `Model.csv` (`OncotreeLineage == Lung`, `ModelType == Cell Line`): **76** lines have any protein; **44** have CLDN4 protein (LUAD 22 / LUSC 7 / other NSCLC 6 / SCLC+NET 7 / other lung 2). NSCLC with CLDN4 protein **n = 35**. CD274 protein is missing in 6 of the 44 CLDN4+ lung lines, so the CD274 row is **n = 38**, not 44. These are cultured lines. There is **no T-cell infiltrate** and **no IFN treatment**. Values are TMT abundance, not surface FACS.

The older `_LUNG` suffix complete-case count for CLDN4+TACSTD2 was n=45. One Gygi `_LUNG` column does not map to a 24Q4 lung cell line, so the Model.csv-restricted n is **44**. Do not cite 45 and 44 as two discoveries.

## Finding (protein extra)

In cancer-cell **basal protein**, CLDN4 is **positively** associated with a Hallmark IFN-γ protein score and with CD274 protein on all-lung complete cases. The sign matches the RNA extra (PR #304) and is the opposite of “CLDN4-high = IFN-cold” at the cell-intrinsic level. MHC-I protein is weaker; its 95% CI includes 0. NSCLC-only CD274 protein is also compatible with 0.

| cohort | n (CLDN4 protein) | CLDN4 vs Hallmark IFN-γ protein ρ (p, q) | CLDN4 vs MHC-I protein ρ (p, q) | CLDN4 vs CD274 protein ρ (p, q) |
|---|---:|---|---|---|
| **all lung (primary)** | **44** | **+0.372** (n=44, p=0.013, q=0.033); CI [+0.07, +0.62] | **+0.253** (n=44, p=0.098, q=0.098); CI [−0.07, +0.55] | **+0.371** (n=**38**, p=0.022, q=0.033); CI [+0.06, +0.63] |
| NSCLC | 35 | +0.388 (n=35, p=0.021, q=0.064); CI [+0.08, +0.64] | +0.239 (n=35, p=0.17, q=0.18); CI [−0.12, +0.55] | +0.241 (n=**32**, p=0.18, q=0.18); CI [−0.11, +0.54] |
| LUAD | 22 | +0.374 (n=22, p=0.086) | +0.247 (n=22, p=0.27) | +0.340 (n=20, p=0.14) |
| LUSC | 7 | not tested (n<8) | not tested (n<8) | not tested (n<8) |
| SCLC+NET | 7 | not tested (n<8) | not tested (n<8) | not tested (n<8) |

Spearman, two-sided, pairwise complete. Hallmark IFN-γ protein and MHC-I protein (HLA-A, HLA-B, HLA-C, B2M; 4/4) are the mean of **within-cohort** protein-wise z-scores. A sample needs ≥5 Hallmark members (IFN score) or ≥2 MHC-I members. BH *q* is within the three new protein axes (IFN-γ score, MHC-I score, CD274) on that cohort. Bootstrap 95% CI: 5,000 resamples, seed 0, all-lung and NSCLC only. LUSC / SCLC+NET n=7 is underpowered; do not read those blanks as a histology-specific rule.

Do not cite the RNA n=214 as the protein n. Do not cite RPPA n=118. Those tables cannot score this pair.

### Sensitivity: compact ISG protein score

Same CLDN4+ lung lines. Compact ISG cassette = STAT1/2, IRF1/9, MX1, ISG15, IFIT1/3, OAS1/2/3, EIF2AK2, IFI35, BST2, SAMHD1, TRIM25, ADAR, TAP1/2, PSMB8/9/10 (proteins with ≥8 values in the cohort). Not in the three-axis FDR.

| cohort | n | CLDN4 vs compact ISG protein ρ |
|---|---:|---|
| all lung | 44 | +0.396 (p=0.0078); CI [+0.08, +0.67] |
| NSCLC | 35 | +0.478 (p=0.0037); CI [+0.19, +0.69] |

### Strongest single IFN proteins (all-lung, exploratory)

BH among the per-protein IFN/MHC list (not the scores). Companion TACSTD2 protein is listed only so it is not mistaken for a new IFN finding.

| protein | n | ρ vs CLDN4 protein | p | q (per-protein) |
|---|---:|---:|---:|---:|
| TACSTD2 (companion, not IFN) | 44 | +0.704 | 9.8×10⁻⁸ | — |
| MX1 | 44 | +0.561 | 7.4×10⁻⁵ | 0.0020 |
| OAS2 | 44 | +0.524 | 2.6×10⁻⁴ | 0.0035 |
| IFIT1 | 44 | +0.406 | 0.0063 | 0.057 |
| ISG15 | 44 | +0.381 | 0.011 | 0.069 |
| STAT1 | 44 | +0.204 | 0.18 | 0.28 |
| HLA-A | 44 | +0.150 | 0.33 | 0.45 |
| B2M | 44 | +0.220 | 0.15 | 0.25 |

MX1 / OAS2 are exploratory. Primary inference is the Hallmark IFN score, MHC-I score, and CD274 protein rows above.

Scatter: `figures/fig_cldn4_vs_ifn_mhc_protein.png`. Forest: `figures/fig_cldn4_protein_forest.png`.

## What this is not

- **Not** the RNA ρ. That is PR #304 (n=214). Protein n is 44 / 38.
- **Not** a hot/cold tumor, ICI-response, or IFN-stimulation result. No immune cells are in the dish.
- **Not** surface PD-L1 or MHC-I FACS. CD274 here is TMT abundance.
- **Not** RPPA n=118. RPPA has neither CLDN4 nor PD-L1.
- **Not** a claim that CLDN4 *induces* IFN or PD-L1. Association only.
- Hallmark IFN-γ includes proteins a pure cancer line may not use. The MHC-I cassette is cleaner and weaker (CI includes 0).

## How to rerun

```bash
python3 -m pip install -r methods/ccle_cldn4_protein_ifn/requirements.txt
python3 methods/ccle_cldn4_protein_ifn/download.py --outdir data/ccle_cldn4_protein_ifn
python3 methods/ccle_cldn4_protein_ifn/analyze.py --data data/ccle_cldn4_protein_ifn --outdir methods/ccle_cldn4_protein_ifn
```

Tables: `tables/correlations.tsv`, `tables/protein_coverage.tsv`, `tables/cohort_counts.tsv`, `tables/lung_protein_lines.tsv`, `tables/inventory.json`, `tables/key_stats.json`.
