# Finding — pair GSE123902+GSE205335, CLDN4-only REAL Palantir destinies

ADDITIVE. **CLDN4 only.** No TACSTD2∩CLDN4 dual-high gate. No GSE148071.
This is the pair that **differs**: malignant-cell IFN family DE is already
**logFC = −1.05** (Q4 vs Q1, n=9/9, p=3.16e-4; pair IFN DE folder). That DE
is taken as given and is not re-run. This folder asks whether CLDN4, a
CLDN4-excluded barrier score, and IFN travel together along **real Palantir
destinies** (Setty et al. 2019).

Engine: **palantir 1.4.5** (`anndata_sanitized_multiscale` API). Not DPT-as-Palantir. Inferential unit = **donor (GSE123902) / patient (GSE205335)**. Cell-level ρ is descriptive. Barrier score **excludes CLDN4**. Root rule: GSE123902 NORMAL, CLDN4-low tertile, AT2 ≥ pool median, farthest from malignant centroid in Harmony/PCA. Not CLDN4-high. (cell `GSE123902:LX684:NORMAL:195622403159966`, CLDN4=0.000, AT2=0.300). Root is **not** CLDN4-high.

**What this folder tests (tumor units n=35).** CLDN4 vs IFN n=35, ρ=-0.136, p=0.437; CLDN4 vs barrier n=35, ρ=0.396, p=0.0187; IFN vs Palantir PT n=35, ρ=-0.377, p=0.0257. Within-unit CLDN4-high vs low IFN n=32, W=74.0, Δmed=0.049, p=0.000178. The given pair IFN DE (−1.05, 9 vs 9) is a different contrast (malignant pseudobulk Q4 vs Q1) and is not re-audited here.

## Verdict

REAL Palantir v1.4.5 on GSE123902+GSE205335 epithelium. Root not CLDN4-high (GSE123902 NORMAL, CLDN4-low tertile, AT2 ≥ pool median, farthest from malignant centroid in Harmony/PCA. Not CLDN4-high.; CLDN4=0.000). Tumor-unit CLDN4 vs Palantir PT: n=35, ρ=-0.188, p=0.28. CLDN4 vs barrier (CLDN4 excluded): n=35, ρ=0.396, p=0.0187. CLDN4 vs IFN: n=35, ρ=-0.136, p=0.437. barrier vs IFN: n=35, ρ=0.356, p=0.0357. IFN vs Palantir PT: n=35, ρ=-0.377, p=0.0257. Paired CLDN4-high vs low IFN: n=32, W=74.0, Δmed=0.049, p=0.000178. Paired CLDN4-high vs low barrier: n=32, W=0.0, Δmed=0.585, p=4.66e-10. 2 specified terminal(s): ['AT1', 'malignant']. Pair IFN DE −1.05 is given and not re-run. No dual-high. No GSE148071. Donor/patient is the unit. Destiny tables written.

## Honest n

- Analysis cells after QC (capped ≤250/unit): **n_cells = 8113** (GSE123902 2942, GSE205335 5171).
- Units (GSE123902 donor + GSE205335 patient; NORMAL root-pool units separate): **n_units = 39** (tumor 35, NORMAL-root 4; GSE123902 17, GSE205335 22).
- Tumor units with ≥10 epithelial cells used for Spearman: **n = 35**.
- GSE123902 NORMAL epithelial cells in the object (root pool): **739**.
- Root cell CLDN4 = 0.000 vs object CLDN4-high tertile cut 1.660 (root below the cut).
- Palantir terminals (specified malignant + AT1, not CLDN4-high): **2** — ['AT1', 'malignant'].
- Units with ≥8 cells in both CLDN4-high and CLDN4-low arms: **n = 32**.
- Genes absent from locked sets: {'IFN': ['MARCHF1', 'TENT5A', 'WARS1']}.
- GSE148071 not used. Dual-high not used. GSE205335 normal-tissue epithelium dropped.
- Palantir: available=True; version=1.4.5; waypoints=500.

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30; diffusion components 10; waypoints 500.
- Batch: harmonypy on PCA, batch=dataset.
- DPT is **not** the clock. Palantir destinies are.
- Barrier genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).
- IFN score: Hallmark IFNα ∪ IFNγ (same family as the given −1.05 DE).
- Terminals are pre-specified (malignant farthest from root; max AT1 among non-malignant). **Not** the CLDN4-high quantile.

## Primary (tumor-unit Spearman, BH inside this list)

| Contrast | n_units | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs Palantir pseudotime | 35 | -0.188 | 0.28 | 0.391 |
| CLDN4 vs barrier (no CLDN4) | 35 | 0.396 | 0.0187 | 0.12 |
| CLDN4 vs IFN | 35 | -0.136 | 0.437 | 0.51 |
| barrier vs IFN | 35 | 0.356 | 0.0357 | 0.125 |
| CLDN4 vs AT2 score | 35 | 0.078 | 0.657 | 0.707 |
| IFN vs Palantir pseudotime | 35 | -0.377 | 0.0257 | 0.12 |
| barrier vs Palantir pseudotime | 35 | -0.298 | 0.0817 | 0.229 |
| CLDN4 vs entropy | 35 | 0.155 | 0.375 | 0.478 |
| CLDN4 vs fate_AT1 | 35 | 0.461 | 0.00536 | 0.075 |
| IFN vs fate_AT1 | 35 | -0.043 | 0.807 | 0.807 |
| barrier vs fate_AT1 | 35 | 0.282 | 0.1 | 0.234 |
| CLDN4 vs fate_malignant | 35 | 0.204 | 0.24 | 0.374 |
| IFN vs fate_malignant | 35 | -0.255 | 0.14 | 0.28 |
| barrier vs fate_malignant | 35 | -0.238 | 0.169 | 0.296 |

## Destinies (CLDN4, barrier, IFN along each terminal)

Cell assignment = argmax Palantir fate probability. Trends are bin means of
unsmoothed log-normalized scores vs Palantir pseudotime on that destiny.
Unit-level Spearman is the claim; cell-level is descriptive.

- CLDN4 vs barrier (no CLDN4): n=35, ρ=0.396, p=0.0187
- CLDN4 vs IFN: n=35, ρ=-0.136, p=0.437
- barrier vs IFN: n=35, ρ=0.356, p=0.0357
- IFN vs Palantir pseudotime: n=35, ρ=-0.377, p=0.0257
- barrier vs Palantir pseudotime: n=35, ρ=-0.298, p=0.0817
- CLDN4 vs fate_AT1: n=35, ρ=0.461, p=0.00536
- IFN vs fate_AT1: n=35, ρ=-0.043, p=0.807
- barrier vs fate_AT1: n=35, ρ=0.282, p=0.1
- CLDN4 vs fate_malignant: n=35, ρ=0.204, p=0.24
- IFN vs fate_malignant: n=35, ρ=-0.255, p=0.14
- barrier vs fate_malignant: n=35, ρ=-0.238, p=0.169

## Sensitivity (not in the BH family)

| Contrast | n_units | ρ | p |
| --- | ---: | ---: | ---: |
| GSE123902-only CLDN4 vs pseudotime | 13 | -0.126 | 0.681 |
| GSE123902-only CLDN4 vs IFN | 13 | -0.220 | 0.471 |
| GSE123902-only CLDN4 vs barrier | 13 | 0.659 | 0.0142 |
| GSE123902-only IFN vs pseudotime | 13 | 0.214 | 0.482 |
| GSE205335-only CLDN4 vs pseudotime | 22 | -0.223 | 0.318 |
| GSE205335-only CLDN4 vs IFN | 22 | 0.012 | 0.958 |
| GSE205335-only CLDN4 vs barrier | 22 | 0.299 | 0.177 |
| GSE205335-only IFN vs pseudotime | 22 | -0.606 | 0.0028 |
| all-units incl NORMAL-root CLDN4 vs IFN | 39 | -0.132 | 0.422 |
| drop SCLC/NUT CLDN4 vs IFN | 30 | -0.050 | 0.793 |
| drop SCLC/NUT CLDN4 vs barrier | 30 | 0.560 | 0.00128 |
| drop SCLC/NUT CLDN4 vs pseudotime | 30 | -0.455 | 0.0116 |

## Extra figure — CLDN4-high vs CLDN4-low (unit-paired)

Emitted: **True**. Observed tumor-unit Spearman(CLDN4, IFN) n=35, ρ=-0.136, p=0.437.

| Paired contrast (high − low) | n_units | Δ median | p |
| --- | ---: | ---: | ---: |
| barrier (no CLDN4) high vs low | 32 | 0.585 | 4.66e-10 |
| IFN high vs low | 32 | 0.049 | 0.000178 |
| Palantir PT high vs low | 32 | -0.016 | 0.0156 |
| AT2 high vs low | 32 | 0.133 | 8.1e-07 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- This is not a redo of the pair IFN DE (−1.05 is given).
- Palantir destinies on a two-cohort merge mix protocol and tissue; they are not a clock.
- The kNN graph had unreachable cells; one destiny can absorb most cells. That is reported, not hidden.
- Terminals were not defined as CLDN4-high. Post-hoc labels are descriptive.
- Barrier score excludes CLDN4. No TACSTD2∩CLDN4 both-high gate.
- GSE148071 is not in this object.
- Not evidence that CLDN4 *causes* IFN or barrier change.
- GSE123902 malignant/epithelium is marker-gated (no author AT2/malignant labels).

## Outputs

- `results/tables/destiny_cell.tsv` — **done criterion (cell destinies)**
- `results/tables/destiny_unit.tsv` — **done criterion (donor/patient destinies)**
- `results/tables/destiny_terminals.tsv`
- `results/tables/destiny_trends.tsv` — CLDN4 / barrier / IFN along destinies
- `results/tables/destiny_unit_spearman.tsv`
- `results/tables/honest_n.tsv`
- `results/figures/fig_destiny_trends.png`
- `results/figures/fig_unit_destiny.png`
- `results/figures/fig_extra_cldn4_tertile.png`
- `results/figures/fig_extra_root_not_cldn4high.png`
- `results/figures/fig_honest_n.png`
- `results/summary.json`

## Reproduce

```bash
pip install -r methods/pair_123902_205335_palantir_cldn4/requirements.txt
python3 methods/pair_123902_205335_palantir_cldn4/scripts/download.py \
  --outdir /tmp/pair_123902_205335_palantir
python3 methods/pair_123902_205335_palantir_cldn4/scripts/extract.py \
  --data /tmp/pair_123902_205335_palantir \
  --out /tmp/pair_123902_205335_palantir/epithelium.h5ad
python3 methods/pair_123902_205335_palantir_cldn4/scripts/analyze.py \
  --input /tmp/pair_123902_205335_palantir/epithelium.h5ad \
  --outdir methods/pair_123902_205335_palantir_cldn4/results \
  --finding methods/pair_123902_205335_palantir_cldn4/FINDING.md
```

