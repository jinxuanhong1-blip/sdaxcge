# Finding — pair GSE123902+GSE189357, CLDN4-only REAL Slingshot/PAGA

ADDITIVE. **CLDN4 only.** No dual-high. No CellChat.

The pair that already **differs** is **taken as given** and is not re-audited (PR #459):

- **GSE123902 + GSE189357 %pos · n=22 · ρ=−0.638 · p=0.003 · I²=0%**
- Members: GSE123902 marker-malignant donors n=13 + GSE189357 marker-malignant patients n=9
- Between-patient Q4 vs Q1 on that same vector is **r=−1.000 on tails 7/5 — thin**. That is not the trajectory n.

This folder asks a **different** question: where do **CLDN4**, a CLDN4-excluded **barrier/keratin** score, and an **IFN** score sit on a **real Slingshot** lineage (PAGA for geometry) of marker epithelium from this pair.

Primary clock: **Bioconductor slingshot 2.10.0** (Street et al. 2018). PAGA is scanpy (Wolf et al. 2019). Inferential unit = **patient/donor** on the given n=22. Cell-level ρ is descriptive. Barrier/keratin **excludes CLDN4**. Root = GSE123902 NORMAL, not CLDN4-high, median AT2 score. Root is **not** CLDN4-high.

**What holds (given units, Spearman n=22).** CLDN4 tracks a CLDN4-excluded barrier/keratin score (n=22, ρ=0.649, p=0.00109, q=0.009). Within-unit CLDN4-high vs low barrier Δmed = +0.568, p=1.91e-06 (n=20). The Slingshot root is real in the sense that it is GSE123902 NORMAL, CLDN4-low (expr=0), Leiden 8.

**What does not hold / is thin.** Pooled given-unit CLDN4 vs Slingshot PT is null (n=22, ρ=0.292, p=0.187). CLDN4 vs IFN is null (n=22, ρ=0.045, p=0.844). Barrier vs PT and IFN vs PT are null. AT2 vs PT is **positive** (n=22, ρ=0.494, p=0.0194) — this primary lineage is **not** an AT2-declining clock. SFTPC/SFTPA1 were absent from the gene intersection, so the AT2 score is incomplete (SFTPB/NAPSA/LAMP3/ABCA3 only). Q4 vs Q1 tails 7/5 from PR #459 stay thin and are not re-used as a trajectory test. Cell-level lineage ρ uses thousands of cells and is descriptive.

## Verdict

Real Slingshot (2.10.0) produced 12 lineage(s) rooted on Leiden 8 (GSE123902 NORMAL, not CLDN4-high, median AT2 score; root tertile=low, CLDN4=0.000). Primary lineage Lineage1. Given-unit CLDN4 vs Slingshot PT: n=22, ρ=0.292, p=0.187. CLDN4 vs barrier/keratin (no CLDN4): n=22, ρ=0.649, p=0.00109. CLDN4 vs IFN (Hallmark IFNα): n=22, ρ=0.045, p=0.844. Barrier vs PT: n=22, ρ=-0.077, p=0.732. IFN vs PT: n=22, ρ=-0.270, p=0.223. AT2 vs PT (control): n=22, ρ=0.494, p=0.0194. PAGA has 1 component(s) at connectivity>0 among 18 Leiden vertices. Given combo n=22 is not re-audited. Q4 vs Q1 tails 7/5 are thin. No dual-high gate. Not a T/NK redo.

## Honest n

- Given combo (do not re-audit): **n=22** (13+9). Q4 vs Q1 tails **7/5 are thin**.
- Analysis cells after QC (capped ≤350/unit): **n_cells = 6514** (GSE123902 3419, GSE189357 3095).
- Given tumor/met units in the object: **n_given_units = 22** (GSE123902 13, GSE189357 9).
- Units with ≥10 cells used for Spearman: **n = 22**.
- GSE123902 NORMAL root-pool cells / units: 755 / 4.
- CLDN4 tertile cells: low 3132, mid 1211, high 2171.
- Units with ≥8 cells in both CLDN4-high and CLDN4-low arms: **n = 20**.
- Genes absent from locked sets: {'AT2': ['SFTPC', 'SFTPA1'], 'AT1': [], 'club': ['SCGB1A1', 'SCGB3A2', 'SCGB3A1'], 'basal': ['KRT5', 'TP63', 'NGFR'], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': ['CEACAM5', 'CEACAM6'], 'ifn_isg': [], 'IFN_hallmark_IFNa': ['BATF2', 'CXCL11', 'RNF31', 'TENT5A', 'WARS1'], 'single_genes': ['SFTPC', 'SCGB1A1', 'KRT5', 'PTPRC', 'PTPRC']}.
- Slingshot: available=True; version=2.10.0; n_lineages=12; start_cluster=8.
- PAGA components at connectivity>0: **1** among 18 Leiden vertices.
- Catalog before cap is in `results/tables/extract_inventory.json`. Thin libraries (e.g. LX699 46 malignant) stay in the given n=22 and are flagged, not dropped.

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- Batch: harmonypy on PCA, batch=dataset.
- Slingshot / DPT root: GSE123902 NORMAL, not CLDN4-high, median AT2 score (root cell GSE123902:LX682_NORMAL, Leiden 8). Never CLDN4-high.
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).
- IFN score: mean of Hallmark interferon-alpha response genes present (CLDN4 not in the set). Compact ISG panel is extra.
- Marker-epithelial gate (same as PR #459): (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0.

## Lineage table (done criterion)

Real Slingshot lineages. Start cluster is the root Leiden (NORMAL / not CLDN4-high). Cell-level ρ along each lineage is descriptive. Unit-level tests are in the Spearman table.

| lineage | start | end | n_clusters | n_cells | path | ρ CLDN4 | ρ barrier | ρ IFN |
| --- | --- | --- | ---: | ---: | --- | ---: | ---: | ---: |
| Lineage1 | 8 | 2 | 5 | 2698 | `8->1->13->0->2` | 0.422 | 0.279 | 0.145 |
| Lineage2 | 8 | 3 | 5 | 2542 | `8->1->13->0->3` | 0.612 | 0.523 | 0.176 |
| Lineage3 | 8 | 11 | 4 | 1751 | `8->1->12->11` | 0.456 | 0.663 | 0.323 |
| Lineage4 | 8 | 6 | 3 | 951 | `8->5->6` | -0.018 | 0.042 | 0.366 |
| Lineage5 | 8 | 15 | 3 | 1713 | `8->1->15` | 0.417 | 0.663 | -0.08 |
| Lineage6 | 8 | 16 | 3 | 1466 | `8->1->16` | 0.395 | 0.669 | 0.021 |
| Lineage7 | 8 | 17 | 3 | 1467 | `8->1->17` | 0.468 | 0.701 | 0.163 |
| Lineage8 | 8 | 4 | 2 | 962 | `8->4` | -0.11 | 0.251 | -0.293 |
| Lineage9 | 8 | 7 | 2 | 1025 | `8->7` | 0.08 | 0.289 | 0.291 |
| Lineage10 | 8 | 9 | 2 | 857 | `8->9` | 0.682 | 0.635 | 0.049 |
| Lineage11 | 8 | 10 | 2 | 692 | `8->10` | 0.021 | 0.094 | -0.026 |
| Lineage12 | 8 | 14 | 2 | 678 | `8->14` | 0.049 | 0.072 | -0.404 |

Table: `results/tables/slingshot_lineages.tsv`. Primary lineage for unit-level tests: **Lineage1**.

## Primary (given-unit Spearman, BH inside this list)

Patient/donor is the unit. Given n=22; Spearman n is units with ≥10 cells after QC/cap. Q4 vs Q1 tails 7/5 remain **thin** and are not this table.

| Contrast | n_units | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs Slingshot PT | 22 | 0.292 | 0.187 | 0.447 |
| CLDN4 vs DPT (companion) | 22 | 0.016 | 0.942 | 1 |
| CLDN4 vs AT2 score | 22 | 0.518 | 0.0136 | 0.0543 |
| CLDN4 vs barrier/keratin (no CLDN4) | 22 | 0.649 | 0.00109 | 0.00917 |
| CLDN4 vs IFN (Hallmark IFNα) | 22 | 0.045 | 0.844 | 1 |
| CLDN4 vs compact ISG (extra) | 22 | -0.233 | 0.296 | 0.466 |
| CLDN4 vs malignant-like score | 22 | -0.226 | 0.311 | 0.466 |
| CLDN4 vs TACSTD2 (comparator) | 22 | 0.634 | 0.00153 | 0.00917 |
| barrier vs Slingshot PT | 22 | -0.077 | 0.732 | 0.976 |
| IFN vs Slingshot PT | 22 | -0.270 | 0.223 | 0.447 |
| SFTPC vs Slingshot PT (control) | 0 | NA | NA | NA |
| AT2 score vs Slingshot PT (control) | 22 | 0.494 | 0.0194 | 0.0583 |

## Sensitivity (not in the BH family)

| Contrast | n_units | ρ | p |
| --- | ---: | ---: | ---: |
| GSE123902-only CLDN4 vs Slingshot PT | 13 | 0.198 | 0.517 |
| GSE189357-only CLDN4 vs Slingshot PT | 9 | 0.167 | 0.668 |
| GSE123902-only CLDN4 vs barrier/keratin (no CLDN4) | 13 | 0.643 | 0.0178 |
| GSE189357-only CLDN4 vs barrier/keratin (no CLDN4) | 9 | 0.700 | 0.0358 |
| GSE123902-only CLDN4 vs IFN | 13 | 0.033 | 0.915 |
| GSE189357-only CLDN4 vs IFN | 9 | 0.083 | 0.831 |
| GSE123902-only barrier vs Slingshot PT | 13 | -0.264 | 0.384 |
| GSE189357-only barrier vs Slingshot PT | 9 | 0.250 | 0.516 |
| GSE123902-only IFN vs Slingshot PT | 13 | -0.093 | 0.762 |
| GSE189357-only IFN vs Slingshot PT | 9 | -0.350 | 0.356 |

## Extra figure — CLDN4-high vs CLDN4-low (unit-paired)

Paired tertile n=20. Thin arms are possible on GSE123902 mets (LX699, LX701).

| Paired contrast (high − low) | n_units | Δ median | p |
| --- | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 20 | 0.568 | 1.91e-06 |
| IFN (Hallmark IFNα) high vs low | 20 | 0.059 | 3.62e-05 |
| AT2 high vs low | 20 | 0.692 | 1.91e-06 |
| Slingshot PT high vs low | 19 | 17.084 | 3.81e-06 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- The given T/NK %pos Spearman (n=22, ρ=−0.638) is **not re-audited**.
- Q4 vs Q1 tails **7/5 are thin**. Do not treat r=−1 as a trajectory result.
- Marker-epithelial is **not** author malignant and **not** CNV.
- GSE123902 NORMAL cells are a root pool, not extra inferential units.
- No TACSTD2∩CLDN4 both-high gate.
- Do not write “AT2 differentiates into LUAD because PAGA is connected.”
- Do not write “CLDN4 marks the malignant terminal.”
- Slingshot is an ordering on this object, not a developmental clock. Twelve lineages; several are two-cluster stubs.
- AT2 vs primary-lineage PT is positive. Do not write that AT2 falls along this clock.
- SFTPC was missing from the cohort gene intersection. Club genes were also missing.
- Not ICI / MPR / RECIST. Not CellChat.

## Outputs

- `results/tables/slingshot_lineages.tsv` — **done criterion**
- `results/tables/sample_level_spearman.tsv`
- `results/tables/sample_means.tsv`
- `results/tables/extract_inventory.json`
- `results/figures/fig_trajectory_cldn4.png`
- `results/figures/fig_extra_along_pseudotime.png`
- `results/figures/fig_extra_cldn4_tertile.png`
- `results/figures/fig_extra_umap_batch_at2.png`
- `results/figures/fig_honest_n.png`
- `results/summary.json`

## Reproduce

```bash
pip install -r methods/pair_123902_189357_slingshot_cldn4/requirements.txt
bash methods/pair_123902_189357_slingshot_cldn4/scripts/install_r_slingshot.sh
python3 methods/pair_123902_189357_slingshot_cldn4/scripts/download.py \
  --out /tmp/geo_pair_123902_189357
python3 methods/pair_123902_189357_slingshot_cldn4/scripts/extract.py \
  --tars /tmp/geo_pair_123902_189357 \
  --out /tmp/geo_pair_123902_189357/epithelium.h5ad
python3 methods/pair_123902_189357_slingshot_cldn4/scripts/analyze.py \
  --input /tmp/geo_pair_123902_189357/epithelium.h5ad \
  --outdir methods/pair_123902_189357_slingshot_cldn4/results \
  --finding methods/pair_123902_189357_slingshot_cldn4/FINDING.md
```

