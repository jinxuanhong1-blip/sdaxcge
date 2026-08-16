# Public lung scRNA: B-cell / TLS-like vs malignant TACSTD2 / CLDN4

**Slice:** `methods|scripts|results/scrna_tls_meta/` only.
**Additive.** User A6 / B6 immune-cold is taken as given and is not re-cut.
**ICI labels:** not the estimand. This slice cannot test ICI benefit or histologic TLS.

### TL;DR (per-feature, not a pooled overclaim)

- TACSTD2 vs B fraction RE-meta: k=7 n=162 ρ=-0.110 (95% CI -0.389 to +0.188) p=0.47 τ²=0.102. Heterogeneous; not a general inverse.
- TACSTD2 vs TLS12 RE-meta: k=7 n=162 ρ=+0.078 (95% CI -0.175 to +0.322) p=0.549 τ²=0.059. Not a general inverse.
- CLDN4 vs B fraction RE-meta: k=7 n=162 ρ=-0.092 (95% CI -0.329 to +0.157) p=0.472 τ²=0.055.
- CLDN4 vs TLS12 RE-meta: k=7 n=162 ρ=-0.129 (95% CI -0.379 to +0.138) p=0.344 τ²=0.072.
- A6/B6 T-cell immune-cold is taken as given and is not tested here.

### What was tested

Patient-level Spearman of malignant (or epithelial / author gate) TACSTD2 and CLDN4 versus B-cell fraction and the 12-chemokine TLS score, then DerSimonian–Laird random-effects meta. Public processed matrices only.

### Coverage (honest missingness)

| Cohort | Status | Patients eligible | Cells (tumor) | TACSTD2 | CLDN4 | TLS12 k/12 | B cells | Note |
|---|---|---:|---:|---|---|---:|---:|---|
| GSE131907 | included | 32 | 128065 | yes | yes | 12 | 16439 | epithelial; malignant |
| GSE148071 | included | 41 | 89887 | yes | yes | 12 | 662 | epithelial; malignant |
| GSE154826 | included | 30 | 193699 | yes | yes | 12 | 20989 | gate |
| GSE207422 | included | 15 | 92330 | yes | yes | 12 | 7648 | epithelial; malignant |
| GSE241934_IIT | included | 11 | 78691 | yes | yes | 12 | 13056 | epithelial |
| GSE241934_RWC | included | 24 | 229505 | yes | yes | 12 | 44681 | epithelial; malignant |
| GSE253013 | included | 9 | 177205 | yes | yes | 12 | 23774 | epithelial; malignant |

### Primary patient-level Spearman (honest n / ρ / p)

**TACSTD2 vs B fraction**

| Cohort | n | ρ | p | note |
|---|---:|---:|---:|---|
| GSE131907 | 32 | +0.045 | 0.805 |  |
| GSE148071 | 41 | -0.002 | 0.992 |  |
| GSE154826 | 30 | -0.229 | 0.224 |  |
| GSE207422 | 15 | -0.636 | 0.0109 |  |
| GSE241934_IIT | 11 | -0.555 | 0.0767 |  |
| GSE241934_RWC | 24 | +0.510 | 0.0108 |  |
| GSE253013 | 9 | -0.083 | 0.831 |  |
| RE meta k=7 n=162 | 162 | -0.110 (95% CI -0.389 to +0.188) | 0.47 | τ²=0.102  |

**TACSTD2 vs TLS12**

| Cohort | n | ρ | p | note |
|---|---:|---:|---:|---|
| GSE131907 | 32 | +0.438 | 0.0121 |  |
| GSE148071 | 41 | -0.217 | 0.174 |  |
| GSE154826 | 30 | +0.334 | 0.0713 |  |
| GSE207422 | 15 | +0.279 | 0.315 |  |
| GSE241934_IIT | 11 | -0.291 | 0.385 |  |
| GSE241934_RWC | 24 | -0.022 | 0.92 |  |
| GSE253013 | 9 | -0.350 | 0.356 |  |
| RE meta k=7 n=162 | 162 | +0.078 (95% CI -0.175 to +0.322) | 0.549 | τ²=0.059  |

**CLDN4 vs B fraction**

| Cohort | n | ρ | p | note |
|---|---:|---:|---:|---|
| GSE131907 | 32 | -0.303 | 0.0921 |  |
| GSE148071 | 41 | +0.215 | 0.176 |  |
| GSE154826 | 30 | -0.252 | 0.179 |  |
| GSE207422 | 15 | -0.239 | 0.39 |  |
| GSE241934_IIT | 11 | -0.609 | 0.0467 |  |
| GSE241934_RWC | 24 | +0.244 | 0.25 |  |
| GSE253013 | 9 | +0.267 | 0.488 |  |
| RE meta k=7 n=162 | 162 | -0.092 (95% CI -0.329 to +0.157) | 0.472 | τ²=0.055  |

**CLDN4 vs TLS12**

| Cohort | n | ρ | p | note |
|---|---:|---:|---:|---|
| GSE131907 | 32 | -0.227 | 0.211 |  |
| GSE148071 | 41 | -0.292 | 0.0641 |  |
| GSE154826 | 30 | +0.397 | 0.0298 |  |
| GSE207422 | 15 | +0.054 | 0.85 |  |
| GSE241934_IIT | 11 | -0.173 | 0.612 |  |
| GSE241934_RWC | 24 | -0.104 | 0.628 |  |
| GSE253013 | 9 | -0.750 | 0.0199 |  |
| RE meta k=7 n=162 | 162 | -0.129 (95% CI -0.379 to +0.138) | 0.344 | τ²=0.072  |

### Sensitivity: partial Spearman residualizing epithelial fraction

Composition check. `frac_B` can fall when more epithelial cells are captured. Ranks of x and y are residualized on ranks of epithelial/gate fraction.

| Contrast | Cohort | n | ρ | p |
|---|---|---:|---:|---:|
| TACSTD2 vs B fraction given epi fraction | GSE131907 | 32 | -0.064 | 0.727 |
| CLDN4 vs TLS12 given epi fraction | GSE131907 | 32 | +0.036 | 0.845 |
| TACSTD2 vs B fraction given epi fraction | GSE148071 | 41 | +0.046 | 0.776 |
| CLDN4 vs TLS12 given epi fraction | GSE148071 | 41 | -0.168 | 0.292 |
| TACSTD2 vs B fraction given epi fraction | GSE154826 | 30 | -0.087 | 0.649 |
| CLDN4 vs TLS12 given epi fraction | GSE154826 | 30 | +0.341 | 0.0648 |
| TACSTD2 vs B fraction given epi fraction | GSE207422 | 15 | -0.637 | 0.0107 |
| CLDN4 vs TLS12 given epi fraction | GSE207422 | 15 | +0.067 | 0.813 |
| TACSTD2 vs B fraction given epi fraction | GSE241934_IIT | 11 | -0.326 | 0.328 |
| CLDN4 vs TLS12 given epi fraction | GSE241934_IIT | 11 | -0.199 | 0.557 |
| TACSTD2 vs B fraction given epi fraction | GSE241934_RWC | 24 | +0.572 | 0.00346 |
| CLDN4 vs TLS12 given epi fraction | GSE241934_RWC | 24 | -0.105 | 0.626 |
| TACSTD2 vs B fraction given epi fraction | GSE253013 | 9 | -0.068 | 0.863 |
| CLDN4 vs TLS12 given epi fraction | GSE253013 | 9 | -0.515 | 0.156 |

### Extra: TLS at single-cell

CXCL13+ / MS4A1+ detection by lineage. This is not a follicle or spatial TLS call.

See `tables/tls_single_cell.tsv`. CXCL13+ is T-enriched in mixed TME objects; MS4A1+ is B-restricted. GSE154826 TLS12 is scored inside a CD45-bead library (immune-internal), not a dissociated whole-tumor fraction.

| Cohort | lineage | n_cells | CXCL13 %pos | MS4A1 %pos |
|---|---|---:|---:|---:|
| GSE131907 | B | 16439 | 0.85 | 68.43 |
| GSE131907 | T | 43285 | 5.71 | 2.18 |
| GSE131907 | ALL | 128065 | 2.49 | 10.12 |
| GSE148071 | B | 662 | 1.21 | 88.82 |
| GSE148071 | T | 3247 | 23.71 | 3.33 |
| GSE148071 | ALL | 89887 | 2.03 | 1.34 |
| GSE154826 | B | 20989 | 0.91 | 82.52 |
| GSE154826 | T | 88049 | 7.81 | 4.89 |
| GSE154826 | gate | 16221 | 4.94 | 5.09 |
| GSE154826 | ALL | 193699 | 4.86 | 12.46 |
| GSE207422 | B | 7648 | 1.61 | 97.57 |
| GSE207422 | T | 32977 | 22.29 | 5.28 |
| GSE207422 | ALL | 92330 | 10.71 | 11.02 |
| GSE241934_IIT | B | 13056 | 1.52 | 97.37 |
| GSE241934_IIT | T | 38600 | 9.34 | 4.77 |
| GSE241934_IIT | ALL | 78691 | 5.49 | 20.08 |
| GSE241934_RWC | B | 44681 | 1.38 | 96.48 |
| GSE241934_RWC | T | 99004 | 16.72 | 4.26 |
| GSE241934_RWC | ALL | 229505 | 8.83 | 22.19 |
| GSE253013 | B | 23774 | 2.08 | 89.19 |
| GSE253013 | T | 65977 | 13.15 | 3.38 |
| GSE253013 | ALL | 177205 | 7.08 | 14.76 |

### Caveats

1. Dissociated 10x / Rhapsody is not histologic TLS.
2. GSE154826 is CD45-bead CITE-seq; epithelium sits in the author gate. Sort bias is recorded.
3. GSE253013 and GSE148071 malignant calls are marker proxies, not CopyKAT.
4. GSE241934 IIT and RWC are non-overlapping patients and are two strata of one paper.
5. Underpowered cohorts are labeled inconclusive. A non-significant ρ is not evidence of no association.
6. A6/B6 T-cell immune-cold is not re-tested.

### Reproduce

```bash
python3 methods/scrna_tls_meta/download.py
python3 methods/scrna_tls_meta/extract.py
python3 methods/scrna_tls_meta/analyze.py
```

