# EXTRA: immune-library TACSTD2/CLDN4 leak vs MPR

**Compartment: CD45+ / immune scRNA. Not malignant RNA.**
User A3 epithelial claim is taken as given.

## Extra table (patient-level detection fraction)

| Gene | Dataset | n MPR vs non-MPR | Median frac MPR | Median frac non-MPR | MWU p | Hedges’ g (MPR−non) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| TACSTD2 | GSE243013 immune atlas | 130 vs 112 | 0.00618 | 0.01085 | **9.76×10⁻⁴** | **−0.44** (−0.70, −0.19) |
| TACSTD2 | GSE229353 CD45+ beads | 2 vs 4 | 0.0130 | 0.0245 | 0.53 | −0.55 (−2.28, 1.18) |
| TACSTD2 | RE meta | 132 vs 116 | — | — | **5.85×10⁻⁴** | **−0.44** (−0.70, −0.19); I²=0 |
| CLDN4 | GSE243013 immune atlas | 130 vs 112 | 0.00204 | 0.00523 | **4.79×10⁻⁵** | **−0.52** (−0.78, −0.26) |
| CLDN4 | GSE229353 CD45+ beads | 2 vs 4 | 0.00956 | 0.0210 | 0.80 | −0.56 (−2.28, 1.17) |
| CLDN4 | RE meta | 132 vs 116 | — | — | **6.2×10⁻⁵** | **−0.52** (−0.77, −0.26); I²=0 |

Direction is **lower leak in MPR**. The meta is GSE243013 plus a tiny
same-direction second series. GSE229353 cannot move the pooled estimate
(2 complete MPR libraries after GEO dropped P03).

## RECIST (GSE243013 only)

Deposited `radiological_response` CR+PR (n=152) vs SD+PD (n=73):
TACSTD2 p=0.54; CLDN4 p=0.11. Not significant. Several labels are
unevaluable (including non-English strings) and were left out, not recoded.

## Must-try exclusions (honest)

- **GSE154826**: CD45+ / CITE-seq, 35 early-stage lesions. GEO family SOFT
  has no ICI / MPR / RECIST labels. The paper’s ICI claim is a separate bulk
  trial, not these scRNA patients.
- **Hui 2022 *Cell Death Dis***: CD45+ neoadjuvant, 12 patients. No GEO
  accession. Not public.
- **2024–2026 leftovers**: GSE280232 (T-cell GEX+TCR, not CD45-all);
  GSE303680 (PD1-IL2v T cells). Wrong library for epithelial leak.
- **GSE207422 / GSE241934 / GSE291670**: unsorted TME. A3 / other extras.

## Interpretation

This extra does **not** say immune cells biologically express TROP2 or
claudin-4 as an ICI-resistance program. Non-MPR tumors have more residual
viable epithelium by definition. Sparse TACSTD2/CLDN4 in a CD45+ matrix is
compatible with ambient RNA, doublet/misassignment, and residual tumor leak.
The large GSE243013 immune atlas is the only public series that can test
this at honest n.
