# Harvested GSE148071 patient scores (no new scRNA download)

Wu *et al.* *Nat Commun* 2021 (GSE148071): 42 advanced NSCLC tumors, one
sample per patient. These copies are already-computed scores from this
repository. The 10x matrices were not re-downloaded.

| Path | Source | Used for |
|---|---|---|
| `tisch_GSE148071_units.tsv` | PR TISCH NSCLC pool (`methods/tisch_nsclc_pool`) | Malignant/epithelial CLDN4 vs T/NK fraction |
| `tisch_GSE148071_given_spearman.tsv` | same | Locked continuous ρ (taken as given) |
| `tls_GSE148071_patients.tsv` | PR #274 / #320 TLS extract | Malignant-compartment CLDN4 vs CXCL13+ among T |

TISCH eligible unit: ≥20 scored epithelial/malignant cells and ≥20 T/NK.
TLS CXCL13+ row: `compartment=malignant` and finite `frac_CXCL13pos_T`.
Do not write n=42 for either test.
