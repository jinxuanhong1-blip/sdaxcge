# METHODS — markermalig15_cldn4_hiend

Additive **CLDN4-only** CellChat-style + NicheNet-style analysis on the marker-malignant combo that already differs.

## Combo (given, not re-audited)

PR #290 / #279 highlighted cut:

- Marker-malignant **%pos** vs T/NK
- Cohorts: **GSE253013** (9 tumor patients) + **GSE291670** (6 tumors)
- **N=15 · CLDN4 ρ=−0.714 · p=0.0217 · I²=23%**

Sample list: `data/combo/samples_n15.tsv`. Adjacent-normal GSE253013 rows are not in the combo. n is not inflated.

## Sender / receiver

- **Marker-malignant** matches the combo notes:
  - GSE253013: marker epithelium with normal-lung score ≤ 0.05
  - GSE291670: marker epithelium below the 75th percentile of the normal-lung program
- **CLDN4-high / low**: median split of malignant `log1p(CP10k)` CLDN4. **TACSTD2 is not used.** Dual-high is scored only as a companion count.
- **T/NK**: marker-argmax T + NK.

## CellChat-style

CellChatDB v2 protein pairs. Probability is the published Hill / mass-action form on 10% truncated means, with `expr_prop ≥ 0.10`. High/low labels among malignant cells are permuted 100 times. The CellChat R package is not run. Smallest p = 1/101 = 0.0099.

Cell-pooled truncated means are **not** a 15-patient mixed model. n=15 is reported as the sample unit.

## NicheNet-style

NicheNet-v2 ligand–target prior (Browaeys et al., Zenodo 7074291). Ligand activity is Pearson / AUROC / AUPR of the prior column against an **a priori** T/NK IFN + cytotoxicity gene set. Ligands must be expressed in ≥10% of CLDN4-high marker-malignant cells. Direction is not assumed from the unsigned prior.

## Honest n

n=15 is small. MRC004 and MRC007 have 16 and 18 marker-malignant cells in the combo notes. Those samples stay in the n=15 list; they are not dropped to dress the p-value, and they are not replaced.
