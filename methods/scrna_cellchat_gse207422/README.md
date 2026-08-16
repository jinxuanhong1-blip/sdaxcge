# methods/scrna_cellchat_gse207422

Additive ligand–receptor analysis on the **public** GSE207422 processed UMI matrix
(Hu et al., *Genome Med* 2023; PMID 36869384). User A3 lineage and malignant-like
rules are taken as given.

Question: do **TACSTD2-high** malignant-like cells show weaker outgoing
T-recruiting / IFN / MHC-I communication toward T/NK, and weaker incoming
signals from T/NK, than TACSTD2-low malignant-like cells? The same split is
scored for **CLDN4**.

## What this is / is not

- **Is:** CellPhoneDB-pair scoring plus LIANA’s CellPhoneDB method when the
  Python package imports. Honest pathway ranks with n cells and n patients.
- **Is not:** CellChat output. R is not used here; no CellChat probability
  tables are written or renamed.

## Reproduce

```bash
cd methods/scrna_cellchat_gse207422
python3 -m pip install -r requirements.txt
python3 scripts/00_download.py
python3 scripts/01_build_pairs.py
python3 scripts/02_run_ccc.py
```

Outputs land in `results/` (`RESULTS.md`, ranks, LIANA CSVs if run, figures).

## Methods in one paragraph

Malignant-like = A3 epithelial gate AND zero UMI of SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3.
High/low = global median of log1p(CP10k) TACSTD2 or CLDN4 among malignant-like cells.
Pair list = CellPhoneDB v5 `interactors` plus a small curated overlay for the
three axes. Score = mean of partner means on log1p(CP10k), complexes = min
subunit (Efremova 2020; Garcia-Alonso 2022). Inferential unit = patient
(paired Wilcoxon). LIANA is secondary and downsampled; its p-values are not
patient-level.
