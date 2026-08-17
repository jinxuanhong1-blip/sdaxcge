# Methods — CellChat-style malignant CLDN4-high vs T/NK (GSE131907)

Public data only. **GSE207422 CellChat is a different analysis** (`methods/scrna_cellchat` / `methods/scrna_cellchat_cldn4`) and is not re-run here. CellChat R and LIANA are not installed.

## Dataset

GEO `GSE131907` (Kim et al., *Nat Commun* 2020, PMID 32385277). Treatment-naive LUAD atlas. Processed files used:

- `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz` (~0.38 GB)
- `GSE131907_Lung_Cancer_cell_annotation.txt.gz`
- `GSE131907_series_matrix.txt.gz` (patient / stage / origin)

**Not used:** 2.86 GB `normalized_log2TPM_matrix.txt.gz`; EGA raw FASTQ `EGAD00001005054`; GSE207422.

There are **no ICI / MPR / RECIST labels** on GEO. This is not an A3 NMPR>MPR test.

## Compartments (author labels)

| Compartment | Rule |
| --- | --- |
| Malignant | `Cell_type == Epithelial cells` and `Cell_subtype ∈ {Malignant cells, tS1, tS2, tS3}` |
| T/NK | `Cell_type ∈ {T lymphocytes, NK cells}`, merged as `TNK` |

Honest label limits:

- Primary tLung tumor epithelium is annotated **tS1 / tS2 / tS3**, not `Malignant cells`.
- Author `Malignant cells` sit in **tL/B, mLN, mBrain**.
- PE epithelial barcodes are mostly **unlabeled** and are **excluded** from the malignant set.
- CopyKAT was not re-run. Normal lung / normal LN are excluded from the kept split.

## CLDN4 split

Among malignant cells in the sample set, `log1p(CP10k)` CLDN4 is split by **tertiles** (high = top third, low = bottom third; middle dropped). TACSTD2 is not used to define groups.

Kept split: tumor origins **tLung + tL/B + mLN + mBrain**. Sensitivities: tLung only; mets (tL/B + mLN + mBrain) only.

## CellChat-style probability

Communication probability follows Jin et al. 2021 on **CellChatDB v2** protein pairs (Secreted / Cell–Cell Contact / ECM–Receptor; non-protein pairs dropped):

1. Keep a pair only if every ligand and receptor subunit is present in the UMI matrix.
2. Per cell group: 10% truncated mean of `log1p(CP10k)` for each gene; complexes = geometric mean of subunits (0 if any subunit mean is 0).
3. Detected if each complex has expressing-cell fraction ≥ 0.10 (AND rule).
4. \(P = (L \cdot R) / (K_h + L \cdot R)\) with \(K_h = 0.5\).
5. Groups with fewer than 25 cells are not scored.

## Significance

Null for the CLDN4 contrast: permute high/low labels among malignant cells (T/NK labels fixed), 100 permutations, \(p = (\#\{P_{\mathrm{perm}} \ge P_{\mathrm{obs}}\} + 1) / 101\). A pair is significant if detected, \(P > 0\), and \(p < 0.05\). The smallest possible p is 1/101 = 0.0099.

Outgoing = Mal → TNK (the ligand table). Incoming = TNK → Mal (supplement). Differential pairs are those significant on at least one arm with \(P_{\mathrm{high}} \ne P_{\mathrm{low}}\).

Means are **cell-pooled**. That is not a patient-level mixed model. Dominant samples are written to `results/dominance_*.tsv`.

## Software

Python (numpy / pandas / scipy / matplotlib). CellChatDB v2 parsed from the public `CellChatDB.human.rda` (`jinworks/CellChat`). Raw FASTQ was not used.
