# Methods — CellChat-style TACSTD2-high vs low vs T/NK (GSE207422)

Public data only. Author Seurat / CopyKAT objects are not deposited; A3/B6 cell-state framing is taken as given and is not re-tested here.

## Dataset

GEO `GSE207422` (Hu et al., *Genome Medicine* 2023, PMID 36869384). One processed UMI matrix (`GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz`; 24,292 genes × 92,330 barcodes) plus the 15-row sample sheet. Stage IIIA NSCLC; 3 pre-treatment biopsies and 12 post-treatment resections after neoadjuvant PD-1 + chemotherapy. Pathologic response on the sheet: MPR (including pCR) n=4, NMPR n=8 among the 12 resections.

GSE241934 and GSE253013 were not used: GSE207422 already has malignant epithelium, T/NK, and MPR labels in a single matrix.

## Lineage and compartments

Canonical marker scores on `log1p(CP10k)` (same panel as the A3 GSE207422 reconstruction). Assigned lineage = argmax; T vs NK broken by CD3E when scores are close. **Malignant compartment** = epithelial lineage (CopyKAT barcode IDs are not public). **T/NK** = T or NK, merged as `TNK` for ligand–receptor tests.

Primary splits use the 12 post-treatment samples only.

## TACSTD2 split

Among epithelial cells in the sample set, `log1p(CP10k)` TACSTD2 is split by **tertiles** (high = top third, low = bottom third; middle dropped) or by the **median**. Combinatorial split: the same tertiles crossed with MPR vs NMPR (communication scored within response arm).

## CellChat-style probability

CellChat R was not installed. Communication probability follows Jin et al. 2021 on **CellChatDB v2** protein pairs (Secreted / Cell–Cell Contact / ECM–Receptor; non-protein pairs dropped):

1. Keep a pair only if every ligand and receptor subunit is present in the UMI matrix.
2. Per cell group: 10% truncated mean of `log1p(CP10k)` for each gene; complexes = geometric mean of subunits (0 if any subunit mean is 0).
3. Detected if each complex has expressing-cell fraction ≥ 0.10 (AND rule).
4. \(P = (L \cdot R) / (K_h + L \cdot R)\) with \(K_h = 0.5\).
5. Groups with fewer than 25 cells are not scored.

## Significance

Null for the TACSTD2 contrast: permute high/low labels among malignant cells (T/NK labels fixed), 100 permutations, \(p = (\#\{P_{\mathrm{perm}} \ge P_{\mathrm{obs}}\} + 1) / 101\). A pair is significant if detected, \(P > 0\), and \(p < 0.05\).

Outgoing = Mal → TNK. Incoming = TNK → Mal. Differential pairs are those significant on at least one arm with \(P_{\mathrm{high}} \ne P_{\mathrm{low}}\).

No network is drawn for a non-significant pair.

## Which split is kept

Each split gets an immune-cold / barrier score from significant differential pairs: plus one for barrier or inhibitory ligands higher in TACSTD2-high (outgoing), recruiting chemokines/cytokines lower in high (outgoing), and attack ligands (IFNG, TNF, FASLG, TRAIL, LTA) lower into high (incoming); minus the opposite. The kept split is the highest score, then more contributing patients, then the simple post-treatment tertile.

## Software

Python (numpy / pandas / scipy). CellChatDB v2 parsed from the public `CellChatDB.human.rda` (`jinworks/CellChat`). Raw FASTQ was not used.
