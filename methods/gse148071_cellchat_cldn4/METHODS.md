# Methods — CellChat-style CLDN4-high vs low vs T/NK (GSE148071)

Additive folder. Prior TACSTD2 / GSE207422 CellChat analyses are taken as given and are not re-run. **This folder splits putative malignant cells on CLDN4 only.**

Public data only. Author Seurat / inferCNV objects are not deposited on GEO.

## Dataset

GEO `GSE148071` (Wu et al., *Nature Communications* 2021, PMID 33953163). Forty-two diagnostic biopsies from stage III/IV NSCLC (Singleron GEXSCOPE; HiSeq X Ten). The GEO supplementary archive is `GSE148071_RAW.tar` (42 `GSM*_P*_exp.txt.gz` raw-count matrices). The series matrix has age and sex; it does **not** have histology, ICI response, or cell-type labels.

GSE207422, GSE131907, and GSE253013 were not used here.

## Lineage and compartments

No published barcode labels are in the GEO archive. Canonical marker scores on `log1p(CP10k)` (same panel as the GSE207422 CellChat folders). Assigned lineage = argmax; T vs NK broken by CD3E when scores are close. **Malignant compartment** = epithelial lineage (**putative**; CopyKAT IDs are not public). **T/NK** = T or NK, merged as `TNK` for ligand–receptor tests.

A patient is **eligible** for the communication contrast if it has ≥25 epithelial cells **and** ≥25 T/NK cells. Patients below that floor are reported in `results/per_sample.tsv` and are not scored.

## CLDN4 split

Among epithelial cells from eligible patients, `log1p(CP10k)` **CLDN4** is split by **tertiles** (high = top third, low = bottom third; middle dropped) or by the **median**. TACSTD2 is recorded but is **not** used to define groups. There is no MPR/NMPR combinatorial split: this series has no pathologic-response labels.

## CellChat-style probability

CellChat R was not installed. Communication probability follows Jin et al. 2021 on **CellChatDB v2** protein pairs (Secreted / Cell–Cell Contact / ECM–Receptor; non-protein pairs dropped):

1. Keep a pair only if every ligand and receptor subunit is present in the first patient’s gene list (shared Singleron universe).
2. Per cell group: 10% truncated mean of `log1p(CP10k)` for each gene; complexes = geometric mean of subunits (0 if any subunit mean is 0).
3. Detected if each complex has expressing-cell fraction ≥ 0.10 (AND rule).
4. \(P = (L \cdot R) / (K_h + L \cdot R)\) with \(K_h = 0.5\).
5. Groups with fewer than 25 cells are not scored.

## Significance

Null for the CLDN4 contrast: permute high/low labels among malignant cells (T/NK labels fixed), 100 permutations, \(p = (\#\{P_{\mathrm{perm}} \ge P_{\mathrm{obs}}\} + 1) / 101\). A pair is significant if detected, \(P > 0\), and \(p < 0.05\).

Outgoing = Mal → TNK. Incoming = TNK → Mal. Differential pairs are those significant on at least one arm with \(P_{\mathrm{high}} \ne P_{\mathrm{low}}\).

No network is drawn for a non-significant pair.

## Which split is kept

Tertile vs median: keep the single contrast with the higher immune-cold / barrier score (barrier or inhibitory ligands higher in CLDN4-high outgoing; recruiting chemokines/cytokines lower in high outgoing; attack ligands lower into high incoming; minus the opposite). Ties go to fewer recruit-up pairs, then more malignant cells.

## Honest n

Cell-pooled truncated means are not a patient-level mixed model. Per-sample epithelial and T/NK counts are written to `results/per_sample.tsv` and plotted (`fig_n_per_sample.png`). One biopsy can dominate the cell-pooled means; that is reported, not hidden. Eligible n is the number of patients meeting the 25/25 floor, not 42.

## Software

Python (numpy / pandas). CellChatDB v2 parsed from the public `CellChatDB.human.rda` (`jinworks/CellChat`). Raw FASTQ was not used. CellChat R and LIANA were not run.
