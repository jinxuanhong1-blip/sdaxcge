# Methods

1. **Cohort.** Public GSE207422 processed UMI matrix and sample sheet (Hu et al., *Genome Med* 2023, PMID 36869384). Post-treatment surgery samples only. Pathologic response from the GEO sheet: pCR counted as MPR (P06). GSE253013 was listed and not downloaded (9.3 GB RDS; no MPR/NMPR).

2. **Cell types.** DRMref public annotations (`GSE207422_Tor` + `GSE207422_Sin`; Liu et al., *NAR* 2024). 30,877 barcodes match GEO. Senders = `Malignant cells`. Receivers = `CD8+ T cells` + `CD4+ T cells` + `NK cells`. These are marker-based labels, not the undeposited Hu CopyKAT calls.

3. **Normalization.** Gene UMIs from the GEO matrix. Size factor = DRMref `nCount_RNA`. Score = `log1p(1e4 * UMI / nCount)`.

4. **High-barrier malignant (pre-registered).** Among malignant cells, TACSTD2 and CLDN4 `log1p(CP10k)` both ≥ the malignant-cell median.

5. **Patient scores.** Per patient: mean malignant TACSTD2, CLDN4, and their average; fraction high-barrier; mean T/NK cytotoxicity and exhaustion (a priori lists in `scripts/gene_sets.py`); T/NK fraction. Wilcoxon NMPR vs MPR and Spearman vs barrier use **n=12 patients**. Cells are never the independent unit for those tests.

6. **Expressed genes.** Detected in ≥ 10% of cells in the named compartment (NicheNet 10x-style rule). Potential ligand = ligand expressed in the sender set, ≥1 receptor expressed in the receiver set, and a column in the v2 ligand–target matrix.

7. **Prior.** NicheNet-v2 human `ligand_target_matrix_nsga2r_final.rds` and `lr_network_human_21122021.rds` (Zenodo 10.5281/zenodo.7074291; Browaeys et al.). Converted with Python `rdata` (no R / `nichenetr`). Ligand activity = Pearson, AUROC, and AUPR of the ligand’s prior target scores versus gene-set membership on background genes. Rank by Pearson.

8. **Background.** Intersection of T/NK-expressed extracted panel genes with prior target rows (256 genes). This is narrower than a full-transcriptome NicheNet run and is reported as such.

9. **Gene sets.** A priori cytotoxicity and exhaustion lists (playbook). Empirical NMPR-up / NMPR-down: patient-mean T/NK `log1p(CP10k)`, Wilcoxon p < 0.15, split by sign. That cutoff is honest about n=4 vs 8; it is not an FDR<0.05 call.

10. **Combinatorial.** Repeat potential-ligand construction and activity ranking with NMPR-only high-barrier senders → NMPR T/NK, and MPR-only high-barrier senders → MPR T/NK. Shared ligands receive the same activity number (prior × set unchanged). Report set difference and patient-level ligand means (Wilcoxon, unadjusted).

11. **Software.** Python 3.12; pandas / numpy / scipy / scikit-learn / matplotlib. Scripts: `scripts/00_download.py` … `03_analyze.py`.
