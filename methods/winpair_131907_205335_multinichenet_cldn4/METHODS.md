# Methods — MultiNicheNet / multi-sample NicheNet on winning pair GSE131907+GSE205335

ADDITIVE. **CLDN4 only.** TACSTD2 is not a gate. Dual-high is not run.
**GSE207422 is not used** (NicheNet PR #334 NS). Patient is the unit.

## Cohorts

| Series | Role | Public input | Skipped |
| --- | --- | --- | --- |
| [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907) | Treatment-naive LUAD atlas (Kim et al., *Nat Commun* 2020, PMID 32385277) | Author cell annotation + processed raw UMI (~0.38 GB gzip) + series matrix | 2.86 GB log2TPM; EGA `EGAD00001005054` |
| [GSE205335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335) | Palliative ICI biopsy / effusion scRNA (Hu / Ahn / Lee) | Author identity + processed UMI `dgCMatrix` + GEO SOFT | EGA `EGAD00001008703` |
| GSE207422 | **Not used** | — | PR #334 already NS; do not redo; do not add |

GSE131907 has **no ICI / MPR**. GSE205335 RECIST is recorded and is **not** substituted for MPR.

## Labels (author, not recomputed)

**GSE131907**

- Tumor-origin: `Sample_Origin` ∈ {tLung, tL/B, mLN, mBrain, PE}.
- Malignant: tumor-origin **and** `Cell_subtype` ∈ {Malignant cells, tS1, tS2, tS3}.
- T/NK: tumor-origin **and** `Cell_type` ∈ {T lymphocytes, NK cells}.
- PE epithelium is unlabeled and is **not** malignant. nLung AT1/AT2/Club/Ciliated are not malignant.

**GSE205335**

- Patient ID is parsed from author `orig.ident` (`EBUS-06-3P` → `P0006`; multi-site libraries with the same number are one patient). The GEO family SOFT on FTP is series-level only (no per-GSM characteristics).
- Normal libraries (`LUNG-N*`, `NS-*`) are dropped. Remaining libraries are treated as tumor/met.
- Malignant: tumor library **and** `lineage.sub == "Malignant cells"`.
- T/NK: tumor library **and** `lineage.total == "T/NK cells"`.

Tumor samples are **pooled per patient**. Cells are not n.

## CLDN4 split (pre-registered, CLDN4 only)

Among malignant cells **within each dataset**, CLDN4 is `log1p(CP10k)` using the cell’s full UMI total. **High** = at or above that dataset’s malignant median; **low** = below. TACSTD2 is not used.

A patient enters the paired sender contrast if it has ≥10 CLDN4-high malignant, ≥10 CLDN4-low malignant, and ≥20 T/NK cells.

## Multi-sample NicheNet (documented prior; no R)

`nichenetr` and `multinichenetr` are not installed. The published NicheNet-v2 human ligand–target matrix and LR network (Zenodo [10.5281/zenodo.7074291](https://doi.org/10.5281/zenodo.7074291); Browaeys et al.) are converted with Python `rdata`.

1. **Potential ligand.** Ligand detected in ≥10% of CLDN4-high malignant cells (cell-weighted across eligible patients), ≥1 receptor detected in ≥10% of T/NK cells, and a column in the v2 matrix.
2. **Ligand activity.** Pearson, AUROC, and AUPR of the ligand’s prior target scores versus gene-set membership on background genes (T/NK-expressed panel ∩ prior rows). Rank by Pearson (NicheNet default). This is unsigned regulatory potential.
3. **Sender DE (MultiNicheNet term).** Per patient, mean ligand `log1p(CP10k)` in CLDN4-high vs CLDN4-low malignant cells. Test = paired Wilcoxon across patients (combined and per dataset). FDR-BH within dataset scope.
4. **Receiver association.** Per-patient T/NK gene means vs malignant CLDN4 (Spearman). Empirical up/down sets use p < 0.15 (honest; not FDR<0.05).
5. **Prioritization.** For each a priori set, min–max scale (i) Pearson activity, (ii) paired high−low ligand Δ, (iii) fraction of patients with a T/NK receptor at ≥10%, then average. This is the MultiNicheNet-like rank written to `ligand_activity_table.tsv`.

Primary gene sets: a priori IFN and cytotoxicity (`scripts/gene_sets.py`). Exhaustion is sensitivity only.

## Software

Python 3; pandas / numpy / scipy / scikit-learn / matplotlib / rdata. No R. Scripts: `00_download.py` … `03_analyze.py`. Matrices and RDS priors are not stored in git.
