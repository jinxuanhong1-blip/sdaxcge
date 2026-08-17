# Methods — merged GSE131907 + GSE205335 CellChat-style CLDN4 → T/NK

ADDITIVE. **CLDN4 only.** TACSTD2 does not define groups. Patient is the unit.
GSE207422 is a different analysis and is not run. CellChat R and LIANA are not run.

## Datasets (public processed UMI only)

| Cohort | Citation | Matrix used | Not used |
| --- | --- | --- | --- |
| GSE131907 | Kim et al., *Nat Commun* 2020, PMID 32385277 | `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz` (~0.38 GB) + author annotation + series matrix | 2.86 GB `normalized_log2TPM_matrix.txt.gz`; EGA FASTQ `EGAD00001005054` |
| GSE205335 | Hu et al., palliative ICI biopsy/effusion scRNA | `GSE205335_Lung_IO_UMI_matrix.rds.gz` (dgCMatrix) + author identity + GEO SOFT | EGA raw `EGAD00001008703` |

No FASTQ. TPM text is skipped because UMI exists.

## Compartments (author labels)

| Cohort | Malignant | T/NK |
| --- | --- | --- |
| GSE131907 | `Cell_type == Epithelial cells` and `Cell_subtype ∈ {Malignant cells, tS1, tS2, tS3}` on tumor origins `tLung, tL/B, mLN, mBrain` | `Cell_type ∈ {T lymphocytes, NK cells}` on the same tumor origins |
| GSE205335 | `lineage.sub == Malignant cells` | `lineage.total == T/NK cells` |

Honest label limits:

- Primary tLung epithelium is annotated **tS1 / tS2 / tS3**, not `Malignant cells`. Those tS* barcodes are included in the primary paired test.
- Author `Malignant cells` in GSE131907 sit in tL/B, mLN, mBrain.
- PE unlabeled epithelium, nLung, and nLN are excluded.
- CopyKAT was not re-run.
- GSE205335 MPR/NMPR is unlabeled; RECIST is not used as MPR.

GEO `patient_id` (GSE131907) or GEO `patient` (GSE205335) is the unit. GSE131907 tumor-origin samples from the same patient are pooled. Same-patient T/NK means T/NK from those same tumor origins, not a normal-lung match.

GSE131907 gene symbols still use pre-HGNC nectin names (`PVRL2` = NECTIN2). Those rows are aliased to CellChatDB v2 symbols before scoring. GSE205335 already uses `NECTIN2`.

## CLDN4 splits (malignant cells, within patient)

Score = `log1p(CP10k)` CLDN4 on that patient's malignant cells.

| Split | High | Low |
| --- | --- | --- |
| Median | ≥ median (if median is 0: CLDN4>0 vs =0) | < median (or =0) |
| Q4 vs Q1 | `pd.qcut` on average ranks, Q4 tail | Q1 tail |

**Honest paired n.** A patient is scored only if:

- malignant CLDN4-high bin ≥ 10 cells
- malignant CLDN4-low bin ≥ 10 cells
- same-patient T/NK ≥ 20 cells

Patients missing either bin or T/NK are out. Cells are not n.

## CellChat-style probability

Jin et al. 2021 on **CellChatDB v2** protein pairs (Secreted / Cell–Cell Contact / ECM–Receptor; non-protein dropped). A pair is kept only if every ligand and receptor subunit is present in **both** UMI matrices.

1. Per cell group: 10% truncated mean of `log1p(CP10k)`; complexes = geometric mean of subunits (0 if any subunit mean is 0).
2. Detected if each complex has expressing-cell fraction ≥ 0.10 (AND rule).
3. \(P = (L \cdot R) / (K_h + L \cdot R)\) with \(K_h = 0.5\).

Outgoing = Mal_high or Mal_low → that patient's T/NK. Incoming is not the primary table.

## Tests

**Primary.** Wilcoxon signed-rank on per-patient \(P_{\mathrm{high}}\) vs \(P_{\mathrm{low}}\) for pairs detected on **both** arms. Minimum paired n after the detect gate = 6. Thin n (<8) is flagged, not hidden. p-values are descriptive.

**Companion (PR #320 slice).** Locked extract: GSE131907 author-`Malignant cells` samples with n_mal≥20 (n=21, sample-level in that extract) + GSE205335 patients with ≥20 malignant and ≥20 T/NK (n=22). Within-cohort median and Q4 vs Q1 of patient/sample malignant CLDN4 %pos. Score all-malignant → same-unit T/NK. Mann–Whitney on per-unit *P* (detected ≥3 per arm). This slice's T/NK association is **given** (n=23 Q4 vs Q1 r=−0.705 p=0.0003; continuous n=43 ρ=−0.479) and is not re-ranked.

Focused pairs are always written if subunits exist: classical MHC-I (HLA-A/B/C–CD8), T-recruit (CXCL9/10/16, CCL4/5), CD274–PDCD1, NECTIN2–TIGIT.

## Software

Python (numpy / pandas / scipy / matplotlib / rdata). CellChatDB v2 parsed from public `CellChatDB.human.rda` (`jinworks/CellChat`). No R CellChat runtime.
