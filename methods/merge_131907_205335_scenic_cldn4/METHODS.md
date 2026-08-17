# Methods — merged GSE131907+GSE205335 CLDN4-only GRN proxy

ADDITIVE. **CLDN4 only.** Patient is the unit. A10 ELF3–CLDN4 is given.
No dual-high. No GSE207422. No invented ChIP. No pySCENIC cisTarget.

## Cohorts

| accession | malignant label | patient key |
|---|---|---|
| GSE131907 (Kim et al. 2020, PMID 32385277) | author `Cell_subtype` ∈ {Malignant cells, tS1, tS2, tS3} | GEO `patient_id` |
| GSE205335 (Hu / Ahn ICI biopsies) | author `lineage.sub` == Malignant cells | GEO SOFT `patient` |

GSE131907 tLung malignant cells are labeled tS1/tS2/tS3, not "Malignant cells".
They are included. PR #320's T/NK extract used `Malignant cells` only and
therefore dropped tLung; that is a T/NK-combo definition, not the GRN gate.

Processed UMI only. Skipped: GSE131907 2.86 GB log2TPM text, EGA FASTQ,
GSE207422.

## CLDN4 split

Within each patient, among author-malignant cells with UMI ≥ 200, split at the
median of CLDN4 log1p(CP10k). High = ≥ median. If the median is 0 (zero-inflated
CLDN4), high = CLDN4 > 0 vs low = 0 — otherwise every cell would be called high
and CLDN4-low patients would be dropped. CLDN4 is held out of every regulon and
program set. TACSTD2 is recorded and does not define groups.

Eligibility: ≥20 malignant cells. Primary paired test: ≥20 high **and** ≥20 low.
Sensitivity (≥10/10) is computed in the patient table but is not the claim.

## Regulons / programs

**Programs** (AUCell, not TF regulons): Hallmark IFN-α / IFN-γ (PR #267 A8
freeze), custom MHC-I antigen presentation, GO tight-junction organization ∪
bicellular TJ assembly, compact KRT + GO keratinization / keratinocyte
differentiation. Hallmark apical junction is a companion set.

**TF priors:** TRRUST v2 + DoRothEA + CollecTRI for STAT1/2, IRF1/7/9, NLRC5,
RFX5, GRHL1/2, OVOL1/2, KLF4/5, TP63, ELF3 (given), NKX2-1, SOX2 (controls).
A prior∩program AUCell is scored when the intersection is non-empty.

**AUCell:** Aibar linear recovery, `auc_threshold=0.05`, ranking inside the
extracted universe (program genes + TFs + prior targets). Not R AUCell binary.
Not motif enrichment.

**Pearson:** TF vs program genes in CLDN4-high cells, per cohort. High-specific
= r_high ≥ 0.15 and r_low < 0.10. Co-expression, not binding.

## Statistics

Patient-mean AUCell in high vs low: two-sided Wilcoxon signed-rank.
Cohorts are scored separately; patient deltas are stacked (k=2). Not Harmony.
Between-patient companion: Spearman of malignant CLDN4 %pos vs patient-mean
program AUCell, Fisher-z pooled. p-values are descriptive. Honest n = patients.

## Software

Python (numpy / pandas / scipy / matplotlib / statsmodels / rdata).
`pyscenic` is optional and unused.
