# GSE207422 demo readout (computed, not fabricated)

Hu et al., *Genome Med* 2023 (PMID 36869384). 92,330 cells × 1,451 streamed
genes, 15 patients, marker-gated lineages (GEO deposited no cell-type column).
Epithelium split by TACSTD2+CLDN4 score, top quartile = high
(threshold 1.285; 3,052 high / 9,153 low). 13 patients had ≥20 cells in
**both** states and enter the paired test.

## Primary — patient-level chemokine pseudobulk

`gse207422_chemokine_pseudobulk.csv` — paired Wilcoxon on log1p(CPM),
FDR-BH across 10 ligands.

| Chemokine | n paired | Δ (high − low) | padj |
|-----------|----------|----------------|------|
| CXCL10 | 13 | −0.64 | 0.54 |
| CCL3 | 13 | −0.52 | 0.54 |
| CXCL13 | 13 | −0.46 | 0.54 |
| CXCL9 | 13 | −0.36 | 0.76 |
| CCL21 | 13 | −0.28 | 0.72 |
| CCL4 | 13 | −0.24 | 0.54 |
| CCL19 | 13 | −0.24 | 0.54 |
| CXCL11 | 13 | −0.05 | 1.00 |
| CCL5 | 13 | +0.14 | 0.54 |
| CXCL16 | 13 | +0.20 | 0.54 |

**8/10 ligands are lower in the TACSTD2/CLDN4-high state; 0/10 reach FDR < 0.05.**
The CXCR3-axis trio (CXCL9/10/11) all point the hypothesized direction.
CXCL16 (and CCL5) point the other way.

## LIANA focus (hypothesis only)

`gse207422_liana_focus.csv` / `gse207422_liana_high_vs_low.csv`.

After `expr_prop=0.10` and restricting to CXCR3/CCR5/CXCR6/CCR1, the only
surviving focus edge is **CXCL16–CXCR6** (to T and NK). Its
`magnitude_rank` is *better* (stronger) from the **high** sender — consistent
with the pseudobulk (CXCL16 slightly higher in high), **not** with a
"reduced recruitment" story. CXCL9/10/11–CXCR3 and CCL5–CCR5 did not clear
the expression-fraction filter in epithelium.

## What we will and will not write

- **Will:** in this 15-patient BD-Rhapsody ICI cohort, TACSTD2/CLDN4-high
  epithelium shows a *non-significant* downward shift in CXCR3-axis ligands
  (8/10 ligands lower; min padj = 0.54). LIANA does not recover a weaker
  CXCR3 edge.
- **Will not:** "TACSTD2/CLDN4-high epithelium reduces T-cell recruitment."
  The primary test is not significant; LIANA's only focus edge goes the
  other way; there is no spatial / IHC triangulation in this demo.
- This is exactly the playbook's `07` joint rule: primary does **not**
  support → do not claim recruitment from LIANA/SCENIC/NicheNet.

## Caveats (copy into any write-up)

- BD Rhapsody, not 10x. Marker gates ≠ published annotation.
- Ambient RNA not re-estimable from the processed matrix.
- n=13 paired; MPR n=4 — do not split on response here.
- Chemokine dropout is high; fraction-expressing should be inspected before
  any follow-up claim.
