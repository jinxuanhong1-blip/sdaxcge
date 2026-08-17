# Winning-pair CLDN4-only patient-pseudobulk DE (GSE131907 + GSE205335)

ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4. Not a bigger merge.
Winning pair is taken from PR #320 (author malignant CLDN4 %pos vs T/NK:
Spearman ρ=−0.479, N=43; Q4 vs Q1 r=−0.705, n=23 compared).

**Method: patient-pseudobulk OLS on log2(TMM-CPM+1).** Not muscat.
muscat `pbDS` is the same idea (sum cells → patient × cell-type count matrix →
bulk DE). Here the patient (GSE131907: sample) is already the unit, so a
direct limma-style OLS is the honest tool. No R/limma/edgeR/muscat in this
environment; TMM + OLS + BH is the implementation. p-values are descriptive.

## Design (locked)

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | GSE131907 + GSE205335 only | Not GSE207422 / 291670 / 253013 / 148071 |
| Split | Within-cohort author malignant CLDN4 **%pos Q4 vs Q1** (PR #320) | Quartile tails, not median; mid quartiles unused in binary DE |
| Continuous | CLDN4 %pos z-scored, all gated units, cohort covariate | Linear; not a causal model |
| Malignant | Author malignant UMI-sum (existing GSEA matrices; n_mal≥30) | P4001 (27 malignant cells) drops from malignant DE |
| T/NK | Same patients; author T/NK UMI-sum; ≥20 T/NK cells | GSE131907 is **sample-level** |
| Model | ~ cohort + CLDN4_Q4 (combined); ~ CLDN4_Q4 (single) | Small n; no muscat mixed model; no voom weights |
| Families | IFN (Hallmark IFNα/γ), MHC-I custom + MHC-II, TJ (KEGG/GO + focal), chemokine panel | Chemokine panel is compact, not MSigDB C2 |

## Honest n

Winning-pair labels (PR #320 author %pos):

- GSE131907: n=21 samples (tumor origins, n_mal≥20). Q1=6 Q4=5 (sample-level).
- GSE205335: n=22 patients. Q1=6 Q4=6.

| contrast | compartment | n_low / n_high | n_genes | note |
|---|---|---:|---:|---|
| Q4 vs Q1 combined | malignant | 11/11 | 15883 | cohort covariate; P4001 out of malignant matrix |
| Q4 vs Q1 combined | T/NK | 12/11 | 12020 | same Q labels; T/NK from those patients |
| continuous combined | malignant | n=42 | 16825 | CLDN4 %pos z |
| continuous combined | T/NK | n=43 | 12735 | CLDN4 %pos z |

Per-cohort Q4 vs Q1 n is in `tables/n_honest.tsv`. Do not quote a pooled n that
ignores the sample-vs-patient unit difference on GSE131907.

## IFN / MHC / TJ / chemokine (combined Q4 vs Q1)

Positive logFC = higher in CLDN4-high (Q4) than CLDN4-low (Q1).

### Malignant cells

| family | n_tested | n p<0.05 (up/down) | n FDR<0.05 | median logFC | top gene (logFC, p, FDR) |
|---|---:|---|---:|---:|---|
| IFN | 214 | 42 (2/40) | 0 | -0.530 | CCL5 (-2.889, 5.02e-05, 0.09751) |
| MHC | 29 | 7 (0/7) | 0 | -0.922 | TAP2 (-1.903, 0.00242, 0.1731) |
| TJ | 191 | 24 (19/5) | 0 | +0.065 | PARD6A (+1.672, 0.001065, 0.1382) |
| chemokine | 21 | 8 (0/8) | 0 | -0.823 | CCL5 (-2.889, 5.02e-05, 0.09751) |

### T/NK cells (same patients)

| family | n_tested | n p<0.05 (up/down) | n FDR<0.05 | median logFC | top gene (logFC, p, FDR) |
|---|---:|---|---:|---:|---|
| IFN | 206 | 39 (0/39) | 0 | -0.452 | JAK2 (-1.348, 0.000934, 0.3402) |
| MHC | 29 | 4 (0/4) | 0 | -0.453 | PSMB10 (-1.099, 0.004764, 0.4229) |
| TJ | 149 | 35 (20/15) | 0 | -0.050 | CLDN4 (+3.487, 0.0003, 0.242) |
| chemokine | 20 | 3 (0/3) | 0 | -0.531 | CCL2 (-1.412, 0.01425, 0.4534) |

T/NK top TJ hits are **epithelial genes** (CLDN4, CDH1, CLDN3, TJP1). That is
ambient / doublet leakage from malignant cells, not a T/NK tight-junction
program. The T/NK IFN arm is one-sided down (39/39 p<0.05 genes down; JAK2
top). Do not quote T/NK CLDN4 logFC as a lymphocyte finding.

Headline genes (combined Q4 vs Q1, family members only, lowest p):

| compartment | family | gene | n_Q1 | n_Q4 | logFC | p | FDR |
|---|---|---|---:|---:|---:|---:|---:|
| malignant | IFN|chemokine | CCL5 | 11 | 11 | -2.889 | 5.02e-05 | 0.09751 |
| malignant | IFN | GZMA | 11 | 11 | -2.524 | 0.0001133 | 0.09751 |
| malignant | IFN | GBP4 | 11 | 11 | -2.760 | 0.0003456 | 0.1038 |
| malignant | IFN | CD69 | 11 | 11 | -2.738 | 0.000444 | 0.1038 |
| malignant | chemokine | CCL4 | 11 | 11 | -2.334 | 0.0007034 | 0.1214 |
| malignant | chemokine | CCL3 | 11 | 11 | -2.333 | 0.0007827 | 0.1259 |
| malignant | TJ | PARD6A | 11 | 11 | +1.672 | 0.001065 | 0.1382 |
| malignant | IFN | GBP2 | 11 | 11 | -2.678 | 0.001173 | 0.1396 |
| T/NK | TJ | CLDN4 | 12 | 11 | +3.487 | 0.0003 | 0.242 |
| T/NK | TJ | CDH1 | 12 | 11 | +2.300 | 0.0003367 | 0.242 |
| T/NK | TJ | CLDN3 | 12 | 11 | +3.253 | 0.00034 | 0.242 |
| T/NK | IFN | JAK2 | 12 | 11 | -1.348 | 0.000934 | 0.3402 |
| T/NK | TJ | MARVELD2 | 12 | 11 | +1.413 | 0.001883 | 0.4166 |
| T/NK | TJ | EPB41L4B | 12 | 11 | +1.404 | 0.002591 | 0.4166 |
| T/NK | TJ | TJP1 | 12 | 11 | +1.706 | 0.003841 | 0.4218 |
| T/NK | TJ | F11R | 12 | 11 | +0.869 | 0.004729 | 0.4229 |

## What this is not

- Not muscat mixed-model DE and not a cell-level Wilcoxon.
- Not a dual-high TACSTD2×CLDN4 score.
- Not a merge beyond GSE131907+GSE205335.
- Not evidence that CLDN4 *causes* IFN/MHC/TJ/chemokine change.
- Genome-wide FDR on n≈11 vs 11 is expected to be thin; family counts are the claim.

## Files

- `tables/de_all.tsv` — all genes, all contrasts
- `tables/de_families.tsv` — IFN/MHC/TJ/chemokine rows
- `tables/de_q4q1_combined_families.tsv` — headline DE table
- `tables/n_honest.tsv`
- `tables/sample_inventory.tsv`
- `figures/` — volcano, family heatmap, key-gene boxes, forest, CLDN4 strip, n bars

Reproduce:

```bash
# T/NK matrices (needs GEO files in /tmp/winpair_geo)
python3 methods/winpair_131907_205335_muscat_cldn4/build_tnk_pseudobulk.py
python3 methods/winpair_131907_205335_muscat_cldn4/analyze.py
```
