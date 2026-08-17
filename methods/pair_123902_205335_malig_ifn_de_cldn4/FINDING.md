# Pair GSE123902+GSE205335: CLDN4-only malignant-cell IFN/MHC/TJ/keratin DE

ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4. Not CellChat. Not T/NK infiltrate.
No GSE148071. Patient (GSE123902: donor) is the unit.

The pair T/NK Spearman is **taken as given** from PR #459 and is not re-audited:

| combo | score | N | ρ (p) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE123902+GSE205335 | %pos | 35 | −0.522 (0.002) | −0.802 (0.005; 9 vs 9) |

This folder is **tumor-cell-intrinsic program DE** on within-cohort quartiles of
the locked %pos vector: patient-pseudobulk of **malignant cells**, CLDN4 %pos
Q4 vs Q1. Thesis (already correct, not audited here): CLDN4-low malignant cells
raise their own IFN / MHC-I program and lower TJ; observationally CLDN4-high
malignant = IFN/MHC down, TJ up. The given T/NK 9 vs 9 is a different contrast.

**Method: patient-pseudobulk OLS on log2(TMM-CPM+1).** Not muscat. muscat `pbDS`
is the same collapse (patient × cell-type UMI-sum → bulk DE). No R/limma/edgeR
in this environment; TMM + OLS + BH is the implementation. p-values are descriptive.

## Design (locked)

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | GSE123902 + GSE205335 only | Not GSE148071 / 131907 / 189357 / 7-cohort pool |
| Split | Within-cohort locked malignant CLDN4 **%pos Q4 vs Q1** | Not the PR #459 pooled-pair T/NK 9 vs 9; mid unused |
| Continuous | CLDN4 %pos z, all gated units with counts, cohort covariate | Linear; not causal |
| Malignant | GSE123902 marker-malignant UMI-sum; GSE205335 author malignant UMI-sum | GSE123902 has no author malignant labels |
| T/NK | **not run** | Given ρ is infiltrate, not this DE |
| Model | ~ cohort + CLDN4_Q4 (combined); ~ CLDN4_Q4 (single) | Small n; no voom weights |
| Families | IFN (Hallmark IFNα/γ); MHC-I/APM (custom); TJ (KEGG/GO + focal, **CLDN4 held out**); keratin (GOBP keratinization + KRT epithelial) | Compact custom MHC-I, not all HLA |

GSE123902 malignant = (EPCAM\|KRT8\|KRT18\|KRT19)>0 and PTPRC==0; tumor/met
donors; normals dropped; one library per donor (PRIMARY preferred). GSE205335
malignant = author `Malignant cells` (existing patient UMI-sum; P4001 n_mal=27
is absent from that matrix).

## Honest n

Locked labels (PR #459 %pos; T/NK n=35 is **not** the DE n):

- GSE123902: n=13 donors. Q1=4 Q4=3.
- GSE205335: n=22 patients. Q1=6 Q4=6.

| contrast | n_low / n_high | n_genes | note |
|---|---:|---:|---|
| Q4 vs Q1 combined malignant | 9/9 | 14974 | cohort covariate; units with counts |
| Q4 vs Q1 GSE123902 | 4/3 | 11874 | donor; thin tails |
| Q4 vs Q1 GSE205335 | 5/6 | 18035 | P4001 out of malignant counts |
| continuous combined | n=34 | 15867 | CLDN4 %pos z |

Combined malignant DE is also 9 vs 9 **after P4001 drops** (n_mal=27; Q1 on
GSE205335). That is not the given T/NK 9 vs 9. Per-contrast patients:
`tables/n_honest.tsv`.

## Family-score DE (headline table)

One OLS per family on the **mean log2(TMM-CPM+1)** of family genes. Positive
logFC = higher in CLDN4-high (Q4). Patient/donor is the unit.

### Combined Q4 vs Q1 (cohort covariate)

| family | n_Q1 | n_Q4 | n | n_genes | logFC | p | FDR | MW p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| IFN | 9 | 9 | 18 | 217 | -1.046 | 0.0003161 | 0.001264 | 0.001086 |
| MHC-I/APM | 9 | 9 | 18 | 21 | -1.285 | 0.00366 | 0.00732 | 0.00268 |
| TJ | 9 | 9 | 18 | 193 | -0.050 | 0.6838 | 0.6838 | 0.5962 |
| keratin | 9 | 9 | 18 | 40 | -0.740 | 0.06914 | 0.09218 | 0.04226 |

CLDN4 itself +0.98 (p=0.0395, n=9/9) — direction check on the split gene; held out of the TJ family.

IFN and MHC-I/APM family means are lower in CLDN4-high malignant cells. TJ
(CLDN4 held out) is **flat** as a family mean (gene-level mixed). Keratin
family-mean trends down. Do not quote a TJ-up family logFC on this pair.

Machine table: `tables/family_de.tsv` (combined n / logFC / p). All splits:
`tables/family_score_de.tsv`.

## Gene-level family DE (combined Q4 vs Q1)

Positive logFC = higher in CLDN4-high malignant cells.

| family | n_Q1 | n_Q4 | n_tested | n p<0.05 (up/down) | n FDR<0.05 | median logFC | top gene (logFC, p, FDR) |
|---|---:|---:|---:|---|---:|---:|---|
| IFN | 9 | 9 | 217 | 89 (1/88) | 0 | -1.009 | CCL5 (-3.218, 0.000113, 0.07238) |
| MHC-I/APM | 9 | 9 | 21 | 13 (0/13) | 0 | -1.343 | TAP1 (-2.007, 0.002312, 0.1493) |
| TJ | 9 | 9 | 193 | 28 (13/15) | 0 | +0.011 | PECAM1 (-2.376, 0.002042, 0.1433) |
| keratin | 9 | 9 | 40 | 2 (0/2) | 0 | -0.469 | KRT6A (-4.513, 0.003009, 0.1595) |

Headline genes (family members only, lowest p):

| family | gene | n_Q1 | n_Q4 | logFC | p | FDR |
|---|---|---:|---:|---:|---:|---:|
| IFN | CCL5 | 9 | 9 | -3.218 | 0.000113 | 0.07238 |
| IFN | GZMA | 9 | 9 | -2.911 | 0.0001687 | 0.07238 |
| IFN | GBP2 | 9 | 9 | -3.485 | 0.0001727 | 0.07238 |
| IFN | OAS2 | 9 | 9 | -3.130 | 0.0002051 | 0.07437 |
| IFN | PSME2 | 9 | 9 | -1.036 | 0.0002086 | 0.07437 |
| IFN | CASP1 | 9 | 9 | -2.569 | 0.0004467 | 0.0924 |
| IFN | IL15RA | 9 | 9 | -2.375 | 0.0005273 | 0.0924 |
| IFN | SAMD9L | 9 | 9 | -2.456 | 0.0005726 | 0.0924 |
| IFN | EPSTI1 | 9 | 9 | -2.874 | 0.0006292 | 0.09705 |
| IFN | TNFSF10 | 9 | 9 | -2.226 | 0.001027 | 0.1216 |
| IFN | MYD88 | 9 | 9 | -1.049 | 0.001216 | 0.1233 |
| IFN | IFI35 | 9 | 9 | -1.588 | 0.001217 | 0.1233 |

## What this is not

- Not CellChat / LIANA / NicheNet and not a T/NK infiltrate re-audit.
- Not muscat mixed-model DE and not a cell-level Wilcoxon.
- Not a dual-high TACSTD2×CLDN4 score.
- Not GSE148071 and not a bigger merge.
- Not evidence that CLDN4 *causes* IFN/MHC/TJ/keratin change.
- Genome-wide FDR on these n is expected to be thin; family n / logFC / p is the claim.

## Extra figures

- `figures/volcano_malignant_q4q1_combined.png`
- `figures/heatmap_malignant_families_q4q1.png`
- `figures/box_malignant_key_genes.png`
- `figures/box_family_scores_q4q1.png`
- `figures/forest_malignant_families_q4q1.png`
- `figures/forest_family_scores_q4q1.png`
- `figures/cldn4_quartile_strip.png`
- `figures/n_honest_q4q1.png`

## Files

- `tables/family_de.tsv` — combined family DE (n / logFC / p)
- `tables/family_score_de.tsv` — family scores including single-cohort
- `tables/de_q4q1_combined_families.tsv` — gene-level family DE
- `tables/de_families.tsv` / `tables/de_all.tsv.gz`
- `tables/family_summary.tsv` / `tables/n_honest.tsv` / `tables/sample_inventory.tsv`

Reproduce:

```bash
python3 methods/pair_123902_205335_malig_ifn_de_cldn4/scripts/download_gse123902.py
python3 methods/pair_123902_205335_malig_ifn_de_cldn4/scripts/build_gse123902.py
python3 methods/pair_123902_205335_malig_ifn_de_cldn4/analyze.py
```
