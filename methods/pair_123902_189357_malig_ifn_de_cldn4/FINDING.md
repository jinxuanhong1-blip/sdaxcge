# Pair GSE123902+GSE189357: CLDN4-only malignant IFN/MHC/TJ DE

ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4. **Not CellChat.**
Tumor-cell-intrinsic program only (marker-malignant / epithelial cells).

Given cut (PR #459, **not re-audited**): GSE123902 + GSE189357 malignant
CLDN4 %pos vs same-unit T/NK, n=22, Spearman ρ=−0.638 (p=0.003).
This extra does **not** re-audit that T/NK ρ.

**Thesis (already correct):** CLDN4 KD / low raises the malignant cell's own
IFN / MHC-I; CLDN4-high should be IFN/MHC down, TJ up.

**Method: patient-pseudobulk OLS on log2(TMM-CPM+1).** Not muscat.
Patient/donor is the unit (GSE123902 donor-level; GSE189357 patient-level).
p-values are descriptive. Genome-wide FDR is thin on these n.

## Design (locked)

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | GSE123902 + GSE189357 only | The PR #459 pair that differs; not a bigger merge |
| Split | Within-cohort marker-malignant CLDN4 **%pos Q4 vs Q1** | Quartile tails, not median; mid quartiles unused in binary DE |
| Tails | Combined Q1 n=7 / Q4 n=5 | **Thin.** GSE189357 Q4 n=2 — single-cohort binary DE skipped |
| Continuous | CLDN4 %pos z, all 22 units, cohort covariate | Linear; not a causal model |
| Malignant | Marker-malignant UMI-sum: (EPCAM\|KRT8\|KRT18\|KRT19)>0 and PTPRC==0 | Not author CNV-malignant; same gate as PR #459 |
| T/NK | not a DE compartment here | T/NK ρ taken as given from PR #459; not re-audited |
| Model | ~ cohort + CLDN4_Q4 (combined); ~ CLDN4_Q4 (GSE123902 only) | Small n; no muscat mixed model; no voom weights |
| Families | IFN (Hallmark IFNα/γ), MHC-I/APM (custom), TJ (KEGG/GO + junction focal; **CLDN4 held out**) | MHC-II and chemokine panels are not the claim |

## Honest n

PR #459 pair labels (marker-malignant %pos; T/NK ρ not re-scored):

- GSE123902: n=13 tumor/met donors (normals dropped). Q1=4 Q4=3.
- GSE189357: n=9 patients. Q1=3 Q4=2.

Combined Q1: LX675, LX682, LX699, LX701, TD2, TD4, TD7 (n=7).
Combined Q4: LX653, LX680, LX684, TD6, TD9 (n=5).

| contrast | n_low / n_high | n_genes | note |
|---|---:|---:|---|
| Q4 vs Q1 combined | 7/5 | 13754 | cohort covariate; **thin tails** |
| Q4 vs Q1 GSE123902 | 4/3 | 11874 | donor-level; thin |
| Q4 vs Q1 GSE189357 | 3/2 | 0 | thin tail n_Q1=3 n_Q4=2 (need ≥3 each) |
| continuous combined | n=22 | 14593 | CLDN4 %pos z; all gated units |

Q1 includes LX699 (46 malignant cells) and LX701 (90). Those units stay in
because they are in the given n=22 pair. Do not quote a cell-level n.

## IFN / MHC-I/APM / TJ (combined Q4 vs Q1)

Positive logFC = higher in CLDN4-high (Q4) than CLDN4-low (Q1).
Expected under the thesis: IFN down, MHC-I/APM down, TJ up.

| family | n_tested | n p<0.05 (up/down) | n FDR<0.05 | median logFC | top gene (logFC, p, FDR) |
|---|---:|---|---:|---:|---|
| IFN | 216 | 45 (4/41) | 0 | -0.463 | B2M (-1.340, 0.0007885, 0.2) |
| MHC-I/APM | 21 | 7 (0/7) | 0 | -0.724 | B2M (-1.340, 0.0007885, 0.2) |
| TJ | 184 | 29 (17/12) | 0 | +0.102 | TBCD (+1.048, 0.0002197, 0.169) |

CLDN4 itself +1.18 (p=0.06075, FDR=0.4372) — direction check on the split gene (held out of the TJ family).

Headline family genes (combined Q4 vs Q1, lowest p):

| family | gene | n_Q1 | n_Q4 | logFC | p | FDR |
|---|---|---:|---:|---:|---:|---:|
| TJ | TBCD | 7 | 5 | +1.048 | 0.0002197 | 0.169 |
| IFN|MHC-I/APM | B2M | 7 | 5 | -1.340 | 0.0007885 | 0.2 |
| IFN | IL15RA | 7 | 5 | -2.098 | 0.000861 | 0.2035 |
| IFN | SAMHD1 | 7 | 5 | -1.258 | 0.0008687 | 0.2035 |
| IFN | ARID5B | 7 | 5 | -1.505 | 0.0008935 | 0.2035 |
| TJ | ARPC2 | 7 | 5 | -0.922 | 0.0009766 | 0.2035 |
| IFN | ELF1 | 7 | 5 | -0.803 | 0.001167 | 0.2127 |
| IFN | CSF2RB | 7 | 5 | -2.241 | 0.001827 | 0.2173 |
| TJ | SCRIB | 7 | 5 | +0.833 | 0.002378 | 0.231 |
| TJ | CLDN3 | 7 | 5 | +1.435 | 0.002612 | 0.231 |
| IFN | MYD88 | 7 | 5 | -0.631 | 0.003732 | 0.2551 |
| TJ | FZD5 | 7 | 5 | +1.375 | 0.003982 | 0.2551 |

### Continuous CLDN4 %pos (n=22, cohort covariate)

Same sign as Q4 vs Q1. This is the better-powered direction check when the
Q4 tail is only 5 donors/patients.

| family | n_tested | n p<0.05 (up/down) | n FDR<0.05 | median logFC | top gene (logFC, p, FDR) |
|---|---:|---|---:|---:|---|
| IFN | 218 | 47 (9/38) | 4 | -0.140 | CSF2RB (-1.048, 1.16e-05, 0.013) |
| MHC-I/APM | 21 | 8 (1/7) | 0 | -0.225 | B2M (-0.498, 0.002641, 0.08963) |
| TJ | 187 | 38 (24/14) | 1 | +0.074 | TBCD (+0.548, 2.58e-05, 0.01695) |

## What this is not

- Not CellChat / LIANA / NicheNet.
- Not a re-audit of the PR #459 T/NK ρ.
- Not muscat mixed-model DE and not a cell-level Wilcoxon.
- Not a dual-high TACSTD2×CLDN4 score.
- Not a merge beyond GSE123902+GSE189357.
- Not evidence that CLDN4 *causes* IFN/MHC/TJ change.
- Q4 n=5 (and GSE189357 Q4 n=2) is a thin tail; genome-wide FDR on 7 vs 5
  is expected to be empty. Family median logFC and sign counts are the claim.

## Files

- `tables/de_q4q1_combined_families.tsv` — headline family DE table
- `tables/family_summary.tsv` — IFN / MHC-I/APM / TJ counts
- `tables/de_families.tsv` — family rows, all contrasts
- `tables/de_all.tsv` — all genes, all contrasts
- `tables/n_honest.tsv`
- `tables/sample_inventory.tsv`
- `figures/` — volcano, family heatmap, key-gene boxes, forest, family-median bar, CLDN4 strip, n bars

Reproduce:

```bash
python3 methods/pair_123902_189357_malig_ifn_de_cldn4/download.py
python3 methods/pair_123902_189357_malig_ifn_de_cldn4/build_malignant_pseudobulk.py
python3 methods/pair_123902_189357_malig_ifn_de_cldn4/analyze.py
```
