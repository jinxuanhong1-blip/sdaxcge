# FINDING — Seurat GSE131907 CLDN4-only (patient unit)

ADDITIVE. **CLDN4-only.** No TACSTD2∩CLDN4 dual-high. Not GSE148071.
Kim et al., *Nat Commun* 2020, PMID 32385277 ([GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907)).
Thesis (already correct; not re-derived): CLDN4-high malignant cells have lower own
IFN/MHC-I, and patients have lower T/NK. CLDN4-low/KD opens IFN/MHC. Matching extras
here are IFN/MHC **DOWN** in CLDN4-high. Do not sell a DoRothEA IFN-up-in-high as the
KD direction.

Primary engine: **Seurat 5.5.1 `CreateSeuratObject`** on the public processed raw UMI
matrix (gene-subset of locked IFN/MHC/TJ/chemokine + anchors). Full-library UMI sizes
were streamed in the same pass and used for `NormalizeData` / CPM. The 2.86 GB log2TPM
text and EGA FASTQ were not used. Python-only is not the primary.

**Honest unit = patient.** Tumor-bearing sites only (tLung, tL/B, mLN, PE, mBrain).
Malignant = author `Cell_subtype == Malignant cells` (same author-malig gate as the
concordant-4 GSE131907 slice). Eligible if n_malignant ≥ 20. nLung has 0 author-malignant
cells and cannot join. p-values are descriptive.

## 1. CLDN4 vs T/NK fraction

Primary (patient, author Malignant, %pos): **n=21**, Spearman ρ=-0.478, p=0.02843.
Mean log-normalized CLDN4 vs T/NK: ρ=-0.532, p=0.01296.
Q4 vs Q1 T/NK (within-patient-set %pos quartiles): rank-biserial r=-0.600, n_Q1/n_Q4=6/5, p=0.1207; median T/NK Q1=0.428 Q4=0.048.
Unique patients in the eligible set: 21. Patients with >1 tumor-bearing sample in this slice: 4
(P1006, P1011, P1012, P1013). Each adds a PE capture with n_malignant=0. Sample-level T/NK
uses only the malignant-bearing capture (concordant-4 given ρ=−0.522). Patient-level T/NK
pools that PE immune into the same patient, so ρ softens to −0.478. That is the honest
patient unit, not a discrepancy.

| contrast | unit | malignant | n | rho %pos (p) | rho mean (p) | Q4 vs Q1 r (n_Q1/n_Q4, p) |
|---|---|---|---:|---|---|---|
| primary | patient | author Malignant cells | 21 | -0.478 (0.02843) | -0.532 (0.01296) | -0.600 (6/5, 0.1207) |
| sensitivity | patient | Malignant + tS1/tS2/tS3 | 31 | -0.343 (0.05877) | -0.392 (0.02939) | -0.406 (8/8, 0.1893) |
| comparison | sample | author Malignant cells | 21 | -0.522 (0.0152) | -0.594 (0.004564) | -0.600 (6/5, 0.1207) |

The sample-level author-Malignant row is the comparison to the prior GSE131907
concordant-4 sample n=21 (given ρ=−0.522). It is not a re-audit of the four-set pool.
tS1/tS2/tS3 are **not** in the primary gate (primary tLung uses those labels; author
`Malignant cells` are mets / tL/B / mLN in this atlas).

## 2. Malignant IFN / MHC-I / TJ — Q4 vs Q1

Method: **patient-pseudobulk log2(CPM+1)** of author-malignant UMI sums, using the
full-library UMI total as the size factor. Not TMM. Not muscat. Not a cell-level Wilcoxon.
Quartiles = the same patient %pos Q labels as the T/NK test. Positive logFC = higher in
CLDN4-high. CLDN4 is held out of TJ. Expected under the thesis: IFN down, MHC-I/APM down, TJ up.

Honest DE n = Q1+Q4 patients in the count collapse: **6 / 5** (Q1/Q4).
Do not quote 208,506 cells as n.

| family | n | n_Q1 / n_Q4 | n_genes | logFC | r (Q4 vs Q1) | p | expected |
|---|---:|---|---:|---:|---:|---:|---|
| IFN | 21 | 6 / 5 | 221 | -0.001 | -0.067 | 0.9273 | DOWN |
| MHC-I/APM | 21 | 6 / 5 | 21 | -0.164 | -0.133 | 0.7842 | DOWN |
| TJ (CLDN4 out) | 21 | 6 / 5 | 204 | 0.182 | 0.400 | 0.3153 | UP |
| chemokine | 21 | 6 / 5 | 26 | -0.390 | -0.667 | 0.08284 | DOWN |
| CLDN4 (split check) | 21 | 6 / 5 | 1 | 2.940 | 0.867 | 0.02248 | UP |

IFN is near-flat in this single cohort (logFC −0.001). That matches the published
GSE131907-only family row in the four-set table (already the near-null IFN set).
The four-set pooled IFN DOWN is not re-derived here. MHC-I/APM and chemokine are
DOWN and TJ is UP; tails are thin (6 vs 5). CLDN4 itself is UP (split check).

Gene families are the locked A8/concordant-4 sets (Hallmark IFN-α∪γ; custom MHC-I/APM;
KEGG TJ ∪ GOBP tight-junction organization minus CLDN4; compact chemokine panel).
Genes in the locked lists but absent from the UMI matrix: AFDN, CLDN24, MARCHF1, MIR105-1, MIR142, PALS1, PATJ, PRKACG, TENT5A, WARS1.

## Honest n

- CreateSeuratObject cells: **208506** (do not quote as n).
- GEO patients / samples: **44 / 58**.
- Author Malignant cells: **24784**. tS1+tS2+tS3: **6352**.
- PRIMARY eligible patients: **21** (Q1=6, Q4=5).
- Author-malignant cells in nLung: **0**.
- Seurat 5.5.1; NormalizeData LogNormalize on full-library size.
- Treatment-naive LUAD atlas: no ICI / MPR / RECIST labels.

## What this is not

- Not a mega-merge and not GSE148071 / GSE207422 / GSE205335 / GSE123902 / GSE189357.
- Not a dual-high TACSTD2×CLDN4 score.
- Not a re-derivation of the thesis and not a re-audit of the concordant-4 pooled n=65.
- Not muscat mixed-model DE and not a cell-level Wilcoxon sold as n.
- Not evidence that CLDN4 *causes* the T/NK or IFN/MHC change.
- Q4 vs Q1 tails are thin in a single cohort; family direction is the claim, not genome-wide FDR.

## Files

- `results/tables/tnk_tests.tsv` — **headline T/NK table**
- `results/tables/family_q4q1.tsv` — **headline IFN/MHC/TJ table**
- `results/tables/patient_units.tsv` — patient-level CLDN4 and T/NK
- `results/tables/n_honest.tsv`
- `results/figures/fig_cldn4_vs_tnk.png`
- `results/figures/fig_tnk_q4q1.png`
- `results/figures/fig_family_q4q1.png`
- `results/figures/fig_honest_n.png`

## Reproduce

```bash
bash methods/seurat_gse131907_cldn4/scripts/download.sh /tmp/gse131907
python3 methods/seurat_gse131907_cldn4/scripts/extract_gene_matrix.py \
  --datadir /tmp/gse131907 --outdir /tmp/gse131907/subset
Rscript methods/seurat_gse131907_cldn4/scripts/analyze_seurat.R \
  /tmp/gse131907/subset /tmp/gse131907 \
  methods/seurat_gse131907_cldn4/results \
  methods/seurat_gse131907_cldn4/FINDING.md \
  methods/seurat_gse131907_cldn4/data/families.json
```
