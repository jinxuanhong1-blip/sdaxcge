# Finding — Seurat GSE189357 CLDN4-only malignant vs T/NK and IFN/MHC/TJ

ADDITIVE. **CLDN4 only.** Zhu et al., *Exp Mol Med* 2022 (DOI 10.1038/s12276-022-00896-9), GEO **GSE189357**: nine treatment-naïve resected LUAD lesions (TD1–TD9; AIS=3, MIA=3, IAC=3). Public processed 10x MTX. **R + Seurat** `CreateSeuratObject`. No TACSTD2∩CLDN4 dual-high. Not GSE148071. Not a Python-only primary. Not a concordant-4 pool redo.

Thesis (already correct; not re-derived): CLDN4-high malignant cells have lower own IFN/MHC-I, and patients have lower T/NK. This folder is the **GSE189357-only Seurat** slice. Matching extras here are T/NK **DOWN** and malignant IFN/MHC **DOWN** as CLDN4 rises. Do not sell a dual-high or GSE148071 row as this answer.

Honest unit = **patient**. One 10x tumor each. **n may be 9 — say so.** Cell counts are descriptive. p-values are descriptive. Marker-malignant ≠ CNV.

## Verdict

Patient-level malignant CLDN4 %pos vs T/NK fraction: n=9, ρ=-0.600, p=0.0876 (n may be 9). CLDN4 mean vs T/NK: n=9, ρ=-0.400, p=0.2861. CLDN4 %pos vs malignant IFN (Hallmark α∪γ, patient-pseudobulk log2 CPM+1): n=9, ρ=-0.367, p=0.3317. vs MHC-I/APM: n=9, ρ=-0.250, p=0.5165. vs TJ (CLDN4 held out): n=9, ρ=-0.400, p=0.2861. Q4 vs Q1 T/NK: skipped: need >=3 patients in each tail. Patient is the unit. Not a TACSTD2 redo. No both-high gate.

## Honest n

- GEO catalog: **n_patients = 9** (TD1–TD9). This is the catalog n and the Spearman n. **Say so.**
- Stage split: AIS=3, MIA=3, IAC=3 (n=3/stage — not a stage test).
- Cells after CreateSeuratObject (no extra QC drop): **n_cells = 122373**.
- Marker-malignant: **n_cells = 14252**. Gate `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`.
- T/NK: **n_cells = 62770**. Gate `(CD3D|CD3E|CD8A|NKG7|GNLY|KLRD1)>0` and not malignant-like.
- Patients with ≥20 malignant cells used for Spearman: **n = 9**.
- Malignant cells per patient: {TD1=2186, TD2=1207, TD3=664, TD4=491, TD5=1736, TD6=1211, TD7=849, TD8=1019, TD9=4889}.
- T/NK cells per patient: {TD1=9915, TD2=11334, TD3=6525, TD4=7096, TD5=11055, TD6=3392, TD7=1851, TD8=7257, TD9=4345}.
- Locked family genes present (first sample, shared features): IFN 224/224 (WARS/MARCH1 stand in for WARS1/MARCHF1); MHC-I/APM 21/21; TJ 211/214 (CLDN4 out).
- Genes absent after alias resolution: IFN []; MHC []; TJ [MIR105-1, MIR142, PALS1].
- Seurat 5.0.1; CreateSeuratObject used. Harmony / clustering not required for this patient-level claim.

## Locked choices

- Input: GEO `GSE189357_RAW.tar` processed MTX/TSV. No FASTQ. No GSE189487 spatial.
- Object: per-patient `ReadMtx` → `CreateSeuratObject` → `NormalizeData` (log1p CP10k).
- Malignant: marker gate, **not CNV**.
- T/NK: locked PR #459 marker gate, not author labels (none published as a compact column on this tar).
- CLDN4 score primary = malignant **%pos**. Mean is the matching extra.
- IFN = Hallmark IFNα ∪ IFNγ (CLDN4 not in the set).
- MHC-I/APM = locked custom classical MHC-I / APM (MHC-II excluded).
- TJ = KEGG tight junction ∪ GOBP TJ organization ∪ focal TJ genes; **CLDN4 held out**.
- Family scores primary = patient-pseudobulk log2(UMI-sum CPM + 1) mean of present genes.
- AddModuleScore and cell-mean scores are extras.
- Q4 vs Q1 requires ≥3 patients in each tail. n=9 tails are thin.

## Primary (patient-level Spearman, BH inside this list)

| contrast | n_patients | rho | p | q |
| --- | --- | --- | --- | --- |
| CLDN4 %pos vs T/NK fraction | 9 | -0.600 | 0.0876 | 0.4147 |
| CLDN4 mean vs T/NK fraction | 9 | -0.400 | 0.2861 | 0.4147 |
| CLDN4 %pos vs malignant IFN (pseudobulk) | 9 | -0.367 | 0.3317 | 0.4147 |
| CLDN4 %pos vs malignant MHC-I/APM (pseudobulk) | 9 | -0.250 | 0.5165 | 0.5165 |
| CLDN4 %pos vs malignant TJ (pseudobulk, CLDN4 out) | 9 | -0.400 | 0.2861 | 0.4147 |

n may be 9. A significant p at n=9 is a small-n result, not a cohort.

## Q4 vs Q1 (not the headline)

Within-cohort CLDN4 %pos quartiles: n_Q1=3 n_Q4=2. skipped: need >=3 patients in each tail.

## Extra — cell-mean / AddModuleScore / TACSTD2 comparator

| contrast | family | n_patients | rho | p |
| --- | --- | --- | --- | --- |
| CLDN4 %pos vs malignant IFN (cell mean) | extra | 9 | 0.317 | 0.4064 |
| CLDN4 %pos vs malignant MHC-I/APM (cell mean) | extra | 9 | -0.200 | 0.6059 |
| CLDN4 %pos vs malignant TJ (cell mean, CLDN4 out) | extra | 9 | 0.050 | 0.8984 |
| CLDN4 %pos vs TACSTD2 %pos (comparator) | comparator | 9 | 0.850 | 0.0037 |
| CLDN4 %pos vs AddModuleScore IFN | extra | 9 | 0.483 | 0.1875 |
| CLDN4 %pos vs AddModuleScore MHC-I/APM | extra | 9 | -0.050 | 0.8984 |
| CLDN4 %pos vs AddModuleScore TJ | extra | 9 | -0.450 | 0.2242 |

## Extra — within-patient CLDN4-high vs low (median split)

Patients with ≥10 high and ≥10 low malignant cells: **n = 9**.

| contrast | n_patients | delta_median | p |
| --- | --- | --- | --- |
| within-patient IFN high-low | 9 | 0.026 | 0.0244 |
| within-patient MHC-I/APM high-low | 9 | 0.024 | 0.8127 |
| within-patient TJ high-low | 9 | 0.064 | 0.0092 |

Delta = high − low cell-mean family score. Observed IFN is slightly **up** in CLDN4-high (Δmed=+0.026, p=0.0244) — opposite the KD / between-patient IFN-down direction, and not sold as the claim. TJ is higher in CLDN4-high (Δmed=+0.064, p=0.0092). Between-patient pseudobulk IFN/MHC remain the matching extras (negative ρ, n=9, not significant). This is not the between-patient n=9 claim.

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- **n may be 9.** This is not the concordant-4 N=65 pool and not a stage-powered test (n=3/stage).
- Malignant is a marker gate, **not CNV** (inferCNV was not run).
- No TACSTD2∩CLDN4 both-high gate.
- Not GSE148071 / GSE127465 / GSE207422 / GSE154826.
- Not evidence that CLDN4 *causes* the T/NK or IFN/MHC change.
- Not a Python-only primary. Seurat `CreateSeuratObject` was run.

## Outputs

- `results/tables/patient_units.tsv` — **headline patient table** (honest n)
- `results/tables/spearman_patient.tsv`
- `results/tables/q4q1_tnk.tsv`
- `results/tables/paired_high_low.tsv`
- `results/tables/gene_coverage.tsv`
- `results/figures/fig_scatter_cldn4_tnk.png`
- `results/figures/fig_scatter_cldn4_programs.png`
- `results/figures/fig_honest_n.png`
- `results/figures/fig_cldn4_by_stage.png`

## Reproduce

```bash
bash methods/seurat_gse189357_cldn4/scripts/download.sh /tmp/gse189357
Rscript methods/seurat_gse189357_cldn4/scripts/analyze.R \
  --tar /tmp/gse189357/GSE189357_RAW.tar \
  --extract /tmp/gse189357/raw \
  --outdir methods/seurat_gse189357_cldn4/results \
  --finding methods/seurat_gse189357_cldn4/FINDING.md
```

