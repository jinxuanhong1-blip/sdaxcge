# Methods — CLDN4-only high-end CellChat on GSE123902 + GSE189357

## Given (not re-audited)

PR #459 `methods/scrna_cldn4_combo_enum`, pair that **differs**:

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE123902+GSE189357 | %pos | 22 | −0.638 (0.003, 0%) | −1.000 (0.003); tails **7/5 thin** |

Members: GSE123902 marker-malignant donors n=13; GSE189357 marker-malignant patients n=9.

This slice does **not** recompute that Spearman. Dual-high TACSTD2∩CLDN4 is not used.
The 7-cohort pool is not the answer (PR #459). GSE148071 / GSE127465 do not differ and are not added.

## Cohorts

Public processed GEO only. No FASTQ. Matrices are **not** concatenated.

- **GSE123902** (Laughney et al., *Nat Med* 2020). Dense CSV in `GSE123902_RAW.tar` (~90 MB). Tumor/met libraries only; NORMAL dropped. Donor is the unit (PRIMARY preferred if both exist; none do). Marker-malignant.
- **GSE189357** (Zhu/Wang AIS–IAC atlas). 10x MTX in `GSE189357_RAW.tar` (~624 MB). One tumor sample per patient (TD1–TD9). Marker-malignant.

Locked marker gates (same as PR #459):

- malignant-like = `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`
- T/NK = `(CD3D|CD3E|CD8A|NKG7|GNLY|KLRD1)>0` and not malignant-like

## Estimand

Within each patient/donor, outgoing ligand–receptor score from **CLDN4-high vs CLDN4-low**
marker-malignant cells toward **same-unit T/NK**.

Patient/donor is the unit. Honest n = units that pass floors, which may be **< 22**.
The given between-patient Q4 vs Q1 (7/5) is thin and is **not** the CellChat n.

## Split and floors

- Score = `log1p(UMI / total × 10⁴)` CLDN4. **CLDN4 only.**
- Primary: within-unit **median**. Extra high-end: **tertile** (middle third dropped) and **Q4 vs Q1**.
- Floors: ≥10 high and ≥10 low malignant; ≥20 T/NK.
- Q4 vs Q1 additionally requires ≥40 malignant cells. Units that pass the floor with thin tails (high or low <20 cells) are kept and **flagged**.

## Scores

**CellChat-style** (Jin et al. 2021; CellChat R not run): 10% truncated mean, Hill \(K_h=0.5\), CellChatDB v2 protein pairs, `expr_prop ≥ 0.10`. Complexes = geometric mean of subunits.

**LIANA-style extra** (Efremova 2020; LIANA package optional): CellPhoneDB mean-of-means on log1p(CP10k), CellPhoneDB v5 pairs, `expr_prop ≥ 0.10`. Complexes = min of subunits.

Inference: paired Wilcoxon on patient high vs low scores. BH-FDR inside each method × rule. Cells are not replicates.

## Cannot test

Author malignant labels (none published as a compact column on these two tars), histologic TLS, CopyKAT, dual-high gates, ICI/MPR on these accessions.
