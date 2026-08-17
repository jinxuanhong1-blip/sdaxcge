# Methods — pair GSE189357+GSE205335 CLDN4 Slingshot/PAGA

## Scope

Public processed scRNA only. GSE189357 RAW 10x UMI tar + GSE205335 author UMI RDS. No FASTQ, no EGA, no GSE131907, no GSE207422. PR #459 combo-enum T/NK Spearman is not re-fit.

## Gates

- GSE189357: marker epithelium `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`. Author cell types are absent. Locked patients TD1–TD9.
- GSE205335: author `lineage.total == Epithelial cells` on non-normal tissue. Locked 22 patients from PR #459. Histology is mixed (ADC/SQ/SCLC/NUT).
- Cap ≤350 cells / patient (seed 1). QC: ≥200 genes, ≥500 UMI.

## Trajectory

- log1p(CP10k); 3000 HVG; 30 PCs; Harmony batch = dataset; 30 neighbors; Leiden r=0.6; PAGA.
- **Slingshot** (Street et al. 2018) on Harmony/PCA + Leiden. Start cluster = highest AT2 score among Leiden clusters **not** in the top tercile of mean CLDN4 (prefer author AT2 if present). Never the CLDN4-high cluster.
- DPT from the same start cluster is a companion clock only.
- Barrier/keratin = KRT8/18/19/7, CDKN1A, PLAUR (**CLDN4 out**).
- IFN = Hallmark IFNα ∩ IFNγ (**CLDN4 out**).
- TACSTD2 is a comparator. No dual-high gate.

## Inference

Patient is the unit. Spearman uses patients with ≥10 analysis cells. Cell-level ρ on a lineage is descriptive. Paired tertile extra requires ≥8 cells in both CLDN4-high and CLDN4-low arms.
