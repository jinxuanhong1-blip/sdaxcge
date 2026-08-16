# w200 / A5 — real public KL lung: GSE165641

**Request:** TROP2-high immune-resistant subset; GSE76628 only if it is actually lung KL;
otherwise find a real public KL lung set. Honest.

**GSE76628 is not KL lung** (flank-skin Ad-VEGF / gastric-stroma arrays in nude mice).
Rejection note: `results/w200/A5_GSE76628/README.md`.

## Dataset used

| Field | Value |
|---|---|
| Accession | [GSE165641](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165641) |
| Paper | Wang & Zhong; PMID [34369094](https://pubmed.ncbi.nlm.nih.gov/34369094/) |
| Model | *KrasLSL-G12D/+; Lkb1fl/fl* (KL) GEMM, Ad-Cre, 10 weeks |
| Tissue | Lung tumor sections, RBC-depleted cells |
| Assay | 10X scRNA-seq, 2 mice (KL1 4,400 cells; KL2 2,780 cells) |
| Matrix | Author-processed log-normalized GEO file (31,053 × 7,180) |

This is the closest **public** series that is actually KL lung **and** contains both
epithelial and lymphocyte compartments, so "Tacstd2-high tumor subset vs CD8/NK"
is computable. No private KL mice were used.

Other public KL lung options that were **not** used as the primary set: GSE69552
(bulk beadchip of KL tumors — no cell-level contrast); GSE180963 (1 Kras + 1 KL
scRNA library, unpublished); GSE165640 (bulk whole-lung KL/KP — composition only).

## Honest verdict

1. **Tacstd2 is an epithelial gene in KL lung, not a CD8/NK program.** Detection
   64.7% in marker-defined epithelium (n=468) vs 5.0% in CD8 T (n=100) vs 2.0%
   in NK (n=50). Epithelial vs CD8: Mann–Whitney p = 2.9e-24, rank-biserial +0.62,
   Fisher OR for detection = 35 (p = 9.4e-31). Epithelial vs NK: p = 5.0e-14,
   r_rb +0.62, OR = 90 (p = 1.4e-19). CD8 vs NK do not differ (p = 0.39).
2. **There is no evidence for a Tacstd2-high immune-resistant epithelial subset.**
   Among epithelial cells, the Tacstd2-high tertile is *Cldn4*/*Epcam*/*Krt8*-high
   (junctional epithelium) and is **not** MHC/IFN-low. *H2-D1*, *Tap1*, *Psmb8*,
   and *Irf1* are higher in the high tertile (FDR < 0.01). *Cd274* (PD-L1) is
   near-zero in both tertiles. *Cxcl9* is undetectable in epithelium.
3. **The KL TME is lymphocyte-poor and neutrophil-rich**, which is the known
   STK11/LKB1-loss phenotype — a **genotype/TME** fact, not a Tacstd2-high
   subset fact. CD8+NK are 150 / 7,180 cells (2.1%). KL1 is neutrophil-dominated
   (2,098 / 4,400); KL2 has more T/NK.

Nothing here supports (or needs) a claim that TROP2-high KL tumor cells are a
distinct immune-resistant subset relative to TROP2-low KL tumor cells.

## Cell typing (limitations stated)

No published cluster labels ship with the GEO matrix. Types are mutually exclusive
marker rules (`scripts/w200_a5_gse165641.py`):

| Type | Rule | n |
|---|---|---|
| Epithelial | (Epcam+ or (Krt8+ and Krt18+) or Nkx2-1+) and Ptprc− | 468 |
| CD8_T | (Cd8a+ or Cd8b1+) and (Cd3e+ or Cd3d+) | 100 |
| NK | (Ncr1+ or Klrb1c+) and not CD8 | 50 |
| Other_T | Cd3e+ or Cd3d+, not CD8 | 165 |
| Neutrophil | Ptprc+ and S100a8 ≥ 2 | 2,339 |
| Myeloid | Ptprc+ and (Cd68+ or Adgre1+), not the above | 1,569 |
| Other_immune / Other | remainder | 583 / 1,906 |

Epithelial requires Ptprc− so Epcam+ immune droplets are not counted as tumor.
This is **not** CNV-based malignant calling: residual AT2/club epithelium can
enter the epithelial bin (`Sftpc` is slightly higher in the Tacstd2-low tertile).

**Ambient / doublet caveat:** 539 Tacstd2+ Ptprc+ droplets exist, almost all
neutrophil/myeloid (84% S100a8+), almost none CD8/NK. In a neutrophil-rich KL
TME this is the expected ambient-epithelial or doublet pattern, not a
lymphocyte Tacstd2 program. Neutrophil Tacstd2 detection (15%) should not be
read as Trop2+ neutrophils.

## Primary contrast

From `tacstd2_epithelial_vs_cd8_nk.csv`:

| Comparison | n | mean Tacstd2 | detection | MW p | r_rb |
|---|---|---|---|---|---|
| Epithelial vs CD8_T | 468 vs 100 | 0.95 vs 0.04 | 65% vs 5% | 2.9e-24 | +0.62 |
| Epithelial vs NK | 468 vs 50 | 0.95 vs 0.03 | 65% vs 2% | 5.0e-14 | +0.62 |
| CD8_T vs NK | 100 vs 50 | 0.04 vs 0.03 | 5% vs 2% | 0.39 | +0.03 |

## Within-epithelium Tacstd2 tertiles (n=156 high vs 156 low)

Selected rows from `epithelial_tacstd2_high_vs_low.csv` (BH-FDR over the
pre-specified panel):

| Gene | mean high | mean low | r_rb | FDR |
|---|---|---|---|---|
| Cldn4 | 1.17 | 0.25 | +0.55 | 5.4e-19 |
| Epcam | 2.44 | 1.42 | +0.52 | 9.8e-15 |
| Tap1 | 0.15 | 0.04 | +0.27 | 1.1e-8 |
| Irf1 | 0.31 | 0.24 | +0.27 | 8.6e-6 |
| H2-D1 | 1.95 | 1.36 | +0.26 | 2.3e-4 |
| Psmb8 | 0.42 | 0.36 | +0.21 | 1.0e-3 |
| Sftpc | 0.89 | 1.53 | −0.15 | 0.019 |
| Cd274 | 0.007 | 0.037 | −0.00 | 0.96 |
| Cxcl9 | 0 | 0 | — | not tested (all zero) |

A Tacstd2-high immune-resistant subset would have predicted **lower** MHC/IFN
in the high tertile. The opposite (or null) is what is observed.

## What this does not do

- Does not identify malignant vs residual normal epithelium by CNV.
- Does not use author clusters (none provided).
- Does not re-align FASTQ; uses the GEO processed matrix.
- Is two mice only; KL1 vs KL2 immune composition differs.
- Has no ICI treatment arm, so "immune-resistant" here means tumor-cell MHC/IFN
  / lymphocyte identity, not checkpoint-blockade outcome.
- Does not analyze private KL mice.

## Files

| File | Contents |
|---|---|
| `tacstd2_by_celltype.csv` | n, detection, mean/median/p90 Tacstd2 per type |
| `tacstd2_epithelial_vs_cd8_nk.csv` | Primary MW + Fisher tests |
| `epithelial_tacstd2_high_vs_low.csv` | Within-epithelium tertile tests + FDR |
| `celltype_composition.csv` | Counts per mouse |
| `per_cell_markers.csv` | Per-cell type + key genes |
| `epithelial_cells.csv` | Epithelial barcodes and Tacstd2 tertile |
| `tacstd2_in_immune_droplets.csv` | Ambient/doublet caveat |
| `fig_tacstd2_by_celltype.png` | Tacstd2 boxplots by type |
| `fig_detection_and_epithelial_highlow.png` | Detection rates + high vs low genes |
| `provenance.json` | URL, MD5, rules, rejection of GSE76628 |

## Reproduce

```bash
pip install pandas numpy scipy matplotlib rdata
python3 scripts/w200_a5_gse165641.py   # downloads the GEO Rdata into data/kl_lung/ if absent
```
