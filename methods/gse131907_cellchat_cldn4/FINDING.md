# FINDING — CellChat-style malignant CLDN4-high → T/NK (GSE131907)

**Additive only.** GSE207422 CellChat (`methods/scrna_cellchat`, `methods/scrna_cellchat_cldn4`) is a different agent and is not re-run. This slice is **GSE131907 + CLDN4**. TACSTD2 is not used to define groups.

Public Kim et al. treatment-naive LUAD atlas (*Nat Commun* 2020, PMID 32385277; GEO [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907)). CellChat R and LIANA were not run. Pairs are Jin et al. 2021 Hill probability on CellChatDB v2 protein pairs plus 100 permutations of CLDN4-high/low among malignant cells. A pair is significant if detected (`expr_prop ≥ 0.10` both sides), `P > 0`, and permutation `p < 0.05`. Smallest possible p = 1/101 = 0.0099.

Primary table: [`results/ligand_table.tsv`](results/ligand_table.tsv) (outgoing Mal → T/NK, kept split).

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| Barcodes × genes in the GEO UMI | yes | **208,506 × 29,634** | 27,163 unique symbols after stripping `.` suffixes |
| Samples / patients in the series | yes | **58 / 44** | GEO series matrix |
| ICI / MPR / RECIST labels | **no** | **0** | treatment-naive atlas; NMPR>MPR is not testable |
| Author malignant / tS* cells | yes | **31,136** | `Epithelial cells` ∩ `{Malignant cells, tS1, tS2, tS3}` |
| — of which tLung tS1/tS2/tS3 | yes | 6,352 | primary tumor epithelium is **not** labeled `Malignant cells` |
| — of which author `Malignant cells` | yes | 24,784 | tL/B 6,400 + mLN 2,961 + mBrain **15,423 (49.5%)** |
| PE epithelial (unlabeled subtype) | yes | 396 | **excluded** from malignant |
| Author T+NK in the whole atlas | yes | 91,227 | includes nLung / nLN / PE |
| **Kept: tumor-origin T+NK** | yes | **34,741** | tLung 19,591 + tL/B 3,515 + mLN 6,866 + mBrain 4,769 |
| **Kept: CLDN4-high malignant** | yes | **10,379 cells / 30 samples / 30 patients** | top tertile of `log1p(CP10k)` CLDN4 |
| **Kept: CLDN4-low malignant** | yes | **10,379 cells / 32 samples / 32 patients** | bottom tertile; middle third dropped |
| **Kept: T/NK partner** | yes | **34,741 cells / 32 samples / 32 patients** | same 32 tumor-origin samples |
| CellChatDB v2 protein pairs in matrix | yes | 2,116 / 2,239 | every subunit present |
| CopyKAT / author Seurat object | **no** | 0 | labels are the GEO annotation file |
| Patient-level mixed model | **no** | — | means are **cell-pooled** |

Do not write n=44 for the ligand test. The computable kept n is **32 tumor-origin samples / 32 patients**, with **10,379 vs 10,379** malignant cells. Two of the 32 samples have no cell in the global CLDN4-high tertile (Mal_high n_samples = 30).

Cell-pool dominance (not a patient mean):

| Group | Top sample | Cells | Fraction | Origin |
|---|---|---:|---:|---|
| Mal_high | EBUS_28 | 2,847 | **27.4%** | tL/B |
| Mal_low | NS_07 | 1,354 | 13.0% | mBrain |
| TNK | LUNG_T31 | 3,713 | 10.7% | tLung |

tLung-only sensitivity is **11 samples / 11 patients** (2,118 vs 2,118 malignant; 19,591 T/NK). Mets-only is **21 samples / 21 patients** (8,262 vs 8,262; 15,150 T/NK). LUNG_T34 is 31.7% of tLung Mal_high and 39.4% of tLung Mal_low.

## Kept split

**Tumor-origin CLDN4 tertile** (tLung + tL/B + mLN + mBrain). Tertile cuts on malignant `log1p(CP10k)` CLDN4: 1.45 / 2.11. Mean CLDN4: high 2.57, low 0.64, T/NK 0.06.

Normal lung, normal LN, and PE are out. PE T/NK (13,271) are out because PE epithelium is not author-malignant.

| Split | Mal high | Mal low | TNK | Samples (pts) | Detected LR | Sig LR (out / in) | Sig outgoing Δ |
|---|---:|---:|---:|---|---:|---|---:|
| **tumor tertile (kept)** | **10,379** | **10,379** | **34,741** | **32 (32)** | **164** | **81 (60 / 21)** | **60** |
| tLung tertile | 2,118 | 2,118 | 19,591 | 11 (11) | 159 | 78 (62 / 16) | 62 |
| mets tertile | 8,262 | 8,262 | 15,150 | 21 (21) | 151 | 72 (50 / 22) | 50 |

## Ligand table (outgoing Mal → T/NK)

Full 2,116-row table with honest n columns: [`results/ligand_table.tsv`](results/ligand_table.tsv). Figure: [`results/fig_top_outgoing_tumor_tertile.png`](results/fig_top_outgoing_tumor_tertile.png).

Among 2,116 pairs × 4 directed tests = 8,464 rows: detected 164; significant 81; outgoing Mal→TNK significant 60 (43 high-arm, 17 low-arm).

### Higher in CLDN4-high (p_high = 0.0099)

Barrier / inhibitory (pre-specified class):

| Pair | Pathway | Class | ΔP (high−low) | P high | P low |
|---|---|---|---:|---:|---:|
| HLA-E–CD8A | MHC-I | inhibitory | +0.057 | 0.352 | 0.295 |
| HLA-E–CD8B | MHC-I | inhibitory | +0.046 | 0.249 | 0.204 |
| LGALS9–CD45 | GALECTIN | inhibitory | +0.043 | 0.223 | 0.180 |
| CDH1–ITGAE/ITGB7 | CDH1 | barrier\|inhibitory | +0.039 | 0.107 | 0.068 |
| LGALS9–CD44 | GALECTIN | inhibitory | +0.038 | 0.188 | 0.150 |
| HLA-F–CD8A | MHC-I | inhibitory | +0.039 | 0.080 | 0.042 |
| JAM1–ITGAL/ITGB2 | JAM | barrier | +0.031 | 0.110 | 0.079 |
| HLA-F–CD8B | MHC-I | inhibitory | +0.025 | 0.051 | 0.026 |
| LGALS9–P4HB | GALECTIN | inhibitory | +0.016 | 0.072 | 0.056 |
| HLA-G–CD8A | MHC-I | inhibitory | +0.006 | 0.006 | 0 |
| CDH1–KLRG1 | CDH1 | barrier\|inhibitory | +0.006 | 0.014 | 0.009 |
| HLA-G–CD8B | MHC-I | inhibitory | +0.004 | 0.004 | 0 |

Recruiting pair **higher** in CLDN4-high (not down):

| Pair | ΔP | P high | P low |
|---|---:|---:|---:|
| CXCL16–CXCR6 | +0.007 | 0.022 | 0.015 |

Largest absolute outgoing deltas are **ECM → CD44**, also higher in CLDN4-high: COL1A1–CD44 (+0.195), LAMB3–CD44 (+0.157), LAMA5–CD44 (+0.117), COL6A1–CD44 (+0.067). MHC-I classical (HLA-A/B/C–CD8A/B) is higher in high (ΔP +0.043 to +0.053). Those rows are in the ligand table; they are not in the barrier/inhibitory class list above.

### Lower in CLDN4-high (significant on the low arm)

| Pair | Pathway | ΔP | P high | P low | p_low |
|---|---|---:|---:|---:|---:|
| SPP1–CD44 | SPP1 | −0.108 | 0.060 | 0.168 | 0.0099 |
| MIF–CD74/CD44 | MIF | −0.035 | 0.735 | 0.769 | 0.0099 |
| MIF–CD74/CXCR4 | MIF | −0.028 | 0.806 | 0.833 | 0.0099 |
| PVR–TIGIT | PVR | −0.003 | 0.0008 | 0.0036 | 0.0099 |
| HLA-DRB1–CD4 / HLA-DRA–CD4 | MHC-II | −0.006 / −0.005 | ~0.007 / 0.013 | ~0.013 / 0.018 | 0.0099 |

### Not detected (do not invent edges)

NECTIN2–TIGIT / NECTIN2–CD226, CD274–PDCD1, CEACAM family, TGFB outgoing, CXCL9/10, CCL5, IL15, MICA: `expr_prop` failed on at least one side. No network is drawn for those pairs.

## Incoming T/NK → Mal (not the ligand table)

21 significant incoming pairs. IFNG–IFNGR1/IFNGR2 is **higher** into CLDN4-high (ΔP +0.0086, p_high = 0.0099). TNFSF10–TNFRSF10B is slightly higher into high. TGFB1–ACVR1B/TGFBR2 is lower into high. Full incoming: `results/contrast_tumor_tertile_incoming.tsv`.

## Sensitivities (same method, different n)

CDH1–ITGAE/ITGB7, JAM1–ITGAL/ITGB2, LGALS9–CD45, HLA-E–CD8A, and CXCL16–CXCR6 stay higher in CLDN4-high in **both** tLung (n=11) and mets (n=21). PVR–TIGIT flips: higher in high on tLung (ΔP +0.008), higher in low on mets (ΔP −0.006). LUNG_T34 dominates the tLung malignant pool.

## What is not claimed

- This is **not** GSE207422 and **not** an ICI / MPR contrast.
- Permutation tests on cell-pooled truncated means are not a 32-patient mixed model. EBUS_28 is 27% of the kept high arm; mBrain is half of author-malignant cells.
- CXCL16–CXCR6 and IFNG incoming are **up**, not down, in CLDN4-high. That is the opposite of an immune-cold / recruit-down story.
- NECTIN2 and PD-L1–PD-1 are not detected here.
- CopyKAT was not re-run. tS1–tS3 are author tumor-epithelial labels, not a new malignant call.

## Reproduce

```bash
python3 methods/gse131907_cellchat_cldn4/scripts/download.py --out /tmp/gse131907
python3 methods/gse131907_cellchat_cldn4/scripts/analyze.py \
  --matrix /tmp/gse131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz \
  --ann /tmp/gse131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz \
  --series /tmp/gse131907/GSE131907_series_matrix.txt.gz \
  --out methods/gse131907_cellchat_cldn4/results
```
