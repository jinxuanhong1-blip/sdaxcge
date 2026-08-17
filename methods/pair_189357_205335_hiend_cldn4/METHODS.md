# Methods — pair GSE189357 + GSE205335, CLDN4-only high-end

ADDITIVE. **CLDN4 only.** No dual-high TACSTD2∩CLDN4 score. Patient is the unit.
PR #459 combo Spearman / Q4 vs Q1 is **given** and is not re-ranked.
CellChat R and LIANA are not run.

## Pair (given)

| Cohort | Malignant | T/NK | Eligible n | Source |
| --- | --- | --- | ---: | --- |
| GSE189357 | marker-malignant: (EPCAM\|KRT8\|KRT18\|KRT19)>0 and PTPRC==0 | (CD3D\|CD3E\|CD8A\|NKG7\|GNLY\|KLRD1)>0 and not malignant | 9 | Zhu/Wang AIS–IAC 10x MTX; PR #459 |
| GSE205335 | author `Malignant cells` | author `T/NK cells` | 22 | Hu et al. UMI RDS + identity; PR #279 / #320 / #459 |

Given combo (PR #459, %pos): **n=31, ρ=−0.478 (p=0.009, I²=0%), Q4 vs Q1 r=−0.750 (p=0.010; 8/8)**.

Malignant definitions are **not the same**. Quartiles for the given row are within-cohort.

## High-end split (this PR)

Score = `log1p(CP10k)` CLDN4 on that patient's malignant cells.

| Split | High | Low |
| --- | --- | --- |
| Median | ≥ median (if median is 0: CLDN4>0 vs =0) | < median (or =0) |
| Q4 vs Q1 | `pd.qcut` on average ranks, Q4 tail | Q1 tail |

**Honest paired n.** A patient is scored only if:

- malignant CLDN4-high bin ≥ 10 cells
- malignant CLDN4-low bin ≥ 10 cells
- same-patient T/NK ≥ 20 cells

Patients missing either bin or T/NK are out. Cells are not n. The given n=31 is the combo Spearman n, not the communication n.

## CellChat-style probability

Jin et al. 2021 on **CellChatDB v2** protein pairs (Secreted / Cell–Cell Contact / ECM–Receptor; non-protein dropped). A pair is kept only if every ligand and receptor subunit is present in **both** matrices.

1. Per cell group: 10% truncated mean of `log1p(CP10k)`; complexes = geometric mean of subunits (0 if any subunit mean is 0).
2. Detected if each complex has expressing-cell fraction ≥ 0.10 (AND).
3. \(P = (L \cdot R) / (K_h + L \cdot R)\) with \(K_h = 0.5\).

Outgoing = Mal_high or Mal_low → that patient's T/NK. Incoming is not the primary table.

## Tests

**Primary.** Wilcoxon signed-rank on per-patient \(P_{\mathrm{high}}\) vs \(P_{\mathrm{low}}\) for pairs detected on **both** arms. Minimum paired n after the detect gate = 6. Thin n (<8) is flagged. p-values are descriptive.

**Companion.** All-malignant → same-patient T/NK on the locked 31-patient combo. Within-cohort Q4 vs Q1 (and median) of malignant CLDN4 %pos, then stacked. Mann–Whitney on per-patient *P* (detected ≥3 per arm). This does **not** re-rank the given Spearman.

Focused pairs are always written if subunits exist: classical MHC-I (HLA-A/B/C–CD8), T-recruit (CXCL9/10/16, CCL4/5), CD274–PDCD1, NECTIN2–TIGIT.

## Software

Python (numpy / pandas / scipy / matplotlib / rdata). CellChatDB v2 from `jinworks/CellChat`. No R CellChat runtime. No FASTQ.
