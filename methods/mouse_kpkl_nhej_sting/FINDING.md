# Public mouse lung KP/KL beyond GSE137244: Cldn4, NHEJ, STING

Public matrices only. Private 8 KL mice were not read and were not merged with GSE165641, GSE180963, GSE154977, or anything else. GSE137244 stays the locked cell-line baseline and was not recomputed: KL n=5 vs KP n=5, Cldn4 Δ=+5.57, Mann–Whitney p=0.00794.

## What was asked

Outside that cell-line series, do open mouse lung KP (Kras/Trp53) and KL (Kras/Lkb1) datasets still show higher Cldn4 in KL, a lower NHEJ module, and a lower STING module?

## Modules

- **Cldn4** and **Sting1** (alias Tmem173) are single genes.
- **NHEJ** is the mean within-contrast z-score of Xrcc5, Xrcc6, Prkdc, Dclre1c, Lig4, Xrcc4, Nhej1, Poll, Polm, Pnkp, Aptx, Paxx, Aplf. Paxx is absent on several older platforms (12/13 genes used).
- **STING_core** is Cgas (Mb21d1), Sting1, Tbk1, Irf3, Ikbke. Downstream ISGs are a separate score and are not called STING.
- **Stk11** is a label check, not an endpoint.

Test: two-sided Mann–Whitney on the pre-declared contrast. Benjamini–Hochberg FDR is within the in-vivo KL-vs-KP family and, separately, within the cell-line / isogenic KL-vs-KP family, across Cldn4, Sting1, NHEJ, and STING_core. n<4 is marked underpowered: the exact two-sided test cannot reach p<0.05. Datasets are not pooled.

## Inventory

Scored, with a real KP and KL (or isogenic Lkb1-loss) arm: GSE137396, GSE274351, GSE274352 (empty vector), GSE164758, GSE6135, GSE69552, GSE244452, GSE322570. Scored as isogenic Stk11 KO, not as a KP GEMM: GSE338923. Scored as KL vs Kras-only or vs Kras/Keap1, not as KP: GSE175479, GSE21581, GSE133895, GSE193895.

Open but not scored: GSE133714 (no KP; no standalone annotation), GSE165640 (compare tables only), GSE118246 (Lkb1 models are not Kras;Lkb1), GSE277929 and GSE253613 and GSE182228 (KL only), GSE194166 (CD45-sorted), GSE180963 / GSE165641 / GSE154977 (scRNA; not merged with private 8 KL). Full table: `results/mouse_kpkl_nhej_sting/tables/inventory.tsv`.

## Cldn4 is not higher in KL than KP once GSE137244 is left behind

Where Stk11 mRNA actually marks the KL arm, Cldn4 does not.

| Series | What | n KL vs KP | Stk11 Δ (p) | Cldn4 Δ | Cldn4 p | FDR |
|---|---|---:|---:|---:|---:|---:|
| GSE137396 | autochthonous nodules | 5 vs 5 | −1.69 (0.0079) | +0.43 | 1.00 | 1 |
| GSE6135 | Ji primary tumors, mouse-level | 7 vs 5 | −0.34 (0.010) | +0.22 | 0.88 | 0.91 |
| GSE69552 | KL papillary + adenosquamous vs KP | 12 vs 5 | −0.94 (3.2×10⁻⁴) | +0.20 | 0.28 | 0.51 |
| GSE69552 | papillary only | 6 vs 5 | −1.02 (0.0043) | −0.09 | 0.93 | — |
| GSE274351 | LCM adenomas | 5 vs 5 | −2.38 (0.20) | −1.04 | 0.67 | 0.87 |
| GSE164758 | primary GEMM tumors | 9 vs 8 | −0.06 (0.89) | −0.35 | 0.42 | 0.73 |

GSE6135 Cldn4 Δ=+0.22 is the same mouse-level contrast already reported for this array. GSE137396 is the in-vivo nodule series from the same study family as GSE137244; Stk11 separates cleanly and Cldn4 does not.

GSE164758 is the largest primary-tumor RNA-seq, and its KL columns are not Stk11-low in bulk RNA (Δ=−0.06). Cldn4 there is slightly lower in the KL-labeled columns, not higher. KL vs Kras-only in that same file is Cldn4 Δ=+1.50 (p=0.006), but Stk11 is not lower in those KL columns versus Kras-only either, so that rise is not a confirmed Lkb1-loss effect.

Two underpowered KL-high directions, neither able to reach p<0.05:

- GSE244452 subcutaneous syngeneic tumors, n=3 vs 3. Stk11 is lower (complete separation). Cldn4 Δ=+8.98 log2. Mann–Whitney p=0.10, which is the floor. Actb is absent from the table. This is the only public series here that looks like the GSE137244 Cldn4 gap, and it cannot be tested.
- GSE274352 empty-vector cell lines, n=3 vs 3. Stk11 lower. Cldn4 Δ=−0.20, p=1. Not the cell-line repeat of GSE137244.

KL versus Kras-only, which is not a KP contrast: GSE175479 bulk lung, n=4 vs 4, Stk11 Δ=−1.30 (p=0.029), Cldn4 Δ=+1.68 (p=0.029). That delta matches the earlier Welch result on this series. GSE21581 primary tumors, n=9 vs 9, Stk11 p=0.0011, Cldn4 Δ=+0.55, p=0.052.

## NHEJ is not a KL-versus-KP result

No in-vivo KL-vs-KP NHEJ contrast passes FDR.

The closest is GSE164758, module Δ=−0.35, p=0.021, FDR=0.17, on the series whose Stk11 mRNA does not confirm the label. The drop is carried by Xrcc4, Xrcc5, Xrcc6 (each p=8.2×10⁻⁵) and Prkdc (p=0.0025). Those single-gene tests were not in the FDR family.

Stk11-confirmed GSE137396 nodules: NHEJ p=0.22. GSE6135 vs KP: p=0.53. GSE69552: p=0.19. GSE274351: p=0.69.

Against Kras-only, not KP, the module is lower in GSE6135 (p=0.018) and GSE21581 (p=0.034). It is not lower in GSE175479 (point estimate higher, p=0.20).

## STING does not come out as Sting1 silencing in KL versus KP

Sting1 itself is not lower in KL than KP in any scored series at FDR<0.05.

The one in-vivo FDR hit on Sting1 is in the other direction. GSE164758 KL-labeled bulk tumors have higher Sting1 than KP (Δ=+0.62, p=0.0037, FDR=0.044). Those columns are not Stk11-low, and the matrix is bulk tumor, so this is not tumor-cell STING.

GSE69552 STING_core is higher in KL (Δ=+0.97 z, p=3.2×10⁻⁴, FDR=0.0078). Sting1 in the same contrast is lower (Δ=−0.39, p=0.082). The module rise is Tbk1 (p=3.2×10⁻⁴), Irf3 (p=6.5×10⁻⁴), and Ikbke (p=3.2×10⁻⁴). STING_ISG is not up (Δ=−0.34, p=0.38). Papillary-only STING_core remains higher (p=0.030) and Sting1 remains lower (p=0.13). This is not a Sting1 result.

GSE137396, the clean nodule contrast: Sting1 Δ=−0.06, p=0.69. GSE6135: Δ=−0.18, p=0.88. GSE274351 adenomas: Sting1 Δ=−2.87, p=0.12, and the ISG score is lower (p=0.032); neither is an FDR call.

Sorted tumor cells, GSE322570 sgLkb1 vs sgNeo on KP, n=3: Sting1 Δ=−0.10, p=0.70. Cldn4 Δ=+0.25, p=0.70.

Lacun3 in vitro, GSE338923, parental vs Stk11 KO, n=4 vs 4. Stk11 Δ=−0.99, p=0.029 (complete separation; that p is the floor). Sting1 Δ=−0.34 and Cldn4 Δ=−0.46, both p=0.029. Actb and Gapdh also separate at the same floor, and Cldn4 is near zero in both arms (0.08 vs 0.54 log2 CPM). This is the only pure tumor-cell Stk11 deletion here. It does not reproduce KL-high Cldn4, and it is not specific to STING.

## Cldn4 does not track the modules inside a contrast

Spearman of Cldn4 vs NHEJ, STING_core, or Sting1 inside each KL-vs-KP contrast: one FDR<0.05 hit. GSE244452 Cldn4 vs NHEJ ρ=−1, n=6, FDR=0. That is the same 3-vs-3 separation (KL high Cldn4, low NHEJ score), not a second cohort. On those same six samples, Cldn4 vs STING_core is ρ=+0.94, FDR=0.077. Every other within-contrast correlation has FDR≥0.30.

## Putting STING back does not move Cldn4 in a testable way

GSE274352, STING-V154M versus empty, paired inside the STING file (not pasted onto the IFN-β file; those two matrices are normalized separately):

- Cldn4 up in 4/4 lines, median Δ=+0.25, Wilcoxon p=0.125. That is the smallest two-sided p at n=4.
- Sting1 up in 3/4 (KL1 Δ=−0.05), median Δ=+0.53, p=0.25.
- NHEJ median Δ=−0.05, p=0.63.

IFN-β, 6 pairs: Sting1 up in 6/6, median Δ=+0.74, p=0.031. Cldn4 up in 4/6, p=0.56.

## Bottom line

Public mouse lung beyond GSE137244 does not carry a KL>KP Cldn4 result. The Stk11-confirmed nodule, microarray, and histology contrasts are flat. Cldn4 is higher for KL than Kras-only in GSE175479, and directionally huge in one n=3 syngeneic series that cannot be tested. NHEJ is not an FDR-supported KL-versus-KP module. Sting1 is not silenced in KL versus KP in these matrices; the FDR hits are a higher bulk Sting1 in GSE164758 (Stk11 mRNA does not confirm that label) and a STING_core score in GSE69552 that is Tbk1/Irf3/Ikbke, with Sting1 itself down. Cldn4 does not correlate with either module once that n=6 syngeneic separation is set aside. Private 8 KL mice were not used.
