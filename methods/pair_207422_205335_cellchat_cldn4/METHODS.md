# Methods — pairwise GSE207422 + GSE205335 CLDN4 CellChat / LIANA

ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. **Include 207422.**
Patient is the unit. Prior TACSTD2 A3 and the multi-cohort Q4 meta are given.

## Datasets

- [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422) — Hu et al., *Genome Medicine* 2023, neoadjuvant PD-1 NSCLC scRNA. Public UMI (24,292 × 92,330). Author CopyKAT IDs are not public. Malignant = DRMref `Malignant cells` (Liu et al., *NAR* 2024); T/NK = CD8+ T + CD4+ T + NK. Locked n=12 post-treatment patients.
- [GSE205335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335) — Hu et al., palliative ICI biopsy/effusion scRNA. Author `lineage.sub == Malignant cells` and `lineage.total == T/NK cells`. Locked n=22 with ≥20 cells in each compartment.

## Combo rho and Q4 vs Q1

Within each cohort, rank malignant CLDN4 (mean log1p(CP10k) is the locked aligned score; %pos is reported where available) and cut quartiles with `pd.qcut` on average ranks. Spearman ρ vs same-patient T/NK fraction (`n_tnk / n_cells`). Q4 vs Q1 is two-sided Mann–Whitney on that fraction; rank-biserial \(r = 2U/(n_4 n_1)-1\).

GSE207422 %pos is computed from the public UMI on DRMref malignant cells. It is **not** on the locked A3 TACSTD2 table. GSE207422 mean vs %pos quartile membership differs (P06 is %pos Q4 with 15 malignant cells and is dropped from LR). Ligand tables therefore use **mean** quartiles so all 9 vs 9 tails are LR-eligible.

**Combo** is DerSimonian–Laird random-effects on Fisher-z(ρ) (and on *r* for the tails). Quartiles are **within-cohort**. A stacked 34-patient Spearman is not the primary estimator.

## CellChat-style ligands

Jin et al. 2021 Hill probability on CellChatDB v2 protein pairs (same parse as `methods/scrna_cellchat_cldn4`):

1. Keep a pair only if every subunit is in the UMI matrix.
2. Per compartment: 10% truncated mean of `log1p(CP10k)`; complexes = geometric mean (0 if any subunit mean is 0).
3. Detected if each complex has expressing-cell fraction ≥ 0.10.
4. \(P = (L \cdot R)/(K_h + L \cdot R)\) with \(K_h = 0.5\).
5. A patient is scored only if it has ≥20 malignant and ≥20 T/NK cells.

Outgoing = malignant → **same-patient** T/NK. Q4 vs Q1 is Mann–Whitney on per-patient *P* (detect gate ≥3 per tail). GSE205335 *P* is reused from the given per-patient table, remapped onto **mean** quartiles. GSE207422 is scored from the public UMI + DRMref labels.

## LIANA-style ligands

Primary combo table: CellPhoneDB-style **mean of partner means** on the same CellChat complexes (comparable across cohorts). GSE207422 also writes a native CellPhoneDB v5 table (min-subunit, arithmetic mean) as a companion. LIANA `mt.cellphonedb` permutations and CellChat R are not run.

## Honest n

- Unit = patient, not cell.
- GSE207422 Q4 vs Q1 is 3 vs 3. Combo tails are 9 vs 9.
- P06 (15 DRMref-malignant cells) stays on the locked mean table and is dropped from LR.
- GSE205335 Q4 mixes SCLC with ADC; that mix is reported.
