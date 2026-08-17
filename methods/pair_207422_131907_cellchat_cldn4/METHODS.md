# Methods — Pair GSE207422 + GSE131907, CLDN4-only

Additive pairwise merge. **CLDN4 only.** No dual-high TACSTD2×CLDN4 score.
GSE207422 is included. Prior single-cohort folders are taken as given and
are not re-run. GEO UMI matrices are not re-downloaded.

## 1. Patient-level malignant CLDN4 vs T/NK

| Cohort | Malignant | Immune | Score | Unit | Source |
| --- | --- | --- | --- | --- | --- |
| GSE207422 | DRMref `Malignant cells` on the locked A3 12-patient table | T+NK fraction | mean log1p(CP10k) CLDN4 | patient | `data/tnk/GSE207422_drmref_patients.tsv` |
| GSE131907 | author `Malignant cells`, `n_malignant ≥ 20` | T+NK fraction | mean log1p(CP10k) CLDN4 | sample | `data/tnk/GSE131907_samples.tsv` |

GSE131907 primary tLung epithelium is labeled tS1–tS3, not `Malignant cells`.
Those samples drop out of the primary T/NK row (`n_malignant = 0`). That is
stated, not patched with all-epithelial CLDN4. A tLung epithelial sensitivity
is computed and labeled as not primary.

Within each cohort:

- Spearman ρ of malignant CLDN4 vs T/NK.
- Quartiles of CLDN4 (average ranks, `qcut`); Mann–Whitney U on T/NK in Q4 vs Q1;
  rank-biserial \(r = 2U/(n_4 n_1) - 1\). \(r < 0\) means Q4 is colder.

Pair combo: DerSimonian–Laird random effects on Fisher-z of ρ (or of *r*),
back-transformed. Stouffer weights \(\sqrt{n-3}\). Quartile tails are **not**
pooled across cohorts before ranking — each study is cut on its own CLDN4
distribution. Honest n = 12 + n_mal≥20 samples for Spearman; n_Q1 + n_Q4
summed across studies for Q4 vs Q1.

p-values are descriptive.

## 2. CellChat-style outgoing CLDN4-high → T/NK

CellChat R and LIANA are not installed. Each cohort already has Jin et al.
2021 Hill probability on CellChatDB v2 protein pairs (Secreted / Cell–Cell
Contact / ECM–Receptor) plus 100 permutations of CLDN4-high/low among
malignant cells. A pair is significant if detected (`expr_prop ≥ 0.10` both
sides), \(P > 0\), and permutation \(p < 0.05\). Smallest p = 1/101 = 0.0099.

Kept splits (not re-chosen here):

- GSE207422: all-post **median** CLDN4 on marker epithelium (12 patients).
- GSE131907: tumor-origin **tertile** on author malignant + tS1–tS3 (32 samples).

The combo LR table is an outer join on `interaction_name`. Concordance:

- `both_sig_same` / `both_sig_discordant`
- `only_207422` / `only_131907` (other cohort detected but not significant)
- `only_*_undetected_*` (other cohort failed detection)
- `both_detected_ns`

Outgoing = Mal → T/NK. Incoming contrasts are stored but are not the ligand
table. Cells are not pooled across studies (different chemistry / labels).

## Software

Python (numpy / pandas / scipy / matplotlib). Same `lib_stats.py` as
`methods/cldn4_malig_q4_tnk`.
