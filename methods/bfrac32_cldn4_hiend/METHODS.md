# Methods — CLDN4-only high-end LR on the given B-frac combo

## Given (not re-audited)

PR #290 `methods/scrna_cldn4_combo`, family `author/b/mean`:

| combo | k | N | CLDN4 ρ | p | I² |
|---|---:|---:|---:|---:|---:|
| GSE131907 + GSE241934 IIT | 2 | 32 | −0.513 | 0.00388 | 0% |

Members: GSE131907 author-malignant n=21; GSE241934 IIT author-epi n=11.

This slice does **not** recompute that Spearman. It does **not** redo GSE207422-only B-frac (PR #400, NS). Dual-high TACSTD2∩CLDN4 is not used.

## Cohorts (from existing methods/notes)

Public processed GEO only. Named in `methods/scrna_tls_meta/playbook.md` and `methods/scrna_cldn4_combo/FINDING.md`.

- **GSE131907** (Kim et al., *Nat Commun* 2020, PMID 32385277). Author labels. Malignant = `Cell_subtype == Malignant cells` on the 21 tL/B + mLN + mBrain samples that already enter the combo. tLung tS1–tS3 is a different row and is not added. B = `B lymphocytes`. T/NK = `T lymphocytes` + `NK cells`.
- **GSE241934 IIT** (NEOTIDE). Author Epi as the malignant proxy (combo `author_epi`). B / T / NK from author major type. Real/RWC is a different stratum and is not downloaded.

Matrices are not concatenated.

## Estimand

Within each patient, outgoing ligand–receptor score from **CLDN4-high vs CLDN4-low** malignant (or IIT Epi) cells toward:

1. B cells
2. TLS-like cells = B + CXCL13+ T/NK (dissociated proxy, not a follicle)
3. T/NK

Patient is the unit. Honest n = patients who pass floors, which may be < 32.

## Split and floors

- Score = `log1p(UMI / total × 10⁴)` CLDN4. **CLDN4 only.**
- Primary: within-patient **median**. Extra high-end: within-patient **tertile** (middle third dropped).
- Floors: ≥10 high and ≥10 low malignant; ≥10 B; ≥10 TLS-like; ≥20 T/NK (same floors as `methods/gse131907_liana_cldn4`).

## Scores

**CellChat-style** (Jin et al. 2021; CellChat R not run): 10% truncated mean, Hill \(K_h=0.5\), CellChatDB v2 protein pairs, `expr_prop ≥ 0.10`. Complexes = geometric mean of subunits.

**LIANA-style** (Efremova 2020 / Garcia-Alonso 2022; LIANA package optional): CellPhoneDB mean-of-means on log1p(CP10k), CellPhoneDB v5 pairs, `expr_prop ≥ 0.10`. Complexes = min of subunits.

Inference: paired Wilcoxon on patient high vs low scores. BH-FDR inside each method × destination × rule. Cells are not replicates.

## Cannot test

Histologic TLS, CopyKAT malignancy, GSE131907 ICI/MPR (none), GSE207422 B-frac, dual-high gates.
