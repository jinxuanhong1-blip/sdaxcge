# FINDING — pairwise CLDN4-only combo GSE127465 + GSE148071

**Additive. CLDN4 only. Not a triple/quad. No dual-high.** Patient is the unit.
Public processed matrices are both **<2 GB** and both carry **CLDN4** plus malignant and T/NK lineages. The pair was therefore scored (not stopped).

This is **not** a cell-level merge. GSE127465 is author-normalized inDrops (Zilionis et al., *Immunity* 2019, PMID 30979687). GSE148071 is 10x UMI per-sample TXT (Wu et al., *Nat Commun* 2021, PMID 33953163). CellChat-style outgoing is run **per cohort** only if the two patient-level CLDN4 vs T/NK associations differ.

Numbers below the gate are written by `scripts/analyze.py` into `results/summary.json` and copied here after the run.

---

## Gate

| Series | Processed expression | Bytes | <2 GB | CLDN4 | Malignant / epithelium | T/NK |
|---|---|---:|:---:|:---:|:---:|:---:|
| [GSE127465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465) | `GSE127465_human_counts_normalized_54773x41861.mtx.gz` | 528,303,938 | yes | present | author `PatientN-specific` (tumor) | author `tT cells` + `tNK cells` |
| [GSE148071](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148071) | `GSE148071_RAW.tar` (42 `*_exp.txt.gz`) | 180,193,280 | yes | present | no GEO label; marker-argmax epithelium (putative) | marker-argmax T + NK |

Mouse MTX on GSE127465 (166.2 Mb) was not used. Blood cells on GSE127465 were not used for malignant CLDN4 or T/NK. GSE127465 `RAW.tar` (168.0 Mb) is per-library TSV and was not needed once the series MTX + metadata were present.

Stop rule (*if matrix missing, stop*) does **not** apply.

---

## Honest n (deposited / lineage)

Do not write n = 54,773 or n = 89,887 as the test n.

| Item | n | Note |
|---|---:|---|
| GSE127465 human patients | **7** | p1–p7 |
| GSE127465 tumor / blood cells | 40,362 / 14,411 | blood excluded |
| GSE127465 author own-patient `PatientN-specific` in tumor | 3,574 | primary malignant (p1 1,289; p2 98; p3 527; p4 172; p5 757; p6 214; p7 517) |
| GSE127465 tumor T / NK / T+NK | 12,776 / 1,116 / 13,892 | all 7 patients have ≥20 malignant and ≥20 T/NK |
| GSE148071 patients deposited | **42** | one biopsy each |
| Dual-high TACSTD2×CLDN4 | **0** | not defined |
| Extra series in this pair | **0** | not a triple/quad |

Eligible n and Spearman ρ are in `results/combo_cldn4_tnk.tsv` after `analyze.py`.

---

## Methods

- Predictor: malignant **CLDN4 only** (mean log1p). TACSTD2 is not a gate.
- GSE127465 malignant = `Major cell type == Patient{N}-specific` inside that patient's **tumor**. Normal Type I/II, club, ciliated, fibroblasts, and endothelium are not malignant.
- GSE127465 T/NK = tumor `tT cells` + `tNK cells`. Fraction = T/NK / tumor cells.
- GSE148071 epithelium / T / NK = marker-argmax on log1p(CP10k) modules (same marker list as `methods/gse148071_cellchat_cldn4`). Epithelium is **putative** malignant.
- Eligible patient: ≥20 scored malignant/epithelial cells **and** ≥20 T/NK.
- Association: Spearman ρ, patient-level. Fisher-z 95% CI. Two-cohort difference = opposite sign (with |ρ| threshold) or Fisher-z two-sample p < 0.05 or one |ρ| ≥ 0.30 and the other |ρ| < 0.15.
- If they differ: CellChat-style **outgoing** only (Mal CLDN4-high vs low → T/NK). Hill / mass-action on 10% truncated means, `expr_prop ≥ 0.10`, CellChatDB v2 protein pairs, 100 permutations of the high/low label among malignant cells. CellChat R is not installed. Cohorts are not stacked into one communication object.
- p-values are descriptive on n=7.

---

## Patient-level malignant CLDN4 vs T/NK

*Filled after `python3 methods/pair_127465_148071_cldn4/scripts/analyze.py`.*

| cohort | assay | n deposited | n eligible | ρ [95% CI] | p | note |
|---|---|---:|---:|---|---:|---|
| GSE127465 | inDrops, author malignant | 7 | — | — | — | pending run |
| GSE148071 | 10x, putative epithelium | 42 | — | — | — | pending run |

Differ? — pending run. CellChat-style outgoing: only if they differ.

---

## What this is not

- Not a triple or quad (no GSE131907 / GSE205335 / GSE207422).
- Not dual-high TACSTD2×CLDN4.
- Not a merged Seurat/AnnData object.
- Not CellChat R / LIANA / NicheNet packages.
- Not ICI response (neither deposit has a usable RECIST/MPR table in the processed files used here).

Reproduce: `python3 methods/pair_127465_148071_cldn4/scripts/analyze.py`
