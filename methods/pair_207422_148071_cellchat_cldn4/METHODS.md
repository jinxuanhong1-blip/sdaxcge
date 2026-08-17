# Methods — pairwise merge GSE207422 + GSE148071 (CLDN4 only)

Additive folder. Prior single-dataset CellChat CLDN4 analyses are **taken as given** and are not re-run:

- `methods/scrna_cellchat_cldn4` (GSE207422; PR #324)
- `methods/gse148071_cellchat_cldn4` (GSE148071; PR #348)

**CLDN4 only.** TACSTD2 is not a gate. No dual-high. GSE207422 is included.

Public processed tables only. Raw FASTQ, CellChat R, and LIANA were not run.

## Datasets

| Series | Citation | What is used here |
| --- | --- | --- |
| GSE148071 | Wu et al., *Nat Commun* 2021, PMID 33953163; 42 stage III/IV NSCLC biopsies (Singleron) | Given CellChat per-sample + kept ligand table. Sensitivity: TISCH unit table from the Q4 folder. |
| GSE207422 | Hu et al., *Genome Med* 2023, PMID 36869384; neoadjuvant PD-1 + chemo | Given CellChat **post-treatment** per-sample (12) + kept ligand table. Three pre-treatment samples stay excluded. Sensitivity: author DRMref 12-patient table from the combo folder. |

GSE131907, GSE205335, and GSE253013 are not used.

## Patient-level CLDN4 vs T/NK

Unit = patient / biopsy. Predictor = mean `log1p(CP10k)` CLDN4 in the malignant/epithelial compartment of that sample. Endpoint = T/NK fraction.

**Primary definition (matched to the CellChat folders):**

- Malignant proxy = marker-argmax epithelial lineage (putative; no CopyKAT on these GEO extracts).
- T/NK = T + NK lineages, merged.
- Eligible if ≥25 epithelial cells **and** ≥25 T/NK cells.
- Fraction = `n_TNK / n_total` (all cells in the given per-sample table).
- GSE148071: 25 / 42 eligible. GSE207422: 12 / 12 post (all pass the floor). Combo n = 37.

Spearman ρ per cohort (`scipy.stats.spearmanr`). Combo = DerSimonian–Laird random-effects on Fisher-z(ρ), back-transformed to ρ, with I² and 95% CI. Signed Stouffer (weight √(n−3)) is reported as a second p. p-values are descriptive.

Sensitivity:

1. Same patients, fraction = `n_TNK / (n_epithelial + n_TNK)`.
2. Given TISCH locked eligible (≥20/20) on 148071 + author DRMref malignant mean vs `frac_tnk` on 207422. Mixed definition; not the primary.

## CellChat-style LR (given, then intersected)

Each given folder already computed Jin et al. 2021 Hill probability on CellChatDB v2 protein pairs (10% truncated mean, \(K_h=0.5\), `expr_prop ≥ 0.10`, groups ≥25 cells) and 100 permutations of CLDN4-high/low labels among malignant cells. Kept splits: 148071 **median**; 207422 **median_post**. Those permutations are not repeated.

This folder joins the two kept ligand tables on `interaction_name` + `direction`. A pair is **consensus** if it is present and significant in both tables and `sign(ΔP)` matches. Outgoing = Mal → T/NK. Incoming = T/NK → Mal. No network is drawn for a non-consensus pair.

LGALS9–PTPRC (148071-only) and CD274–PDCD1 (207422-only / opposite story) are reported as **not** consensus.

## Honest n

Deposited n (42 and 15) is not the test n. Combo patient n is eligible 25+12=37. Cell-pooled LR n is cells, not patients; one biopsy can dominate a given pool (P7 in 148071 T/NK; BD_immune07 in 207422 epithelium). That is inherited from the given folders and is not hidden.

## Software

Python (numpy / pandas / scipy / matplotlib). `lib_stats.py` is the same Spearman / DL / Stouffer helper used in the CLDN4 combo folder. Matrices were not re-downloaded.
