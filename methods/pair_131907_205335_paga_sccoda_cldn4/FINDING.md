# FINDING — pair GSE131907 + GSE205335 PAGA + scCODA, CLDN4 only

ADDITIVE **CLDN4-only** high-end on the combo that already differs
(PR #320: author-malignant CLDN4 %pos vs T/NK, Q4 vs Q1 **n=23**
r=−0.705). This run does **not** re-audit that combo, add GSE148071,
or re-run the 7-cohort / 6-unit pool. No TACSTD2∩CLDN4 dual-high gate.

A previous attempt of this exact task OOM'd mid-Harmony/UMAP on
**53,296** unsampled malignant cells. This retry **subsamples per unit**
(malignant cap 150, T/NK cap 80) and runs **scCODA first** from the
PR #320 tables so composition results do not depend on the embedding.

**Composition holds; trajectory DPT does not.** On the same Q4 vs Q1
units (n=23 = 12/11), T/NK ALR is −2.09 (perm p=0.002, q=0.012) and
B ALR is −2.19 (p=0.006, q=0.018). Five of six primary tests recover
at q<0.10; GSE205335 B is the miss. Unit-mean CLDN4 vs DPT is null
(n=43, ρ=−0.065, p=0.68). CLDN4 still tracks a malignant-like score
(ρ=0.587, q=2.5×10⁻⁴) and, more weakly, barrier/keratin without CLDN4
(ρ=0.328, q=0.074).

## Honest n

Unit = GSE131907 **sample** (not patient) + GSE205335 **patient**.
Cells are library size / embedding weight only.

| cohort | unit_type | n_eligible | n_q1 | n_q4 | n_q4q1_compared | note |
|---|---|---|---|---|---|---|
| GSE131907 | sample | 21 | 6 | 5 | 11 | GSE131907 sample is the unit (not patient) |
| GSE205335 | patient | 22 | 6 | 6 | 12 | GSE205335 B = B+plasma; all 22 patients |
| merge_within_cohort | sample+patient | 43 | 12 | 11 | 23 | PR #320 Q4 vs Q1 n=23; not a 7-cohort pool |

Merged Q4 vs Q1 compared n = **23** (PR #320 12/11). Do not write n=43
patients — GSE131907 contributes 21 **samples**.

Catalog before cap: GSE131907 24,784 malignant + 15,150 T/NK; GSE205335
28,512 malignant + 31,642 T/NK (53,296 malignant if merged unsampled —
the OOM object). After per-unit cap: **6,134 malignant + 3,382 T/NK =
9,516 cells**, 43 units. Harmony ran on that subsample (True). DPT root:
min CLDN4 among the top AT2 tercile (never CLDN4-high; not PR #325
nLung AT2).

## scCODA (primary = Q4 vs Q1)

Engine: ALR of T/NK or B vs Other + unit-level permutation p.
scCODA HMC: `not_available`. Not faked.
Primary family = 6 tests (T/NK and B × GSE131907 / GSE205335 / merge).
Recovery = ALR effect < 0 and BH q < 0.10. Recoveries: **5**.

Merge T/NK Q4 vs Q1: ALR -2.086, perm p=0.002, q=0.012, n_low/n_high=12/11.

### Primary table

| slice | compartment | n | n_low | n_high | alr_effect | alr_p_perm | q_bh_primary | frac_delta_high_minus_low | frac_mwu_p | spearman_rho | recover_q10 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GSE131907 | TNK | 11 | 6 | 5 | -2.275 | 0.026 | 0.039 | -0.342 | 0.126 | -0.327 | True |
| GSE131907 | B | 11 | 6 | 5 | -3.387 | 0.017 | 0.035 | -0.079 | 0.017 | -0.791 | True |
| GSE205335 | TNK | 12 | 6 | 6 | -1.961 | 0.039 | 0.047 | -0.398 | 0.026 | -0.608 | True |
| GSE205335 | B | 12 | 6 | 6 | -1.171 | 0.258 | 0.258 | -0.028 | 0.818 | 0.014 | False |
| merge_within_cohort | TNK | 23 | 12 | 11 | -2.086 | 0.002 | 0.012 | -0.304 | 0.005 | -0.501 | True |
| merge_within_cohort | B | 23 | 12 | 11 | -2.194 | 0.006 | 0.018 | -0.054 | 0.029 | -0.546 | True |

Full tests (Q4 vs Q1 + median sensitivity): [`results/tables/sccoda_tests.tsv`](results/tables/sccoda_tests.tsv).

## PAGA / DPT (malignant, subsampled)

Not a redo of PR #325 (GSE131907 nLung+tLung epithelium, AT2-rooted). This object is author-malignant cells from the winning pair only (GSE131907 metastases with ≥20 malignant cells + GSE205335 tumor libraries). Barrier/keratin excludes CLDN4. TACSTD2 is a comparator Spearman only.

| slice | contrast | n | rho | p | q_bh_merge |
|---|---|---|---|---|---|
| merge | CLDN4 vs DPT | 43 | -0.065 | 0.681 | 0.681 |
| merge | CLDN4 vs AT2 | 43 | 0.252 | 0.103 | 0.152 |
| merge | CLDN4 vs barrier/keratin | 43 | 0.328 | 0.032 | 0.074 |
| merge | CLDN4 vs club | 43 | 0.248 | 0.108 | 0.152 |
| merge | CLDN4 vs basal | 43 | 0.189 | 0.225 | 0.262 |
| merge | CLDN4 vs malignant-like | 43 | 0.587 | 3.56e-05 | 2.49e-04 |
| merge | CLDN4 vs TACSTD2 (comparator) | 43 | 0.404 | 0.007 | 0.026 |

Paired within-unit CLDN4 tertiles (extra; 23 units with ≥6 cells/arm
after the cap): barrier/keratin Δmed=+0.210, p=7.2×10⁻⁷. DPT Δmed=+0.003
(p=1.5×10⁻⁴) is a tiny consistent shift, not a trajectory effect.
AT2 paired Δmed=0, p=0.069.

Per-cohort Spearman rows, Leiden vertices, and T/NK PAGA connectivities:
`results/tables/`.

GSE131907 T/NK fraction MWU on the Q4 vs Q1 tails is p=0.126 — the
companion is weaker than the ALR primary on that slice. That is stated.

## Extra figures

- `results/figures/fig_malignant_trajectory_cldn4.png` — UMAP CLDN4 / DPT / cohort / unit quartile
- `results/figures/fig_malignant_paga.png` — PAGA graph
- `results/figures/fig_malignant_extra_programs.png` — unit-mean CLDN4 vs DPT / barrier / AT2
- `results/figures/fig_tnk_trajectory_cldn4.png` — T/NK UMAP (no AT2 DPT)
- `results/figures/fig_sccoda_fractions.png` — T/NK and B fractions, Q4 vs Q1
- `results/figures/fig_sccoda_alr_forest.png` — primary ALR forest
- `results/figures/fig_honest_n.png` — unit counts

## What this is not

- Not dual-high TACSTD2×CLDN4.
- Not GSE148071.
- Not the 7-cohort / 6-unit malignant pool.
- Not PR #325 GSE131907-only epithelial PAGA and not PR #339 multi-cohort ICI scCODA.

## Reproduce

```bash
python3 methods/pair_131907_205335_paga_sccoda_cldn4/scripts/run_sccoda.py
python3 methods/pair_131907_205335_paga_sccoda_cldn4/scripts/download.py
python3 methods/pair_131907_205335_paga_sccoda_cldn4/scripts/extract.py
python3 methods/pair_131907_205335_paga_sccoda_cldn4/scripts/run_paga.py \
  --input /tmp/pair_131907_205335/pair_subsample.h5ad
python3 methods/pair_131907_205335_paga_sccoda_cldn4/scripts/write_finding.py
```
