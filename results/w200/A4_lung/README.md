# W200-A4 · TISMO NSCLC / lung-only — Tacstd2 and Cldn4 paired ICB

**Task:** claim A4 (TISMO `Tacstd2` up in 49/64 ICI models, p=5.8e-5), restricted
to **lung-carcinoma / NSCLC models**, now for both **Tacstd2** and **Cldn4**.

## Verdict

**Not supported as a lung/NSCLC result.**

TISMO annotates two lung-carcinoma lines: **LLC** and **CMT-167**. The only
paired ICB design is **GSE155972 LLC** (subcutaneous flank; anti-PD1 +
anti-CTLA4; WT vs Setdb1_KO). CMT-167 has no ICB arm.

In that one study both genes sit near the detection floor (mean TPM ≈ 0.1–0.2).
One uncorrected Mann–Whitney U (Tacstd2 WT ICB vs baseline, p=0.045) is
nominally <0.05; it does **not** survive Benjamini–Hochberg FDR across the four
primary tests (FDR=0.18). Welch t is n.s. (p=0.15). Cldn4 is n.s. in both arms.
The Setdb1_KO mean-up is one sample (**SRX8918393**). Immune controls
(Cd8a / Cd274 / Gzmb) **do** rise, so the design can detect ICB effects.

The A4 **49/64** count is a **pan-TISMO mean-direction** tally (49 up / 15 down
of 64 ICB treated-vs-baseline cohorts). Only **2/64** cohorts are lung. Cldn4
has no consistent up-direction even pan-TISMO (34 up / 29 down / 2 tie).

## What A4 actually is

User claim (PPT): mouse `Tacstd2` is up in 49 of 64 TISMO ICI comparisons,
p=5.8×10⁻⁵.

Reproduced here from the live TISMO gene-module export
(`icbList` = all six ICB treatments, `tumorList=All`):

| Rule | Gene | UP | DOWN | TIE | n | two-sided binomial p | matches 49/64? |
|---|---|---:|---:|---:|---:|---:|---|
| mean | Tacstd2 | 49 | 15 | 0 | 64 | 2.44e-5 | yes |
| median | Tacstd2 | 40 | 18 | 6 | 64 | 0.060 | no |
| mean | Cldn4 | 34 | 29 | 2 | 65 | 0.80 | no |
| median | Cldn4 | 35 | 22 | 8 | 65 | 0.62 | no |

- Exact two-sided binomial for 49/64 vs 0.5 is **2.44e-5**; one-sided greater is
  **1.22e-5**. The user-reported **5.8e-5 was not recovered**.
- Wilcoxon signed-rank on the 64 Tacstd2 mean-deltas (greater): p=2.92e-5.
- Two 4T1 GSE130472 pairs share baseline samples, so the 64 tests are not fully
  independent.
- Cldn4 has 65 cohorts because TISMO’s Cldn4 export splits YTN16 n-labels
  differently and adds MOC22.

**Lung slice of the 64:** 2 cohorts, both LLC GSE155972. Mean-direction is UP
in both, but that is the near-floor / outlier pattern below — not a 49/64-style
result.

## Lung inventory (TISMO `cancerType == Lung carcinoma`)

19 design groups, 8 studies, 2 cell lines.

| Study | Line | ICB? | Note |
|---|---|---|---|
| GSE155972 | LLC WT / Setdb1_KO | **yes** — antiPD1+antiCTLA4 | only paired ICB |
| GSE100412 | CMT-167, LLC | no | orthotopic, untreated |
| GSE131271 | LLC | no | orthotopic, untreated |
| GSE115109, GSE71491, GSE80678, GSE84535, GSE98672, GSE148101 | LLC | no | chemo / vehicle / genotype / metastasis |

KPB25L / p53-2225L / p53-2336R appear in the ICB gene export but TISMO labels
them **Mammary cancer, NOS**, not lung. They were not added to the lung set.

GSE155972 response labels are confounded with genotype: WT ICB = Non-responders
(n=6); Setdb1_KO ICB = Responders (n=7). There is no within-genotype R vs NR
contrast. Baseline is untreated (n=10 per genotype).

## Primary lung tests (keep outlier)

TISMO values are log2(TPM+1). TPM = 2^value − 1.

| Gene | Arm | nB | nT | mean B→T (log2p1) | mean TPM B→T | Welch p | MWU p | MWU FDR |
|---|---|---:|---:|---|---|---:|---:|---:|
| Tacstd2 | WT (NR) | 10 | 6 | 0.157 → 0.290 | 0.11 → 0.22 | 0.146 | **0.045** | 0.179 |
| Tacstd2 | Setdb1_KO (R) | 10 | 7 | 0.280 → 1.005 | 0.21 → 1.01 | 0.366 | 0.305 | 0.407 |
| Cldn4 | WT (NR) | 10 | 6 | 0.128 → 0.209 | 0.09 → 0.16 | 0.169 | 0.174 | 0.348 |
| Cldn4 | Setdb1_KO (R) | 10 | 7 | 0.267 → 0.601 | 0.20 → 0.52 | 0.482 | 0.887 | 0.887 |

Drop **SRX8918393** (KO ICB; Tacstd2=5.41, Cldn4=3.23): KO Tacstd2 mean goes
0.280 → 0.271 (slightly **down**); KO Cldn4 0.267 → 0.163 (down). WT numbers
are unchanged. TISMO’s LLC-only Tacstd2 export drops this sample; the
All-models export keeps it. This analysis keeps it in the primary table and
drops it only in the sensitivity rows of `paired_icb_stats.tsv`.

Tacstd2 vs Cldn4 Spearman in LLC: ρ=0.33 (p=0.063, n=33); ρ=0.26 (p=0.15)
without the outlier.

## Controls (same LLC ICB design)

| Gene | WT MWU p | KO MWU p | Direction |
|---|---:|---:|---|
| Cd8a | 0.0017 | 0.014 | up after ICB |
| Cd274 | 0.0010 | 0.0020 | up after ICB |
| Gzmb | 0.0075 | 0.0097 | up after ICB |
| Ifng | 0.011 | 0.27 | up in WT |
| Actb | 0.015 (down) | 0.071 | housekeeping, not up |
| Epcam | 0.87 | 0.56 | near floor, flat |

## Honest limits

- This is **not** 64 lung models. It is 1 study, 1 line, 2 genotype arms.
- LLC is a mouse Lewis lung carcinoma, not a human NSCLC histology series.
- Implantation in the ICB study is subcutaneous, not orthotopic lung.
- A4’s 49/64 is real **only** as a pan-cancer mean-direction count. Most of
  those 49 shifts are tiny and not individually significant in TISMO’s own
  p-value column.
- Cldn4 was not part of the original A4 claim; it does not show a 49/64-style
  bias even pan-TISMO.
- No GEO FASTQ reprocessing. Values are TISMO’s uniformly processed export.
- No in-vitro cytokine arm (out of scope for paired ICB).

## Files

| File | Role |
|---|---|
| `summary.json` / `audit.json` | Verdict and method notes |
| `lung_design_inventory.tsv` | All TISMO lung-carcinoma design groups |
| `icb_cohort_universe.tsv` | 64 Tacstd2 + 65 Cldn4 ICB treated-vs-baseline cohorts |
| `direction_tally.tsv` | UP/DOWN/TIE + binomial p |
| `per_sample_expression.tsv` | LLC samples, targets + controls |
| `primary_lung_tests.tsv` | 4 primary tests + BH FDR |
| `paired_icb_stats.tsv` | Primary, outlier-dropped, and control tests |
| `tacstd2_cldn4_correlation.tsv` | Spearman in LLC |
| `figures/fig_llc_paired_icb.png` | Tacstd2/Cldn4 by arm |
| `figures/fig_direction_tally.png` | A4 49/64 vs Cldn4 |
| `figures/fig_llc_expression_floor.png` | Floor vs Actb/immune genes |
| `raw/` | Live TISMO exports used here |

## Reproduce

```bash
python3 scripts/w200/A4_lung/fetch_tismo.py
python3 scripts/w200/A4_lung/analyze.py
```

Source: [TISMO](https://tismo.pku-genomics.org/) (Zeng et al., *NAR* 2022,
PMID 34534350). GSE155972 = Griffin et al., *Nature* 2021 (Setdb1 / ICB).
