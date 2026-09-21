# FINDING — concordant-4 decoupler / DoRothEA (malignant CLDN4-high vs low)

ADDITIVE. **CLDN4 only.** No TACSTD2∩CLDN4 dual-high.
Cohorts are the locked concordant-4 set: **GSE123902 + GSE131907 + GSE205335 + GSE189357**.
Not GSE148071, GSE127465, GSE207422, or GSE154826.

The unit is the patient / donor / sample malignant pseudobulk, with the **PR #503 quartile labels** (within-cohort malignant CLDN4 %pos). This is not the within-patient cell-tail contrast of PR #454.

p-values are descriptive.

## What the KD direction predicts, and what this table is

External expectation, not estimated here: CLDN4 loss is followed by **IFN / MHC up** (private KD; public lung KO GSE50927, PR #43). If a CLDN4-high vs low observation had that same sign, CLDN4-high malignant cells would show **lower** STAT1 / IRF / NLRC5 / MHC activity.

Junction TFs do not have that prediction. The same public KO did not collapse junctions.

This table is observational. A matching sign is not evidence that CLDN4 causes the IFN change.

## Honest n

| Item | n | Note |
|---|---:|---|
| Concordant-4 units with a CLDN4 %pos | 65 | 13 + 21 + 22 + 9. T/NK n from PR #503. Not the TF n. |
| Units in the malignant UMI matrices | **64** | P4001 (GSE205335) has no malignant pseudobulk |
| **Q4 vs Q1 stacked** | **34** | **18 low / 16 high** |
| GSE189357 inside that 34 | 3 low / 2 high | kept in the stack; single-cohort test skipped |
| DoRothEA A+B+C TFs with ≥5 targets in the matrix | 297 | IRF7 has 2 targets here and is not scored |

Do not quote 65 as the activity n.

## Verdict

On this patient-level pseudobulk, DoRothEA **STAT1, STAT2, and IRF1 activity is lower in CLDN4-high**. That sign **matches** the KD expectation and matches the IFN gene score on the same matrix (Hallmark α∪γ −0.584, p=0.00243; the PR #503 number).

It does **not** match the earlier cell-level DoRothEA result. PR #454 (winning pair, within-patient malignant cell Q4 vs Q1, raw wmean) found STAT1 **higher** in CLDN4-high cells (0.175 vs 0.129, Δ=+0.046, 50/52, p=1.18×10⁻⁹). That sign is **opposite of KD**. This analysis does not re-run that cell split and does not erase it. The two contrasts answer different questions. Do not describe the cell-level IFN-up-in-high as the knockdown direction.

NLRC5 is not a regulator in DoRothEA A–D or in CollecTRI, so there is no NLRC5 activity score. NLRC5 mRNA itself is lower in CLDN4-high (−0.485, p=0.21): the sign matches KD and the estimate is too wide to claim.

DoRothEA junction TFs (GRHL2, ELF3, KLF4, TFAP2A) are not different at n=34. There is no stacked junction-TF program opposing the IFN result.

## Primary table — ULM, cohort-adjusted OLS, Q4 vs Q1 (n=18 / 16)

Activity is the decoupler 2.2.0 univariate linear model t-value (DoRothEA A+B+C, confidence-scaled weights). Positive coefficient = higher activity in CLDN4-high.

| program | TF | targets | coef | p | BH FDR (297 TFs) | vs KD |
|---|---|---:|---:|---:|---:|---|
| IFN | STAT1 | 1062 | −2.132 | 0.00584 | 0.058 | matches (lower in high) |
| IFN | STAT2 | 102 | −2.169 | 0.00443 | 0.058 | matches |
| IFN | IRF1 | 74 | −1.407 | 0.00696 | 0.061 | matches |
| IFN | IRF2 | 73 | −0.547 | 0.086 | 0.30 | same sign, not p<0.05 |
| IFN | IRF3 | 69 | −0.189 | 0.40 | 0.68 | same sign, not p<0.05 |
| IFN | IRF7 | 2 | — | — | — | not scored (tmin=5) |
| IFN | IRF8 | 8 | −0.076 | 0.53 | 0.77 | flat |
| IFN | IRF9 | 45 | −0.539 | 0.11 | 0.34 | same sign, not p<0.05 |
| MHC | RFX5 | 87 | −1.388 | 0.073 | 0.26 | same sign, not p<0.05 |
| junction | GRHL2 | 74 | +0.481 | 0.44 | 0.70 | no KD sign; not different |
| junction | ELF3 | 52 | +0.467 | 0.33 | 0.61 | no KD sign; not different |
| junction | KLF4 | 48 | −0.095 | 0.75 | 0.91 | no KD sign; flat |
| junction | TFAP2A | 123 | +0.009 | 0.98 | 0.99 | no KD sign; flat |

STAT1, STAT2, and IRF1 are pre-specified. Their nominal p-values are below 0.05. Across all 297 scored TFs, Benjamini–Hochberg leaves them at FDR ≈ 0.058. Eleven other TFs reach FDR<0.05; they were not pre-specified and are not a result of this note (see `tf_tests_all.tsv`).

## Same patients, four decoupler statistics

Sign of the Q4 coefficient. “matches KD” means negative for an IFN/MHC regulator.

| TF | ULM | MLM | VIPER | wmean | four-method sign |
|---|---|---|---|---|---|
| STAT1 | −2.132 (0.0058) | −0.849 (0.011) | −0.410 (0.051) | −0.120 (0.098) | all negative |
| STAT2 | −2.169 (0.0044) | −1.721 (0.0039) | −0.866 (0.0064) | −0.508 (0.012) | all negative, all p<0.05 |
| IRF1 | −1.407 (0.0070) | −0.549 (0.063) | −1.063 (0.0048) | −0.399 (0.017) | all negative |
| RFX5 | −1.388 (0.073) | −1.290 (0.073) | −0.720 (0.059) | −0.343 (0.13) | all negative, none p<0.05 |
| GRHL2 | +0.481 (0.44) | +0.212 (0.62) | +0.508 (0.14) | +0.188 (0.31) | all positive, none p<0.05 |
| ELF3 | +0.467 (0.33) | +0.202 (0.48) | +0.177 (0.36) | +0.211 (0.21) | all positive, none p<0.05 |

ULM is the primary statistic (decoupler benchmark). On this matrix the weighted mean does **not** flip STAT1 positive. The IFN-up-in-high pattern is not a general property of wmean. It is what PR #454 saw in within-patient cell tails.

IRF9 is the exception among the weaker IRFs: ULM −0.54 (p=0.11) but MLM +0.25 (p=0.095). Do not treat every IRF as concordant. The agreeing set is STAT1, STAT2, and IRF1.

Replacing confidence-scaled weights with unsigned ±1 does not flip STAT1 (ULM −1.83, p=0.021).

## Why a DoRothEA IFN score can disagree with the IFN gene program

STAT1’s A+B+C regulon in this matrix has 1,062 targets. **83** overlap Hallmark IFN-α∪γ. Those 83 have mean log2FC −0.66 (90% negative in CLDN4-high). The other **979** have mean log2FC −0.03 (49% negative).

ULM is a within-sample regression of expression on the signed weights, so the IFN-overlapping tail can pull STAT1 activity down even though most targets are not IFN genes. A cell-level weighted mean can instead be carried by that non-IFN majority and by expression depth. PR #454’s OXPHOS control rose with CLDN4-high cells (33/52, Δ=+0.037, p=0.015). On this pseudobulk the OXPHOS control does not (coef +0.037, p=0.64).

So: patient-level ULM here agrees with the KD sign; cell-level wmean in PR #454 does not. Both statements stay in the record.

## Gene readouts on the same matrix (not TF activity)

These use the PR #503 normalization, log2(TMM-CPM+1), and the same cohort-adjusted Q4 model. IFN and MHC-I scores match PR #503 exactly.

| readout | genes | coef | p | vs KD |
|---|---:|---:|---:|---|
| Hallmark IFN α∪γ | 221 | −0.584 | 0.00243 | matches |
| MHC-I / APM | 21 | −0.779 | 0.00988 | matches |
| TJ, CLDN4 held out | 191 | +0.059 | 0.39 | no KD prediction; flat |
| OXPHOS control | 184 | +0.037 | 0.64 | depth control, flat |
| STAT1 mRNA | 1 | −0.991 | 0.00109 | matches |
| IRF1 mRNA | 1 | −0.577 | 0.075 | same sign |
| NLRC5 mRNA | 1 | −0.485 | 0.21 | same sign; not a regulon |
| CIITA mRNA | 1 | −0.312 | 0.42 | same sign; CIITA is not in DoRothEA |
| RFX5 mRNA | 1 | +0.129 | 0.64 | transcript flat; activity was negative and n.s. |
| CLDN4 mRNA | 1 | +1.626 | 0.014 | splitter check |

NLRC5 and CIITA have no DoRothEA activity. CollecTRI (a different network, labeled as sensitivity) does contain CIITA, RFXANK, and RFXAP, and still does not contain NLRC5.

CollecTRI ULM, same patients: CIITA −1.21 (p=0.029), RFX5 −1.31 (p=0.062), RFXANK −1.43 (p=0.056), RFXAP −1.39 (p=0.074), STAT1 −2.49 (p=0.00039). MHC-regulator signs that can be scored agree with the gene program. CollecTRI GRHL2 has only 10 targets and is up (+1.46, p=9.6×10⁻⁵). That is not the DoRothEA GRHL2 result (p=0.44) and is not used as the junction claim.

## Where the stacked STAT1 result comes from

Within-cohort ULM, Q4 vs Q1, no cohort covariate:

| cohort | n low/high | STAT1 coef | p |
|---|---|---:|---:|
| GSE123902 | 4/3 | −1.96 | 0.33 |
| GSE131907 | 6/5 | −0.11 | 0.92 |
| GSE205335 | 5/6 | −4.58 | 0.0054 |
| GSE189357 | 3/2 | — | not tested |

Leave-one-cohort-out of the stacked model: dropping GSE205335 leaves STAT1 at −0.95 (n=23, p=0.24). The coefficient stays negative in every leave-one-out, and the p<0.05 call does not. This is the same cohort pattern as the PR #503 IFN gene score (GSE205335 logFC −1.24; GSE131907 −0.011).

Continuous CLDN4 %pos (n=64, cohort-adjusted): STAT1 ULM −0.54 per SD (p=0.073); STAT2 −0.69 (p=0.025). The quartile contrast is the primary one.

## Specification sweep

The pre-specified rows above were not replaced. A second script, `scripts/sweep.py`, re-scored the same 64 malignant pseudobulks under the settings listed below and kept every row. Single-cohort specs are in the table and were not eligible to be the reported panel. No cohort was removed to force a sign. All-epithelial pseudobulks are not in the locked count matrices, so that contrast was not invented.

Varied: DoRothEA A, A+B, A+B+C, A–D; CollecTRI; PROGENy top 500; unsigned Hallmark IFN and tight-junction gene sets. Statistics: ULM, MLM, wmean, wsum, and wmean normalized to 100 permutations (NES) for A+B+C and CollecTRI. Cutoffs, within cohort: locked Q4/Q1, equal-count quartile, tertile, median, top/bottom 30%, top/bottom 20%. Scopes: stacked OLS with cohort indicators, and inverse-variance meta of cohorts that had at least 3 units per arm.

A spec is counted as jointly nominal when the IFN-arm feature is lower in CLDN4-high at p<0.05 and at least one barrier feature is higher in CLDN4-high at p<0.05. The IFN arm is STAT1 when the network has it, otherwise the IFN gene set or PROGENy JAK-STAT. The barrier arm is the best of GRHL2, ELF3, KLF4, TFAP2A, the 191-gene TJ set, and Hallmark apical junction.

| network | specs (stacked or meta) | both signs | both nominal p<0.05 | both within-spec FDR<0.05 | of which stacked |
|---|---:|---:|---:|---:|---:|
| CollecTRI | 60 | 60 | 47 | 42 | 24 |
| DoRothEA A | 48 | 43 | 6 | 3 | 3 |
| DoRothEA A+B | 48 | 36 | 6 | 1 | 2 |
| DoRothEA A+B+C | 60 | 60 | 6 | 0 | 1 |
| DoRothEA A–D | 36 | 8 | 0 | 0 | 0 |
| IFN / TJ gene sets | 24 | 23 | 3 | 1 | 0 |

276 stacked or meta specs. 68 are jointly nominal. Those 68 are overlapping re-cuts of the same patients, not 68 independent cohorts. The minimum p after that search is not a confirmatory p-value. Within-spec FDR is only across the handful of features in that spec.

Of the 68, the barrier feature is GRHL2 in 49, TFAP2A in 14, the 191-gene TJ set in 3, and ELF3 in 2. The 191-gene TJ hits are inverse-variance metas, not the stacked model. Apical junction never clears the joint bar.

### Strongest joint panel

Smallest max(p) among stacked or meta specs with both nominal p<0.05: CollecTRI, wmean NES (100 permutations), locked quartiles, inverse-variance meta of the three cohorts with ≥3 per arm (n=15/14; GSE189357 tails are 3/2 and stay out).

| feature | targets | coef | p | vs thesis |
|---|---:|---:|---:|---|
| STAT1 | 284 | −2.520 | 3.9×10⁻⁵ | matches |
| GRHL2 | 10 | +1.535 | 1.9×10⁻⁶ | matches |
| IRF1 | 147 | −1.664 | 0.0078 | matches |
| TFAP2A | 285 | −0.933 | 0.014 | wrong sign |
| ELF3 | 31 | −0.115 | 0.76 | wrong sign, not significant |

The same network and the same locked quartiles, scored by ULM and tested with the original stacked model (n=18/16, GSE189357 included), is the panel on the original sample size:

| feature | targets | coef | p |
|---|---:|---:|---:|
| STAT1 | 284 | −2.486 | 3.9×10⁻⁴ |
| IRF1 | 147 | −1.896 | 0.0047 |
| GRHL2 | 10 | +1.461 | 9.6×10⁻⁵ |
| TFAP2A | 285 | −1.027 | 0.0082 |

GRHL2 is positive in each testable cohort (GSE123902 +2.04, p=0.0080, n=4/3; GSE131907 +1.89, p=0.029, n=6/5; GSE205335 +0.90, p=0.15, n=5/6). STAT1 is not: GSE205335 −4.72 (p=0.0018), GSE123902 −3.26 (p=0.060), GSE131907 −0.31 (p=0.73). The two arms are not carried by the same cohort. CollecTRI GRHL2 has 10 targets. TFAP2A, the large CollecTRI junction TF, is significantly the wrong way. This is not a general barrier-program result.

### What the sweep did not rescue

- Pre-specified DoRothEA A+B+C ULM, locked quartiles, stacked: GRHL2 stays +0.48, p=0.44. In that network GSE205335 GRHL2 is negative (−1.45, p=0.19).
- DoRothEA A+B+C has one stacked joint nominal hit: MLM, tertiles, STAT1 −0.63 (p=0.037) and TFAP2A +0.33 (p=0.041). Within-spec FDR does not pass. The other five DoRothEA A+B+C joint hits are metas, and none pass within-spec FDR.
- DoRothEA A–D produces no joint nominal spec. Adding the D-level edges does not create the barrier result.
- The 191-gene TJ set is significant only in the inverse-variance meta (tail 30% wmean: IFN −0.57, p=0.0023; TJ +0.13, p=0.015). The stacked test of that same cutoff is TJ +0.076, p=0.25. GSE205335’s TJ coefficient is negative.
- PROGENy has no tight-junction pathway. On the locked quartile stack, JAK-STAT ULM is −3.85 (p=0.0032, n=18/16), which matches the IFN arm only.
- Four single-cohort specs are jointly nominal. All four are CollecTRI GRHL2, with tail sizes 4/4, 6/6, or 7/7. They were not eligible as the panel.

Full rows: `results/sweep_panel.tsv`, `results/sweep_long.tsv`, `results/sweep_counts.tsv`.

## What this is not

- Not a re-analysis of within-patient cell Q4 vs Q1, and not a retraction of PR #454.
- Not a causal KD result. Sign agreement is a label against the external IFN-up-after-loss expectation.
- Not an NLRC5 regulon. That TF is absent from both networks used here.
- Not a significant DoRothEA junction-TF program.
- Not a genome-wide TF discovery (pre-specified STAT1/STAT2/IRF1 sit at FDR ≈ 0.058).
- Not cell-level n, and not N=65.
- Not a confirmatory test of the specification search. The joint CollecTRI p-values were selected after seeing 276 stacked or meta specs.

## Figures

- `figures/fig_ulm_forest.png` — primary ULM coefficients
- `figures/fig_method_coefs.png` — ULM / MLM / VIPER / wmean
- `figures/fig_ulm_boxes.png` — patient activities, Q1 vs Q4
- `figures/fig_ulm_by_cohort.png` — which cohort carries STAT1
- `figures/fig_stat1_targets.png` — IFN-overlap vs the rest of the STAT1 regulon
- `figures/fig_sweep_signs.png` — IFN-arm vs barrier-arm coefficients across the sweep
- `figures/fig_sweep_best.png` — smallest joint nominal spec (CollecTRI wmean NES, locked quartiles, meta)

## Reproduce

```bash
python3 methods/concordant4_dorothea_decoupler_cldn4/scripts/analyze.py
```

Inputs are the PR #503 malignant pseudobulk UMI sums and quartile labels. DoRothEA A+B+C is the committed OmniPath table (`resources/dorothea_hs_ABC.tsv`).
