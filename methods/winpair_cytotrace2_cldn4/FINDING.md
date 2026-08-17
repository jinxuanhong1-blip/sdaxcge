# Finding — winning-pair GSE131907+GSE205335 malignant cells, CLDN4 vs CytoTRACE2

ADDITIVE. **CLDN4 only.** Winning pair from the CLDN4-first combinatorial search (author malignant CLDN4 %pos GSE131907+GSE205335 vs T/NK). This folder does **not** use GSE148071. GSE207422 is not added. No TACSTD2∩CLDN4 dual-high gate.

Primary potency = **CytoTRACE2** (`cytotrace2-py` 1.1.0.4, Kang et al. *Nat Methods* 2025). Score 0 = differentiated, 1 = totipotent. This is not a residual-n_genes dump. Inferential unit = **patient** (GSE131907 `Sample`, GSE205335 `patient`). Cell-level ρ is descriptive. Barrier/keratin **excludes CLDN4**.

**Question.** Are CLDN4-high author-malignant cells more differentiated / barrier-locked (lower potency) than CLDN4-low cells? **Answer:** at the **patient** unit, higher mean CLDN4 tracks **lower** CytoTRACE2 (n=43, ρ=-0.364, p=0.0165). The association is GSE131907-driven; GSE205335-only is null. Within-patient potency Δ is significant but tiny. CLDN4-high cells are strongly barrier/keratin-high in the same patient.

## Verdict

Patient-level CLDN4 vs CytoTRACE2/potency: n=43, ρ=-0.364, p=0.0165. CLDN4 vs barrier/keratin (CLDN4 excluded): n=43, ρ=0.288, p=0.0615. Potency vs barrier/keratin: n=43, ρ=-0.525, p=3.02e-04. Paired within-patient CLDN4-high vs low potency: n=43, W=193.0, Δmed=-0.006, p=0.0012. Paired barrier/keratin: n=43, W=0.0, Δmed=+0.357, p=1.65e-08. GSE131907-only: n=21, ρ=-0.499, p=0.0214. GSE205335-only: n=22, ρ=-0.149, p=0.5095. Gulati 2020 gene-count CytoTRACE does **not** agree with CytoTRACE2 (patient-level n=43, ρ=0.259, p=0.0930; paired n=43, W=9.0, Δmed=+0.246, p=3.15e-08). tLung contributes **0** author-malignant cells (Kim labels those cells tS1/tS2/tS3, not `Malignant cells`). Not a TACSTD2 redo. No both-high gate. GSE148071 not used.

**What holds.** Pooled patient-mean CLDN4 vs CytoTRACE2 is negative (n=43, ρ=-0.364, p=0.0165; BH q=0.0231). GSE131907-only is the same direction (n=21, ρ=-0.499, p=0.0214). CLDN4-high cells are barrier/keratin-high in the same patient (n=43, W=0.0, Δmed=+0.357, p=1.65e-08). Higher barrier tracks lower CytoTRACE2 (n=43, ρ=-0.525, p=3.02e-04). CLDN4 tracks the Differentiated fraction (n=43, ρ=0.443, p=0.0029).

**What does not hold.** GSE205335-only CLDN4 vs CytoTRACE2 is null (n=22, ρ=-0.149, p=0.5095). ADC-only n=14 is also null. Between-patient CLDN4 vs barrier is NS (n=43, ρ=0.288, p=0.0615). The paired potency shift is tiny (Δmed=-0.006) — do not quote it as a large within-tumor stemness drop. Gulati 2020 gene-count CytoTRACE is **opposite** CytoTRACE2 at the paired test; that is why this folder is not a residual-n_genes dump.

## Honest n

- Analysis cells after QC (capped ≤200/unit): **n_cells = 7999** (GSE131907 3857, GSE205335 4142).
- Units (GSE131907 Sample + GSE205335 patient): **n_units = 43** (GSE131907 21, GSE205335 22).
- Units with ≥20 malignant cells used for Spearman: **n = 43**.
- Units with ≥8 cells in both CLDN4-high and CLDN4-low arms: **n = 43**.
- Author malignant rule: GSE131907 `Cell_subtype==Malignant cells` (catalog malignant 24784; tLung author-malignant = 0). GSE205335 `lineage.sub==Malignant cells` on non-normal tissues (catalog malignant 28512; analyzed 4142).
- Histology (cells): {'LUAD': 3857, 'ADC': 2542, 'SCLC': 800, 'SQ': 600, 'NUT': 200}.
- CytoTRACE2 potency categories (cells): {'Differentiated': 6155, 'Unipotent': 1292, 'Oligopotent': 386, 'Multipotent': 166}.
- CLDN4 tertile cells: {'low': 2667, 'high': 2667, 'mid': 2665}.
- Genes absent from locked sets: {'barrier_keratin': [], 'malignant_like': [], 'AT2': [], 'focal': [], 'comparator': []}.
- GSE148071 not used. GSE207422 not used. Dual-high not used.
- CytoTRACE2 ran: **True**. Fallback/sensitivity Gulati 2020 ran: **True**.

## Locked choices

- Malignant = author label only. No CopyKAT / inferCNV re-call.
- Cap ≤200 cells / unit after QC (honest n reports catalog vs analysis).
- CLDN4 = log1p(CP10k) from raw UMI.
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).
- Inferential n: GSE131907 `Sample` + GSE205335 `patient`.
- Extra figure: within-unit CLDN4-high vs low potency and barrier (min 8 cells/arm).
- Unused: GSE148071; GSE207422; dual-high TACSTD2∩CLDN4; ICI / MPR / RECIST as a potency label.

## Primary (patient-level Spearman, BH inside this list)

| Contrast | n_units | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs CytoTRACE2 | 43 | -0.364 | 0.0165 | 0.0231 |
| CLDN4 vs barrier/keratin (no CLDN4) | 43 | 0.288 | 0.0615 | 0.0718 |
| CytoTRACE2 vs barrier/keratin (no CLDN4) | 43 | -0.525 | 3.02e-04 | 0.0011 |
| CLDN4 vs frac Differentiated | 43 | 0.443 | 0.0029 | 0.0068 |
| CLDN4 vs malignant-like | 43 | 0.556 | 1.08e-04 | 7.54e-04 |
| CLDN4 vs TACSTD2 (comparator) | 43 | 0.381 | 0.0116 | 0.0203 |
| CLDN4 vs Gulati2020 CytoTRACE | 43 | 0.259 | 0.0930 | 0.0930 |

## Sensitivity (not in the BH family)

| Contrast | n_units | ρ | p |
| --- | ---: | ---: | ---: |
| GSE131907-only CLDN4 vs CytoTRACE2 | 21 | -0.499 | 0.0214 |
| GSE131907-only CLDN4 vs barrier/keratin (no CLDN4) | 21 | 0.199 | 0.3879 |
| GSE205335-only CLDN4 vs CytoTRACE2 | 22 | -0.149 | 0.5095 |
| GSE205335-only CLDN4 vs barrier/keratin (no CLDN4) | 22 | 0.248 | 0.2660 |
| GSE205335 ADC-only CLDN4 vs CytoTRACE2 | 14 | -0.310 | 0.2809 |
| GSE205335 ADC+SQ CLDN4 vs CytoTRACE2 | 17 | -0.326 | 0.2016 |
| GSE131907 mBrain CLDN4 vs CytoTRACE2 | 10 | -0.442 | 0.2004 |
| GSE131907 mLN CLDN4 vs CytoTRACE2 | 7 | -0.607 | 0.1482 |
| GSE131907 tL/B CLDN4 vs CytoTRACE2 | 4 | 0.000 | 1.0000 |

## Extra figure — CLDN4-high vs CLDN4-low (unit-paired)

Emitted: **True**. Rule: patient-level Spearman(CLDN4, potency) or Spearman(CLDN4, barrier) p<0.05, or any paired Wilcoxon p<0.05, or n_paired≥4.

| Paired contrast (high − low) | n_units | Δ median | p |
| --- | ---: | ---: | ---: |
| potency high vs low | 43 | -0.006 | 0.0012 |
| barrier/keratin (no CLDN4) high vs low | 43 | 0.357 | 1.65e-08 |
| Gulati2020 high vs low | 43 | 0.246 | 3.15e-08 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- tLung is not in the author-malignant universe. Do not quote this as a tLung AT2 trajectory.
- GSE205335 mixes ADC / SQ / SCLC / NUT. ADC-only n is in the sensitivity table.
- Malignant is the author label, not CNV.
- This is not an ICI / MPR / RECIST test.
- No TACSTD2∩CLDN4 both-high gate.
- GSE148071 is not used.
- Do not write “CLDN4 marks a stem-like malignant state” unless the patient-level potency ρ is negative and significant.

## Outputs

- `results/tables/patient_cldn4_vs_potency.tsv` — **done criterion**
- `results/tables/patient_means.tsv`
- `results/tables/paired_high_vs_low.tsv`
- `results/tables/stats.tsv`
- `results/figures/fig_patient_cldn4_vs_potency.png`
- `results/figures/fig_extra_paired_potency.png`
- `results/figures/fig_extra_cldn4_vs_barrier.png`
- `results/figures/fig_extra_potency_category.png`
- `results/figures/fig_honest_n.png`
- `results/figures/fig_umap_cldn4.png`
- `results/figures/fig_umap_potency.png`
- `results/summary.json`

## Reproduce

```bash
python3 -m venv /tmp/winpair_ct2_venv
/tmp/winpair_ct2_venv/bin/pip install -r methods/winpair_cytotrace2_cldn4/requirements.txt
/tmp/winpair_ct2_venv/bin/python methods/winpair_cytotrace2_cldn4/scripts/download.py \
  --out /tmp/winpair_cytotrace2_cldn4
/tmp/winpair_ct2_venv/bin/python methods/winpair_cytotrace2_cldn4/scripts/extract_malignant.py \
  --data /tmp/winpair_cytotrace2_cldn4 \
  --out /tmp/winpair_cytotrace2_cldn4/malignant.h5ad
/tmp/winpair_ct2_venv/bin/python methods/winpair_cytotrace2_cldn4/scripts/analyze.py \
  --input /tmp/winpair_cytotrace2_cldn4/malignant.h5ad \
  --outdir methods/winpair_cytotrace2_cldn4/results \
  --finding methods/winpair_cytotrace2_cldn4/FINDING.md
```
