# IMpower133 — what is and is not public ("the leftovers")

**Short version: IMpower133 per-patient RNA-seq is NOT public.** I did not, and
cannot, run a real IMpower133 expression analysis from open data, and I have
**not fabricated one**.

## What IMpower133 is
Phase 3 RCT (Horn et al. 2018, NEJM) in treatment-naïve extensive-stage SCLC:
carboplatin/etoposide ± atezolizumab (anti-PD-L1). Established EP+atezo as a
standard of care. Baseline tumor RNA-seq (n≈276 treatment-naïve) was later used
by Gay et al. 2021 (Cancer Cell) to validate SCLC-A/N/P/I and to show SCLC-I
derives the greatest benefit from adding immunotherapy.

## What is controlled access
- **RNA-seq (subtype-classification gene subset)**: EGA `EGAD00001006928` —
  access by application, 1-year publication embargo.
- **RNA-seq (full transcriptome)**: EGA `EGAD00001006927` — not granted until a
  follow-up publication.
- **Subtype assignments (per sample)**: EGA `EGAD00001006926` — by application.
- Study umbrella: EGA `EGAS00001004888` (DAC `EGAC00001001932`). Roche/Genentech
  controlled; individual patient records cannot be linked across sources.

Because there is no PD-L1×gene-expression×outcome table released at the patient
level, **TACSTD2/CLDN4 vs immune-subtype in IMpower133 cannot be computed from
public data.**

## The genuinely public "leftovers" (from the Gay 2021 paper text/figures)
These are qualitative, cohort-level statements — usable as context, not as
re-analyzable data:
- IMpower133 subtype distribution: SCLC-A 51%, SCLC-N 23%, SCLC-I 18%, SCLC-P 7%.
- The Ayers 18-gene IFN-γ T-cell GEP is high specifically in SCLC-I (their Fig 3G).
- OS benefit of adding atezolizumab is numerically greatest in SCLC-I; SCLC-I vs
  all-other in the atezo arm HR≈0.566 (95% CI 0.321–0.998), not seen in placebo —
  i.e. SCLC-I looks predictive, not merely prognostic (their Fig 3H–J, Sup Fig 4).

## How to run the analysis IF you obtain EGA access
`src/hunt_sclc_ici/impower133_template.py` is a ready-to-run template that
consumes an expression matrix + subtype table in the same shape as the George
inputs and produces the identical TACSTD2/CLDN4-vs-immune-subtype outputs. It
**refuses to run on placeholder data** and prints acquisition instructions,
precisely so that no fake IMpower133 result can be produced by accident.
