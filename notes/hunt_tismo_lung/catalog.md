# TISMO lung-model catalog (complete, not ICB-only)

Source: TISMO `cellLineMeta` + `vivoMeta` + official `TISMO_vivosample_annotations.csv` (1,518 in-vivo samples) and `TISMO_vitrosample_annotations.csv` (605 in-vitro samples). Retrieved 2026-08-16 from `https://tismo.pku-genomics.org/`.

## Cell lines labeled lung in TISMO

| Cell line | Cancer_type | In vivo? | In vitro? | ICB pairing? |
|---|---|---|---|---|
| LLC | Lung carcinoma | yes (68 samples, 9 studies) | yes (55 samples) | **yes — GSE155972 only** |
| CMT-167 | Lung carcinoma | yes (3 samples, GSE100412) | no | **no** |
| MLE12 | Lung adenocarcinoma | **no** | yes (4 untreated, GSE103548) | **no** |

No other TISMO cell line is annotated `Lung carcinoma` or `Lung adenocarcinoma`. KP, 344SQ, LKR13, 393P, and other common murine NSCLC lines are **not in TISMO**.

## In-vivo lung studies (n=71 samples)

| Study | Line | n | Implantation | ICB | Usable for claim A4 pairing? |
|---|---|---|---|---|---|
| **GSE155972** | LLC WT | 10 baseline + 6 anti-PD1+anti-CTLA4 (NR) | subcutaneous flank | yes | **yes** (1 of 2 official groups) |
| **GSE155972** | LLC Setdb1-KO | 10 baseline + 7 anti-PD1+anti-CTLA4 (R) | subcutaneous flank | yes | **yes** (2 of 2). Response is confounded with genotype. |
| GSE100412 | CMT-167, LLC | 3+5 | orthotopic left lung | no | no |
| GSE115109 | LLC | 8 | subcutaneous | no | no |
| GSE131271 | LLC ± shSocs1 | 6 | orthotopic left lung | no | no |
| GSE148101 | LLC | 2 | intravenous (metastasis) | no | no |
| GSE71491 | LLC ± cyclophosphamide | 4 | subcutaneous | no | no (chemo, not ICB) |
| GSE80678 | LLC | 2 | subcutaneous | no | no |
| GSE84535 | LLC ± host p110γ-KO | 6 | subcutaneous | no | no |
| GSE98672 | LLC | 2 | subcutaneous back | no | no |

**Honest ICB universe:** 1 cell line, 1 study, 2 official TISMO groups, 33 samples. The ICB tumors are **subcutaneous**, not orthotopic lung.

## In-vitro lung

- **LLC** (RTM28723893 / Manguso in-house): parental + clones 7/8/9, IFNβ / IFNγ / TNFα / no_treatment. Tacstd2 and Cldn4 are at the floor (almost all exact 0). Cd274 induces with IFN (positive control).
- **MLE12** (GSE103548): untreated only. Tacstd2 is high (mean log2(TPM+1) ≈ 5.88); Cldn4 is ~0. No cytokine and no ICB arm.

## What 49/64 is

Official TISMO Gene-module ICB CSV groups (`cell_line` field, e.g. `LLC_GSE155972_antiCTLA4&antiPD1(n=16)`). Pairing is `Baseline==1` vs `Baseline==0` (R and NR pooled). Tacstd2 is missing the in-house `MOC22_RU31562203_antiPD1` group, which is why Tacstd2 has 64 groups and Cldn4 has 65.
