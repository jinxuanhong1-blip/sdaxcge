# methods/scrna_scvi_cldn4

Additive public-only scRNA. No claim-audit. **Not a TACSTD2 redo.**

## Question

After a real scVI (else Harmony) latent on **GSE207422** (batch = 10x sample), what is the sample-level Spearman of **malignant CLDN4** vs T/NK? Honest n / ρ / p.

## Why this is additive

| Existing | Not this |
|---|---|
| scVI combo GSE207422 + GSE241934 IIT scored TACSTD2 (n=26) | GSE207422 only; gene = CLDN4 |
| Marker-only GSE207422 malignant CLDN4 (no latent) | Same UMI, malignant also from latent Leiden |
| A3 TACSTD2 NMPR/MPR / T/NK | Taken as given |

TACSTD2 is a companion column only. Dual-high is not run.

## Unit

12 post-treatment patients (MPR n=4 including pCR P06; NMPR n=8). The three pre-treatment biopsies stay in the latent and in the 15-sample sensitivity table. Cell-level p-values are not reported. Eligible: ≥10 malignant-like and ≥20 T/NK cells.

## Reproduce

```bash
python3 methods/scrna_scvi_cldn4/download.py
python3 methods/scrna_scvi_cldn4/prepare.py
python3 methods/scrna_scvi_cldn4/integrate_score.py
```

scVI is preferred. If training fails, the script falls back to Harmony and records the error.

## Result (public run)

scVI trained (no Harmony fallback). scvi-tools 1.5.0.post1, 40 epochs, `n_latent=10`, 2,000 HVGs, batch = 10x sample.

Primary (post, Leiden-cluster malignant CLDN4 vs T/NK): **n=10/12, ρ=−0.03, p=0.93**. Dropped P06 (5 cells) and P13 (3 cells). Not significant.

Write-up: `FINDING.md`.
