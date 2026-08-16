# A10 KLF4 vs TACSTD2 / CLDN4 — public lung

**Claim:** KLF4 is associated in the same direction with TACSTD2 and with CLDN4 in public lung.

**Honest verdict (TCGA-LUAD primary): `PARTIAL`**

TCGA-LUAD primary tumors: KLF4 is significant for one target only (TACSTD2=POSITIVE, CLDN4=NULL). Same-direction tracking of both genes is not supported. TCGA-LUSC (pre-specified second histology): TACSTD2=POSITIVE, CLDN4=POSITIVE.

Filters were not tuned to force concordance.

## Primary pairs (BH-FDR within this list)

| Cohort | Layer | Pair | n | Spearman ρ | p | FDR | 95% CI | Direction |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- |
| TCGA-LUAD | RNA_log2TPM | KLF4 vs TACSTD2 | 516 | 0.207 | 2.22e-06 | 6.2e-06 | [0.119, 0.290] | POSITIVE |
| TCGA-LUAD | RNA_log2TPM | KLF4 vs CLDN4 | 516 | -0.004 | 9.33e-01 | 0.933 | [-0.095, 0.082] | NULL |
| TCGA-LUSC | RNA_log2TPM | KLF4 vs TACSTD2 | 501 | 0.369 | 1.20e-17 | 1.68e-16 | [0.285, 0.447] | POSITIVE |
| TCGA-LUSC | RNA_log2TPM | KLF4 vs CLDN4 | 501 | 0.198 | 7.99e-06 | 1.86e-05 | [0.112, 0.288] | POSITIVE |
| CPTAC-LUAD | RNA_RSEM_UQ_log2 | KLF4 vs TACSTD2 | 110 | 0.049 | 6.10e-01 | 0.657 | [-0.160, 0.242] | NULL |
| CPTAC-LUAD | RNA_RSEM_UQ_log2 | KLF4 vs CLDN4 | 110 | -0.312 | 9.18e-04 | 0.00161 | [-0.480, -0.126] | INVERSE |
| CPTAC-LUAD | protein_TMT_log2 | KLF4 vs TACSTD2 | 85 | -0.084 | 4.43e-01 | 0.564 | [-0.306, 0.141] | NULL |
| CPTAC-LUAD | protein_TMT_log2 | KLF4 vs CLDN4 | 60 | 0.088 | 5.05e-01 | 0.589 | [-0.176, 0.354] | NULL |
| CPTAC-LSCC | RNA_RSEM_UQ_log2 | KLF4 vs TACSTD2 | 108 | 0.403 | 1.50e-05 | 3e-05 | [0.231, 0.553] | POSITIVE |
| CPTAC-LSCC | RNA_RSEM_UQ_log2 | KLF4 vs CLDN4 | 108 | 0.171 | 7.76e-02 | 0.109 | [-0.016, 0.353] | NULL |
| CPTAC-LSCC | protein_TMT_log2 | KLF4 vs TACSTD2 | 103 | 0.494 | 1.16e-07 | 4.06e-07 | [0.314, 0.635] | POSITIVE |
| CPTAC-LSCC | protein_TMT_log2 | KLF4 vs CLDN4 | 73 | 0.208 | 7.80e-02 | 0.109 | [-0.028, 0.422] | NULL |
| DepMap24Q4-lung | RNA_log2TPM | KLF4 vs TACSTD2 | 214 | 0.401 | 1.14e-09 | 7.99e-09 | [0.279, 0.512] | POSITIVE |
| DepMap24Q4-lung | RNA_log2TPM | KLF4 vs CLDN4 | 214 | 0.392 | 2.98e-09 | 1.39e-08 | [0.266, 0.501] | POSITIVE |

## Exploratory (not used for the verdict)

| Cohort | Pair | n | Spearman ρ | p | Direction |
| --- | --- | ---: | ---: | ---: | --- |
| TCGA-LUAD | TACSTD2 vs CLDN4 | 516 | 0.531 | 7.19e-39 | POSITIVE |
| TCGA-LUAD | KLF4 vs KRT5 | 516 | 0.251 | 7.68e-09 | POSITIVE |
| TCGA-LUAD | KLF4 vs CLDN7 | 516 | -0.088 | 4.54e-02 | INVERSE |
| TCGA-LUAD | KLF4 vs PECAM1 | 516 | 0.332 | 9.43e-15 | POSITIVE |
| TCGA-LUSC | TACSTD2 vs CLDN4 | 501 | 0.391 | 8.59e-20 | POSITIVE |
| TCGA-LUSC | KLF4 vs KRT5 | 501 | 0.323 | 1.23e-13 | POSITIVE |
| TCGA-LUSC | KLF4 vs CLDN7 | 501 | 0.165 | 2.03e-04 | POSITIVE |
| TCGA-LUSC | KLF4 vs PECAM1 | 501 | -0.045 | 3.15e-01 | NULL |
| CPTAC-LUAD | TACSTD2 vs CLDN4 | 110 | 0.392 | 2.30e-05 | POSITIVE |
| CPTAC-LUAD | TACSTD2 vs CLDN4 | 79 | 0.275 | 1.43e-02 | POSITIVE |
| CPTAC-LSCC | TACSTD2 vs CLDN4 | 108 | 0.446 | 1.30e-06 | POSITIVE |
| CPTAC-LSCC | TACSTD2 vs CLDN4 | 78 | 0.081 | 4.83e-01 | NULL |
| DepMap24Q4-lung | TACSTD2 vs CLDN4 | 214 | 0.607 | 5.64e-23 | POSITIVE |
| DepMap24Q4-LUAD | KLF4 vs TACSTD2 | 80 | 0.384 | 4.35e-04 | POSITIVE |
| DepMap24Q4-LUAD | KLF4 vs CLDN4 | 80 | 0.403 | 2.12e-04 | POSITIVE |
| DepMap24Q4-LUAD | TACSTD2 vs CLDN4 | 80 | 0.756 | 5.11e-16 | POSITIVE |

## Coverage (honest missingness)

- TCGA-LUAD primary tumors: n=516; ABSOLUTE purity available for 503.
- TCGA-LUSC primary tumors: n=501; ABSOLUTE purity available for 493.
- CPTAC LUAD RNA: n=110; non-NA KLF4/TACSTD2/CLDN4 = 110/110/110.
- CPTAC LUAD protein: n=110; non-NA KLF4/TACSTD2/CLDN4 = 85/110/79.
- CPTAC LSCC RNA: n=108; non-NA KLF4/TACSTD2/CLDN4 = 108/108/108.
- CPTAC LSCC protein: n=108; non-NA KLF4/TACSTD2/CLDN4 = 103/108/78.
- DepMap 24Q4 lung cell lines: n=214; LUAD subset n=80.

## Secondary: purity-adjusted TCGA (does not change the verdict)

- TCGA-LUAD KLF4 vs TACSTD2: partial ρ = 0.238, p = 6.82e-08, n = 503
- TCGA-LUAD KLF4 vs CLDN4: partial ρ = 0.029, p = 5.17e-01, n = 503
- TCGA-LUAD KLF4 vs KRT5: partial ρ = 0.236, p = 8.40e-08, n = 503
- TCGA-LUAD KLF4 vs CLDN7: partial ρ = -0.045, p = 3.16e-01, n = 503
- TCGA-LUAD KLF4 vs PECAM1: partial ρ = 0.303, p = 4.30e-12, n = 503
- TCGA-LUSC KLF4 vs TACSTD2: partial ρ = 0.374, p = 8.30e-18, n = 493
- TCGA-LUSC KLF4 vs CLDN4: partial ρ = 0.172, p = 1.29e-04, n = 493
- TCGA-LUSC KLF4 vs KRT5: partial ρ = 0.316, p = 7.18e-13, n = 493
- TCGA-LUSC KLF4 vs CLDN7: partial ρ = 0.144, p = 1.39e-03, n = 493
- TCGA-LUSC KLF4 vs PECAM1: partial ρ = 0.041, p = 3.58e-01, n = 493

## What this is not

- Direct KLF4 binding or transcriptional control at TACSTD2/CLDN4 (no ChIP/perturbation in this slice).
- ICI response or TROP2-ADC outcome.
- KLF4 IHC or protein from RNA (except the CPTAC protein rows, which are reported separately).
- A reason to drop CPTAC protein or DepMap because the sign disagrees with TCGA RNA.
- Evidence that KLF4 transcriptionally activates or represses TACSTD2 or CLDN4.

## Caveats

1. Bulk tumor RNA mixes epithelium, stroma, endothelium, and immune cells. In TCGA-LUAD, KLF4 tracks PECAM1 (ρ = 0.33) as well as TACSTD2 (ρ = 0.21). That does not prove the TACSTD2 association is endothelial, but it forbids reading KLF4 as a purely epithelial/TJ transcription factor from bulk RNA.
2. TCGA-LUAD KLF4–CLDN4 is a true null (ρ = −0.004, CI includes 0), not a weak inverse that we rounded away. Purity adjustment does not create an inverse.
3. CPTAC-LUAD RNA is not a silent non-replication: KLF4–CLDN4 is significantly inverse (ρ = −0.31, n=110) while KLF4–TACSTD2 is null. That is the opposite of a same-sign pair. CPTAC-LUAD protein is null for both, with KLF4 missing in 25/110 and CLDN4 missing in 31/110 tumors.
4. TCGA-LUSC and DepMap lung lines are concordant-positive. They are reported; they do not override the pre-specified LUAD verdict. DepMap lung n=214 includes NSCLC, neuroendocrine, and a few non-cancerous lines; the LUAD-only subset (n=80) is also positive for both pairs.
5. The TACSTD2–CLDN4 pair itself is robustly positive in TCGA LUAD/LUSC and DepMap. The assay is not broken; the KLF4 same-sign claim is.
6. TCGA/CPTAC/DepMap are not ICI cohorts.

## Rerun

```bash
python3 scripts/w200/A10_KLF4/download.py
python3 scripts/w200/A10_KLF4/analyze.py
```

See `summary.json` and `correlations.csv`.
