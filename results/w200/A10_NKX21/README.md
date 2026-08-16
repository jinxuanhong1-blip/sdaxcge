# A10 NKX2-1 inverse with TACSTD2 / CLDN4 — public LUAD

**Claim:** NKX2-1 is inversely correlated with TACSTD2 and with CLDN4 in public LUAD.

**Honest verdict (TCGA-LUAD primary): `OPPOSITE`**

TCGA-LUAD primary tumors: NKX2-1 is significantly positively correlated with CLDN4 (TACSTD2 is NULL, not inverse). The inverse claim is not supported.

Filters were not tuned to force an inverse.

## Primary pairs (BH-FDR within this list)

| Cohort | Layer | Pair | n | Spearman ρ | p | FDR | 95% CI | Direction |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- |
| TCGA-LUAD | RNA_log2TPM | NKX2-1 vs TACSTD2 | 516 | -0.020 | 6.52e-01 | 0.768 | [-0.112, 0.071] | NULL |
| TCGA-LUAD | RNA_log2TPM | NKX2-1 vs CLDN4 | 516 | 0.240 | 3.55e-08 | 2.84e-07 | [0.156, 0.321] | POSITIVE |
| CPTAC-LUAD | RNA_RSEM_UQ_log2 | NKX2-1 vs TACSTD2 | 110 | -0.009 | 9.28e-01 | 0.928 | [-0.206, 0.192] | NULL |
| CPTAC-LUAD | RNA_RSEM_UQ_log2 | NKX2-1 vs CLDN4 | 110 | 0.273 | 3.89e-03 | 0.00778 | [0.078, 0.448] | POSITIVE |
| CPTAC-LUAD | protein_TMT_log2 | NKX2-1 vs TACSTD2 | 110 | 0.041 | 6.72e-01 | 0.768 | [-0.173, 0.247] | NULL |
| CPTAC-LUAD | protein_TMT_log2 | NKX2-1 vs CLDN4 | 79 | -0.216 | 5.62e-02 | 0.0898 | [-0.432, 0.009] | NULL |
| DepMap24Q4-LUAD | RNA_log2TPM | NKX2-1 vs TACSTD2 | 80 | 0.445 | 3.54e-05 | 9.44e-05 | [0.230, 0.634] | POSITIVE |
| DepMap24Q4-LUAD | RNA_log2TPM | NKX2-1 vs CLDN4 | 80 | 0.532 | 3.78e-07 | 1.51e-06 | [0.342, 0.687] | POSITIVE |

## Exploratory (not used for the verdict)

| Cohort | Pair | n | Spearman ρ | p | Direction |
| --- | --- | ---: | ---: | ---: | --- |
| TCGA-LUAD | TACSTD2 vs CLDN4 | 516 | 0.531 | 7.19e-39 | POSITIVE |
| TCGA-LUAD | NKX2-1 vs SFTPB | 516 | 0.571 | 6.49e-46 | POSITIVE |
| TCGA-LUAD | NKX2-1 vs NAPSA | 516 | 0.661 | 3.76e-66 | POSITIVE |
| TCGA-LUAD | NKX2-1 vs KRT5 | 516 | -0.179 | 4.37e-05 | INVERSE |
| TCGA-LUAD | NKX2-1 vs CLDN7 | 516 | 0.307 | 1.07e-12 | POSITIVE |
| CPTAC-LUAD | TACSTD2 vs CLDN4 | 110 | 0.392 | 2.30e-05 | POSITIVE |
| CPTAC-LUAD | TACSTD2 vs CLDN4 | 79 | 0.275 | 1.43e-02 | POSITIVE |
| DepMap24Q4-LUAD | TACSTD2 vs CLDN4 | 80 | 0.756 | 5.11e-16 | POSITIVE |

## Coverage (honest missingness)

- TCGA-LUAD primary tumors: n=516; ABSOLUTE purity available for 503.
- CPTAC LUAD RNA: n=110; all three genes complete.
- CPTAC LUAD protein: NKX2-1 and TACSTD2 complete (n=110); CLDN4 non-NA = 79 / 110 (pairwise n=79 for NKX2-1 vs CLDN4).
- DepMap 24Q4 LUAD cell lines: n=80 (OncotreeSubtype Lung Adenocarcinoma).

## Secondary: purity-adjusted TCGA (does not change the verdict)

- NKX2-1 vs TACSTD2: partial ρ = -0.053, p = 2.32e-01, n = 503
- NKX2-1 vs CLDN4: partial ρ = 0.212, p = 1.64e-06, n = 503
- NKX2-1 vs SFTPB: partial ρ = 0.585, p = 2.46e-47, n = 503
- NKX2-1 vs NAPSA: partial ρ = 0.671, p = 7.12e-67, n = 503
- NKX2-1 vs KRT5: partial ρ = -0.100, p = 2.45e-02, n = 503
- NKX2-1 vs CLDN7: partial ρ = 0.244, p = 3.12e-08, n = 503

## What this is not

- Direct NKX2-1 binding or repression at TACSTD2/CLDN4 (no ChIP/perturbation in this slice).
- ICI response or TROP2-ADC outcome.
- TTF-1 IHC (this is continuous RNA/protein, not a clinical IHC call).
- A reason to drop CPTAC protein or DepMap because the sign disagrees with a hoped-for inverse.
- Evidence that NKX2-1 represses TACSTD2. The TACSTD2 RNA association is consistent with zero.

## Caveats

1. Bulk tumor RNA mixes epithelium, stroma, and immune cells. NKX2-1 is lineage-restricted; TACSTD2/CLDN4 are not.
2. CPTAC protein CLDN4 is missing in 31/110 tumors. The protein ρ = −0.22 (n=79) has FDR = 0.09 and a CI that includes 0. That is not an inverse call.
3. DepMap cell lines are not tumors. Both pairs are **positive** in vitro (ρ = 0.45 and 0.53).
4. Exploratory sanity checks behave as expected: NKX2-1 tracks SFTPB/NAPSA (alveolar) and is weakly inverse with KRT5 (basal). The assay is not broken; the TACSTD2/CLDN4 inverse claim is.
5. TCGA/CPTAC are not ICI cohorts.

## Rerun

```bash
python3 scripts/w200/A10_NKX21/download.py
python3 scripts/w200/A10_NKX21/analyze.py
```

See `summary.json` and `correlations.csv`.
