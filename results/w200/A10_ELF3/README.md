# A10 ELF3 public ChIP / co-expression vs TACSTD2 / CLDN4 in lung

**Claim:** ELF3 binds TACSTD2 and CLDN4 in lung (public ChIP) and is co-expressed with them in public lung tumors.

**Honest verdict: `CHIP_NOT_TESTABLE_IN_LUNG`**

Public ELF3 ChIP does not exist in lung, so a lung binding claim at TACSTD2/CLDN4 is not testable. Non-lung ChIP (pancreas/liver/esophagus/biliary) is mixed and is not a lung result. TCGA-LUAD primary tumors: ELF3 vs TACSTD2 POSITIVE (ρ=0.369 (n=516, FDR=9.8e-18)); ELF3 vs CLDN4 POSITIVE (ρ=0.528 (n=516, FDR=2.3e-37)). TCGA-LUSC: TACSTD2 POSITIVE (ρ=0.464 (n=501, FDR=1.4e-27)); CLDN4 POSITIVE (ρ=0.500 (n=501, FDR=2.7e-32)). Positive bulk co-expression is expected for epithelial genes and is not evidence that ELF3 binds or regulates TACSTD2/CLDN4 in lung.

Filters were not tuned to force a sign. Motif scanning was not run.

## 1. Public ELF3 ChIP (lung is absent)

No public human lung / NSCLC ELF3 ChIP-seq was found in ENCODE, ChIP-Atlas (ELF3.5 target matrix), or a GEO ELF3+ChIP-seq survey. ELF3 binding at TACSTD2/CLDN4 therefore cannot be tested in lung from public ChIP. Available non-lung ChIP (CFPAC-1 pancreas, HepG2 liver, ESO-26 esophagus, HBDEC2 biliary) is mixed and is not a lung result.

Catalog (all `is_lung=false` except the mouse false-positive GEO text hit): `chip_experiments.csv`. Locus overlaps: `chip_locus_overlap.csv`.

Non-lung overlap calls (not a lung result):

| Call | Result |
| --- | --- |
| TACSTD2_coding_TSS_pm2kb_CFPAC1 | True |
| TACSTD2_coding_TSS_pm2kb_HepG2_ENCODE_IDR | False |
| TACSTD2_coding_TSS_pm2kb_ESO26_or_HBDEC2 | False |
| CLDN4_coding_TSS_pm2kb_CFPAC1 | True |
| CLDN4_alt_TSS_pm2kb_CFPAC1 | True |
| CLDN4_alt_TSS_pm2kb_HepG2_ENCODE_IDR | True |
| CLDN4_coding_TSS_pm2kb_HepG2_ENCODE_IDR | False |

What the non-lung peaks actually are:

- **CFPAC-1 (PDAC):** MACS2 bed05 peaks recur near the TACSTD2 coding TSS (±2 kb) in several replicates, and along the CLDN4 gene span / alt TSS. A peak at the CLDN4 *coding* TSS (±2 kb) is rare (one CFPAC-1 library).
- **HepG2 (ENCODE optimal IDR ENCFF080FAU):** no peak at TACSTD2 gene ±10 kb. One peak at chr7:73799582–73799998, which is the CLDN4 *alt* TSS (~31 kb upstream of the protein-coding TSS), not the coding promoter.
- **ESO-26 and HBDEC2:** no bed05 peak at either locus. HBDEC2 is biliary epithelium (GSE156165), not lung.
- **HEK293 (GSE280165):** additional 2024 ELF3 ChIP, not lung; peak bed not pulled.

## 2. Public lung co-expression (primary pairs, BH-FDR within this list)

| Cohort | Layer | Pair | n | Spearman ρ | p | FDR | 95% CI | Direction |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- |
| TCGA-LUAD | RNA_log2TPM | ELF3 vs TACSTD2 | 516 | 0.369 | 3.93e-18 | 9.84e-18 | [0.292, 0.443] | POSITIVE |
| TCGA-LUAD | RNA_log2TPM | ELF3 vs CLDN4 | 516 | 0.528 | 2.27e-38 | 2.27e-37 | [0.456, 0.595] | POSITIVE |
| TCGA-LUSC | RNA_log2TPM | ELF3 vs TACSTD2 | 501 | 0.464 | 4.16e-28 | 1.39e-27 | [0.388, 0.533] | POSITIVE |
| TCGA-LUSC | RNA_log2TPM | ELF3 vs CLDN4 | 501 | 0.500 | 5.49e-33 | 2.74e-32 | [0.423, 0.568] | POSITIVE |
| CPTAC-LUAD | RNA_RSEM_UQ_log2 | ELF3 vs TACSTD2 | 110 | 0.358 | 1.23e-04 | 0.000175 | [0.187, 0.510] | POSITIVE |
| CPTAC-LUAD | RNA_RSEM_UQ_log2 | ELF3 vs CLDN4 | 110 | 0.276 | 3.58e-03 | 0.00447 | [0.084, 0.452] | POSITIVE |
| CPTAC-LUAD | protein_TMT_log2 | ELF3 vs TACSTD2 | 110 | 0.177 | 6.41e-02 | 0.0712 | [-0.020, 0.363] | NULL |
| CPTAC-LUAD | protein_TMT_log2 | ELF3 vs CLDN4 | 79 | 0.132 | 2.45e-01 | 0.245 | [-0.093, 0.334] | NULL |
| DepMap24Q4-LUAD | RNA_log2TPM | ELF3 vs TACSTD2 | 80 | 0.636 | 2.38e-10 | 4.76e-10 | [0.481, 0.752] | POSITIVE |
| DepMap24Q4-LUAD | RNA_log2TPM | ELF3 vs CLDN4 | 80 | 0.609 | 2.12e-09 | 3.53e-09 | [0.429, 0.750] | POSITIVE |

## Exploratory (not used for the ChIP verdict)

| Cohort | Pair | n | Spearman ρ | p | Direction |
| --- | --- | ---: | ---: | ---: | --- |
| TCGA-LUAD | TACSTD2 vs CLDN4 | 516 | 0.531 | 7.19e-39 | POSITIVE |
| TCGA-LUAD | ELF3 vs EPCAM | 516 | 0.277 | 1.53e-10 | POSITIVE |
| TCGA-LUAD | ELF3 vs KRT8 | 516 | 0.246 | 1.49e-08 | POSITIVE |
| TCGA-LUAD | ELF3 vs KRT5 | 516 | 0.021 | 6.37e-01 | NULL |
| TCGA-LUAD | ELF3 vs NKX2-1 | 516 | 0.198 | 6.08e-06 | POSITIVE |
| TCGA-LUSC | TACSTD2 vs CLDN4 | 501 | 0.391 | 8.59e-20 | POSITIVE |
| TCGA-LUSC | ELF3 vs EPCAM | 501 | 0.249 | 1.68e-08 | POSITIVE |
| TCGA-LUSC | ELF3 vs KRT8 | 501 | 0.326 | 6.60e-14 | POSITIVE |
| TCGA-LUSC | ELF3 vs KRT5 | 501 | 0.083 | 6.43e-02 | NULL |
| TCGA-LUSC | ELF3 vs NKX2-1 | 501 | 0.092 | 4.00e-02 | POSITIVE |
| CPTAC-LUAD | TACSTD2 vs CLDN4 | 110 | 0.392 | 2.30e-05 | POSITIVE |
| CPTAC-LUAD | TACSTD2 vs CLDN4 | 79 | 0.275 | 1.43e-02 | POSITIVE |
| DepMap24Q4-LUAD | TACSTD2 vs CLDN4 | 80 | 0.756 | 5.11e-16 | POSITIVE |
| DepMap24Q4-Lung | ELF3 vs TACSTD2 | 214 | 0.541 | 1.10e-17 | POSITIVE |
| DepMap24Q4-Lung | ELF3 vs CLDN4 | 214 | 0.696 | 2.35e-32 | POSITIVE |
| DepMap24Q4-Lung | TACSTD2 vs CLDN4 | 214 | 0.607 | 5.64e-23 | POSITIVE |

## Coverage

- TCGA-LUAD primary tumors: n=516; ABSOLUTE purity for 503.
- TCGA-LUSC primary tumors: n=501; ABSOLUTE purity for 493.
- CPTAC LUAD RNA: n=110.
- CPTAC LUAD protein non-NA: ELF3=110, TACSTD2=110, CLDN4=79 (CLDN4 pairwise n=79).
- DepMap 24Q4 LUAD cell lines: n=80.
- DepMap 24Q4 all lung cell lines (exploratory): n=214.

## Secondary: purity-adjusted TCGA (does not create lung ChIP)

- TCGA-LUAD ELF3 vs TACSTD2: partial ρ = 0.358, p = 1.30e-16, n = 503
- TCGA-LUAD ELF3 vs CLDN4: partial ρ = 0.521, p = 2.72e-36, n = 503
- TCGA-LUAD ELF3 vs EPCAM: partial ρ = 0.198, p = 7.47e-06, n = 503
- TCGA-LUAD ELF3 vs KRT8: partial ρ = 0.231, p = 1.71e-07, n = 503
- TCGA-LUAD ELF3 vs KRT5: partial ρ = 0.103, p = 2.11e-02, n = 503
- TCGA-LUAD ELF3 vs NKX2-1: partial ρ = 0.133, p = 2.76e-03, n = 503
- TCGA-LUSC ELF3 vs TACSTD2: partial ρ = 0.476, p = 3.98e-29, n = 493
- TCGA-LUSC ELF3 vs CLDN4: partial ρ = 0.489, p = 5.75e-31, n = 493
- TCGA-LUSC ELF3 vs EPCAM: partial ρ = 0.183, p = 4.44e-05, n = 493
- TCGA-LUSC ELF3 vs KRT8: partial ρ = 0.274, p = 6.20e-10, n = 493
- TCGA-LUSC ELF3 vs KRT5: partial ρ = 0.069, p = 1.25e-01, n = 493
- TCGA-LUSC ELF3 vs NKX2-1: partial ρ = 0.155, p = 5.69e-04, n = 493

## What this is not

- ELF3 ChIP-seq in human lung / NSCLC / A549 (none public in ENCODE, ChIP-Atlas, or GEO survey).
- Whether ELF3 directly transactivates TACSTD2 or CLDN4 in lung (no lung ChIP, no ELF3 KD/KO lung RNA in this slice).
- ICI response or TROP2-ADC outcome.
- A motif presence/absence argument (not run; user asked ChIP/co-expression only).
- Evidence that CFPAC-1 or HepG2 peaks transfer to lung epithelium.

## Caveats

1. Bulk tumor RNA mixes epithelium, stroma, and immune cells. ELF3, TACSTD2, and CLDN4 are all epithelial-leaning; a positive ρ can be purity / epithelial fraction.
2. Partial Spearman on ABSOLUTE purity is a check, not a cell-type deconvolution.
3. DepMap lines are not tumors.
4. ChIP-Atlas gene scores are MACS2 peak scores assigned to a 5 kb gene window; a non-zero score is not the same as a coding-TSS peak.
5. CLDN4's Ensembl gene span starts ~31 kb upstream of the protein-coding TSS. Peaks at the alt TSS are not automatically coding-promoter binding.
6. TCGA/CPTAC are not ICI cohorts.

## Rerun

```bash
python3 scripts/w200/A10_ELF3/chip.py
python3 scripts/w200/A10_ELF3/download.py
python3 scripts/w200/A10_ELF3/analyze.py
```

See `summary.json`, `chip_verdict.json`, and `correlations.csv`.
