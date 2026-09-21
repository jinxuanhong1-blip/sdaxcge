# Concordant-4 malignant NHEJ, STING/IFN, and MHC-I

ADDITIVE. Observational only. This is not a causal mediation and not a
re-derivation of the locked concordant-4 result (malignant CLDN4 %pos vs
T/NK, ρ = −0.531, n = 65). Same four cohorts and the same patient / donor /
sample units: **GSE123902 + GSE131907 + GSE205335 + GSE189357**. Not
GSE148071, GSE127465, GSE154826, GSE200563, or E-MTAB-13526.

Malignant cells only. GSE123902 and GSE189357 use the marker gate
(EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0. GSE131907 and GSE205335 use author
malignant labels. Scores are patient pseudobulk means of log2(CPM+1), CPM from
that unit's malignant UMI total. p-values are descriptive.

## Modules

- **NHEJ**: MSigDB KEGG_NON_HOMOLOGOUS_END_JOINING (hsa03450):
  DCLRE1C, DNTT, FEN1, LIG4, MRE11, NHEJ1, POLL, POLM, PRKDC, RAD50, XRCC4, XRCC5, XRCC6.
  Sensitivity drops DNTT and FEN1.
- **STING core** (sensor/adapter, not the ISG program): CGAS, STING1, TBK1, IRF3, IFI16, DDX41.
  NHEJ genes that sit in Reactome STING lists (XRCC5, XRCC6, PRKDC, MRE11) are excluded.
  Negative regulators (TREX1, NLRC3, NLRP4) are excluded. IRF7 and ZBP1 are left in the IFN arm.
- **IFN**: Hallmark IFNα ∪ IFNγ, the same union as the concordant-4 Seurat run.
- **STING/IFN module**: equal-weight mean of the within-cohort z-scored STING core and IFN scores.
- **MHC-I**: the 21-gene MHC-I/APM list used in that Seurat run. Twelve of those genes also sit in Hallmark IFN-γ; a sensitivity IFN score drops them.

Aliases collapsed before scoring: MRE11A→MRE11, TMEM173→STING1, MB21D1→CGAS, XLF→NHEJ1, MARCH1→MARCHF1, WARS→WARS1.
GSE205335 stores STING1 as TMEM173; that alias is applied before the STING score.
Genes absent from a dataset are left out of that dataset's mean. Zeros stay in.

## Honest n

- Units scored: **n = 65** (123902 13, 131907 21, 205335 22, 189357 9).
- Malignant-cell gate check against the locked table: CLDN4 %pos Spearman ρ = 1.000, max |Δ| = 0.0000 percentage points, n matched = 65.
- Locked IFN score vs this Hallmark IFN score: Spearman ρ = 0.960 (same units, different normalization: theirs is TMM-CPM).
- Do not quote malignant-cell or total-cell counts as n.

## 1. Patient-level association of CLDN4 with each module

Primary CLDN4 measure is malignant **% positive**, the locked concordant-4 patient score.
Module scores stay on log2(CPM+1). Pseudobulk CLDN4 abundance is a sensitivity, because
it does not rank the same patients (cohort-adjusted partial ρ of %pos vs log2(CPM+1) = 0.525).
Two summaries of each association:

1. **Partial Spearman controlling for cohort**: rank CLDN4 and the module across the 65 units, residualize both ranks on cohort indicators, then Pearson-correlate the residuals. Degrees of freedom subtract the three cohort indicators.
2. **Within-cohort Spearman**, then DerSimonian–Laird on the Fisher z values (same pooling as the locked T/NK meta).

| module | partial ρ (cohort) | p | DL ρ | DL p | I² | DL 95% CI |
|---|---:|---:|---:|---:|---:|---|
| CLDN4 %pos vs NHEJ | -0.032 | 0.8058 | 0.038 | 0.7794 | 0% | -0.227 to 0.298 |
| CLDN4 %pos vs STING/IFN | -0.244 | 0.0557 | -0.275 | 0.1400 | 44% | -0.576 to 0.092 |
| CLDN4 %pos vs MHC-I | -0.316 | 0.0123 | -0.229 | 0.0903 | 0% | -0.464 to 0.037 |

Benjamini–Hochberg q across the three primary partial tests: NHEJ q=0.8058, STING/IFN q=0.0836, MHC-I q=0.0368.

Within-cohort Spearmans (context, not a second n):

| module | 123902 | 131907 | 205335 | 189357 |
|---|---:|---:|---:|---:|
| NHEJ | 0.302 (n=13) | -0.040 (n=21) | 0.034 (n=22) | -0.167 (n=9) |
| STING/IFN | -0.374 (n=13) | 0.171 (n=21) | -0.496 (n=22) | -0.450 (n=9) |
| MHC-I | -0.110 (n=13) | 0.006 (n=21) | -0.475 (n=22) | -0.250 (n=9) |

Supporting rows (not the three primary tests):

| contrast | partial ρ | p | DL ρ | DL p |
|---|---:|---:|---:|---:|
| CLDN4 %pos vs STING core | -0.196 | 0.1268 | -0.247 | 0.2436 |
| CLDN4 %pos vs Hallmark IFN | -0.280 | 0.0276 | -0.254 | 0.0586 |
| CLDN4 %pos vs IFN minus MHC genes | -0.258 | 0.0426 | -0.250 | 0.0631 |
| CLDN4 %pos vs NHEJ core | -0.099 | 0.4430 | -0.010 | 0.9441 |
| CLDN4 logCPM vs NHEJ | 0.071 | 0.5826 | -0.008 | 0.9562 |
| CLDN4 logCPM vs STING/IFN | 0.028 | 0.8263 | 0.015 | 0.9142 |
| CLDN4 logCPM vs MHC-I | 0.059 | 0.6486 | 0.077 | 0.5734 |
| CLDN4 logCPM vs Hallmark IFN | 0.053 | 0.6844 | 0.049 | 0.7604 |

## 2. Does NHEJ-low go with IFN-high inside CLDN4 strata?

This is the mediation-style question. It asks whether the NHEJ–IFN association
is still there after CLDN4 is held fixed. It does not identify a mechanism.

- Partial Spearman, NHEJ vs Hallmark IFN, controlling for CLDN4 %pos and cohort: **ρ = 0.190**, p = 0.1434, n = 65.
- Same partial, controlling for CLDN4 log2(CPM+1) instead of %pos: ρ = 0.188, p = 0.1472.
- Same partial, additionally controlling for a 6-gene proliferation score: ρ = 0.119, p = 0.3631.
- Same partial using the NHEJ-core sensitivity set: ρ = 0.189, p = 0.1437.
- Partial Spearman, NHEJ vs STING core, controlling for CLDN4 and cohort: ρ = 0.146, p = 0.2607.

Within-cohort median split on CLDN4 %pos. NHEJ and IFN are within-cohort z-scores, then pooled. Negative ρ means NHEJ-low with IFN-high inside that stratum.

| CLDN4 stratum | n | Spearman NHEJ vs IFN | p |
|---|---:|---:|---:|
| CLDN4 low (within-cohort bottom half) | 31 | 0.168 | 0.3659 |
| CLDN4 high (within-cohort top half) | 34 | 0.199 | 0.2582 |
| CLDN4 %pos T1 | 21 | 0.345 | 0.1251 |
| CLDN4 %pos T2 | 21 | -0.099 | 0.6704 |
| CLDN4 %pos T3 | 23 | 0.278 | 0.1996 |

Within-cohort partial Spearman (NHEJ vs IFN | CLDN4), DL meta: ρ = 0.138, p = 0.3113, I² = 0%, k = 4.

Cohort-specific partials (small n; GSE189357 is 9 units):

| cohort | n | partial ρ | p |
|---|---:|---:|---:|
| GSE123902 | 13 | 0.229 | 0.4739 |
| GSE131907 | 21 | -0.089 | 0.7103 |
| GSE205335 | 22 | 0.303 | 0.1818 |
| GSE189357 | 9 | 0.115 | 0.7861 |

## 3. Observational path split (not an effect)

Within-cohort z-scores of CLDN4 %pos, NHEJ, and IFN. OLS. Path a is CLDN4 with NHEJ. Path b is NHEJ with IFN when CLDN4 is in the model. Path c is CLDN4 with IFN. Path c′ is CLDN4 with IFN when NHEJ is in the model. The product a×b is a descriptive split of the CLDN4–IFN association, not a causal indirect effect. No sequential ignorability, no intervention.

| path | beta | p | n |
|---|---:|---:|---:|
| a. NHEJ ~ CLDN4 %pos | -0.024 | 0.8495 | 65 |
| c. IFN ~ CLDN4 %pos | -0.259 | 0.0370 | 65 |
| c′. IFN ~ CLDN4 %pos, NHEJ also in the model | -0.257 | 0.0393 | 65 |
| b. IFN ~ NHEJ, CLDN4 %pos also in the model | 0.107 | 0.3851 | 65 |
| b, proliferation also in the model | 0.028 | 0.8574 | 65 |
| STING core ~ NHEJ, CLDN4 %pos also in the model | 0.170 | 0.1692 | 65 |

Stratified bootstrap of a×b (2000 resamples within cohort): -0.003 (percentile interval -0.066 to 0.063; fraction of draws below 0 = 0.53). The interval is uncertainty of this decomposition, not a causal confidence interval.

## Reading

Cohort-adjusted partial correlations on CLDN4 %pos: CLDN4 %pos vs NHEJ is near null (ρ = -0.032); CLDN4 %pos vs STING/IFN is negative (ρ = -0.244); CLDN4 %pos vs MHC-I is negative (ρ = -0.316). DerSimonian–Laird intervals for STING/IFN and MHC-I still include 0. STING/IFN I² is 44%: GSE131907 is positive and the other three cohorts are negative. NHEJ vs IFN | CLDN4, cohort is near null (ρ = 0.190). Inside the CLDN4-low half, NHEJ vs IFN ρ = 0.168 (n = 31); inside the CLDN4-high half, ρ = 0.199 (n = 34). Both stratum estimates are positive. The cytosolic-DNA pattern would be a negative partial (NHEJ-low with IFN-high at a fixed CLDN4 rank). That pattern is not what these 65 units show.

Abundance sensitivity, same partial Spearman on malignant CLDN4 log2(CPM+1): CLDN4 logCPM vs NHEJ ρ = 0.071; CLDN4 logCPM vs STING/IFN ρ = 0.028; CLDN4 logCPM vs MHC-I ρ = 0.059; CLDN4 logCPM vs Hallmark IFN ρ = 0.053. The locked Q4-versus-Q1 IFN/MHC shift used CLDN4 %pos, not this abundance rank.

Cross-sectional patient pseudobulks. Cytosolic-DNA logic is why NHEJ and STING are on the table. The split above is an association decomposition. STING1 is on the matrix in every cohort after the TMEM173 alias, but the sensor arm is a 6-gene mean next to a 224-gene IFN mean, so the composite is not a pure STING score (the STING-core row is separate). Twelve MHC-I genes also sit in Hallmark IFN-γ, so the MHC-I and IFN rows share genes. The IFN-minus-MHC sensitivity stays negative on CLDN4 %pos.

## Reproduce

```
bash methods/nhej_sting_concordant4/scripts/download.sh /tmp/geo_nhej
python3 methods/nhej_sting_concordant4/scripts/extract_sums.py
Rscript methods/nhej_sting_concordant4/scripts/export_gse205335.R /tmp/geo_nhej /tmp/geo_nhej/gse205335_export
python3 methods/nhej_sting_concordant4/scripts/analyze_associations.py
```

