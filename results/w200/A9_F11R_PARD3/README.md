# A9 — F11R and PARD3 vs TACSTD2-high public lung

## Honest result

**F11R:** TCGA_CALL_NO_INDEPENDENT_REPLICATION. F11R meets the pre-specified TCGA call in both histologies (TCGA_LUAD ρ=0.182, TCGA_LUSC ρ=0.248), but the independent OncoSG LUAD cohort does not replicate (ρ=0.117, n=169, no-call). Treat the TCGA call as histology-detectable, not as a confirmed public-lung partner. Also: OncoSG LUAD ρ=0.117 (n=169), no-call; DepMap lung lines ρ=0.437 (n=214); pooled NSCLC ρ=0.296 (n=1017) — sensitivity only.

**PARD3:** NOT_A_HISTOLOGY_INDEPENDENT_PARTNER. PARD3 meets the call only in TCGA_LUSC and fails in TCGA_LUAD (TCGA_LUAD ρ=-0.020). A pooled LUAD+LUSC correlation is not evidence of a shared partner: TACSTD2-high is LUSC-enriched. Do not treat this gene as histology-independent. Also: OncoSG LUAD ρ=0.062 (n=169), no-call; DepMap lung lines ρ=0.187 (n=214); pooled NSCLC ρ=0.164 (n=1017) — sensitivity only.

These are bulk RNA associations in public tumors. They do not show that TROP2
binds F11R or PARD3, that either gene is required for a TROP2-high state, or
that the private cohort's effect sizes transfer.

### Primary public lung (TCGA, within-histology tertiles)

- F11R TCGA_LUAD: ρ=0.182 (n=515), tertile log2FC=+0.223, Welch FDR=0.0007, CALL, partial-CPE ρ=0.184
- F11R TCGA_LUSC: ρ=0.248 (n=502), tertile log2FC=+0.429, Welch FDR=2.7e-09, CALL, partial-CPE ρ=0.254
- PARD3 TCGA_LUAD: ρ=-0.020 (n=515), tertile log2FC=-0.046, Welch FDR=0.51, no-call, partial-CPE ρ=-0.025
- PARD3 TCGA_LUSC: ρ=0.151 (n=502), tertile log2FC=+0.255, Welch FDR=0.0024, CALL, partial-CPE ρ=0.152

### Independent LUAD (OncoSG 2020)

- F11R OncoSG_LUAD: ρ=0.117 (n=169), tertile log2FC=+0.332, Welch FDR=0.17, no-call, partial-CPE ρ=0.051
- PARD3 OncoSG_LUAD: ρ=0.062 (n=169), tertile log2FC=+0.237, Welch FDR=0.2, no-call, partial-CPE ρ=0.037

### Stroma-free sensitivity (DepMap 24Q4 lung cell lines)

- F11R DepMap24Q4_lung_cell_lines: ρ=0.437 (n=214), tertile log2FC=+1.586, Welch FDR=1e-09, CALL
- PARD3 DepMap24Q4_lung_cell_lines: ρ=0.187 (n=214), tertile log2FC=+0.326, Welch FDR=0.012, CALL

### Context (not the claim)

CLDN1 and CLDN4 are stronger TACSTD2 correlates than F11R or PARD3 in every
bulk cohort tested (e.g. TCGA-LUAD CLDN4 ρ=0.46 vs F11R ρ=0.18; OncoSG CLDN4
ρ=0.50 vs F11R ρ=0.12). A TROP2-high tight-junction story that needs F11R or
PARD3 as load-bearing public evidence is not supported at the same strength
as CLDN1/CLDN4.

DepMap lung cell lines (no stroma) give F11R ρ=0.44, so the weak LUAD bulk
signal is not proof that F11R is only a stromal artifact. It still does not
make F11R a replicated bulk-tumor partner.

### What this is not

- Not a protein / IHC / spatial result.
- Not a malignant-cell-only result (except the optional DepMap cell-line slice).
- Not evidence that F11R or PARD3 belong to a TROP2-high *program* on their
  own: CLDN1/CLDN4 are stronger public correlates; PARD3 is in the weak tail.
- Pooled LUAD+LUSC is **not** the primary test. TACSTD2-high tumors are
  LUSC-enriched, so an unstratified pool can inflate a junction-gene signal.

## Pre-specified methods

| Item | Choice |
|------|--------|
| Claim genes | F11R, PARD3 |
| Primary cohorts | TCGA-LUAD and TCGA-LUSC primary tumors (`-01`), separate |
| TACSTD2-high | Within-cohort upper tertile; lower tertile = low |
| Primary association | Spearman ρ vs continuous TACSTD2 |
| High-vs-low test | Welch t-test on log2 expression |
| Call | ρ>0 AND log2FC>0 AND Welch BH-FDR<0.05, FDR across the 2 claim genes |
| Purity | Partial Spearman vs Aran CPE / ESTIMATE / ABSOLUTE |
| Independent | OncoSG LUAD 2020 (cBioPortal z-scored RNA; Spearman is invariant) |
| Stroma-free | DepMap 24Q4 lung cell lines, if downloaded |
| Comparators | CLDN1/4/7 (other user TJ genes), EPCAM/KRT8, PTPRC/CD8A |

Median split, pooled NSCLC, and DepMap are labeled sensitivity. They were not
used to flip a primary call.

## Sources

- UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` and `TCGA.LUSC.sampleMap/HiSeqV2` (log2(norm_count+1))
- Aran, Sirota & Butte, *Nat Commun* 2015, Supplementary Data 1 (purity)
- OncoSG LUAD (Chen et al. 2020) via [cBioPortal](https://www.cbioportal.org/study/summary?id=luad_oncosg_2020)
- DepMap Public 24Q4 (optional): https://doi.org/10.25452/figshare.plus.27993248.v1

## Reproduce

```bash
python3 scripts/w200/A9_F11R_PARD3/download.py
python3 scripts/w200/A9_F11R_PARD3/analyze.py
```
