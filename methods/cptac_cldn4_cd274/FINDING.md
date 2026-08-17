# Finding — CPTAC LUAD / LSCC protein: CLDN4 vs CD274 and IFN proteins

**Additive only.** Prior CPTAC CLDN4 **protein vs RNA-layer** ImmuneScore / GEP18 / CD8A (PR [#315](https://github.com/jinxuanhong1-blip/sdaxcge/pull/315) Q4; PR [#289](https://github.com/jinxuanhong1-blip/sdaxcge/pull/289) LSCC; PR [#245](https://github.com/jinxuanhong1-blip/sdaxcge/pull/245) LUAD TJ) is taken as given and is **not** re-estimated. This page adds **protein–protein** CLDN4 versus **CD274** and a locked IFN list on the same public TMT freeze.

Public CPTAC TMT freeze v1.2. **LUAD** (Gillette *Cell* 2020) and **LSCC** (Satpathy *Cell* 2021) kept separate. Tumor protein matrix only. Pairwise-complete Spearman ρ, two-sided p, 2,000-resample bootstrap 95% CI (seed `20260817`). Honest n. Do not write n=110 / n=108 for these tests.

Treatment-naive surgical proteomes. **No ICI labels.** CD274 here is **protein**, not PD-L1 IHC and not RNA.

## Verdict

| Endpoint | LUAD | LSCC |
|---|---|---|
| **CLDN4 vs CD274 protein (primary)** | n=**79** ρ=**+0.056** p=**0.63** | n=**78** ρ=**−0.098** p=**0.39** |
| IFN ligands (IFNG, IFNA1, IFNB1) | **absent** | **absent** |
| IFN-core protein score (10/10 signaling+ISG) | n=79 ρ=−0.192 p=0.089 | n=78 ρ=−0.175 p=0.12 |
| MHC-I protein score (HLA-A/B/C; B2M absent) | n=79 ρ=+0.049 p=0.67 | n=78 ρ=−0.410 p=**2.0×10⁻⁴** |

**What holds:** CD274 protein is quantified in every tumor (LUAD 110/110, LSCC 108/108). The pairwise n is the CLDN4-complete set (**79** / **78**). CLDN4 protein and CD274 protein **do not co-vary** in either histology. Both 95% CIs include 0.

**What does not hold:** There is no CPTAC protein support that high CLDN4 protein is high (or low) PD-L1 protein. IFN **ligands** are not in this TMT table — do not invent IFNG / IFNA / IFNB protein correlations. The 10-gene IFN-core score is a weak negative trend and is **not** significant in either cohort.

LSCC has a separate, extra inverse between CLDN4 protein and HLA-A/B/C and some IFN receptors / JAKs. That is **not** a CD274 result and is **not** an ISG result (ISG15, MX1, OAS1, IFI35 are null in LSCC). Do not upgrade the LSCC MHC-I extra to a PD-L1 claim.

Do not pool LUAD+LSCC TMT values. Histology is not interchangeable here.

## Honest n

| | LUAD | LSCC |
|---|---:|---:|
| Tumors in protein matrix | 110 | 108 |
| CLDN4 protein quantified | **79** | **78** |
| CLDN4 protein NA (TMT dropout) | **31 (28%)** | **30 (28%)** |
| CD274 protein quantified | **110 / 110** | **108 / 108** |
| CLDN4 ∩ CD274 (this test) | **79** | **78** |
| IFNG / IFNA1 / IFNB1 protein | **0 / 0 / 0** | **0 / 0 / 0** |
| B2M protein | **0** | **0** |
| IFNAR2 protein ∩ CLDN4 | 0 (row absent) | **0** (row in 10 tumors, none with CLDN4) |

CLDN4 is `ENSG00000189143.9`. CD274 is `ENSG00000120217.14`. Missing CLDN4 is real TMT dropout, not a join error. CD274 is not the limiting n.

## Primary: CLDN4 protein vs CD274 protein

| Cohort | n | ρ | 95% CI | p |
|---|---:|---:|---|---:|
| LUAD | 79 | +0.056 | −0.176 to +0.282 | 0.63 |
| LSCC | 78 | −0.098 | −0.312 to +0.140 | 0.39 |

Extra scatter: `figures/fig_cldn4_vs_cd274.png`.

## IFN proteins if present

Locked list. Not a search. Ligands first, then receptors, signaling, ISGs. Absent or pairwise n<8 → no ρ.

### Ligands — not quantified

| Protein | LUAD | LSCC |
|---|---|---|
| IFNG | absent | absent |
| IFNA1 | absent | absent |
| IFNB1 | absent | absent |

### Receptors

| Protein | LUAD n, ρ, p | LSCC n, ρ, p |
|---|---|---|
| IFNAR1 | 61, −0.156, 0.23 | 45, −0.422, **0.0039** |
| IFNAR2 | absent | present in 10 tumors; pairwise n=**0** |
| IFNGR1 | 79, −0.007, 0.95 | 78, −0.346, **0.0019** |
| IFNGR2 | absent | absent |

### Signaling

| Protein | LUAD n, ρ, p | LSCC n, ρ, p |
|---|---|---|
| STAT1 | 79, −0.218, 0.054 | 78, −0.254, **0.025** |
| STAT2 | 79, −0.136, 0.23 | 78, −0.094, 0.41 |
| IRF1 | 67, −0.107, 0.39 | 65, −0.219, 0.079 |
| IRF9 | 79, −0.283, **0.011** | 78, −0.281, **0.013** |
| JAK1 | 79, −0.036, 0.75 | 78, −0.369, **8.8×10⁻⁴** |
| JAK2 | 79, −0.119, 0.30 | 78, −0.419, **1.3×10⁻⁴** |

### ISGs

| Protein | LUAD n, ρ, p | LSCC n, ρ, p |
|---|---|---|
| ISG15 | 79, −0.055, 0.63 | 78, +0.010, 0.93 |
| MX1 | 79, −0.037, 0.75 | 78, −0.036, 0.76 |
| OAS1 | 79, −0.183, 0.11 | 78, +0.056, 0.62 |
| IFI35 | 79, −0.128, 0.26 | 78, +0.015, 0.90 |

IFN-core score = mean of per-gene z-scores on the 10 signaling+ISG members that pass n≥8 (all 10 in both cohorts). LUAD ρ=−0.192 (n=79, p=0.089, CI −0.423 to +0.032). LSCC ρ=−0.175 (n=78, p=0.12, CI −0.394 to +0.055). The LSCC JAK/STAT inverses are real on those rows; the ISGs pull the mean-z score back to a null.

## Extra: MHC-I proteins (IFN axis, not IFN ligands)

B2M is **absent**. Score uses HLA-A, HLA-B, HLA-C (3/4).

| Protein | LUAD n, ρ, p | LSCC n, ρ, p |
|---|---|---|
| HLA-A | 79, −0.101, 0.37 | 78, −0.332, **0.0029** |
| HLA-B | 79, +0.039, 0.73 | 78, −0.227, **0.045** |
| HLA-C | 79, +0.071, 0.54 | 78, −0.344, **0.0020** |
| B2M | absent | absent |
| MHC-I score (A/B/C) | 79, +0.049, 0.67 | 78, −0.410, **2.0×10⁻⁴** |

LSCC MHC-I inverse is the strongest extra on this page. It is still not CD274 and not an IFN ligand. LUAD MHC-I is null.

Extra IFN-axis scatter: `figures/fig_cldn4_vs_ifn_proteins.png`.

## What this does not claim

- It does not re-fit CLDN4 protein vs ImmuneScore / GEP18 / CD8A RNA.
- It does not test PD-L1 IHC, RNA, or ICI response.
- It does not test IFN-γ **ligand** protein (row absent).
- It does not treat the LSCC JAK / IFNGR1 / HLA inverse as a CD274 or ISG finding.
- It does not pool LUAD and LSCC.
- It does not treat n=110 / n=108 as the CLDN4–CD274 n.
- DepMap lung **RNA** CLDN4 vs CD274 was positive (PR [#304](https://github.com/jinxuanhong1-blip/sdaxcge/pull/304)). That is cell-line RNA, not this tumor TMT protein pair.

## Methods (this slice)

- **Source:** open S3 `cptac-pancancer-data / data_freeze_v1.2_reorganized/{LUAD,LSCC}/`. Filenames from the LinkedOmics CPTAC-pancan index; HEAD HTTP 200.
- **Files:** `{LUAD,LSCC}_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt`. Already log2 vs reference; not logged again.
- **IDs:** Ensembl prefix match (`ENSG…` ± version). Locked map in `genes.py`.
- **Predictor:** CLDN4 protein only. No TACSTD2 gate.
- **Primary endpoint:** CD274 protein.
- **IFN list (locked):** ligands IFNG, IFNA1, IFNB1; receptors IFNAR1, IFNAR2, IFNGR1, IFNGR2; signaling STAT1, STAT2, IRF1, IRF9, JAK1, JAK2; ISGs ISG15, MX1, OAS1, IFI35. Extra MHC-I: HLA-A, HLA-B, HLA-C, B2M.
- **Scores:** mean of per-gene z-scores; IFN-core requires ≥4 members; MHC-I requires ≥3 (B2M is missing).
- **Test:** Spearman, two-sided, pairwise-complete. Bootstrap 95% CI: 2,000 resamples, seed 20260817. Minimum usable n = 8.
- **Not done:** no RNA endpoints; no ImmuneScore re-fit; no ICI model; no gene-set fishing; no LUAD+LSCC pool.

## Outputs

- `tables/n_table.tsv` — honest n
- `tables/coverage.tsv` — present / absent / pairwise n
- `tables/correlations.tsv` — n, ρ, p, CI
- `tables/sample_scores.tsv`
- `tables/summary.json`
- `figures/fig_cldn4_vs_cd274.png` — extra scatter (primary)
- `figures/fig_cldn4_vs_ifn_proteins.png` — extra scatter (present IFN-axis proteins)

```bash
python3 -m pip install -r methods/cptac_cldn4_cd274/requirements.txt
python3 methods/cptac_cldn4_cd274/download.py --outdir data/cptac_cldn4_cd274
python3 methods/cptac_cldn4_cd274/analyze.py --data data/cptac_cldn4_cd274 --outdir methods/cptac_cldn4_cd274
```
