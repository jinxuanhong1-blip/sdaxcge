# Finding — CPTAC LUAD/LSCC CLDN4 vs CXCL9 / CXCL10 protein (if present)

**Protein vs protein.** Public CPTAC TMT freeze v1.2, **LUAD** and **LSCC** separately. **No TACSTD2 gate.** CXCL9 and CXCL10 are used only when a protein row exists. Honest pairwise-complete n (tumors with both proteins quantified). LinkedOmics NArm CCT is a presence check only: it drops CLDN4, so it cannot test these pairs.

These are treatment-naive surgical cohorts (Gillette *Cell* 2020; Satpathy *Cell* 2021). **No ICI response labels.** Do not read these correlations as immunotherapy outcomes.

## Verdict

| Pair | LUAD (freeze v1.2) | LSCC (freeze v1.2) |
|---|---|---|
| CLDN4 vs CXCL9 protein | ρ=**−0.203**; p=**0.23**; n=**37** | ρ=**−0.278**; p=**0.027**; n=**63** |
| CLDN4 vs CXCL10 protein | ρ=**0.014**; p=**0.91**; n=**71** | ρ=**−0.108**; p=**0.36**; n=**73** |

**What is present:** Both chemokines have a TMT row in freeze v1.2. This is not an absence result on the primary table.

**What holds:** In **LSCC**, CLDN4 protein is nominally inverse to CXCL9 protein (n=63). Same sign as the existing LSCC CLDN4-protein vs ImmuneScore / GEP18 / CD8A RNA slice, but this is a different layer (chemokine protein, not an RNA score).

**What does not hold:** **LUAD** is null on both pairs. **CXCL10 protein** is null in both histologies. Do not write a pan-lung “CLDN4 vs CXCL9/10 protein” claim. Do not upgrade the LSCC CXCL9 p=0.027 to a CXCL10 result.

Do not pool LUAD+LSCC TMT values. Do not use LinkedOmics NArm n (CLDN4 row is gone; pairwise n=0).

## Honest n

| | LUAD | LSCC |
|---|---:|---:|
| Tumors in freeze protein table | 110 | 108 |
| CLDN4 protein quantified | **79** (31 NA, 28%) | **78** (30 NA, 28%) |
| CXCL9 protein quantified | **45** (65 NA, 59%) | **88** (20 NA, 19%) |
| CXCL10 protein quantified | **80** (30 NA, 27%) | **98** (10 NA, 9%) |
| CLDN4 ∩ CXCL9 (tested n) | **37** | **63** |
| CLDN4 ∩ CXCL10 (tested n) | **71** | **73** |
| TACSTD2 filter | **none** | **none** |

Do not write n=110 / n=108 for the Spearman tests. LUAD CXCL9 is the thin pair: 65/110 tumors have no CXCL9 protein (TMT dropout), so the tested n is 37, not 79 and not 110.

### LinkedOmics NArm (cannot test)

| | LUAD | LSCC |
|---|---:|---:|
| CLDN4 protein row | **absent** (dropped) | **absent** (dropped) |
| CXCL9 protein | absent | 59 / 108 |
| CXCL10 protein | 55 / 110 | 88 / 108 |
| CLDN4 ∩ CXCL9 / ∩ CXCL10 | **0 / 0** | **0 / 0** |

NArm is why CLDN4 vanishes as a *row*, not as a row of NAs. Presence of CXCL10 on that table does not rescue a CLDN4–CXCL10 test.

## Methods (this slice)

- Primary source: open S3 `cptac-pancancer-data / data_freeze_v1.2_reorganized/{LUAD,LSCC}/` gene-abundance tumor TMT (already log2 reference-intensity). Filenames from the LinkedOmics CPTAC-pancan index; HEAD HTTP 200.
- Secondary source: LinkedOmics NArm CCT (gene symbols). Used only to document which proteins survive NArm.
- Predictor: CLDN4 protein `ENSG00000189143.9`. TACSTD2 is not a gate, not a covariate, and not a co-filter.
- Endpoints: CXCL9 protein `ENSG00000138755.6`; CXCL10 protein `ENSG00000169245.6`. Symbol aliases checked (`CXCL9`/`MIG`/`SCYB9`; `CXCL10`/`IP10`/`SCYB10`). No RNA substitution if protein is missing.
- Test: two-sided Spearman on pairwise-complete tumors. No quartile contrast (CXCL9 LUAD n=37 is already the limiting n).
- Cohorts kept separate. No meta-analysis. No ICI labels.

## LUAD (Gillette 2020; freeze v1.2)

110 tumors. CLDN4 protein n=**79**. CXCL9 protein n=**45**. CXCL10 protein n=**80**.

| Pair | n pairwise | ρ | p |
|---|---:|---:|---:|
| CLDN4 vs CXCL9 protein | **37** | −0.203 | 0.23 |
| CLDN4 vs CXCL10 protein | **71** | 0.014 | 0.91 |

Both null. The CXCL9 test is underpowered relative to the CLDN4-complete set (37 vs 79) because CXCL9 is missing in 59% of tumors. A null at n=37 is not proof of a zero effect, and it is not a license to borrow LSCC.

## LSCC (Satpathy 2021; freeze folder LSCC / TCGA synonym LUSC)

108 tumors. CLDN4 protein n=**78**. CXCL9 protein n=**88**. CXCL10 protein n=**98**.

| Pair | n pairwise | ρ | p |
|---|---:|---:|---:|
| CLDN4 vs CXCL9 protein | **63** | −0.278 | **0.027** |
| CLDN4 vs CXCL10 protein | **73** | −0.108 | 0.36 |

CXCL9 protein is the only nominal hit. CXCL10 protein does not follow. This is untreated resected squamous TMT, not an ICI-resistance result.

## What this does not claim

- It does not use a TACSTD2-high, TACSTD2-low, or TACSTD2-complete gate.
- It does not test CLDN4 RNA or CXCL9/CXCL10 RNA as stand-ins for protein.
- It does not test ICI response, PFS, or OS.
- It does not pool LUAD and LSCC.
- It does not treat n=110 / n=108 as the tested n.
- It does not treat LinkedOmics NArm as a negative protein result for CLDN4 (that table dropped the row).
- LUAD CXCL9 n=37 is modest. The LUAD null is underpowered relative to LSCC n=63, not a histology-proof zero.

## Outputs

- `results/n_table.tsv` — honest n by cohort and source
- `results/presence.tsv` — row present / n observed / n NA
- `results/spearman.tsv` — pairwise Spearman (or n=0 if protein missing)
- `results/sample_scores.tsv` — per-tumor protein values
- `results/summary.json`
- `results/fig_cldn4_vs_cxcl_protein.png`

```bash
python3 methods/cptac_cldn4_cxcl/download.py --outdir data/cptac_cldn4_cxcl
python3 methods/cptac_cldn4_cxcl/analyze.py --data data/cptac_cldn4_cxcl --outdir methods/cptac_cldn4_cxcl
```
