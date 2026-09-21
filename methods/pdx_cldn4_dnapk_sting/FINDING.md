# Finding — public NSCLC PDX: CLDN4 protein vs DNA-PK subunits and STING-pathway proteins

**Protein vs protein, human tumor rows only.** Mirhadi et al., *Nat Commun* 2022 (doi:10.1038/s41467-022-29444-9), Supplementary Data 1, sheet `Normalized log2 Prot Quant`. PRIDE **PXD016579**. This is the public TMT proteome of NSCLC patient-derived xenografts (human tumor vs mouse stroma separated by species-unique peptides). It is not the private PDX proteome, not CPTAC primary-tumor TMT, and not an ICI cohort.

LUAD and LUSC are the primary cohorts. Technical replicate pairs were averaged to one model before testing. Spearman uses pairwise-complete log2 values. CLDN4 missingness is not filled in as zero.

## Verdict

Human CLDN4 protein is quantified in a minority of models. Among those models, the three DNA-PK subunits are **positive in LUAD and near zero in LUSC**, and none clear a nominal 0.05 or the within-histology BH threshold. The four STING-core proteins are **mixed in sign** in both histologies. A LUAD+LUSC pool (n=47) does not create a correlation that the separate histologies lack.

**What is present:** Human CLDN4 (O14493), PRKDC (P78527), XRCC5/Ku80 (P13010), XRCC6/Ku70 (P12956), STING1 (Q86WV6), CGAS (Q8N884), TBK1 (Q9UHD2), and IRF3 (Q14653) all have a human row. DNA-PK subunits and TBK1 are quantified in every model. The limiting count is CLDN4.

**What the coefficients are:** LUAD DNA-PK score ρ=0.357, p=0.087, n=24 (Ku80 ρ=0.347, Ku70 ρ=0.341, PRKDC ρ=0.245; all q≥0.36). LUSC DNA-PK score ρ=0.124, p=0.57, n=23. STING scores are ρ=0.145 (LUAD) and ρ=−0.112 (LUSC). Do not write a public “CLDN4 protein tracks DNA-PK subunits and STING-pathway proteins” sentence from this matrix. The LUAD DNA-PK sign is consistent across the three subunits, at n=24, and it is not a significant correlation.

## Honest n

133 TMT columns. Three QA replicate pairs (PHLC113, PHLC116, PHLC277) averaged to one model each → **130 models**. PHLC113 CLDN4 is observed in only one of the two replicates; the model value is that single log2.

| | LUAD | LUSC | LUAD+LUSC (sensitivity) | other histology |
|---|---:|---:|---:|---:|
| Models in the protein table | 56 | 61 | 117 | 13 |
| CLDN4 protein quantified | **24** | **23** | **47** | 2 |
| CLDN4 missing | 32 | 38 | 70 | 11 |
| PRKDC, XRCC5, XRCC6, TBK1 quantified | 56 | 61 | 117 | 13 |
| STING1 quantified | 50 | 52 | 102 | 9 |
| CGAS quantified | 51 | 51 | 102 | 12 |
| IRF3 quantified | 46 | 47 | 93 | 9 |
| Tested n, CLDN4 ∩ DNA-PK subunit | **24** | **23** | **47** | not tested |
| Tested n, CLDN4 ∩ STING1 or CGAS or TBK1 | **24** | **23** | **47** | not tested |
| Tested n, CLDN4 ∩ IRF3 | **18** | **17** | **35** | not tested |

Do not write n=133, n=130, or n=56/61 for the Spearman tests. Other histologies (LCNEC, LCC, pleomorphic, combined) contribute 2 CLDN4-quantified models and are left out of both the primary tests and the pool.

## Primary Spearman (human protein)

BH q is within histology, across the 7 pre-specified proteins. Scores are composites of those proteins and are not given a second q.

### LUAD (24 models with CLDN4)

| Partner | n | ρ | p | q |
|---|---:|---:|---:|---:|
| PRKDC (DNA-PKcs) | 24 | 0.245 | 0.25 | 0.49 |
| XRCC5 (Ku80) | 24 | 0.347 | 0.097 | 0.36 |
| XRCC6 (Ku70) | 24 | 0.341 | 0.10 | 0.36 |
| DNA-PK score | 24 | 0.357 | 0.087 | — |
| STING1 | 24 | 0.216 | 0.31 | 0.49 |
| CGAS | 24 | −0.201 | 0.35 | 0.49 |
| TBK1 | 24 | 0.174 | 0.42 | 0.49 |
| IRF3 | 18 | −0.092 | 0.72 | 0.72 |
| STING score | 24 | 0.145 | 0.50 | — |

### LUSC (23 models with CLDN4)

| Partner | n | ρ | p | q |
|---|---:|---:|---:|---:|
| PRKDC | 23 | 0.134 | 0.54 | 0.76 |
| XRCC5 | 23 | 0.045 | 0.84 | 0.91 |
| XRCC6 | 23 | 0.024 | 0.91 | 0.91 |
| DNA-PK score | 23 | 0.124 | 0.57 | — |
| STING1 | 23 | −0.190 | 0.39 | 0.70 |
| CGAS | 23 | −0.229 | 0.29 | 0.70 |
| TBK1 | 23 | 0.225 | 0.30 | 0.70 |
| IRF3 | 17 | −0.218 | 0.40 | 0.70 |
| STING score | 23 | −0.112 | 0.61 | — |

IRF3 is the thin pair (missing in part of the CLDN4-quantified set). A null at n=17–18 is underpowered relative to n=24/23. It is still the tested n.

## Sensitivities (not a second primary)

**Pooled LUAD+LUSC, n=47** (IRF3 n=35). DNA-PK ρ=0.141–0.167 (score 0.168, p=0.26). STING1 ρ=0.094, CGAS ρ=−0.156, TBK1 ρ=0.196 (p=0.19), IRF3 ρ=−0.157. STING score ρ=0.115, p=0.44. All 7 protein q≥0.43. The pool is reported so a combined n is on the record; it is not the primary test, because LUAD and LUSC DNA-PK coefficients are not the same.

**Partial Spearman** on the median log2 of all human proteins (global abundance) leaves the same pattern: LUAD Ku80/Ku70 stay near ρ=0.33–0.34 (p≈0.10–0.12); LUSC DNA-PK stays ~0; STING stays mixed. The LUAD DNA-PK sign is not a global-intensity artifact, and adjusting for intensity does not make it significant.

**CLDN4 quantified vs missing** (Mann–Whitney; Glass rank-biserial > 0 means higher protein when CLDN4 was quantified). This is not a correlation among quantified tumors, and NA is still not imputed. DNA-PK subunits do not differ by CLDN4 detection (LUAD |r|≤0.13, p≥0.43; LUSC r≈0.16–0.17, p≥0.26). Human STING1 is nominally lower when CLDN4 is quantified (LUAD r=−0.30, p=0.069; LUSC r=−0.32, p=0.051) and does not survive BH. The only detection contrast with p<0.05 is **LUSC IFI16** (r=−0.41, p=0.007, q=0.065 within the 9 LUSC detection tests). Among LUSC tumors where both proteins were quantified, the IFI16 Spearman is ρ=−0.157, p=0.47, n=23. LUAD STAT1 is higher in the CLDN4-quantified group (r=0.37, p=0.020, q=0.18) while the quantified-only Spearman is ρ=−0.047, p=0.83, n=24. Neither detection contrast is a DNA-PK result or a STING-core correlation.

**Extended IFN / NF-κB proteins** (IFI16, STAT1, STAT2, IRF9, MAVS, DDX58, IFIH1, IRF7, NFKB1, RELA, JAK1, TYK2), Spearman, BH within histology across these 12: no q<0.45. IRF7 is sparse (LUAD n=9, LUSC n=6). Human TREX1 has no human row (the table’s TREX1 is mouse) and was not tested.

**Mouse stroma rows** exist for PRKDC, XRCC5, XRCC6, STING1, and STAT1. They were not substituted for the human protein. The only stroma coefficient with p<0.05 is LUSC mouse STING1 vs human CLDN4 (ρ=−0.438, p=0.037, n=23). That row is murine stroma, one histology, and outside the human FDR family.

## Methods (this slice)

- Matrix: author-normalized log2 TMT intensities after the paper’s human vs mouse cellularity conversion. No second normalization.
- Species: `H/M == H`. Ambiguous `H/M` rows (IKBKE, ATR) are excluded. Mouse rows are the stroma contrast only.
- Predictor: human CLDN4. No TACSTD2 gate.
- DNA-PK subunits: PRKDC, XRCC5, XRCC6. Score = mean of within-cohort z-scores, all three required.
- STING core: STING1, CGAS, TBK1, IRF3. Score = mean of within-cohort z-scores of members present, at least 2 of 4.
- Test: two-sided Spearman. BH within cohort and analysis, proteins only.
- GSE166999 is the matched RNA matrix and is not used as a protein stand-in.
- No ICI labels, no survival test, no proteotype fishing.

## What this does not claim

- It does not report a significant public correlation between CLDN4 protein and DNA-PK subunits.
- It does not report a positive correlation between CLDN4 protein and STING1, CGAS, TBK1, or IRF3.
- It does not impute the 32 LUAD / 38 LUSC models with no CLDN4 protein as zero, and it does not call that missingness a correlation.
- It does not use mouse stroma protein, RNA, or CPTAC primary tumors as a stand-in.
- It does not pool LUAD and LUSC in the primary test.
- It does not treat n=133 as the tested n.
- LUAD n=24 and LUSC n=23 are modest. A non-significant LUAD DNA-PK ρ≈0.3 is not a precise zero, and it is not a license to borrow a private PDX result.

## Outputs

- `results/n_table.tsv` — models and pairwise n
- `results/presence.tsv` — human row, UniProt entry, observed vs missing models
- `results/spearman.tsv` — primary, partial, extended, and pooled Spearman
- `results/detection.tsv` — CLDN4 quantified vs missing
- `results/stroma_contrast.tsv` — mouse protein vs human CLDN4
- `results/replicate_audit.tsv` — the three averaged pairs
- `results/sample_scores.tsv` — per-model values
- `results/summary.json`
- `results/fig_cldn4_vs_dnapk_sting.png`

```bash
python3 methods/pdx_cldn4_dnapk_sting/download.py --outdir data/pdx_cldn4_dnapk_sting
python3 methods/pdx_cldn4_dnapk_sting/analyze.py --data data/pdx_cldn4_dnapk_sting --outdir methods/pdx_cldn4_dnapk_sting
```
