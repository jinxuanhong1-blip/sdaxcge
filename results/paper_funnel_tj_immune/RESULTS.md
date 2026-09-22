# RESULTS — PAPER FUNNEL TJ / CLDN4 vs T/NK / CD8 / IFN

Public-only PPT corroboration. All numbers copied from cited PRs. No re-fit. No fabrication.

## Slide verdict

| Limb | Supports TJ/CLDN4 → fewer T/NK / IFN? | Lead number |
|---|---|---|
| **Concordant-4** patient CLDN4 %pos vs T/NK | **YES — STRONG** | **ρ=−0.531**; I²=0%; N=65 |
| Largest \|ρ\| (I²=0) | YES (sensitivity) | **ρ=−0.533** |
| Largest \|ρ\| (sign-only) | heterogeneous — do not promote | ρ=−0.552; I²=67% |
| **GSE131907** patient | **YES** | ρ=−0.478 (n=21) |
| Pathway fgsea NES | **YES — IFN/MHC DOWN** | IFNγ **−3.85**; MHC **−2.75** |
| GSVA Spearman | YES (IFN) | IFNγ **−0.40**; IFNα **−0.32** |
| CosMx CLDN4 spatial | **YES** | ratio **0.36 / 0.52** |
| TCGA structural TJ vs CD8 | YES (bulk orthogonal) | ρ=−0.29 |
| Public mouse GEMM forest | **NO** | ρ=**+0.48** vs T/NK |
| Mouse max-\|ρ\| T frac | exploratory only | ρ=−1 on n=4; p=0.083 |

## Files

- `PPT_EVIDENCE_TABLE.md` — paste into slide
- `evidence_table.tsv` — full rows + caveats
- `provenance.json` — PR #503/#539/#712/#541/#691/#698/#90/#677/#700/#200/#715

## Honest skips

- Public mouse does not corroborate the human inverse at the pre-specified forest.
- KEGG TJ NES crosses zero — not “TJ pathway up” on this ranking.
- Junction gene modules (#715) are weaker than CLDN4 %pos alone.
- CD8/NK effector/exhaustion state vs CLDN4 is null (#649); fraction exclusion stays.
