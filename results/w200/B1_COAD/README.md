# B1 analog — TCGA-COAD TACSTD2 × CLDN4 surface rank

Honest ranking of every in-silico surfaceome gene by co-expression with
`TACSTD2` (TROP2) in TCGA colon adenocarcinoma, then reporting where
`CLDN4` actually lands.

**Answer:** CLDN4 is **not** a TACSTD2 surface-gene partner in COAD.
Spearman ρ = 0.007 (FDR q = 0.93), rank **#1876 of 2672**. The same null
reproduces on the independent Xena HiSeqV2 matrix (ρ = −0.058, rank #2136
of 2510).

This is the colon analog of `results/w200/B1_BRCA` (where CLDN4 is #4,
ρ ≈ 0.35). The BRCA result does not transfer.

## Reproduce

```bash
pip install -r requirements.txt
python3 scripts/download_data.py    # GDC STAR-Counts + surfaceome + Xena
python3 scripts/coexpression.py     # writes this directory
```

Raw matrices stay in git-ignored `data/`.
