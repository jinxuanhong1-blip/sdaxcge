# CPTAC LSCC protein TACSTD2/CLDN4 vs xCell CD8/immune (purity residual)

Focused recompute of the PR23 xCell slice. Outputs only under
`results/w200/CPTAC_LSCC_xcell/`.

```bash
pip install -r scripts/w200/CPTAC_LSCC_xcell/requirements.txt
python3 scripts/w200/CPTAC_LSCC_xcell/download.py
python3 scripts/w200/CPTAC_LSCC_xcell/analyze.py
```

- Proteins: TACSTD2 `ENSG00000184292`, CLDN4 `ENSG00000189143`
- xCell: `xCell_T_cell_CD8+`, `xCell_immune_score`
- Primary purity residual: `WES_purity` (DNA). Sensitivity: `WGS_purity`, ESTIMATE cosine purity.
- Cohort is treatment-naive (no ICI labels). See `results/w200/CPTAC_LSCC_xcell/WRITEUP.md`.
