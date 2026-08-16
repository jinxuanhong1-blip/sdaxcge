# Claims B3 + B4 reproduction

Independent, honest reproduction of two tight-junction (TJ) signature claims:

- **B3** — TCGA-LUAD: TJ signature vs CD8 and GEP signatures.
- **B4** — GSE126044 (NSCLC anti-PD-1): non-responders vs responders TJ signature.

## Quick start

```bash
pip install -r requirements.txt
bash scripts/download_data.sh          # fetch TCGA-LUAD, GSE126044, KEGG hsa04530
python3 scripts/analyze_B3B4.py        # core B3 + B4 statistics
python3 scripts/sensitivity_and_figures.py  # B4 robustness sweep + figures
```

Results, figures, and a detailed honest write-up are in
[`results/claim_B3B4/`](results/claim_B3B4/README.md).

## Summary

- **B3 is supported** (Spearman `p ≈ 1e-9`, well below `p < 1e-6`) — but the
  correlation with CD8/GEP is **negative** and depends on the TJ gene set.
- **B4 is not reproduced as stated** — the NR>R direction is robust across 36
  configurations, but the effect is not significant (best one-sided `p ≈ 0.06`);
  the claimed `p = 0.019` could not be recovered.
