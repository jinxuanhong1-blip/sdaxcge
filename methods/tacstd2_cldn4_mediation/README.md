# TACSTD2–immune partial on CLDN4, and mediation by the locked 221-gene score

OncoSG LUAD plus the PR 590 GEO LUAD bulks (GSE273377 discovery, GSE273377 validation, GSE282774, GSE233774 tumors). The 221 genes are the locked PR 590 signature (`data/signature_genes.tsv`, SHA-256 `9d398e71f8a2569d3c5c7bb940577a24eaca1f19db0a1d589c84e13a0bac6a3a`). CLDN4 and TACSTD2 are not in that list.

The script recomputes partial Spearman correlations and a rank-regression mediation model. It does not refit the signature. Write-up: [FINDING.md](FINDING.md).

```bash
pip install -r methods/tacstd2_cldn4_mediation/requirements.txt
python3 methods/tacstd2_cldn4_mediation/analyze.py
```
