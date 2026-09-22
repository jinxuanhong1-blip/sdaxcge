# TACSTD2–CLDN4 same-cell coexpression BRIDGE (not mediation)

Public evidence that TACSTD2 and CLDN4 **co-occur in the same malignant/tumor cells**, mark an **overlapping cold niche**, and should be framed as a **bridge** — not as TACSTD2→CLDN4→immune mediation (null on CosMx + concordant-4, PR #718).

## PPT numbers (computed here)

| Source | Number |
|---|---|
| Concordant-4 within-cell Spearman | median **ρ = 0.746** (63 units ≥30 malignant; 100% ρ>0) |
| Concordant-4 detection lift | median **1.31×** vs independence; median OR **13.65** |
| CosMx within-cell Spearman | median **ρ = 0.180** (8/8 >0; matches PR #726 inventory) |
| CosMx detection lift | median **1.17×**; **8/8** sections >1; median OR **2.87** |
| DepMap Gygi lung protein | **ρ = 0.693** (n=45) → rounds to **0.69** |
| CosMx cold niche @10 µm | both-high immune frac **0.033** vs both-low **0.091**; HH<LL **8/8** (PR #726) |

## Reproduce

```bash
bash scripts/download.sh /tmp/geo_c4 /workspace/data/cosmx_nsclc
export GEO_C4=/tmp/geo_c4
export COSMX_H5AD=/workspace/data/cosmx_nsclc/cosmx_human_nsclc_clustered.h5ad
python3 -m pip install -r requirements.txt
# R + Matrix required for GSE205335 focus export
python3 scripts/analyze_bridge.py
```

Write-up: `FINDING.md`. PPT paste table: `results/PPT_EVIDENCE_TABLE.md`.

## Not recomputed / not claimed

- Locked CosMx CLDN4 cytotoxic ratios **0.36 / 0.52**
- Locked concordant-4 CLDN4 %pos vs T/NK **ρ = −0.531**
- Mediation / “CLDN4 explains TACSTD2 exclusion”
- Private 8KL matrices
