# B6_radius — tumor-only CLDN4 vs immune neighborhoods at hex rings 1/2/3

**Honest verdict: DOES_NOT_SUPPORT_CLAIM.**

Reworked Visium (E-MTAB-13530, 20 tumour sections): CLDN4 vs broad-immune hex-ring neighborhoods, tumor spots only, immune-rich spots excluded from the CLDN4 score.

| Ring | median partial ρ | Wilcoxon p |
|---|---:|---:|
| 1 | −0.038 | 0.22 |
| 2 | +0.016 | 0.78 |
| 3 | +0.013 | 0.60 |

|median ρ| ≤ 0.038, still ≤ the PR67 mismatch band of 0.06. The user spatial-exclusion claim is not supported. Details: `WRITEUP.md`, `summary.json`.

## Rework rules (pre-specified, not tuned)

1. Tumor-spot only (E-MTAB-13530 `P*_T*`; epithelial-high).
2. Visium hex rings 1 / 2 / 3 (self excluded).
3. Immune-rich spots excluded from the CLDN4 score (immune ≥ section Q3).
4. GeoMx GSE271689 = tumor / CK compartment only (same-AOI correlation; no rings).

## Run

```bash
pip install -r results/rework/B6_radius/requirements.txt
python3 results/rework/B6_radius/download.py
python3 results/rework/B6_radius/analyze.py
```

Processed inputs cache under `/tmp/b6_radius_data` (not committed).
