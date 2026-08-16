# Claims B3 + B4 reproduction

Independent, honest reproduction of two tight-junction (TJ) signature claims
from user PPT 2026-08-17 (`claims/B3.md`, `claims/B4.md`):

- **B3** — TCGA-LUAD: TJ-high tumors have low CD8 / low GEP (`p < 1e-6`).
- **B4** — GSE126044 (NSCLC anti-PD-1): non-responders have higher TJ (`p = 0.019`).

Full write-up, tables, and figures: [`results/claim_B3B4/`](results/claim_B3B4/README.md).

## Verdict

- **B3 is supported.** Official 15-gene TJ vs CD8 Spearman ρ = −0.29, p = 1.6e-11
  (GEP p = 9.0e-11). Survives Aran ESTIMATE purity adjustment (partial ρ = −0.26,
  p = 1.7e-9). The broad KEGG tight-junction pathway does **not** support the
  claim vs CD8 (ρ ≈ +0.07, p ≈ 0.10).
- **B4 p = 0.019 is recovered only for the 7-gene module named on the B3 claim
  page** (`CLDN1, CLDN4, CLDN7, F11R, TJP1, TJP2, OCLN`) scored as z-mean +
  Mann-Whitney U. The 15-gene set, KEGG/GOCC, ssGSEA of the same 7 genes, and
  CLDN4 alone do not reach 0.019. Small n (5 R / 11 NR); definition-dependent.

## Quick start

```bash
pip install -r requirements.txt
bash scripts/download_data.sh
python3 scripts/analyze_B3B4.py
python3 scripts/sensitivity_and_figures.py
python3 scripts/analyze_B3B4_continue.py
```
