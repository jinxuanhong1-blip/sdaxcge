"""Independent sanity checks on the committed B1_STAD ranking.

Recomputes the TACSTD2–CLDN4 Spearman/Pearson coefficients with SciPy
(not the vectorised ranker) and checks that the written tables match.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
EXPR = ROOT / "data" / "TCGA-STAD.HiSeqV2.gz"
OUT = ROOT / "results" / "w200" / "B1_STAD"


def main() -> int:
    if not EXPR.exists():
        print("[skip] expression matrix not present; run scripts/download_data.py")
        return 0

    expr = pd.read_csv(EXPR, sep="\t", index_col=0)
    samples = [
        c
        for c in expr.columns
        if len(c.split("-")) >= 4 and c.split("-")[3][:2] == "01"
    ]
    x = expr.loc["TACSTD2", samples].astype(float)
    y = expr.loc["CLDN4", samples].astype(float)
    rho, p_s = stats.spearmanr(x, y)
    r, p_p = stats.pearsonr(x, y)

    full = pd.read_csv(OUT / "coexpression_TACSTD2_surfaceome.csv")
    cl = full.loc[full["gene"] == "CLDN4"].iloc[0]
    summary = json.loads((OUT / "summary.json").read_text())

    errors: list[str] = []
    if abs(rho - cl["spearman_r"]) > 1e-10:
        errors.append(f"Spearman mismatch: scipy={rho} table={cl['spearman_r']}")
    if abs(r - cl["pearson_r"]) > 1e-10:
        errors.append(f"Pearson mismatch: scipy={r} table={cl['pearson_r']}")
    if int(cl["spearman_rank"]) != 33:
        errors.append(f"unexpected Spearman rank {cl['spearman_rank']}")
    if int(cl["pearson_rank"]) != 26:
        errors.append(f"unexpected Pearson rank {cl['pearson_rank']}")
    if summary["n_samples"] != 415:
        errors.append(f"unexpected n_samples {summary['n_samples']}")
    if summary["focus_result"]["is_top_spearman"] is not False:
        errors.append("focus_result.is_top_spearman should be false")
    if full.iloc[0]["gene"] != "PVRL4":
        errors.append(f"expected #1 PVRL4, got {full.iloc[0]['gene']}")
    if len(full) != 2618:
        errors.append(f"expected 2618 surface genes, got {len(full)}")
    if int(full["spearman_r"].notna().sum()) != 2617:
        errors.append("expected 2617 ranked genes")
    if set(full.loc[full["spearman_r"].isna(), "gene"]) != {"OR2T29"}:
        errors.append("expected only OR2T29 unrankable")
    top200 = pd.read_csv(OUT / "top200.csv")
    if len(top200) != 200:
        errors.append(f"top200 has {len(top200)} rows")
    if "CLDN4" not in set(top200["gene"]):
        errors.append("CLDN4 missing from top200")

    if errors:
        print("FAIL")
        for e in errors:
            print(" -", e)
        return 1

    print("OK")
    print(f"  n_primary={len(samples)}")
    print(f"  scipy spearman rho={rho:.6f} p={p_s:.3e}")
    print(f"  scipy pearson  r={r:.6f} p={p_p:.3e}")
    print(f"  table rank Spearman #{int(cl['spearman_rank'])} Pearson #{int(cl['pearson_rank'])}")
    print(f"  #1={full.iloc[0]['gene']} rho={full.iloc[0]['spearman_r']:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
