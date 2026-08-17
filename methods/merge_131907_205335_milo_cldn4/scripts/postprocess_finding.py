#!/usr/bin/env python3
"""Refresh FINDING + sensitivity tables from an already-finished run."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_merge import extra_summary_figures, fdr_counts, min_present_sensitivity, write_finding  # noqa: E402


def main() -> None:
    root = Path("methods/merge_131907_205335_milo_cldn4")
    tables = root / "tables"
    figures = root / "figures"
    summary_path = tables / "summary.json"
    overall = json.loads(summary_path.read_text())
    hit_rows = []
    for name, graph in overall["graphs"].items():
        da_path = tables / name / "da_malignant_cldn4.tsv"
        nh_path = tables / name / "nhoods.tsv"
        if not da_path.exists():
            continue
        da = pd.read_csv(da_path, sep="\t")
        graph["da_malignant_cldn4"] = fdr_counts(da)
        graph["sensitivity_min_present"] = min_present_sensitivity(da)
        t = da[da["testable"]]
        sig = t[t["SpatialFDR"] < 0.1]
        nh = pd.read_csv(nh_path, sep="\t") if nh_path.exists() else pd.DataFrame()
        for _, r in sig.iterrows():
            rec = {
                "graph": name,
                "nhood": int(r["nhood"]),
                "n_samples_present": int(r["n_samples_present"]),
                "spearman_rho": float(r["spearman_rho"]),
                "p": float(r["p"]),
                "SpatialFDR": float(r["SpatialFDR"]),
                "perfect_rank": abs(float(r["spearman_rho"])) >= 0.999,
            }
            if len(nh):
                row = nh.loc[nh["nhood"] == r["nhood"]]
                if len(row):
                    rec["dominant_lineage"] = row.iloc[0].get("dominant_lineage")
                    rec["n_malig"] = int(row.iloc[0].get("n_malig", -1))
                    rec["n_tnk"] = int(row.iloc[0].get("n_tnk", -1))
            hit_rows.append(rec)
        fig, ax = plt.subplots(figsize=(6.2, 4.4))
        ax.scatter(t["n_samples_present"], t["spearman_rho"], s=10, c="#4c72b0", alpha=0.5)
        if len(sig):
            ax.scatter(sig["n_samples_present"], sig["spearman_rho"], s=28, c="#c44e52", label="SpatialFDR<0.1")
        ax.axhline(1, c="0.7", lw=0.6)
        ax.axhline(-1, c="0.7", lw=0.6)
        ax.set_xlabel("n units present")
        ax.set_ylabel("Spearman ρ (nhood prop vs malignant CLDN4)")
        ax.set_title(f"{name}: |ρ|=1 hits sit on the n=5 floor")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(figures / f"fig_rho_vs_n_present_{name}.png", dpi=160)
        plt.close(fig)
    if hit_rows:
        pd.DataFrame(hit_rows).to_csv(tables / "spatialfdr_hits.tsv", sep="\t", index=False)
    extra_summary_figures(overall["graphs"], figures, tables)
    tables.joinpath("summary.json").write_text(json.dumps(overall, indent=2, default=str))
    write_finding(overall, root / "FINDING.md")
    print("updated FINDING + sensitivity", flush=True)


if __name__ == "__main__":
    main()
