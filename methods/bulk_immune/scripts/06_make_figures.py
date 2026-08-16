#!/usr/bin/env python3
"""Optional figures for the playbook. Fail soft if matplotlib is missing."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _plot_cohort(indir: Path, outdir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    corr = pd.read_csv(indir / "correlation_spearman.tsv", sep="\t")
    highlight = corr[corr["family"].isin(["exclusion", "inflamed"])].copy()
    if highlight.empty:
        return
    for target, sub in highlight.groupby("target"):
        sub = sub.dropna(subset=["spearman_r"]).sort_values("spearman_r")
        if sub.empty:
            continue
        fig, ax = plt.subplots(figsize=(8, max(3, 0.28 * len(sub))))
        colors = ["#b2182b" if f == "exclusion" else "#2166ac" for f in sub["family"]]
        ax.barh(sub["score"], sub["spearman_r"], color=colors)
        ax.axvline(0, color="black", lw=0.6)
        ax.set_xlabel(f"Spearman r vs {target}")
        ax.set_title(f"{indir.name}: {target} vs exclusion (red) / inflamed (blue)")
        fig.tight_layout()
        fig.savefig(outdir / f"{indir.name}_{target}_highlight.png", dpi=140)
        plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    args = p.parse_args()
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        print("[figures] matplotlib not installed; skipping")
        return
    results = args.root / "results"
    for name in ("GSE126044", "GSE135222", "TCGA_NSCLC"):
        indir = results / name
        if (indir / "correlation_spearman.tsv").exists():
            _plot_cohort(indir, indir)
            print("[figures]", name)


if __name__ == "__main__":
    main()
