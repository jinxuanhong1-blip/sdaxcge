#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path("results/stereo_seq_luad")


def _order(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values("patient", key=lambda s: s.str.replace("LUAD_P", "").astype(int))


def fig_rho(df: pd.DataFrame, title: str, path: Path) -> None:
    df = _order(df)
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    x = np.arange(len(df))
    ax.axhline(0, color="#666", lw=0.8)
    ax.plot(x, df["rho_cldn4_cd8a"], "o-", color="#1f4e79", label="all units CLDN4 vs CD8A")
    ax.plot(x, df["rho_krt8_residual_cd8a"], "s--", color="#c45c26", label="KRT8 residual vs CD8A")
    ax.plot(x, df["rho_cldn4_cd8a_epi"], "d:", color="#2a7f62", label="epithelial-only CLDN4 vs CD8A")
    ax.set_xticks(x)
    ax.set_xticklabels(df["patient"], rotation=45, ha="right")
    ax.set_ylabel("Spearman ρ")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    ax.set_ylim(-0.55, 0.25)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_distance(df: pd.DataFrame, title: str, path: Path) -> None:
    df = _order(df)
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    x = np.arange(len(df))
    w = 0.35
    ax.bar(x - w / 2, df["median_nn_um_high"], w, color="#1f4e79", label="CLDN4-high epithelial")
    ax.bar(x + w / 2, df["median_nn_um_low"], w, color="#9aa5b1", label="CLDN4-low epithelial")
    ax.set_xticks(x)
    ax.set_xticklabels(df["patient"], rotation=45, ha="right")
    ax.set_ylabel("Median nearest CD8A+ distance (µm)")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main():
    both = pd.read_csv(OUT / "tables" / "per_sample_metrics.csv")
    cell = both[both["unit"] == "cell_bin"]
    bin50 = both[both["unit"] == "bin50"]
    fig_rho(cell, "GSE328481 cell-bin Spearman (CLDN4-only)", OUT / "figures" / "cellbin_spearman")
    fig_rho(bin50, "GSE328481 bin50 Spearman (CLDN4-only)", OUT / "figures" / "bin50_spearman")
    fig_distance(cell, "GSE328481 cell-bin nearest CD8A+ from epithelial bins", OUT / "figures" / "cellbin_nn_cd8")
    fig_distance(bin50, "GSE328481 bin50 nearest CD8A+ from epithelial bins", OUT / "figures" / "bin50_nn_cd8")


if __name__ == "__main__":
    main()
