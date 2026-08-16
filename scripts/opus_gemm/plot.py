#!/usr/bin/env python3
"""Figures for the Tacstd2 / Cldn4 NSCLC GEMM ICI slice."""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = "results/opus_gemm/figures"
os.makedirs(OUT, exist_ok=True)

contr = pd.read_csv("results/opus_gemm/focus_contrasts.tsv", sep="\t")
vals = pd.read_csv("results/opus_gemm/sample_values.tsv", sep="\t")
corr = pd.read_csv("results/opus_gemm/signature_correlations.tsv", sep="\t")


def _save(fig, name: str) -> None:
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)


# ---- 1. forest of ICI + genotype contrasts for both genes ----
keep_fam = {"ici", "genotype", "resistance"}
sub = contr[contr["family"].isin(keep_fam)].copy()
sub = sub[sub["n_test"] >= 2]
sub["label"] = sub["dataset"].str.replace(r"_aPD1|_total|_nodules|_adenoma|_AMPK|_ICBrelapse", "", regex=True)
sub["label"] = sub["accession"] + " " + sub["contrast"] + "  " + sub["symbol"]
sub = sub.sort_values(["family", "symbol", "delta_log2"])

fig, ax = plt.subplots(figsize=(9.5, max(3.5, 0.32 * len(sub))))
y = np.arange(len(sub))
colors = {"ici": "#1f4e79", "genotype": "#7b2d26", "resistance": "#6b4c9a"}
for fam, chunk in sub.groupby("family", sort=False):
    yy = y[sub["family"].to_numpy() == fam]
    ax.errorbar(
        chunk["delta_log2"], yy, xerr=None, fmt="o", color=colors[fam],
        label=fam, markersize=5,
    )
    for yi, (_, r) in zip(yy, chunk.iterrows()):
        q = r.get("fdr_bh_within_family")
        mark = ""
        if pd.notna(r["p_value"]):
            if pd.notna(q) and q < 0.05:
                mark = " *"
            elif r["p_value"] < 0.05:
                mark = " ·"
        ax.text(r["delta_log2"] + 0.04, yi, f"n={int(r['n_test'])}/{int(r['n_ref'])}{mark}",
                va="center", fontsize=7, color="#333")
ax.axvline(0, color="#888", lw=0.8)
ax.set_yticks(y)
ax.set_yticklabels(sub["label"].tolist(), fontsize=7)
ax.set_xlabel("Δ log2 (test − reference)")
ax.set_title("Tacstd2 / Cldn4 contrasts  (* FDR<0.05 within family; · raw p<0.05)")
ax.legend(loc="lower right", fontsize=8)
ax.grid(axis="x", ls=":", alpha=0.4)
_save(fig, "forest_contrasts.png")

# ---- 2. KL aPD1 strip plots (primary ICI series) ----
kl = vals[vals["dataset"] == "GSE182228_KL_aPD1"].copy()
order = ["vehicle", "aPD1", "palbociclib", "aPD1_palbociclib"]
fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6), sharex=True)
rng = np.random.default_rng(1)
for ax, gene in zip(axes, ["Tacstd2", "Cldn4"]):
    for i, g in enumerate(order):
        yv = kl.loc[kl["group"] == g, gene].to_numpy(dtype=float)
        x = np.full(len(yv), i) + rng.uniform(-0.08, 0.08, len(yv))
        ax.scatter(x, yv, s=28, c="#1f4e79", zorder=3)
        ax.hlines(np.mean(yv), i - 0.22, i + 0.22, color="#111", lw=1.4)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(["vehicle", "aPD-1", "palbo", "aPD-1+palbo"], rotation=20, ha="right")
    ax.set_ylabel(f"{gene}  log2(FPKM+1)")
    ax.set_title(gene)
    ax.grid(axis="y", ls=":", alpha=0.4)
fig.suptitle("GSE182228  LKB1-deficient LUAD tumours", fontsize=11)
fig.tight_layout()
_save(fig, "gse182228_kl_apd1.png")

# ---- 3. KL vs KP genotype (two independent GEMM nodule/adenoma series) ----
fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8))
panels = [
    ("GSE137396_KL_vs_KP_nodules", ["KL_nodule", "KP_nodule"], ["KL", "KP"]),
    ("GSE274351_K_KP_KL_adenoma", ["K_adenoma", "KP_adenoma", "KL_adenoma", "normal_lung"],
     ["K", "KP", "KL", "NL"]),
]
for ax, (key, groups, labels) in zip(axes, panels):
    d = vals[vals["dataset"] == key]
    rng = np.random.default_rng(2)
    for gene, col, off in [("Tacstd2", "#1f4e79", -0.15), ("Cldn4", "#b85c38", 0.15)]:
        for i, g in enumerate(groups):
            yv = d.loc[d["group"] == g, gene].to_numpy(dtype=float)
            x = np.full(len(yv), i) + off + rng.uniform(-0.05, 0.05, len(yv))
            ax.scatter(x, yv, s=22, c=col, label=gene if i == 0 else None, zorder=3)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_ylabel("log2 expression")
    ax.set_title(key.split("_", 1)[0])
    ax.grid(axis="y", ls=":", alpha=0.4)
    ax.legend(fontsize=7)
fig.suptitle("Tacstd2 / Cldn4 by KRAS-GEMM genotype", fontsize=11)
fig.tight_layout()
_save(fig, "genotype_kl_kp.png")

# ---- 4. correlation heatmap (in-vivo) ----
if not corr.empty:
    corr["p"] = pd.to_numeric(corr["p_value"], errors="coerce")
    corr["rho"] = pd.to_numeric(corr["spearman_rho"], errors="coerce")
    # one cell per dataset x gene x signature
    pivot_idx = corr["dataset"].str.replace(r"_.*", "", regex=True) + " " + corr["symbol"]
    mat = corr.assign(idx=pivot_idx).pivot_table(index="idx", columns="signature", values="rho")
    fig, ax = plt.subplots(figsize=(9.2, max(3.2, 0.38 * len(mat))))
    im = ax.imshow(mat.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.03, label="Spearman ρ")
    ax.set_title("Tacstd2 / Cldn4 vs immune signatures (in-vivo / sorted)")
    fig.tight_layout()
    _save(fig, "signature_correlations.png")

print("done")
