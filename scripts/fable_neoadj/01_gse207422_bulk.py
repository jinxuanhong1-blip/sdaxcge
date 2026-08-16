"""GSE207422 -- primary analysis.

Bulk RNA-seq (log2 TPM) of 24 PRE-TREATMENT biopsies from resectable NSCLC
patients later treated with neoadjuvant anti-PD-1 + platinum chemotherapy
(Zhang et al., PMID 36869384). Baseline TACSTD2 / CLDN4 expression is tested
against pathologic response (MPR/pCR vs non-MPR) and against residual-tumor %.

This is a baseline-biomarker design: does tumor-intrinsic expression of the
ADC targets TROP2 (TACSTD2) and Claudin-4 (CLDN4) before therapy associate
with major pathologic response to neoadjuvant chemo-immunotherapy?
"""
from __future__ import annotations

import gzip
import json
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import GENES, DATA, FIG, TAB, mwu_stats, spearman

EXPR = os.path.join(DATA, "GSE207422_bulk_log2TPM.txt.gz")
META = os.path.join(DATA, "GSE207422_bulk_metadata.xlsx")


def load_expr(genes):
    with gzip.open(EXPR, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        samples = [s.strip('"') for s in header[1:]]
        want = set(genes)
        rows = {}
        for line in fh:
            g = line.split("\t", 1)[0].strip('"')
            if g in want:
                vals = np.array(line.rstrip("\n").split("\t")[1:], float)
                rows[g] = vals
    return pd.DataFrame(rows, index=samples)


def load_meta():
    df = pd.read_excel(META, sheet_name="sheet1")
    df = df[df["Sample"].notna()].copy()
    df["Sample"] = df["Sample"].astype(str).str.strip()
    resp = df["Pathologic Response"].astype(str).str.strip()
    # Binary MPR: MPR or MPR (pCR) => responder; NMPR => non-responder
    df["MPR"] = resp.str.upper().str.startswith("MPR").astype(int)
    df["pCR"] = resp.str.contains("pCR", case=False).astype(int)
    df["residual"] = pd.to_numeric(df["Residual Tumor"], errors="coerce")
    df["response_label"] = resp
    return df.set_index("Sample")


def main():
    expr = load_expr(GENES + ["EPCAM"])  # EPCAM for tumor-content context
    meta = load_meta()
    common = [s for s in expr.index if s in meta.index]
    expr = expr.loc[common]
    meta = meta.loc[common]

    df = meta[["Patient", "response_label", "MPR", "pCR", "residual",
               "Pathology", "PD1 Antibody", "Chemotherapy"]].copy()
    for g in GENES + ["EPCAM"]:
        df[g] = expr[g].values
    df.index.name = "Sample"
    df.to_csv(os.path.join(TAB, "GSE207422_bulk_per_sample.csv"))

    results = []
    for g in GENES:
        pos = df.loc[df["MPR"] == 1, g].values
        neg = df.loc[df["MPR"] == 0, g].values
        st = mwu_stats(pos, neg)
        # Spearman vs residual tumor fraction (higher residual = worse response)
        sp = spearman(df[g].values, df["residual"].values)
        row = {"dataset": "GSE207422_bulk", "gene": g,
               "comparison": "MPR_vs_nonMPR", "unit": "log2TPM_baseline"}
        row.update(st)
        row.update({"spearman_rho_vs_residual": sp["rho"],
                    "spearman_p_vs_residual": sp["p"], "spearman_n": sp["n"]})
        results.append(row)

    res = pd.DataFrame(results)
    res.to_csv(os.path.join(TAB, "GSE207422_bulk_stats.csv"), index=False)
    print(res.to_string())

    # ---- Figures ----
    _boxplots(df)
    _scatter_residual(df)

    summary = {
        "dataset": "GSE207422_bulk",
        "n_samples": int(len(df)),
        "n_MPR": int((df["MPR"] == 1).sum()),
        "n_nonMPR": int((df["MPR"] == 0).sum()),
        "n_pCR": int((df["pCR"] == 1).sum()),
        "sample_type": "pre-treatment biopsy (baseline)",
        "regimen": "neoadjuvant anti-PD-1 + platinum chemotherapy",
    }
    with open(os.path.join(TAB, "GSE207422_bulk_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


def _boxplots(df):
    fig, axes = plt.subplots(1, len(GENES), figsize=(4.2 * len(GENES), 4.2))
    for ax, g in zip(axes, GENES):
        groups = [df.loc[df["MPR"] == 0, g].values, df.loc[df["MPR"] == 1, g].values]
        bp = ax.boxplot(groups, widths=0.6, patch_artist=True, showfliers=False)
        for patch, c in zip(bp["boxes"], ["#8aa1b1", "#d1495b"]):
            patch.set_facecolor(c)
            patch.set_alpha(0.55)
        rng = np.random.default_rng(0)
        for i, vals in enumerate(groups, 1):
            x = rng.normal(i, 0.06, len(vals))
            ax.scatter(x, vals, s=22, color="#22303c", zorder=3, alpha=0.8)
        st = mwu_stats(groups[1], groups[0])
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"non-MPR\n(n={len(groups[0])})", f"MPR/pCR\n(n={len(groups[1])})"])
        ax.set_ylabel("log2 TPM (baseline)")
        ax.set_title(f"{g}\nMWU p={st['p']:.3f}, AUC={st['auc']:.2f}")
    fig.suptitle("GSE207422 baseline bulk RNA-seq: ADC targets vs MPR", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "GSE207422_bulk_boxplots.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def _scatter_residual(df):
    fig, axes = plt.subplots(1, len(GENES), figsize=(4.2 * len(GENES), 4.2))
    for ax, g in zip(axes, GENES):
        ax.scatter(df["residual"].values, df[g].values, s=30, color="#2a6f97")
        sp = spearman(df[g].values, df["residual"].values)
        ax.set_xlabel("residual tumor fraction")
        ax.set_ylabel(f"{g} log2 TPM")
        ax.set_title(f"{g}: Spearman rho={sp['rho']:.2f}, p={sp['p']:.3f}")
    fig.suptitle("GSE207422 baseline expression vs residual tumor", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "GSE207422_bulk_scatter_residual.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
