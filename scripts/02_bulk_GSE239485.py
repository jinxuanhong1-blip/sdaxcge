#!/usr/bin/env python3
"""GSE239485 (mouse LLC tumors, bulk RNA-seq, n=24): Tacstd2-high vs -low GSEA.

Input: GSE239485_Processed_data.xlsx, sheet "DataNorm" (log2-scale normalized
expression; values approx. log2 CPM, range ~[-5.3, 14.4]).

Primary analysis : median split of samples on Tacstd2; per-gene Welch t
                   statistic (high vs low) as the ranking metric; preranked
                   GSEA against mouse keratin / tight-junction gene sets.
Sensitivity      : same contrast with treatment arm (C/D/T) as covariate in a
                   per-gene OLS (expr ~ group + arm), because Tacstd2 is
                   confounded with treatment in this study.
"""
from pathlib import Path

import gseapy as gp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

XLSX = Path("/workspace/data/GSE239485/GSE239485_Processed_data.xlsx")
GMT = Path("/workspace/results/w200/A8_mouse/genesets/mouse_keratin_tj.gmt")
OUT = Path("/workspace/results/w200/A8_mouse/GSE239485_bulk")
OUT.mkdir(parents=True, exist_ok=True)
SEED = 42


def load_expression():
    df = pd.read_excel(XLSX, sheet_name="DataNorm")
    samples = [c for c in df.columns if c[:2] in ("C_", "D_", "T_")]
    # collapse duplicate symbols: keep the row with highest mean expression
    df["_mean"] = df[samples].mean(axis=1)
    df = (
        df.sort_values("_mean", ascending=False)
        .drop_duplicates("gene_name")
        .set_index("gene_name")
    )
    expr = df[samples]
    # drop genes never rising above ~1 CPM (log2 <= 0) in any sample
    expr = expr[(expr > 0).any(axis=1)]
    return expr, samples


def welch_t(expr, hi, lo):
    a, b = expr[hi].to_numpy(), expr[lo].to_numpy()
    t, _ = stats.ttest_ind(a, b, axis=1, equal_var=False)
    return pd.Series(t, index=expr.index)


def adjusted_t(expr, samples, hi):
    """t statistic of the Tacstd2-group coefficient in expr ~ group + arm."""
    group = np.array([1.0 if s in hi else 0.0 for s in samples])
    arm_d = np.array([1.0 if s.startswith("D_") else 0.0 for s in samples])
    arm_t = np.array([1.0 if s.startswith("T_") else 0.0 for s in samples])
    X = np.column_stack([np.ones_like(group), group, arm_d, arm_t])
    Y = expr[samples].to_numpy().T  # samples x genes
    beta, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ beta
    dof = X.shape[0] - X.shape[1]
    sigma2 = (resid**2).sum(axis=0) / dof
    xtx_inv = np.linalg.inv(X.T @ X)
    se = np.sqrt(sigma2 * xtx_inv[1, 1])
    with np.errstate(divide="ignore", invalid="ignore"):
        t = beta[1] / se
    return pd.Series(t, index=expr.index).replace([np.inf, -np.inf], np.nan).dropna()


def run_prerank(rnk, label):
    res = gp.prerank(
        rnk=rnk.sort_values(ascending=False).reset_index(),
        gene_sets=str(GMT),
        min_size=5,
        max_size=1000,
        permutation_num=10000,
        seed=SEED,
        threads=4,
        outdir=None,
        no_plot=True,
    )
    df = res.res2d.copy()
    df.insert(0, "ranking", label)
    return res, df


def nes_barplot(df, path, title):
    df = df.sort_values("NES")
    colors = ["#c0392b" if f < 0.05 else "#f1948a" if f < 0.25 else "#95a5a6"
              for f in df["FDR q-val"].astype(float)]
    fig, ax = plt.subplots(figsize=(7, 0.42 * len(df) + 1.4))
    ax.barh(df["Term"], df["NES"].astype(float), color=colors)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("NES (positive = enriched in Tacstd2-high)")
    ax.set_title(title, fontsize=10)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in ("#c0392b", "#f1948a", "#95a5a6")]
    ax.legend(handles, ["FDR<0.05", "FDR<0.25", "n.s."], fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


FOCAL = [
    "Tacstd2", "Cldn4", "Cldn1", "Cldn7", "Ocln", "Tjp1", "Tjp2",
    "Krt5", "Krt8", "Krt17", "Krt18", "Krt19", "Cdh1", "Epcam",
]


def focal_table(expr, samples, anchor):
    rows = []
    a = expr.loc[anchor, samples].astype(float)
    for g in FOCAL:
        if g not in expr.index or g == anchor:
            continue
        b = expr.loc[g, samples].astype(float)
        r, p = stats.spearmanr(a, b)
        rows.append({
            "anchor": anchor, "gene": g,
            "spearman_r": r, "spearman_p": p,
            "gene_median": float(b.median()),
            "gene_max": float(b.max()),
            "n_samples_gt0": int((b > 0).sum()),
        })
    return pd.DataFrame(rows)


def median_split(series):
    order = series.sort_values(ascending=False)
    n = len(order)
    return list(order.index[: n // 2]), list(order.index[n // 2 :])


def main():
    expr, samples = load_expression()
    tac = expr.loc["Tacstd2", samples].astype(float)
    cld = expr.loc["Cldn4", samples].astype(float)
    hi, lo = median_split(tac)
    order = tac.sort_values(ascending=False)

    cld_hi, cld_lo = median_split(cld)
    meta = pd.DataFrame(
        {
            "sample": samples,
            "arm": [{"C": "Control", "D": "PolyIC+antiPD1", "T": "PolyIC+antiPD1+antiC5aR1"}[s[0]] for s in samples],
            "Tacstd2_log2norm": tac[samples].to_numpy(),
            "Cldn4_log2norm": cld[samples].to_numpy(),
            "tac_group": ["Tacstd2_high" if s in hi else "Tacstd2_low" for s in samples],
            "cld_group": ["Cldn4_high" if s in cld_hi else "Cldn4_low" for s in samples],
        }
    )
    meta.to_csv(OUT / "sample_groups.tsv", sep="\t", index=False)
    xtab = pd.crosstab(meta["arm"], meta["tac_group"])
    xtab.to_csv(OUT / "group_by_arm_crosstab.tsv", sep="\t")
    print(xtab, "\n")
    print("Cldn4 vs arm:\n", pd.crosstab(meta["arm"], meta["cld_group"]), "\n")

    focal_table(expr, samples, "Tacstd2").to_csv(OUT / "focal_vs_Tacstd2.tsv", sep="\t", index=False)
    focal_table(expr, samples, "Cldn4").to_csv(OUT / "focal_vs_Cldn4.tsv", sep="\t", index=False)
    r_tc, p_tc = stats.spearmanr(tac, cld)
    print(f"Tacstd2 vs Cldn4 Spearman r={r_tc:.3f} p={p_tc:.4g}")
    print(f"Tacstd2>0: {(tac>0).sum()}/24  Cldn4>0: {(cld>0).sum()}/24  Cldn4 max={cld.max():.3f}")

    fig, ax = plt.subplots(figsize=(6, 3.5))
    for i, arm in enumerate(["Control", "PolyIC+antiPD1", "PolyIC+antiPD1+antiC5aR1"]):
        sub = meta[meta.arm == arm]
        ax.scatter(
            np.full(len(sub), i) + np.random.default_rng(SEED).uniform(-0.12, 0.12, len(sub)),
            sub.Tacstd2_log2norm,
            c=["#c0392b" if g == "Tacstd2_high" else "#2980b9" for g in sub.tac_group],
            s=28,
        )
    ax.axhline(order.iloc[11], ls="--", c="grey", lw=0.8)
    ax.set_xticks(range(3), ["Control", "PolyIC\n+antiPD1", "PolyIC+antiPD1\n+antiC5aR1"])
    ax.set_ylabel("Tacstd2 log2 normalized expr")
    ax.set_title("GSE239485: Tacstd2 by arm (red = 'high' group; dashed = split)", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "tacstd2_by_arm.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    ax.scatter(tac, cld, c=["#c0392b" if s[0] == "C" else "#27ae60" if s[0] == "D" else "#8e44ad" for s in samples], s=28)
    ax.set_xlabel("Tacstd2"); ax.set_ylabel("Cldn4")
    ax.set_title(f"GSE239485 LLC: Tacstd2 vs Cldn4\nSpearman r={r_tc:.2f} p={p_tc:.2g}", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "tacstd2_vs_cldn4.png", dpi=200)
    plt.close(fig)

    all_res = []
    rnk1 = welch_t(expr, hi, lo).dropna().drop(index="Tacstd2", errors="ignore")
    rnk1.sort_values(ascending=False).to_csv(OUT / "rank_median_split_welch_t.rnk", sep="\t", header=False)
    res1, df1 = run_prerank(rnk1, "Tacstd2_median_split_welch_t")
    all_res.append(df1)

    rnk2 = adjusted_t(expr, samples, hi).drop(index="Tacstd2", errors="ignore")
    rnk2.sort_values(ascending=False).to_csv(OUT / "rank_arm_adjusted_t.rnk", sep="\t", header=False)
    _, df2 = run_prerank(rnk2, "Tacstd2_arm_adjusted_t")
    all_res.append(df2)

    rnk3 = welch_t(expr, cld_hi, cld_lo).dropna().drop(index="Cldn4", errors="ignore")
    rnk3.sort_values(ascending=False).to_csv(OUT / "rank_Cldn4_median_split_welch_t.rnk", sep="\t", header=False)
    _, df3 = run_prerank(rnk3, "Cldn4_median_split_welch_t")
    all_res.append(df3)

    gsea = pd.concat(all_res, ignore_index=True)
    gsea.to_csv(OUT / "gsea_results.tsv", sep="\t", index=False)
    print(gsea[["ranking", "Term", "NES", "NOM p-val", "FDR q-val"]].to_string(index=False))

    nes_barplot(df1, OUT / "nes_primary_median_split.png",
                "GSE239485 bulk: Tacstd2-high vs -low (median split, Welch t)")
    nes_barplot(df2, OUT / "nes_sensitivity_arm_adjusted.png",
                "GSE239485 bulk: Tacstd2 group adjusted for treatment arm")
    nes_barplot(df3, OUT / "nes_Cldn4_median_split.png",
                "GSE239485 bulk: Cldn4-high vs -low (median split; Cldn4 is at floor)")

    for term in ("CURATED_KERATIN_FAMILY", "CURATED_TIGHT_JUNCTION_CORE",
                 "GO_KERATIN_FILAMENT", "GO_BICELLULAR_TIGHT_JUNCTION",
                 "BARRIER_CORE_TACSTD2_CLDN4"):
        try:
            ax = res1.plot(terms=term)
            fig = ax.figure if hasattr(ax, "figure") else plt.gcf()
            fig.savefig(OUT / f"enrichment_{term}.png", dpi=200, bbox_inches="tight")
            plt.close("all")
        except Exception as e:  # plotting must not kill the stats
            print(f"plot failed for {term}: {e}")


if __name__ == "__main__":
    main()
