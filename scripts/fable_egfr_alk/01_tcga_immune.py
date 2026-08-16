"""Immune axis in TCGA lung cohorts.

Question: within EGFR-mutant / ALK-driven vs wild-type lung tumours, how do the
ADC/TJ targets TACSTD2 (TROP2) and CLDN4 (Claudin-4) relate to the anti-tumour
immune microenvironment?

Outputs (results/fable_egfr_alk/):
  tables/tcga_<cohort>_target_immune_correlations.csv
  tables/tcga_<cohort>_target_by_driver_group.csv
  tables/tcga_<cohort>_signature_scores.csv (per-sample scores + targets)
  figures/tcga_<cohort>_target_immune_heatmap.png
  figures/tcga_<cohort>_target_by_driver.png
"""
from __future__ import annotations

import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as C

TARGET_NAMES = list(C.TARGETS.keys())              # TACSTD2, CLDN4
IMMUNE_MARKERS = list(C.IMMUNE_GENES.keys())
DRIVER_ORDER = ["EGFR_ALK_wt", "EGFR_mut", "ALK_driven"]


def load(cohort: str) -> pd.DataFrame:
    df = pd.read_csv(f"{C.PROCESSED}/tcga_{cohort.lower()}_samples.csv")
    return df


def add_signatures(df: pd.DataFrame) -> pd.DataFrame:
    expr = df.set_index("sample")[[g for g in IMMUNE_MARKERS if g in df.columns]]
    for name, genes in C.SIGNATURES.items():
        df[f"sig_{name}"] = C.signature_score(expr, genes).values
    return df


def correlations(df: pd.DataFrame, cohort: str, group_label: str) -> pd.DataFrame:
    """Spearman correlation of each target vs each immune feature."""
    feats = IMMUNE_MARKERS + [f"sig_{s}" for s in C.SIGNATURES]
    rows = []
    for tgt in TARGET_NAMES:
        for f in feats:
            x = df[tgt].astype(float)
            y = df[f].astype(float)
            ok = x.notna() & y.notna()
            if ok.sum() < 10:
                continue
            rho, p = stats.spearmanr(x[ok], y[ok])
            rows.append({
                "cohort": cohort, "group": group_label, "n": int(ok.sum()),
                "target": tgt, "feature": f, "spearman_rho": rho, "p": p,
            })
    res = pd.DataFrame(rows)
    # FDR within each target x group
    res["q_BH"] = np.nan
    for (t, g), idx in res.groupby(["target", "group"]).groups.items():
        res.loc[idx, "q_BH"] = C.benjamini_hochberg(res.loc[idx, "p"].values)
    return res


def target_by_driver(df: pd.DataFrame, cohort: str) -> pd.DataFrame:
    rows = []
    sub = df[df["driver_group"].isin(DRIVER_ORDER)]
    for tgt in TARGET_NAMES:
        groups = [sub.loc[sub["driver_group"] == g, tgt].dropna() for g in DRIVER_ORDER]
        # global Kruskal-Wallis
        kw_h, kw_p = stats.kruskal(*[g for g in groups if len(g) > 0])
        rows.append({
            "cohort": cohort, "target": tgt, "test": "kruskal_all_groups",
            "stat": kw_h, "p": kw_p,
            "n_wt": len(groups[0]), "n_egfr": len(groups[1]), "n_alk": len(groups[2]),
            "median_wt": groups[0].median(), "median_egfr": groups[1].median(),
            "median_alk": groups[2].median(),
        })
        # pairwise vs WT
        for gi, gname in zip((1, 2), ("EGFR_mut", "ALK_driven")):
            if len(groups[gi]) >= 3:
                u, pu = stats.mannwhitneyu(groups[gi], groups[0], alternative="two-sided")
                # rank-biserial effect size
                n1, n2 = len(groups[gi]), len(groups[0])
                rb = 1 - 2 * u / (n1 * n2)
                rows.append({
                    "cohort": cohort, "target": tgt,
                    "test": f"MWU_{gname}_vs_wt", "stat": u, "p": pu,
                    "n_wt": n2, "n_egfr": len(groups[1]), "n_alk": len(groups[2]),
                    "median_wt": groups[0].median(), "median_egfr": groups[1].median(),
                    "median_alk": groups[2].median(), "rank_biserial": rb,
                })
    return pd.DataFrame(rows)


def plot_heatmap(res_overall: pd.DataFrame, cohort: str) -> None:
    feats = IMMUNE_MARKERS + [f"sig_{s}" for s in C.SIGNATURES]
    mat = res_overall.pivot(index="feature", columns="target", values="spearman_rho")
    mat = mat.reindex(feats)
    qmat = res_overall.pivot(index="feature", columns="target", values="q_BH").reindex(feats)
    fig, ax = plt.subplots(figsize=(4.2, 0.28 * len(feats) + 1))
    im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-0.8, vmax=0.8, aspect="auto")
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, fontsize=9)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=7)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            q = qmat.values[i, j]
            star = "*" if (pd.notna(q) and q < 0.05) else ""
            ax.text(j, i, f"{v:.2f}{star}", ha="center", va="center",
                    fontsize=6, color="black")
    ax.set_title(f"TCGA-{cohort}: target vs immune (Spearman rho)\n* BH-FDR<0.05", fontsize=9)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=7)
    fig.tight_layout()
    fig.savefig(f"{C.FIGURES}/tcga_{cohort.lower()}_target_immune_heatmap.png", dpi=150)
    plt.close(fig)


def plot_target_by_driver(df: pd.DataFrame, cohort: str) -> None:
    sub = df[df["driver_group"].isin(DRIVER_ORDER)]
    fig, axes = plt.subplots(1, len(TARGET_NAMES), figsize=(3.2 * len(TARGET_NAMES), 3.4))
    for ax, tgt in zip(np.atleast_1d(axes), TARGET_NAMES):
        data = [sub.loc[sub["driver_group"] == g, tgt].dropna().values for g in DRIVER_ORDER]
        bp = ax.boxplot(data, tick_labels=["WT", "EGFR", "ALK"], showfliers=False,
                        patch_artist=True, widths=0.6)
        for patch, col in zip(bp["boxes"], ["#bdbdbd", "#d95f02", "#1b9e77"]):
            patch.set_facecolor(col); patch.set_alpha(0.6)
        for i, d in enumerate(data):
            x = np.random.normal(i + 1, 0.06, size=len(d))
            ax.scatter(x, d, s=6, color="k", alpha=0.35, zorder=3)
        ax.set_title(f"{tgt}", fontsize=10)
        ax.set_ylabel("log2(TPM+1)", fontsize=8)
        ax.tick_params(labelsize=8)
    fig.suptitle(f"TCGA-{cohort}: target expression by driver group", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{C.FIGURES}/tcga_{cohort.lower()}_target_by_driver.png", dpi=150)
    plt.close(fig)


def main(cohort: str) -> None:
    np.random.seed(C.RANDOM_SEED)
    df = add_signatures(load(cohort))

    # per-sample signature table (small, committed)
    keep = ["sample", "driver_group", "EGFR_mut", "ALK_driven", "OS", "OS.time"] \
        + TARGET_NAMES + [f"sig_{s}" for s in C.SIGNATURES]
    df[keep].to_csv(f"{C.TABLES}/tcga_{cohort.lower()}_signature_scores.csv", index=False)

    # correlations: overall + per driver group
    all_res = [correlations(df, cohort, "ALL")]
    for g in DRIVER_ORDER:
        sub = df[df["driver_group"] == g]
        if len(sub) >= 15:
            all_res.append(correlations(sub, cohort, g))
    res = pd.concat(all_res, ignore_index=True)
    res.to_csv(f"{C.TABLES}/tcga_{cohort.lower()}_target_immune_correlations.csv", index=False)

    # target expression across driver groups
    tbd = target_by_driver(df, cohort)
    tbd.to_csv(f"{C.TABLES}/tcga_{cohort.lower()}_target_by_driver_group.csv", index=False)

    # figures
    plot_heatmap(res[res["group"] == "ALL"], cohort)
    plot_target_by_driver(df, cohort)

    # console highlights
    print(f"\n=== TCGA-{cohort} immune-axis highlights (overall) ===")
    ov = res[res["group"] == "ALL"]
    for tgt in TARGET_NAMES:
        key = ov[(ov["target"] == tgt) & (ov["feature"].isin(
            ["sig_CD8_effector", "sig_IFNG_6gene", "sig_CYT", "CD274", "sig_Checkpoint"]))]
        print(f"[{tgt}]")
        print(key[["feature", "spearman_rho", "p", "q_BH"]].to_string(index=False))
    print(f"\nwrote tables + figures for {cohort}")


if __name__ == "__main__":
    cohort = (sys.argv[1] if len(sys.argv) > 1 else "LUAD").upper()
    main(cohort)
