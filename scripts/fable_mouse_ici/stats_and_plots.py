#!/usr/bin/env python3
"""Statistics and figures for Tacstd2/Cldn4 by treatment.

Reads results/fable_mouse_ici/per_sample_expression.csv. For every dataset with
biological replicates in both the control and a treatment arm, runs a Welch
two-sample t-test and a Mann-Whitney U test (two-sided) for each target gene and
each treatment-vs-control contrast, and applies Benjamini-Hochberg FDR across all
contrasts within a dataset. Datasets with a single sample per condition are
summarised descriptively (no inferential p-values). Produces per-dataset figures
and a combined overview.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

RESULTS = Path(__file__).resolve().parents[2] / "results" / "fable_mouse_ici"
FIG = RESULTS / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# control label per dataset (everything else = treatment arm)
CONTROL = {
    "GSE239485": "Control vehicle",
    "GSE297630": "Control",
    "GSE241978": "Control",
    "GSE330658": "Control",
    "E-MTAB-13704": "vehicle",
    "GSE197260": None,   # 1/arm -> descriptive
    "GSE133604": None,
    "GSE129297": None,
    "GSE297632": None,
    "GSE222158": None,
}


def bh_fdr(pvals):
    p = np.asarray(pvals, float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    idx = np.where(ok)[0]
    if len(idx) == 0:
        return q
    pv = p[idx]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(pv)
    adj = ranked * n / (np.arange(1, n + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    q[idx] = out
    return q


def run_stats(df):
    rows = []
    for ds, dsub in df.groupby("dataset"):
        ctrl = CONTROL[ds]
        for gene, gsub in dsub.groupby("gene"):
            if ctrl is None or ctrl not in set(gsub.group):
                continue
            cvals = gsub.loc[gsub.group == ctrl, "value_log2"].to_numpy()
            for grp, arm in gsub.groupby("group"):
                if grp == ctrl:
                    continue
                tvals = arm["value_log2"].to_numpy()
                n1, n2 = len(cvals), len(tvals)
                rec = dict(dataset=ds, gene=gene, comparison=f"{grp} vs {ctrl}",
                           n_control=n1, n_treat=n2,
                           mean_control=round(float(np.mean(cvals)), 4),
                           mean_treat=round(float(np.mean(tvals)), 4),
                           sd_control=round(float(np.std(cvals, ddof=1)), 4) if n1 > 1 else np.nan,
                           sd_treat=round(float(np.std(tvals, ddof=1)), 4) if n2 > 1 else np.nan,
                           log2FC_treat_minus_control=round(float(np.mean(tvals) - np.mean(cvals)), 4),
                           value_type=arm["value_type"].iloc[0])
                if n1 >= 2 and n2 >= 2:
                    t = stats.ttest_ind(tvals, cvals, equal_var=False)
                    rec["welch_t_stat"] = round(float(t.statistic), 4)
                    rec["welch_p"] = float(t.pvalue)
                    try:
                        u = stats.mannwhitneyu(tvals, cvals, alternative="two-sided")
                        rec["mwu_p"] = float(u.pvalue)
                    except ValueError:
                        rec["mwu_p"] = np.nan
                else:
                    rec["welch_t_stat"] = np.nan
                    rec["welch_p"] = np.nan
                    rec["mwu_p"] = np.nan
                rows.append(rec)
    res = pd.DataFrame(rows)
    # BH-FDR on Welch p within each dataset (across genes x arms)
    res["welch_FDR_within_dataset"] = np.nan
    for ds, sub in res.groupby("dataset"):
        res.loc[sub.index, "welch_FDR_within_dataset"] = bh_fdr(sub["welch_p"].to_numpy())
    return res


def descriptive(df):
    rows = []
    for ds in [d for d, c in CONTROL.items() if c is None]:
        dsub = df[df.dataset == ds]
        for (gene, grp), s in dsub.groupby(["gene", "group"]):
            rows.append(dict(dataset=ds, gene=gene, group=grp,
                             n=len(s),
                             mean_log2=round(float(s.value_log2.mean()), 4),
                             value_type=s.value_type.iloc[0]))
    return pd.DataFrame(rows)


def plot_dataset(df, ds):
    dsub = df[df.dataset == ds]
    genes = ["Tacstd2", "Cldn4"]
    groups = list(dict.fromkeys(dsub.group))  # preserve order
    ctrl = CONTROL[ds]
    if ctrl in groups:
        groups.remove(ctrl)
        groups = [ctrl] + groups
    fig, axes = plt.subplots(1, 2, figsize=(max(7, 1.2 * len(groups) * 2), 4.2))
    vt = dsub.value_type.iloc[0]
    for ax, gene in zip(axes, genes):
        g = dsub[dsub.gene == gene]
        data = [g.loc[g.group == grp, "value_log2"].to_numpy() for grp in groups]
        positions = np.arange(len(groups))
        # box only where >1 point
        for i, (grp, vals) in enumerate(zip(groups, data)):
            if len(vals) > 1:
                ax.boxplot(vals, positions=[i], widths=0.55,
                           showfliers=False, patch_artist=True,
                           boxprops=dict(facecolor="#cfe3f2", alpha=0.8),
                           medianprops=dict(color="#08519c"))
            jitter = (np.random.RandomState(i).rand(len(vals)) - 0.5) * 0.18
            ax.scatter(np.full(len(vals), i) + jitter, vals, s=26,
                       color="#c0392b", zorder=3, edgecolor="white", linewidth=0.4)
        ax.set_xticks(positions)
        ax.set_xticklabels(groups, rotation=30, ha="right", fontsize=8)
        ax.set_title(gene, fontsize=11)
        ax.set_ylabel(f"expression ({vt})", fontsize=8)
        ax.grid(axis="y", alpha=0.25)
    fig.suptitle(f"{ds}  —  Tacstd2 / Cldn4 by treatment", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG / f"{ds}_tacstd2_cldn4.png", dpi=130)
    plt.close(fig)


def plot_overview(res):
    """Forest-style overview of log2FC for replicate datasets."""
    r = res.dropna(subset=["welch_p"]).copy()
    if r.empty:
        return
    r["label"] = r.dataset + " | " + r.comparison
    fig, axes = plt.subplots(1, 2, figsize=(13, max(4, 0.42 * len(r))))
    for ax, gene in zip(axes, ["Tacstd2", "Cldn4"]):
        gg = r[r.gene == gene].reset_index(drop=True)
        y = np.arange(len(gg))
        colors = ["#c0392b" if p < 0.05 else "#7f8c8d" for p in gg.welch_p]
        ax.scatter(gg.log2FC_treat_minus_control, y, color=colors, s=40, zorder=3)
        for yi, (_, row) in zip(y, gg.iterrows()):
            ax.text(row.log2FC_treat_minus_control, yi + 0.15,
                    f"p={row.welch_p:.3g}", fontsize=6, ha="center")
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(gg.label, fontsize=7)
        ax.set_xlabel("log2 fold-change (treatment − control)", fontsize=9)
        ax.set_title(gene, fontsize=11)
        ax.grid(axis="x", alpha=0.25)
    fig.suptitle("Tacstd2 / Cldn4: treatment vs control (replicate datasets; red = Welch p<0.05)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(FIG / "overview_log2FC.png", dpi=140)
    plt.close(fig)


def main():
    df = pd.read_csv(RESULTS / "per_sample_expression.csv")
    res = run_stats(df)
    res.to_csv(RESULTS / "stats_results.csv", index=False)
    desc = descriptive(df)
    desc.to_csv(RESULTS / "descriptive_single_replicate.csv", index=False)

    for ds in CONTROL:
        plot_dataset(df, ds)
    plot_overview(res)

    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 30)
    print("=== inferential stats (replicate datasets) ===")
    show = ["dataset", "gene", "comparison", "n_control", "n_treat",
            "mean_control", "mean_treat", "log2FC_treat_minus_control",
            "welch_p", "mwu_p", "welch_FDR_within_dataset"]
    print(res[show].to_string(index=False))
    print("\n=== descriptive (single sample per condition) ===")
    print(desc.to_string(index=False))


if __name__ == "__main__":
    main()
