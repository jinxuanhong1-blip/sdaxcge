#!/usr/bin/env python3
"""Sensitivity combined tests + per-dataset figures for the hunt extension.

Reads already-computed tables (no re-download). Writes:
  results/fable_mouse_ici/sensitivity_combined.json
  results/fable_mouse_ici/primary_icb_cldn4.csv
  results/fable_mouse_ici/figures/GSE114601|157880|309199|169196_tacstd2_cldn4.png
  results/fable_mouse_ici/figures/sensitivity_tacstd2.png
"""
import json
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

PRIMARY = [
    ("GSE239485", "PolyIC+anti-PD1 vs Control vehicle", "polyIC+ICB"),
    ("GSE297630", "anti-PD-1 vs Control", "monotherapy"),
    ("E-MTAB-13704", "aPD-L1 vs vehicle", "monotherapy"),
    ("GSE114601", "anti-PD1 vs Vehicle", "monotherapy"),
    ("GSE157880", "PD-1 0Gy vs IgG 0Gy", "monotherapy"),
    ("GSE309199", "aPD1 vs Ctrl", "monotherapy"),
    ("GSE169196", "A2V+aPD1 vs IgG", "combo"),
]


def stouffer(p_up):
    p = np.clip(np.asarray(p_up, float), 1e-15, 1 - 1e-15)
    z = stats.norm.isf(p)
    zc = float(z.sum() / np.sqrt(len(z)))
    return zc, float(stats.norm.sf(zc))


def summarise(df, label):
    n = len(df)
    n_up = int((df.log2FC > 0).sum())
    sign_p = float(stats.binomtest(n_up, n, 0.5, alternative="greater").pvalue) if n else np.nan
    with_p = df.dropna(subset=["welch_p_up"])
    if len(with_p):
        zc, pz = stouffer(with_p.welch_p_up.to_numpy())
    else:
        zc, pz = np.nan, np.nan
    return dict(
        set=label, n=n, n_up=n_up, n_down=n - n_up,
        sign_test_p_greater=sign_p, stouffer_Z=zc, stouffer_p_up=pz,
        datasets=", ".join(df.dataset.tolist()),
    )


def pick(stats_df, gene):
    rows = []
    for ds, comp, kind in PRIMARY:
        hit = stats_df[(stats_df.dataset == ds) & (stats_df.gene == gene)
                       & (stats_df.comparison == comp)]
        if hit.empty:
            continue
        rec = hit.iloc[0].to_dict()
        rec["kind"] = kind
        rows.append(rec)
    return pd.DataFrame(rows)


def plot_gene_primary(df, gene, out):
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    y = np.arange(len(df))
    colors = ["#c0392b" if (pd.notna(p) and p < 0.05) else "#7f8c8d"
              for p in df.welch_p_two]
    ax.axvline(0, color="k", lw=0.8)
    ax.scatter(df.log2FC, y, c=colors, s=55, zorder=3)
    for yi, (_, r) in enumerate(df.iterrows()):
        lab = f"p={r.welch_p_two:.3g}" if pd.notna(r.welch_p_two) else ""
        ax.text(r.log2FC, yi + 0.18, lab, fontsize=7, ha="center")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.dataset} | {r.comparison}  [{r.kind}]"
                        for _, r in df.iterrows()], fontsize=7.5)
    ax.set_xlabel(f"{gene} log2FC (ICB − control)")
    ax.set_title(f"Primary ICB-vs-control {gene} (red = two-sided Welch p<0.05)")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def plot_dataset(long_df, ds):
    dsub = long_df[(long_df.dataset == ds) & (long_df.gene.isin(["Tacstd2", "Cldn4"]))]
    if dsub.empty:
        return
    groups = list(dict.fromkeys(dsub.group))
    ctrl = dsub.loc[dsub.is_control, "group"]
    if len(ctrl):
        c = ctrl.iloc[0]
        if c in groups:
            groups.remove(c)
            groups = [c] + groups
    fig, axes = plt.subplots(1, 2, figsize=(max(7, 1.15 * len(groups) * 2), 4.1))
    vt = dsub.value_type.iloc[0]
    for ax, gene in zip(axes, ["Tacstd2", "Cldn4"]):
        g = dsub[dsub.gene == gene]
        for i, grp in enumerate(groups):
            vals = g.loc[g.group == grp, "value_log2"].to_numpy()
            if len(vals) > 1:
                ax.boxplot(vals, positions=[i], widths=0.55, showfliers=False,
                           patch_artist=True,
                           boxprops=dict(facecolor="#cfe3f2", alpha=0.8),
                           medianprops=dict(color="#08519c"))
            jitter = (np.random.RandomState(i + 3).rand(len(vals)) - 0.5) * 0.18
            ax.scatter(np.full(len(vals), i) + jitter, vals, s=26,
                       color="#c0392b", zorder=3, edgecolor="white", lw=0.4)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, rotation=28, ha="right", fontsize=8)
        ax.set_title(gene, fontsize=11)
        ax.set_ylabel(f"expression ({vt})", fontsize=8)
        ax.grid(axis="y", alpha=0.25)
    fig.suptitle(f"{ds}  —  Tacstd2 / Cldn4 by treatment", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG / f"{ds}_tacstd2_cldn4.png", dpi=130)
    plt.close(fig)


def main():
    ext = pd.read_csv(RESULTS / "extended_stats.csv")
    old = pd.read_csv(RESULTS / "stats_results.csv")
    # unify
    old = old.rename(columns={
        "log2FC_treat_minus_control": "log2FC",
        "welch_p": "welch_p_two",
        "welch_t_stat": "welch_t",
    })
    if "welch_p_up" not in old.columns:
        df_approx = old["n_control"] + old["n_treat"] - 2
        old["welch_p_up"] = [
            float(stats.t.sf(t, d)) if pd.notna(t) else np.nan
            for t, d in zip(old["welch_t"], df_approx)
        ]
    keep = ["dataset", "gene", "comparison", "n_control", "n_treat",
            "mean_control", "mean_treat", "log2FC", "welch_t",
            "welch_p_two", "welch_p_up", "mwu_p"]
    stats_df = pd.concat([ext[keep], old[keep]], ignore_index=True)
    stats_df = stats_df.drop_duplicates(["dataset", "gene", "comparison"], keep="first")

    tac = pick(stats_df, "Tacstd2")
    cld = pick(stats_df, "Cldn4")
    tac.to_csv(RESULTS / "primary_icb_tacstd2.csv", index=False)
    cld.to_csv(RESULTS / "primary_icb_cldn4.csv", index=False)

    sets = {
        "all_7_primary": tac,
        "monotherapy_only": tac[tac.kind == "monotherapy"],
        "exclude_polyIC_GSE239485": tac[tac.dataset != "GSE239485"],
        "monotherapy_n_ge_3": tac[(tac.kind == "monotherapy")
                                  & (tac.n_control + tac.n_treat >= 6)],
    }
    summary = {k: summarise(v, k) for k, v in sets.items()}
    summary["cldn4_all_7_primary"] = summarise(cld, "cldn4_all_7_primary")
    summary["notes"] = {
        "GSE239485": "Poly I:C + anti-PD-1 (Poly I:C confounder); only well-powered UP",
        "GSE169196": "A2V+aPD1 combo; no aPD1-monotherapy total-tumour arm",
        "GSE157880_4Gy": ("PD-1 4Gy vs IgG 0Gy Tacstd2 +1.42 p=0.003 is radiation-"
                          "confounded: IgG 4Gy vs IgG 0Gy is +1.78 p=0.013. "
                          "Primary contrast correctly uses the 0 Gy pair."),
        "TISMO": ("49/64 p=5.8e-5 remains an external all-cancer prior; "
                  "expression matrix not independently re-downloadable."),
    }
    (RESULTS / "sensitivity_combined.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    plot_gene_primary(tac, "Tacstd2", FIG / "primary_icb_tacstd2.png")
    plot_gene_primary(cld, "Cldn4", FIG / "primary_icb_cldn4.png")

    # sensitivity forest: all-7 vs monotherapy
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
    for ax, (title, sub) in zip(axes, [
        ("All 7 primary (incl. PolyIC / combo)", tac),
        ("Monotherapy only (n=5)", tac[tac.kind == "monotherapy"]),
    ]):
        sub = sub.reset_index(drop=True)
        y = np.arange(len(sub))
        colors = ["#c0392b" if (pd.notna(p) and p < 0.05) else "#7f8c8d"
                  for p in sub.welch_p_two]
        ax.axvline(0, color="k", lw=0.8)
        ax.scatter(sub.log2FC, y, c=colors, s=50, zorder=3)
        for yi, (_, r) in enumerate(sub.iterrows()):
            ax.text(r.log2FC, yi + 0.18, f"p={r.welch_p_two:.3g}", fontsize=7, ha="center")
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r.dataset}" for _, r in sub.iterrows()], fontsize=8)
        ax.set_xlabel("Tacstd2 log2FC (ICB − control)")
        ax.set_title(title, fontsize=10)
        ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "sensitivity_tacstd2.png", dpi=140)
    plt.close(fig)

    long_df = pd.read_csv(RESULTS / "extended_per_sample.csv")
    for ds in ["GSE114601", "GSE157880", "GSE309199", "GSE169196"]:
        plot_dataset(long_df, ds)


if __name__ == "__main__":
    main()
