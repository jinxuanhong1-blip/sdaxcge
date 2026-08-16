#!/usr/bin/env python3
"""Compositional tests: TACSTD2-high vs T/NK drop; NMPR vs higher epi / lower T/NK.

Primary method: CLR + ILR (Aitchison) on 6-part lineage counts, plus a
Dirichlet-multinomial two-group LRT. scCODA was not installable (rpy2/R).

Patient is the unit. All n and p are computed here from public matrices.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compositional import (  # noqa: E402
    LINEAGE_COLS,
    attach_transforms,
    cohort_center,
    dm_two_group,
    mwu,
    spearman,
)

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = RES / "figures"
MIN_EPI_TAC = 20
MIN_CELLS = 200

# User trend
# H1: T/NK drops when malignant TACSTD2 is high  (rho < 0; high < low)
# H2: NMPR has higher epithelial and lower T/NK


def load_table() -> pd.DataFrame:
    df = pd.read_csv(RES / "composition_all.tsv", sep="\t")
    df = attach_transforms(df)
    df["group"] = df["group"].fillna("").astype(str)
    df["tacstd2_high"] = np.nan
    return df


def eligible_tac(df: pd.DataFrame) -> pd.DataFrame:
    return df[(df["n_cells"] >= MIN_CELLS) & (df["n_epi"] >= MIN_EPI_TAC)].copy()


def eligible_resp(df: pd.DataFrame) -> pd.DataFrame:
    return df[(df["n_cells"] >= MIN_CELLS) & df["group"].isin(["MPR", "NMPR"])].copy()


def median_split(df: pd.DataFrame, col: str = "malignant_TACSTD2_mean_log1p") -> pd.DataFrame:
    out = df.copy()
    med = float(out[col].median())
    out["tacstd2_high"] = np.where(out[col] >= med, "high", "low")
    out["tacstd2_median"] = med
    return out


def h1_tests(df: pd.DataFrame, label: str) -> dict:
    """T/NK vs malignant TACSTD2 (continuous + median split)."""
    a = eligible_tac(df)
    out = {"label": label, "n": int(len(a)), "test": "H1_TACSTD2_vs_TNK"}
    if len(a) < 4:
        out["note"] = "n<4 after min_epi=20"
        return out
    a = median_split(a)
    out["n_high"] = int((a["tacstd2_high"] == "high").sum())
    out["n_low"] = int((a["tacstd2_high"] == "low").sum())
    out["tacstd2_median"] = float(a["tacstd2_median"].iloc[0])
    out["spearman_frac_TNK"] = spearman(a["malignant_TACSTD2_mean_log1p"], a["frac_T_NK"])
    out["spearman_clr_TNK"] = spearman(a["malignant_TACSTD2_mean_log1p"], a["clr_T_NK"])
    out["spearman_ilr_TNK"] = spearman(a["malignant_TACSTD2_mean_log1p"], a["ilr_T_NK"])
    hi = a[a.tacstd2_high == "high"]
    lo = a[a.tacstd2_high == "low"]
    # a = high, b = low; expected high < low for T/NK  -> alternative less
    out["mwu_frac_TNK_high_vs_low"] = mwu(hi["frac_T_NK"], lo["frac_T_NK"], "less")
    out["mwu_clr_TNK_high_vs_low"] = mwu(hi["clr_T_NK"], lo["clr_T_NK"], "less")
    out["mwu_ilr_TNK_high_vs_low"] = mwu(hi["ilr_T_NK"], lo["ilr_T_NK"], "less")
    out["mwu_frac_TNK_two_sided"] = mwu(hi["frac_T_NK"], lo["frac_T_NK"], "two-sided")
    out["direction_rho_negative"] = bool(out["spearman_clr_TNK"]["rho"] < 0)
    out["direction_high_lower_TNK"] = bool(hi["frac_T_NK"].median() < lo["frac_T_NK"].median())
    out["matches_H1"] = bool(out["direction_rho_negative"] and out["direction_high_lower_TNK"])
    # DM: high vs low
    if len(hi) >= 2 and len(lo) >= 2:
        out["dm_high_vs_low"] = dm_two_group(
            a[LINEAGE_COLS].to_numpy(),
            a["tacstd2_high"].to_numpy(),
            n_perm=500 if len(a) <= 20 else 200,
            seed=7,
        )
    return out


def h2_tests(df: pd.DataFrame, label: str) -> dict:
    """NMPR vs MPR: higher epi, lower T/NK."""
    a = eligible_resp(df)
    out = {"label": label, "n": int(len(a)), "test": "H2_NMPR_vs_MPR_composition"}
    nm = a[a.group == "NMPR"]
    mp = a[a.group == "MPR"]
    out["n_NMPR"] = int(len(nm))
    out["n_MPR"] = int(len(mp))
    if len(nm) < 2 or len(mp) < 2:
        out["note"] = "need >=2 per response group"
        return out
    # a = NMPR, b = MPR
    out["mwu_frac_Epi_NMPR_gt_MPR"] = mwu(nm["frac_Epithelial"], mp["frac_Epithelial"], "greater")
    out["mwu_clr_Epi_NMPR_gt_MPR"] = mwu(nm["clr_Epithelial"], mp["clr_Epithelial"], "greater")
    out["mwu_ilr_Epi_NMPR_gt_MPR"] = mwu(nm["ilr_Epithelial"], mp["ilr_Epithelial"], "greater")
    out["mwu_frac_TNK_NMPR_lt_MPR"] = mwu(nm["frac_T_NK"], mp["frac_T_NK"], "less")
    out["mwu_clr_TNK_NMPR_lt_MPR"] = mwu(nm["clr_T_NK"], mp["clr_T_NK"], "less")
    out["mwu_ilr_TNK_NMPR_lt_MPR"] = mwu(nm["ilr_T_NK"], mp["ilr_T_NK"], "less")
    out["mwu_frac_Epi_two_sided"] = mwu(nm["frac_Epithelial"], mp["frac_Epithelial"], "two-sided")
    out["mwu_frac_TNK_two_sided"] = mwu(nm["frac_T_NK"], mp["frac_T_NK"], "two-sided")
    out["direction_NMPR_higher_epi"] = bool(nm["frac_Epithelial"].median() > mp["frac_Epithelial"].median())
    out["direction_NMPR_lower_TNK"] = bool(nm["frac_T_NK"].median() < mp["frac_T_NK"].median())
    out["matches_H2"] = bool(out["direction_NMPR_higher_epi"] and out["direction_NMPR_lower_TNK"])
    out["dm_NMPR_vs_MPR"] = dm_two_group(
        a[LINEAGE_COLS].to_numpy(),
        a["group"].to_numpy(),
        n_perm=500 if len(a) <= 20 else 200,
        seed=11,
    )
    # residual-tumor confound note
    out["median_n_epi_NMPR"] = float(nm["n_epi"].median())
    out["median_n_epi_MPR"] = float(mp["n_epi"].median())
    return out


def combo_ols(df: pd.DataFrame, label: str) -> dict:
    """CLR(T/NK) ~ TACSTD2 + C(cohort) on eligible TACSTD2 rows."""
    a = eligible_tac(df)
    out = {"label": label, "n": int(len(a))}
    if len(a) < 6 or a["cohort"].nunique() < 2:
        out["note"] = "skip OLS (need multi-cohort n>=6)"
        return out
    try:
        fit = smf.ols("clr_T_NK ~ malignant_TACSTD2_mean_log1p + C(cohort)", data=a).fit()
        out["beta_TACSTD2"] = float(fit.params["malignant_TACSTD2_mean_log1p"])
        out["p_TACSTD2"] = float(fit.pvalues["malignant_TACSTD2_mean_log1p"])
        out["r2"] = float(fit.rsquared)
        out["direction_negative"] = bool(out["beta_TACSTD2"] < 0)
    except Exception as e:
        out["note"] = str(e)
    a2 = eligible_resp(df)
    if len(a2) >= 6 and a2["cohort"].nunique() >= 2 and set(a2["group"]) >= {"MPR", "NMPR"}:
        try:
            fit2 = smf.ols("clr_Epithelial ~ C(group, Treatment('MPR')) + C(cohort)", data=a2).fit()
            key = [k for k in fit2.params.index if "group" in k][0]
            out["beta_NMPR_vs_MPR_clrEpi"] = float(fit2.params[key])
            out["p_NMPR_vs_MPR_clrEpi"] = float(fit2.pvalues[key])
            fit3 = smf.ols("clr_T_NK ~ C(group, Treatment('MPR')) + C(cohort)", data=a2).fit()
            key3 = [k for k in fit3.params.index if "group" in k][0]
            out["beta_NMPR_vs_MPR_clrTNK"] = float(fit3.params[key3])
            out["p_NMPR_vs_MPR_clrTNK"] = float(fit3.pvalues[key3])
        except Exception as e:
            out["ols_response_note"] = str(e)
    # cohort-centered Spearman
    a = cohort_center(a, ["malignant_TACSTD2_mean_log1p", "clr_T_NK", "frac_T_NK"])
    out["spearman_cc_clr_TNK"] = spearman(a["malignant_TACSTD2_mean_log1p__cc"], a["clr_T_NK__cc"])
    return out


def slice_sets(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    post = df["timepoint"].astype(str).str.contains("Post", case=False, na=False)
    return {
        "GSE207422_post": df[(df.dataset == "GSE207422") & post],
        "GSE207422_all_labeled": df[(df.dataset == "GSE207422") & df.group.isin(["MPR", "NMPR"])],
        "GSE241934_IIT": df[df.cohort == "GSE241934_IIT"],
        "GSE241934_REAL": df[df.cohort == "GSE241934_REAL"],
        "GSE241934": df[df.dataset == "GSE241934"],
        "GSE291670": df[df.dataset == "GSE291670"],
        "GSE207422_post+GSE241934": df[((df.dataset == "GSE207422") & post) | (df.dataset == "GSE241934")],
        "GSE207422_post+GSE241934_REAL": df[((df.dataset == "GSE207422") & post) | (df.cohort == "GSE241934_REAL")],
        "GSE207422_post+GSE241934_IIT": df[((df.dataset == "GSE207422") & post) | (df.cohort == "GSE241934_IIT")],
        "GSE207422_post+GSE291670": df[((df.dataset == "GSE207422") & post) | (df.dataset == "GSE291670")],
        "GSE241934+GSE291670": df[df.dataset.isin(["GSE241934", "GSE291670"])],
        "all_three_post": df[((df.dataset == "GSE207422") & post) | (df.dataset.isin(["GSE241934", "GSE291670"]))],
    }


def summarize_row(name: str, h1: dict, h2: dict, ols: dict) -> dict:
    rho = h1.get("spearman_clr_TNK", {}).get("rho", float("nan"))
    p_rho = h1.get("spearman_clr_TNK", {}).get("p", float("nan"))
    p_h1 = h1.get("mwu_clr_TNK_high_vs_low", {}).get("p", float("nan"))
    p_epi = h2.get("mwu_frac_Epi_NMPR_gt_MPR", {}).get("p", float("nan"))
    p_tnk = h2.get("mwu_frac_TNK_NMPR_lt_MPR", {}).get("p", float("nan"))
    match = bool(h1.get("matches_H1") and h2.get("matches_H2"))
    return {
        "combo": name,
        "n_H1": h1.get("n"),
        "n_H2_NMPR": h2.get("n_NMPR"),
        "n_H2_MPR": h2.get("n_MPR"),
        "H1_clr_rho": rho,
        "H1_clr_rho_p": p_rho,
        "H1_high_lt_low_clr_p": p_h1,
        "H1_matches": h1.get("matches_H1"),
        "H2_NMPR_higher_epi": h2.get("direction_NMPR_higher_epi"),
        "H2_NMPR_lower_TNK": h2.get("direction_NMPR_lower_TNK"),
        "H2_epi_p_greater": p_epi,
        "H2_tnk_p_less": p_tnk,
        "H2_matches": h2.get("matches_H2"),
        "both_match_user_trend": match,
        "ols_beta_TACSTD2": ols.get("beta_TACSTD2"),
        "ols_p_TACSTD2": ols.get("p_TACSTD2"),
        "ols_cc_rho": ols.get("spearman_cc_clr_TNK", {}).get("rho") if isinstance(ols.get("spearman_cc_clr_TNK"), dict) else None,
        "ols_cc_p": ols.get("spearman_cc_clr_TNK", {}).get("p") if isinstance(ols.get("spearman_cc_clr_TNK"), dict) else None,
    }


def pick_matching_combo(summary: pd.DataFrame) -> dict:
    """Among combos matching both directions, prefer largest H1 n; else none."""
    hit = summary[summary["both_match_user_trend"] == True]  # noqa: E712
    if hit.empty:
        # partial: H1 match only, then H2
        h1 = summary[summary["H1_matches"] == True]  # noqa: E712
        return {
            "combo": None,
            "note": "no combination matches both H1 and H2 directions",
            "H1_only": h1["combo"].tolist(),
            "H2_only": summary[summary["H2_matches"] == True]["combo"].tolist(),  # noqa: E712
        }
    hit = hit.sort_values(["n_H1", "n_H2_NMPR"], ascending=False)
    top = hit.iloc[0]
    return {
        "combo": top["combo"],
        "n_H1": int(top["n_H1"]),
        "n_H2": f"{int(top['n_H2_NMPR'])} vs {int(top['n_H2_MPR'])}",
        "H1_clr_rho": float(top["H1_clr_rho"]),
        "H1_clr_rho_p": float(top["H1_clr_rho_p"]),
        "H2_epi_p": float(top["H2_epi_p_greater"]),
        "H2_tnk_p": float(top["H2_tnk_p_less"]),
        "all_matching_combos": hit["combo"].tolist(),
        "note": "largest-n combo whose H1 and H2 directions both match the user trend",
    }


def plot_all(df: pd.DataFrame, summary: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    # Figure 1: TACSTD2 vs frac T/NK, post-tx / all labeled
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6), sharey=False)
    panels = [
        ("GSE207422", df[(df.dataset == "GSE207422") & df.timepoint.astype(str).str.contains("Post", na=False)]),
        ("GSE241934", df[df.dataset == "GSE241934"]),
        ("GSE291670", df[df.dataset == "GSE291670"]),
    ]
    for ax, (title, sub) in zip(axes, panels):
        sub = eligible_tac(sub)
        for g, col, m in [("NMPR", "#b23a48", "o"), ("MPR", "#2b6cb0", "s"), ("", "#888888", "x")]:
            ss = sub[sub.group == g] if g else sub[~sub.group.isin(["MPR", "NMPR"])]
            if ss.empty:
                continue
            ax.scatter(
                ss["malignant_TACSTD2_mean_log1p"],
                ss["frac_T_NK"],
                c=col,
                marker=m,
                s=36,
                label=g or "unlabeled",
                edgecolors="k",
                linewidths=0.3,
            )
        sp = spearman(sub["malignant_TACSTD2_mean_log1p"], sub["frac_T_NK"])
        ax.set_title(f"{title}\nn={sp['n']}  ρ={sp['rho']:.2f}  p={sp['p']:.3g}")
        ax.set_xlabel("malignant TACSTD2  mean log1p(UMI)")
        ax.set_ylabel("T/NK fraction")
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_tacstd2_vs_tnk_fraction.png", dpi=160)
    fig.savefig(FIG / "fig1_tacstd2_vs_tnk_fraction.pdf")
    plt.close(fig)

    # Figure 2: response boxplots
    fig, axes = plt.subplots(2, 3, figsize=(11.2, 6.4))
    metrics = [("frac_Epithelial", "epithelial fraction"), ("frac_T_NK", "T/NK fraction")]
    for row, (col, ylab) in enumerate(metrics):
        for ax, (title, sub) in zip(axes[row], panels):
            sub = eligible_resp(sub)
            data, labs, cols = [], [], []
            for g, c in [("MPR", "#2b6cb0"), ("NMPR", "#b23a48")]:
                v = sub.loc[sub.group == g, col].dropna().values
                data.append(v)
                labs.append(f"{g}\nn={len(v)}")
                cols.append(c)
            bp = ax.boxplot(data, tick_labels=labs, patch_artist=True, widths=0.55)
            for patch, c in zip(bp["boxes"], cols):
                patch.set_facecolor(c)
                patch.set_alpha(0.45)
            ax.set_ylabel(ylab)
            ax.set_title(title if row == 0 else "")
    fig.suptitle("H2: NMPR vs MPR composition (patient-level)", y=1.01)
    fig.tight_layout()
    fig.savefig(FIG / "fig2_nmpr_vs_mpr_fractions.png", dpi=160)
    fig.savefig(FIG / "fig2_nmpr_vs_mpr_fractions.pdf")
    plt.close(fig)

    # Figure 3: stacked composition for GSE241934 (author labels, best n)
    sub = df[df.dataset == "GSE241934"].sort_values(["group", "frac_Epithelial"], ascending=[True, False])
    fig, ax = plt.subplots(figsize=(12.5, 4.2))
    x = np.arange(len(sub))
    bottom = np.zeros(len(sub))
    colors = {
        "Epithelial": "#d4a017",
        "T_NK": "#2b6cb0",
        "B": "#6b4c9a",
        "Myeloid": "#3d8b6e",
        "Stromal": "#c46b3a",
        "Other": "#bbbbbb",
    }
    for name in LINEAGE_COLS:
        ax.bar(x, sub[f"frac_{name}"], bottom=bottom, color=colors[name], width=0.9, label=name)
        bottom = bottom + sub[f"frac_{name}"].to_numpy()
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{r.sample_id}\n{r.group}" for r in sub.itertuples()],
        fontsize=6,
        rotation=90,
    )
    ax.set_ylabel("fraction")
    ax.set_title("GSE241934 author major.cell.type (IIT + REAL)")
    ax.legend(ncol=6, frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "fig3_gse241934_stacked.png", dpi=160)
    fig.savefig(FIG / "fig3_gse241934_stacked.pdf")
    plt.close(fig)

    # Figure 4: combo match heatmap-ish table figure
    fig, ax = plt.subplots(figsize=(10.8, 4.8))
    ax.axis("off")
    show = summary.copy()
    cells = []
    for _, r in show.iterrows():
        cells.append(
            [
                r["combo"],
                f"{r['n_H1']}",
                f"{r['n_H2_NMPR']} vs {r['n_H2_MPR']}",
                f"{r['H1_clr_rho']:.2f}" if pd.notna(r["H1_clr_rho"]) else "NA",
                f"{r['H1_clr_rho_p']:.3g}" if pd.notna(r["H1_clr_rho_p"]) else "NA",
                "yes" if r["H1_matches"] else "no",
                "yes" if r["H2_matches"] else "no",
                "YES" if r["both_match_user_trend"] else "no",
            ]
        )
    tbl = ax.table(
        cellText=cells,
        colLabels=["combo", "n H1", "NMPR vs MPR", "CLR ρ", "ρ p", "H1 dir", "H2 dir", "both"],
        loc="center",
        cellLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7)
    tbl.scale(1.05, 1.25)
    ax.set_title("Direction match to user trend (not a significance claim)")
    fig.tight_layout()
    fig.savefig(FIG / "fig4_combo_direction_table.png", dpi=160)
    fig.savefig(FIG / "fig4_combo_direction_table.pdf")
    plt.close(fig)


def write_finding(summary: pd.DataFrame, picked: dict, payload: dict) -> None:
    lines = []
    lines.append("# FINDING — compositional analysis (scCODA fallback)")
    lines.append("")
    lines.append("Public only. Patient is the unit. scCODA was **not installable**")
    lines.append("(pip `sccoda` requires rpy2 → system R). Fallback: CLR + ILR on")
    lines.append("6-part lineage fractions, plus Dirichlet-multinomial two-group LRT")
    lines.append("with permutation p.")
    lines.append("")
    lines.append("User trend tested")
    lines.append("")
    lines.append("1. **H1** T/NK fraction drops when patient malignant TACSTD2 is high.")
    lines.append("2. **H2** NMPR has higher epithelial / lower T/NK.")
    lines.append("")
    lines.append("## Combination that matches the trend")
    lines.append("")
    if picked.get("combo"):
        lines.append(f"**{picked['combo']}** is the only slice whose *directions* match both H1 and H2.")
        lines.append("")
        lines.append(f"- H1 (CLR T/NK vs TACSTD2): n={picked['n_H1']}, ρ={picked['H1_clr_rho']:.3f}, p={picked['H1_clr_rho_p']:.3g}")
        lines.append(f"- H2: NMPR vs MPR {picked['n_H2']}; epi one-sided p={picked['H2_epi_p']:.3g}; T/NK one-sided p={picked['H2_tnk_p']:.3g}")
        lines.append(f"- All direction-matching slices: {', '.join(picked['all_matching_combos'])}")
        lines.append("")
        lines.append("Matching means **sign of the effect**, not p<0.05. On this winning")
        lines.append("slice H1 is a null (ρ ≈ 0). No two- or three-series combination")
        lines.append("keeps both directions: adding GSE207422 flips H1 to positive.")
    else:
        lines.append("**No combination matches both H1 and H2 directions.**")
        lines.append(f"- H1-only: {picked.get('H1_only')}")
        lines.append(f"- H2-only: {picked.get('H2_only')}")
    lines.append("")
    lines.append("## Per-cohort honest n / p (primary)")
    lines.append("")
    lines.append("| Cohort | H1 n | CLR ρ (TACSTD2 vs T/NK) | ρ p | H1 dir | H2 n (NMPR vs MPR) | epi p (NMPR>) | T/NK p (NMPR<) | H2 dir |")
    lines.append("|---|---:|---:|---:|---|---:|---:|---:|---|")
    primary = [
        "GSE207422_post",
        "GSE241934_IIT",
        "GSE241934_REAL",
        "GSE241934",
        "GSE291670",
        "GSE207422_post+GSE241934",
        "GSE207422_post+GSE241934_REAL",
        "GSE207422_post+GSE291670",
        "GSE241934+GSE291670",
        "all_three_post",
    ]
    for name in primary:
        r = summary[summary.combo == name]
        if r.empty:
            continue
        r = r.iloc[0]
        def fmt(x, nd=3):
            return "NA" if pd.isna(x) else (f"{x:.{nd}g}" if abs(x) < 0.01 or abs(x) >= 100 else f"{x:.3f}")
        lines.append(
            f"| {name} | {r['n_H1']} | {fmt(r['H1_clr_rho'])} | {fmt(r['H1_clr_rho_p'], 3)} | "
            f"{'yes' if r['H1_matches'] else 'no'} | {r['n_H2_NMPR']} vs {r['n_H2_MPR']} | "
            f"{fmt(r['H2_epi_p_greater'])} | {fmt(r['H2_tnk_p_less'])} | "
            f"{'yes' if r['H2_matches'] else 'no'} |"
        )
    lines.append("")
    lines.append("## Confounds that stay in the write-up")
    lines.append("")
    lines.append("- Post-neoadjuvant **MPR/pCR residuals have fewer epithelial cells by definition**.")
    lines.append("  H2 epithelial enrichment in NMPR is partly that residual-tumor fact,")
    lines.append("  especially in GSE241934.")
    lines.append("- GSE207422 / GSE291670 use a marker hierarchy (CopyKAT / author barcodes")
    lines.append("  are not on GEO). GSE241934 uses author `major.cell.type`.")
    lines.append("- GSE291670 is n=3 vs 3. Exact Wilcoxon cannot go below p=0.10 two-sided.")
    lines.append("- Cell-level p-values are not reported (pseudoreplication).")
    lines.append("- scCODA credible effects are **not** claimed.")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `composition_all.tsv` — patient-level counts + TACSTD2")
    lines.append("- `tests_summary.tsv` / `tests_full.json`")
    lines.append("- `figures/`")
    (RES / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    df = load_table()
    sets = slice_sets(df)
    full = {}
    rows = []
    for name, sub in sets.items():
        print("testing", name, "n=", len(sub), flush=True)
        h1 = h1_tests(sub, name)
        h2 = h2_tests(sub, name)
        ols = combo_ols(sub, name)
        full[name] = {"H1": h1, "H2": h2, "OLS": ols}
        rows.append(summarize_row(name, h1, h2, ols))
    summary = pd.DataFrame(rows)
    summary.to_csv(RES / "tests_summary.tsv", sep="\t", index=False)
    picked = pick_matching_combo(summary)
    payload = {
        "sccoda_installable": False,
        "sccoda_reason": "pip sccoda requires rpy2, which requires system R; not present",
        "method": "CLR + ILR (6-part) + Dirichlet-multinomial two-group LRT (permutation p)",
        "min_epi_for_TACSTD2": MIN_EPI_TAC,
        "min_cells": MIN_CELLS,
        "pseudocount": 0.5,
        "picked_combo": picked,
        "tests": full,
    }
    (RES / "tests_full.json").write_text(json.dumps(payload, indent=2, default=str))
    (RES / "picked_combo.json").write_text(json.dumps(picked, indent=2))
    plot_all(df, summary)
    write_finding(summary, picked, payload)
    print(json.dumps(picked, indent=2))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
