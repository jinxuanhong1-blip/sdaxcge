#!/usr/bin/env python3
"""ADDITIVE CLDN4-only ADC+ICI quadrant counts.

No dual-high. TACSTD2 is never a gate.
On GSE285029, GSE218989, GSE126044, GSE166449: count
  CLDN4 Q4 ∩ (IFN-γ Q4 ∪ CD274 Q4)   vs   CLDN4 Q4 ∩ IFN-low
Honest n. Extra scatter.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
HARVEST = HERE / "harvested"
DATA = HERE / "data"
FIG = HERE / "figures"
TAB = HERE / "tables"

# A11 a-priori IFN-γ compact (not mined here).
IFN_COMPACT = ["IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "IDO1", "GBP1"]
AYERS6 = ["IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG"]


def zmean(log_expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in log_expr.index]
    if not present:
        return pd.Series(np.nan, index=log_expr.columns), present
    block = log_expr.loc[present]
    mu = block.mean(axis=1)
    sd = block.std(axis=1, ddof=0).replace(0, np.nan)
    z = block.sub(mu, axis=0).div(sd, axis=0)
    return z.mean(axis=0), present


def quartiles(s: pd.Series) -> pd.Series:
    """Equal-count Q1–Q4 via rank(method='first') then qcut. Ties still fill four bins."""
    r = s.rank(method="first")
    return pd.qcut(r, 4, labels=["Q1", "Q2", "Q3", "Q4"])


def q_cut_value(s: pd.Series, q: float) -> float:
    return float(s.quantile(q))


def fisher_q4_overlap(a: pd.Series, b: pd.Series) -> dict:
    d = pd.concat([a, b], axis=1).dropna()
    d.columns = ["a", "b"]
    qa = quartiles(d["a"])
    qb = quartiles(d["b"])
    both = int(((qa == "Q4") & (qb == "Q4")).sum())
    a4 = int((qa == "Q4").sum())
    b4 = int((qb == "Q4").sum())
    n = int(len(d))
    table = np.array(
        [
            [both, a4 - both],
            [b4 - both, n - a4 - b4 + both],
        ]
    )
    if table.min() < 0:
        or_, p = np.nan, np.nan
    else:
        or_, p = stats.fisher_exact(table, alternative="two-sided")
    return {
        "n": n,
        "n_anchor_q4": a4,
        "n_feature_q4": b4,
        "n_both_q4": both,
        "frac_anchor_q4_in_feature_q4": both / a4 if a4 else np.nan,
        "expected_frac": a4 * b4 / n / a4 if a4 and n else np.nan,
        "expected_n_both": a4 * b4 / n if n else np.nan,
        "odds_ratio": float(or_) if or_ == or_ else np.nan,
        "fisher_p": float(p) if p == p else np.nan,
    }


def load_gse285029() -> pd.DataFrame:
    df = pd.read_csv(HARVEST / "GSE285029_sample_scores.csv", index_col=0)
    out = pd.DataFrame(
        {
            "sample": df.index.astype(str),
            "CLDN4": df["CLDN4"].astype(float),
            "CD274": df["CD274"].astype(float),
            "IFN": df["IFN_compact"].astype(float),
            "IFNG": df["IFNG"].astype(float),
            "IFN_sens": df["HALLMARK_IFNG"].astype(float),
        }
    )
    out["response"] = np.nan
    out["cohort"] = "GSE285029"
    return out


def load_gse218989() -> pd.DataFrame:
    h = pd.read_csv(HARVEST / "GSE218989_per_patient.tsv", sep="\t")
    ext = pd.read_csv(DATA / "GSE218989_gene_extract.tsv", sep="\t", index_col=0)
    ext.index = ext.index.astype(str)
    logx = np.log2(ext.clip(lower=0) + 1.0)
    ifn, used = zmean(logx, IFN_COMPACT)
    cd274 = logx.loc["CD274"]
    # Join on patient IDs (TPM columns).
    h = h.set_index("patient")
    out = pd.DataFrame(
        {
            "sample": h.index.astype(str),
            "CLDN4": h["CLDN4_log2tpm1"].astype(float).values,
            "CD274": cd274.reindex(h.index).astype(float).values,
            "IFN": ifn.reindex(h.index).astype(float).values,
            "IFNG": h["IFNG_log2tpm1"].astype(float).values,
            "IFN_sens": h["IFNG_Ayers6_meanz"].astype(float).values,
            "response": h["geo_response"].astype(str).values,
        }
    )
    out["cohort"] = "GSE218989"
    out.attrs["ifn_compact_genes"] = used
    # Sanity: harvested CLDN4 must match log2(TPM+1) of the extract.
    if "CLDN4" in logx.index:
        rec = logx.loc["CLDN4"].reindex(h.index)
        delta = float(np.nanmax(np.abs(rec.values - out["CLDN4"].values)))
        if delta > 1e-8:
            raise SystemExit(f"GSE218989 CLDN4 extract mismatch max|Δ|={delta}")
    return out


def load_gse126044() -> pd.DataFrame:
    counts = pd.read_csv(DATA / "GSE126044_counts.txt.gz", sep="\t", index_col=0)
    counts = counts.apply(pd.to_numeric, errors="coerce")
    lib = counts.sum(axis=0)
    logcpm = np.log2(1e6 * counts.div(lib, axis=1) + 1.0)
    ifn, used = zmean(logcpm, IFN_COMPACT)
    scores = pd.read_csv(HARVEST / "GSE126044_per_sample_scores.csv")
    scores = scores.set_index("sample")
    # CLDN4 from the already-published log2 CPM column (same transform).
    cldn4 = scores["CLDN4_log2cpm"].astype(float)
    rec = logcpm.loc["CLDN4"].reindex(cldn4.index)
    delta = float(np.nanmax(np.abs(rec.values - cldn4.values)))
    if delta > 1e-8:
        raise SystemExit(f"GSE126044 CLDN4 recompute mismatch max|Δ|={delta}")
    out = pd.DataFrame(
        {
            "sample": cldn4.index.astype(str),
            "CLDN4": cldn4.values,
            "CD274": logcpm.loc["CD274"].reindex(cldn4.index).astype(float).values,
            "IFN": ifn.reindex(cldn4.index).astype(float).values,
            "IFNG": logcpm.loc["IFNG"].reindex(cldn4.index).astype(float).values,
            "IFN_sens": np.nan,
            "response": scores["response"].astype(str).values,
        }
    )
    out["cohort"] = "GSE126044"
    out.attrs["ifn_compact_genes"] = used
    return out


def load_gse166449() -> pd.DataFrame:
    # Deposited "Raw gene TPM" matches the harvested CLDN4 column as-is.
    # Prior writeup: use deposited matrix; do not re-log.
    expr = pd.read_csv(DATA / "GSE166449_Raw_gene_TPM_matrix.txt.gz", sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    ifn, used = zmean(expr, IFN_COMPACT)
    samp = pd.read_csv(HARVEST / "GSE166449_sample_level_expression.csv").set_index("sample_id")
    rec = expr.loc["CLDN4"].reindex(samp.index)
    delta = float(np.nanmax(np.abs(rec.values - samp["CLDN4"].values)))
    if delta > 1e-8:
        raise SystemExit(f"GSE166449 CLDN4 extract mismatch max|Δ|={delta}")
    out = pd.DataFrame(
        {
            "sample": samp.index.astype(str),
            "CLDN4": samp["CLDN4"].astype(float).values,
            "CD274": expr.loc["CD274"].reindex(samp.index).astype(float).values,
            "IFN": ifn.reindex(samp.index).astype(float).values,
            "IFNG": expr.loc["IFNG"].reindex(samp.index).astype(float).values,
            "IFN_sens": np.nan,
            "response": samp["response"].astype(str).values,
        }
    )
    out["cohort"] = "GSE166449"
    out.attrs["ifn_compact_genes"] = used
    return out


def assign_quadrants(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d["CLDN4_q"] = quartiles(d["CLDN4"])
    d["IFN_q"] = quartiles(d["IFN"])
    d["CD274_q"] = quartiles(d["CD274"])
    d["IFNG_q"] = quartiles(d["IFNG"])
    d["cldn4_q4"] = d["CLDN4_q"].eq("Q4")
    d["ifn_q4"] = d["IFN_q"].eq("Q4")
    d["cd274_q4"] = d["CD274_q"].eq("Q4")
    d["ifn_q1"] = d["IFN_q"].eq("Q1")
    d["ifn_below_median"] = d["IFN"] <= d["IFN"].median()
    d["combo"] = d["cldn4_q4"] & (d["ifn_q4"] | d["cd274_q4"])
    d["ifn_low"] = d["cldn4_q4"] & d["ifn_q1"]
    d["combo_and_ifn_low"] = d["combo"] & d["ifn_low"]
    d["cldn4_q4_ifn_q4"] = d["cldn4_q4"] & d["ifn_q4"]
    d["cldn4_q4_cd274_q4"] = d["cldn4_q4"] & d["cd274_q4"]
    d["cldn4_q4_ifn_below_med"] = d["cldn4_q4"] & d["ifn_below_median"]
    # Exclusive labels among CLDN4 Q4 for the extra scatter.
    lab = np.full(len(d), "other", dtype=object)
    lab = np.where(d["cldn4_q4"], "CLDN4 Q4 mid-IFN", lab)
    lab = np.where(d["ifn_low"] & ~d["combo"], "CLDN4 Q4 ∩ IFN-low", lab)
    lab = np.where(d["combo"] & ~d["ifn_low"], "CLDN4 Q4 ∩ IFN/CD274 Q4", lab)
    lab = np.where(d["combo_and_ifn_low"], "overlap (IFN Q1 ∩ CD274 Q4)", lab)
    d["quadrant_label"] = lab
    return d


def cohort_counts(d: pd.DataFrame) -> dict:
    n = int(len(d))
    n_q4 = int(d["cldn4_q4"].sum())
    n_combo = int(d["combo"].sum())
    n_low = int(d["ifn_low"].sum())
    n_overlap = int(d["combo_and_ifn_low"].sum())
    n_ifn_q4 = int(d["cldn4_q4_ifn_q4"].sum())
    n_cd274_q4 = int(d["cldn4_q4_cd274_q4"].sum())
    n_below = int(d["cldn4_q4_ifn_below_med"].sum())
    n_mid = n_q4 - n_combo - n_low + n_overlap
    q_counts = d["CLDN4_q"].value_counts().to_dict()
    fish_ifn = fisher_q4_overlap(d["CLDN4"], d["IFN"])
    fish_cd274 = fisher_q4_overlap(d["CLDN4"], d["CD274"])
    rho_ifn = stats.spearmanr(d["CLDN4"], d["IFN"])
    rho_cd274 = stats.spearmanr(d["CLDN4"], d["CD274"])
    # Response among CLDN4 Q4, if labels exist (extra; not the claim).
    resp = d["response"].astype(str)
    has_resp = resp.notna() & ~resp.isin(["nan", "NaN", ""])
    extra_resp = {}
    if has_resp.any() and resp[has_resp].nunique() >= 2:
        sub = d.loc[d["cldn4_q4"] & has_resp]
        extra_resp = {
            "n_cldn4_q4_with_label": int(len(sub)),
            "response_counts_cldn4_q4": sub["response"].value_counts().to_dict(),
            "combo_response": d.loc[d["combo"], "response"].value_counts().to_dict(),
            "ifn_low_response": d.loc[d["ifn_low"], "response"].value_counts().to_dict(),
        }
    return {
        "cohort": d["cohort"].iloc[0],
        "n_matrix": n,
        "n_complete": n,
        "n_cldn4_q1": int(q_counts.get("Q1", 0)),
        "n_cldn4_q2": int(q_counts.get("Q2", 0)),
        "n_cldn4_q3": int(q_counts.get("Q3", 0)),
        "n_cldn4_q4": n_q4,
        "n_combo_cldn4_q4_and_ifn_or_cd274_q4": n_combo,
        "n_cldn4_q4_and_ifn_q4": n_ifn_q4,
        "n_cldn4_q4_and_cd274_q4": n_cd274_q4,
        "n_cldn4_q4_and_ifn_q1": n_low,
        "n_cldn4_q4_and_ifn_below_median": n_below,
        "n_overlap_combo_and_ifn_q1": n_overlap,
        "n_cldn4_q4_neither_combo_nor_ifn_q1": n_mid,
        "frac_cldn4_q4_combo": n_combo / n_q4 if n_q4 else np.nan,
        "frac_cldn4_q4_ifn_low": n_low / n_q4 if n_q4 else np.nan,
        "cldn4_vs_ifn_rho": float(rho_ifn.statistic),
        "cldn4_vs_ifn_p": float(rho_ifn.pvalue),
        "cldn4_vs_cd274_rho": float(rho_cd274.statistic),
        "cldn4_vs_cd274_p": float(rho_cd274.pvalue),
        "ifn_q4_overlap": fish_ifn,
        "cd274_q4_overlap": fish_cd274,
        "cldn4_q4_cut": q_cut_value(d["CLDN4"], 0.75),
        "ifn_q4_cut": q_cut_value(d["IFN"], 0.75),
        "cd274_q4_cut": q_cut_value(d["CD274"], 0.75),
        "ifn_q1_cut": q_cut_value(d["IFN"], 0.25),
        "extra_response": extra_resp,
        "tacstd2_gate": False,
        "dual_high": False,
    }


def plot_extra_scatter(d: pd.DataFrame, cohort: str) -> None:
    colors = {
        "other": "#bdbdbd",
        "CLDN4 Q4 mid-IFN": "#9ecae1",
        "CLDN4 Q4 ∩ IFN-low": "#3182bd",
        "CLDN4 Q4 ∩ IFN/CD274 Q4": "#e6550d",
        "overlap (IFN Q1 ∩ CD274 Q4)": "#756bb1",
    }
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.4))
    for ax, ycol, ylab, yq in (
        (axes[0], "IFN", "IFN-γ compact (mean-z)", "IFN_q"),
        (axes[1], "CD274", "CD274", "CD274_q"),
    ):
        for lab, col in colors.items():
            m = d["quadrant_label"] == lab
            if not m.any():
                continue
            ax.scatter(
                d.loc[m, "CLDN4"],
                d.loc[m, ycol],
                s=28 if d["cohort"].iloc[0] in {"GSE126044", "GSE166449"} else 16,
                c=col,
                alpha=0.85,
                edgecolors="none",
                label=lab if ax is axes[0] else None,
            )
        xq4 = d.loc[d["CLDN4_q"].eq("Q4"), "CLDN4"].min()
        ax.axvline(xq4, color="#636363", ls="--", lw=0.8)
        yq4 = d.loc[d[yq].eq("Q4"), ycol].min()
        ax.axhline(yq4, color="#636363", ls="--", lw=0.8)
        if ycol == "IFN":
            yq1 = d.loc[d["IFN_q"].eq("Q1"), "IFN"].max()
            ax.axhline(yq1, color="#9ecae1", ls=":", lw=0.8)
        ax.set_xlabel("CLDN4")
        ax.set_ylabel(ylab)
        ax.set_title(f"{cohort}  n={len(d)}")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(FIG / f"{cohort}_extra_scatter.png", dpi=160)
    fig.savefig(FIG / f"{cohort}_extra_scatter.pdf")
    plt.close(fig)


def plot_count_bars(rows: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    names = [r["cohort"] for r in rows]
    x = np.arange(len(names))
    w = 0.36
    combo = [r["n_combo_cldn4_q4_and_ifn_or_cd274_q4"] for r in rows]
    low = [r["n_cldn4_q4_and_ifn_q1"] for r in rows]
    ax.bar(x - w / 2, combo, w, color="#e6550d", label="CLDN4 Q4 ∩ (IFN-γ Q4 ∪ CD274 Q4)")
    ax.bar(x + w / 2, low, w, color="#3182bd", label="CLDN4 Q4 ∩ IFN Q1")
    for i, r in enumerate(rows):
        ax.text(i, max(combo[i], low[i]) + 0.4, f"Q4 n={r['n_cldn4_q4']}", ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("patients (honest count)")
    ax.set_title("CLDN4-only ADC+ICI quadrant  (no dual-high)")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "quadrant_counts.png", dpi=160)
    fig.savefig(FIG / "quadrant_counts.pdf")
    plt.close(fig)


def main() -> None:
    FIG.mkdir(exist_ok=True)
    TAB.mkdir(exist_ok=True)

    loaders = {
        "GSE285029": load_gse285029,
        "GSE218989": load_gse218989,
        "GSE126044": load_gse126044,
        "GSE166449": load_gse166449,
    }
    assigned = []
    rows = []
    for name, fn in loaders.items():
        raw = fn()
        d = assign_quadrants(raw)
        assigned.append(d)
        rows.append(cohort_counts(d))
        plot_extra_scatter(d, name)
        d.to_csv(TAB / f"{name}_sample_quadrants.tsv", sep="\t", index=False)

    plot_count_bars(rows)

    # Flat quadrant table (the PR deliverable).
    flat = []
    for r in rows:
        flat.append(
            {
                "cohort": r["cohort"],
                "n_matrix": r["n_matrix"],
                "n_CLDN4_Q4": r["n_cldn4_q4"],
                "n_CLDN4_Q4_and_IFNG_Q4": r["n_cldn4_q4_and_ifn_q4"],
                "n_CLDN4_Q4_and_CD274_Q4": r["n_cldn4_q4_and_cd274_q4"],
                "n_CLDN4_Q4_and_IFNG_or_CD274_Q4": r["n_combo_cldn4_q4_and_ifn_or_cd274_q4"],
                "n_CLDN4_Q4_and_IFN_Q1": r["n_cldn4_q4_and_ifn_q1"],
                "n_CLDN4_Q4_and_IFN_below_median": r["n_cldn4_q4_and_ifn_below_median"],
                "n_overlap_combo_and_IFN_Q1": r["n_overlap_combo_and_ifn_q1"],
                "n_CLDN4_Q4_neither": r["n_cldn4_q4_neither_combo_nor_ifn_q1"],
                "frac_Q4_combo": r["frac_cldn4_q4_combo"],
                "frac_Q4_IFN_low": r["frac_cldn4_q4_ifn_low"],
                "CLDN4_vs_IFN_rho": r["cldn4_vs_ifn_rho"],
                "CLDN4_vs_IFN_p": r["cldn4_vs_ifn_p"],
                "CLDN4_vs_CD274_rho": r["cldn4_vs_cd274_rho"],
                "CLDN4_vs_CD274_p": r["cldn4_vs_cd274_p"],
                "IFN_Q4_overlap_OR": r["ifn_q4_overlap"]["odds_ratio"],
                "IFN_Q4_overlap_p": r["ifn_q4_overlap"]["fisher_p"],
                "CD274_Q4_overlap_OR": r["cd274_q4_overlap"]["odds_ratio"],
                "CD274_Q4_overlap_p": r["cd274_q4_overlap"]["fisher_p"],
                "TACSTD2_gate": False,
                "dual_high": False,
            }
        )
    flat_df = pd.DataFrame(flat)
    flat_df.to_csv(TAB / "quadrant_table.tsv", sep="\t", index=False)

    all_s = pd.concat(assigned, ignore_index=True)
    all_s.to_csv(TAB / "all_sample_quadrants.tsv", sep="\t", index=False)

    summary = {
        "generated": "CLDN4-only ADC+ICI quadrant",
        "dual_high": False,
        "tacstd2_gate": False,
        "ifn_primary": "IFN-compact mean-z (A11 list)",
        "ifn_low": "IFN Q1 (equal-count quartile)",
        "combo": "CLDN4 Q4 AND (IFN Q4 OR CD274 Q4)",
        "quartile_rule": "rank(method='first') then qcut(4)",
        "cohorts": rows,
        "quadrant_table": flat,
    }
    (TAB / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(flat_df.to_string(index=False))


if __name__ == "__main__":
    main()
