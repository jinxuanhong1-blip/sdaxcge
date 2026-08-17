#!/usr/bin/env python3
"""ADDITIVE CLDN4-only pairwise merge: GSE207422 + GSE131907.

1) Patient-level malignant CLDN4 vs T/NK (honest n, Spearman + Q4 vs Q1).
2) CellChat-style outgoing CLDN4-high → T/NK combo LR table.

No dual-high TACSTD2×CLDN4. Prior single-cohort CellChat folders are taken
as given and are not re-run. Matrices are not re-downloaded.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from lib_stats import random_effects_dl, spearman, stouffer

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TNK = DATA / "tnk"
CCC = DATA / "cellchat"
TABLES = HERE / "results"
FIGS = HERE / "figures"
NEG = "#7a2d0b"
POS = "#4a4a4a"
SHARE = "#1f4e79"
DISC = "#8a4b08"

# Pre-specified outgoing pairs to always show (barrier / checkpoint / recruit).
KEY_PAIRS = [
    "CDH1_ITGAE_ITGB7",
    "CDH1_KLRG1",
    "JAM1_ITGAL_ITGB2",
    "NECTIN2_TIGIT",
    "NECTIN2_CD226",
    "PVR_TIGIT",
    "CD274_PDCD1",
    "HLA-E_CD8A",
    "HLA-E_CD8B",
    "HLA-G_CD8A",
    "HLA-G_CD8B",
    "LGALS9_CD45",
    "LGALS9_CD44",
    "CXCL16_CXCR6",
    "CXCL9_CXCR3",
    "CXCL10_CXCR3",
    "CCL5_CCR5",
    "SPP1_CD44",
    "MIF_CD74_CD44",
    "MIF_CD74_CXCR4",
]


def q4_vs_q1(cldn4, immune) -> dict | None:
    s = pd.DataFrame({"c": np.asarray(cldn4, dtype=float), "i": np.asarray(immune, dtype=float)})
    s = s[np.isfinite(s["c"]) & np.isfinite(s["i"])].copy()
    n = int(len(s))
    if n < 6:
        return None
    ranks = s["c"].rank(method="average")
    try:
        qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    except ValueError:
        return None
    if qs.nunique() < 4:
        return None
    q1 = s.loc[qs == "Q1", "i"]
    q4 = s.loc[qs == "Q4", "i"]
    n1, n4 = int(len(q1)), int(len(q4))
    if n1 < 2 or n4 < 2:
        return None
    u, p = stats.mannwhitneyu(q4.values, q1.values, alternative="two-sided")
    r_rb = (2.0 * float(u)) / (n4 * n1) - 1.0
    return {
        "n": n,
        "n_q1": n1,
        "n_q4": n4,
        "n_compared": n1 + n4,
        "median_q1": float(q1.median()),
        "median_q4": float(q4.median()),
        "delta_median": float(q4.median() - q1.median()),
        "mwu_u": float(u),
        "p": float(p),
        "r_rb": float(r_rb),
        "thin": n < 8 or n1 < 3 or n4 < 3,
        "q1": q1.to_numpy(),
        "q4": q4.to_numpy(),
        "cldn4": s["c"].to_numpy(),
        "immune": s["i"].to_numpy(),
        "quartile": qs.astype(str).to_numpy(),
    }


def atomic(cohort, malig_def, immune_def, score, cldn4, immune, **note) -> dict:
    rho, p, n = spearman(cldn4, immune)
    rec = {
        "cohort": cohort,
        "malig_def": malig_def,
        "immune_def": immune_def,
        "score": score,
        "n": n,
        "rho": rho,
        "p": p,
    }
    rec.update(note)
    q = q4_vs_q1(cldn4, immune)
    if q is None:
        rec.update(
            {
                "n_q1": np.nan,
                "n_q4": np.nan,
                "n_compared": np.nan,
                "median_q1": np.nan,
                "median_q4": np.nan,
                "delta_median": np.nan,
                "r_rb": np.nan,
                "p_q4q1": np.nan,
                "thin_q4q1": True,
                "poolable_q4q1": False,
            }
        )
        rec["_q"] = None
    else:
        rec.update(
            {
                "n_q1": q["n_q1"],
                "n_q4": q["n_q4"],
                "n_compared": q["n_compared"],
                "median_q1": q["median_q1"],
                "median_q4": q["median_q4"],
                "delta_median": q["delta_median"],
                "r_rb": q["r_rb"],
                "p_q4q1": q["p"],
                "thin_q4q1": q["thin"],
                "poolable_q4q1": (not q["thin"]) and abs(q["r_rb"]) < 0.999,
            }
        )
        rec["_q"] = q
    rec["_x"] = np.asarray(cldn4, dtype=float)
    rec["_y"] = np.asarray(immune, dtype=float)
    return rec


def pool_spearman(rhos, ns, ps):
    rhos = [float(r) for r in rhos]
    ns = [int(n) for n in ns]
    ps = [float(p) for p in ps]
    re = random_effects_dl(rhos, ns) if len(rhos) > 1 else {
        "k": 1, "n_patients_total": ns[0], "pooled_rho": rhos[0], "p": ps[0], "I2": 0.0,
        "ci95_rho": [float("nan"), float("nan")],
    }
    return re, stouffer(rhos, ps, ns)


def pool_q4q1(rr, ns, ps):
    rr = [float(r) for r in rr]
    ns = [int(n) for n in ns]
    ps = [float(p) for p in ps]
    re = random_effects_dl(rr, ns) if len(rr) > 1 else {
        "k": 1, "n_patients_total": ns[0], "pooled_rho": rr[0], "p": ps[0], "I2": 0.0,
        "ci95_rho": [float("nan"), float("nan")],
    }
    return re, stouffer(rr, ps, ns)


def load_pair_atomics() -> list[dict]:
    a3 = pd.read_csv(TNK / "GSE207422_drmref_patients.tsv", sep="\t")
    d = pd.read_csv(TNK / "GSE131907_samples.tsv", sep="\t")
    tumor = d[d["origin"].isin(["tLung", "tL/B", "mLN", "PE", "mBrain"])]
    mal = tumor[tumor["n_malignant"] >= 20].copy()
    tlung_epi = tumor[(tumor["origin"] == "tLung") & (tumor["n_epithelial"] >= 20)]

    rows = [
        atomic(
            "GSE207422",
            "author_DRMref",
            "tnk",
            "mean",
            a3["malig_CLDN4_mean"],
            a3["frac_tnk"],
            unit="patient",
            note="locked A3 12-patient DRMref table; CLDN4 only; TACSTD2 slide not re-cut",
            n_malignant_cells=int(a3["n_malignant"].sum()),
            n_tnk_cells=int(a3["n_tnk"].sum()),
            origins="post-treatment resection (12/12)",
        ),
        atomic(
            "GSE131907",
            "author_malig",
            "tnk",
            "mean",
            mal["mal_CLDN4_mean"],
            mal["frac_tnk"],
            unit="sample",
            note="author Malignant cells, n_mal>=20; sample-level; tLung tS* not in this column",
            n_malignant_cells=int(mal["n_malignant"].sum()),
            n_tnk_cells=int(mal["n_tnk"].sum()),
            origins=",".join(sorted(mal["origin"].unique())),
        ),
        atomic(
            "GSE131907",
            "author_malig",
            "tnk",
            "pct",
            mal["mal_CLDN4_pct"],
            mal["frac_tnk"],
            unit="sample",
            note="sensitivity: malignant CLDN4 %pos vs T/NK fraction",
            n_malignant_cells=int(mal["n_malignant"].sum()),
            n_tnk_cells=int(mal["n_tnk"].sum()),
            origins=",".join(sorted(mal["origin"].unique())),
        ),
    ]
    if len(tlung_epi) >= 4:
        rows.append(
            atomic(
                "GSE131907_tLung",
                "author_epi",
                "tnk",
                "mean",
                tlung_epi["epi_CLDN4_mean"],
                tlung_epi["frac_tnk"],
                unit="sample",
                note="sensitivity only: tLung epithelium (tS1–tS3 live here; not primary)",
                n_malignant_cells=int(tlung_epi["n_epithelial"].sum()),
                n_tnk_cells=int(tlung_epi["n_tnk"].sum()),
                origins="tLung",
            )
        )
    return rows


def _bool_series(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    return s.astype(str).str.lower().isin(["true", "1", "yes"])


def load_outgoing(path: Path, cohort: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df[df["direction"].astype(str).str.lower() == "outgoing"].copy()
    for col in ("sig_high", "sig_low", "sig_either", "sig_diff"):
        if col in df.columns:
            df[col] = _bool_series(df[col])
    df["cohort"] = cohort
    df["detected"] = (df["prob_high"] > 0) | (df["prob_low"] > 0)
    df["sig"] = df.get("sig_diff", df["sig_high"] | df["sig_low"])
    if "sig_diff" in df.columns:
        df["sig"] = df["sig_diff"] & df["detected"]
    else:
        df["sig"] = (df["sig_high"] | df["sig_low"]) & df["detected"]
    return df


def combo_lr(out_a: pd.DataFrame, out_b: pd.DataFrame) -> pd.DataFrame:
    """Join outgoing contrasts on interaction_name. A=207422, B=131907."""
    keep = [
        "interaction_name",
        "pathway_name",
        "annotation",
        "ligand",
        "receptor",
        "ligand_genes",
        "receptor_genes",
        "ligand_class",
        "prob_high",
        "prob_low",
        "delta_prob",
        "pval_high",
        "pval_low",
        "sig_high",
        "sig_low",
        "sig",
        "detected",
    ]
    a = out_a[keep].copy()
    b = out_b[keep].copy()
    m = a.merge(b, on="interaction_name", how="outer", suffixes=("_207422", "_131907"))
    for col in (
        "pathway_name",
        "annotation",
        "ligand",
        "receptor",
        "ligand_genes",
        "receptor_genes",
        "ligand_class",
    ):
        m[col] = m[f"{col}_207422"].fillna(m[f"{col}_131907"])
        m.drop(columns=[f"{col}_207422", f"{col}_131907"], inplace=True)

    def _arm(sig, delta) -> str:
        if pd.isna(sig) or not bool(sig):
            return "ns"
        if float(delta) > 0:
            return "high"
        if float(delta) < 0:
            return "low"
        return "ns"

    m["arm_207422"] = [
        _arm(s, d) for s, d in zip(m["sig_207422"], m["delta_prob_207422"])
    ]
    m["arm_131907"] = [
        _arm(s, d) for s, d in zip(m["sig_131907"], m["delta_prob_131907"])
    ]

    def _conc(r) -> str:
        a_arm, b_arm = r["arm_207422"], r["arm_131907"]
        a_det = bool(r["detected_207422"]) if pd.notna(r["detected_207422"]) else False
        b_det = bool(r["detected_131907"]) if pd.notna(r["detected_131907"]) else False
        if a_arm != "ns" and b_arm != "ns":
            return "both_sig_same" if a_arm == b_arm else "both_sig_discordant"
        if a_arm != "ns" and b_arm == "ns":
            return "only_207422" if b_det else "only_207422_undetected_131907"
        if b_arm != "ns" and a_arm == "ns":
            return "only_131907" if a_det else "only_131907_undetected_207422"
        if a_det and b_det:
            return "both_detected_ns"
        return "not_differential"

    m["concordance"] = m.apply(_conc, axis=1)
    m["sig_either"] = (m["arm_207422"] != "ns") | (m["arm_131907"] != "ns")
    m["mean_delta"] = m[["delta_prob_207422", "delta_prob_131907"]].mean(axis=1)
    m["abs_mean_delta"] = m["mean_delta"].abs()
    m["key_pair"] = m["interaction_name"].isin(KEY_PAIRS)
    # Honest CellChat n (from the source FINDINGs; cell-pooled, not patient mixed model).
    m["n_mal_high_207422"] = 4203
    m["n_mal_low_207422"] = 4204
    m["n_tnk_207422"] = 33760
    m["n_patients_207422"] = 12
    m["n_mal_high_131907"] = 10379
    m["n_mal_low_131907"] = 10379
    m["n_tnk_131907"] = 34741
    m["n_patients_131907"] = 32
    m["n_patients_pair"] = 12 + 32
    return m.sort_values(
        ["sig_either", "key_pair", "abs_mean_delta"],
        ascending=[False, False, False],
    )


def forest(rows: list[dict], title: str, path: Path, xlabel: str, effect: str, pcol: str) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 1.15 + 0.48 * max(len(rows), 1)))
    y = np.arange(len(rows))
    for i, r in enumerate(rows):
        eff = float(r[effect])
        p = float(r[pcol])
        color = NEG if eff < 0 else POS
        if r.get("pooled"):
            color = SHARE if eff < 0 else POS
            ax.plot(eff, i, "D", color=color, ms=8, zorder=3)
        else:
            ax.plot(eff, i, "o", color=color, ms=7, zorder=3)
        ax.text(
            0.98,
            i,
            f"{eff:+.3f} p={p:.3g}",
            va="center",
            ha="right",
            fontsize=7,
            family="monospace",
            transform=ax.get_yaxis_transform(),
        )
    ax.axvline(0, color="#888", lw=0.8, ls="--")
    ax.set_xlim(-1.05, 1.05)
    ax.set_yticks(y)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def scatter(path: Path, x, y, xlabel, ylabel, title) -> None:
    fig, ax = plt.subplots(figsize=(4.5, 3.7))
    ax.scatter(list(x), list(y), c=NEG, s=38, zorder=3)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def box_q4q1(path: Path, q: dict, title: str, ylabel: str) -> None:
    q1, q4 = q["q1"], q["q4"]
    fig, ax = plt.subplots(figsize=(4.3, 3.8))
    bp = ax.boxplot(
        [q1, q4],
        tick_labels=[f"Q1\nn={len(q1)}", f"Q4\nn={len(q4)}"],
        patch_artist=True,
        widths=0.55,
    )
    for patch, color in zip(bp["boxes"], ["#6a8aaa", NEG]):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for i, vals in enumerate((q1, q4), start=1):
        ax.scatter(i + rng.uniform(-0.08, 0.08, size=len(vals)), vals, c="black", s=16, zorder=3)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_honest_n(path: Path, a207, a131) -> None:
    labels = [
        "GSE207422 patients\n(DRMref malig)",
        "GSE131907 samples\n(author malig ≥20)",
        "Pair N (Spearman)",
        "Q4+Q1 compared",
    ]
    vals = [
        a207["n"],
        a131["n"],
        a207["n"] + a131["n"],
        int(a207["n_compared"] + a131["n_compared"]),
    ]
    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    bars = ax.bar(labels, vals, color=[NEG, NEG, SHARE, "#6a8aaa"])
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.4, str(v), ha="center", fontsize=9)
    ax.set_ylabel("Honest n")
    ax.set_title("Pair GSE207422 + GSE131907 · patient/sample n (not cell n)", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_lr_concordance(path: Path, lr: pd.DataFrame) -> None:
    sig = lr[lr["sig_either"]].copy()
    counts = (
        sig["concordance"]
        .value_counts()
        .reindex(
            [
                "both_sig_same",
                "both_sig_discordant",
                "only_207422",
                "only_207422_undetected_131907",
                "only_131907",
                "only_131907_undetected_207422",
            ]
        )
        .fillna(0)
        .astype(int)
    )
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    colors = [SHARE, DISC, NEG, "#b07a5a", "#3d6b8a", "#7aa0b8"]
    ax.barh(counts.index.astype(str), counts.values, color=colors)
    for i, v in enumerate(counts.values):
        ax.text(v + 0.3, i, str(int(v)), va="center", fontsize=8)
    ax.set_xlabel("Outgoing Mal → T/NK pairs (significant in ≥1 cohort)")
    ax.set_title("Combo LR concordance · CLDN4-high vs low", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_shared_delta(path: Path, lr: pd.DataFrame) -> None:
    both = lr[lr["concordance"].isin(["both_sig_same", "both_sig_discordant"])].copy()
    if both.empty:
        return
    both = both.sort_values("abs_mean_delta", ascending=True).tail(16)
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    y = np.arange(len(both))
    ax.scatter(both["delta_prob_207422"], y, c=NEG, s=36, label="GSE207422 ΔP", zorder=3)
    ax.scatter(both["delta_prob_131907"], y, c=SHARE, s=36, label="GSE131907 ΔP", zorder=3)
    for i, r in enumerate(both.itertuples()):
        ax.plot(
            [r.delta_prob_207422, r.delta_prob_131907],
            [i, i],
            color="#bbb",
            lw=0.8,
            zorder=1,
        )
    ax.axvline(0, color="#888", lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(both["interaction_name"], fontsize=7)
    ax.set_xlabel("ΔP (CLDN4-high − CLDN4-low), outgoing Mal → T/NK")
    ax.set_title("Shared significant outgoing pairs", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_extra_ligand_table(path: Path, lr: pd.DataFrame) -> None:
    show = lr[lr["key_pair"] | (lr["concordance"] == "both_sig_same")].copy()
    show = show.drop_duplicates("interaction_name")
    show = show.sort_values("abs_mean_delta", ascending=True)
    if show.empty:
        return
    fig, ax = plt.subplots(figsize=(9.2, 0.42 * len(show) + 1.6))
    y = np.arange(len(show))
    ax.scatter(show["delta_prob_207422"], y - 0.12, c=NEG, s=28, label="GSE207422")
    ax.scatter(show["delta_prob_131907"], y + 0.12, c=SHARE, s=28, label="GSE131907")
    ax.axvline(0, color="#888", lw=0.8, ls="--")
    labels = []
    for r in show.itertuples():
        tags = []
        if r.arm_207422 != "ns":
            tags.append(f"422:{r.arm_207422}")
        else:
            tags.append("422:ns")
        if r.arm_131907 != "ns":
            tags.append(f"907:{r.arm_131907}")
        else:
            tags.append("907:ns")
        labels.append(f"{r.interaction_name}  ({', '.join(tags)})")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Outgoing ΔP (high − low)")
    ax.set_title("Extra: key + shared outgoing CLDN4-high → T/NK pairs", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def _fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return f"{float(x):+.{nd}f}" if nd else f"{x}"


def _fp(p):
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "NA"
    return f"{float(p):.3g}"


def write_finding(
    a207: dict,
    a131: dict,
    a131_pct: dict,
    a_tlung: dict | None,
    re_s: dict,
    st_s: dict,
    re_q: dict,
    st_q: dict,
    re_s_pct: dict,
    re_q_pct: dict,
    lr: pd.DataFrame,
    path: Path,
) -> None:
    both_same = lr[lr["concordance"] == "both_sig_same"]
    both_disc = lr[lr["concordance"] == "both_sig_discordant"]
    only_a = lr[lr["concordance"].str.startswith("only_207422")]
    only_b = lr[lr["concordance"].str.startswith("only_131907")]
    shared_up = both_same[both_same["arm_207422"] == "high"].sort_values(
        "abs_mean_delta", ascending=False
    )
    shared_dn = both_same[both_same["arm_207422"] == "low"].sort_values(
        "abs_mean_delta", ascending=False
    )
    key = lr[lr["key_pair"]].copy()

    def _row(r):
        return (
            f"| {r.interaction_name} | {r.pathway_name} | {r.ligand_class} | "
            f"{r.concordance} | {_fmt(r.delta_prob_207422)} | {_fmt(r.delta_prob_131907)} | "
            f"{r.arm_207422} | {r.arm_131907} |"
        )

    n_pair = a207["n"] + a131["n"]
    n_cmp = int(a207["n_compared"] + a131["n_compared"])
    ci = re_s.get("ci95_rho") or [float("nan"), float("nan")]
    ci_q = re_q.get("ci95_rho") or [float("nan"), float("nan")]

    lines = [
        "# FINDING — Pair GSE207422 + GSE131907, CLDN4-only (patient T/NK then CellChat outgoing)",
        "",
        "**Additive. CLDN4 only. No dual-high.** GSE207422 is included. Prior single-cohort",
        "CellChat folders (`methods/scrna_cellchat_cldn4`, `methods/gse131907_cellchat_cldn4`)",
        "and the malignant Q4 T/NK extract (`methods/cldn4_malig_q4_tnk`) are taken as given",
        "and were **not** re-run. Matrices were not re-downloaded.",
        "",
        "Primary tables: [`results/combo_rho.tsv`](results/combo_rho.tsv) and",
        "[`results/combo_lr_table.tsv`](results/combo_lr_table.tsv)",
        "(also `ligand_table.tsv`). Extra figures under [`figures/`](figures/).",
        "",
        "## 1. Patient-level malignant CLDN4 vs T/NK (honest n)",
        "",
        "Unit is the **patient** on GSE207422 and the **sample** on GSE131907 (Kim atlas",
        "has one tumor-origin sample per patient in the author-malignant extract).",
        "Quartiles are **within cohort**, then Q4 vs Q1 T/NK is pooled as rank-biserial *r*",
        "on Fisher-z (DerSimonian–Laird). Spearman is the same RE on Fisher-z(ρ).",
        "p-values are descriptive.",
        "",
        "| Item | Public? | n | Note |",
        "|---|---|---:|---|",
        f"| GSE207422 patients (DRMref malignant) | yes | **{a207['n']}** | locked A3 table; post-treatment resections |",
        f"| GSE207422 Q4 vs Q1 tails | yes | **{int(a207['n_q1'])} vs {int(a207['n_q4'])}** | n=12 so each tail is 3; poolable but thin |",
        f"| GSE131907 author-malignant samples | yes | **{a131['n']}** | `n_malignant ≥ 20`; origins {a131['origins']} |",
        f"| GSE131907 Q4 vs Q1 tails | yes | **{int(a131['n_q1'])} vs {int(a131['n_q4'])}** | sample-level |",
        f"| **Pair Spearman N** | yes | **{n_pair}** | 12 patients + {a131['n']} samples |",
        f"| **Pair Q4+Q1 compared** | yes | **{n_cmp}** | {int(a207['n_q1']+a131['n_q1'])} Q1 + {int(a207['n_q4']+a131['n_q4'])} Q4 |",
        "| tLung tS1–tS3 in the primary T/NK row | **no** | 0 | primary tLung epithelium is not `Malignant cells` |",
        "| Dual-high TACSTD2×CLDN4 | **no** | 0 | CLDN4-only |",
        "| Cell-pooled CellChat n (below) | yes | see §2 | not a 44-patient mixed model |",
        "",
        "Do not write n=44 (Kim series patients) or n=15 (Hu series samples) for the pair T/NK test.",
        f"The computable pair n is **{n_pair}** for Spearman and **{n_cmp}** compared tails for Q4 vs Q1.",
        "",
        "### Singles",
        "",
        "| cohort | malig | score | unit | n | Spearman ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) |",
        "|---|---|---|---|---:|---|---|",
        f"| GSE207422 | author_DRMref | mean | patient | {a207['n']} | {_fmt(a207['rho'])} ({_fp(a207['p'])}) | {_fmt(a207['r_rb'])} ({_fp(a207['p_q4q1'])}; {int(a207['n_q1'])}/{int(a207['n_q4'])}) |",
        f"| GSE131907 | author_malig | mean | sample | {a131['n']} | {_fmt(a131['rho'])} ({_fp(a131['p'])}) | {_fmt(a131['r_rb'])} ({_fp(a131['p_q4q1'])}; {int(a131['n_q1'])}/{int(a131['n_q4'])}) |",
        f"| GSE131907 (sensitivity) | author_malig | pct | sample | {a131_pct['n']} | {_fmt(a131_pct['rho'])} ({_fp(a131_pct['p'])}) | {_fmt(a131_pct['r_rb'])} ({_fp(a131_pct['p_q4q1'])}; {int(a131_pct['n_q1'])}/{int(a131_pct['n_q4'])}) |",
    ]
    if a_tlung:
        lines.append(
            f"| GSE131907 tLung (not primary) | author_epi | mean | sample | {a_tlung['n']} | {_fmt(a_tlung['rho'])} ({_fp(a_tlung['p'])}) | {_fmt(a_tlung.get('r_rb'))} ({_fp(a_tlung.get('p_q4q1'))}; {a_tlung.get('n_q1')}/{a_tlung.get('n_q4')}) |"
        )
    lines += [
        "",
        "### Combo (primary = mean / mean)",
        "",
        "| analysis | k | N | effect | p | I² | 95% CI | Stouffer p |",
        "|---|---:|---:|---|---|---:|---|---|",
        f"| Spearman ρ | 2 | {n_pair} | **{_fmt(re_s.get('pooled_rho'))}** | {_fp(re_s.get('p'))} | {re_s.get('I2', float('nan')):.0f}% | [{_fmt(ci[0])}, {_fmt(ci[1])}] | {_fp(st_s.get('p'))} |",
        f"| Q4 vs Q1 r | 2 | {n_cmp} | **{_fmt(re_q.get('pooled_rho'))}** | {_fp(re_q.get('p'))} | {re_q.get('I2', float('nan')):.0f}% | [{_fmt(ci_q[0])}, {_fmt(ci_q[1])}] | {_fp(st_q.get('p'))} |",
        f"| Spearman ρ (131907 %pos sensitivity) | 2 | {a207['n']+a131_pct['n']} | {_fmt(re_s_pct.get('pooled_rho'))} | {_fp(re_s_pct.get('p'))} | {re_s_pct.get('I2', float('nan')):.0f}% | — | — |",
        f"| Q4 vs Q1 r (131907 %pos sensitivity) | 2 | {int(a207['n_compared']+a131_pct['n_compared'])} | {_fmt(re_q_pct.get('pooled_rho'))} | {_fp(re_q_pct.get('p'))} | {re_q_pct.get('I2', float('nan')):.0f}% | — | — |",
        "",
        "GSE207422 alone is near-null (ρ = −0.09, n=12). GSE131907 author-malignant mean is",
        "negative and larger. The pair stays **CLDN4-negative vs T/NK** with I² = 0, but the",
        "Q4 vs Q1 p is above 0.05 on the primary mean/mean cut. The %pos sensitivity on",
        "GSE131907 is the stronger single (see table). That is the honest pair, not a hidden n.",
        "",
        "## 2. CellChat-style outgoing CLDN4-high → T/NK (combo LR)",
        "",
        "CellChat R was not run. Each cohort already has Jin et al. 2021 Hill probability on",
        "CellChatDB v2 protein pairs + 100 permutations of CLDN4-high/low among malignant",
        "cells (`expr_prop ≥ 0.10`, p < 0.05; smallest p = 1/101 = 0.0099).",
        "",
        "| Cohort | Kept split | Mal high / low | T/NK | Patients/samples | Detected outgoing | Sig outgoing |",
        "|---|---|---:|---:|---:|---:|---:|",
        f"| GSE207422 | median_post (epithelial proxy) | 4,203 / 4,204 | 33,760 | 12 | {int((lr['detected_207422']==True).sum())} | {int((lr['arm_207422']!='ns').sum())} |",
        f"| GSE131907 | tumor tertile (author malig + tS*) | 10,379 / 10,379 | 34,741 | 32 | {int((lr['detected_131907']==True).sum())} | {int((lr['arm_131907']!='ns').sum())} |",
        "",
        "Means are **cell-pooled**. BD_immune07 is ~61% of GSE207422 post epithelium;",
        "EBUS_28 is ~27% of GSE131907 Mal_high. That is not a 12+32 patient mixed model.",
        "",
        "### Concordance (outgoing, significant in ≥1 cohort)",
        "",
        f"| Concordance | n pairs |",
        f"|---|---:|",
        f"| both significant, same direction | **{len(both_same)}** |",
        f"| both significant, discordant | {len(both_disc)} |",
        f"| only GSE207422 | {len(only_a)} |",
        f"| only GSE131907 | {len(only_b)} |",
        "",
        "### Shared same-direction outgoing (higher in CLDN4-high)",
        "",
        "| Pair | Pathway | Class | Concordance | ΔP 207422 | ΔP 131907 | arm 422 | arm 907 |",
        "|---|---|---|---|---:|---:|---|---|",
    ]
    for r in shared_up.itertuples():
        lines.append(_row(r))
    if shared_up.empty:
        lines.append("| — | — | — | — | — | — | — | — |")
    lines += [
        "",
        "### Shared same-direction outgoing (higher in CLDN4-low)",
        "",
        "| Pair | Pathway | Class | Concordance | ΔP 207422 | ΔP 131907 | arm 422 | arm 907 |",
        "|---|---|---|---|---:|---:|---|---|",
    ]
    for r in shared_dn.itertuples():
        lines.append(_row(r))
    if shared_dn.empty:
        lines.append("| — | — | — | — | — | — | — | — |")
    lines += [
        "",
        "### Pre-specified key pairs (always shown, even if undetected)",
        "",
        "| Pair | Pathway | Class | Concordance | ΔP 207422 | ΔP 131907 | arm 422 | arm 907 |",
        "|---|---|---|---|---:|---:|---|---|",
    ]
    for r in key.sort_values("interaction_name").itertuples():
        lines.append(_row(r))
    lines += [
        "",
        "Shared barrier / inhibitory outgoing that stay **higher in CLDN4-high** on both",
        "cohorts (when detected): CDH1–ITGAE/ITGB7, JAM1–ITGAL/ITGB2, HLA-E–CD8A/B.",
        "CXCL16–CXCR6 is also **higher** in CLDN4-high on both kept splits (recruit-up,",
        "not recruit-down). NECTIN2–TIGIT and CD274–PDCD1 are GSE207422-only (not detected",
        "on GSE131907). Do not invent those edges on the Kim atlas.",
        "",
        "Full combo table: [`results/combo_lr_table.tsv`](results/combo_lr_table.tsv).",
        "Significant-in-either outgoing slice: [`results/ligand_table.tsv`](results/ligand_table.tsv).",
        "",
        "## Extra figures",
        "",
        "- [`figures/fig_honest_n.png`](figures/fig_honest_n.png) — pair n vs compared tails",
        "- [`figures/fig_combo_rho_forest.png`](figures/fig_combo_rho_forest.png)",
        "- [`figures/fig_combo_q4q1_forest.png`](figures/fig_combo_q4q1_forest.png)",
        "- [`figures/fig_scatter_GSE207422.png`](figures/fig_scatter_GSE207422.png)",
        "- [`figures/fig_scatter_GSE131907.png`](figures/fig_scatter_GSE131907.png)",
        "- [`figures/fig_q4q1_box_GSE207422.png`](figures/fig_q4q1_box_GSE207422.png)",
        "- [`figures/fig_q4q1_box_GSE131907.png`](figures/fig_q4q1_box_GSE131907.png)",
        "- [`figures/fig_extra_lr_concordance.png`](figures/fig_extra_lr_concordance.png)",
        "- [`figures/fig_extra_shared_outgoing.png`](figures/fig_extra_shared_outgoing.png)",
        "- [`figures/fig_extra_ligand_table.png`](figures/fig_extra_ligand_table.png)",
        "",
        "## What is not claimed",
        "",
        "- This is **not** dual-high TACSTD2×CLDN4 and **not** a 6-unit / 4-unit search.",
        "- GSE207422 CellChat uses marker epithelium (CopyKAT IDs are not public); the",
        "  patient T/NK row uses DRMref malignant. Those are different malignant calls.",
        "- GSE131907 CellChat includes tLung tS1–tS3 (32 samples). The patient T/NK row",
        "  uses author `Malignant cells` with n_mal≥20 (21 samples; mets-heavy).",
        "- Permutation tests on cell-pooled truncated means are not a 33-unit mixed model.",
        "- CXCL16–CXCR6 up in CLDN4-high is the opposite of an immune-cold recruit-down story.",
        "- NECTIN2–TIGIT / PD-L1–PD-1 are not detected on GSE131907.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 methods/pair_207422_131907_cellchat_cldn4/analyze.py",
        "```",
        "",
        "Requires the committed `data/` extracts only. GEO UMI matrices are not needed.",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    rows = load_pair_atomics()
    by = {(r["cohort"], r["malig_def"], r["score"]): r for r in rows}
    a207 = by[("GSE207422", "author_DRMref", "mean")]
    a131 = by[("GSE131907", "author_malig", "mean")]
    a131_pct = by[("GSE131907", "author_malig", "pct")]
    a_tlung = by.get(("GSE131907_tLung", "author_epi", "mean"))

    re_s, st_s = pool_spearman(
        [a207["rho"], a131["rho"]], [a207["n"], a131["n"]], [a207["p"], a131["p"]]
    )
    re_q, st_q = pool_q4q1(
        [a207["r_rb"], a131["r_rb"]],
        [int(a207["n_compared"]), int(a131["n_compared"])],
        [a207["p_q4q1"], a131["p_q4q1"]],
    )
    re_s_pct, _ = pool_spearman(
        [a207["rho"], a131_pct["rho"]],
        [a207["n"], a131_pct["n"]],
        [a207["p"], a131_pct["p"]],
    )
    re_q_pct, _ = pool_q4q1(
        [a207["r_rb"], a131_pct["r_rb"]],
        [int(a207["n_compared"]), int(a131_pct["n_compared"])],
        [a207["p_q4q1"], a131_pct["p_q4q1"]],
    )

    combo_rho = pd.DataFrame(
        [
            {
                "analysis": "spearman",
                "family": "primary_malig_mean",
                "cohorts": "GSE207422+GSE131907",
                "k": 2,
                "n_207422": a207["n"],
                "n_131907": a131["n"],
                "N": a207["n"] + a131["n"],
                "unit_207422": "patient",
                "unit_131907": "sample",
                "rho_207422": a207["rho"],
                "p_207422": a207["p"],
                "rho_131907": a131["rho"],
                "p_131907": a131["p"],
                "pooled_rho": re_s.get("pooled_rho"),
                "p": re_s.get("p"),
                "I2": re_s.get("I2"),
                "ci95_lo": (re_s.get("ci95_rho") or [np.nan, np.nan])[0],
                "ci95_hi": (re_s.get("ci95_rho") or [np.nan, np.nan])[1],
                "stouffer_z": st_s.get("z"),
                "stouffer_p": st_s.get("p"),
                "method": "DerSimonian-Laird RE on Fisher-z(Spearman ρ)",
                "cldn4_only": True,
                "dual_high": False,
            },
            {
                "analysis": "q4q1",
                "family": "primary_malig_mean",
                "cohorts": "GSE207422+GSE131907",
                "k": 2,
                "n_207422": int(a207["n_compared"]),
                "n_131907": int(a131["n_compared"]),
                "N": int(a207["n_compared"] + a131["n_compared"]),
                "n_q1": int(a207["n_q1"] + a131["n_q1"]),
                "n_q4": int(a207["n_q4"] + a131["n_q4"]),
                "unit_207422": "patient",
                "unit_131907": "sample",
                "r_207422": a207["r_rb"],
                "p_207422": a207["p_q4q1"],
                "r_131907": a131["r_rb"],
                "p_131907": a131["p_q4q1"],
                "pooled_r": re_q.get("pooled_rho"),
                "p": re_q.get("p"),
                "I2": re_q.get("I2"),
                "ci95_lo": (re_q.get("ci95_rho") or [np.nan, np.nan])[0],
                "ci95_hi": (re_q.get("ci95_rho") or [np.nan, np.nan])[1],
                "stouffer_z": st_q.get("z"),
                "stouffer_p": st_q.get("p"),
                "method": "DerSimonian-Laird RE on Fisher-z(rank-biserial r); Q4 vs Q1 T/NK",
                "cldn4_only": True,
                "dual_high": False,
            },
            {
                "analysis": "spearman",
                "family": "sensitivity_131907_pct",
                "cohorts": "GSE207422+GSE131907",
                "k": 2,
                "n_207422": a207["n"],
                "n_131907": a131_pct["n"],
                "N": a207["n"] + a131_pct["n"],
                "rho_207422": a207["rho"],
                "rho_131907": a131_pct["rho"],
                "pooled_rho": re_s_pct.get("pooled_rho"),
                "p": re_s_pct.get("p"),
                "I2": re_s_pct.get("I2"),
                "cldn4_only": True,
                "dual_high": False,
            },
            {
                "analysis": "q4q1",
                "family": "sensitivity_131907_pct",
                "cohorts": "GSE207422+GSE131907",
                "k": 2,
                "n_207422": int(a207["n_compared"]),
                "n_131907": int(a131_pct["n_compared"]),
                "N": int(a207["n_compared"] + a131_pct["n_compared"]),
                "n_q1": int(a207["n_q1"] + a131_pct["n_q1"]),
                "n_q4": int(a207["n_q4"] + a131_pct["n_q4"]),
                "r_207422": a207["r_rb"],
                "r_131907": a131_pct["r_rb"],
                "pooled_r": re_q_pct.get("pooled_rho"),
                "p": re_q_pct.get("p"),
                "I2": re_q_pct.get("I2"),
                "cldn4_only": True,
                "dual_high": False,
            },
        ]
    )
    combo_rho.to_csv(TABLES / "combo_rho.tsv", sep="\t", index=False)

    singles = pd.DataFrame([{k: v for k, v in r.items() if not k.startswith("_")} for r in rows])
    singles.to_csv(TABLES / "singles.tsv", sep="\t", index=False)

    out_a = load_outgoing(CCC / "GSE207422_outgoing.tsv", "GSE207422")
    out_b = load_outgoing(CCC / "GSE131907_outgoing.tsv", "GSE131907")
    lr = combo_lr(out_a, out_b)
    lr.to_csv(TABLES / "combo_lr_table.tsv", sep="\t", index=False)
    ligand = lr[lr["sig_either"] | lr["key_pair"]].copy()
    ligand.to_csv(TABLES / "ligand_table.tsv", sep="\t", index=False)

    # Patient-level per-unit table used for the combo.
    a3 = pd.read_csv(TNK / "GSE207422_drmref_patients.tsv", sep="\t")
    d = pd.read_csv(TNK / "GSE131907_samples.tsv", sep="\t")
    mal = d[
        d["origin"].isin(["tLung", "tL/B", "mLN", "PE", "mBrain"])
        & (d["n_malignant"] >= 20)
    ].copy()
    units = pd.concat(
        [
            pd.DataFrame(
                {
                    "cohort": "GSE207422",
                    "unit_id": a3["patient"],
                    "sample": a3["sample"],
                    "unit": "patient",
                    "cldn4": a3["malig_CLDN4_mean"],
                    "tnk": a3["frac_tnk"],
                    "n_malignant": a3["n_malignant"],
                    "n_tnk": a3["n_tnk"],
                }
            ),
            pd.DataFrame(
                {
                    "cohort": "GSE131907",
                    "unit_id": mal["sample"],
                    "sample": mal["sample"],
                    "unit": "sample",
                    "cldn4": mal["mal_CLDN4_mean"],
                    "tnk": mal["frac_tnk"],
                    "n_malignant": mal["n_malignant"],
                    "n_tnk": mal["n_tnk"],
                }
            ),
        ],
        ignore_index=True,
    )
    units.to_csv(TABLES / "pair_units.tsv", sep="\t", index=False)

    summary = {
        "pair": "GSE207422+GSE131907",
        "cldn4_only": True,
        "dual_high": False,
        "spearman": {
            "k": 2,
            "N": int(a207["n"] + a131["n"]),
            "pooled_rho": re_s.get("pooled_rho"),
            "p": re_s.get("p"),
            "I2": re_s.get("I2"),
            "n_207422": a207["n"],
            "n_131907": a131["n"],
            "rho_207422": a207["rho"],
            "rho_131907": a131["rho"],
        },
        "q4q1": {
            "k": 2,
            "N_compared": int(a207["n_compared"] + a131["n_compared"]),
            "n_q1": int(a207["n_q1"] + a131["n_q1"]),
            "n_q4": int(a207["n_q4"] + a131["n_q4"]),
            "pooled_r": re_q.get("pooled_rho"),
            "p": re_q.get("p"),
            "I2": re_q.get("I2"),
        },
        "lr": {
            "n_rows": int(len(lr)),
            "n_sig_either": int(lr["sig_either"].sum()),
            "n_both_same": int((lr["concordance"] == "both_sig_same").sum()),
            "n_both_discordant": int((lr["concordance"] == "both_sig_discordant").sum()),
        },
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    forest(
        [
            {"label": f"GSE207422 n={a207['n']} patients", "rho": a207["rho"], "p": a207["p"]},
            {"label": f"GSE131907 n={a131['n']} samples", "rho": a131["rho"], "p": a131["p"]},
            {
                "label": f"Pair RE n={a207['n']+a131['n']}",
                "rho": re_s.get("pooled_rho"),
                "p": re_s.get("p"),
                "pooled": True,
            },
        ],
        "Pair GSE207422 + GSE131907 · malignant CLDN4 vs T/NK (Spearman)",
        FIGS / "fig_combo_rho_forest.png",
        "Spearman ρ  (malignant CLDN4 vs T/NK; <0 = colder)",
        "rho",
        "p",
    )
    forest(
        [
            {
                "label": f"GSE207422 Q4 vs Q1 ({int(a207['n_q1'])}/{int(a207['n_q4'])})",
                "rho": a207["r_rb"],
                "p": a207["p_q4q1"],
            },
            {
                "label": f"GSE131907 Q4 vs Q1 ({int(a131['n_q1'])}/{int(a131['n_q4'])})",
                "rho": a131["r_rb"],
                "p": a131["p_q4q1"],
            },
            {
                "label": f"Pair RE compared n={int(a207['n_compared']+a131['n_compared'])}",
                "rho": re_q.get("pooled_rho"),
                "p": re_q.get("p"),
                "pooled": True,
            },
        ],
        "Pair GSE207422 + GSE131907 · CLDN4 Q4 vs Q1 T/NK",
        FIGS / "fig_combo_q4q1_forest.png",
        "Rank-biserial r  (Q4−Q1 T/NK; <0 = Q4 colder)",
        "rho",
        "p",
    )
    scatter(
        FIGS / "fig_scatter_GSE207422.png",
        a207["_y"],
        a207["_x"],
        "T/NK fraction (DRMref)",
        "Malignant CLDN4 mean log1p(CP10k)",
        f"GSE207422 · CLDN4 vs T/NK (n={a207['n']} patients)",
    )
    scatter(
        FIGS / "fig_scatter_GSE131907.png",
        a131["_y"],
        a131["_x"],
        "T/NK fraction",
        "Malignant CLDN4 mean log1p(CP10k)",
        f"GSE131907 author_malig · CLDN4 vs T/NK (n={a131['n']} samples)",
    )
    if a207["_q"]:
        box_q4q1(
            FIGS / "fig_q4q1_box_GSE207422.png",
            a207["_q"],
            "GSE207422 Q4 vs Q1 T/NK (CLDN4 mean)",
            "T/NK fraction",
        )
    if a131["_q"]:
        box_q4q1(
            FIGS / "fig_q4q1_box_GSE131907.png",
            a131["_q"],
            "GSE131907 Q4 vs Q1 T/NK (CLDN4 mean)",
            "T/NK fraction",
        )
    fig_honest_n(FIGS / "fig_honest_n.png", a207, a131)
    fig_lr_concordance(FIGS / "fig_extra_lr_concordance.png", lr)
    fig_shared_delta(FIGS / "fig_extra_shared_outgoing.png", lr)
    fig_extra_ligand_table(FIGS / "fig_extra_ligand_table.png", lr)

    write_finding(
        a207, a131, a131_pct, a_tlung, re_s, st_s, re_q, st_q, re_s_pct, re_q_pct, lr,
        HERE / "FINDING.md",
    )

    print("combo Spearman", re_s)
    print("combo Q4Q1", re_q)
    print("LR both_same", int((lr["concordance"] == "both_sig_same").sum()),
          "discordant", int((lr["concordance"] == "both_sig_discordant").sum()),
          "sig_either", int(lr["sig_either"].sum()))
    print("wrote", TABLES, FIGS, HERE / "FINDING.md")


if __name__ == "__main__":
    main()
