#!/usr/bin/env python3
"""Assemble public per-patient malignant/epithelial TACSTD2 and CLDN4, then IV meta + forest."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
INP = ROOT / "inputs"
OUT = ROOT / "results"
FIG = ROOT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

MIN_CELLS = 20
# NMPR minus MPR (or NR minus R): positive = higher in non-responders


def hedges_g(a: np.ndarray, b: np.ndarray) -> dict:
    """Hedges g for a vs b (a = NMPR/NR, b = MPR/R)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return {"n1": n1, "n2": n2, "g": np.nan, "se": np.nan, "var": np.nan}
    s1, s2 = float(np.std(a, ddof=1)), float(np.std(b, ddof=1))
    df = n1 + n2 - 2
    sp = np.sqrt(((n1 - 1) * s1**2 + (n2 - 1) * s2**2) / df) if df > 0 else np.nan
    d = (float(np.mean(a)) - float(np.mean(b))) / sp if sp and sp > 0 else 0.0
    j = 1.0 - 3.0 / (4.0 * df - 1.0) if (4.0 * df - 1.0) > 0 else 1.0
    g = j * d
    var = (n1 + n2) / (n1 * n2) + (g**2) / (2.0 * (n1 + n2))
    return {
        "n1": n1,
        "n2": n2,
        "mean_a": float(np.mean(a)),
        "mean_b": float(np.mean(b)),
        "sd_a": s1,
        "sd_b": s2,
        "g": float(g),
        "var": float(var),
        "se": float(np.sqrt(var)),
    }


def rank_biserial(a: np.ndarray, b: np.ndarray) -> dict:
    """Glass rank-biserial from MWU U (a vs b). r = 2U/(n1 n2) - 1; + = a stochastically larger."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n1, n2 = len(a), len(b)
    if n1 < 1 or n2 < 1:
        return {"U": np.nan, "p": np.nan, "r_rb": np.nan}
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided", method="auto")
    r = 2.0 * float(u) / (n1 * n2) - 1.0
    # SE for rank-biserial via U variance under null, then delta method
    var_u = n1 * n2 * (n1 + n2 + 1) / 12.0
    se_r = 2.0 * np.sqrt(var_u) / (n1 * n2)
    return {"U": float(u), "p": float(p), "r_rb": float(r), "se_rb": float(se_r), "var_rb": float(se_r**2)}


def iv_meta(effects: list[dict], key_g="g", key_var="var") -> dict:
    rows = [e for e in effects if np.isfinite(e.get(key_g, np.nan)) and np.isfinite(e.get(key_var, np.nan)) and e[key_var] > 0]
    if not rows:
        return {"k": 0}
    g = np.array([e[key_g] for e in rows], dtype=float)
    v = np.array([e[key_var] for e in rows], dtype=float)
    w = 1.0 / v
    g_fe = float(np.sum(w * g) / np.sum(w))
    se_fe = float(1.0 / np.sqrt(np.sum(w)))
    q = float(np.sum(w * (g - g_fe) ** 2))
    k = len(rows)
    df = k - 1
    p_q = float(1.0 - stats.chi2.cdf(q, df)) if df > 0 else np.nan
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w)) if df > 0 else np.nan
    tau2 = max(0.0, (q - df) / c) if c and c > 0 else 0.0
    w_re = 1.0 / (v + tau2)
    g_re = float(np.sum(w_re * g) / np.sum(w_re))
    se_re = float(1.0 / np.sqrt(np.sum(w_re)))
    i2 = max(0.0, (q - df) / q) if q > 0 and df > 0 else 0.0
    z_fe = g_fe / se_fe
    z_re = g_re / se_re
    return {
        "k": k,
        "n_patients": int(sum(e.get("n1", 0) + e.get("n2", 0) for e in rows)),
        "n_NMPR": int(sum(e.get("n1", 0) for e in rows)),
        "n_MPR": int(sum(e.get("n2", 0) for e in rows)),
        "g_fe": g_fe,
        "se_fe": se_fe,
        "ci95_fe": [g_fe - 1.96 * se_fe, g_fe + 1.96 * se_fe],
        "p_fe": float(2.0 * (1.0 - stats.norm.cdf(abs(z_fe)))),
        "g_re": g_re,
        "se_re": se_re,
        "ci95_re": [g_re - 1.96 * se_re, g_re + 1.96 * se_re],
        "p_re": float(2.0 * (1.0 - stats.norm.cdf(abs(z_re)))),
        "Q": q,
        "p_Q": p_q,
        "tau2": tau2,
        "I2": i2,
    }


def contrast_row(cohort: str, gene: str, endpoint: str, metric: str, a: np.ndarray, b: np.ndarray, note: str) -> dict:
    hg = hedges_g(a, b)
    rb = rank_biserial(a, b)
    return {
        "cohort": cohort,
        "gene": gene,
        "endpoint": endpoint,
        "metric": metric,
        "n_NMPR_or_NR": hg["n1"],
        "n_MPR_or_R": hg["n2"],
        "n_patients": hg["n1"] + hg["n2"],
        "mean_NMPR_or_NR": hg.get("mean_a", np.nan),
        "mean_MPR_or_R": hg.get("mean_b", np.nan),
        "sd_NMPR_or_NR": hg.get("sd_a", np.nan),
        "sd_MPR_or_R": hg.get("sd_b", np.nan),
        "hedges_g": hg.get("g", np.nan),
        "g_se": hg.get("se", np.nan),
        "g_var": hg.get("var", np.nan),
        "mwu_U": rb.get("U", np.nan),
        "mwu_p": rb.get("p", np.nan),
        "rank_biserial": rb.get("r_rb", np.nan),
        "rb_se": rb.get("se_rb", np.nan),
        "rb_var": rb.get("var_rb", np.nan),
        "note": note,
        "n1": hg["n1"],
        "n2": hg["n2"],
        "g": hg.get("g", np.nan),
        "var": hg.get("var", np.nan),
    }


def load_gse207422() -> tuple[pd.DataFrame, list[dict]]:
    df = pd.read_csv(INP / "GSE207422_per_sample.tsv", sep="\t")
    post = df[df["timing"] == "post"].copy()
    rows = []
    # Primary A3-given public contrast: all epithelial, post-treatment, min 20 epi cells (keeps 4 MPR)
    a = post[post["n_epithelial"] >= MIN_CELLS]
    for gene, col in [
        ("TACSTD2", "epi_TACSTD2_mean_log1p"),
        ("CLDN4", "epi_CLDN4_mean_log1p"),
    ]:
        nmpr = a.loc[a["path_response"] == "NMPR", col].astype(float).values
        mpr = a.loc[a["path_response"] == "MPR", col].astype(float).values
        rows.append(
            contrast_row(
                "GSE207422 (A3, epithelial)",
                gene,
                "MPR",
                col,
                nmpr,
                mpr,
                "Hu 2023 public UMI; post-tx only; all epithelial; pCR with MPR; A3 MPR labels taken as given",
            )
        )
    # Sensitivity: malignant-like with >=1 cell (2 MPR have 0 malignant-like)
    b = post[post["n_malig_like"] >= 1]
    for gene, col in [
        ("TACSTD2", "malig_TACSTD2_mean_log1p"),
        ("CLDN4", "malig_CLDN4_mean_log1p"),
    ]:
        nmpr = b.loc[b["path_response"] == "NMPR", col].astype(float).values
        mpr = b.loc[b["path_response"] == "MPR", col].astype(float).values
        rows.append(
            contrast_row(
                "GSE207422 malignant-like (sensitivity)",
                gene,
                "MPR",
                col,
                nmpr,
                mpr,
                "marker epithelial AND zero SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3; 2/4 MPR have 0 malignant-like cells",
            )
        )
    keep = post[["Sample", "Patient", "path_response", "n_epithelial", "n_malig_like",
                 "epi_TACSTD2_mean_log1p", "epi_CLDN4_mean_log1p",
                 "malig_TACSTD2_mean_log1p", "malig_CLDN4_mean_log1p"]].copy()
    keep.insert(0, "cohort", "GSE207422")
    return keep, rows


def load_gse241934() -> tuple[pd.DataFrame, list[dict]]:
    tac = pd.read_csv(INP / "GSE241934_scrna_per_sample.csv")
    cldn_path = INP / "GSE241934_cldn4_per_sample.csv"
    if cldn_path.exists():
        cldn = pd.read_csv(cldn_path)
        tac = tac.merge(
            cldn[["cohort", "sampleID", "malignant_CLDN4_mean_log1p", "malignant_CLDN4_frac_pos"]],
            on=["cohort", "sampleID"],
            how="left",
        )
        has_cldn = True
    else:
        tac["malignant_CLDN4_mean_log1p"] = np.nan
        has_cldn = False
    rows = []
    parts = [
        ("GSE241934 IIT", tac[tac["cohort"] == "IIT_EGFRmut"]),
        ("GSE241934 Real", tac[tac["cohort"] == "REAL_WT"]),
        ("GSE241934 IIT+Real (pooled patients)", tac),
    ]
    for label, sub in parts:
        a = sub[(sub["n_epi"] >= MIN_CELLS) & (sub["group"].isin(["MPR", "NMPR"]))]
        for gene, col in [("TACSTD2", "malignant_TACSTD2_mean_log1p"), ("CLDN4", "malignant_CLDN4_mean_log1p")]:
            if gene == "CLDN4" and not has_cldn:
                continue
            nmpr = a.loc[a["group"] == "NMPR", col].astype(float).dropna().values
            mpr = a.loc[a["group"] == "MPR", col].astype(float).dropna().values
            rows.append(
                contrast_row(
                    label,
                    gene,
                    "MPR",
                    col,
                    nmpr,
                    mpr,
                    "author Epi labels; pCR with MPR; min 20 Epi cells; mean log1p(UMI)",
                )
            )
    keep = tac.copy()
    return keep, rows


def load_gse291670() -> tuple[pd.DataFrame, list[dict]]:
    df = pd.read_csv(INP / "GSE291670_per_sample.tsv", sep="\t")
    rows = []
    for gene, col in [
        ("TACSTD2", "mal_TACSTD2_mean_log1p_cp10k"),
        ("CLDN4", "mal_CLDN4_mean_log1p_cp10k"),
    ]:
        nmpr = df.loc[df["response"] == "NMPR", col].astype(float).values
        mpr = df.loc[df["response"] == "MPR", col].astype(float).values
        rows.append(
            contrast_row(
                "GSE291670",
                gene,
                "MPR",
                col,
                nmpr,
                mpr,
                "Xia 2025; 3 vs 3; marker malignant; mean log1p(CP10k); exact MWU min p at 3v3 is 0.10",
            )
        )
    keep = df.copy()
    keep.insert(0, "cohort", "GSE291670")
    return keep, rows


def load_gse205335() -> tuple[pd.DataFrame, list[dict]]:
    df = pd.read_csv(INP / "GSE205335_per_patient_tacstd2.csv")
    mal = df[(df["compartment"] == "Malignant") & (df["n_cells"] >= MIN_CELLS)].copy()
    rows = []
    # RECIST: NR = SD/PD, R = PR
    a = mal[mal["response"].isin(["R", "NR"])]
    nmpr = a.loc[a["response"] == "NR", "mean_log1p_cp10k"].astype(float).values
    mpr = a.loc[a["response"] == "R", "mean_log1p_cp10k"].astype(float).values
    rows.append(
        contrast_row(
            "GSE205335",
            "TACSTD2",
            "RECIST",
            "mean_log1p_cp10k",
            nmpr,
            mpr,
            "Kim/Ahn 2024; author malignant; R=PR, NR=SD/PD; NE dropped; 3 R + 1 NR had 0 malignant cells",
        )
    )
    cldn_path = INP / "GSE205335_per_patient_cldn4.csv"
    if cldn_path.exists():
        c = pd.read_csv(cldn_path)
        cm = c[(c["compartment"] == "Malignant") & (c["n_cells"] >= MIN_CELLS) & (c["response"].isin(["R", "NR"]))]
        rows.append(
            contrast_row(
                "GSE205335",
                "CLDN4",
                "RECIST",
                "mean_log1p_cp10k",
                cm.loc[cm["response"] == "NR", "mean_log1p_cp10k"].astype(float).values,
                cm.loc[cm["response"] == "R", "mean_log1p_cp10k"].astype(float).values,
                "author malignant; R=PR, NR=SD/PD",
            )
        )
    return mal, rows


def forest(rows: list[dict], title: str, out_png: Path, meta: dict | None, effect_key="hedges_g", se_key="g_se") -> None:
    plot_rows = [r for r in rows if np.isfinite(r.get(effect_key, np.nan))]
    if not plot_rows:
        return
    labels = []
    xs, xse, ns = [], [], []
    for r in plot_rows:
        labels.append(f"{r['cohort']}  {r['n_NMPR_or_NR']} vs {r['n_MPR_or_R']}")
        xs.append(r[effect_key])
        xse.append(r[se_key])
        ns.append(r["n_patients"])
    if meta and meta.get("k", 0) >= 1:
        labels.append(f"IV fixed  n={meta['n_patients']}")
        xs.append(meta["g_fe"])
        xse.append(meta["se_fe"])
        labels.append(f"IV random  n={meta['n_patients']}")
        xs.append(meta["g_re"])
        xse.append(meta["se_re"])
    y = np.arange(len(labels))[::-1]
    fig, ax = plt.subplots(figsize=(8.2, 0.55 * len(labels) + 1.8))
    for i, (xi, sei, lab) in enumerate(zip(xs, xse, labels)):
        color = "#222222" if i < len(plot_rows) else "#1f4e79"
        marker = "o" if i < len(plot_rows) else "D"
        ax.errorbar(xi, y[i], xerr=1.96 * sei, fmt=marker, color=color, ecolor=color, capsize=3, ms=6)
        ax.text(xi + 1.96 * sei + 0.08, y[i], f"{xi:.2f}", va="center", fontsize=8, color=color)
    ax.axvline(0, color="#888888", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Hedges g  (NMPR/NR − MPR/R)")
    ax.set_title(title, fontsize=11)
    ax.set_xlim(min(-2.2, min(np.array(xs) - 1.96 * np.array(xse)) - 0.3),
                max(2.2, max(np.array(xs) + 1.96 * np.array(xse)) + 0.8))
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


def main() -> None:
    p207, r207 = load_gse207422()
    p241, r241 = load_gse241934()
    p291, r291 = load_gse291670()
    p205, r205 = load_gse205335()

    p207.to_csv(OUT / "per_patient_GSE207422.tsv", sep="\t", index=False)
    p241.to_csv(OUT / "per_patient_GSE241934.tsv", sep="\t", index=False)
    p291.to_csv(OUT / "per_patient_GSE291670.tsv", sep="\t", index=False)
    p205.to_csv(OUT / "per_patient_GSE205335.tsv", sep="\t", index=False)

    all_rows = r207 + r241 + r291 + r205
    table = pd.DataFrame(all_rows)
    table.to_csv(OUT / "cohort_effects.tsv", sep="\t", index=False)

    # Primary MPR meta: independent cohorts (IIT and Real separate; exclude combined and sensitivity)
    primary_names = {
        "GSE207422 (A3, epithelial)",
        "GSE241934 IIT",
        "GSE241934 Real",
        "GSE291670",
    }
    summary = {"primary_rule": "independent MPR cohorts; GSE241934 IIT and Real kept separate; GSE205335 RECIST held out"}
    for gene in ["TACSTD2", "CLDN4"]:
        prim = [r for r in all_rows if r["gene"] == gene and r["endpoint"] == "MPR" and r["cohort"] in primary_names]
        meta_g = iv_meta(prim, "hedges_g", "g_var")
        meta_rb = iv_meta(
            [{**r, "g": r["rank_biserial"], "var": r["rb_var"]} for r in prim],
            "g",
            "var",
        )
        summary[f"MPR_{gene}_hedges_g"] = meta_g
        summary[f"MPR_{gene}_rank_biserial"] = meta_rb
        forest(
            prim,
            f"{gene}  malignant/epithelial  NMPR vs MPR",
            FIG / f"forest_mpr_{gene.lower()}.png",
            meta_g,
        )
        # also IIT+Real as one row instead of two (sensitivity, no double-count)
        alt_names = {"GSE207422 (A3, epithelial)", "GSE241934 IIT+Real (pooled patients)", "GSE291670"}
        alt = [r for r in all_rows if r["gene"] == gene and r["endpoint"] == "MPR" and r["cohort"] in alt_names]
        summary[f"MPR_{gene}_hedges_g_IITRealPooled"] = iv_meta(alt, "hedges_g", "g_var")

    recist = [r for r in all_rows if r["endpoint"] == "RECIST"]
    pd.DataFrame(recist).to_csv(OUT / "recist_GSE205335.tsv", sep="\t", index=False)
    if recist:
        forest(
            [r for r in recist if r["gene"] == "TACSTD2"],
            "GSE205335  malignant TACSTD2  NR vs R (RECIST)",
            FIG / "forest_recist_gse205335.png",
            None,
        )

    # patient totals (no double-count of GSE241934)
    n_map = {}
    for r in all_rows:
        if r["cohort"] in primary_names and r["gene"] == "TACSTD2" and r["endpoint"] == "MPR":
            n_map[r["cohort"]] = r["n_patients"]
    summary["n_patients_primary_MPR"] = int(sum(n_map.values()))
    summary["n_patients_by_cohort"] = n_map
    summary["footnote_GSE243013"] = (
        "GSE243013 (Liu 2025, Cell): public MPR labels, n=243 chemo-IO, CD45-sorted immune MTX only. "
        "0 malignant/epithelial cells in the public matrix. Not entered in the malignant table."
    )
    summary["leftover_raw_only"] = (
        "Molecular Cancer 2025 neoadjuvant scRNA (3 NMPR / 2 MPR / 1 pCR) is PRJNA1068179 raw SRA only; "
        "no public processed count matrix. Not entered."
    )

    (OUT / "meta_summary.json").write_text(json.dumps(summary, indent=2, default=float))
    # compact human table
    show = table[
        [
            "cohort",
            "gene",
            "endpoint",
            "n_NMPR_or_NR",
            "n_MPR_or_R",
            "n_patients",
            "mean_NMPR_or_NR",
            "mean_MPR_or_R",
            "hedges_g",
            "g_se",
            "mwu_p",
            "rank_biserial",
            "note",
        ]
    ]
    show.to_csv(OUT / "cohort_effects_compact.tsv", sep="\t", index=False)
    print(show.to_string(index=False))
    print(json.dumps({k: summary[k] for k in summary if k.startswith("n_") or k.startswith("MPR_TACSTD2_hedges")}, indent=2, default=float))


if __name__ == "__main__":
    main()
