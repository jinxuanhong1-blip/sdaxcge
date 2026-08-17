#!/usr/bin/env python3
"""CLDN4-first gates (median; Q4 vs Q1) vs T/NK, CD8A, ImmuneScore, CXCL13.

Uses existing public patient-level scores only. TACSTD2 is a companion column.
Does not require TACSTD2-high. Does not re-audit GSE207422 A3.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lib_stats import cliffs_delta, mwu, random_effects_dl, spearman

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "existing"
OUT = HERE / "results"
FIG = HERE / "figures"
TAB = HERE / "tables"
for d in (OUT, FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

# Barrier / immune-low thesis: CLDN4-high should sit lower on immune axes.
HIGH_C = "#C0392B"
LOW_C = "#1A7A6D"
MID_C = "#7F8C8D"


def _gate(cldn4: pd.Series) -> pd.DataFrame:
    x = pd.to_numeric(cldn4, errors="coerce")
    med = float(x.median())
    q1 = float(x.quantile(0.25))
    q3 = float(x.quantile(0.75))
    out = pd.DataFrame({"cldn4": x})
    out["median_gate"] = np.where(x.isna(), None, np.where(x >= med, "high", "low"))
    q = np.full(len(x), None, dtype=object)
    ok = x.notna().to_numpy()
    xv = x.to_numpy()
    q[ok & (xv >= q3)] = "Q4"
    q[ok & (xv <= q1)] = "Q1"
    q[ok & (xv > q1) & (xv < q3)] = "mid"
    out["q_gate"] = q
    out.attrs["median"] = med
    out.attrs["q1"] = q1
    out.attrs["q3"] = q3
    return out


def _pair(cohort, modality, cldn4_def, endpoint, endpoint_col, gate, n_total, high, low, rho, rp, rn, tacstd2_rho=np.nan):
    _, p_med, nh, nl = mwu(high, low)
    delta = cliffs_delta(high, low)
    return {
        "cohort": cohort,
        "modality": modality,
        "cldn4_def": cldn4_def,
        "endpoint": endpoint,
        "endpoint_col": endpoint_col,
        "gate": gate,
        "n_total": int(n_total),
        "n_high": int(nh) if nh == nh else 0,
        "n_low": int(nl) if nl == nl else 0,
        "cliff_delta": delta,
        "mwu_p": p_med,
        "spearman_rho": rho,
        "spearman_p": rp,
        "spearman_n": rn,
        "tacstd2_companion_rho": tacstd2_rho,
        "median_high": float(np.nanmedian(high)) if len(high) else float("nan"),
        "median_low": float(np.nanmedian(low)) if len(low) else float("nan"),
    }


def _endpoint_rows(cohort, modality, cldn4_def, df, cldn4_col, endpoints, tacstd2_col=None):
    rows = []
    g = _gate(df[cldn4_col])
    work = df.copy()
    work["_median_gate"] = g["median_gate"].values
    work["_q_gate"] = g["q_gate"].values
    work["_cldn4"] = g["cldn4"].values
    tac_rho = float("nan")
    if tacstd2_col and tacstd2_col in work.columns:
        tac_rho, _, _ = spearman(work["_cldn4"], work[tacstd2_col])
    for endpoint, col in endpoints:
        if col not in work.columns:
            continue
        y = pd.to_numeric(work[col], errors="coerce")
        pair = work["_cldn4"].notna() & y.notna()
        rho, rp, rn = spearman(work.loc[pair, "_cldn4"], y[pair])
        for gate, mask_h, mask_l in (
            ("median", work["_median_gate"] == "high", work["_median_gate"] == "low"),
            ("Q4_vs_Q1", work["_q_gate"] == "Q4", work["_q_gate"] == "Q1"),
        ):
            high = y[mask_h & pair].to_numpy()
            low = y[mask_l & pair].to_numpy()
            rows.append(
                _pair(
                    cohort, modality, cldn4_def, endpoint, col, gate,
                    int(pair.sum()), high, low, rho, rp, rn, tac_rho,
                )
            )
    return rows, work


def load_cohorts() -> tuple[list[dict], dict[str, pd.DataFrame]]:
    rows = []
    frames = {}

    # GSE207422 A3 given 12-patient table. TACSTD2 slide taken as given; CLDN4 from same table.
    a3 = pd.read_csv(DATA / "GSE207422_drmref_patients.tsv", sep="\t")
    r, f = _endpoint_rows(
        "GSE207422", "scRNA", "malig_CLDN4_mean (A3 given table)",
        a3, "malig_CLDN4_mean",
        [("T/NK", "frac_tnk")],
        tacstd2_col="malig_TACSTD2_mean",
    )
    rows += r
    frames["GSE207422"] = f.assign(patient=a3["patient"], endpoint_TNK=a3["frac_tnk"], cldn4=a3["malig_CLDN4_mean"])

    # GSE253013 existing 9-patient tumor extract
    g253 = pd.read_csv(DATA / "GSE253013_patients.tsv", sep="\t")
    g253 = g253[(g253["tissue"] == "Tumor") & (g253["eligible_malig"])].copy()
    r, f = _endpoint_rows(
        "GSE253013", "scRNA", "malig_CLDN4_mean_log1p (extract)",
        g253, "CLDN4_mean_log1p",
        [("T/NK", "tnk_fraction")],
        tacstd2_col="TACSTD2_mean_log1p",
    )
    rows += r
    frames["GSE253013"] = f.assign(patient=g253["patient"], endpoint_TNK=g253["tnk_fraction"], cldn4=g253["CLDN4_mean_log1p"])

    # GSE131907 tumor-site epithelium (same filter as scrna_meta_tnk: n_epi>=20)
    g131 = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    tumor = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
    g131 = g131[g131["origin"].isin(tumor) & (g131["n_epithelial"] >= 20)].copy()
    r, f = _endpoint_rows(
        "GSE131907", "scRNA", "epi_CLDN4_mean (tumor sites)",
        g131, "epi_CLDN4_mean",
        [("T/NK", "frac_tnk"), ("CD8A", "frac_cd8"), ("ImmuneScore", "frac_immune")],
        tacstd2_col="epi_TACSTD2_mean",
    )
    rows += r
    frames["GSE131907"] = f.assign(
        patient=g131["sample"], endpoint_TNK=g131["frac_tnk"], endpoint_CD8A=g131["frac_cd8"],
        endpoint_ImmuneScore=g131["frac_immune"], cldn4=g131["epi_CLDN4_mean"],
    )

    # GSE241934 IIT / Real (eligible Real n=24 as in scrna_meta_tnk)
    iit = pd.read_csv(DATA / "GSE241934_IIT_patients.tsv", sep="\t")
    r, f = _endpoint_rows(
        "GSE241934_IIT", "scRNA", "epi_CLDN4_log1p_cp10k",
        iit, "CLDN4_log1p_cp10k",
        [("T/NK", "frac_tnk"), ("CD8A", "frac_cd8")],
        tacstd2_col="TACSTD2_log1p_cp10k",
    )
    rows += r
    frames["GSE241934_IIT"] = f.assign(patient=iit["patient"], endpoint_TNK=iit["frac_tnk"], endpoint_CD8A=iit["frac_cd8"], cldn4=iit["CLDN4_log1p_cp10k"])

    real = pd.read_csv(DATA / "GSE241934_Real_patients.tsv", sep="\t")
    real = real[real["eligible"] == 1].copy()
    r, f = _endpoint_rows(
        "GSE241934_Real", "scRNA", "epi_CLDN4_log1p_cp10k",
        real, "CLDN4_log1p_cp10k",
        [("T/NK", "frac_tnk"), ("CD8A", "frac_cd8")],
        tacstd2_col="TACSTD2_log1p_cp10k",
    )
    rows += r
    frames["GSE241934_Real"] = f.assign(patient=real["patient"], endpoint_TNK=real["frac_tnk"], endpoint_CD8A=real["frac_cd8"], cldn4=real["CLDN4_log1p_cp10k"])

    # GSE148071 TISCH immune fraction (ImmuneScore analog) + TLS CXCL13
    tisch = pd.read_csv(DATA / "GSE148071_tisch_immune.csv")
    r, f = _endpoint_rows(
        "GSE148071", "scRNA", "malig_CLDN4_mean (TISCH)",
        tisch, "CLDN4_mean_malignant",
        [("ImmuneScore", "immune_fraction")],
        tacstd2_col="TACSTD2_mean_malignant",
    )
    rows += r
    frames["GSE148071"] = f.assign(patient=tisch["patient"], endpoint_ImmuneScore=tisch["immune_fraction"], cldn4=tisch["CLDN4_mean_malignant"])

    tls = pd.read_csv(DATA / "scrna_tls_patient_scores.tsv", sep="\t")
    tls = tls[tls["eligible"]].copy()
    for cohort, label in (
        ("GSE148071", "GSE148071_tls"),
        ("GSE207422", "GSE207422_tls"),
        ("GSE253013", "GSE253013_tls"),
        ("GSE131907", "GSE131907_tls"),
        ("GSE241934_IIT", "GSE241934_IIT_tls"),
    ):
        sub = tls[tls["cohort"] == cohort].copy()
        if sub.empty:
            continue
        r, f = _endpoint_rows(
            label, "scRNA", "mal_CLDN4_mean_log1p_cp10k (TLS extract)",
            sub, "cldn4_mal_mean_log1p_cp10k",
            [("CXCL13", "frac_CXCL13pos_T"), ("CXCL13_mean", "cxcl13_mean_log1p_cp10k")],
            tacstd2_col="tacstd2_mal_mean_log1p_cp10k",
        )
        rows += r
        frames[label] = f.assign(
            patient=sub["patient"], endpoint_CXCL13=sub["frac_CXCL13pos_T"],
            cldn4=sub["cldn4_mal_mean_log1p_cp10k"],
        )

    # GSE218989 bulk ICI
    g218 = pd.read_csv(DATA / "GSE218989_per_patient.tsv", sep="\t")
    r, f = _endpoint_rows(
        "GSE218989", "bulk", "CLDN4_log2tpm1",
        g218, "CLDN4_log2tpm1",
        [("CD8A", "CD8A_log2tpm1"), ("ImmuneScore", "ESTIMATE_Immune_meanz")],
        tacstd2_col="TACSTD2_log2tpm1",
    )
    rows += r
    frames["GSE218989"] = f.assign(
        patient=g218["patient"], endpoint_CD8A=g218["CD8A_log2tpm1"],
        endpoint_ImmuneScore=g218["ESTIMATE_Immune_meanz"], cldn4=g218["CLDN4_log2tpm1"],
    )

    # GSE285029 bulk PD-L1
    g285 = pd.read_csv(DATA / "GSE285029_sample_scores.csv")
    if "Unnamed: 0" in g285.columns:
        g285 = g285.rename(columns={"Unnamed: 0": "patient"})
    r, f = _endpoint_rows(
        "GSE285029", "bulk", "CLDN4_log2p1",
        g285, "CLDN4",
        [("CD8A", "CD8A"), ("ImmuneScore", "immune8")],
        tacstd2_col="TACSTD2",
    )
    rows += r
    frames["GSE285029"] = f.assign(
        patient=g285["patient"], endpoint_CD8A=g285["CD8A"],
        endpoint_ImmuneScore=g285["immune8"], cldn4=g285["CLDN4"],
    )

    # CPTAC LUAD RNA (primary) + protein companion
    luad = pd.read_csv(DATA / "CPTAC_LUAD_sample_table.tsv", sep="\t")
    r, f = _endpoint_rows(
        "CPTAC_LUAD", "bulk_RNA", "CLDN4_RNA",
        luad, "CLDN4_RNA",
        [("CD8A", "CIBERSORT_T_cell_CD8+"), ("ImmuneScore", "ESTIMATE_ImmuneScore"), ("CXCL13", "TLS_chemokine_RNA")],
        tacstd2_col="TACSTD2_RNA",
    )
    rows += r
    frames["CPTAC_LUAD"] = f.assign(
        patient=luad["case_id"], endpoint_CD8A=luad["CIBERSORT_T_cell_CD8+"],
        endpoint_ImmuneScore=luad["ESTIMATE_ImmuneScore"], endpoint_CXCL13=luad["TLS_chemokine_RNA"],
        cldn4=luad["CLDN4_RNA"],
    )
    r, _ = _endpoint_rows(
        "CPTAC_LUAD_protein", "protein", "CLDN4_protein",
        luad, "CLDN4_protein",
        [("CD8A", "CIBERSORT_T_cell_CD8+"), ("ImmuneScore", "ESTIMATE_ImmuneScore"), ("CXCL13", "TLS_chemokine_RNA")],
        tacstd2_col="TACSTD2_protein",
    )
    rows += r

    lscc = pd.read_csv(DATA / "CPTAC_LSCC_sample_level_features.tsv", sep="\t")
    if "Unnamed: 0" in lscc.columns:
        lscc = lscc.rename(columns={"Unnamed: 0": "patient"})
    r, f = _endpoint_rows(
        "CPTAC_LSCC", "bulk_RNA", "CLDN4_RNA",
        lscc, "rna_CLDN4",
        [("CD8A", "CIBERSORT_T_cell_CD8+"), ("ImmuneScore", "ESTIMATE_ImmuneScore")],
        tacstd2_col="rna_TACSTD2",
    )
    rows += r
    frames["CPTAC_LSCC"] = f.assign(
        patient=lscc["patient"], endpoint_CD8A=lscc["CIBERSORT_T_cell_CD8+"],
        endpoint_ImmuneScore=lscc["ESTIMATE_ImmuneScore"], cldn4=lscc["rna_CLDN4"],
    )
    r, _ = _endpoint_rows(
        "CPTAC_LSCC_protein", "protein", "CLDN4_protein",
        lscc, "protein_CLDN4",
        [("CD8A", "CIBERSORT_T_cell_CD8+"), ("ImmuneScore", "ESTIMATE_ImmuneScore")],
        tacstd2_col="protein_TACSTD2",
    )
    rows += r

    return rows, frames


def _fmt_p(p):
    if p != p:
        return "NA"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.3g}"


def _fmt_d(x, nd=3):
    if x != x:
        return "NA"
    return f"{x:+.{nd}f}"


def plot_forest(df: pd.DataFrame, gate: str, path: Path, title: str, drop_protein_tls: bool = True):
    sub = df[(df["gate"] == gate) & (df["n_high"] >= 3) & (df["n_low"] >= 3)].copy()
    # drop CXCL13_mean companion to keep the forest readable
    sub = sub[sub["endpoint"] != "CXCL13_mean"]
    # primary rows only (drop protein companion from the main forest)
    if drop_protein_tls:
        sub = sub[~sub["cohort"].str.contains("protein|_tls")]
    sub = sub.sort_values(["endpoint", "cliff_delta"])
    if sub.empty:
        return
    labels = [
        f"{r.cohort}  {r.endpoint}  n={r.n_high}/{r.n_low}"
        for r in sub.itertuples()
    ]
    y = np.arange(len(sub))
    fig_h = max(4.2, 0.38 * len(sub) + 1.6)
    fig, ax = plt.subplots(figsize=(9.2, fig_h))
    colors = [HIGH_C if d < 0 else "#2C3E50" for d in sub["cliff_delta"]]
    ax.axvline(0, color="#222", lw=1.0, zorder=0)
    ax.axvline(-0.33, color="#bbb", ls="--", lw=0.7, zorder=0)
    ax.axvline(-0.474, color="#bbb", ls=":", lw=0.7, zorder=0)
    ax.scatter(sub["cliff_delta"], y, s=64, c=colors, zorder=3, edgecolors="white", linewidths=0.4)
    for i, (d, p) in enumerate(zip(sub["cliff_delta"], sub["mwu_p"])):
        star = " **" if p == p and p < 0.01 else (" *" if p == p and p < 0.05 else "")
        ax.text(d + (0.03 if d >= 0 else -0.03), i, f"{d:+.2f}{star}", va="center",
                ha="left" if d >= 0 else "right", fontsize=7.5, color="#222")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Cliff δ  (CLDN4-high − CLDN4-low)   negative = immune-low")
    ax.set_title(title)
    ax.set_xlim(-1.15, 1.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=180)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_violins(frames: dict[str, pd.DataFrame], specs: list[tuple], path: Path, title: str):
    """specs: (frame_key, endpoint_col, panel_title)"""
    n = len(specs)
    fig, axes = plt.subplots(1, n, figsize=(3.15 * n, 4.4), sharey=False)
    if n == 1:
        axes = [axes]
    rng = np.random.default_rng(7)
    for ax, (key, ecol, lab) in zip(axes, specs):
        df = frames[key]
        y = pd.to_numeric(df[ecol], errors="coerce")
        q = df["_q_gate"]
        low = y[q == "Q1"].dropna().to_numpy()
        high = y[q == "Q4"].dropna().to_numpy()
        parts = ax.violinplot([low, high], positions=[0, 1], showmeans=False, showmedians=True, widths=0.72)
        for i, body in enumerate(parts["bodies"]):
            body.set_facecolor(LOW_C if i == 0 else HIGH_C)
            body.set_alpha(0.78)
            body.set_edgecolor("white")
        for k in ("cbars", "cmins", "cmaxes", "cmedians"):
            parts[k].set_color("#222")
            parts[k].set_linewidth(1.0)
        ax.scatter(rng.normal(0, 0.045, size=low.size), low, s=10, c="#0E3D38", zorder=4, alpha=0.7)
        ax.scatter(rng.normal(1, 0.045, size=high.size), high, s=10, c="#5B130C", zorder=4, alpha=0.7)
        d = cliffs_delta(high, low)
        _, p, nh, nl = mwu(high, low)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"CLDN4 Q1\nn={nl}", f"CLDN4 Q4\nn={nh}"], fontsize=8)
        ax.set_title(f"{lab}\nδ={d:+.2f}  p={_fmt_p(p)}", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(title, y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=180, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_median_violins(frames, specs, path, title):
    n = len(specs)
    fig, axes = plt.subplots(1, n, figsize=(3.15 * n, 4.4))
    if n == 1:
        axes = [axes]
    rng = np.random.default_rng(11)
    for ax, (key, ecol, lab) in zip(axes, specs):
        df = frames[key]
        y = pd.to_numeric(df[ecol], errors="coerce")
        g = df["_median_gate"]
        low = y[g == "low"].dropna().to_numpy()
        high = y[g == "high"].dropna().to_numpy()
        parts = ax.violinplot([low, high], positions=[0, 1], showmeans=False, showmedians=True, widths=0.72)
        for i, body in enumerate(parts["bodies"]):
            body.set_facecolor(LOW_C if i == 0 else HIGH_C)
            body.set_alpha(0.78)
            body.set_edgecolor("white")
        for k in ("cbars", "cmins", "cmaxes", "cmedians"):
            parts[k].set_color("#222")
            parts[k].set_linewidth(1.0)
        ax.scatter(rng.normal(0, 0.045, size=low.size), low, s=9, c="#0E3D38", zorder=4, alpha=0.65)
        ax.scatter(rng.normal(1, 0.045, size=high.size), high, s=9, c="#5B130C", zorder=4, alpha=0.65)
        d = cliffs_delta(high, low)
        _, p, nh, nl = mwu(high, low)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"CLDN4-low\nn={nl}", f"CLDN4-high\nn={nh}"], fontsize=8)
        ax.set_title(f"{lab}\nδ={d:+.2f}  p={_fmt_p(p)}", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(title, y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=180, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    rows, frames = load_cohorts()
    stats = pd.DataFrame(rows)
    stats.to_csv(TAB / "per_cohort_gates.tsv", sep="\t", index=False)

    # Full-pool continuous Spearman: one row per endpoint (RE of per-cohort ρ). Not the answer.
    pool_rows = []
    # One continuous-Spearman row per endpoint. TLS CXCL13 is the CXCL13 family.
    # Protein companions stay out of the pool. This row is not the answer.
    primary = stats[
        (stats["gate"] == "median")
        & (~stats["cohort"].str.contains("protein"))
        & (stats["endpoint"].isin(["T/NK", "CD8A", "ImmuneScore", "CXCL13"]))
    ].copy()
    primary = primary[~((primary["endpoint"] == "CXCL13") & (~primary["cohort"].str.contains("tls|CPTAC_LUAD")))]
    primary = primary.drop_duplicates(["cohort", "endpoint"])
    for endpoint, g in primary.groupby("endpoint"):
        g = g[g["spearman_n"] >= 4]
        re = random_effects_dl(g["spearman_rho"].tolist(), g["spearman_n"].tolist())
        pool_rows.append({
            "row": "full_pool_continuous_spearman",
            "endpoint": endpoint,
            "k": re.get("k", 0),
            "N": re.get("n_patients_total", 0),
            "pooled_rho": re.get("pooled_rho", float("nan")),
            "p": re.get("p", float("nan")),
            "I2": re.get("I2", float("nan")),
            "ci95_lo": re.get("ci95_rho", [float("nan"), float("nan")])[0] if re.get("k") else float("nan"),
            "ci95_hi": re.get("ci95_rho", [float("nan"), float("nan")])[1] if re.get("k") else float("nan"),
            "cohorts": "+".join(g["cohort"].tolist()),
            "note": "one row, not the answer",
        })
    pool = pd.DataFrame(pool_rows)
    pool.to_csv(TAB / "full_pool_spearman_one_row.tsv", sep="\t", index=False)

    # Figures — Q4 vs Q1 is the visually large cut
    plot_forest(
        stats, "Q4_vs_Q1", FIG / "forest_q4q1_cliff",
        "CLDN4 Q4 vs Q1  ·  Cliff δ on immune endpoints\n(negative = CLDN4-high is immune-low)",
    )
    plot_forest(
        stats, "median", FIG / "forest_median_cliff",
        "CLDN4 median-high vs median-low  ·  Cliff δ on immune endpoints",
    )
    # CXCL13 forest keeps TLS extracts (strongest CLDN4-first cut)
    cx = stats[(stats["endpoint"] == "CXCL13") & (stats["n_high"] >= 3) & (stats["n_low"] >= 3)].copy()
    plot_forest(
        cx.assign(cohort=cx["cohort"].str.replace("_tls", "", regex=False)),
        "Q4_vs_Q1", FIG / "forest_q4q1_cxcl13",
        "CLDN4 Q4 vs Q1  ·  CXCL13 (TLS extract + CPTAC LUAD TLS chemokine)",
    )
    plot_forest(
        cx.assign(cohort=cx["cohort"].str.replace("_tls", "", regex=False)),
        "median", FIG / "forest_median_cxcl13",
        "CLDN4 median-high vs median-low  ·  CXCL13",
    )
    prot = stats[stats["cohort"].str.contains("protein")].copy()
    prot["cohort"] = prot["cohort"].str.replace("_protein", " protein")
    plot_forest(
        prot, "Q4_vs_Q1", FIG / "forest_q4q1_protein",
        "CLDN4 protein Q4 vs Q1  ·  companion (not the RNA gate)",
        drop_protein_tls=False,
    )
    plot_violins(
        frames,
        [
            ("GSE218989", "endpoint_CD8A", "GSE218989  CD8A"),
            ("GSE218989", "endpoint_ImmuneScore", "GSE218989  ImmuneScore"),
            ("GSE285029", "endpoint_CD8A", "GSE285029  CD8A"),
            ("CPTAC_LUAD", "endpoint_CD8A", "CPTAC LUAD  CD8"),
        ],
        FIG / "violin_q4q1_bulk_cd8_immune",
        "CLDN4 Q4 vs Q1  ·  bulk CD8 / ImmuneScore",
    )
    plot_violins(
        frames,
        [
            ("GSE131907", "endpoint_TNK", "GSE131907  T/NK"),
            ("GSE148071_tls", "endpoint_CXCL13", "GSE148071  CXCL13+ T"),
            ("GSE207422_tls", "endpoint_CXCL13", "GSE207422  CXCL13+ T"),
            ("CPTAC_LUAD", "endpoint_CXCL13", "CPTAC LUAD  TLS/CXCL13"),
        ],
        FIG / "violin_q4q1_tnk_cxcl13",
        "CLDN4 Q4 vs Q1  ·  T/NK and CXCL13",
    )
    plot_median_violins(
        frames,
        [
            ("GSE218989", "endpoint_CD8A", "GSE218989  CD8A"),
            ("GSE131907", "endpoint_TNK", "GSE131907  T/NK"),
            ("GSE148071_tls", "endpoint_CXCL13", "GSE148071  CXCL13+ T"),
            ("CPTAC_LUAD", "endpoint_CD8A", "CPTAC LUAD  CD8"),
        ],
        FIG / "violin_median_selected",
        "CLDN4 median-high vs median-low",
    )

    # compact JSON for FINDING
    summary = {
        "policy": "CLDN4-first. Gates are CLDN4-high vs CLDN4-low (median and Q4 vs Q1). TACSTD2 is companion only.",
        "a3": "GSE207422 A3 TACSTD2 slide taken as given. CLDN4 from the same 12-patient table.",
        "n_stat_rows": int(len(stats)),
        "full_pool": pool.to_dict(orient="records"),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(stats.to_string(index=False))
    print("\nFULL-POOL (one row, not the answer)")
    print(pool.to_string(index=False))


if __name__ == "__main__":
    main()
