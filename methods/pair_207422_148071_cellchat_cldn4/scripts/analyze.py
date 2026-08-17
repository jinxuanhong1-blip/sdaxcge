#!/usr/bin/env python3
"""ADDITIVE pairwise merge: GSE207422 + GSE148071, CLDN4 only.

1. Patient-level malignant CLDN4 vs T/NK (per cohort + combo rho).
2. CellChat-style outgoing CLDN4-high → T/NK (given ligand tables; consensus).

Prior CellChat folders are taken as given. Matrices are not re-downloaded.
TACSTD2 is not used. No dual-high.
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from lib_stats import random_effects_dl, spearman, stouffer  # noqa: E402

DATA = ROOT / "data"
OUT = ROOT / "results"
FIG = ROOT / "figures"

MIN_EPI = 25
MIN_TNK = 25
TISCH_MIN = 20

KEY_PAIRS = [
    "NECTIN2_TIGIT",
    "LGALS9_PTPRC",
    "LGALS9_CD44",
    "HLA-E_CD8A",
    "HLA-E_CD8B",
    "F11R_ITGAL_ITGB2",
    "JAM1_ITGAL_ITGB2",
    "CDH1_ITGAE_ITGB7",
    "CXCL16_CXCR6",
    "CD274_PDCD1",
    "PVR_TIGIT",
    "MDK_NCL",
    "CD8A_CEACAM5",
    "IFNG_IFNGR1_IFNGR2",
    "TNF_TNFRSF1A",
]


def _rho_row(family: str, cohort: str, definition: str, x, y, extra: dict | None = None) -> dict:
    rho, p, n = spearman(x, y)
    if n > 3 and np.isfinite(rho):
        z = float(np.arctanh(np.clip(rho, -0.999999, 0.999999)))
        se = 1.0 / np.sqrt(n - 3)
        ci95_lo = float(np.tanh(z - 1.96 * se))
        ci95_hi = float(np.tanh(z + 1.96 * se))
    else:
        ci95_lo = ci95_hi = float("nan")
    rec = {
        "family": family,
        "cohort": cohort,
        "definition": definition,
        "n": n,
        "rho": rho,
        "p": p,
        "k": 1,
        "I2": 0.0,
        "ci95_lo": ci95_lo,
        "ci95_hi": ci95_hi,
        "method": "Spearman (patient unit)",
    }
    if extra:
        rec.update(extra)
    return rec


def _combo_row(family: str, definition: str, singles: list[dict], extra: dict | None = None) -> dict:
    rhos = [r["rho"] for r in singles]
    ns = [r["n"] for r in singles]
    ps = [r["p"] for r in singles]
    re = random_effects_dl(rhos, ns)
    st = stouffer(rhos, ps, ns)
    rec = {
        "family": family,
        "cohort": "+".join(r["cohort"] for r in singles),
        "definition": definition,
        "n": int(re.get("n_patients_total", sum(ns))),
        "rho": re.get("pooled_rho", float("nan")),
        "p": re.get("p", float("nan")),
        "k": re.get("k", len(singles)),
        "I2": re.get("I2", float("nan")),
        "ci95_lo": (re.get("ci95_rho") or [np.nan, np.nan])[0],
        "ci95_hi": (re.get("ci95_rho") or [np.nan, np.nan])[1],
        "method": re.get("method", "DL RE Fisher-z"),
        "stouffer_p": st.get("p", float("nan")),
        "stouffer_z": st.get("z", float("nan")),
        "fixed_rho": re.get("fixed_rho", float("nan")),
        "Q": re.get("Q", float("nan")),
        "tau2": re.get("tau2", float("nan")),
    }
    if extra:
        rec.update(extra)
    return rec


def cellchat_patients() -> pd.DataFrame:
    a = pd.read_csv(DATA / "gse148071_cellchat_per_sample.tsv", sep="\t")
    a["dataset"] = "GSE148071"
    a["response"] = ""
    a["eligible"] = a["eligible"].astype(bool)
    a["frac_tnk"] = a["n_TNK"] / a["n_total"]
    a["frac_tnk_among_epi_tnk"] = a["n_TNK"] / (a["n_epithelial"] + a["n_TNK"]).replace(0, np.nan)
    a["cldn4"] = a["mean_CLDN4_epithelial"]
    a["patient"] = a["sample"]

    b = pd.read_csv(DATA / "gse207422_cellchat_per_sample.tsv", sep="\t")
    b["dataset"] = "GSE207422"
    b["eligible"] = (b["n_epithelial"] >= MIN_EPI) & (b["n_TNK"] >= MIN_TNK)
    b["frac_tnk"] = b["n_TNK"] / b["n_total"]
    b["frac_tnk_among_epi_tnk"] = b["n_TNK"] / (b["n_epithelial"] + b["n_TNK"]).replace(0, np.nan)
    b["cldn4"] = b["mean_CLDN4_epithelial"]
    b["patient"] = b["sample"]
    cols = [
        "dataset",
        "patient",
        "sample",
        "response",
        "n_epithelial",
        "n_TNK",
        "n_total",
        "eligible",
        "cldn4",
        "frac_tnk",
        "frac_tnk_among_epi_tnk",
        "frac_CLDN4_pos_epithelial",
    ]
    for df in (a, b):
        for c in cols:
            if c not in df.columns:
                df[c] = np.nan
    return pd.concat([a[cols], b[cols]], ignore_index=True)


def given_patients() -> tuple[pd.DataFrame, pd.DataFrame]:
    tisch = pd.read_csv(DATA / "gse148071_tisch_units.tsv", sep="\t")
    tisch["dataset"] = "GSE148071"
    tisch["patient"] = tisch["patient"]
    tisch["cldn4"] = tisch["CLDN4_epi_mean"]
    tisch["frac_tnk"] = tisch["frac_tnk"]
    tisch["n_malignant"] = tisch["n_malignant"]
    tisch["n_scored_epi"] = tisch["n_scored_epi"]
    tisch["n_tnk"] = tisch["n_tnk"]
    tisch["eligible_locked"] = tisch["eligible"].astype(bool)
    tisch["malignant_only"] = tisch["epi_definition"].astype(str).eq("Malignant") & tisch["eligible_locked"]

    drm = pd.read_csv(DATA / "gse207422_drmref_patients.tsv", sep="\t")
    drm["dataset"] = "GSE207422"
    drm["cldn4"] = drm["malig_CLDN4_mean"]
    drm["n_tnk"] = drm["n_tnk"]
    drm["n_malignant"] = drm["n_malignant"]
    drm["eligible"] = (drm["n_malignant"] >= MIN_EPI) & (drm["n_tnk"] >= MIN_TNK)
    return tisch, drm


def build_rho(patients: pd.DataFrame, tisch: pd.DataFrame, drm: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    rows: list[dict] = []

    # Primary: CellChat marker-argmax, eligible ≥25/25, T/NK / all cells
    p148 = patients[(patients.dataset == "GSE148071") & patients.eligible]
    p422 = patients[(patients.dataset == "GSE207422") & patients.eligible]
    r148 = _rho_row(
        "primary_cellchat",
        "GSE148071",
        "marker_epi_mean vs TNK/all; eligible ≥25/25",
        p148.cldn4,
        p148.frac_tnk,
        {"n_deposited": 42, "n_eligible": int(p148.shape[0]), "note": "putative malignant = epithelium"},
    )
    r422 = _rho_row(
        "primary_cellchat",
        "GSE207422",
        "marker_epi_mean vs TNK/all; post ≥25/25",
        p422.cldn4,
        p422.frac_tnk,
        {"n_deposited": 15, "n_eligible": int(p422.shape[0]), "note": "12 post; 3 pre excluded in given CellChat"},
    )
    rows.extend([r148, r422])
    combo_primary = _combo_row(
        "primary_cellchat",
        "marker_epi_mean vs TNK/all; eligible ≥25/25",
        [r148, r422],
        {"n_deposited": 42 + 15, "n_eligible": int(p148.shape[0] + p422.shape[0]), "note": "DL RE; include 207422"},
    )
    rows.append(combo_primary)

    # Same patients, T/NK among epi+TNK only (composition of the two scored compartments)
    r148b = _rho_row(
        "cellchat_among_epi_tnk",
        "GSE148071",
        "marker_epi_mean vs TNK/(epi+TNK); ≥25/25",
        p148.cldn4,
        p148.frac_tnk_among_epi_tnk,
    )
    r422b = _rho_row(
        "cellchat_among_epi_tnk",
        "GSE207422",
        "marker_epi_mean vs TNK/(epi+TNK); post ≥25/25",
        p422.cldn4,
        p422.frac_tnk_among_epi_tnk,
    )
    rows.extend([r148b, r422b, _combo_row("cellchat_among_epi_tnk", r148b["definition"], [r148b, r422b])])

    # Sensitivity: given TISCH + author DRMref (taken as given)
    t_lock = tisch[tisch.eligible_locked]
    r_tisch = _rho_row(
        "given_tisch_drmref",
        "GSE148071",
        "TISCH locked eligible ≥20/20; CLDN4_epi_mean vs n_TNK/n_cells",
        t_lock.cldn4,
        t_lock.frac_tnk,
        {"n_deposited": 42, "n_eligible": int(t_lock.shape[0]), "note": "given Q4 extract; 3 epithelial-like inside 25"},
    )
    t_mal = tisch[tisch.malignant_only]
    r_tisch_mal = _rho_row(
        "given_tisch_malignant",
        "GSE148071",
        "TISCH eligible and Malignant only",
        t_mal.cldn4,
        t_mal.frac_tnk,
        {"n_eligible": int(t_mal.shape[0])},
    )
    r_drm = _rho_row(
        "given_tisch_drmref",
        "GSE207422",
        "author DRMref malignant mean vs frac_tnk; n=12 post",
        drm.cldn4,
        drm.frac_tnk,
        {"n_deposited": 15, "n_eligible": int(drm.shape[0]), "note": "given combo table; author malignant"},
    )
    rows.extend([r_tisch, r_tisch_mal, r_drm])
    rows.append(
        _combo_row(
            "given_tisch_drmref",
            "TISCH locked + author DRMref (given extracts)",
            [r_tisch, r_drm],
            {"note": "sensitivity; mixed malignant definition"},
        )
    )

    rho_df = pd.DataFrame(rows)
    return rho_df, {
        "primary": combo_primary,
        "gse148071_cellchat": r148,
        "gse207422_cellchat": r422,
        "gse148071_tisch": r_tisch,
        "gse207422_drmref": r_drm,
        "n_148071_eligible": int(p148.shape[0]),
        "n_207422_eligible": int(p422.shape[0]),
    }


def build_lr() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    a = pd.read_csv(DATA / "gse148071_cellchat_ligand_table.tsv", sep="\t")
    b = pd.read_csv(DATA / "gse207422_cellchat_ligand_table.tsv", sep="\t")
    a["dataset"] = "GSE148071"
    b["dataset"] = "GSE207422"
    a["kept_split"] = "median (eligible 25)"
    b["kept_split"] = "median_post (12 post)"
    both = pd.concat([a, b], ignore_index=True)

    keys = ["interaction_name", "direction"]
    merged = a.merge(b, on=keys, how="outer", suffixes=("_148071", "_207422"))
    merged["in_148071"] = merged["delta_prob_148071"].notna()
    merged["in_207422"] = merged["delta_prob_207422"].notna()
    merged["same_sign"] = np.sign(merged["delta_prob_148071"].fillna(0)) == np.sign(
        merged["delta_prob_207422"].fillna(0)
    )
    merged["both_sig_same_dir"] = (
        merged["in_148071"]
        & merged["in_207422"]
        & merged["same_sign"]
        & (merged["delta_prob_148071"].fillna(0) != 0)
        & (merged["delta_prob_207422"].fillna(0) != 0)
    )
    # Prefer 207422 annotation columns when 148071 missing
    for col in ("pathway_name", "annotation", "ligand", "receptor", "ligand_class", "contrast"):
        left, right = f"{col}_148071", f"{col}_207422"
        if left in merged.columns and right in merged.columns:
            merged[col] = merged[left].combine_first(merged[right])
    merged["mean_delta"] = merged[["delta_prob_148071", "delta_prob_207422"]].mean(axis=1)
    merged["min_abs_delta"] = merged[["delta_prob_148071", "delta_prob_207422"]].abs().min(axis=1)
    merged = merged.sort_values(
        ["both_sig_same_dir", "direction", "mean_delta"],
        ascending=[False, True, False],
    )
    return both, merged, a, b


def honest_n(patients: pd.DataFrame, tisch: pd.DataFrame, drm: pd.DataFrame, lr_merged: pd.DataFrame) -> pd.DataFrame:
    p148 = patients[patients.dataset == "GSE148071"]
    p422 = patients[patients.dataset == "GSE207422"]
    rows = [
        {"item": "GSE148071_deposited", "n": 42, "note": "Wu 2021; one biopsy each; not the test n"},
        {"item": "GSE148071_cellchat_eligible_25_25", "n": int(p148.eligible.sum()), "note": "≥25 epi and ≥25 T/NK; marker-argmax"},
        {"item": "GSE148071_cellchat_excluded", "n": int((~p148.eligible).sum()), "note": "mostly epithelium with almost no T/NK"},
        {"item": "GSE148071_TISCH_locked_eligible", "n": int(tisch.eligible_locked.sum()), "note": "given ≥20/20; 3 epithelial-like inside"},
        {"item": "GSE207422_deposited_samples", "n": 15, "note": "Hu 2023; 3 pre + 12 post"},
        {"item": "GSE207422_post_used", "n": 12, "note": "pre-treatment excluded in given CellChat"},
        {"item": "GSE207422_cellchat_eligible_25_25", "n": int(p422.eligible.sum()), "note": "all 12 post pass ≥25/25"},
        {"item": "GSE207422_drmref_patients", "n": int(drm.shape[0]), "note": "author malignant; given combo table"},
        {
            "item": "combo_primary_n_patients",
            "n": int(p148.eligible.sum() + p422.eligible.sum()),
            "note": "25 + 12; do not write 42+15",
        },
        {"item": "GSE148071_sig_LR_pairs", "n": 52, "note": "kept median; 39 out / 13 in"},
        {"item": "GSE207422_sig_LR_pairs", "n": 100, "note": "kept median_post; 68 out / 32 in"},
        {
            "item": "consensus_both_sig_same_dir",
            "n": int(lr_merged.both_sig_same_dir.sum()),
            "note": "intersection of kept ligand tables",
        },
        {
            "item": "consensus_outgoing_both_sig_same_dir",
            "n": int(((lr_merged.direction == "outgoing") & lr_merged.both_sig_same_dir).sum()),
            "note": "Mal CLDN4-high → T/NK",
        },
    ]
    return pd.DataFrame(rows)


def plot_scatter(patients: pd.DataFrame, rho_info: dict, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), sharey=False)
    colors = {"GSE148071": "#2a6f97", "GSE207422": "#c44536"}
    for ax, ds, rec in (
        (axes[0], "GSE148071", rho_info["gse148071_cellchat"]),
        (axes[1], "GSE207422", rho_info["gse207422_cellchat"]),
    ):
        d = patients[(patients.dataset == ds) & patients.eligible]
        ax.scatter(d.cldn4, d.frac_tnk, s=36, c=colors[ds], edgecolors="white", linewidths=0.4, zorder=3)
        ax.set_title(
            f"{ds}\nn={rec['n']}  ρ={rec['rho']:+.3f}  p={rec['p']:.3g}",
            fontsize=10,
        )
        ax.set_xlabel("Malignant/epithelial CLDN4 mean (log1p CP10k)")
        ax.set_ylabel("T/NK fraction (of all cells)")
        ax.axhline(d.frac_tnk.median(), color="#999", lw=0.6, ls="--")
        ax.axvline(d.cldn4.median(), color="#999", lw=0.6, ls="--")
        for _, r in d.iterrows():
            if r.n_epithelial >= 2000 or r.n_TNK >= 1000 or r.frac_tnk > 0.4:
                ax.annotate(r.patient, (r.cldn4, r.frac_tnk), fontsize=6, alpha=0.75)
    combo = rho_info["primary"]
    fig.suptitle(
        f"Patient-level CLDN4 vs T/NK  ·  combo DL ρ={combo['rho']:+.3f} "
        f"(p={combo['p']:.3g}, I²={combo['I2']:.0f}%, n={combo['n']})",
        fontsize=11,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_forest(rho_df: pd.DataFrame, path: Path) -> None:
    show = rho_df[rho_df.family.isin(["primary_cellchat", "given_tisch_drmref"])].copy()
    show = show.reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8.4, 3.6))
    y = np.arange(len(show))[::-1]
    for i, rec in show.iterrows():
        color = "#1b4f72" if rec.family == "primary_cellchat" else "#7f8c8d"
        if rec.k > 1:
            color = "#6c3483" if rec.family == "primary_cellchat" else "#616a6b"
        ax.plot([rec.ci95_lo, rec.ci95_hi], [y[i], y[i]], color=color, lw=1.6)
        ax.plot(rec.rho, y[i], "o", color=color, ms=7 if rec.k > 1 else 5)
        tag = "primary" if rec.family == "primary_cellchat" else "given TISCH/DRMref"
        label = f"{tag} · {rec.cohort}  n={int(rec.n)}  ρ={rec.rho:+.3f}  p={rec.p:.3g}"
        ax.text(0.98, y[i], label, transform=ax.get_yaxis_transform(), va="center", ha="right", fontsize=7)
    ax.axvline(0, color="#333", lw=0.7)
    ax.set_yticks([])
    ax.set_xlabel("Spearman ρ (CLDN4 vs T/NK)  ·  combo = DL RE on Fisher-z")
    ax.set_title("Honest combo rho — include GSE207422; CLDN4 only")
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_n(patients: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.4))
    for ax, ds, title in (
        (axes[0], "GSE148071", "GSE148071 (42 biopsies; 25 eligible)"),
        (axes[1], "GSE207422", "GSE207422 (12 post; all eligible)"),
    ):
        d = patients[patients.dataset == ds].sort_values("n_epithelial", ascending=False)
        x = np.arange(len(d))
        ax.bar(x - 0.18, d.n_epithelial, 0.36, label="epithelial", color="#d4a373")
        ax.bar(x + 0.18, d.n_TNK, 0.36, label="T/NK", color="#2a6f97")
        for i, rec in enumerate(d.itertuples()):
            if not rec.eligible:
                ax.axvline(i, color="#c0392b", lw=0.4, alpha=0.35)
        ax.set_xticks(x)
        ax.set_xticklabels(d.patient, rotation=90, fontsize=6)
        ax.set_ylabel("cells")
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7, loc="upper right")
        ax.set_yscale("log")
    fig.suptitle("Honest n — red ticks = below 25/25 floor (not scored for rho or CellChat)", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_consensus_lr(merged: pd.DataFrame, path: Path) -> None:
    d = merged[(merged.direction == "outgoing") & merged.both_sig_same_dir].copy()
    d = d.sort_values("mean_delta", ascending=False)
    if d.empty:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.text(0.5, 0.5, "No consensus outgoing pairs", ha="center")
        fig.savefig(path, dpi=140)
        plt.close(fig)
        return
    fig, ax = plt.subplots(figsize=(8.8, max(3.2, 0.28 * len(d) + 1.4)))
    y = np.arange(len(d))[::-1]
    ax.barh(y + 0.16, d.delta_prob_148071, 0.30, label="GSE148071 ΔP", color="#2a6f97")
    ax.barh(y - 0.16, d.delta_prob_207422, 0.30, label="GSE207422 ΔP", color="#c44536")
    ax.axvline(0, color="#333", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(d.interaction_name.str.replace("_", "–"), fontsize=7.5)
    ax.set_xlabel("ΔP (CLDN4-high − CLDN4-low)  outgoing Mal → T/NK")
    ax.set_title(f"Consensus outgoing pairs (sig in both, same sign)  n={len(d)}")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_key_pairs(merged: pd.DataFrame, path: Path) -> None:
    rows = []
    for name in KEY_PAIRS:
        sub = merged[merged.interaction_name == name]
        if sub.empty:
            rows.append({"pair": name, "d148": np.nan, "d422": np.nan, "dir": ""})
            continue
        # prefer outgoing if both exist
        hit = sub[sub.direction == "outgoing"]
        rec = (hit if not hit.empty else sub).iloc[0]
        rows.append(
            {
                "pair": name,
                "d148": rec.get("delta_prob_148071", np.nan),
                "d422": rec.get("delta_prob_207422", np.nan),
                "dir": rec.get("direction", ""),
            }
        )
    d = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    y = np.arange(len(d))[::-1]
    ax.barh(y + 0.16, d.d148, 0.30, label="GSE148071", color="#2a6f97")
    ax.barh(y - 0.16, d.d422, 0.30, label="GSE207422", color="#c44536")
    ax.axvline(0, color="#333", lw=0.7)
    labels = [f"{p.replace('_', '–')} ({r})" if r else p.replace("_", "–") for p, r in zip(d.pair, d.dir)]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("ΔP (high − low); missing bar = not in that kept ligand table")
    ax.set_title("Extra figure — key CLDN4 pairs (not dual-high)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_ligand_table(merged: pd.DataFrame, path: Path) -> None:
    d = merged[merged.both_sig_same_dir].copy()
    d = d.sort_values(["direction", "mean_delta"], ascending=[True, False])
    fig, ax = plt.subplots(figsize=(11.2, max(3.5, 0.22 * len(d) + 1.6)))
    ax.axis("off")
    cols = [
        "interaction_name",
        "direction",
        "pathway_name",
        "delta_prob_148071",
        "delta_prob_207422",
        "mean_delta",
    ]
    tab = d[cols].copy()
    tab["interaction_name"] = tab["interaction_name"].str.replace("_", "–")
    tab["delta_prob_148071"] = tab["delta_prob_148071"].map(lambda x: f"{x:+.3f}")
    tab["delta_prob_207422"] = tab["delta_prob_207422"].map(lambda x: f"{x:+.3f}")
    tab["mean_delta"] = tab["mean_delta"].map(lambda x: f"{x:+.3f}")
    table = ax.table(
        cellText=tab.values,
        colLabels=["pair", "dir", "pathway", "ΔP 148071", "ΔP 207422", "mean ΔP"],
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1, 1.15)
    ax.set_title("Extra figure — consensus LR table (both sig, same direction)", pad=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    patients = cellchat_patients()
    tisch, drm = given_patients()
    rho_df, rho_info = build_rho(patients, tisch, drm)
    both_lr, merged, lig148, lig422 = build_lr()
    n_df = honest_n(patients, tisch, drm, merged)

    patients.to_csv(OUT / "patient_scores.tsv", sep="\t", index=False)
    rho_df.to_csv(OUT / "combo_rho.tsv", sep="\t", index=False)
    n_df.to_csv(OUT / "honest_n.tsv", sep="\t", index=False)
    both_lr.to_csv(OUT / "lr_pairs_by_dataset.tsv", sep="\t", index=False)

    consensus = merged[merged.both_sig_same_dir].copy()
    keep_cols = [
        c
        for c in [
            "interaction_name",
            "direction",
            "pathway_name",
            "annotation",
            "ligand",
            "receptor",
            "ligand_class",
            "delta_prob_148071",
            "delta_prob_207422",
            "mean_delta",
            "min_abs_delta",
            "prob_high_148071",
            "prob_low_148071",
            "prob_high_207422",
            "prob_low_207422",
            "pval_high_148071",
            "pval_high_207422",
            "kept_split_148071",
            "kept_split_207422",
            "both_sig_same_dir",
        ]
        if c in consensus.columns or c in merged.columns
    ]
    # some kept_split cols may be missing after merge
    for c in ("kept_split_148071", "kept_split_207422"):
        if c not in merged.columns:
            pass
    out_cols = [c for c in keep_cols if c in merged.columns]
    consensus[out_cols].to_csv(OUT / "lr_table.tsv", sep="\t", index=False)
    merged.to_csv(OUT / "lr_join.tsv", sep="\t", index=False)
    lig148.to_csv(OUT / "lr_table_gse148071.tsv", sep="\t", index=False)
    lig422.to_csv(OUT / "lr_table_gse207422.tsv", sep="\t", index=False)

    outgoing = consensus[consensus.direction == "outgoing"].sort_values("mean_delta", ascending=False)
    outgoing.to_csv(OUT / "lr_table_outgoing.tsv", sep="\t", index=False)

    plot_scatter(patients, rho_info, FIG / "fig_combo_rho_scatter.png")
    plot_forest(rho_df, FIG / "fig_combo_rho_forest.png")
    plot_n(patients, FIG / "fig_n_honest.png")
    plot_consensus_lr(merged, FIG / "fig_extra_consensus_outgoing.png")
    plot_key_pairs(merged, FIG / "fig_extra_key_pairs.png")
    plot_ligand_table(merged, FIG / "fig_extra_ligand_table.png")
    # also copy extra figures next to results for the PR convention used in single-dataset folders
    for src in FIG.glob("*.png"):
        dest = OUT / src.name
        dest.write_bytes(src.read_bytes())

    summary = {
        "additive": True,
        "marker": "CLDN4",
        "dual_high": False,
        "include_207422": True,
        "datasets": ["GSE207422", "GSE148071"],
        "not_used": ["TACSTD2 gate", "dual-high", "CellChat R", "LIANA", "raw FASTQ"],
        "primary_combo_rho": rho_info["primary"],
        "gse148071_cellchat_rho": rho_info["gse148071_cellchat"],
        "gse207422_cellchat_rho": rho_info["gse207422_cellchat"],
        "gse148071_tisch_rho": rho_info["gse148071_tisch"],
        "gse207422_drmref_rho": rho_info["gse207422_drmref"],
        "n_148071_eligible": rho_info["n_148071_eligible"],
        "n_207422_eligible": rho_info["n_207422_eligible"],
        "n_combo_primary": int(rho_info["n_148071_eligible"] + rho_info["n_207422_eligible"]),
        "n_consensus_lr": int(consensus.shape[0]),
        "n_consensus_outgoing": int(outgoing.shape[0]),
        "n_consensus_outgoing_high": int((outgoing.mean_delta > 0).sum()),
        "n_consensus_outgoing_low": int((outgoing.mean_delta < 0).sum()),
        "consensus_outgoing_pairs": outgoing.interaction_name.tolist(),
        "honest_n": n_df.to_dict(orient="records"),
        "algorithm": (
            "Patient Spearman on given CellChat per-sample tables; "
            "combo = DerSimonian–Laird RE on Fisher-z(Spearman). "
            "LR tables taken as given from PR #348 and PR #324; "
            "consensus = significant in both kept splits, same ΔP sign."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps({k: summary[k] for k in (
        "primary_combo_rho",
        "gse148071_cellchat_rho",
        "gse207422_cellchat_rho",
        "n_combo_primary",
        "n_consensus_lr",
        "n_consensus_outgoing",
        "consensus_outgoing_pairs",
    )}, indent=2, default=str))


if __name__ == "__main__":
    main()
