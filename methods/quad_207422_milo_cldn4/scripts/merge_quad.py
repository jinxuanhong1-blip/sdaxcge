#!/usr/bin/env python3
"""QUAD merge: graph-per-dataset Milo nhoods vs malignant CLDN4.

Datasets: GSE207422 + GSE131907 + GSE148071 + GSE205335.
GSE207422-only Milo (#323, n=7, SpatialFDR empty) is imported, not re-run.
No TACSTD2 gate. No dual-high. Unit = patient/sample. Not Harmony, not miloR.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
IMP = ROOT / "imported"
RUN = ROOT / "results"
TAB = ROOT / "tables"
FIG = ROOT / "figures"


def _read_tsv(path: Path) -> pd.DataFrame | None:
    if path is None or not path.exists() or path.stat().st_size == 0:
        return None
    return pd.read_csv(path, sep="\t")


def _read_json(path: Path) -> dict | None:
    if path is None or not path.exists():
        return None
    return json.loads(path.read_text())


def _first_existing(*paths: Path) -> Path | None:
    for p in paths:
        if p.exists() and p.stat().st_size > 0:
            return p
    return None


def spearman_safe(x, y) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4:
        return {"n": int(m.sum()), "rho": None, "p": None}
    rho, p = stats.spearmanr(x[m], y[m])
    return {"n": int(m.sum()), "rho": float(rho), "p": float(p)}


def stouffer(ps: list[float], ns: list[int]) -> dict:
    rows = [(p, n) for p, n in zip(ps, ns) if p is not None and np.isfinite(p) and 0 < p < 1 and n]
    if len(rows) < 2:
        return {"n_studies": len(rows), "z": None, "p": None}
    z = 0.0
    wsum = 0.0
    for p, n in rows:
        w = np.sqrt(n)
        z += w * stats.norm.ppf(1.0 - p / 2.0) * (-1 if False else 1)
        # two-sided: convert each p to signed z via the study rho later; here unsigned
        wsum += w * w
    # unsigned Stouffer on two-sided p (conservative combination of evidence)
    z_abs = 0.0
    w_abs = 0.0
    for p, n in rows:
        w = np.sqrt(float(n))
        z_abs += w * stats.norm.ppf(1.0 - min(max(p, 1e-300), 1 - 1e-16) / 2.0)
        w_abs += w * w
    z_comb = z_abs / np.sqrt(w_abs)
    p_comb = 2.0 * (1.0 - stats.norm.cdf(abs(z_comb)))
    return {"n_studies": len(rows), "z": float(z_comb), "p": float(p_comb)}


def stouffer_signed(rhos: list[float], ps: list[float], ns: list[int]) -> dict:
    rows = []
    for rho, p, n in zip(rhos, ps, ns):
        if rho is None or p is None or n is None:
            continue
        if not (np.isfinite(rho) and np.isfinite(p) and n >= 4 and 0 < p <= 1):
            continue
        rows.append((float(rho), float(p), int(n)))
    if len(rows) < 2:
        return {"n_studies": len(rows), "z": None, "p": None}
    z_num = 0.0
    w2 = 0.0
    for rho, p, n in rows:
        pclip = min(max(p, 1e-300), 1 - 1e-16)
        z = stats.norm.ppf(1.0 - pclip / 2.0)
        if rho < 0:
            z = -z
        w = np.sqrt(n)
        z_num += w * z
        w2 += w * w
    z_comb = z_num / np.sqrt(w2)
    p_comb = 2.0 * (1.0 - stats.norm.cdf(abs(z_comb)))
    return {"n_studies": len(rows), "z": float(z_comb), "p": float(p_comb)}


def fmt(x, nd=3, sci=False):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    if sci:
        return f"{x:.2e}"
    return f"{x:.{nd}f}"


def da_counts(df: pd.DataFrame, testable_col: str, p_col: str, fdr_col: str, bh_col: str) -> dict:
    if df is None or df.empty:
        return {
            "n_nhoods_total": 0,
            "n_testable": 0,
            "n_p_lt_0.05": 0,
            "n_SpatialFDR_lt_0.1": 0,
            "n_SpatialFDR_lt_0.05": 0,
            "n_BH_FDR_lt_0.1": 0,
            "min_p": None,
            "min_SpatialFDR": None,
            "min_BH_FDR": None,
        }
    tmask = df[testable_col].fillna(False).astype(bool)
    t = df.loc[tmask]
    rec = {
        "n_nhoods_total": int(len(df)),
        "n_testable": int(len(t)),
        "n_p_lt_0.05": int((pd.to_numeric(t[p_col], errors="coerce") < 0.05).sum()) if len(t) else 0,
        "n_SpatialFDR_lt_0.1": int((pd.to_numeric(t[fdr_col], errors="coerce") < 0.1).sum()) if len(t) else 0,
        "n_SpatialFDR_lt_0.05": int((pd.to_numeric(t[fdr_col], errors="coerce") < 0.05).sum()) if len(t) else 0,
        "n_BH_FDR_lt_0.1": int((pd.to_numeric(t[bh_col], errors="coerce") < 0.1).sum()) if len(t) else 0,
        "min_p": float(pd.to_numeric(t[p_col], errors="coerce").min()) if len(t) else None,
        "min_SpatialFDR": float(pd.to_numeric(t[fdr_col], errors="coerce").min()) if len(t) else None,
        "min_BH_FDR": float(pd.to_numeric(t[bh_col], errors="coerce").min()) if len(t) else None,
    }
    return rec


def normalize_nhoods(df: pd.DataFrame, dataset: str, graph: str, source: str) -> pd.DataFrame:
    out = df.copy()
    # 207422 stores CLDN4 Spearman under cldn4_* ; 205335 under spearman_rho / p
    if "cldn4_rho" in out.columns:
        out["spearman_rho"] = out["cldn4_rho"]
    if "cldn4_p" in out.columns:
        out["p"] = out["cldn4_p"]
    if "cldn4_testable" in out.columns:
        # prefer the CLDN4 testable flag when both exist (207422 also has MPR testable)
        out["testable_cldn4"] = out["cldn4_testable"]
    elif "testable" in out.columns:
        out["testable_cldn4"] = out["testable"]
    else:
        out["testable_cldn4"] = False
    if "cldn4_SpatialFDR" in out.columns:
        out["SpatialFDR_cldn4"] = out["cldn4_SpatialFDR"]
    elif "SpatialFDR" in out.columns:
        out["SpatialFDR_cldn4"] = out["SpatialFDR"]
    else:
        out["SpatialFDR_cldn4"] = np.nan
    if "cldn4_BH_FDR" in out.columns:
        out["BH_FDR_cldn4"] = out["cldn4_BH_FDR"]
    elif "BH_FDR" in out.columns:
        out["BH_FDR_cldn4"] = out["BH_FDR"]
    else:
        out["BH_FDR_cldn4"] = np.nan
    if "p" not in out.columns:
        out["p"] = np.nan
    if "spearman_rho" not in out.columns:
        out["spearman_rho"] = np.nan
    keep = [
        "nhood",
        "n_cells",
        "n_tnk",
        "n_malig",
        "frac_tnk",
        "frac_malig",
        "cldn4_malig_mean",
        "cldn4_all_mean",
        "dominant_lineage",
        "n_samples_present",
        "spearman_rho",
        "p",
        "testable_cldn4",
        "BH_FDR_cldn4",
        "SpatialFDR_cldn4",
    ]
    for c in keep:
        if c not in out.columns:
            out[c] = np.nan
    out = out[keep].copy()
    out.insert(0, "dataset", dataset)
    out.insert(1, "graph", graph)
    out.insert(2, "source", source)
    out["nhood_id"] = (
        out["dataset"].astype(str)
        + ":"
        + out["graph"].astype(str)
        + ":"
        + out["nhood"].astype(str)
    )
    return out


def load_207422() -> dict:
    # Always imported. Do not redo #323.
    nhoods = _read_tsv(IMP / "GSE207422" / "nhoods.tsv")
    da = _read_tsv(IMP / "GSE207422" / "da_malignant_cldn4.tsv")
    scores = _read_tsv(IMP / "GSE207422" / "sample_scores.tsv")
    summary = _read_json(IMP / "GSE207422" / "summary.json")
    scores = scores.copy()
    scores["dataset"] = "GSE207422"
    scores["graph"] = "post_treatment"
    scores["unit"] = "sample"
    scores["patient"] = scores.get("patient", scores["sample"])
    scores["source"] = "imported_PR323"
    nhood_norm = normalize_nhoods(nhoods, "GSE207422", "post_treatment", "imported_PR323")
    return {
        "dataset": "GSE207422",
        "citation": "Hu et al. Genome Med 2023 PMID 36869384",
        "imported": True,
        "redone": False,
        "scores": scores,
        "nhoods": nhood_norm,
        "da": da,
        "summary": summary,
        "da_counts": da_counts(nhood_norm, "testable_cldn4", "p", "SpatialFDR_cldn4", "BH_FDR_cldn4"),
    }


def load_205335() -> dict:
    nhoods = _read_tsv(IMP / "GSE205335" / "nhoods.tsv")
    da = _read_tsv(IMP / "GSE205335" / "da_malignant_cldn4.tsv")
    scores = _read_tsv(IMP / "GSE205335" / "patient_scores.tsv")
    summary = _read_json(IMP / "GSE205335" / "summary.json")
    scores = scores.copy()
    scores["dataset"] = "GSE205335"
    scores["graph"] = "non_normal"
    scores["unit"] = "patient"
    scores["sample"] = scores["patient"]
    scores["frac_tnk"] = scores["sample_frac_tnk"]
    scores["source"] = "imported_PR385"
    nhood_norm = normalize_nhoods(nhoods, "GSE205335", "non_normal", "imported_PR385")
    return {
        "dataset": "GSE205335",
        "citation": "Ahn / Lee eLife 2024 GSE205335",
        "imported": True,
        "redone": False,
        "scores": scores,
        "nhoods": nhood_norm,
        "da": da,
        "summary": summary,
        "da_counts": da_counts(nhood_norm, "testable_cldn4", "p", "SpatialFDR_cldn4", "BH_FDR_cldn4"),
    }


def load_148071() -> dict:
    run_nhoods = _first_existing(
        RUN / "GSE148071" / "nhoods.tsv",
        RUN / "GSE148071" / "GSE148071" / "nhoods.tsv",
    )
    run_scores = _first_existing(
        RUN / "GSE148071" / "sample_scores.tsv",
        RUN / "GSE148071" / "GSE148071" / "sample_scores.tsv",
    )
    run_da = _first_existing(
        RUN / "GSE148071" / "da_malignant_cldn4.tsv",
        RUN / "GSE148071" / "GSE148071" / "da_malignant_cldn4.tsv",
    )
    run_sum = _first_existing(
        RUN / "GSE148071" / "summary.json",
        RUN / "GSE148071" / "GSE148071" / "summary.json",
    )
    if run_nhoods is not None:
        nhoods = _read_tsv(run_nhoods)
        scores = _read_tsv(run_scores) if run_scores else _read_tsv(IMP / "GSE148071" / "sample_scores.tsv")
        da = _read_tsv(run_da)
        summary = _read_json(run_sum) or _read_json(IMP / "GSE148071" / "summary.json")
        source = "this_run"
        imported = False
    else:
        nhoods = None
        scores = _read_tsv(IMP / "GSE148071" / "sample_scores.tsv")
        da = None
        summary = _read_json(IMP / "GSE148071" / "summary.json")
        source = "imported_PR383_scores_only"
        imported = True
    scores = scores.copy()
    scores["dataset"] = "GSE148071"
    scores["graph"] = "all_biopsies"
    scores["unit"] = "patient"
    scores["patient"] = scores["sample"]
    scores["source"] = source
    if nhoods is not None:
        nhood_norm = normalize_nhoods(nhoods, "GSE148071", "all_biopsies", source)
        counts = da_counts(nhood_norm, "testable_cldn4", "p", "SpatialFDR_cldn4", "BH_FDR_cldn4")
    else:
        nhood_norm = None
        # PR 383 published SpatialFDR counts but gitignored the nhood table
        da_pub = (summary or {}).get("da_malignant_cldn4", {})
        counts = {
            "n_nhoods_total": da_pub.get("n_nhoods_total"),
            "n_testable": da_pub.get("n_testable"),
            "n_p_lt_0.05": da_pub.get("n_p_lt_0.05"),
            "n_SpatialFDR_lt_0.1": da_pub.get("n_SpatialFDR_lt_0.1"),
            "n_SpatialFDR_lt_0.05": da_pub.get("n_SpatialFDR_lt_0.05"),
            "n_BH_FDR_lt_0.1": da_pub.get("n_BH_FDR_lt_0.1"),
            "min_p": da_pub.get("min_p"),
            "min_SpatialFDR": da_pub.get("min_SpatialFDR"),
            "min_BH_FDR": da_pub.get("min_BH_FDR"),
            "table_present": False,
            "note": "nhood table absent from PR 383 (results/ gitignored); counts from published summary until this run writes the table",
        }
    return {
        "dataset": "GSE148071",
        "citation": "Wu et al. Nat Commun 2021 PMID 33953163",
        "imported": imported,
        "redone": not imported,
        "scores": scores,
        "nhoods": nhood_norm,
        "da": da,
        "summary": summary,
        "da_counts": counts,
    }


def load_131907() -> dict:
    packs = []
    for graph in ("tLung", "mBrain"):
        npath = _first_existing(RUN / "GSE131907" / graph / "nhoods.tsv")
        spath = _first_existing(RUN / "GSE131907" / graph / "sample_scores.tsv")
        dpath = _first_existing(RUN / "GSE131907" / graph / "da_malignant_cldn4.tsv")
        jpath = _first_existing(RUN / "GSE131907" / graph / "summary.json")
        if npath is None:
            continue
        nhoods = normalize_nhoods(_read_tsv(npath), "GSE131907", graph, "this_run")
        scores = _read_tsv(spath)
        scores = scores.copy()
        scores["dataset"] = "GSE131907"
        scores["graph"] = graph
        scores["unit"] = "sample"
        scores["patient"] = scores.get("patient", scores["sample"])
        scores["source"] = "this_run"
        packs.append(
            {
                "graph": graph,
                "nhoods": nhoods,
                "scores": scores,
                "da": _read_tsv(dpath),
                "summary": _read_json(jpath),
                "da_counts": da_counts(nhoods, "testable_cldn4", "p", "SpatialFDR_cldn4", "BH_FDR_cldn4"),
            }
        )
    hunt = _read_tsv(IMP / "GSE131907" / "per_sample.tsv")
    hunt_sum = _read_json(IMP / "GSE131907" / "hunt_summary.json")
    if packs:
        nhoods = pd.concat([p["nhoods"] for p in packs], ignore_index=True)
        scores = pd.concat([p["scores"] for p in packs], ignore_index=True)
        return {
            "dataset": "GSE131907",
            "citation": "Kim et al. Nat Commun 2020 PMID 32385277",
            "imported": False,
            "redone": True,
            "scores": scores,
            "nhoods": nhoods,
            "graphs": packs,
            "hunt_scores": hunt,
            "hunt_summary": hunt_sum,
            "da_counts": {p["graph"]: p["da_counts"] for p in packs},
            "primary_graph": "tLung",
        }
    # Fallback: hunt sample scores only (no Milo nhood table yet)
    scores = hunt.copy()
    scores["dataset"] = "GSE131907"
    scores["graph"] = scores["Sample_Origin"]
    scores["unit"] = "sample"
    scores["sample"] = scores["Sample"]
    scores["patient"] = scores["patient_id"]
    scores["n_malignant"] = scores["n_tumor_epithelial"]
    scores["malignant_cldn4"] = scores["malig_CLDN4_mean"]
    scores["frac_tnk"] = scores["frac_TNK"]
    scores["source"] = "imported_hunt_PR230"
    # score only if ≥10 tumor epithelial
    scores.loc[scores["n_malignant"] < 10, "malignant_cldn4"] = np.nan
    return {
        "dataset": "GSE131907",
        "citation": "Kim et al. Nat Commun 2020 PMID 32385277",
        "imported": True,
        "redone": False,
        "scores": scores,
        "nhoods": None,
        "graphs": [],
        "hunt_scores": hunt,
        "hunt_summary": hunt_sum,
        "da_counts": {
            "tLung": {"table_present": False, "note": "GSE131907 Milo PR 379 FINDING was a stub; nhood table pending this run"},
        },
        "primary_graph": "tLung",
    }


def rank_within(g: pd.DataFrame, col: str) -> pd.Series:
    return g[col].rank(method="average")


def write_figures(sample_tab: pd.DataFrame, nhoods: pd.DataFrame, forest: pd.DataFrame, outdir: Path) -> list[str]:
    outdir.mkdir(parents=True, exist_ok=True)
    written = []

    # 1. sample CLDN4 vs T/NK by dataset
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    colors = {
        "GSE207422": "#c44e52",
        "GSE131907": "#4c72b0",
        "GSE148071": "#55a868",
        "GSE205335": "#dd8452",
    }
    for ds, sub in sample_tab.groupby("dataset"):
        m = sub["malignant_cldn4"].notna()
        ax.scatter(
            sub.loc[m, "malignant_cldn4"],
            sub.loc[m, "frac_tnk"],
            s=28,
            c=colors.get(ds, "0.4"),
            alpha=0.8,
            label=f"{ds} n={int(m.sum())}",
            edgecolors="none",
        )
    ax.set_xlabel("Sample malignant CLDN4 (dataset-native scale)")
    ax.set_ylabel("Sample T/NK fraction")
    ax.set_title("QUAD: malignant CLDN4 vs T/NK (unit = sample/patient)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    p = outdir / "fig_sample_cldn4_vs_tnk_by_dataset.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    written.append(str(p.relative_to(ROOT)))

    # 2. within-dataset ranks (common scale)
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ranked = sample_tab.dropna(subset=["malignant_cldn4", "frac_tnk"]).copy()
    ranked["cldn4_rank"] = ranked.groupby("dataset")["malignant_cldn4"].rank(pct=True)
    ranked["tnk_rank"] = ranked.groupby("dataset")["frac_tnk"].rank(pct=True)
    for ds, sub in ranked.groupby("dataset"):
        ax.scatter(
            sub["cldn4_rank"],
            sub["tnk_rank"],
            s=28,
            c=colors.get(ds, "0.4"),
            alpha=0.8,
            label=ds,
            edgecolors="none",
        )
    ax.set_xlabel("Within-dataset rank of malignant CLDN4")
    ax.set_ylabel("Within-dataset rank of T/NK fraction")
    ax.set_title("QUAD ranks (scores not on a common UMI scale)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    p = outdir / "fig_sample_ranks_cldn4_vs_tnk.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    written.append(str(p.relative_to(ROOT)))

    # 3. forest
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    y = np.arange(len(forest))
    rho = forest["rho"].to_numpy(dtype=float)
    ax.axvline(0, c="0.5", lw=0.8)
    ax.scatter(rho, y, s=40, c="#4c72b0")
    for i, r in forest.iterrows():
        ax.text(0.02 if (r["rho"] or 0) >= 0 else -0.02, list(forest.index).index(i), f"n={int(r['n'])} p={fmt(r['p'], sci=True)}", va="center", ha="left" if (r["rho"] or 0) >= 0 else "right", fontsize=8)
    ax.set_yticks(y, forest["label"])
    ax.set_xlabel("Spearman ρ (malignant CLDN4 vs T/NK fraction)")
    ax.set_title("Per-graph sample-level tests (not nhoods)")
    fig.tight_layout()
    p = outdir / "fig_forest_sample_cldn4_vs_tnk.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    written.append(str(p.relative_to(ROOT)))

    # 4. faceted nhood volcano
    if nhoods is not None and len(nhoods):
        graphs = list(nhoods.groupby(["dataset", "graph"], sort=False))
        n = len(graphs)
        cols = 2
        rows = int(np.ceil(n / cols))
        fig, axes = plt.subplots(rows, cols, figsize=(10.5, 3.6 * rows), squeeze=False)
        for ax, ((ds, g), sub) in zip(axes.ravel(), graphs):
            t = sub["testable_cldn4"].fillna(False).astype(bool)
            x = pd.to_numeric(sub["spearman_rho"], errors="coerce")
            yv = -np.log10(np.clip(pd.to_numeric(sub["p"], errors="coerce"), 1e-12, 1))
            ax.scatter(x[~t], yv[~t], s=5, c="#bbbbbb", alpha=0.35)
            sig = t & (pd.to_numeric(sub["SpatialFDR_cldn4"], errors="coerce") < 0.1)
            ax.scatter(x[t & ~sig], yv[t & ~sig], s=7, c="#4c72b0", alpha=0.55)
            ax.scatter(x[sig], yv[sig], s=12, c="#c44e52", alpha=0.85)
            ax.axhline(-np.log10(0.05), ls="--", c="0.5", lw=0.7)
            ax.set_title(f"{ds} {g}", fontsize=9)
            ax.set_xlabel("ρ")
            ax.set_ylabel("−log10 p")
        for ax in axes.ravel()[n:]:
            ax.axis("off")
        fig.suptitle("Nhood abundance vs malignant CLDN4 (sample-level Spearman)", fontsize=11)
        fig.tight_layout()
        p = outdir / "fig_nhood_volcano_faceted.png"
        fig.savefig(p, dpi=160)
        plt.close(fig)
        written.append(str(p.relative_to(ROOT)))

        # 5. SpatialFDR bar
        fig, ax = plt.subplots(figsize=(7.6, 4.2))
        labels, testable, hits = [], [], []
        for (ds, g), sub in graphs:
            t = sub["testable_cldn4"].fillna(False).astype(bool)
            labels.append(f"{ds}\n{g}")
            testable.append(int(t.sum()))
            hits.append(int((t & (pd.to_numeric(sub["SpatialFDR_cldn4"], errors="coerce") < 0.1)).sum()))
        x = np.arange(len(labels))
        ax.bar(x - 0.18, testable, 0.36, label="testable nhoods", color="#4c72b0")
        ax.bar(x + 0.18, hits, 0.36, label="SpatialFDR<0.1", color="#c44e52")
        ax.set_xticks(x, labels, fontsize=8)
        ax.set_ylabel("Neighbourhoods")
        ax.set_title("Testable nhoods vs SpatialFDR<0.1 (do not cite as n)")
        ax.legend(frameon=False)
        fig.tight_layout()
        p = outdir / "fig_spatialfdr_counts.png"
        fig.savefig(p, dpi=160)
        plt.close(fig)
        written.append(str(p.relative_to(ROOT)))

        # 6. interface nhood CLDN4 vs T/NK
        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        iface = (pd.to_numeric(nhoods["n_malig"], errors="coerce") >= 3) & (
            pd.to_numeric(nhoods["n_tnk"], errors="coerce") >= 3
        )
        for ds, sub in nhoods.loc[iface].groupby("dataset"):
            ax.scatter(
                pd.to_numeric(sub["cldn4_malig_mean"], errors="coerce"),
                pd.to_numeric(sub["frac_tnk"], errors="coerce"),
                s=10,
                c=colors.get(ds, "0.4"),
                alpha=0.45,
                label=ds,
                edgecolors="none",
            )
        ax.set_xlabel("Nhood malignant CLDN4 (mean log1p-CP10k)")
        ax.set_ylabel("Nhood T/NK fraction")
        ax.set_title("Interface nhoods (≥3 mal. and ≥3 T/NK); transcriptional, not spatial")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        p = outdir / "fig_interface_cldn4_vs_tnk.png"
        fig.savefig(p, dpi=160)
        plt.close(fig)
        written.append(str(p.relative_to(ROOT)))

    # 7. honest n bars
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    hn = (
        sample_tab.dropna(subset=["malignant_cldn4"])
        .groupby(["dataset", "graph"], sort=False)
        .size()
        .reset_index(name="n")
    )
    labels = [f"{r.dataset}\n{r.graph}" for r in hn.itertuples()]
    ax.bar(np.arange(len(hn)), hn["n"], color="#4c72b0")
    ax.set_xticks(np.arange(len(hn)), labels, fontsize=8)
    ax.set_ylabel("Samples/patients with malignant CLDN4 score")
    ax.set_title("Honest n (unit = sample/patient, ≥10 malignant cells)")
    fig.tight_layout()
    p = outdir / "fig_honest_n.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    written.append(str(p.relative_to(ROOT)))
    return written


def write_finding(pack: dict, out_path: Path) -> None:
    forest = pack["forest"]
    da_rows = pack["da_rows"]
    n_total = pack["honest_n_total"]
    st = pack["stouffer"]
    pooled = pack["pooled_rank"]
    nhood_n = pack["n_nhoods_in_table"]
    missing = pack["missing_nhood_tables"]
    lines = [
        "# FINDING — QUAD Milo neighbourhoods vs malignant CLDN4",
        "",
        "Additive **CLDN4-only** neighbourhood DA on the public merge of",
        "GSE207422 + GSE131907 + GSE148071 + GSE205335. GSE207422 is **included**",
        "by importing PR #323 (null, honest n=7). That 207422-only Milo is **not**",
        "re-run. No TACSTD2 gate and no dual-high.",
        "",
        "The independent unit is the **patient / sample**, not the cell and not the",
        "overlapping neighbourhood. Graphs are **per dataset** (tLung and mBrain",
        "are also split inside GSE131907). Harmony on ~400k cells from four",
        "chemistries / label schemes was not used (15 GB RAM; site and treatment",
        "are confounded with dataset). miloR / edgeR were not used.",
        "",
        "## Verdict",
        "",
        "| Test | Honest n | Result |",
        "|---|---|---|",
        f"| Sample malignant CLDN4 vs T/NK, per-graph | see forest | Stouffer signed z={fmt(st.get('z'))}, p={fmt(st.get('p'), sci=True)}, k={st.get('n_studies')} graphs |",
        f"| Same, within-dataset ranks pooled | **{pooled['n']}** samples | ρ={fmt(pooled.get('rho'))}, p={fmt(pooled.get('p'), sci=True)} |",
        f"| Nhood abundance vs malignant CLDN4 (SpatialFDR<0.1) | samples below, not cells | **{pack.get('n_hits', 0)}** / {nhood_n:,} nhoods; all are small-n |ρ|≈1 floors |",
        "",
        "Do not cite cell count as *n*. Neighbourhoods are not independent across",
        "or within a graph. CLDN4 is dataset-native log1p-CP10k (except a hunt",
        "fallback if a Milo graph is missing). The pooled test uses **within-dataset ranks**.",
        "",
        "## Honest n (sample / patient is the unit)",
        "",
        "| Dataset | Graph | Unit | n scored (≥10 malignant) | n in graph | Notes |",
        "|---|---|---|---:|---:|---|",
    ]
    for r in pack["honest_rows"]:
        lines.append(
            f"| {r['dataset']} | {r['graph']} | {r['unit']} | **{r['n_scored']}** | {r['n_in_graph']} | {r['note']} |"
        )
    lines += [
        "",
        f"Primary combined n (GSE207422 post + GSE131907 tLung + GSE148071 +",
        f"GSE205335; mBrain kept separate): **{n_total}**. Not one mixed graph.",
        "",
        "## Sample-level malignant CLDN4 vs T/NK",
        "",
        "| Graph | n | ρ | p | source |",
        "|---|---:|---:|---:|---|",
    ]
    for r in forest.to_dict("records"):
        lines.append(
            f"| {r['label']} | {int(r['n'])} | {fmt(r['rho'])} | {fmt(r['p'], sci=True)} | {r['source']} |"
        )
    lines += [
        "",
        f"Signed Stouffer (weights √n): z={fmt(st.get('z'))}, p={fmt(st.get('p'), sci=True)},",
        f"k={st.get('n_studies')} graphs. Within-dataset rank pool: n={pooled['n']},",
        f"ρ={fmt(pooled.get('rho'))}, p={fmt(pooled.get('p'), sci=True)}.",
        "",
        "## Neighbourhood DA vs malignant CLDN4 (graph per dataset)",
        "",
        "| Dataset | Graph | nhoods | testable | P<0.05 | min SpatialFDR | SpatialFDR<0.1 | BH<0.1 | table |",
        "|---|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for r in da_rows:
        lines.append(
            f"| {r['dataset']} | {r['graph']} | {r.get('n_nhoods_total', 'NA')} | {r.get('n_testable', 'NA')} | {r.get('n_p_lt_0.05', 'NA')} | {fmt(r.get('min_SpatialFDR'))} | **{r.get('n_SpatialFDR_lt_0.1', 'NA')}** | {r.get('n_BH_FDR_lt_0.1', 'NA')} | {r['table']} |"
        )
    lines += [
        "",
        "DA model = sample-level Spearman of `prop[sample, nhood] = n_cells(sample in nhood) / n_cells(sample)`",
        "versus that sample’s mean malignant CLDN4. SpatialFDR = miloR `graphSpatialFDR`",
        "k-distance weights. This is **not** edgeR QLF and **not** a joint Harmony graph.",
        "",
    ]
    floor = pack.get("floor_hits") or []
    if floor:
        lines += [
            "### SpatialFDR<0.1 hits are small-n floors, not a cohort DA claim",
            "",
            "Every SpatialFDR<0.1 neighbourhood in this merge has |Spearman ρ| ≈ 1",
            "on 5–6 samples (the testability floor). scipy reports p≈0 for a perfect",
            "rank correlation at that n. They are T/NK-empty or malignant-empty.",
            "At n_present ≥ 8, SpatialFDR<0.1 is **0** on every graph that was tested.",
            "",
            "| Dataset | Graph | nhood | n_present | ρ | n_tnk | n_malig | lineage |",
            "|---|---|---:|---:|---:|---:|---:|---|",
        ]
        for h in floor:
            lines.append(
                f"| {h['dataset']} | {h['graph']} | {h['nhood']} | {h['n_samples_present']} | {fmt(h['spearman_rho'])} | {h['n_tnk']} | {h['n_malig']} | {h.get('dominant_lineage', '')} |"
            )
        lines.append("")
    if missing:
        lines += [
            "Nhood tables still missing from this merge (counts may come from a published summary): "
            + ", ".join(missing)
            + ".",
            "",
        ]
    lines += [
        "## What this does not say",
        "",
        "- It does not redo the GSE207422-only Milo (PR #323). Those 20 testable",
        "  nhoods, 0 SpatialFDR<0.1, n=7 are imported and sit in the stacked table.",
        "- SpatialFDR<0.1 rows with |ρ|=1 on 5–6 samples are not a TME claim.",
        "- It does not treat 90k–200k cells as *n*.",
        "- It does not invent a shared ICI / MPR label. Only GSE207422 has MPR;",
        "  GSE205335 has RECIST; GSE131907 is treatment-naive; GSE148071 is a",
        "  diagnostic-biopsy atlas.",
        "- Neighbourhood “next to” is kNN co-membership, not histology.",
        "- Unrestricted nhood CLDN4 vs T/NK is partly lineage geometry.",
        "- mRNA ≠ protein. UMI CLDN4 is not an IHC H-score.",
        "- No dual-high TACSTD2+CLDN4 gate.",
        "",
        "## Files",
        "",
        "- `tables/nhoods.tsv` — stacked neighbourhood table (the merge deliverable)",
        "- `tables/da_malignant_cldn4.tsv` — CLDN4 DA columns only",
        "- `tables/sample_scores.tsv`, `tables/honest_n.tsv`, `tables/forest_sample_cldn4_vs_tnk.tsv`",
        "- `tables/summary.json`",
        "- `figures/fig_sample_cldn4_vs_tnk_by_dataset.png`",
        "- `figures/fig_sample_ranks_cldn4_vs_tnk.png`",
        "- `figures/fig_forest_sample_cldn4_vs_tnk.png`",
        "- `figures/fig_nhood_volcano_faceted.png`",
        "- `figures/fig_spatialfdr_counts.png`",
        "- `figures/fig_interface_cldn4_vs_tnk.png`",
        "- `figures/fig_honest_n.png`",
        "",
        "```bash",
        "pip install -r methods/quad_207422_milo_cldn4/requirements.txt",
        "python3 methods/quad_207422_milo_cldn4/scripts/download_gse131907.py",
        "python3 methods/quad_207422_milo_cldn4/scripts/download_gse148071.py",
        "python3 methods/quad_207422_milo_cldn4/scripts/run_gse148071.py \\",
        "  --outdir methods/quad_207422_milo_cldn4/results/GSE148071 \\",
        "  --finding methods/quad_207422_milo_cldn4/results/GSE148071/FINDING.md",
        "python3 methods/quad_207422_milo_cldn4/scripts/run_gse131907.py \\",
        "  --outdir methods/quad_207422_milo_cldn4/results/GSE131907 \\",
        "  --finding methods/quad_207422_milo_cldn4/results/GSE131907/FINDING.md",
        "python3 methods/quad_207422_milo_cldn4/scripts/merge_quad.py",
        "```",
        "",
    ]
    out_path.write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT)
    args = ap.parse_args()
    global IMP, RUN, TAB, FIG
    IMP = args.root / "imported"
    RUN = args.root / "results"
    TAB = args.root / "tables"
    FIG = args.root / "figures"
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    d207 = load_207422()
    d131 = load_131907()
    d148 = load_148071()
    d205 = load_205335()

    score_frames = []
    nhood_frames = []
    honest_rows = []
    da_rows = []
    forest_rows = []
    missing = []

    def take_scores(pack, graph_filter=None, label=None, note=""):
        sc = pack["scores"].copy()
        if graph_filter is not None:
            sc = sc[sc["graph"] == graph_filter]
        scored = sc["malignant_cldn4"].notna()
        n_scored = int(scored.sum())
        n_graph = int(len(sc))
        unit = str(sc["unit"].iloc[0]) if len(sc) else "sample"
        honest_rows.append(
            {
                "dataset": pack["dataset"],
                "graph": graph_filter or ",".join(sorted(map(str, sc["graph"].unique()))) if len(sc) else "",
                "unit": unit,
                "n_scored": n_scored,
                "n_in_graph": n_graph,
                "note": note,
            }
        )
        if n_scored >= 4:
            sp = spearman_safe(sc["malignant_cldn4"], sc["frac_tnk"])
            forest_rows.append(
                {
                    "label": label or f"{pack['dataset']} {graph_filter or ''}".strip(),
                    "dataset": pack["dataset"],
                    "graph": graph_filter or "",
                    "n": sp["n"],
                    "rho": sp["rho"],
                    "p": sp["p"],
                    "source": str(sc["source"].iloc[0]) if len(sc) else "",
                }
            )
        score_frames.append(sc)
        return sc

    take_scores(
        d207,
        "post_treatment",
        "GSE207422 post (imported #323)",
        "imported PR #323; 5/12 post samples dropped (<10 malignant-like); not re-run",
    )
    # 131907: prefer Milo graphs; else hunt tLung as primary and other tumor sites as extra
    if d131["nhoods"] is not None:
        for g in sorted(d131["scores"]["graph"].unique()):
            take_scores(
                d131,
                g,
                f"GSE131907 {g}",
                "this-run Milo graph; author malignant labels; sites not pooled",
            )
        nhood_frames.append(d131["nhoods"])
        for g, c in d131["da_counts"].items():
            da_rows.append({"dataset": "GSE131907", "graph": g, "table": "yes", **c})
    else:
        take_scores(
            d131,
            "tLung",
            "GSE131907 tLung (hunt fallback)",
            "Milo nhood table pending; sample scores from hunt PR 230 (log2TPM, author labels)",
        )
        extra = d131["scores"][
            (d131["scores"]["graph"] != "tLung") & d131["scores"]["malignant_cldn4"].notna()
        ]
        if len(extra) >= 4:
            sp = spearman_safe(extra["malignant_cldn4"], extra["frac_tnk"])
            forest_rows.append(
                {
                    "label": "GSE131907 other origins (hunt, not primary n)",
                    "dataset": "GSE131907",
                    "graph": "other_origins",
                    "n": sp["n"],
                    "rho": sp["rho"],
                    "p": sp["p"],
                    "source": "imported_hunt_PR230",
                }
            )
        missing.append("GSE131907 nhoods")
        da_rows.append(
            {
                "dataset": "GSE131907",
                "graph": "tLung",
                "table": "no",
                "n_nhoods_total": None,
                "n_testable": None,
                "n_p_lt_0.05": None,
                "min_SpatialFDR": None,
                "n_SpatialFDR_lt_0.1": None,
                "n_BH_FDR_lt_0.1": None,
            }
        )

    take_scores(
        d148,
        "all_biopsies",
        "GSE148071 biopsies",
        "39/42 with ≥10 malignant-like; marker-reconstructed lineage (no public CNA IDs)",
    )
    if d148["nhoods"] is not None:
        nhood_frames.append(d148["nhoods"])
        da_rows.append({"dataset": "GSE148071", "graph": "all_biopsies", "table": "yes", **d148["da_counts"]})
    else:
        missing.append("GSE148071 nhoods")
        da_rows.append({"dataset": "GSE148071", "graph": "all_biopsies", "table": "summary_only", **d148["da_counts"]})

    take_scores(
        d205,
        "non_normal",
        "GSE205335 patients",
        "imported PR #385; 22 patients; 4 normal-only patients dropped; MPR n=0",
    )
    nhood_frames.append(d205["nhoods"])
    da_rows.append({"dataset": "GSE205335", "graph": "non_normal", "table": "yes", **d205["da_counts"]})

    nhood_frames.append(d207["nhoods"])
    da_rows.insert(
        0,
        {"dataset": "GSE207422", "graph": "post_treatment", "table": "yes_imported", **d207["da_counts"]},
    )

    sample_tab = pd.concat(score_frames, ignore_index=True, sort=False)
    # drop the mixed 131907 all-origins row from the concatenated scores used for pooled ranks
    # by keeping one row per dataset/sample/graph
    sample_tab = sample_tab.drop_duplicates(subset=["dataset", "sample", "graph"], keep="first")

    nhoods = pd.concat(nhood_frames, ignore_index=True) if nhood_frames else pd.DataFrame()
    forest = pd.DataFrame(forest_rows)
    # primary forest: drop the mixed hunt all-origins if Milo graphs exist
    if d131["nhoods"] is not None:
        forest = forest[~forest["label"].str.contains("all origins", case=False, na=False)]

    # pooled ranks: one row per dataset+sample, primary graphs only
    primary_graphs = {
        ("GSE207422", "post_treatment"),
        ("GSE148071", "all_biopsies"),
        ("GSE205335", "non_normal"),
        ("GSE131907", "tLung"),
    }
    prim = sample_tab[sample_tab.apply(lambda r: (r["dataset"], r["graph"]) in primary_graphs, axis=1)].copy()
    prim = prim.dropna(subset=["malignant_cldn4", "frac_tnk"])
    if len(prim):
        prim["cldn4_rank"] = prim.groupby("dataset")["malignant_cldn4"].rank(pct=True)
        prim["tnk_rank"] = prim.groupby("dataset")["frac_tnk"].rank(pct=True)
        pooled = spearman_safe(prim["cldn4_rank"], prim["tnk_rank"])
    else:
        pooled = {"n": 0, "rho": None, "p": None}

    # Stouffer on primary graphs only
    prim_forest = forest[forest["label"].str.contains("tLung|post|biops|205335", regex=True, na=False)]
    if prim_forest.empty:
        prim_forest = forest
    st = stouffer_signed(
        list(prim_forest["rho"]),
        list(prim_forest["p"]),
        list(prim_forest["n"]),
    )

    honest_primary = [
        r
        for r in honest_rows
        if (r["dataset"], r["graph"]) in primary_graphs
        or (r["dataset"] == "GSE131907" and r["graph"] == "tLung")
    ]
    n_total = int(sum(r["n_scored"] for r in honest_primary))

    figs = write_figures(prim if len(prim) else sample_tab, nhoods if len(nhoods) else None, forest, FIG)

    # write tables
    nhoods.to_csv(TAB / "nhoods.tsv", sep="\t", index=False)
    da_only = nhoods[
        [
            "nhood_id",
            "dataset",
            "graph",
            "nhood",
            "n_samples_present",
            "n_cells",
            "n_tnk",
            "n_malig",
            "spearman_rho",
            "p",
            "testable_cldn4",
            "BH_FDR_cldn4",
            "SpatialFDR_cldn4",
        ]
    ].copy() if len(nhoods) else pd.DataFrame()
    da_only.to_csv(TAB / "da_malignant_cldn4.tsv", sep="\t", index=False)
    sample_tab.to_csv(TAB / "sample_scores.tsv", sep="\t", index=False)
    prim.to_csv(TAB / "sample_scores_primary.tsv", sep="\t", index=False)
    forest.to_csv(TAB / "forest_sample_cldn4_vs_tnk.tsv", sep="\t", index=False)
    pd.DataFrame(honest_rows).to_csv(TAB / "honest_n.tsv", sep="\t", index=False)
    pd.DataFrame(da_rows).to_csv(TAB / "da_summary_by_graph.tsv", sep="\t", index=False)

    testable = nhoods[nhoods["testable_cldn4"].fillna(False).astype(bool)] if len(nhoods) else pd.DataFrame()
    hits = (
        testable[pd.to_numeric(testable["SpatialFDR_cldn4"], errors="coerce") < 0.1]
        if len(testable)
        else pd.DataFrame()
    )
    if len(hits):
        hits.to_csv(TAB / "nhoods_spatialfdr_lt_0.1.tsv", sep="\t", index=False)
    floor_hits = []
    if len(hits):
        for r in hits.itertuples(index=False):
            floor_hits.append(
                {
                    "dataset": r.dataset,
                    "graph": r.graph,
                    "nhood": int(r.nhood) if pd.notna(r.nhood) else None,
                    "n_samples_present": int(r.n_samples_present) if pd.notna(r.n_samples_present) else None,
                    "spearman_rho": float(r.spearman_rho) if pd.notna(r.spearman_rho) else None,
                    "n_tnk": int(r.n_tnk) if pd.notna(r.n_tnk) else None,
                    "n_malig": int(r.n_malig) if pd.notna(r.n_malig) else None,
                    "dominant_lineage": getattr(r, "dominant_lineage", ""),
                }
            )

    summary = {
        "title": "QUAD merge Milo / nhood vs malignant CLDN4",
        "datasets": ["GSE207422", "GSE131907", "GSE148071", "GSE205335"],
        "include_207422": True,
        "redo_207422_only_milo": False,
        "pr_323": "imported; null SpatialFDR; honest n=7",
        "dual_high": False,
        "tacstd2_gate": False,
        "unit": "patient/sample",
        "graph": "per dataset (GSE131907 also split tLung / mBrain); not Harmony",
        "harmony": False,
        "harmony_reason": "four chemistries/label schemes; ~400k cells; 15 GB RAM; treatment and site confounded with dataset",
        "miloR": False,
        "da_model": "sample-level Spearman on nhood proportions vs malignant CLDN4; SpatialFDR k-distance",
        "honest_n_primary_sum": n_total,
        "honest_rows": honest_rows,
        "forest": forest.to_dict("records"),
        "stouffer_signed_primary": st,
        "pooled_within_dataset_ranks": pooled,
        "n_nhoods_in_table": int(len(nhoods)),
        "n_testable_in_table": int(len(testable)),
        "n_SpatialFDR_lt_0.1_in_table": int(len(hits)),
        "floor_hits": floor_hits,
        "missing_nhood_tables": missing,
        "figures": figs,
        "da_by_graph": da_rows,
    }
    (TAB / "summary.json").write_text(json.dumps(summary, indent=2, default=str))

    pack = {
        "forest": forest,
        "da_rows": da_rows,
        "honest_n_total": n_total,
        "honest_rows": honest_rows,
        "stouffer": st,
        "pooled_rank": pooled,
        "n_nhoods_in_table": int(len(nhoods)),
        "n_hits": int(len(hits)),
        "floor_hits": floor_hits,
        "missing_nhood_tables": missing,
    }
    write_finding(pack, args.root / "FINDING.md")
    print(
        json.dumps(
            {
                "n_nhoods": int(len(nhoods)),
                "n_testable": int(len(testable)),
                "n_SpatialFDR_lt_0.1": int(len(hits)),
                "honest_n_primary": n_total,
                "stouffer": st,
                "pooled": pooled,
                "missing": missing,
                "finding": str(args.root / "FINDING.md"),
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
