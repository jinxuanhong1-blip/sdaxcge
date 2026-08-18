#!/usr/bin/env python3
"""ADDITIVE mouse, Cldn4-only: GSE179502 Lkb1 restore in KRAS tumors.

Public Cell Ranger aggr mtx + features + barcodes on GEO. Sorted
neoplastic epithelium (FACS Lin−). Author GEO cohort is Restored vs
NonRestored (n=3 mice vs 3 mice). No dual-high. No human. No GSE179501.

Questions
1. When Lkb1 is restored vs not, does Cldn4 go down and/or IFN / AT2
   go up in neoplastic epithelium?
2. In Cldn4-high vs Cldn4-low cells, is IFN/MHC down and TJ/barrier up?

Honest unit for the restore contrast is the mouse (n=3 vs 3). Cell-level
tests are descriptive (barcodes are not independent).
"""
from __future__ import annotations

import gzip
import json
import math
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import io, stats
from scipy.sparse import csc_matrix

from gene_sets import AT2, IFN_CORE, MHC_CORE, TJ_CORE, load_a8_mouse, present_subset

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TABLES = HERE / "tables"
FIGS = HERE / "figures"

GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179502"
FILES = {
    "barcodes.tsv.gz": f"{GEO_BASE}/suppl/GSE179502_XTR_sorted_scRNAseq_barcodes.tsv.gz",
    "features.tsv.gz": f"{GEO_BASE}/suppl/GSE179502_XTR_sorted_scRNAseq_features.tsv.gz",
    "matrix.mtx.gz": f"{GEO_BASE}/suppl/GSE179502_XTR_sorted_scRNAseq_matrix.mtx.gz",
    "soft.gz": f"{GEO_BASE}/soft/GSE179502_family.soft.gz",
}

MIN_GENES = 200
MIN_UMI = 500
MAX_MITO = 0.20
MIN_N_Q4 = 8
MIN_N_SPEARMAN = 4


def _fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "—"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    return f"{float(x):.{nd}f}"


def ensure_files() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        dest = DATA / name
        if dest.exists() and dest.stat().st_size > 1000:
            continue
        print("DOWNLOAD", url)
        urllib.request.urlretrieve(url, dest)


def _maybe_text(path: Path) -> str:
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        return gzip.decompress(raw).decode()
    return raw.decode()


def parse_soft() -> pd.DataFrame:
    text = _maybe_text(DATA / "soft.gz")
    rows = []
    cur: dict | None = None
    for line in text.splitlines():
        if line.startswith("^SAMPLE"):
            if cur:
                rows.append(cur)
            cur = {"gsm": line.split("=", 1)[1].strip()}
        elif cur is not None and line.startswith("!Sample_"):
            k, v = line.split("=", 1)
            k = k.replace("!Sample_", "").strip()
            v = v.strip()
            if k == "characteristics_ch1":
                if ":" in v:
                    ck, cv = v.split(":", 1)
                    cur[ck.strip().lower()] = cv.strip()
            elif k == "title":
                cur["mouse"] = v
            elif k == "geo_accession":
                cur["gsm"] = v
    if cur:
        rows.append(cur)
    meta = pd.DataFrame(rows)
    meta["restorable"] = meta["genotype"].str.contains("FLPo", na=False)
    meta["lkb1_restored"] = meta["cohort"].eq("Restored")
    return meta


def load_barcodes() -> pd.DataFrame:
    text = _maybe_text(DATA / "barcodes.tsv.gz")
    bcs = [ln.strip() for ln in text.splitlines() if ln.strip()]
    mice = [b.split("_", 1)[0] for b in bcs]
    return pd.DataFrame({"barcode": bcs, "mouse": mice})


def load_features() -> pd.DataFrame:
    raw = (DATA / "features.tsv.gz").read_bytes()
    if raw[:2] == b"\x1f\x8b":
        text = gzip.decompress(raw).decode()
    else:
        text = raw.decode()
    rows = [ln.split("\t") for ln in text.splitlines() if ln.strip()]
    feat = pd.DataFrame(rows, columns=["ensembl", "symbol", "feature_type"][: len(rows[0])])
    feat["idx"] = np.arange(len(feat), dtype=int)
    return feat


def load_matrix() -> csc_matrix:
    print("READ mtx")
    mat = io.mmread(DATA / "matrix.mtx.gz")
    if not hasattr(mat, "tocsc"):
        raise RuntimeError("expected sparse Matrix Market counts")
    return mat.tocsc()


def qc_mask(mat: csc_matrix, symbols: np.ndarray) -> tuple[np.ndarray, pd.DataFrame]:
    n_umi = np.asarray(mat.sum(axis=0)).ravel()
    n_genes = np.asarray((mat > 0).sum(axis=0)).ravel()
    mito = np.array([str(s).startswith("mt-") for s in symbols])
    mito_umi = np.asarray(mat[mito, :].sum(axis=0)).ravel() if mito.any() else np.zeros_like(n_umi)
    mito_frac = np.divide(mito_umi, n_umi, out=np.zeros_like(n_umi, dtype=float), where=n_umi > 0)
    keep = (n_genes >= MIN_GENES) & (n_umi >= MIN_UMI) & (mito_frac <= MAX_MITO)
    qc = pd.DataFrame(
        {
            "n_umi": n_umi,
            "n_genes": n_genes,
            "mito_frac": mito_frac,
            "pass_qc": keep,
        }
    )
    return keep, qc


def log_cp10k(mat: csc_matrix) -> csc_matrix:
    n_umi = np.asarray(mat.sum(axis=0)).ravel().astype(float)
    n_umi[n_umi == 0] = 1.0
    # scale columns then log1p; keep CSC
    mat = mat.astype(np.float32)
    mat = mat.multiply(10000.0 / n_umi)
    mat.data = np.log1p(mat.data)
    return mat.tocsc()


def gene_vector(mat: csc_matrix, symbols: np.ndarray, gene: str) -> np.ndarray:
    hits = np.where(symbols == gene)[0]
    if hits.size == 0:
        return np.full(mat.shape[1], np.nan)
    return np.asarray(mat[hits[0], :].todense()).ravel()


def set_score(mat: csc_matrix, symbols: np.ndarray, genes: list[str]) -> np.ndarray:
    idx = [int(np.where(symbols == g)[0][0]) for g in genes if g in set(symbols)]
    if not idx:
        return np.full(mat.shape[1], np.nan)
    return np.asarray(mat[idx, :].mean(axis=0)).ravel()


def detect_frac(counts: csc_matrix, symbols: np.ndarray, gene: str, keep: np.ndarray) -> float:
    hits = np.where(symbols == gene)[0]
    if hits.size == 0:
        return float("nan")
    x = np.asarray(counts[hits[0], :].todense()).ravel()
    x = x[keep]
    return float(np.mean(x > 0)) if x.size else float("nan")


def welch_mwu(a, b) -> dict:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    out = {
        "n_a": int(a.size),
        "n_b": int(b.size),
        "mean_a": float(a.mean()) if a.size else float("nan"),
        "mean_b": float(b.mean()) if b.size else float("nan"),
        "median_a": float(np.median(a)) if a.size else float("nan"),
        "median_b": float(np.median(b)) if b.size else float("nan"),
        "delta_mean": float("nan"),
        "delta_median": float("nan"),
        "welch_p": float("nan"),
        "mwu_p": float("nan"),
        "r_rb": float("nan"),
    }
    if a.size == 0 or b.size == 0:
        return out
    out["delta_mean"] = out["mean_b"] - out["mean_a"]
    out["delta_median"] = out["median_b"] - out["median_a"]
    if a.size >= 2 and b.size >= 2 and a.std(ddof=1) + b.std(ddof=1) > 0:
        _, p = stats.ttest_ind(b, a, equal_var=False)
        out["welch_p"] = float(p)
    if a.size >= 1 and b.size >= 1:
        try:
            u, p = stats.mannwhitneyu(b, a, alternative="two-sided")
            out["mwu_p"] = float(p)
            out["r_rb"] = (2.0 * float(u)) / (b.size * a.size) - 1.0
        except ValueError:
            pass
    return out


def q4_vs_q1(rank_on, y, min_n: int = MIN_N_Q4) -> dict:
    """Top vs bottom quartile by rank. Ties (Cldn4 zeros) stay in Q1.

    pd.qcut(..., duplicates='drop') fails when >25% of cells are zero.
    Explicit rank tails always return two groups when n>=min_n.
    """
    c = np.asarray(rank_on, float)
    y = np.asarray(y, float)
    m = np.isfinite(c) & np.isfinite(y)
    c, y = c[m], y[m]
    n = int(c.size)
    out = {
        "n": n,
        "n_q1": np.nan,
        "n_q4": np.nan,
        "median_q1": np.nan,
        "median_q4": np.nan,
        "delta_median": np.nan,
        "r_rb": np.nan,
        "p": np.nan,
        "usable": False,
    }
    n_tail = n // 4
    if n < min_n or n_tail < 2:
        return out
    order = np.argsort(c, kind="mergesort")
    q1 = y[order[:n_tail]]
    q4 = y[order[-n_tail:]]
    n1, n4 = int(q1.size), int(q4.size)
    u, p = stats.mannwhitneyu(q4, q1, alternative="two-sided")
    out.update(
        {
            "n_q1": n1,
            "n_q4": n4,
            "median_q1": float(np.median(q1)),
            "median_q4": float(np.median(q4)),
            "delta_median": float(np.median(q4) - np.median(q1)),
            "r_rb": (2.0 * float(u)) / (n4 * n1) - 1.0,
            "p": float(p),
            "usable": True,
        }
    )
    return out


def pos_vs_neg(counts, y) -> dict:
    c = np.asarray(counts, float)
    y = np.asarray(y, float)
    m = np.isfinite(c) & np.isfinite(y)
    c, y = c[m], y[m]
    neg = y[c <= 0]
    pos = y[c > 0]
    d = welch_mwu(neg, pos)
    return {
        "n": int(c.size),
        "n_q1": d["n_a"],
        "n_q4": d["n_b"],
        "median_q1": d["median_a"],
        "median_q4": d["median_b"],
        "delta_median": d["delta_median"],
        "r_rb": d["r_rb"],
        "p": d["mwu_p"],
        "usable": d["n_a"] >= 2 and d["n_b"] >= 2,
    }


def partial_spearman(x, y, z) -> tuple[float, float, int]:
    """Spearman of residuals after ranking out z (depth)."""
    xa = np.asarray(x, float)
    ya = np.asarray(y, float)
    za = np.asarray(z, float)
    m = np.isfinite(xa) & np.isfinite(ya) & np.isfinite(za)
    xa, ya, za = xa[m], ya[m], za[m]
    n = int(xa.size)
    if n < MIN_N_SPEARMAN:
        return float("nan"), float("nan"), n
    rx, ry, rz = stats.rankdata(xa), stats.rankdata(ya), stats.rankdata(za)
    rx = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    ry = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    rho, p = stats.spearmanr(rx, ry)
    return float(rho), float(p), n


def spearman(x, y) -> tuple[float, float, int]:
    xa = np.asarray(x, float)
    ya = np.asarray(y, float)
    m = np.isfinite(xa) & np.isfinite(ya)
    xa, ya = xa[m], ya[m]
    n = int(xa.size)
    if n < MIN_N_SPEARMAN:
        return float("nan"), float("nan"), n
    rho, p = stats.spearmanr(xa, ya)
    return float(rho), float(p), n


def coverage_row(name: str, wanted: list[str], present: list[str]) -> dict:
    return {
        "set": name,
        "n_wanted": len(wanted),
        "n_present": len(present),
        "genes_present": ",".join(present),
        "genes_missing": ",".join([g for g in wanted if g not in set(present)]),
    }


def draw_restore_figure(cells: pd.DataFrame, mice: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6))
    order = ["NonRestored", "Restored"]
    colors = {"NonRestored": "#b56576", "Restored": "#2a9d8f"}
    ax = axes[0]
    data = [cells.loc[cells["cohort"] == c, "Cldn4"].to_numpy() for c in order]
    parts = ax.violinplot(data, positions=[0, 1], showextrema=False, widths=0.8)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(colors[order[i]])
        pc.set_alpha(0.35)
    for i, c in enumerate(order):
        sub = mice[mice["cohort"] == c]
        ax.scatter(
            np.full(len(sub), i),
            sub["Cldn4"],
            c=colors[c],
            s=46,
            zorder=3,
            edgecolors="k",
            linewidths=0.4,
        )
        for _, r in sub.iterrows():
            ax.annotate(r["mouse"], (i + 0.08, r["Cldn4"]), fontsize=7, va="center")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["NonRestored\n(n=3 mice)", "Restored\n(n=3 mice)"])
    ax.set_ylabel("Cldn4 log1p(CP10k)")
    ax.set_title("Cldn4 by Lkb1-restore cohort")
    ax = axes[1]
    for lab, col in [("IFN_A8", "C0"), ("AT2", "C2"), ("MHC_A8", "C1")]:
        ax.plot(
            [0, 1],
            [
                mice.loc[mice["cohort"] == "NonRestored", lab].mean(),
                mice.loc[mice["cohort"] == "Restored", lab].mean(),
            ],
            color=col,
            lw=1.2,
            alpha=0.5,
        )
        for i, c in enumerate(order):
            sub = mice[mice["cohort"] == c]
            ax.scatter(
                np.full(len(sub), i) + (0.04 if lab == "AT2" else (-0.04 if lab == "IFN_A8" else 0)),
                sub[lab],
                label=lab if i == 0 else None,
                s=40,
                color=col,
                edgecolors="k",
                linewidths=0.4,
            )
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["NonRestored", "Restored"])
    ax.set_ylabel("mouse-mean set score")
    ax.set_title("IFN / AT2 / MHC-I by restore")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def draw_highlow_figure(rows: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    fams = ["IFN_A8", "MHC_A8", "TJ_A8", "AT2"]
    sub = rows[
        rows["stratum"].eq("all_cells")
        & rows["gate"].eq("Q4_vs_Q1_rank_tails")
        & rows["family"].isin(fams)
    ].copy()
    sub = sub.set_index("family").loc[fams]
    x = np.arange(len(fams))
    colors = ["#2a9d8f" if d < 0 else "#b56576" for d in sub["delta_median"]]
    ax.bar(x, sub["delta_median"], color=colors, edgecolor="k", linewidth=0.4)
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(["IFN", "MHC-I/APM", "TJ (Cldn4 out)", "AT2"])
    ax.set_ylabel("Δ median (Cldn4 Q4 − Q1)")
    ax.set_title("Cldn4-high vs low neoplastic cells")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    ensure_files()

    meta = parse_soft()
    bc = load_barcodes()
    feat = load_features()
    symbols = feat["symbol"].to_numpy()
    present = set(symbols.tolist())
    a8 = load_a8_mouse(present)
    sets = {
        "AT2": present_subset(AT2, present),
        "IFN_CORE": present_subset(IFN_CORE, present),
        "MHC_CORE": present_subset(MHC_CORE, present),
        "TJ_CORE": present_subset(TJ_CORE, present),
        **a8,
    }
    cov = pd.DataFrame(
        [
            coverage_row("AT2", AT2, sets["AT2"]),
            coverage_row("IFN_CORE", IFN_CORE, sets["IFN_CORE"]),
            coverage_row("MHC_CORE", MHC_CORE, sets["MHC_CORE"]),
            coverage_row("TJ_CORE", TJ_CORE, sets["TJ_CORE"]),
            coverage_row("IFN_A8", [], sets["IFN_A8"]),
            coverage_row("MHC_A8", [], sets["MHC_A8"]),
            coverage_row("TJ_A8", [], sets["TJ_A8"]),
        ]
    )
    cov.to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)

    counts = load_matrix()
    if counts.shape != (len(feat), len(bc)):
        raise RuntimeError(f"shape mismatch {counts.shape} vs {len(feat)} x {len(bc)}")

    keep, qc = qc_mask(counts, symbols)
    print(f"QC keep {int(keep.sum())} / {keep.size}")
    mat = log_cp10k(counts[:, keep])
    cells = bc.loc[keep].reset_index(drop=True)
    cells = cells.merge(
        meta[
            ["mouse", "gsm", "cohort", "treatment", "sex", "genotype", "restorable", "lkb1_restored"]
        ],
        on="mouse",
        how="left",
    )
    cells = pd.concat([cells, qc.loc[keep].reset_index(drop=True)], axis=1)
    raw_keep = counts[:, keep]

    cells["Cldn4"] = gene_vector(mat, symbols, "Cldn4")
    cells["Stk11"] = gene_vector(mat, symbols, "Stk11")
    cells["Tacstd2"] = gene_vector(mat, symbols, "Tacstd2")
    cells["Epcam"] = gene_vector(mat, symbols, "Epcam")
    cells["Cldn4_count"] = np.asarray(raw_keep[int(np.where(symbols == "Cldn4")[0][0]), :].todense()).ravel()
    for name, genes in sets.items():
        cells[name] = set_score(mat, symbols, genes)

    # mouse-level means (honest unit)
    mouse_rows = []
    for mouse, sub in cells.groupby("mouse", sort=False):
        row = {
            "mouse": mouse,
            "gsm": sub["gsm"].iloc[0],
            "cohort": sub["cohort"].iloc[0],
            "treatment": sub["treatment"].iloc[0],
            "sex": sub["sex"].iloc[0],
            "restorable": bool(sub["restorable"].iloc[0]),
            "lkb1_restored": bool(sub["lkb1_restored"].iloc[0]),
            "n_cells_qc": int(len(sub)),
            "n_cells_raw": int((bc["mouse"] == mouse).sum()),
            "frac_Cldn4_pos": float((sub["Cldn4_count"] > 0).mean()),
            "median_n_umi": float(sub["n_umi"].median()),
        }
        for col in ["Cldn4", "Stk11", "Tacstd2", "Epcam"] + list(sets):
            row[col] = float(sub[col].mean())
        mouse_rows.append(row)
    mice = pd.DataFrame(mouse_rows)
    mice.to_csv(TABLES / "mouse_scores.tsv", sep="\t", index=False)

    samples = meta.merge(
        mice[["mouse", "n_cells_qc", "n_cells_raw", "frac_Cldn4_pos"]],
        on="mouse",
        how="left",
    )
    samples.to_csv(TABLES / "samples.tsv", sep="\t", index=False)

    rest = mice[mice["lkb1_restored"]]
    nonr = mice[~mice["lkb1_restored"]]
    contrast_rows = []
    for feature in ["Cldn4", "Stk11", "Tacstd2", "Epcam", "frac_Cldn4_pos"] + list(sets):
        d = welch_mwu(nonr[feature], rest[feature])
        contrast_rows.append(
            {
                "contrast": "Restored_vs_NonRestored",
                "unit": "mouse",
                "feature": feature,
                **d,
            }
        )
        # cell-level descriptive
        dcell = welch_mwu(
            cells.loc[~cells["lkb1_restored"], feature if feature != "frac_Cldn4_pos" else "Cldn4"],
            cells.loc[cells["lkb1_restored"], feature if feature != "frac_Cldn4_pos" else "Cldn4"],
        )
        if feature == "frac_Cldn4_pos":
            continue
        contrast_rows.append(
            {
                "contrast": "Restored_vs_NonRestored",
                "unit": "cell_descriptive",
                "feature": feature,
                **dcell,
            }
        )
    contrasts = pd.DataFrame(contrast_rows)
    contrasts.to_csv(TABLES / "restore_contrasts.tsv", sep="\t", index=False)

    # Cldn4-high vs low
    hl_rows = []
    families = ["IFN_A8", "MHC_A8", "TJ_A8", "AT2", "IFN_CORE", "MHC_CORE", "TJ_CORE", "Stk11"]
    strata = {
        "all_cells": cells,
        "NonRestored": cells[~cells["lkb1_restored"]],
        "Restored": cells[cells["lkb1_restored"]],
    }
    for stratum, sub in strata.items():
        for fam in families:
            q = q4_vs_q1(sub["Cldn4"], sub[fam])
            det = pos_vs_neg(sub["Cldn4_count"], sub[fam])
            rho, p, n = spearman(sub["Cldn4"], sub[fam])
            prho, pp, pn = partial_spearman(sub["Cldn4"], sub[fam], np.log1p(sub["n_umi"]))
            hl_rows.append(
                {
                    "stratum": stratum,
                    "gate": "Q4_vs_Q1_rank_tails",
                    "unit": "cell_descriptive",
                    "family": fam,
                    **q,
                    "spearman_rho": rho,
                    "spearman_p": p,
                    "spearman_n": n,
                    "partial_rho_umi": prho,
                    "partial_p_umi": pp,
                }
            )
            hl_rows.append(
                {
                    "stratum": stratum,
                    "gate": "Cldn4pos_vs_neg",
                    "unit": "cell_descriptive",
                    "family": fam,
                    **det,
                    "spearman_rho": rho,
                    "spearman_p": p,
                    "spearman_n": n,
                    "partial_rho_umi": prho,
                    "partial_p_umi": pp,
                }
            )
    # within-mouse Q4 vs Q1, then mouse as unit
    per_mouse_hl = []
    for mouse, sub in cells.groupby("mouse", sort=False):
        for fam in families:
            q = q4_vs_q1(sub["Cldn4"], sub[fam], min_n=MIN_N_Q4)
            det = pos_vs_neg(sub["Cldn4_count"], sub[fam])
            per_mouse_hl.append(
                {
                    "mouse": mouse,
                    "cohort": sub["cohort"].iloc[0],
                    "lkb1_restored": bool(sub["lkb1_restored"].iloc[0]),
                    "gate": "Q4_vs_Q1_rank_tails",
                    "family": fam,
                    **q,
                }
            )
            per_mouse_hl.append(
                {
                    "mouse": mouse,
                    "cohort": sub["cohort"].iloc[0],
                    "lkb1_restored": bool(sub["lkb1_restored"].iloc[0]),
                    "gate": "Cldn4pos_vs_neg",
                    "family": fam,
                    **det,
                }
            )
    per_mouse_hl_df = pd.DataFrame(per_mouse_hl)
    per_mouse_hl_df.to_csv(TABLES / "cldn4_highlow_by_mouse.tsv", sep="\t", index=False)
    for gate in ("Q4_vs_Q1_rank_tails", "Cldn4pos_vs_neg"):
        for fam in families:
            dlt = per_mouse_hl_df.loc[
                per_mouse_hl_df["family"].eq(fam)
                & per_mouse_hl_df["gate"].eq(gate)
                & per_mouse_hl_df["usable"],
                "delta_median",
            ]
            if len(dlt) >= 4:
                w, p = stats.wilcoxon(dlt.to_numpy(), alternative="two-sided")
            else:
                w, p = float("nan"), float("nan")
            hl_rows.append(
                {
                    "stratum": "within_mouse_then_wilcoxon",
                    "gate": gate,
                    "unit": "mouse",
                    "family": fam,
                    "n": int(dlt.size),
                    "n_q1": np.nan,
                    "n_q4": np.nan,
                    "median_q1": np.nan,
                    "median_q4": np.nan,
                    "delta_median": float(dlt.median()) if len(dlt) else float("nan"),
                    "r_rb": float((dlt > 0).mean() - (dlt < 0).mean()) if len(dlt) else float("nan"),
                    "p": float(p) if math.isfinite(p) else float("nan"),
                    "usable": len(dlt) >= 4,
                    "spearman_rho": float("nan"),
                    "spearman_p": float("nan"),
                    "spearman_n": int(dlt.size),
                    "wilcoxon_W": float(w) if math.isfinite(w) else float("nan"),
                    "n_mice_delta_neg": int((dlt < 0).sum()) if len(dlt) else 0,
                    "n_mice_delta_pos": int((dlt > 0).sum()) if len(dlt) else 0,
                }
            )
    # mouse-mean Spearman (n=6)
    for fam in families:
        rho, p, n = spearman(mice["Cldn4"], mice[fam] if fam in mice.columns else mice["Stk11"])
        hl_rows.append(
            {
                "stratum": "mouse_mean",
                "gate": "mouse_mean_spearman",
                "unit": "mouse",
                "family": fam,
                "n": n,
                "n_q1": np.nan,
                "n_q4": np.nan,
                "median_q1": np.nan,
                "median_q4": np.nan,
                "delta_median": float("nan"),
                "r_rb": float("nan"),
                "p": p,
                "usable": n >= MIN_N_SPEARMAN,
                "spearman_rho": rho,
                "spearman_p": p,
                "spearman_n": n,
            }
        )
    hl = pd.DataFrame(hl_rows)
    hl.to_csv(TABLES / "cldn4_highlow.tsv", sep="\t", index=False)
    ifn_tab = hl[
        hl["family"].isin(["IFN_A8", "IFN_CORE", "MHC_A8", "MHC_CORE", "TJ_A8", "TJ_CORE", "AT2"])
        & hl["gate"].isin(["Q4_vs_Q1_rank_tails", "Cldn4pos_vs_neg", "mouse_mean_spearman"])
    ].copy()
    ifn_tab.to_csv(TABLES / "cldn4_highlow_ifn.tsv", sep="\t", index=False)

    # one-row restore table
    def grab(feature, unit="mouse"):
        r = contrasts[
            contrasts["feature"].eq(feature) & contrasts["unit"].eq(unit)
        ].iloc[0]
        return r

    one = pd.DataFrame(
        [
            {
                "dataset": "GSE179502",
                "model": "KT;Lkb1XTR KRAS lung, FACS neoplastic scRNA",
                "contrast": "Restored vs NonRestored (author GEO cohort)",
                "n_mice_nonrestored": int((~mice["lkb1_restored"]).sum()),
                "n_mice_restored": int(mice["lkb1_restored"].sum()),
                "n_cells_qc": int(len(cells)),
                "Cldn4_delta_mouse": grab("Cldn4")["delta_mean"],
                "Cldn4_welch_p_mouse": grab("Cldn4")["welch_p"],
                "Cldn4_mwu_p_mouse": grab("Cldn4")["mwu_p"],
                "AT2_delta_mouse": grab("AT2")["delta_mean"],
                "AT2_welch_p_mouse": grab("AT2")["welch_p"],
                "IFN_A8_delta_mouse": grab("IFN_A8")["delta_mean"],
                "IFN_A8_welch_p_mouse": grab("IFN_A8")["welch_p"],
                "MHC_A8_delta_mouse": grab("MHC_A8")["delta_mean"],
                "MHC_A8_welch_p_mouse": grab("MHC_A8")["welch_p"],
                "Stk11_delta_mouse": grab("Stk11")["delta_mean"],
                "Stk11_welch_p_mouse": grab("Stk11")["welch_p"],
            }
        ]
    )
    one.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    draw_restore_figure(cells, mice, FIGS / "cldn4_by_restore")
    draw_highlow_figure(hl, FIGS / "cldn4_highlow_ifn_tj")

    summary = {
        "accession": "GSE179502",
        "pmid": 35228570,
        "species": "mouse",
        "compartment": "FACS-sorted neoplastic epithelium",
        "n_mice": int(mice["mouse"].nunique()),
        "n_restored": int(mice["lkb1_restored"].sum()),
        "n_nonrestored": int((~mice["lkb1_restored"]).sum()),
        "n_cells_raw": int(len(bc)),
        "n_cells_qc": int(len(cells)),
        "qc": {"min_genes": MIN_GENES, "min_umi": MIN_UMI, "max_mito": MAX_MITO},
        "mice": mice.to_dict(orient="records"),
        "restore_mouse": {
            k: contrasts[
                contrasts["feature"].eq(k) & contrasts["unit"].eq("mouse")
            ]
            .iloc[0]
            .to_dict()
            for k in ["Cldn4", "Stk11", "AT2", "IFN_A8", "MHC_A8", "TJ_A8", "Tacstd2"]
        },
        "highlow_all_cells": hl[
            hl["stratum"].eq("all_cells") & hl["gate"].eq("Q4_vs_Q1_rank_tails")
        ][
            [
                "family",
                "n",
                "n_q1",
                "n_q4",
                "delta_median",
                "r_rb",
                "p",
                "spearman_rho",
                "spearman_p",
                "partial_rho_umi",
                "partial_p_umi",
            ]
        ].to_dict(orient="records"),
        "highlow_within_mouse": hl[
            hl["stratum"].eq("within_mouse_then_wilcoxon")
            & hl["gate"].eq("Q4_vs_Q1_rank_tails")
        ][["family", "n", "delta_median", "p", "n_mice_delta_neg", "n_mice_delta_pos"]].to_dict(
            orient="records"
        ),
        "gene_coverage": cov[["set", "n_present"]].to_dict(orient="records"),
        "notes": [
            "Cldn4 only. Tacstd2 audit. No dual-high. No human. No GSE179501.",
            "Author GEO cohort Restored = restorable + tamoxifen (3 mice).",
            "NonRestored = non-restorable ± tamoxifen/vehicle or restorable + vehicle (3 mice).",
            "MWU at 3 vs 3 cannot go below p=0.1.",
        ],
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps({k: summary[k] for k in ["n_mice", "n_cells_qc", "restore_mouse"]}, indent=2, default=str))
    print("WROTE", TABLES, FIGS)


if __name__ == "__main__":
    main()
