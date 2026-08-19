#!/usr/bin/env python3
"""ADDITIVE public Visium, CLDN4-only.

Sources (both run if downloadable):
  A) GSE277206 — 10x Visium CytAssist FFPE, never-smoker LUAD MIA (n=2 sections)
  B) Zenodo 13337961 — Frontiers 2024 (10.3389/fimmu.2024.1430163) lepidic vs solid LUAD

No private 8-KL. No fabricated stats. Section is the unit.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io, stats
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "visium_public_luad"
COORD_DIR = Path(__file__).resolve().parent / "coords"
OUT = ROOT / "results" / "visium_public_luad"
TABLES = OUT / "tables"
MAPS = OUT / "maps"
TABLES.mkdir(parents=True, exist_ok=True)
MAPS.mkdir(parents=True, exist_ok=True)

MIN_GENES = 200
MIN_NEI = 3
MIN_N_SPEARMAN = 8
MIN_EPI = 20
PITCH_UM = 100.0  # standard Visium (not HD) center-to-center

EPI_GENES = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]


def dump_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=_json_default) + "\n")


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def fmt_p(p):
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "NA"
    if p == 0 or p < 1e-300:
        return "<1e-300"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4g}"


def fmt_rho(r):
    if r is None or (isinstance(r, float) and not np.isfinite(r)):
        return "NA"
    return f"{r:.3f}"


def spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < MIN_N_SPEARMAN:
        return {"n": n, "rho": None, "p": None}
    rho, p = stats.spearmanr(x[m], y[m])
    return {"n": n, "rho": float(rho), "p": float(p)}


def residualize(y, z):
    """Linear residual of y on z (both finite). Returns y-sized array with NaNs preserved."""
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    out = np.full(y.shape, np.nan, dtype=float)
    m = np.isfinite(y) & np.isfinite(z)
    if int(m.sum()) < 8:
        return out
    zz = z[m]
    yy = y[m]
    if np.nanstd(zz) < 1e-12:
        out[m] = yy - yy.mean()
        return out
    b = np.polyfit(zz, yy, 1)
    out[m] = yy - (b[0] * zz + b[1])
    return out


def partial_spearman(x, y, z):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(m.sum())
    if n < 12:
        return {"n": n, "rho": None, "p": None}
    rx = stats.rankdata(x[m])
    ry = stats.rankdata(y[m])
    rz = stats.rankdata(z[m])
    rz = (rz - rz.mean()) / (rz.std() + 1e-12)
    bx = np.polyfit(rz, rx, 1)
    by = np.polyfit(rz, ry, 1)
    ex = rx - (bx[0] * rz + bx[1])
    ey = ry - (by[0] * rz + by[1])
    rho, p = stats.spearmanr(ex, ey)
    return {"n": n, "rho": float(rho), "p": float(p)}


def mw_q4_q1(q4, q1):
    a = np.asarray(q4, dtype=float)
    b = np.asarray(q1, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    out = {
        "n_q4": int(a.size),
        "n_q1": int(b.size),
        "median_q4": float(np.median(a)) if a.size else None,
        "median_q1": float(np.median(b)) if b.size else None,
        "delta_q4_minus_q1": None,
        "p": None,
    }
    if a.size and b.size:
        out["delta_q4_minus_q1"] = float(np.median(a) - np.median(b))
    if a.size >= 5 and b.size >= 5:
        try:
            out["p"] = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
        except ValueError:
            out["p"] = None
    return out


def hex_neighbors(row, col):
    return [
        (row, col - 2),
        (row, col + 2),
        (row - 1, col - 1),
        (row - 1, col + 1),
        (row + 1, col - 1),
        (row + 1, col + 1),
    ]


def hex_xy_um(row, col, pitch=PITCH_UM):
    x = (col / 2.0) * pitch
    y = row * (np.sqrt(3.0) / 2.0) * pitch
    return x, y


def log_cp10k(mtx):
    lib = np.asarray(mtx.sum(axis=1)).ravel().astype(float)
    lib[lib == 0] = 1.0
    return np.log1p(mtx.multiply(1e4 / lib[:, None]).toarray())


def gene_index(genes):
    idx = {}
    for i, g in enumerate(genes):
        idx.setdefault(str(g).upper(), i)
    return idx


def score_from_log(logx, gidx, names):
    ii = [gidx[g] for g in names if g in gidx]
    if not ii:
        return np.full(logx.shape[0], np.nan)
    return logx[:, ii].mean(axis=1)


def load_coord_map(kind: str) -> pd.DataFrame:
    if kind == "v5_11mm":
        path = COORD_DIR / "visium_v5_cytassist_11mm.csv"
    elif kind == "v4_6p5mm":
        path = COORD_DIR / "visium_v4_cytassist_6p5mm.csv"
    else:
        raise ValueError(kind)
    df = pd.read_csv(path)
    df["barcode"] = df["barcode"].astype(str)
    return df.set_index("barcode")[["array_row", "array_col"]]


def load_10x_h5(path):
    import h5py

    with h5py.File(path, "r") as f:
        mat = f["matrix"]
        data = mat["data"][:]
        indices = mat["indices"][:]
        indptr = mat["indptr"][:]
        shape = tuple(int(x) for x in mat["shape"][:])
        barcodes = [x.decode() if isinstance(x, bytes) else str(x) for x in mat["barcodes"][:]]
        genes = [x.decode() if isinstance(x, bytes) else str(x) for x in mat["features"]["name"][:]]
        mtx = csr_matrix((data, indices, indptr), shape=(shape[1], shape[0]))
    return pd.Series(barcodes, dtype=str), genes, mtx


def load_mtx_dir(d):
    d = Path(d)
    barcodes = pd.read_csv(d / "barcodes.tsv.gz", header=None)[0].astype(str)
    feat = pd.read_csv(d / "features.tsv.gz", sep="\t", header=None)
    genes = feat[1].astype(str).tolist() if feat.shape[1] > 1 else feat[0].astype(str).tolist()
    mtx = io.mmread(d / "matrix.mtx.gz").tocsr()
    if mtx.shape[0] == len(barcodes) and mtx.shape[1] == len(genes):
        pass
    elif mtx.shape[1] == len(barcodes) and mtx.shape[0] == len(genes):
        mtx = mtx.T.tocsr()
    else:
        raise ValueError(f"shape mismatch mtx={mtx.shape} bc={len(barcodes)} genes={len(genes)}")
    return barcodes, genes, mtx


def quartiles(x):
    x = np.asarray(x, dtype=float)
    q1, q3 = np.nanquantile(x, [0.25, 0.75])
    return q1, q3


def analyze_section(name, source, histology, barcodes, genes, mtx, pos):
    gidx = gene_index(genes)
    rec = {
        "section": name,
        "source": source,
        "histology": histology,
        "n_genes_matrix": len(gidx),
        "CLDN4": "CLDN4" in gidx,
        "CD8A": "CD8A" in gidx,
        "KRT8": "KRT8" in gidx,
        "epi_genes_present": [g for g in EPI_GENES if g in gidx],
    }
    if "CLDN4" not in gidx or "CD8A" not in gidx:
        rec["skip"] = "missing_CLDN4_or_CD8A"
        return rec

    n_genes = np.asarray((mtx > 0).sum(axis=1)).ravel()
    pos_all = pos.reindex(barcodes)
    have_xy = pos_all["array_row"].notna() & pos_all["array_col"].notna()
    keep = (n_genes >= MIN_GENES) & have_xy.to_numpy()
    rec["n_raw"] = int(len(barcodes))
    rec["n_with_coords"] = int(have_xy.sum())
    rec["n_qc"] = int(keep.sum())
    rec["coord_match_frac"] = float(have_xy.mean()) if len(barcodes) else 0.0
    if int(keep.sum()) < 50:
        rec["skip"] = "too_few_qc_spots"
        return rec

    mtx_k = mtx[keep]
    bc_k = np.asarray(barcodes)[keep]
    pos_k = pos.reindex(bc_k)
    rows = pos_k["array_row"].to_numpy(dtype=float)
    cols = pos_k["array_col"].to_numpy(dtype=float)
    logx = log_cp10k(mtx_k)
    cl4 = score_from_log(logx, gidx, ["CLDN4"])
    cd8 = score_from_log(logx, gidx, ["CD8A"])
    krt8 = score_from_log(logx, gidx, ["KRT8"])
    epi = score_from_log(logx, gidx, rec["epi_genes_present"])
    cl4_krt8_resid = residualize(cl4, krt8)

    rec["same_CLDN4_CD8A"] = spearman(cl4, cd8)
    rec["same_KRT8_CD8A"] = spearman(krt8, cd8)
    rec["same_CLDN4_KRT8"] = spearman(cl4, krt8)
    rec["partial_CLDN4_CD8A_ctrl_KRT8"] = partial_spearman(cl4, cd8, krt8)
    rec["same_CLDN4residKRT8_CD8A"] = spearman(cl4_krt8_resid, cd8)

    # epithelial-like = section Q3+ of epi score
    epi_cut = float(np.nanquantile(epi, 0.75))
    epi_like = np.isfinite(epi) & (epi >= epi_cut)
    rec["n_epithelial_like"] = int(epi_like.sum())
    rec["epi_q3_cut"] = epi_cut
    if int(epi_like.sum()) < MIN_EPI:
        rec["skip"] = "too_few_epithelial_like"
        rec["note"] = f"epithelial-like n={int(epi_like.sum())} < {MIN_EPI}"
        return rec

    cl4_e = cl4[epi_like]
    q1c, q3c = quartiles(cl4_e)
    q1 = epi_like & (cl4 <= q1c)
    q4 = epi_like & (cl4 >= q3c)
    rec["n_cldn4_q1"] = int(q1.sum())
    rec["n_cldn4_q4"] = int(q4.sum())
    rec["cldn4_q1_cut"] = float(q1c)
    rec["cldn4_q4_cut"] = float(q3c)

    r1c, r3c = quartiles(cl4_krt8_resid[epi_like])
    rq1 = epi_like & (cl4_krt8_resid <= r1c)
    rq4 = epi_like & (cl4_krt8_resid >= r3c)
    rec["n_resid_q1"] = int(rq1.sum())
    rec["n_resid_q4"] = int(rq4.sum())

    cd8_q1, cd8_q3 = quartiles(cd8)
    cd8_high = np.isfinite(cd8) & (cd8 >= cd8_q3)
    rec["n_cd8a_high"] = int(cd8_high.sum())
    rec["cd8a_q4_cut"] = float(cd8_q3)

    # hex ring-1 neighbor CD8A
    coord = {
        (int(r), int(c)): i
        for i, (r, c) in enumerate(zip(rows, cols))
        if np.isfinite(r) and np.isfinite(c)
    }
    nei_cd8 = np.full(len(bc_k), np.nan)
    for (r, c), i in coord.items():
        acc = [cd8[j] for nb in hex_neighbors(r, c) if (j := coord.get(nb)) is not None]
        if len(acc) >= MIN_NEI:
            nei_cd8[i] = float(np.nanmean(acc))
    rec["CLDN4_vs_neiCD8A"] = spearman(cl4, nei_cd8)
    rec["partial_CLDN4_vs_neiCD8A_ctrl_KRT8"] = partial_spearman(cl4, nei_cd8, krt8)
    rec["q4q1_neighbor_CD8A"] = mw_q4_q1(nei_cd8[q4], nei_cd8[q1])
    rec["resid_q4q1_neighbor_CD8A"] = mw_q4_q1(nei_cd8[rq4], nei_cd8[rq1])

    # nearest CD8A-high distance (exclude self)
    xs, ys = hex_xy_um(rows, cols)
    xy = np.column_stack([xs, ys])
    hi_idx = np.where(cd8_high)[0]
    dist = np.full(len(bc_k), np.nan)
    if hi_idx.size >= 1:
        tree = cKDTree(xy[hi_idx])
        kq = 2 if hi_idx.size >= 2 else 1
        dd, _jj = tree.query(xy, k=kq)
        if kq == 1:
            dist = np.where(cd8_high, np.nan, np.asarray(dd, dtype=float))
        else:
            dist = np.where(cd8_high, dd[:, 1], dd[:, 0])
    rec["q4q1_nearest_CD8Ahigh_um"] = mw_q4_q1(dist[q4], dist[q1])
    rec["resid_q4q1_nearest_CD8Ahigh_um"] = mw_q4_q1(dist[rq4], dist[rq1])

    rec["spot_table"] = pd.DataFrame(
        {
            "barcode": bc_k,
            "array_row": rows,
            "array_col": cols,
            "x_um": xs,
            "y_um": ys,
            "CLDN4": cl4,
            "CD8A": cd8,
            "KRT8": krt8,
            "epi": epi,
            "CLDN4_resid_KRT8": cl4_krt8_resid,
            "nei_CD8A": nei_cd8,
            "dist_CD8Ahigh_um": dist,
            "epithelial_like": epi_like.astype(int),
            "cldn4_q": np.where(q4, "Q4", np.where(q1, "Q1", "")),
            "cd8a_high": cd8_high.astype(int),
        }
    )
    return rec


def plot_section(rec, dest):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = rec["spot_table"]
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 9.2), constrained_layout=True)
    fig.suptitle(
        f"{rec['section']}  ({rec['source']}; {rec['histology']}; n={rec['n_qc']} QC spots)",
        fontsize=11,
    )

    def scatter(ax, val, title, cmap="viridis"):
        sc = ax.scatter(df["x_um"], df["y_um"], c=val, s=6, cmap=cmap, linewidths=0)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title, fontsize=10)
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)

    scatter(axes[0, 0], df["CLDN4"], "CLDN4 log1p CP10K", "magma")
    scatter(axes[0, 1], df["CD8A"], "CD8A log1p CP10K", "cividis")
    scatter(axes[1, 0], df["CLDN4_resid_KRT8"], "CLDN4 residual | KRT8", "coolwarm")

    ax = axes[1, 1]
    ax.scatter(df["x_um"], df["y_um"], c="#d9d9d9", s=5, linewidths=0, label="_bg")
    m4 = df["cldn4_q"] == "Q4"
    m1 = df["cldn4_q"] == "Q1"
    hi = df["cd8a_high"] == 1
    ax.scatter(df.loc[hi, "x_um"], df.loc[hi, "y_um"], c="#1f77b4", s=8, linewidths=0, label="CD8A-high (Q4)")
    ax.scatter(df.loc[m1, "x_um"], df.loc[m1, "y_um"], c="#2ca02c", s=10, linewidths=0, label="epi CLDN4 Q1")
    ax.scatter(df.loc[m4, "x_um"], df.loc[m4, "y_um"], c="#d62728", s=10, linewidths=0, label="epi CLDN4 Q4")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Epithelial-like CLDN4 Q4/Q1 vs CD8A-high", fontsize=10)
    ax.legend(loc="upper right", fontsize=7, frameon=False, markerscale=1.6)

    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=160)
    plt.close(fig)


def section_row(rec):
    s = rec["same_CLDN4_CD8A"]
    p = rec["partial_CLDN4_CD8A_ctrl_KRT8"]
    n1 = rec["q4q1_nearest_CD8Ahigh_um"]
    nb = rec["q4q1_neighbor_CD8A"]
    rn = rec["resid_q4q1_nearest_CD8Ahigh_um"]
    rb = rec["resid_q4q1_neighbor_CD8A"]
    return {
        "source": rec["source"],
        "section": rec["section"],
        "histology": rec["histology"],
        "n_raw": rec["n_raw"],
        "n_qc": rec["n_qc"],
        "n_epithelial_like": rec["n_epithelial_like"],
        "n_cldn4_q4": rec["n_cldn4_q4"],
        "n_cldn4_q1": rec["n_cldn4_q1"],
        "n_cd8a_high": rec["n_cd8a_high"],
        "spearman_CLDN4_CD8A_rho": s["rho"],
        "spearman_CLDN4_CD8A_p": s["p"],
        "spearman_CLDN4_CD8A_n": s["n"],
        "partial_CLDN4_CD8A_ctrl_KRT8_rho": p["rho"],
        "partial_CLDN4_CD8A_ctrl_KRT8_p": p["p"],
        "nearest_um_median_Q4": n1["median_q4"],
        "nearest_um_median_Q1": n1["median_q1"],
        "nearest_um_delta_Q4_minus_Q1": n1["delta_q4_minus_q1"],
        "nearest_um_p": n1["p"],
        "neighbor_CD8A_median_Q4": nb["median_q4"],
        "neighbor_CD8A_median_Q1": nb["median_q1"],
        "neighbor_CD8A_delta_Q4_minus_Q1": nb["delta_q4_minus_q1"],
        "neighbor_CD8A_p": nb["p"],
        "resid_nearest_um_delta_Q4_minus_Q1": rn["delta_q4_minus_q1"],
        "resid_nearest_um_p": rn["p"],
        "resid_neighbor_CD8A_delta_Q4_minus_Q1": rb["delta_q4_minus_q1"],
        "resid_neighbor_CD8A_p": rb["p"],
        "CLDN4_vs_neiCD8A_rho": rec["CLDN4_vs_neiCD8A"]["rho"],
        "CLDN4_vs_neiCD8A_p": rec["CLDN4_vs_neiCD8A"]["p"],
    }


def write_results_md(rows, skips, status):
    lines = []
    lines.append("# RESULTS — additive public Visium, CLDN4-only")
    lines.append("")
    lines.append("Public processed Visium only. **CLDN4-only** (no dual-high / TACSTD2 gate). No private 8-KL. No invented accessions or statistics.")
    lines.append("")
    lines.append("## Sources")
    lines.append("")
    lines.append("| ID | Accession | Design | Slide chemistry used for coordinates | Status |")
    lines.append("|---|---|---|---|---|")
    a = status["GSE277206"]
    b = status["zenodo_13337961"]
    lines.append(
        f"| A | [GSE277206](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277206) | 10x Visium CytAssist FFPE, never-smoker LUAD progression (MIA-034, MIA-039; small n) | visium-v5 / 11 mm (barcode whitelist; GEO uploaded H5 only) | {a} |"
    )
    lines.append(
        f"| B | [Zenodo 13337961](https://zenodo.org/records/13337961) / [10.3389/fimmu.2024.1430163](https://doi.org/10.3389/fimmu.2024.1430163) | lepidic vs solid LUAD spatial transcriptome | visium-v4 / CytAssist 6.5 mm (record text says 11 mm; uploaded barcodes are the 4,992-spot 6.5 mm set) | {b} |"
    )
    lines.append("")
    lines.append("## Locked methods")
    lines.append("")
    lines.append("- QC: spots with ≥200 detected genes and a matched Visium array coordinate.")
    lines.append("- Expression: log1p(CP10K).")
    lines.append("- Same-spot test: Spearman **CLDN4 vs CD8A** on all QC spots.")
    lines.append("- Epithelial-like: section Q3+ of the mean of available `EPCAM`, `KRT8`, `KRT18`, `KRT19`, `CDH1`, `KRT7` (KRT18 is absent from the Visium Human Transcriptome Probe Set v2.0 matrices used here).")
    lines.append("- CLDN4-high / low: Q4 / Q1 of CLDN4 **among epithelial-like spots**.")
    lines.append("- CD8A-high: section Q4 of CD8A among all QC spots.")
    lines.append("- Nearest CD8A-high distance: hex-aware Euclidean distance in µm (100 µm Visium pitch). Self is excluded when a query spot is itself CD8A-high.")
    lines.append("- Neighbor CD8A: mean CD8A of hex ring-1 neighbors `(row, col±2)` and `(row±1, col±1)`; require ≥3 neighbors.")
    lines.append("- KRT8 residualization: (i) partial Spearman of CLDN4 vs CD8A controlling for KRT8 ranks; (ii) linear residual of log CLDN4 on log KRT8, then the same Q4/Q1 distance and neighbor tests on residual quartiles among epithelial-like spots.")
    lines.append("- Q4 vs Q1: two-sided Mann–Whitney U. Section is the unit; n=2 per source is too small for a signed-rank across sections.")
    lines.append("- Maps use array coordinates recovered from the public 10x barcode inclusion list. Section-specific H&E / `tissue_positions` pixel maps were **not** deposited with these count matrices.")
    lines.append("")
    lines.append("## Per-section results (computed)")
    lines.append("")
    if not rows:
        lines.append("No section passed CLDN4 + CD8A + QC.")
    else:
        lines.append("| Source | Section | Histology | n QC | n epi-like | Spearman CLDN4 vs CD8A | Partial \\| KRT8 | Nearest CD8A-high µm (Q4 vs Q1) | Neighbor CD8A (Q4 vs Q1) |")
        lines.append("|---|---|---|---:|---:|---|---|---|---|")
        for r in rows:
            s = f"ρ={fmt_rho(r['spearman_CLDN4_CD8A_rho'])}, p={fmt_p(r['spearman_CLDN4_CD8A_p'])} (n={r['spearman_CLDN4_CD8A_n']})"
            pr = f"ρ={fmt_rho(r['partial_CLDN4_CD8A_ctrl_KRT8_rho'])}, p={fmt_p(r['partial_CLDN4_CD8A_ctrl_KRT8_p'])}"
            nd = (
                f"med {r['nearest_um_median_Q4']:.1f} vs {r['nearest_um_median_Q1']:.1f} "
                f"(Δ={r['nearest_um_delta_Q4_minus_Q1']:.1f}, p={fmt_p(r['nearest_um_p'])})"
                if r["nearest_um_median_Q4"] is not None and r["nearest_um_median_Q1"] is not None
                else "NA"
            )
            nb = (
                f"med {r['neighbor_CD8A_median_Q4']:.3f} vs {r['neighbor_CD8A_median_Q1']:.3f} "
                f"(Δ={r['neighbor_CD8A_delta_Q4_minus_Q1']:.3f}, p={fmt_p(r['neighbor_CD8A_p'])})"
                if r["neighbor_CD8A_median_Q4"] is not None and r["neighbor_CD8A_median_Q1"] is not None
                else "NA"
            )
            lines.append(
                f"| {r['source']} | {r['section']} | {r['histology']} | {r['n_qc']} | {r['n_epithelial_like']} | {s} | {pr} | {nd} | {nb} |"
            )
    lines.append("")
    lines.append("### KRT8-residual Q4 vs Q1 (epithelial-like residual quartiles)")
    lines.append("")
    if rows:
        lines.append("| Source | Section | Resid nearest µm Δ (Q4−Q1) | p | Resid neighbor CD8A Δ (Q4−Q1) | p |")
        lines.append("|---|---|---:|---:|---:|---:|")
        for r in rows:
            rd = r["resid_nearest_um_delta_Q4_minus_Q1"]
            rb = r["resid_neighbor_CD8A_delta_Q4_minus_Q1"]
            rd_s = "NA" if rd is None else f"{rd:.1f}"
            rb_s = "NA" if rb is None else f"{rb:.3f}"
            lines.append(
                f"| {r['source']} | {r['section']} | {rd_s} | {fmt_p(r['resid_nearest_um_p'])} | {rb_s} | {fmt_p(r['resid_neighbor_CD8A_p'])} |"
            )
    lines.append("")
    if skips:
        lines.append("## Skips")
        lines.append("")
        for s in skips:
            lines.append(f"- `{s.get('section', '?')}`: {s.get('skip')} {s.get('note', '')}".rstrip())
        lines.append("")
    lines.append("## What this does not claim")
    lines.append("")
    lines.append("- It does not use private 8-KL or any controlled/unpublished matrix.")
    lines.append("- It does not treat n=2 sections as a multi-cohort meta-analysis.")
    lines.append("- Same-spot CLDN4 vs CD8A is composition-sensitive; the KRT8 residual / partial is the control, not a proof of a CLDN4-specific barrier.")
    lines.append("- Array maps are barcode-whitelist coordinates, not the unpublished section H&E.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 methods/visium_public_luad/download.py")
    lines.append("python3 methods/visium_public_luad/analyze.py")
    lines.append("```")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `methods/visium_public_luad/analyze.py`")
    lines.append("- `methods/visium_public_luad/download.py`")
    lines.append("- `methods/visium_public_luad/coords/visium_v4_cytassist_6p5mm.csv`")
    lines.append("- `methods/visium_public_luad/coords/visium_v5_cytassist_11mm.csv`")
    lines.append("- `results/visium_public_luad/tables/section_stats.tsv`")
    lines.append("- `results/visium_public_luad/tables/summary.json`")
    lines.append("- `results/visium_public_luad/maps/`")
    lines.append("")
    (OUT / "RESULTS.md").write_text("\n".join(lines) + "\n")


def run_gse277206():
    d = DATA / "GSE277206"
    pos = load_coord_map("v5_11mm")
    specs = [
        ("GSE277206_MIA-034", "MIA", d / "GSM8516529_V0801_filtered_feature_bc_matrix.h5"),
        ("GSE277206_MIA-039", "MIA", d / "GSM8516530_V0901_filtered_feature_bc_matrix.h5"),
    ]
    if not d.exists() or not all(p.exists() for _, _, p in specs):
        return [], [{"section": "GSE277206", "skip": "not_downloaded", "note": str(d)}]
    out, skips = [], []
    for name, hist, path in specs:
        print("GSE277206", name, flush=True)
        barcodes, genes, mtx = load_10x_h5(path)
        rec = analyze_section(name, "GSE277206", hist, barcodes, genes, mtx, pos)
        if "skip" in rec:
            skips.append(rec)
            continue
        plot_section(rec, MAPS / f"{name}.png")
        rec["spot_table"].to_csv(TABLES / f"{name}_spots.tsv.gz", sep="\t", index=False)
        rec = {k: v for k, v in rec.items() if k != "spot_table"}
        out.append(rec)
    return out, skips


def run_zenodo():
    d = DATA / "zenodo_13337961"
    pos = load_coord_map("v4_6p5mm")
    specs = [
        ("Zenodo13337961_lepidic", "lepidic", d / "lepidic" / "filtered_feature_bc_matrix"),
        ("Zenodo13337961_solid", "solid", d / "solid" / "filtered_feature_bc_matrix"),
    ]
    if not all((p / "matrix.mtx.gz").exists() for _, _, p in specs):
        return [], [{"section": "zenodo_13337961", "skip": "not_downloaded", "note": str(d)}]
    out, skips = [], []
    for name, hist, path in specs:
        print("Zenodo", name, flush=True)
        barcodes, genes, mtx = load_mtx_dir(path)
        rec = analyze_section(name, "zenodo_13337961", hist, barcodes, genes, mtx, pos)
        if "skip" in rec:
            skips.append(rec)
            continue
        plot_section(rec, MAPS / f"{name}.png")
        rec["spot_table"].to_csv(TABLES / f"{name}_spots.tsv.gz", sep="\t", index=False)
        rec = {k: v for k, v in rec.items() if k != "spot_table"}
        out.append(rec)
    return out, skips


def main():
    # Prefer already-fetched /tmp copies if present (cloud workspace).
    gse_src = Path("/tmp/visium_dl/gse277206")
    zen_src = Path("/tmp/visium_dl/zenodo")
    gse_dst = DATA / "GSE277206"
    zen_dst = DATA / "zenodo_13337961"
    if gse_src.exists() and not (gse_dst / "GSM8516529_V0801_filtered_feature_bc_matrix.h5").exists():
        gse_dst.mkdir(parents=True, exist_ok=True)
        for p in gse_src.glob("*.h5"):
            dest = gse_dst / p.name
            if not dest.exists():
                dest.symlink_to(p)
    if (zen_src / "lepidic" / "filtered_feature_bc_matrix" / "matrix.mtx.gz").exists():
        for label in ("lepidic", "solid"):
            src = zen_src / label / "filtered_feature_bc_matrix"
            dst = zen_dst / label / "filtered_feature_bc_matrix"
            if src.exists() and not (dst / "matrix.mtx.gz").exists():
                dst.parent.mkdir(parents=True, exist_ok=True)
                if dst.exists() or dst.is_symlink():
                    pass
                else:
                    dst.symlink_to(src)

    gse_recs, gse_skips = run_gse277206()
    zen_recs, zen_skips = run_zenodo()
    recs = gse_recs + zen_recs
    skips = gse_skips + zen_skips
    rows = [section_row(r) for r in recs]
    pd.DataFrame(rows).to_csv(TABLES / "section_stats.tsv", sep="\t", index=False)
    status = {
        "GSE277206": "ran" if gse_recs else ("missing/controlled" if gse_skips else "empty"),
        "zenodo_13337961": "ran" if zen_recs else ("missing/controlled" if zen_skips else "empty"),
    }
    if gse_recs:
        status["GSE277206"] = f"ran ({len(gse_recs)} sections)"
    if zen_recs:
        status["zenodo_13337961"] = f"ran ({len(zen_recs)} sections)"
    dump_json(
        TABLES / "summary.json",
        {"status": status, "sections": recs, "skips": skips, "methods": {
            "min_genes": MIN_GENES,
            "pitch_um": PITCH_UM,
            "epithelial_like": "section Q3+ of mean(EPCAM,KRT8,KRT18,KRT19,CDH1,KRT7 if present)",
            "cldn4_high": "Q4 of CLDN4 among epithelial-like",
            "cd8a_high": "Q4 of CD8A among QC spots",
            "krt8_residual": "linear residual of log CLDN4 on log KRT8; partial Spearman ranks",
        }},
    )
    write_results_md(rows, skips, status)
    print("wrote", OUT / "RESULTS.md")
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
