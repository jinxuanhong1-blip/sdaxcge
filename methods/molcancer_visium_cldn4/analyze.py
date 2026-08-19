#!/usr/bin/env python3
"""CLDN4-only Visium stats for Molecular Cancer 2025 post-chemoIO (PRJNA1139087)
and open pre-ICI NSCLC Visium GSE292299.

Question: are CLDN4-high tumor spots farther from CD8 / lower CD8 neighbor?
If a pre cohort exists, does post-ICI CLDN4-high still exclude?

Locked metrics: spot Spearman, nearest CD8A distance, hex-radius neighbors,
KRT8 residual. Section is the unit. No TACSTD2 primary. No private 8-KL.
"""

from __future__ import annotations

import json
import tarfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
OUT = ROOT / "results" / "molcancer_visium_cldn4"
TABLES = OUT / "tables"
FIGS = OUT / "maps"
TABLES.mkdir(parents=True, exist_ok=True)
FIGS.mkdir(parents=True, exist_ok=True)

BROAD_EPI = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]
MIN_N_SPEARMAN = 8
MIN_NEI = 3
# Visium center-to-center is 100 um; hex uses (row, col±2) / (row±1, col±1).
UM_PER_COL = 50.0  # col step 2 = 100 um
UM_PER_ROW = 100.0 * np.sqrt(3) / 2.0

META_PRJ = {
    "PA08": {"response": "NMPR", "timepoint": "post_chemoIO", "note": "cold / unresponsive"},
    "PA09": {"response": "NMPR", "timepoint": "post_chemoIO", "note": "cold / unresponsive"},
    "PA10": {"response": "pCR", "timepoint": "post_chemoIO", "note": "pathologic complete response"},
    "PA12": {"response": "NMPR_responsive_like", "timepoint": "post_chemoIO", "note": "NMPR but therapy-responsive histology"},
}


def dump_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, default=_json_default) + "\n")


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < MIN_N_SPEARMAN:
        return {"n": n, "rho": None, "p": None}
    rho, p = stats.spearmanr(x[m], y[m])
    if not np.isfinite(rho):
        return {"n": n, "rho": None, "p": None}
    return {"n": n, "rho": float(rho), "p": float(p)}


def residual_vs(x, z):
    """Rank-residual of x after linear fit on ranks of z (KRT8 control)."""
    x = np.asarray(x, dtype=float)
    z = np.asarray(z, dtype=float)
    m = np.isfinite(x) & np.isfinite(z)
    out = np.full(x.shape, np.nan)
    if m.sum() < 12:
        return out
    rx = stats.rankdata(x[m])
    rz = stats.rankdata(z[m])
    rz_s = (rz - rz.mean()) / (rz.std() + 1e-12)
    b = np.polyfit(rz_s, rx, 1)
    out[m] = rx - (b[0] * rz_s + b[1])
    return out


def partial_spearman(x, y, z):
    ex = residual_vs(x, z)
    ey = residual_vs(y, z)
    return spearman(ex, ey)


def wilcoxon_signed(values):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    n = int(v.size)
    out = {
        "n": n,
        "median": float(np.median(v)) if n else None,
        "n_neg": int((v < 0).sum()) if n else 0,
        "n_pos": int((v > 0).sum()) if n else 0,
        "p": None,
    }
    if n < 3:
        return out
    try:
        out["p"] = float(stats.wilcoxon(v, alternative="two-sided").pvalue)
    except ValueError:
        out["p"] = None
    return out


def mwu(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 5 or len(b) < 5:
        return {
            "n_high": int(len(a)),
            "n_low": int(len(b)),
            "median_high": float(np.median(a)) if len(a) else None,
            "median_low": float(np.median(b)) if len(b) else None,
            "p": None,
        }
    p = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    return {
        "n_high": int(len(a)),
        "n_low": int(len(b)),
        "median_high": float(np.median(a)),
        "median_low": float(np.median(b)),
        "delta_high_minus_low": float(np.median(a) - np.median(b)),
        "p": p,
    }


def hex_neighbors(row, col):
    return [
        (row, col - 2),
        (row, col + 2),
        (row - 1, col - 1),
        (row - 1, col + 1),
        (row + 1, col - 1),
        (row + 1, col + 1),
    ]


def xy_um(row, col):
    return np.asarray(col, dtype=float) * UM_PER_COL, np.asarray(row, dtype=float) * UM_PER_ROW


def log1p_cp10k_from_reads(count, n_reads):
    lib = np.asarray(n_reads, dtype=float)
    lib[lib <= 0] = 1.0
    return np.log1p(np.asarray(count, dtype=float) * (1e4 / lib))


def log1p_cp10k_matrix(counts, lib):
    lib = np.asarray(lib, dtype=float)
    lib[lib <= 0] = 1.0
    return np.log1p(np.asarray(counts, dtype=float) * (1e4 / lib[:, None]))


def neighbor_mean(coord, values, k):
    out = np.full(len(values), np.nan)
    for (r, c), i in coord.items():
        if k == 1:
            nbs = hex_neighbors(r, c)
        else:
            seen = {(r, c)}
            frontier = {(r, c)}
            for _ in range(k):
                nxt = set()
                for rr, cc in frontier:
                    for nb in hex_neighbors(rr, cc):
                        if nb not in seen:
                            nxt.add(nb)
                            seen.add(nb)
                frontier = nxt
            nbs = list(frontier)
        acc = [values[j] for nb in nbs if (j := coord.get(nb)) is not None]
        if len(acc) >= MIN_NEI:
            out[i] = float(np.nanmean(acc))
    return out


def nearest_cd8_um(row, col, cd8_pos, cd8_xy):
    x, y = xy_um(row, col)
    pts = np.column_stack([x, y])
    if cd8_xy.shape[0] == 0:
        return np.full(len(row), np.nan)
    tree = cKDTree(cd8_xy)
    dist, _ = tree.query(pts, k=1)
    # self-hit: if a tumor spot is also CD8+, ignore distance 0 by using 2nd neighbor when possible
    if cd8_xy.shape[0] >= 2:
        dist2, _ = tree.query(pts, k=2)
        is_cd8 = cd8_pos
        dist = np.where(is_cd8, dist2[:, 1], dist)
    return dist


def score_cols(df, names):
    have = [g for g in names if g in df.columns]
    if not have:
        return np.full(len(df), np.nan)
    return df[have].mean(axis=1).to_numpy()


def analyze_section(name, df, timepoint, response, lib_col="n_reads"):
    row = df["array_row"].to_numpy(dtype=float)
    col = df["array_col"].to_numpy(dtype=float)
    lib = df[lib_col].to_numpy(dtype=float)
    cl4 = log1p_cp10k_from_reads(df["CLDN4"], lib)
    cd8 = log1p_cp10k_from_reads(df["CD8A"], lib)
    krt8 = log1p_cp10k_from_reads(df["KRT8"], lib) if "KRT8" in df.columns else np.full(len(df), np.nan)
    epi = np.nanmean(
        np.vstack(
            [
                log1p_cp10k_from_reads(df[g], lib)
                for g in BROAD_EPI
                if g in df.columns
            ]
        ),
        axis=0,
    ) if any(g in df.columns for g in BROAD_EPI) else krt8

    # Tissue already filtered. Tumor = epithelial-high (top half of epi among spots).
    epi_ok = np.isfinite(epi)
    epi_cut = np.nanmedian(epi[epi_ok]) if epi_ok.any() else np.nan
    tumor = epi_ok & (epi >= epi_cut)
    # CD8+ spots: CD8A UMI > 0
    raw_cd8 = df["CD8A"].to_numpy(dtype=float)
    cd8_pos = raw_cd8 > 0
    # CLDN4-high / low among tumor spots: top vs bottom quartile
    cl4_t = cl4[tumor]
    if cl4_t.size >= 20:
        q1, q3 = np.nanpercentile(cl4_t, [25, 75])
        high = tumor & (cl4 >= q3)
        low = tumor & (cl4 <= q1)
    else:
        high = tumor & False
        low = tumor & False

    coord = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(row, col))}
    nei1 = neighbor_mean(coord, cd8, 1)
    nei2 = neighbor_mean(coord, cd8, 2)
    nei3 = neighbor_mean(coord, cd8, 3)
    x_cd8, y_cd8 = xy_um(row[cd8_pos], col[cd8_pos])
    cd8_xy = np.column_stack([x_cd8, y_cd8]) if cd8_pos.any() else np.zeros((0, 2))
    dist = nearest_cd8_um(row, col, cd8_pos, cd8_xy)
    cl4_krt8_resid = residual_vs(cl4, krt8)

    rec = {
        "section": name,
        "timepoint": timepoint,
        "response": response,
        "n_spots": int(len(df)),
        "n_tumor": int(tumor.sum()),
        "n_cd8pos": int(cd8_pos.sum()),
        "n_cldn4_high_tumor": int(high.sum()),
        "n_cldn4_low_tumor": int(low.sum()),
        "frac_cd8pos": float(cd8_pos.mean()),
        "mean_CLDN4_log": float(np.nanmean(cl4)),
        "mean_CD8A_log": float(np.nanmean(cd8)),
        "same_CLDN4_CD8A_all": spearman(cl4, cd8),
        "same_CLDN4_CD8A_tumor": spearman(cl4[tumor], cd8[tumor]),
        "same_KRT8_CD8A_tumor": spearman(krt8[tumor], cd8[tumor]),
        "partial_CLDN4_CD8A_ctrl_KRT8_all": partial_spearman(cl4, cd8, krt8),
        "partial_CLDN4_CD8A_ctrl_KRT8_tumor": partial_spearman(cl4[tumor], cd8[tumor], krt8[tumor]),
        "spearman_CLDN4_vs_nnCD8_tumor": spearman(cl4[tumor], dist[tumor]),
        "spearman_CLDN4residKRT8_vs_nnCD8_tumor": spearman(cl4_krt8_resid[tumor], dist[tumor]),
        "spearman_CLDN4_vs_ring1CD8_tumor": spearman(cl4[tumor], nei1[tumor]),
        "spearman_CLDN4_vs_ring2CD8_tumor": spearman(cl4[tumor], nei2[tumor]),
        "spearman_CLDN4_vs_ring3CD8_tumor": spearman(cl4[tumor], nei3[tumor]),
        "partial_CLDN4_vs_ring1CD8_ctrl_KRT8_tumor": partial_spearman(cl4[tumor], nei1[tumor], krt8[tumor]),
        "nnCD8_high_vs_low_tumor_um": mwu(dist[high], dist[low]),
        "ring1CD8_high_vs_low_tumor": mwu(nei1[high], nei1[low]),
        "ring2CD8_high_vs_low_tumor": mwu(nei2[high], nei2[low]),
    }
    maps = {
        "section": name,
        "row": row,
        "col": col,
        "cl4": cl4,
        "cd8": cd8,
        "krt8": krt8,
        "tumor": tumor,
        "high": high,
        "dist": dist,
        "nei1": nei1,
    }
    return rec, maps


def plot_maps(maps, out_png):
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    row, col = maps["row"], maps["col"]
    specs = [
        (axes[0, 0], maps["cl4"], "CLDN4 (log1p CP10k-proxy)", "viridis"),
        (axes[0, 1], maps["cd8"], "CD8A (log1p CP10k-proxy)", "magma"),
        (axes[0, 2], maps["krt8"], "KRT8 (log1p CP10k-proxy)", "cividis"),
        (axes[1, 0], maps["tumor"].astype(float), "Tumor (epi ≥ median)", "coolwarm"),
        (axes[1, 1], maps["dist"], "Nearest CD8A+ (µm)", "plasma"),
        (axes[1, 2], maps["nei1"], "Ring-1 neighbor CD8A", "magma"),
    ]
    for ax, val, title, cmap in specs:
        sc = ax.scatter(col, -row, c=val, s=6, cmap=cmap, linewidths=0)
        ax.set_title(title, fontsize=9)
        ax.set_aspect("equal")
        ax.axis("off")
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    # overlay CLDN4-high tumor
    h = maps["high"]
    axes[1, 0].scatter(col[h], -row[h], s=10, facecolors="none", edgecolors="k", linewidths=0.4)
    fig.suptitle(maps["section"], fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    fig.savefig(str(out_png).replace(".png", ".pdf"))
    plt.close(fig)


def summarize(rows, key):
    rhos = []
    for r in rows:
        v = r.get(key) or {}
        if isinstance(v, dict) and v.get("rho") is not None:
            rhos.append(v["rho"])
    return wilcoxon_signed(rhos)


def summarize_delta(rows, key, field="delta_high_minus_low"):
    vals = []
    for r in rows:
        v = r.get(key) or {}
        if isinstance(v, dict) and v.get(field) is not None:
            vals.append(v[field])
    return wilcoxon_signed(vals)


def run_prjna():
    rows = []
    count_dir = DATA / "PRJNA1139087" / "counts"
    for patient, meta in META_PRJ.items():
        path = count_dir / f"{patient}_spot_counts.csv"
        if not path.exists():
            rows.append({"section": patient, "skip": "counts_missing"})
            continue
        df = pd.read_csv(path)
        # tissue: keep barcodes with enough total reads (knee ~ paper ~3k/section)
        n = df["n_reads"].to_numpy()
        if len(df) == 0:
            rows.append({"section": patient, "skip": "empty"})
            continue
        cut = max(200, float(np.quantile(n, 0.35)))
        df = df[df["n_reads"] >= cut].copy()
        rec, maps = analyze_section(patient, df, meta["timepoint"], meta["response"])
        rec["n_reads_cut"] = float(cut)
        rec["note"] = meta["note"]
        plot_maps(maps, FIGS / f"{patient}_maps.png")
        rows.append(rec)
    ok = [r for r in rows if "skip" not in r]
    summary = {
        "cohort": "PRJNA1139087_MolCancer2025",
        "timepoint": "post_neoadjuvant_chemoIO",
        "n_sections": len(ok),
        "n_spots": int(sum(r["n_spots"] for r in ok)),
        "same_CLDN4_CD8A_tumor": summarize(ok, "same_CLDN4_CD8A_tumor"),
        "partial_CLDN4_CD8A_ctrl_KRT8_tumor": summarize(ok, "partial_CLDN4_CD8A_ctrl_KRT8_tumor"),
        "spearman_CLDN4_vs_nnCD8_tumor": summarize(ok, "spearman_CLDN4_vs_nnCD8_tumor"),
        "spearman_CLDN4_vs_ring1CD8_tumor": summarize(ok, "spearman_CLDN4_vs_ring1CD8_tumor"),
        "partial_CLDN4_vs_ring1CD8_ctrl_KRT8_tumor": summarize(ok, "partial_CLDN4_vs_ring1CD8_ctrl_KRT8_tumor"),
        "nnCD8_high_minus_low_um": summarize_delta(ok, "nnCD8_high_vs_low_tumor_um"),
        "ring1CD8_high_minus_low": summarize_delta(ok, "ring1CD8_high_vs_low_tumor"),
        "sections": ok,
        "skips": [r for r in rows if "skip" in r],
    }
    pd.json_normalize(ok).to_csv(TABLES / "prjna1139087_sections.csv", index=False)
    dump_json(TABLES / "prjna1139087_summary.json", summary)
    return summary


def load_10x_h5(path):
    import h5py
    from scipy.sparse import csr_matrix

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


def read_positions_csv(path):
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    if "barcode" not in cols:
        df = pd.read_csv(path, header=None)
        if str(df.iloc[0, 0]).lower() in {"barcode", "barcodes"}:
            df = df.iloc[1:].reset_index(drop=True)
        df.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"][: df.shape[1]]
    else:
        rename = {}
        for want in ("barcode", "in_tissue", "array_row", "array_col"):
            if want in cols:
                rename[cols[want]] = want
        df = df.rename(columns=rename)
    df["barcode"] = df["barcode"].astype(str)
    return df.set_index("barcode")


def extract_tissue_positions(tar_path, dest):
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    out = dest / "tissue_positions.csv"
    if out.exists():
        return out
    with tarfile.open(tar_path) as t:
        members = [m for m in t.getmembers() if m.isfile() and "tissue_positions" in Path(m.name).name.lower()]
        if not members:
            raise FileNotFoundError(f"no tissue_positions in {tar_path}")
        members.sort(key=lambda m: m.size)
        out.write_bytes(t.extractfile(members[-1]).read())
    return out


def run_gse292299():
    d = DATA / "GSE292299"
    if not d.exists():
        return {"cohort": "GSE292299", "skip": "data_missing"}
    pos_dir = d / "positions"
    rows = []
    for i in range(1, 17):
        h5s = list(d.glob(f"GSM*_NSCLC_P{i}_filtered_feature_bc_matrix.h5"))
        spas = list(d.glob(f"GSM*_NSCLC_P{i}_spatial.tar.gz"))
        name = f"NSCLC_P{i}"
        if not h5s or not spas:
            rows.append({"section": name, "skip": "missing_files"})
            continue
        barcodes, genes, mtx = load_10x_h5(h5s[0])
        # First Gene Expression / ENSG hit wins. GSE292299 H5 also has antibody
        # features that reuse gene symbols (empty CD8A protein column).
        gidx = {}
        for j, g in enumerate(genes):
            gidx.setdefault(str(g).upper(), j)
        need = ["CLDN4", "CD8A", "KRT8"]
        if any(g not in gidx for g in need):
            rows.append({"section": name, "skip": "missing_gene", "have": [g for g in need if g in gidx]})
            continue
        pos = read_positions_csv(extract_tissue_positions(spas[0], pos_dir / name))
        n_genes = np.asarray((mtx > 0).sum(axis=1)).ravel()
        in_t = pos.reindex(barcodes)["in_tissue"].fillna(0).astype(int).to_numpy()
        keep = (in_t == 1) & (n_genes >= 200)
        if keep.sum() < 50:
            rows.append({"section": name, "skip": "too_few_spots", "n_keep": int(keep.sum())})
            continue
        mtx_k = mtx[keep]
        bc_k = np.asarray(barcodes)[keep]
        pos_k = pos.reindex(bc_k)
        lib = np.asarray(mtx_k.sum(axis=1)).ravel()
        df = pd.DataFrame(
            {
                "barcode": bc_k,
                "array_row": pos_k["array_row"].to_numpy(),
                "array_col": pos_k["array_col"].to_numpy(),
                "n_reads": lib,
                "CLDN4": np.asarray(mtx_k[:, gidx["CLDN4"]].todense()).ravel(),
                "CD8A": np.asarray(mtx_k[:, gidx["CD8A"]].todense()).ravel(),
                "KRT8": np.asarray(mtx_k[:, gidx["KRT8"]].todense()).ravel(),
            }
        )
        for g in BROAD_EPI:
            if g in gidx and g not in df.columns:
                df[g] = np.asarray(mtx_k[:, gidx[g]].todense()).ravel()
        rec, maps = analyze_section(name, df, "pre_ICI_biopsy", "pre_ICI_NSCLC")
        if i <= 4:
            plot_maps(maps, FIGS / f"GSE292299_{name}_maps.png")
        rows.append(rec)
        print("GSE292299", name, rec["n_spots"], flush=True)
    ok = [r for r in rows if "skip" not in r]
    summary = {
        "cohort": "GSE292299",
        "timepoint": "pre_ICI_biopsy",
        "n_sections": len(ok),
        "n_spots": int(sum(r["n_spots"] for r in ok)),
        "same_CLDN4_CD8A_tumor": summarize(ok, "same_CLDN4_CD8A_tumor"),
        "partial_CLDN4_CD8A_ctrl_KRT8_tumor": summarize(ok, "partial_CLDN4_CD8A_ctrl_KRT8_tumor"),
        "spearman_CLDN4_vs_nnCD8_tumor": summarize(ok, "spearman_CLDN4_vs_nnCD8_tumor"),
        "spearman_CLDN4_vs_ring1CD8_tumor": summarize(ok, "spearman_CLDN4_vs_ring1CD8_tumor"),
        "partial_CLDN4_vs_ring1CD8_ctrl_KRT8_tumor": summarize(ok, "partial_CLDN4_vs_ring1CD8_ctrl_KRT8_tumor"),
        "nnCD8_high_minus_low_um": summarize_delta(ok, "nnCD8_high_vs_low_tumor_um"),
        "ring1CD8_high_minus_low": summarize_delta(ok, "ring1CD8_high_vs_low_tumor"),
        "sections": ok,
        "skips": [r for r in rows if "skip" in r],
        "note": "Open pre-treatment NSCLC Visium (immunotherapy spatial biology). Not matched to PRJNA1139087 patients.",
    }
    pd.json_normalize(ok).to_csv(TABLES / "gse292299_sections.csv", index=False)
    dump_json(TABLES / "gse292299_summary.json", summary)
    return summary


def main():
    post = run_prjna()
    pre = run_gse292299()
    dump_json(TABLES / "summary.json", {"post_PRJNA1139087": post, "pre_GSE292299": pre})
    print(json.dumps({"post": {k: post.get(k) for k in post if k != "sections"}, "pre": {k: pre.get(k) for k in pre if k != "sections"}}, indent=2, default=_json_default))


if __name__ == "__main__":
    main()
