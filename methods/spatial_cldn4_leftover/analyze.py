#!/usr/bin/env python3
"""Additive leftover public spatial: CLDN4 vs T/B neighborhood.

Skip already-reported GSE265899, GSE289483, GSE273378 (and prior B6 series).
Honest n / empty: a series is analyzed only if CLDN4 is in the public matrix.
"""

from __future__ import annotations

import json
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io, stats
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "spatial_cldn4_leftover"
OUT = ROOT / "results" / "spatial_cldn4_leftover"
TABLES = OUT / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

T_GENES = ["CD3D", "CD3E", "CD3G", "CD2", "CD8A", "CD8B", "TRAC", "CD247", "IL7R"]
B_GENES = ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "CD22"]
BROAD_EPI = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]
MIN_GENES = 200
MIN_NEI = 3
MIN_N_SPEARMAN = 8


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


def spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < MIN_N_SPEARMAN:
        return {"n": n, "rho": None, "p": None}
    rho, p = stats.spearmanr(x[m], y[m])
    return {"n": n, "rho": float(rho), "p": float(p)}


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


def hex_neighbors(row, col):
    return [(row, col - 2), (row, col + 2), (row - 1, col - 1), (row - 1, col + 1), (row + 1, col - 1), (row + 1, col + 1)]


def log_cp10k(mtx):
    lib = np.asarray(mtx.sum(axis=1)).ravel()
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


def present(names, universe):
    u = {str(x).upper() for x in universe}
    return [g for g in names if g.upper() in u]


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
        for want, opts in {
            "barcode": ["barcode"],
            "in_tissue": ["in_tissue"],
            "array_row": ["array_row"],
            "array_col": ["array_col"],
        }.items():
            for o in opts:
                if o in cols:
                    rename[cols[o]] = want
        df = df.rename(columns=rename)
    df["barcode"] = df["barcode"].astype(str)
    df["in_tissue"] = pd.to_numeric(df["in_tissue"], errors="coerce").fillna(0).astype(int)
    df["array_row"] = pd.to_numeric(df["array_row"], errors="coerce")
    df["array_col"] = pd.to_numeric(df["array_col"], errors="coerce")
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
        src = t.extractfile(members[-1])
        out.write_bytes(src.read())
    return out


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


def load_mtx(barcodes_path, features_path, mtx_path):
    barcodes = pd.read_csv(barcodes_path, header=None)[0].astype(str)
    feat = pd.read_csv(features_path, sep="\t", header=None)
    genes = feat[1].astype(str).tolist() if feat.shape[1] > 1 else feat[0].astype(str).tolist()
    mtx = io.mmread(mtx_path).tocsr()
    if mtx.shape[0] == len(barcodes) and mtx.shape[1] == len(genes):
        pass
    elif mtx.shape[1] == len(barcodes) and mtx.shape[0] == len(genes):
        mtx = mtx.T.tocsr()
    else:
        raise ValueError(f"shape mismatch mtx={mtx.shape} bc={len(barcodes)} genes={len(genes)}")
    return barcodes, genes, mtx


def analyze_visium_section(name, barcodes, genes, mtx, pos):
    gidx = gene_index(genes)
    if "CLDN4" not in gidx:
        return {"section": name, "skip": "CLDN4_absent", "n_genes": len(gidx)}
    n_genes = np.asarray((mtx > 0).sum(axis=1)).ravel()
    in_t = pos.reindex(barcodes)["in_tissue"].fillna(0).astype(int).to_numpy()
    keep = (in_t == 1) & (n_genes >= MIN_GENES)
    if int(keep.sum()) < 50:
        return {"section": name, "skip": "too_few_spots", "n_keep": int(keep.sum())}
    mtx_k = mtx[keep]
    bc_k = np.asarray(barcodes)[keep]
    pos_k = pos.reindex(bc_k)
    logx = log_cp10k(mtx_k)
    cl4 = score_from_log(logx, gidx, ["CLDN4"])
    tac = score_from_log(logx, gidx, ["TACSTD2"])
    broad = score_from_log(logx, gidx, present(BROAD_EPI, gidx))
    tsc = score_from_log(logx, gidx, present(T_GENES, gidx))
    bsc = score_from_log(logx, gidx, present(B_GENES, gidx))
    tb = np.nanmean(np.vstack([tsc, bsc]), axis=0)
    same = {
        "section": name,
        "n_spots": int(keep.sum()),
        "CLDN4": True,
        "TACSTD2": "TACSTD2" in gidx,
        "same_CLDN4_TB": spearman(cl4, tb),
        "same_CLDN4_T": spearman(cl4, tsc),
        "same_CLDN4_B": spearman(cl4, bsc),
        "same_TACSTD2_TB": spearman(tac, tb),
        "same_broad_TB": spearman(broad, tb),
        "partial_same_CLDN4_TB_ctrl_broad": partial_spearman(cl4, tb, broad),
    }
    coord = {
        (int(r), int(c)): i
        for i, (r, c) in enumerate(zip(pos_k["array_row"], pos_k["array_col"]))
        if pd.notna(r) and pd.notna(c)
    }
    rings = []
    for k in (1, 2, 3):
        nei_tb = np.full(len(bc_k), np.nan)
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
            acc = [tb[j] for nb in nbs if (j := coord.get(nb)) is not None]
            if len(acc) >= MIN_NEI:
                nei_tb[i] = float(np.nanmean(acc))
        ok = np.isfinite(nei_tb)
        rings.append(
            {
                "section": name,
                "ring": k,
                "n": int(ok.sum()),
                "CLDN4_vs_neiTB": spearman(cl4[ok], nei_tb[ok]),
                "TACSTD2_vs_neiTB": spearman(tac[ok], nei_tb[ok]),
                "partial_CLDN4_vs_neiTB_ctrl_broad": partial_spearman(cl4[ok], nei_tb[ok], broad[ok]),
            }
        )
    return {"same": same, "rings": rings}


def summarize_visium(same_rows, ring_rows, accession):
    def collect(key):
        rhos = []
        for r in same_rows:
            v = r.get(key) or {}
            if v.get("rho") is not None:
                rhos.append(v["rho"])
        return wilcoxon_signed(rhos)

    out = {
        "accession": accession,
        "n_sections": len(same_rows),
        "same_CLDN4_TB": collect("same_CLDN4_TB"),
        "same_CLDN4_T": collect("same_CLDN4_T"),
        "same_CLDN4_B": collect("same_CLDN4_B"),
        "same_TACSTD2_TB": collect("same_TACSTD2_TB"),
        "partial_same_CLDN4_TB_ctrl_broad": collect("partial_same_CLDN4_TB_ctrl_broad"),
        "rings": {},
    }
    for k in (1, 2, 3):
        rhos = []
        prhos = []
        for r in ring_rows:
            if r.get("ring") != k:
                continue
            v = r.get("CLDN4_vs_neiTB") or {}
            if v.get("rho") is not None:
                rhos.append(v["rho"])
            pv = r.get("partial_CLDN4_vs_neiTB_ctrl_broad") or {}
            if pv.get("rho") is not None:
                prhos.append(pv["rho"])
        out["rings"][str(k)] = {
            "CLDN4_vs_neiTB": wilcoxon_signed(rhos),
            "partial_CLDN4_vs_neiTB_ctrl_broad": wilcoxon_signed(prhos),
        }
    return out


def run_gse292299():
    d = DATA / "GSE292299"
    pos_dir = d / "positions"
    same_rows, ring_rows, skips = [], [], []
    for i in range(1, 17):
        h5s = list(d.glob(f"GSM*_NSCLC_P{i}_filtered_feature_bc_matrix.h5"))
        spas = list(d.glob(f"GSM*_NSCLC_P{i}_spatial.tar.gz"))
        name = f"NSCLC_P{i}"
        print("GSE292299", name, flush=True)
        if not h5s or not spas:
            skips.append({"section": name, "skip": "missing_files"})
            continue
        barcodes, genes, mtx = load_10x_h5(h5s[0])
        pos_path = extract_tissue_positions(spas[0], pos_dir / name)
        pos = read_positions_csv(pos_path)
        res = analyze_visium_section(name, barcodes, genes, mtx, pos)
        if "skip" in res:
            skips.append(res)
            continue
        same_rows.append(res["same"])
        ring_rows.extend(res["rings"])
    summary = summarize_visium(same_rows, ring_rows, "GSE292299")
    summary["skips"] = skips
    pd.json_normalize(same_rows).to_csv(TABLES / "gse292299_samespot.csv", index=False)
    pd.json_normalize(ring_rows).to_csv(TABLES / "gse292299_rings.csv", index=False)
    dump_json(TABLES / "gse292299_summary.json", summary)
    return summary


def run_gse318867():
    ex = DATA / "GSE318867" / "extracted"
    pairs = [
        ("VA", ex / "VA" / "VA_barcodes.tsv.gz", ex / "VA" / "VA_features.tsv.gz", ex / "VA" / "VA_matrix.mtx.gz", ex / "VA2" / "tissue_positions.csv"),
        ("VB", ex / "VB" / "VB_barcodes.tsv.gz", ex / "VB" / "VB_features.tsv.gz", ex / "VB" / "VB_matrix.mtx.gz", ex / "VB1" / "tissue_positions.csv"),
        ("VISA", ex / "VISA" / "VISA_barcodes.tsv.gz", ex / "VISA" / "VISA_features.tsv.gz", ex / "VISA" / "VISA_matrix.mtx.gz", ex / "VA1" / "tissue_positions.csv"),
    ]
    same_rows, ring_rows, skips = [], [], []
    for name, bc, feat, mtxp, posp in pairs:
        print("GSE318867", name, flush=True)
        barcodes, genes, mtx = load_mtx(bc, feat, mtxp)
        pos = read_positions_csv(posp)
        res = analyze_visium_section(name, barcodes, genes, mtx, pos)
        if "skip" in res:
            skips.append(res)
            continue
        same_rows.append(res["same"])
        ring_rows.extend(res["rings"])
    summary = summarize_visium(same_rows, ring_rows, "GSE318867")
    summary["skips"] = skips
    pd.json_normalize(same_rows).to_csv(TABLES / "gse318867_samespot.csv", index=False)
    pd.json_normalize(ring_rows).to_csv(TABLES / "gse318867_rings.csv", index=False)
    dump_json(TABLES / "gse318867_summary.json", summary)
    return summary


def load_cosmx(expr_path, meta_path):
    expr = pd.read_csv(expr_path)
    meta = pd.read_csv(meta_path)
    expr = expr[expr["cell_ID"] > 0].copy()
    key = expr["fov"].astype(str) + "_" + expr["cell_ID"].astype(str)
    meta_key = meta["fov"].astype(str) + "_" + meta["cell_ID"].astype(str)
    meta = meta.set_index(meta_key)
    expr = expr.set_index(key)
    both = expr.index.intersection(meta.index)
    expr = expr.loc[both]
    meta = meta.loc[both]
    gene_cols = [c for c in expr.columns if c not in {"fov", "cell_ID"} and not str(c).startswith("NegPrb")]
    counts = expr[gene_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    ncount = counts.sum(axis=1)
    keep = ncount >= 20
    counts = counts.loc[keep]
    meta = meta.loc[keep]
    ncount = ncount.loc[keep]
    genes = [str(c).upper() for c in counts.columns]
    counts.columns = genes
    return counts, meta, ncount


def run_gse299786():
    d = DATA / "GSE299786"
    tmas = [
        (
            "TMA1",
            d / "GSM9046088_Lung_Adenocarcinoma_TMA1_CosMx_exprMat_file.csv.gz",
            d / "GSM9046088_Lung_Adenocarcinoma_TMA1_CosMx_metadata_file.csv.gz",
        ),
        (
            "TMA2",
            d / "GSM9046089_Lung_Adenocarcinoma_TMA2_CosMx_exprMat_file.csv.gz",
            d / "GSM9046089_Lung_Adenocarcinoma_TMA2_CosMx_metadata_file.csv.gz",
        ),
    ]
    recs = []
    fov_rows = []
    knn_rows = []
    for name, ep, mp in tmas:
        print("GSE299786", name, flush=True)
        counts, meta, ncount = load_cosmx(ep, mp)
        genes = list(counts.columns)
        if "CLDN4" not in genes:
            recs.append({"tma": name, "skip": "CLDN4_absent", "n_genes": len(genes)})
            continue
        lib = np.asarray(ncount.to_numpy(), dtype=float).copy()
        lib[lib == 0] = 1.0
        logx = np.log1p(counts.to_numpy() * (1e4 / lib[:, None]))
        gidx = {g: i for i, g in enumerate(genes)}
        cl4 = logx[:, gidx["CLDN4"]]
        tac = logx[:, gidx["TACSTD2"]] if "TACSTD2" in gidx else np.full(len(cl4), np.nan)
        tsc = score_from_log(logx, gidx, present(T_GENES, genes))
        bsc = score_from_log(logx, gidx, present(B_GENES, genes))
        broad = score_from_log(logx, gidx, present(BROAD_EPI, genes))
        tb = np.nanmean(np.vstack([tsc, bsc]), axis=0)
        rec = {
            "tma": name,
            "n_cells": int(len(counts)),
            "n_fov": int(meta["fov"].nunique()),
            "n_genes": len(genes),
            "CLDN4_pct_pos": float((counts["CLDN4"] > 0).mean()),
            "same_CLDN4_TB": spearman(cl4, tb),
            "same_CLDN4_T": spearman(cl4, tsc),
            "same_CLDN4_B": spearman(cl4, bsc),
            "same_TACSTD2_TB": spearman(tac, tb),
            "partial_same_CLDN4_TB_ctrl_broad": partial_spearman(cl4, tb, broad),
        }
        # FOV means
        tmp = pd.DataFrame({"fov": meta["fov"].to_numpy(), "CLDN4": cl4, "TACSTD2": tac, "TB": tb, "T": tsc, "B": bsc, "broad": broad})
        fov = tmp.groupby("fov").mean(numeric_only=True)
        rec["fov_CLDN4_TB"] = spearman(fov["CLDN4"], fov["TB"])
        rec["fov_CLDN4_T"] = spearman(fov["CLDN4"], fov["T"])
        rec["fov_CLDN4_B"] = spearman(fov["CLDN4"], fov["B"])
        rec["fov_TACSTD2_TB"] = spearman(fov["TACSTD2"], fov["TB"])
        rec["fov_partial_CLDN4_TB_ctrl_broad"] = partial_spearman(fov["CLDN4"], fov["TB"], fov["broad"])
        fov.assign(tma=name).to_csv(TABLES / f"gse299786_{name}_fov_means.csv")
        fov_rows.append(fov.assign(tma=name))
        # within-FOV kNN-20 in pixel space (no invented micron scale)
        xs = meta["CenterX_global_px"].to_numpy(dtype=float)
        ys = meta["CenterY_global_px"].to_numpy(dtype=float)
        fovs = meta["fov"].to_numpy()
        nei = np.full(len(cl4), np.nan)
        for fv in np.unique(fovs):
            ii = np.where(fovs == fv)[0]
            if ii.size < 12:
                continue
            xy = np.column_stack([xs[ii], ys[ii]])
            d2 = ((xy[:, None, :] - xy[None, :, :]) ** 2).sum(axis=2)
            np.fill_diagonal(d2, np.inf)
            k = min(20, ii.size - 1)
            nn = np.argpartition(d2, kth=k - 1, axis=1)[:, :k]
            nei[ii] = tb[ii][nn].mean(axis=1)
        ok = np.isfinite(nei)
        rec["knn20_CLDN4_vs_neiTB"] = spearman(cl4[ok], nei[ok])
        rec["knn20_partial_CLDN4_vs_neiTB_ctrl_broad"] = partial_spearman(cl4[ok], nei[ok], broad[ok])
        knn_rows.append({"tma": name, **rec["knn20_CLDN4_vs_neiTB"]})
        recs.append(rec)
    dump_json(TABLES / "gse299786_summary.json", recs)
    return recs


def run_gse329813_gene_check():
    gx = pd.read_csv(DATA / "GSE329813_processed_data_file_normalized_data.csv.gz", usecols=[0])
    genes = set(gx.iloc[:, 0].astype(str).str.upper())
    rec = {
        "accession": "GSE329813",
        "platform": "GeoMx WTA-ish processed",
        "n_genes": len(genes),
        "CLDN4": "CLDN4" in genes,
        "TACSTD2": "TACSTD2" in genes,
        "T_genes": present(T_GENES, genes),
        "B_genes": present(B_GENES, genes),
        "decision": "skip_no_CLDN4",
        "note": "Neoadjuvant NSCLC GeoMx (22 patients, MPR/NMPR). Public normalized matrix has TACSTD2 and T/B genes but not CLDN4.",
    }
    dump_json(TABLES / "gse329813_gene_check.json", rec)
    return rec


def main():
    hunt = pd.read_csv(Path(__file__).with_name("hunt_catalog.tsv"), sep="\t")
    hunt.to_csv(TABLES / "hunt_catalog.tsv", sep="\t", index=False)
    out = {
        "GSE329813": run_gse329813_gene_check(),
        "GSE299786": run_gse299786(),
        "GSE318867": run_gse318867(),
        "GSE292299": run_gse292299(),
    }
    dump_json(TABLES / "summary.json", out)
    print("wrote", TABLES / "summary.json")


if __name__ == "__main__":
    main()
