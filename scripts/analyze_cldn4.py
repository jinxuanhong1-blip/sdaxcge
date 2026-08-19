#!/usr/bin/env python3
"""CLDN4-only Visium spot analysis: Spearman, Q4/Q1 nearest-CD8 distance, KRT8 residual."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import sys
from dataclasses import dataclass
from typing import Iterable

import numpy as np
from scipy import sparse, stats
from scipy.spatial import cKDTree

GENES = ("CLDN4", "CD8A", "KRT8")
MIN_UMI = 100


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def read_tsv_col(path: str, col: int = 0) -> list[str]:
    out = []
    with _open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            out.append(parts[col])
    return out


def read_features(path: str) -> list[str]:
    names = []
    with _open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                names.append(parts[1])
            else:
                names.append(parts[0])
    return names


def read_mtx(path: str) -> sparse.csr_matrix:
    with _open(path) as f:
        header = f.readline()
        if not header.startswith("%%MatrixMarket"):
            raise ValueError(f"not MTX: {path}")
        line = f.readline()
        while line.startswith("%"):
            line = f.readline()
        n_row, n_col, n_nz = map(int, line.split())
        rows = np.empty(n_nz, dtype=np.int32)
        cols = np.empty(n_nz, dtype=np.int32)
        data = np.empty(n_nz, dtype=np.float64)
        i = 0
        for line in f:
            if not line.strip():
                continue
            r, c, v = line.split()
            rows[i] = int(r) - 1
            cols[i] = int(c) - 1
            data[i] = float(v)
            i += 1
        if i != n_nz:
            rows, cols, data = rows[:i], cols[:i], data[:i]
    return sparse.csr_matrix((data, (rows, cols)), shape=(n_row, n_col))


def read_h5(path: str):
    import h5py

    with h5py.File(path, "r") as h:
        root = "matrix" if "matrix" in h else ""
        grp = h[root] if root else h
        data = grp["data"][:]
        indices = grp["indices"][:]
        indptr = grp["indptr"][:]
        shape = tuple(int(x) for x in grp["shape"][:])
        mat = sparse.csc_matrix((data, indices, indptr), shape=shape).tocsr()
        feats = grp["features"]
        if "name" in feats:
            names = [x.decode() if isinstance(x, bytes) else str(x) for x in feats["name"][:]]
        elif "id" in feats:
            names = [x.decode() if isinstance(x, bytes) else str(x) for x in feats["id"][:]]
        else:
            raise KeyError("no feature names in h5")
        barcodes = [x.decode() if isinstance(x, bytes) else str(x) for x in grp["barcodes"][:]]
    return names, barcodes, mat


def read_positions(path: str) -> dict[str, tuple[int, float, float]]:
    """barcode -> (in_tissue, pxl_row, pxl_col)."""
    out: dict[str, tuple[int, float, float]] = {}
    with _open(path) as f:
        first = f.readline().strip()
        has_header = "barcode" in first.lower() or "in_tissue" in first.lower()
        rows = []
        if has_header:
            reader = csv.DictReader([first] + f.readlines())
            for rec in reader:
                bc = rec.get("barcode") or rec.get("barcodes")
                if bc is None:
                    # first unnamed column
                    bc = list(rec.values())[0]
                in_t = int(float(rec.get("in_tissue", 1)))
                # 10x: pxl_row_in_fullres, pxl_col_in_fullres
                pr = rec.get("pxl_row_in_fullres") or rec.get("pxl_row") or rec.get("imagerow")
                pc = rec.get("pxl_col_in_fullres") or rec.get("pxl_col") or rec.get("imagecol")
                if pr is None or pc is None:
                    vals = list(rec.values())
                    # barcode, in_tissue, array_row, array_col, pxl_row, pxl_col
                    pr, pc = vals[4], vals[5]
                out[bc] = (in_t, float(pr), float(pc))
        else:
            # barcode,in_tissue,array_row,array_col,pxl_row,pxl_col
            f.seek(0)
            for line in f:
                parts = line.strip().split(",")
                if len(parts) < 6:
                    continue
                out[parts[0]] = (int(float(parts[1])), float(parts[4]), float(parts[5]))
    return out


def find_gene(names: list[str], symbol: str) -> int | None:
    up = [n.upper() for n in names]
    if symbol in up:
        return up.index(symbol)
    # sometimes "CLDN4.1" or "GRCh38_CLDN4"
    for i, n in enumerate(up):
        if n == symbol or n.endswith("_" + symbol) or n.split(".")[0] == symbol:
            return i
    return None


def log_norm(counts: np.ndarray) -> np.ndarray:
    lib = counts.sum(axis=0) if counts.ndim == 2 else None
    # counts is (3, n) gene x spot
    lib = counts.sum(axis=0)
    lib = np.where(lib > 0, lib, 1.0)
    return np.log1p(1e4 * counts / lib)


def ols_residual(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    x1 = np.column_stack([np.ones(len(x)), x])
    beta, *_ = np.linalg.lstsq(x1, y, rcond=None)
    return y - x1 @ beta


def mannwhitney(a: np.ndarray, b: np.ndarray):
    if len(a) < 3 or len(b) < 3:
        return np.nan, np.nan
    res = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(res.statistic), float(res.pvalue)


@dataclass
class SlideResult:
    series: str
    sample: str
    n_spots: int
    n_cd8pos: int
    genes_present: str
    spearman_r: float
    spearman_p: float
    krt8_resid_r: float
    krt8_resid_p: float
    q1_n: int
    q4_n: int
    q1_med_dist: float
    q4_med_dist: float
    q4_minus_q1: float
    dist_u: float
    dist_p: float
    note: str


def analyze_slide(
    series: str,
    sample: str,
    names: list[str],
    barcodes: list[str],
    mat: sparse.csr_matrix,
    positions: dict[str, tuple[int, float, float]] | None,
) -> SlideResult:
    idx = {g: find_gene(names, g) for g in GENES}
    present = [g for g, i in idx.items() if i is not None]
    missing = [g for g, i in idx.items() if i is None]
    note = ""
    if missing:
        note = "missing:" + ",".join(missing)
        return SlideResult(
            series, sample, 0, 0, ",".join(present),
            np.nan, np.nan, np.nan, np.nan, 0, 0,
            np.nan, np.nan, np.nan, np.nan, np.nan, note,
        )

    # restrict to in-tissue if positions available
    keep = np.ones(len(barcodes), dtype=bool)
    xy = np.full((len(barcodes), 2), np.nan)
    if positions:
        for i, bc in enumerate(barcodes):
            rec = positions.get(bc)
            if rec is None:
                # try without extra suffix
                rec = positions.get(bc.split("-")[0])
            if rec is None:
                keep[i] = False
                continue
            in_t, pr, pc = rec
            keep[i] = in_t == 1
            xy[i] = (pc, pr)
    umi = np.asarray(mat.sum(axis=0)).ravel()
    keep &= umi >= MIN_UMI
    if keep.sum() < 50:
        return SlideResult(
            series, sample, int(keep.sum()), 0, ",".join(present),
            np.nan, np.nan, np.nan, np.nan, 0, 0,
            np.nan, np.nan, np.nan, np.nan, np.nan, "too_few_spots",
        )

    sub = mat[:, keep]
    xy = xy[keep]
    raw = np.vstack([
        np.asarray(sub[idx["CLDN4"], :].todense()).ravel(),
        np.asarray(sub[idx["CD8A"], :].todense()).ravel(),
        np.asarray(sub[idx["KRT8"], :].todense()).ravel(),
    ])
    norm = log_norm(raw)
    cldn4, cd8a, krt8 = norm

    r, p = stats.spearmanr(cldn4, cd8a)
    resid = ols_residual(cldn4, krt8)
    rr, rp = stats.spearmanr(resid, cd8a)

    cd8_pos = raw[1] > 0
    n_cd8 = int(cd8_pos.sum())
    q1_n = q4_n = 0
    q1_med = q4_med = delta = u = dp = np.nan
    if positions is None or np.isnan(xy).all():
        note = (note + ";no_coords").strip(";")
    elif n_cd8 < 5:
        note = (note + ";too_few_cd8").strip(";")
    else:
        tree = cKDTree(xy[cd8_pos])
        dist, _ = tree.query(xy, k=1)
        q1_cut, q3_cut = np.quantile(cldn4, [0.25, 0.75])
        q1 = cldn4 <= q1_cut
        q4 = cldn4 >= q3_cut
        d1 = dist[q1]
        d4 = dist[q4]
        q1_n, q4_n = int(q1.sum()), int(q4.sum())
        q1_med = float(np.median(d1))
        q4_med = float(np.median(d4))
        delta = q4_med - q1_med
        u, dp = mannwhitney(d4, d1)

    return SlideResult(
        series, sample, int(keep.sum()), n_cd8, ",".join(present),
        float(r), float(p), float(rr), float(rp),
        q1_n, q4_n, q1_med, q4_med, float(delta), float(u), float(dp), note,
    )


def discover_slides(series_dir: str) -> list[dict]:
    """Pair matrix/h5 files with barcodes/features/positions in a series folder."""
    files = sorted(os.listdir(series_dir))
    slides = []
    # MTX triples
    mtx = [f for f in files if f.endswith("matrix.mtx.gz") or f.endswith("matrix.mtx")]
    for m in mtx:
        prefix = m.replace("_matrix.mtx.gz", "").replace("_matrix.mtx", "")
        bc = next((f for f in files if f.startswith(prefix) and "barcodes" in f), None)
        ft = next((f for f in files if f.startswith(prefix) and "features" in f), None)
        pos = next((f for f in files if f.startswith(prefix) and "tissue_positions" in f), None)
        if bc and ft:
            slides.append({
                "sample": prefix,
                "kind": "mtx",
                "matrix": os.path.join(series_dir, m),
                "barcodes": os.path.join(series_dir, bc),
                "features": os.path.join(series_dir, ft),
                "positions": os.path.join(series_dir, pos) if pos else None,
            })
    h5_by_prefix: dict[str, str] = {}
    for f in files:
        if "filtered_feature_bc_matrix.h5" not in f:
            continue
        prefix = f.split("_filtered_feature_bc_matrix")[0]
        prev = h5_by_prefix.get(prefix)
        # Prefer a decompressed .h5 over .h5.gz when both exist.
        if prev is None or (f.endswith(".h5") and not f.endswith(".h5.gz")):
            h5_by_prefix[prefix] = f
    for prefix, h in h5_by_prefix.items():
        pos = next((f for f in files if f.startswith(prefix) and "tissue_positions" in f), None)
        slides.append({
            "sample": prefix,
            "kind": "h5",
            "matrix": os.path.join(series_dir, h),
            "barcodes": None,
            "features": None,
            "positions": os.path.join(series_dir, pos) if pos else None,
        })
    return slides


def load_slide(slide: dict):
    if slide["kind"] == "h5":
        path = slide["matrix"]
        # h5py can read gzipped? usually not — decompress if needed
        if path.endswith(".gz"):
            import gzip as gz
            import tempfile
            raw = path[:-3]
            if not os.path.exists(raw):
                with gz.open(path, "rb") as src, open(raw, "wb") as dst:
                    dst.write(src.read())
            names, barcodes, mat = read_h5(raw)
        else:
            names, barcodes, mat = read_h5(path)
    else:
        names = read_features(slide["features"])
        barcodes = read_tsv_col(slide["barcodes"], 0)
        mat = read_mtx(slide["matrix"])
    pos = read_positions(slide["positions"]) if slide["positions"] else None
    return names, barcodes, mat, pos


def write_csv(path: str, rows: list[SlideResult]) -> None:
    fields = [
        "series", "sample", "n_spots", "n_cd8pos", "genes_present",
        "spearman_r", "spearman_p", "krt8_resid_r", "krt8_resid_p",
        "q1_n", "q4_n", "q1_med_dist", "q4_med_dist", "q4_minus_q1",
        "dist_u", "dist_p", "note",
    ]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r.__dict__)


def summarize(rows: list[SlideResult]) -> dict:
    ok = [r for r in rows if not np.isnan(r.spearman_r)]
    dist_ok = [r for r in ok if not np.isnan(r.q4_minus_q1)]
    if not ok:
        return {"n_slides": 0}
    rs = np.array([r.spearman_r for r in ok])
    rrs = np.array([r.krt8_resid_r for r in ok])
    out = {
        "n_slides": len(ok),
        "n_slides_total": len(rows),
        "median_spearman_r": float(np.median(rs)),
        "mean_spearman_r": float(np.mean(rs)),
        "n_neg_spearman": int((rs < 0).sum()),
        "median_krt8_resid_r": float(np.median(rrs)),
        "mean_krt8_resid_r": float(np.mean(rrs)),
        "n_neg_krt8_resid": int((rrs < 0).sum()),
        "n_dist_slides": len(dist_ok),
    }
    if dist_ok:
        dlt = np.array([r.q4_minus_q1 for r in dist_ok])
        out["median_q4_minus_q1"] = float(np.median(dlt))
        out["mean_q4_minus_q1"] = float(np.mean(dlt))
        out["n_q4_farther"] = int((dlt > 0).sum())
        # Stouffer combine of distance p-values with sign of delta
        zs = []
        for r in dist_ok:
            if r.dist_p > 0 and not np.isnan(r.dist_p):
                z = stats.norm.isf(r.dist_p / 2.0)
                zs.append(np.sign(r.q4_minus_q1) * z)
        if zs:
            zc = float(np.sum(zs) / np.sqrt(len(zs)))
            out["stouffer_z_q4_minus_q1"] = zc
            out["stouffer_p_q4_minus_q1"] = float(2 * stats.norm.sf(abs(zc)))
    return out


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("series_dirs", nargs="+")
    p.add_argument("--out", default="results/cldn4_visium_leftover.csv")
    args = p.parse_args()
    rows: list[SlideResult] = []
    for d in args.series_dirs:
        series = os.path.basename(os.path.abspath(d))
        slides = discover_slides(d)
        print(f"=== {series}: {len(slides)} slides ===")
        for sl in slides:
            print(f"  {sl['sample']} ...", flush=True)
            try:
                names, barcodes, mat, pos = load_slide(sl)
                res = analyze_slide(series, sl["sample"], names, barcodes, mat, pos)
            except Exception as e:
                res = SlideResult(
                    series, sl["sample"], 0, 0, "",
                    np.nan, np.nan, np.nan, np.nan, 0, 0,
                    np.nan, np.nan, np.nan, np.nan, np.nan, f"error:{e}",
                )
            rows.append(res)
            print(
                f"    n={res.n_spots} cd8+={res.n_cd8pos} r={res.spearman_r:.4f} "
                f"resid_r={res.krt8_resid_r:.4f} dQ4-Q1={res.q4_minus_q1:.1f} "
                f"p_dist={res.dist_p:.3g} {res.note}",
                flush=True,
            )
    write_csv(args.out, rows)
    summary_path = os.path.splitext(args.out)[0] + "_summary.json"
    by = {}
    for series in sorted({r.series for r in rows}):
        by[series] = summarize([r for r in rows if r.series == series])
    by["ALL"] = summarize(rows)
    with open(summary_path, "w") as f:
        json.dump(by, f, indent=2)
    print("wrote", args.out)
    print("wrote", summary_path)
    print(json.dumps(by, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
