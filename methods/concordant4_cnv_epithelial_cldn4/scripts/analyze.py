#!/usr/bin/env python3
"""Concordant-4: CopyKAT + InferCNV malignant calls, then CLDN4 vs T/NK.

Annotation-only malignant is recomputed on the same units as the locked
concordant-4 (GSE123902, GSE131907, GSE205335, GSE189357; n = 65).

Both callers build a 5 Mb autosomal profile (chr6 excluded) of each
epithelial cell relative to that unit's own T/NK cells, then subtract the
median profile of normal epithelial cells from the same cohort. GSE189357
has no normal epithelium on GEO and borrows the GSE123902 normal-lung null.
The CNV+ threshold is the normal-epithelium 95th percentile of residual
mean squared error, fixed before the CLDN4 test.

CopyKAT features (Gao et al. 2021): log(sqrt(x)+sqrt(x+1)), cell-mean
center, local-level Kalman smooth (dV=0.16, dW=0.001), T/NK median
baseline. Ward clustering on the normal-corrected 5 Mb residual. The
cluster enriched for the lower-residual half is diploid. Other clusters
are labeled by 1D Wasserstein distance. Consensus correlation >= 0.6
calls every cell diploid (package v1.2). A cell is CopyKAT-aneuploid only
if that label is aneuploid and its residual exceeds the normal-epithelium
threshold.

InferCNV features (Trinity preliminary mode, HMM off): log2(CP10k+1),
subtract the T/NK mean, within-chromosome moving mean (window 101), cell
center. InferCNV+ = residual MSE above the normal-epithelium 95th percentile.

Consensus = CopyKAT-aneuploid AND InferCNV+. CLDN4 %pos uses only those
cells. T/NK fraction is the full-unit fraction from the locked test.

Not run: copykat's MCMC breakpoint sampler, inferCNV's HMM, and Numbat.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import subprocess
import time
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread
from sklearn.cluster import AgglomerativeClustering

ROOT = Path(__file__).resolve().parents[1]
EPI_CAP = 800
REF_CAP = 250
MIN_POS = 20
MIN_REF = 30
MIN_GENES = 800
COPYKAT_WIN = 25
INFER_WIN = 101
MIN_BREADTH = 4
COR_ALL_DIPLOID = 0.6
COR_LOW_CONF = 0.4
DR_MIN = 0.05

LOCKED_RHO = {
    "GSE123902": -0.659,
    "GSE131907": -0.522,
    "GSE205335": -0.435,
    "GSE189357": -0.600,
}
COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
COLORS = {
    "GSE123902": "#1b4f72",
    "GSE131907": "#b9770e",
    "GSE205335": "#196f3d",
    "GSE189357": "#6c3483",
}


def say(*parts) -> None:
    print(time.strftime("%H:%M:%S"), *parts, flush=True)


def read_tsv(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            out = {}
            for col in columns:
                val = row.get(col, "")
                if isinstance(val, float):
                    out[col] = "" if not math.isfinite(val) else f"{val:.6g}"
                else:
                    out[col] = val
            writer.writerow(out)


def cldn4_stats(counts: np.ndarray) -> tuple[float, float]:
    if len(counts) == 0:
        return float("nan"), float("nan")
    x = np.asarray(counts, dtype=float)
    return 100.0 * float(np.mean(x > 0)), float(np.mean(np.log1p(x)))


def pct_on(cldn4: np.ndarray, mask: np.ndarray) -> tuple[float, int]:
    n = int(np.sum(mask))
    if n < MIN_POS:
        return float("nan"), n
    return 100.0 * float(np.mean(cldn4[mask] > 0)), n


# ---------------------------------------------------------------------------
# Gene order (hg38 refGene, longest transcript per symbol, autosomes)
# ---------------------------------------------------------------------------

def load_gene_order(refgene: Path, cache: Path) -> dict[str, tuple[int, int]]:
    if cache.exists():
        order = {}
        with cache.open() as handle:
            next(handle)
            for line in handle:
                sym, chrom, start = line.rstrip("\n").split("\t")
                order[sym] = (int(chrom), int(start))
        say("gene order", len(order), "from", cache.name)
        return order
    best: dict[str, tuple[int, int, int]] = {}
    with gzip.open(refgene, "rt") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[2]
            if chrom.startswith("chr"):
                chrom = chrom[3:]
            if chrom not in {str(i) for i in range(1, 23)}:
                continue
            start = int(parts[4])
            end = int(parts[5])
            sym = parts[12].upper()
            if not sym or sym == "NA":
                continue
            span = end - start
            prev = best.get(sym)
            if prev is None or span > prev[2]:
                best[sym] = (int(chrom), start, span)
    cache.parent.mkdir(parents=True, exist_ok=True)
    symbols = sorted(best, key=lambda s: (best[s][0], best[s][1], s))
    with cache.open("w") as handle:
        handle.write("symbol\tchrom\tstart\n")
        for sym in symbols:
            chrom, start, _span = best[sym]
            handle.write(f"{sym}\t{chrom}\t{start}\n")
    order = {sym: (best[sym][0], best[sym][1]) for sym in symbols}
    say("gene order", len(order), "wrote", cache)
    return order


def order_indices(genes: list[str], order_map: dict[str, tuple[int, int]]):
    hits = []
    seen = set()
    for i, gene in enumerate(genes):
        gene = gene.upper()
        if gene in seen or gene not in order_map:
            continue
        seen.add(gene)
        chrom, start = order_map[gene]
        hits.append((chrom, start, gene, i))
    hits.sort()
    names = [h[2] for h in hits]
    chroms = [h[0] for h in hits]
    idx = [h[3] for h in hits]
    return names, chroms, idx


def rng_for(dataset: str, unit_id: str) -> np.random.Generator:
    import zlib

    return np.random.default_rng(zlib.crc32(f"{dataset}|{unit_id}".encode()) & 0xFFFFFFFF)


def draw(idx: np.ndarray, cap: int, rng: np.random.Generator) -> np.ndarray:
    idx = np.asarray(idx)
    if len(idx) <= cap:
        return np.sort(idx)
    return np.sort(rng.choice(idx, size=cap, replace=False))


def first_occurrence(genes: list[str]) -> tuple[list[str], np.ndarray]:
    seen = set()
    keep = []
    names = []
    for i, gene in enumerate(genes):
        gene = gene.upper()
        if gene in seen:
            continue
        seen.add(gene)
        keep.append(i)
        names.append(gene)
    return names, np.asarray(keep, dtype=int)


def marker_masks(genes: list[str], counts_gxc: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """counts_gxc is genes x cells, already unique-uppercased."""
    index = {g: i for i, g in enumerate(genes)}

    def gt0(name: str) -> np.ndarray:
        i = index.get(name)
        if i is None:
            return np.zeros(counts_gxc.shape[1], dtype=bool)
        return counts_gxc[i] > 0

    mal = (gt0("EPCAM") | gt0("KRT8") | gt0("KRT18") | gt0("KRT19")) & ~gt0("PTPRC")
    tnk = (gt0("CD3D") | gt0("CD3E") | gt0("CD8A") | gt0("NKG7") | gt0("GNLY") | gt0("KLRD1")) & ~mal
    return mal, tnk


def save_prepared(cache: Path, meta: dict, payload: dict | None) -> None:
    dest = cache / meta["dataset"]
    dest.mkdir(parents=True, exist_ok=True)
    uid = meta["unit_id"]
    (dest / f"{uid}.json").write_text(json.dumps(meta))
    if payload is None:
        np.savez(dest / f"{uid}.npz", counts=np.zeros((0, 0), np.float32))
        return
    np.savez(
        dest / f"{uid}.npz",
        counts=payload["counts"],
        chrom=np.asarray(payload["chrom"], dtype=np.int16),
        genes=np.asarray(payload["genes"]),
        libsize=np.asarray(payload["libsize"], dtype=np.float32),
        is_epi=np.asarray(payload["is_epi"], dtype=np.uint8),
        is_ref=np.asarray(payload["is_ref"], dtype=np.uint8),
        is_annot=np.asarray(payload["is_annot"], dtype=np.uint8),
        cldn4=np.asarray(payload["cldn4"], dtype=np.float32),
    )


def prepared_ok(cache: Path, dataset: str, unit_id: str) -> bool:
    return (cache / dataset / f"{unit_id}.json").exists() and (cache / dataset / f"{unit_id}.npz").exists()


def pack_unit(
    genes: list[str],
    counts_gxc: np.ndarray,
    libsize_all: np.ndarray,
    epi_idx: np.ndarray,
    ref_idx: np.ndarray,
    annot_idx_in_matrix: np.ndarray,
    order_map: dict[str, tuple[int, int]],
    cldn4_all_cells: np.ndarray | None = None,
) -> dict | None:
    """Build a CNV matrix on the drawn epithelial + T/NK cells."""
    names, chroms, row_idx = order_indices(genes, order_map)
    if len(row_idx) < MIN_GENES:
        say("  too few ordered genes", len(row_idx))
        return None
    cells = np.concatenate([epi_idx, ref_idx])
    mat = np.asarray(counts_gxc[row_idx][:, cells], dtype=np.float32)
    lib = np.asarray(libsize_all[cells], dtype=np.float32)
    n_epi = len(epi_idx)
    is_epi = np.zeros(len(cells), dtype=np.uint8)
    is_ref = np.zeros(len(cells), dtype=np.uint8)
    is_epi[:n_epi] = 1
    is_ref[n_epi:] = 1
    annot_set = set(int(i) for i in annot_idx_in_matrix)
    is_annot = np.array([1 if int(c) in annot_set else 0 for c in cells], dtype=np.uint8)
    if cldn4_all_cells is None:
        if "CLDN4" in genes:
            cldn4 = np.asarray(counts_gxc[genes.index("CLDN4"), cells], dtype=np.float32)
        else:
            cldn4 = np.zeros(len(cells), dtype=np.float32)
    else:
        cldn4 = np.asarray(cldn4_all_cells[cells], dtype=np.float32)
    # Drop genes that are zero in this draw.
    keep = mat.sum(axis=1) > 0
    kept_i = np.where(keep)[0]
    return {
        "counts": mat[keep],
        "chrom": [chroms[i] for i in kept_i],
        "genes": [names[i] for i in kept_i],
        "libsize": lib,
        "is_epi": is_epi,
        "is_ref": is_ref,
        "is_annot": is_annot,
        "cldn4": cldn4,
    }


# ---------------------------------------------------------------------------
# CopyKAT label + InferCNV score
# ---------------------------------------------------------------------------

def moving_mean(block: np.ndarray, window: int) -> np.ndarray:
    n, k = block.shape
    if n < 5:
        return block.copy()
    if window % 2 == 0:
        window -= 1
    window = min(window, n if n % 2 else n - 1)
    if window < 3:
        return block.copy()
    pad = window // 2
    padded = np.pad(block, ((pad, pad), (0, 0)), mode="edge")
    csum = np.cumsum(padded, axis=0)
    csum = np.vstack([np.zeros((1, k), dtype=csum.dtype), csum])
    return (csum[window : window + n] - csum[:n]) / window


def smooth_by_chrom(values: np.ndarray, chrom: np.ndarray, window: int) -> np.ndarray:
    out = np.empty_like(values)
    for c in np.unique(chrom):
        ix = np.where(chrom == c)[0]
        out[ix] = moving_mean(values[ix], window)
    return out


def kalman_smooth(values: np.ndarray, d_v: float = 0.16, d_w: float = 0.001) -> np.ndarray:
    """Local-level RTS smoother, copykat's dlmModPoly(order=1, dV, dW)."""
    t_steps, _n = values.shape
    mean = np.empty_like(values, dtype=np.float64)
    var = np.empty(t_steps, dtype=np.float64)
    pred = np.empty_like(mean)
    pred_var = np.empty(t_steps, dtype=np.float64)
    mean[0] = values[0]
    var[0] = 1.0
    for t in range(1, t_steps):
        pred[t] = mean[t - 1]
        pred_var[t] = var[t - 1] + d_w
        gain = pred_var[t] / (pred_var[t] + d_v)
        mean[t] = pred[t] + gain * (values[t] - pred[t])
        var[t] = (1.0 - gain) * pred_var[t]
    smooth = np.empty_like(mean)
    smooth[-1] = mean[-1]
    for t in range(t_steps - 2, -1, -1):
        c = var[t] / pred_var[t + 1]
        smooth[t] = mean[t] + c * (smooth[t + 1] - pred[t + 1])
    return smooth


def window_means(values: np.ndarray, window: int) -> np.ndarray | None:
    n_win = values.shape[0] // window
    if n_win < 20:
        return None
    trimmed = values[: n_win * window]
    return trimmed.reshape(n_win, window, values.shape[1]).mean(axis=1)


def copykat_aneuploid(raw: np.ndarray, chrom: np.ndarray, is_ref: np.ndarray) -> tuple[np.ndarray, dict]:
    """Return a boolean per cell (True = aneuploid) and a status dict."""
    n_cells = raw.shape[1]
    diploid = np.zeros(n_cells, dtype=bool)
    info = {"copykat_status": "too_few_cells", "copykat_cor": float("nan")}
    if int(is_ref.sum()) < MIN_REF or n_cells < MIN_REF + MIN_POS:
        return diploid, info
    transformed = np.log(np.sqrt(raw) + np.sqrt(raw + 1.0))
    transformed -= transformed.mean(axis=0, keepdims=True)
    smoothed = kalman_smooth(transformed.astype(np.float64))
    baseline = np.median(smoothed[:, is_ref], axis=1, keepdims=True)
    relative = smoothed - baseline
    relative -= relative.mean(axis=0, keepdims=True)
    windows = window_means(relative, COPYKAT_WIN)
    if windows is None:
        info["copykat_status"] = "too_few_windows"
        return diploid, info
    features = np.asarray(windows.T, dtype=np.float64)
    labels = None
    for k in (4, 3, 2):
        if n_cells <= k:
            continue
        trial = AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(features)
        counts = np.bincount(trial)
        if counts.min() > 10:
            labels = trial
            break
    if labels is None:
        info["copykat_status"] = "cluster_failed"
        return diploid, info
    ref_total = float(is_ref.sum())
    fractions = []
    for cluster in range(labels.max() + 1):
        members = labels == cluster
        fractions.append(float(is_ref[members].sum()) / ref_total)
    fractions = np.asarray(fractions)
    diploid_id = int(np.argmax(fractions))
    aneuploid_id = int(np.argmin(fractions))
    dip_profile = np.median(windows[:, labels == diploid_id], axis=1)
    ane_profile = np.median(windows[:, labels == aneuploid_id], axis=1)
    if np.std(dip_profile) == 0 or np.std(ane_profile) == 0:
        corr = 1.0
    else:
        corr = float(np.corrcoef(dip_profile, ane_profile)[0, 1])
    info["copykat_cor"] = corr
    if not math.isfinite(corr) or corr >= COR_ALL_DIPLOID:
        info["copykat_status"] = "all_diploid_cor>=0.6"
        return diploid, info
    aneuploid = np.zeros(n_cells, dtype=bool)
    for cluster in range(labels.max() + 1):
        profile = np.median(windows[:, labels == cluster], axis=1)
        dist_dip = float(stats.wasserstein_distance(profile, dip_profile))
        dist_ane = float(stats.wasserstein_distance(profile, ane_profile))
        if dist_ane <= dist_dip and cluster != diploid_id:
            aneuploid[labels == cluster] = True
        elif dist_ane < dist_dip:
            aneuploid[labels == cluster] = True
    # Known T/NK cells are the diploid anchor; do not call them malignant.
    aneuploid[is_ref] = False
    if COR_LOW_CONF <= corr < COR_ALL_DIPLOID:
        info["copykat_status"] = "low_conf"
    else:
        info["copykat_status"] = "separated"
    return aneuploid, info


def infercnv_calls(raw: np.ndarray, chrom: np.ndarray, libsize: np.ndarray, is_ref: np.ndarray) -> dict:
    scale = 1e4 / np.maximum(libsize.astype(np.float64), 1.0)
    lognorm = np.log2(raw.astype(np.float64) * scale + 1.0)
    ref = lognorm[:, is_ref]
    resid = lognorm - ref.mean(axis=1, keepdims=True)
    smooth = smooth_by_chrom(resid, chrom, INFER_WIN)
    smooth = smooth - smooth.mean(axis=0, keepdims=True)
    ref_s = smooth[:, is_ref]
    sd = np.maximum(ref_s.std(axis=1, keepdims=True), 1e-8)
    denoised = np.where(np.abs(smooth) <= sd, 0.0, smooth)
    use = chrom != 6
    if int(use.sum()) < 100:
        use = np.ones(len(chrom), dtype=bool)
    score = np.mean(denoised[use] ** 2, axis=0)
    thr = float(np.quantile(score[is_ref], 0.95))
    p95 = score > thr
    n_alt = np.zeros(raw.shape[1], dtype=int)
    for c in np.unique(chrom):
        if int(c) == 6:
            continue
        ix = np.where(chrom == c)[0]
        if len(ix) < 20:
            continue
        stat = denoised[ix].mean(axis=0)
        cut = float(np.quantile(np.abs(stat[is_ref]), 0.95))
        n_alt += (np.abs(stat) > cut).astype(int)
    breadth = p95 & (n_alt >= MIN_BREADTH)
    return {"score": score, "p95": p95, "breadth": breadth, "n_alt": n_alt, "thr": thr}


BIN_BP = 5_000_000


def genome_bins(order_map: dict[str, tuple[int, int]]) -> list[tuple[int, int]]:
    counts = defaultdict(int)
    for _sym, (chrom, start) in order_map.items():
        if chrom == 6:
            continue
        counts[(chrom, start // BIN_BP)] += 1
    return sorted(b for b, n in counts.items() if n >= 8)


def _bin_matrix(values: np.ndarray, genes: list[str], chrom: np.ndarray, bin_keys: list[tuple[int, int]], order_map) -> np.ndarray:
    """Mean of gene values inside each 5 Mb bin. values is genes x cells."""
    n_bins = len(bin_keys)
    n_cells = values.shape[1]
    acc = np.zeros((n_bins, n_cells), dtype=np.float64)
    cnt = np.zeros(n_bins, dtype=np.int32)
    index = {key: i for i, key in enumerate(bin_keys)}
    for i, gene in enumerate(genes):
        if int(chrom[i]) == 6:
            continue
        loc = order_map.get(gene)
        if loc is None:
            continue
        key = (loc[0], loc[1] // BIN_BP)
        b = index.get(key)
        if b is None:
            continue
        acc[b] += values[i]
        cnt[b] += 1
    out = np.full((n_cells, n_bins), np.nan, dtype=np.float64)
    ok = cnt >= 3
    out[:, ok] = (acc[ok] / cnt[ok, None]).T
    return out


def epithelial_profiles(path: Path, order_map, bin_keys) -> dict | None:
    with np.load(path, allow_pickle=False) as data:
        if "genes" not in data.files or "is_epi" not in data.files or data["counts"].size == 0:
            return None
        counts = data["counts"]
        chrom = data["chrom"].astype(int)
        genes = [str(g) for g in data["genes"].tolist()]
        libsize = data["libsize"].astype(np.float64)
        is_epi = data["is_epi"].astype(bool)
        is_ref = data["is_ref"].astype(bool)
        is_annot = data["is_annot"].astype(bool)
        cldn4 = data["cldn4"].astype(np.float64)
    if int(is_ref.sum()) < MIN_REF or int(is_epi.sum()) < 10:
        return None
    dr = (counts > 0).mean(axis=1)
    keep = (dr >= DR_MIN) & (chrom != 6)
    keep_chr = np.zeros(len(chrom), dtype=bool)
    for c in np.unique(chrom):
        if int(c) == 6:
            continue
        mask = (chrom == c) & keep
        if int(mask.sum()) >= 5:
            keep_chr[mask] = True
    if int(keep_chr.sum()) < MIN_GENES:
        return None
    raw = counts[keep_chr].astype(np.float64)
    chrom_k = chrom[keep_chr]
    genes_k = [genes[i] for i in np.where(keep_chr)[0]]
    # InferCNV preliminary: log2 CP10k, T/NK mean, within-chromosome smooth.
    scale = 1e4 / np.maximum(libsize, 1.0)
    lognorm = np.log2(raw * scale + 1.0)
    resid = lognorm - lognorm[:, is_ref].mean(axis=1, keepdims=True)
    smooth = smooth_by_chrom(resid, chrom_k, INFER_WIN)
    smooth -= smooth.mean(axis=0, keepdims=True)
    infer_bins = _bin_matrix(smooth, genes_k, chrom_k, bin_keys, order_map)
    # CopyKAT: package transform, Kalman smooth, T/NK median baseline, then the same bins.
    transformed = np.log(np.sqrt(raw) + np.sqrt(raw + 1.0))
    transformed -= transformed.mean(axis=0, keepdims=True)
    smoothed = kalman_smooth(transformed)
    baseline = np.median(smoothed[:, is_ref], axis=1, keepdims=True)
    relative = smoothed - baseline
    relative -= relative.mean(axis=0, keepdims=True)
    copy_bins = _bin_matrix(relative, genes_k, chrom_k, bin_keys, order_map)
    epi = np.where(is_epi)[0]
    return {
        "infer": infer_bins[epi],
        "copy": copy_bins[epi],
        "cldn4": cldn4[epi],
        "is_annot": is_annot[epi],
        "n_genes": int(keep_chr.sum()),
        "n_epi": int(len(epi)),
        "n_ref": int(is_ref.sum()),
    }


def _mse(profile: np.ndarray, null: np.ndarray) -> np.ndarray:
    delta = profile - null.reshape(1, -1)
    n_obs = np.sum(np.isfinite(delta), axis=1)
    mse = np.full(profile.shape[0], np.nan)
    ok = n_obs >= 30
    mse[ok] = np.nanmean(delta[ok] ** 2, axis=1)
    return mse


def _copykat_from_residual(resid: np.ndarray, mse: np.ndarray, thr: float) -> tuple[np.ndarray, dict]:
    n_cells = resid.shape[0]
    call = np.zeros(n_cells, dtype=bool)
    info = {"copykat_status": "too_few_cells", "copykat_cor": float("nan")}
    finite = np.isfinite(mse)
    if int(finite.sum()) < MIN_POS + 10:
        return call, info
    features = np.nan_to_num(resid[finite], nan=0.0)
    keep = np.abs(features).sum(axis=0) > 0
    features = features[:, keep]
    if features.shape[1] < 20:
        info["copykat_status"] = "too_few_windows"
        return call, info
    labels_full = None
    for k in (4, 3, 2):
        if features.shape[0] <= k:
            continue
        trial = AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(features)
        if np.bincount(trial).min() > 10:
            labels_full = trial
            break
    if labels_full is None:
        info["copykat_status"] = "cluster_failed"
        return call, info
    mse_f = mse[finite]
    low = mse_f <= np.nanmedian(mse_f)
    if int(low.sum()) < 10:
        info["copykat_status"] = "cluster_failed"
        return call, info
    fractions = []
    for cluster in range(int(labels_full.max()) + 1):
        members = labels_full == cluster
        fractions.append(float(low[members].sum()) / float(low.sum()))
    fractions = np.asarray(fractions)
    diploid_id = int(np.argmax(fractions))
    aneuploid_id = int(np.argmin(fractions))
    windows = features.T
    dip_profile = np.median(windows[:, labels_full == diploid_id], axis=1)
    ane_profile = np.median(windows[:, labels_full == aneuploid_id], axis=1)
    if np.std(dip_profile) == 0 or np.std(ane_profile) == 0:
        corr = 1.0
    else:
        corr = float(np.corrcoef(dip_profile, ane_profile)[0, 1])
    info["copykat_cor"] = corr
    if not math.isfinite(corr) or corr >= COR_ALL_DIPLOID:
        info["copykat_status"] = "all_diploid_cor>=0.6"
        return call, info
    ane_local = np.zeros(features.shape[0], dtype=bool)
    for cluster in range(int(labels_full.max()) + 1):
        profile = np.median(windows[:, labels_full == cluster], axis=1)
        dist_dip = float(stats.wasserstein_distance(profile, dip_profile))
        dist_ane = float(stats.wasserstein_distance(profile, ane_profile))
        if dist_ane < dist_dip and cluster != diploid_id:
            ane_local[labels_full == cluster] = True
    # Magnitude gate: above the normal-epithelium null. Homogeneous tumors can
    # still pass InferCNV when this correlation rule calls every cell diploid.
    ane_local &= np.isfinite(mse_f) & (mse_f > thr)
    call[np.where(finite)[0][ane_local]] = True
    info["copykat_status"] = "low_conf" if COR_LOW_CONF <= corr < COR_ALL_DIPLOID else "separated"
    return call, info


def assign_cnv_calls(records: list[dict]) -> None:
    """Null = median 5 Mb profile of normal epithelial cells in the same cohort.

    GSE189357 has no normal epithelium on GEO; it borrows the GSE123902
    normal-lung null and threshold. The threshold is the normal-epithelium
    95th percentile of residual MSE, fixed before the CLDN4 test.
    """
    bin_n = None
    for rec in records:
        prof = rec.get("profile")
        if prof is not None:
            bin_n = prof["infer"].shape[1]
            break
    if bin_n is None:
        return
    by_ds = defaultdict(list)
    for rec in records:
        by_ds[rec["meta"]["dataset"]].append(rec)

    def stack(rows, key):
        mats = [r["profile"][key] for r in rows if r.get("profile") is not None]
        if not mats:
            return None
        return np.vstack(mats)

    nulls = {}
    for dataset, rows in by_ds.items():
        controls = [r for r in rows if r["meta"]["role"] == "control" and r.get("profile") is not None]
        inf = stack(controls, "infer")
        cop = stack(controls, "copy")
        if inf is None or inf.shape[0] < 50:
            nulls[dataset] = None
            continue
        nulls[dataset] = {
            "infer": np.nanmedian(inf, axis=0),
            "copy": np.nanmedian(cop, axis=0),
            "n": int(inf.shape[0]),
        }
    if nulls.get("GSE123902") and nulls.get("GSE189357") is None and "GSE189357" in by_ds:
        nulls["GSE189357"] = dict(nulls["GSE123902"])
        nulls["GSE189357"]["borrowed_from"] = "GSE123902"
        say("GSE189357 borrows GSE123902 normal-lung CNV null", "n", nulls["GSE123902"]["n"])

    # Thresholds from control residuals, after the null is fixed.
    for dataset, rows in by_ds.items():
        null = nulls.get(dataset)
        if not null:
            say("no normal-epithelium null for", dataset)
            continue
        borrowed = null.get("borrowed_from")
        if borrowed:
            # Threshold computed below from the donor cohort's own controls.
            donor_controls = [r for r in by_ds[borrowed] if r["meta"]["role"] == "control" and r.get("profile") is not None]
        else:
            donor_controls = [r for r in rows if r["meta"]["role"] == "control" and r.get("profile") is not None]
        inf_c = stack(donor_controls, "infer")
        cop_c = stack(donor_controls, "copy")
        thr_inf = float(np.nanquantile(_mse(inf_c, null["infer"]), 0.95))
        thr_cop = float(np.nanquantile(_mse(cop_c, null["copy"]), 0.95))
        say("threshold", dataset, "infer", round(thr_inf, 5), "copykat", round(thr_cop, 5), "null_n", null["n"], "borrowed", borrowed or "")
        for rec in rows:
            prof = rec.get("profile")
            meta = rec["meta"]
            if prof is None:
                meta["cnv_ran"] = 0
                continue
            mse_inf = _mse(prof["infer"], null["infer"])
            mse_cop = _mse(prof["copy"], null["copy"])
            infer_call = np.isfinite(mse_inf) & (mse_inf > thr_inf)
            resid = prof["copy"] - null["copy"].reshape(1, -1)
            copy_call, info = _copykat_from_residual(resid, mse_cop, thr_cop)
            consensus = infer_call & copy_call
            meta.update(
                {
                    "cnv_ran": 1,
                    "n_genes": prof["n_genes"],
                    "n_epi_in_matrix": prof["n_epi"],
                    "n_ref_in_matrix": prof["n_ref"],
                    "copykat_status": info["copykat_status"],
                    "copykat_cor": info["copykat_cor"],
                    "infercnv_thr": thr_inf,
                    "copykat_thr": thr_cop,
                    "null_n": null["n"],
                    "null_borrowed_from": borrowed or "",
                    "median_infer_mse": float(np.nanmedian(mse_inf)),
                }
            )
            epi_annot = prof["is_annot"].astype(bool)
            for name, mask in (("copykat", copy_call), ("infercnv", infer_call), ("consensus", consensus)):
                pct, n = pct_on(prof["cldn4"], mask)
                meta[f"n_{name}"] = n
                meta[f"{name}_cldn4_pct"] = pct
                meta[f"frac_epi_{name}"] = float(mask.mean()) if len(mask) else float("nan")
                meta[f"frac_annot_{name}"] = float(mask[epi_annot].mean()) if epi_annot.any() else float("nan")
                other = ~epi_annot
                meta[f"frac_other_epi_{name}"] = float(mask[other].mean()) if other.any() else float("nan")


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def spearman(x, y) -> tuple[float, float]:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = len(x)
    if n < 4:
        return float("nan"), float("nan")
    rho = float(stats.spearmanr(x, y).statistic)
    if not math.isfinite(rho):
        return float("nan"), float("nan")
    if abs(rho) >= 1:
        return rho, 0.0
    t = rho * math.sqrt((n - 2) / (1 - rho * rho))
    p = float(2 * stats.t.sf(abs(t), n - 2))
    return rho, p


def dl_spearman(rhos, ns) -> dict:
    pairs = [(float(r), int(n)) for r, n in zip(rhos, ns) if int(n) > 3 and math.isfinite(float(r))]
    if len(pairs) < 2:
        return {"rho": float("nan"), "p": float("nan"), "I2": float("nan"), "ci_lo": float("nan"), "ci_hi": float("nan"), "k": len(pairs), "N": sum(n for _, n in pairs)}
    rhos = [r for r, _n in pairs]
    ns = [n for _r, n in pairs]
    z = [math.atanh(min(max(r, -0.999999), 0.999999)) for r in rhos]
    var = [1 / (n - 3) for n in ns]
    w = [1 / v for v in var]
    zbar = sum(wi * zi for wi, zi in zip(w, z)) / sum(w)
    q = sum(wi * (zi - zbar) ** 2 for wi, zi in zip(w, z))
    k = len(rhos)
    dfree = k - 1
    cdenom = sum(w) - sum(wi ** 2 for wi in w) / sum(w)
    tau2 = max(0.0, (q - dfree) / cdenom) if dfree > 0 and cdenom > 0 else 0.0
    wstar = [1 / (v + tau2) for v in var]
    zre = sum(wi * zi for wi, zi in zip(wstar, z)) / sum(wstar)
    se = math.sqrt(1 / sum(wstar))
    p = float(2 * stats.norm.sf(abs(zre / se)))
    i2 = max(0.0, (q - dfree) / q) if q > 0 else 0.0
    return {
        "rho": math.tanh(zre),
        "p": p,
        "I2": i2,
        "ci_lo": math.tanh(zre - 1.96 * se),
        "ci_hi": math.tanh(zre + 1.96 * se),
        "k": k,
        "N": int(sum(ns)),
    }


def within_quartile(values: np.ndarray) -> np.ndarray:
    ranks = pd.Series(values).rank(method="average").to_numpy()
    breaks = np.quantile(ranks, [0, 0.25, 0.5, 0.75, 1.0])
    labels = np.array(["NA"] * len(values), dtype=object)
    if len(np.unique(np.round(breaks, 6))) < 5:
        return labels
    cats = pd.cut(ranks, bins=breaks, include_lowest=True, labels=["Q1", "Q2", "Q3", "Q4"])
    return np.array(cats.astype(str))


def rank_biserial(x_q4, x_q1) -> dict:
    x = np.asarray(x_q4, float)
    y = np.asarray(x_q1, float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2:
        return {"r": float("nan"), "p": float("nan"), "n_q1": len(y), "n_q4": len(x)}
    both = np.concatenate([x, y])
    order = np.argsort(both, kind="mergesort")
    sorted_vals = both[order]
    ranks_sorted = np.empty(len(both), float)
    i = 0
    while i < len(both):
        j = i
        while j + 1 < len(both) and sorted_vals[j + 1] == sorted_vals[i]:
            j += 1
        ranks_sorted[i : j + 1] = 0.5 * ((i + 1) + (j + 1))
        i = j + 1
    ranks = np.empty(len(both), float)
    ranks[order] = ranks_sorted
    n4 = len(x)
    u = float(ranks[:n4].sum() - n4 * (n4 + 1) / 2.0)
    r = 2 * u / (n4 * len(y)) - 1
    p = float(stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic").pvalue)
    return {"r": r, "p": p, "n_q1": len(y), "n_q4": n4}


# ---------------------------------------------------------------------------
# Dataset prep
# ---------------------------------------------------------------------------

def select_laughney(rows: list[dict], tissues: set[str]) -> list[dict]:
    chosen = [r for r in rows if r["tissue"] in tissues]
    chosen = sorted(chosen, key=lambda r: (r["patient"], r["tissue"]))
    seen = set()
    out = []
    for row in chosen:
        if row["patient"] in seen:
            continue
        seen.add(row["patient"])
        if tissues == {"PRIMARY", "METASTASIS"}:
            if float(row["n_malignant"]) < 20 or float(row["n_tnk"]) < 20:
                continue
        out.append(row)
    return out


def prepare_dense_csv(path: Path, dataset: str, unit_id: str, role: str, order_map, cache: Path) -> None:
    if prepared_ok(cache, dataset, unit_id):
        say("  cached", dataset, unit_id)
        return
    say("  load", unit_id, path.name)
    frame = pd.read_csv(path, index_col=0)
    genes_raw = [str(c) for c in frame.columns]
    names, keep = first_occurrence(genes_raw)
    mat = frame.to_numpy(dtype=np.float32).T[keep]
    lib = mat.sum(axis=0).astype(np.float64)
    mal, tnk = marker_masks(names, mat)
    epi_all = np.where(mal)[0]
    ref_all = np.where(tnk)[0]
    annot_counts = mat[names.index("CLDN4"), epi_all] if "CLDN4" in names else np.zeros(len(epi_all))
    pct, mean = cldn4_stats(annot_counts)
    n_cells = int(mat.shape[1])
    n_tnk = int(tnk.sum())
    meta = {
        "dataset": dataset,
        "unit_id": unit_id,
        "role": role,
        "n_cells": n_cells,
        "n_tnk": n_tnk,
        "frac_tnk": n_tnk / n_cells if n_cells else float("nan"),
        "n_annot_mal": int(mal.sum()),
        "annot_cldn4_pct": pct,
        "annot_cldn4_mean": mean,
        "n_epi_available": int(mal.sum()),
        "n_ref_available": n_tnk,
    }
    rng = rng_for(dataset, unit_id)
    epi_idx = draw(epi_all, EPI_CAP, rng)
    ref_idx = draw(ref_all, REF_CAP, rng)
    meta["n_epi_drawn"] = int(len(epi_idx))
    meta["n_ref_drawn"] = int(len(ref_idx))
    payload = None
    if len(epi_idx) >= MIN_POS and len(ref_idx) >= MIN_REF:
        payload = pack_unit(names, mat, lib, epi_idx, ref_idx, epi_all, order_map)
    save_prepared(cache, meta, payload)
    say(
        "  ",
        role,
        unit_id,
        "cells",
        n_cells,
        "annot",
        meta["n_annot_mal"],
        "pct",
        round(pct, 2) if math.isfinite(pct) else None,
        "tnk",
        round(meta["frac_tnk"], 3),
    )


def prepare_gse123902(raw: Path, order_map, cache: Path) -> None:
    rows = read_tsv(ROOT / "data" / "GSE123902_marker_units.tsv")
    dest = raw / "gse123902"
    if not any(dest.glob("*dense.csv.gz")):
        dest.mkdir(parents=True, exist_ok=True)
        subprocess.check_call(["tar", "-C", str(dest), "-xf", str(raw / "GSE123902_RAW.tar")])
    tests = select_laughney(rows, {"PRIMARY", "METASTASIS"})
    controls = [r for r in rows if r["tissue"] == "NORMAL"]
    say("GSE123902 test", len(tests), "normal-controls", len(controls))
    if len(tests) != 13:
        raise SystemExit(f"GSE123902 test units {len(tests)} != 13")
    for row in tests:
        prepare_dense_csv(dest / row["file"], "GSE123902", row["patient"], "test", order_map, cache)
    for row in controls:
        prepare_dense_csv(dest / row["file"], "GSE123902", row["patient"] + "_NORMAL", "control", order_map, cache)


def prepare_gse189357(raw: Path, order_map, cache: Path) -> None:
    rows = read_tsv(ROOT / "data" / "GSE189357_marker_units.tsv")
    dest = raw / "gse189357"
    if not any(dest.glob("*_matrix.mtx.gz")):
        dest.mkdir(parents=True, exist_ok=True)
        subprocess.check_call(["tar", "-C", str(dest), "-xf", str(raw / "GSE189357_RAW.tar")])
    say("GSE189357 patients", len(rows))
    if len(rows) != 9:
        raise SystemExit(f"GSE189357 units {len(rows)} != 9")
    for row in rows:
        pat = row["patient"]
        if prepared_ok(cache, "GSE189357", pat):
            say("  cached", pat)
            continue
        mtx = next(dest.glob(f"*_{pat}_matrix.mtx.gz"))
        feat = next(dest.glob(f"*_{pat}_features.tsv.gz"))
        say("  load", pat)
        features = pd.read_csv(feat, sep="\t", header=None)
        genes_raw = features.iloc[:, 1 if features.shape[1] > 1 else 0].astype(str).tolist()
        names, keep = first_occurrence(genes_raw)
        with gzip.open(mtx, "rb") as handle:
            sparse = mmread(handle).tocsr()
        if sparse.shape[0] != len(genes_raw):
            raise SystemExit(f"{pat} mtx shape {sparse.shape} vs features {len(genes_raw)}")
        sparse = sparse[keep]
        mal, tnk = marker_masks_sparse(names, sparse)
        lib = np.asarray(sparse.sum(axis=0)).ravel().astype(np.float64)
        epi_all = np.where(mal)[0]
        ref_all = np.where(tnk)[0]
        if "CLDN4" in names:
            annot_counts = np.asarray(sparse[names.index("CLDN4"), epi_all].todense()).ravel()
        else:
            annot_counts = np.zeros(len(epi_all))
        pct, mean = cldn4_stats(annot_counts)
        n_cells = int(sparse.shape[1])
        n_tnk = int(tnk.sum())
        meta = {
            "dataset": "GSE189357",
            "unit_id": pat,
            "role": "test",
            "n_cells": n_cells,
            "n_tnk": n_tnk,
            "frac_tnk": n_tnk / n_cells if n_cells else float("nan"),
            "n_annot_mal": int(mal.sum()),
            "annot_cldn4_pct": pct,
            "annot_cldn4_mean": mean,
            "n_epi_available": int(mal.sum()),
            "n_ref_available": n_tnk,
        }
        rng = rng_for("GSE189357", pat)
        epi_idx = draw(epi_all, EPI_CAP, rng)
        ref_idx = draw(ref_all, REF_CAP, rng)
        meta["n_epi_drawn"] = int(len(epi_idx))
        meta["n_ref_drawn"] = int(len(ref_idx))
        payload = None
        if len(epi_idx) >= MIN_POS and len(ref_idx) >= MIN_REF:
            # Densify only the ordered genes x drawn cells.
            _names, chroms, row_idx = order_indices(names, order_map)
            cells = np.concatenate([epi_idx, ref_idx])
            dense = np.asarray(sparse[row_idx][:, cells].todense(), dtype=np.float32)
            keep_nz = dense.sum(axis=1) > 0
            n_epi = len(epi_idx)
            is_epi = np.zeros(len(cells), dtype=np.uint8)
            is_ref = np.zeros(len(cells), dtype=np.uint8)
            is_epi[:n_epi] = 1
            is_ref[n_epi:] = 1
            is_annot = is_epi.copy()
            cldn4 = (
                np.asarray(sparse[names.index("CLDN4"), cells].todense()).ravel().astype(np.float32)
                if "CLDN4" in names
                else np.zeros(len(cells), np.float32)
            )
            kept_i = np.where(keep_nz)[0]
            payload = {
                "counts": dense[keep_nz],
                "chrom": [chroms[i] for i in kept_i],
                "genes": [_names[i] for i in kept_i],
                "libsize": lib[cells],
                "is_epi": is_epi,
                "is_ref": is_ref,
                "is_annot": is_annot,
                "cldn4": cldn4,
            }
        save_prepared(cache, meta, payload)
        say("  test", pat, "cells", n_cells, "annot", meta["n_annot_mal"], "pct", round(pct, 2), "tnk", round(meta["frac_tnk"], 3))


def marker_masks_sparse(genes: list[str], mat) -> tuple[np.ndarray, np.ndarray]:
    index = {g: i for i, g in enumerate(genes)}

    def gt0(name: str) -> np.ndarray:
        i = index.get(name)
        if i is None:
            return np.zeros(mat.shape[1], dtype=bool)
        return np.asarray(mat[i].todense()).ravel() > 0

    mal = (gt0("EPCAM") | gt0("KRT8") | gt0("KRT18") | gt0("KRT19")) & ~gt0("PTPRC")
    tnk = (gt0("CD3D") | gt0("CD3E") | gt0("CD8A") | gt0("NKG7") | gt0("GNLY") | gt0("KLRD1")) & ~mal
    return mal, tnk


def prepare_gse131907(raw: Path, order_map, cache: Path) -> None:
    locked = [r for r in read_tsv(ROOT / "data" / "patient_units.tsv") if r["dataset"] == "GSE131907"]
    test_ids = [r["unit_id"] for r in locked]
    if len(test_ids) != 21:
        raise SystemExit(f"GSE131907 locked units {len(test_ids)} != 21")
    if all(prepared_ok(cache, "GSE131907", uid) for uid in test_ids):
        # Controls may still be missing; require the four largest normal lungs too.
        say("GSE131907 test cache complete; checking controls")
    ann_path = raw / "GSE131907_ann.txt.gz"
    mat_path = raw / "GSE131907_umi.txt.gz"
    say("GSE131907 read annotation")
    cells = []
    with gzip.open(ann_path, "rt") as handle:
        for rec in csv.DictReader(handle, delimiter="\t"):
            sample = rec["Sample"]
            origin = rec["Sample_Origin"]
            if sample not in test_ids and origin != "nLung":
                continue
            cells.append(
                {
                    "index": rec["Index"],
                    "sample": sample,
                    "origin": origin,
                    "epi": rec["Cell_type"] == "Epithelial cells",
                    "ref": rec["Cell_type"] in {"T lymphocytes", "NK cells"},
                    "mal": rec["Cell_subtype"] == "Malignant cells",
                }
            )
    by_sample = defaultdict(list)
    for i, rec in enumerate(cells):
        rec["all_i"] = i
        by_sample[rec["sample"]].append(i)
    # Draw CNV cells.
    drawn_all_i = []
    sample_of_drawn = []
    for sample in sorted(by_sample):
        role = "test" if sample in test_ids else "control"
        group = by_sample[sample]
        epi = np.array([i for i in group if cells[i]["epi"]], dtype=int)
        ref = np.array([i for i in group if cells[i]["ref"]], dtype=int)
        # Smaller draw than the other cohorts: the genes×cells text is held as a
        # dense autosomal matrix while it is streamed (208,506 cells).
        rng = rng_for("GSE131907", sample)
        epi_d = draw(epi, 500, rng)
        ref_d = draw(ref, 150, rng)
        for i in list(epi_d) + list(ref_d):
            drawn_all_i.append(int(i))
            sample_of_drawn.append(sample)
        say("  plan", role, sample, "epi", len(epi), "->", len(epi_d), "ref", len(ref), "->", len(ref_d))
    drawn_all_i = np.asarray(drawn_all_i, dtype=int)
    # Map matrix columns.
    say("GSE131907 map header")
    with gzip.open(mat_path, "rb") as handle:
        header = handle.readline().rstrip(b"\n").split(b"\t")
    col_of = {bc.decode(): i for i, bc in enumerate(header)}
    missing = 0
    for rec in cells:
        col = col_of.get(rec["index"])
        if col is None:
            missing += 1
            rec["col"] = -1
        else:
            rec["col"] = col
    say("  annotation cells kept", len(cells), "missing from matrix", missing)
    cnv_cols = []
    cnv_all_i = []
    for all_i in drawn_all_i:
        col = cells[int(all_i)]["col"]
        if col < 0:
            continue
        cnv_all_i.append(int(all_i))
        cnv_cols.append(col)
    n_cnv = len(cnv_cols)
    n_genes = len(order_map)
    say("  cnv cells", n_cnv, "ordered genes", n_genes, "memmap GB", round(n_cnv * n_genes * 4 / 1e9, 2))
    if n_cnv * n_genes * 4 > 3.2e9:
        raise SystemExit("memmap would exceed 3.2 GB; lower caps")
    order_symbols = sorted(order_map, key=lambda s: (order_map[s][0], order_map[s][1], s))
    gene_row = {sym: i for i, sym in enumerate(order_symbols)}
    chrom_of_row = np.array([order_map[s][0] for s in order_symbols], dtype=np.int16)
    mm_path = cache / "gse131907_memmap.dat"
    cache.mkdir(parents=True, exist_ok=True)
    if mm_path.exists():
        mm_path.unlink()
    mem = np.memmap(mm_path, dtype=np.float32, mode="w+", shape=(n_genes, n_cnv))
    lib = np.zeros(n_cnv, dtype=np.float64)
    cldn4 = np.zeros(len(cells), dtype=np.float32)
    col_to_j = {col: j for j, col in enumerate(cnv_cols)}
    all_cols = [(rec["col"], rec["all_i"]) for rec in cells if rec["col"] >= 0]
    say("GSE131907 stream UMI")
    t0 = time.time()
    n_seen = 0
    seen = set()
    with gzip.open(mat_path, "rb") as handle:
        handle.readline()
        for line in handle:
            parts = line.rstrip(b"\n").split(b"\t")
            gene = parts[0].decode().upper()
            if gene in seen:
                continue
            seen.add(gene)
            n_seen += 1
            row = gene_row.get(gene)
            is_cl = gene == "CLDN4"
            if row is None and not is_cl:
                # Still need library size.
                for col, j in col_to_j.items():
                    token = parts[col]
                    if token != b"0":
                        lib[j] += float(token)
            else:
                if is_cl:
                    for col, all_i in all_cols:
                        token = parts[col]
                        if token != b"0":
                            cldn4[all_i] = float(token)
                if row is not None:
                    for col, j in col_to_j.items():
                        token = parts[col]
                        if token != b"0":
                            val = float(token)
                            mem[row, j] = val
                            lib[j] += val
                elif is_cl:
                    for col, j in col_to_j.items():
                        token = parts[col]
                        if token != b"0":
                            lib[j] += float(token)
            if n_seen % 4000 == 0:
                say("  genes", n_seen, "sec", round(time.time() - t0, 1))
    say("  stream done genes", n_seen, "sec", round(time.time() - t0, 1), "CLDN4 sum", float(cldn4.sum()))
    mem.flush()
    # Split by sample. cnv_all_i order matches columns.
    by_col_sample = defaultdict(list)
    for j, all_i in enumerate(cnv_all_i):
        by_col_sample[cells[all_i]["sample"]].append(j)
    for sample in sorted(by_sample):
        role = "test" if sample in set(test_ids) else "control"
        group = by_sample[sample]
        mal_i = [i for i in group if cells[i]["mal"] and cells[i]["col"] >= 0]
        epi_i = [i for i in group if cells[i]["epi"] and cells[i]["col"] >= 0]
        ref_i = [i for i in group if cells[i]["ref"] and cells[i]["col"] >= 0]
        n_cells = len(group)
        n_tnk = len(ref_i)
        pct, mean = cldn4_stats(cldn4[mal_i] if mal_i else np.zeros(0))
        cols = by_col_sample.get(sample, [])
        meta = {
            "dataset": "GSE131907",
            "unit_id": sample,
            "role": role,
            "n_cells": n_cells,
            "n_tnk": n_tnk,
            "frac_tnk": n_tnk / n_cells if n_cells else float("nan"),
            "n_annot_mal": len(mal_i),
            "annot_cldn4_pct": pct,
            "annot_cldn4_mean": mean,
            "n_epi_available": len(epi_i),
            "n_ref_available": n_tnk,
            "n_epi_drawn": int(sum(1 for j in cols if cells[cnv_all_i[j]]["epi"])),
            "n_ref_drawn": int(sum(1 for j in cols if cells[cnv_all_i[j]]["ref"])),
        }
        payload = None
        if cols and meta["n_epi_drawn"] >= MIN_POS and meta["n_ref_drawn"] >= MIN_REF:
            block = np.array(mem[:, cols], dtype=np.float32)
            keep = block.sum(axis=1) > 0
            drawn_all = [cnv_all_i[j] for j in cols]
            is_epi = np.array([1 if cells[i]["epi"] else 0 for i in drawn_all], dtype=np.uint8)
            is_ref = np.array([1 if cells[i]["ref"] and not cells[i]["epi"] else 0 for i in drawn_all], dtype=np.uint8)
            is_annot = np.array([1 if cells[i]["mal"] else 0 for i in drawn_all], dtype=np.uint8)
            # Epi columns first so nothing downstream assumes order, masks carry the roles.
            payload = {
                "counts": block[keep],
                "chrom": chrom_of_row[keep].tolist(),
                "genes": [order_symbols[i] for i in np.where(keep)[0]],
                "libsize": lib[cols],
                "is_epi": is_epi,
                "is_ref": is_ref,
                "is_annot": is_annot,
                "cldn4": cldn4[drawn_all],
            }
        save_prepared(cache, meta, payload)
        say("  saved", role, sample, "annot_pct", None if not math.isfinite(pct) else round(pct, 2), "tnk", round(meta["frac_tnk"], 3))
    del mem
    mm_path.unlink(missing_ok=True)


def parse_soft(path: Path) -> dict[str, dict]:
    rows = []
    cur = {}
    for line in gzip.open(path, "rt", errors="replace"):
        if line.startswith("^SAMPLE"):
            if cur.get("title"):
                rows.append(cur)
            cur = {}
        elif line.startswith("!Sample_title"):
            cur["title"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_characteristics_ch1"):
            val = line.split("=", 1)[1].strip()
            if ": " in val:
                key, value = val.split(": ", 1)
                cur[key] = value
    if cur.get("title"):
        rows.append(cur)
    out = {}
    for row in rows:
        title = row.get("title", "")
        code = title.split(" ", 1)[1] if " " in title else title
        code = code.upper().replace("-", "_")
        out[code] = {"patient": row.get("patient", ""), "tissue": row.get("tissue", ""), "title": title}
    return out


def norm_code(orig: str) -> str:
    import re

    text = orig.upper().replace("-", "_")
    return re.sub(r"_[35]P$", "", text)


def prepare_gse205335(raw: Path, order_map, cache: Path) -> None:
    patients = [r["patient"] for r in read_tsv(ROOT / "data" / "GSE205335_patients.tsv")]
    if len(patients) != 22:
        raise SystemExit(f"GSE205335 patients {len(patients)} != 22")
    mat_dir = cache / "GSE205335_mats"
    flags_path = cache / "gse205335_flags.tsv"
    need = [mat_dir / f"{pat}.mtx" for pat in patients] + [mat_dir / "CONTROL.mtx", mat_dir / "cldn4_mal.tsv"]
    if not all(p.exists() for p in need):
        say("GSE205335 build flags and export RDS")
        soft = parse_soft(raw / "GSE205335_family.soft.gz")
        ident = pd.read_csv(raw / "GSE205335_ident.txt.gz", sep="\t")
        ident["code"] = ident["orig.ident"].map(norm_code)
        ident = ident.merge(
            pd.DataFrame([{"code": k, **v} for k, v in soft.items()]),
            on="code",
            how="left",
        )
        unmatched = ident.loc[ident["patient"].isna(), "orig.ident"].value_counts()
        if len(unmatched):
            say("  unmatched orig.ident", unmatched.head(8).to_dict())
        ident["is_normal"] = ident["tissue"].fillna("").str.startswith("Normal")
        ident["epi"] = ident["lineage.total"] == "Epithelial cells"
        ident["ref"] = ident["lineage.total"] == "T/NK cells"
        ident["mal"] = ident["lineage.sub"] == "Malignant cells"
        flag_rows = []
        # Control: normal-tissue epithelial + T/NK, capped.
        ctrl = ident[ident["is_normal"]].copy()
        rng = rng_for("GSE205335", "CONTROL")
        epi_c = draw(np.where(ctrl["epi"].to_numpy())[0], EPI_CAP, rng)
        ref_c = draw(np.where(ctrl["ref"].to_numpy())[0], REF_CAP, rng)
        ctrl_idx = np.concatenate([epi_c, ref_c]) if len(epi_c) or len(ref_c) else np.array([], dtype=int)
        for i in ctrl_idx:
            row = ctrl.iloc[int(i)]
            flag_rows.append(
                {
                    "barcode": row["barcode"],
                    "patient": "CONTROL",
                    "epi": int(row["epi"]),
                    "ref": int(row["ref"] and not row["epi"]),
                    "mal": 0,
                }
            )
        say("  control epi", len(epi_c), "ref", len(ref_c), "of normal-tissue", int(ctrl["epi"].sum()))
        for pat in patients:
            sub = ident[(ident["patient"] == pat) & (~ident["is_normal"])]
            if sub.empty:
                raise SystemExit(f"no tumor cells for {pat}")
            epi = np.where(sub["epi"].to_numpy())[0]
            ref = np.where(sub["ref"].to_numpy())[0]
            mal = np.where(sub["mal"].to_numpy())[0]
            rng = rng_for("GSE205335", pat)
            epi_d = set(int(i) for i in draw(epi, EPI_CAP, rng))
            ref_d = set(int(i) for i in draw(ref, REF_CAP, rng))
            mal_set = set(int(i) for i in mal)
            # Matrix cells are the CNV draw. CLDN4 is extracted for every malignant cell.
            for i in range(len(sub)):
                if i not in epi_d and i not in ref_d and i not in mal_set:
                    continue
                row = sub.iloc[i]
                flag_rows.append(
                    {
                        "barcode": row["barcode"],
                        "patient": pat,
                        "epi": int(i in epi_d),
                        "ref": int(i in ref_d and i not in epi_d),
                        "mal": int(i in mal_set),
                    }
                )
        flags_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(flag_rows).to_csv(flags_path, sep="\t", index=False)
        rds_script = ROOT / "scripts" / "export_gse205335.R"
        subprocess.check_call(["Rscript", str(rds_script), str(raw), str(mat_dir), str(flags_path)])
    cldn4_tbl = pd.read_csv(mat_dir / "cldn4_mal.tsv", sep="\t")
    ident = pd.read_csv(raw / "GSE205335_ident.txt.gz", sep="\t")
    soft = parse_soft(raw / "GSE205335_family.soft.gz")
    ident["code"] = ident["orig.ident"].map(norm_code)
    ident = ident.merge(pd.DataFrame([{"code": k, **v} for k, v in soft.items()]), on="code", how="left")
    ident["is_normal"] = ident["tissue"].fillna("").str.startswith("Normal")
    ident["ref"] = ident["lineage.total"] == "T/NK cells"
    ident["mal"] = ident["lineage.sub"] == "Malignant cells"
    ident["epi"] = ident["lineage.total"] == "Epithelial cells"
    for pat in patients + ["CONTROL"]:
        role = "control" if pat == "CONTROL" else "test"
        if prepared_ok(cache, "GSE205335", pat):
            say("  cached", pat)
            continue
        if pat == "CONTROL":
            sub = ident[ident["is_normal"]]
        else:
            sub = ident[(ident["patient"] == pat) & (~ident["is_normal"])]
        n_cells = int(len(sub))
        n_tnk = int(sub["ref"].sum())
        n_mal = int(sub["mal"].sum()) if pat != "CONTROL" else 0
        mal_counts = cldn4_tbl.loc[cldn4_tbl["patient"] == pat, "cldn4"].to_numpy() if pat != "CONTROL" else np.zeros(0)
        pct, mean = cldn4_stats(mal_counts)
        flags = pd.read_csv(mat_dir / f"{pat}.flags.tsv", sep="\t")
        genes = [g.upper() for g in (mat_dir / f"{pat}.genes").read_text().splitlines() if g]
        names, keep_u = first_occurrence(genes)
        with (mat_dir / f"{pat}.mtx").open("rb") as handle:
            sparse = mmread(handle).tocsr()
        if sparse.shape[1] != len(flags) or sparse.shape[0] != len(genes):
            raise SystemExit(f"{pat} mtx {sparse.shape} flags {len(flags)} genes {len(genes)}")
        sparse = sparse[keep_u]
        _names, chroms, row_idx = order_indices(names, order_map)
        dense = np.asarray(sparse[row_idx].todense(), dtype=np.float32)
        keep_nz = dense.sum(axis=1) > 0
        lib = np.asarray(sparse.sum(axis=0)).ravel().astype(np.float64)
        # CLDN4 on the drawn cells from the exported matrix (gene may have been dropped if all zero).
        if "CLDN4" in names:
            cldn4 = np.asarray(sparse[names.index("CLDN4")].todense()).ravel().astype(np.float32)
        else:
            cldn4 = np.zeros(sparse.shape[1], np.float32)
        meta = {
            "dataset": "GSE205335",
            "unit_id": pat,
            "role": role,
            "n_cells": n_cells,
            "n_tnk": n_tnk,
            "frac_tnk": n_tnk / n_cells if n_cells else float("nan"),
            "n_annot_mal": n_mal,
            "annot_cldn4_pct": pct if pat != "CONTROL" else float("nan"),
            "annot_cldn4_mean": mean if pat != "CONTROL" else float("nan"),
            "n_epi_available": int(sub["epi"].sum()),
            "n_ref_available": n_tnk,
            "n_epi_drawn": int(flags["epi"].sum()),
            "n_ref_drawn": int(flags["ref"].sum()),
        }
        kept_i = np.where(keep_nz)[0]
        payload = {
            "counts": dense[keep_nz],
            "chrom": [chroms[i] for i in kept_i],
            "genes": [_names[i] for i in kept_i],
            "libsize": lib,
            "is_epi": flags["epi"].to_numpy(dtype=np.uint8),
            "is_ref": flags["ref"].to_numpy(dtype=np.uint8),
            "is_annot": flags["mal"].to_numpy(dtype=np.uint8) if "mal" in flags.columns else np.zeros(len(flags), np.uint8),
            "cldn4": cldn4,
        }
        # Reference library size must include genes outside the ordered set.
        # sparse was already restricted to unique symbols, which is essentially the whole matrix.
        save_prepared(cache, meta, payload)
        shown = round(pct, 2) if math.isfinite(pct) and pat != "CONTROL" else None
        say("  saved", role, pat, "cells", n_cells, "annot_pct", shown, "tnk", round(meta["frac_tnk"], 3))


# ---------------------------------------------------------------------------
# Summaries and figures
# ---------------------------------------------------------------------------

SCORE_COLUMNS = {
    "annotation": "annot_cldn4_pct",
    "copykat": "copykat_cldn4_pct",
    "infercnv": "infercnv_cldn4_pct",
    "consensus": "consensus_cldn4_pct",
}


def cohort_rows(frame: pd.DataFrame, definition: str, paired_to: str | None = None) -> list[dict]:
    col = SCORE_COLUMNS[definition if paired_to is None else "annotation"]
    if paired_to:
        gate = SCORE_COLUMNS[paired_to]
        use = frame[np.isfinite(frame[gate].to_numpy(dtype=float))].copy()
        label = f"annotation_on_{paired_to}_units"
    else:
        use = frame.copy()
        label = definition
        col = SCORE_COLUMNS[definition]
    rows = []
    rhos, ns = [], []
    q4_all, q1_all = [], []
    for cohort in COHORTS:
        sub = use[use["dataset"] == cohort]
        x = sub[col].to_numpy(dtype=float)
        y = sub["frac_tnk"].to_numpy(dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)
        rho, p = spearman(x[ok], y[ok])
        n = int(ok.sum())
        rhos.append(rho)
        ns.append(n)
        rows.append({"definition": label, "cohort": cohort, "n": n, "rho": rho, "p": p, "I2": "", "ci_lo": "", "ci_hi": "", "N": "", "q4q1_r": "", "q4q1_p": "", "n_q1": "", "n_q4": ""})
        if n >= 4:
            labs = within_quartile(x[ok])
            yy = y[ok]
            q4_all.extend(yy[labs == "Q4"].tolist())
            q1_all.extend(yy[labs == "Q1"].tolist())
    meta = dl_spearman(rhos, ns)
    rb = rank_biserial(q4_all, q1_all)
    rows.append(
        {
            "definition": label,
            "cohort": "META",
            "n": meta["k"],
            "rho": meta["rho"],
            "p": meta["p"],
            "I2": meta["I2"],
            "ci_lo": meta["ci_lo"],
            "ci_hi": meta["ci_hi"],
            "N": meta["N"],
            "q4q1_r": rb["r"],
            "q4q1_p": rb["p"],
            "n_q1": rb["n_q1"],
            "n_q4": rb["n_q4"],
        }
    )
    return rows


def make_figures(units: pd.DataFrame, stats_rows: list[dict], fig_dir: Path) -> None:
    fig_dir.mkdir(parents=True, exist_ok=True)
    test = units[units["role"] == "test"].copy()
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), sharey=True)
    for ax, col, title in (
        (axes[0], "annot_cldn4_pct", "Annotation-only malignant"),
        (axes[1], "copykat_cldn4_pct", "CopyKAT-aneuploid epithelium"),
    ):
        for cohort in COHORTS:
            sub = test[test["dataset"] == cohort]
            ax.scatter(sub[col], sub["frac_tnk"], s=28, color=COLORS[cohort], label=cohort, alpha=0.9)
        ax.set_xlabel("CLDN4 % positive")
        ax.set_title(title)
        ax.set_xlim(left=-2)
    axes[0].set_ylabel("T/NK fraction")
    axes[1].legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig_cldn4_vs_tnk.png", dpi=140)
    fig.savefig(fig_dir / "fig_cldn4_vs_tnk.pdf")
    plt.close(fig)

    defs = ["annotation", "copykat", "infercnv"]
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    xpos = np.arange(len(COHORTS) + 1)
    width = 0.18
    for i, definition in enumerate(defs):
        ys = []
        for cohort in COHORTS + ["META"]:
            hit = [r for r in stats_rows if r["definition"] == definition and r["cohort"] == cohort]
            ys.append(hit[0]["rho"] if hit else np.nan)
        ax.bar(xpos + (i - 1) * width, ys, width=width, label=definition)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(xpos)
    ax.set_xticklabels(COHORTS + ["DL meta"])
    ax.set_ylabel("Spearman ρ  CLDN4 %pos vs T/NK")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig_rho_compare.png", dpi=140)
    fig.savefig(fig_dir / "fig_rho_compare.pdf")
    plt.close(fig)


def qc_annotation(units: pd.DataFrame) -> pd.DataFrame:
    locked = pd.read_csv(ROOT / "data" / "patient_units.tsv", sep="\t")
    test = units[units["role"] == "test"]
    merged = test.merge(locked, on=["dataset", "unit_id"], suffixes=("", "_locked"))
    rows = []
    for cohort in COHORTS:
        sub = merged[merged["dataset"] == cohort]
        rho, p = spearman(sub["annot_cldn4_pct"], sub["frac_tnk"])
        d_pct = np.nanmax(np.abs(sub["annot_cldn4_pct"] - sub["mal_CLDN4_pct"])) if len(sub) else float("nan")
        d_tnk = np.nanmax(np.abs(sub["frac_tnk"] - sub["frac_tnk_locked"])) if len(sub) else float("nan")
        rows.append(
            {
                "cohort": cohort,
                "n": len(sub),
                "rho": rho,
                "p": p,
                "locked_rho": LOCKED_RHO[cohort],
                "delta_rho": rho - LOCKED_RHO[cohort],
                "max_abs_pct": d_pct,
                "max_abs_frac_tnk": d_tnk,
            }
        )
        say(
            "QC",
            cohort,
            "n",
            len(sub),
            "rho",
            None if not math.isfinite(rho) else round(rho, 3),
            "locked",
            LOCKED_RHO[cohort],
            "max|pct|",
            None if not math.isfinite(d_pct) else round(d_pct, 3),
            "max|tnk|",
            None if not math.isfinite(d_tnk) else round(float(d_tnk), 5),
        )
    return pd.DataFrame(rows)


def load_units(cache: Path, datasets: list[str], order_map) -> pd.DataFrame:
    bin_keys = genome_bins(order_map)
    say("5 Mb bins", len(bin_keys), "(chr6 excluded)")
    records = []
    for dataset in datasets:
        folder = cache / dataset
        if not folder.exists():
            continue
        for path in sorted(folder.glob("*.json")):
            meta = json.loads(path.read_text())
            prof = epithelial_profiles(folder / f"{path.stem}.npz", order_map, bin_keys)
            records.append({"meta": meta, "profile": prof})
            say("  profile", dataset, meta["unit_id"], "epi", None if prof is None else prof["n_epi"])
    assign_cnv_calls(records)
    rows = []
    for rec in records:
        meta = rec["meta"]
        rows.append(meta)
        if meta.get("role") == "test" and meta.get("cnv_ran"):
            say(
                "call",
                meta["dataset"],
                meta["unit_id"],
                meta.get("copykat_status"),
                "cor",
                None if not isinstance(meta.get("copykat_cor"), float) else round(meta["copykat_cor"], 3),
                "frac_inf",
                None if not isinstance(meta.get("frac_epi_infercnv"), float) else round(meta["frac_epi_infercnv"], 3),
                "frac_ck",
                None if not isinstance(meta.get("frac_epi_copykat"), float) else round(meta["frac_epi_copykat"], 3),
                "n_inf",
                meta.get("n_infercnv"),
                "n_both",
                meta.get("n_consensus"),
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", default="/tmp/concordant4_raw")
    parser.add_argument("--cache", default="/tmp/concordant4_cnv_cache/v2")
    parser.add_argument("--datasets", default=",".join(COHORTS))
    args = parser.parse_args()
    raw = Path(args.raw)
    cache = Path(args.cache)
    datasets = [d for d in args.datasets.split(",") if d]
    order_map = load_gene_order(raw / "refGene.txt.gz", ROOT / "data" / "gene_order_hg38.tsv")
    if "CLDN4" not in order_map:
        raise SystemExit("CLDN4 missing from gene order")
    say("CLDN4", order_map["CLDN4"])
    if "GSE123902" in datasets:
        prepare_gse123902(raw, order_map, cache)
    if "GSE189357" in datasets:
        prepare_gse189357(raw, order_map, cache)
    if "GSE205335" in datasets:
        prepare_gse205335(raw, order_map, cache)
    if "GSE131907" in datasets:
        prepare_gse131907(raw, order_map, cache)

    say("CNV calls")
    units = load_units(cache, datasets, order_map)
    tab = ROOT / "results" / "tables"
    fig = ROOT / "results" / "figures"
    tab.mkdir(parents=True, exist_ok=True)
    units.to_csv(tab / "unit_calls.tsv", sep="\t", index=False)
    qc = qc_annotation(units)
    qc.to_csv(tab / "annotation_qc.tsv", sep="\t", index=False)

    test = units[units["role"] == "test"].copy()
    stat_rows = []
    for definition in SCORE_COLUMNS:
        stat_rows.extend(cohort_rows(test, definition))
    for definition in ("copykat", "infercnv", "consensus"):
        stat_rows.extend(cohort_rows(test, definition, paired_to=definition))
    columns = ["definition", "cohort", "n", "rho", "p", "I2", "ci_lo", "ci_hi", "N", "q4q1_r", "q4q1_p", "n_q1", "n_q4"]
    write_tsv(tab / "cldn4_tnk_stats.tsv", stat_rows, columns)
    make_figures(units, stat_rows, fig)

    controls = units[units["role"] == "control"]
    if len(controls):
        ctrl_cols = [
            "dataset",
            "unit_id",
            "n_epi_drawn",
            "n_ref_drawn",
            "frac_epi_copykat",
            "frac_epi_infercnv",
            "frac_epi_consensus",
            "median_infer_mse",
            "copykat_status",
            "copykat_cor",
        ]
        controls.to_csv(tab / "normal_epithelium_controls.tsv", sep="\t", index=False, columns=[c for c in ctrl_cols if c in controls.columns])
    summary = {
        "epi_cap": EPI_CAP,
        "ref_cap": REF_CAP,
        "min_pos": MIN_POS,
        "infercnv_min_autosomes": MIN_BREADTH,
        "copykat_cor_all_diploid": COR_ALL_DIPLOID,
        "n_test_units": int((units["role"] == "test").sum()),
        "n_control_units": int((units["role"] == "control").sum()),
    }
    (ROOT / "results" / "summary.json").write_text(json.dumps(summary, indent=2))
    say("wrote", tab)
    for row in stat_rows:
        if row["cohort"] == "META":
            say("META", row["definition"], "N", row["N"], "rho", row["rho"], "p", row["p"], "I2", row["I2"])


if __name__ == "__main__":
    main()
