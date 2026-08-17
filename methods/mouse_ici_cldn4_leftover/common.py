"""Shared scoring for leftover mouse ICI Cldn4 (scRNA + spatial)."""
from __future__ import annotations

import gzip
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "notes/mouse_ici_cldn4_leftover/raw"
OUT = ROOT / "results/mouse_ici_cldn4_leftover"
PANEL_PATH = Path(__file__).with_name("gene_panel.tsv")

# Do not use Sftpc/Scgb1a1 (ambient / normal AT2-club).
EPI_MARKERS = ["Epcam", "Cdh1", "Krt8", "Krt18", "Krt19"]
TNK_MARKERS = ["Cd3d", "Cd3e", "Cd3g", "Cd8a", "Nkg7", "Ncr1"]
T_SCORE_GENES = ["Cd3d", "Cd3e", "Cd3g", "Cd2", "Cd8a", "Cd8b1", "Nkg7", "Gzmb", "Prf1", "Ifng"]
EXCLUSION_T = ["Cd3d", "Cd3e", "Cd8a", "Nkg7", "Gzmb", "Prf1"]

ALIASES = {
    "Epcam": ["Tacstd1"],
    "Nkx2-1": ["Nkx2.1", "Ttf1"],
    "Cd8b1": ["Cd8b"],
    "Cldn4": ["CLDN4"],
    "Tacstd2": ["TACSTD2", "Trop2"],
}


def load_panel():
    df = pd.read_csv(PANEL_PATH, sep="\t")
    roles = defaultdict(list)
    for _, r in df.iterrows():
        roles[r.role].append(r.symbol)
    return df["symbol"].tolist(), dict(roles)


def norm_symbol(s: str) -> str:
    return str(s).replace(".", "-").upper()


def gene_index(genes: list[str], wanted: list[str]) -> dict[str, int]:
    lut = {}
    for i, g in enumerate(genes):
        lut[norm_symbol(g)] = i
    out = {}
    for w in wanted:
        keys = [w] + ALIASES.get(w, [])
        for k in keys:
            i = lut.get(norm_symbol(k))
            if i is not None:
                out[w] = i
                break
    return out


def read_features(path: Path) -> list[str]:
    genes = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and not parts[1].isdigit():
                genes.append(parts[1])
            else:
                genes.append(parts[0])
    return genes


def mtx_shape(path: Path) -> tuple[int, int, int]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("%"):
                continue
            ng, nc, nz = map(int, line.split()[:3])
            return ng, nc, nz
    raise ValueError(path)


def stream_mtx_panel(mtx_path: Path, idx: dict[str, int], n_cells: int, min_umi: int | None):
    """One-pass MTX: per-cell UMI totals + panel counts. 10x is genes x cells, 1-indexed."""
    want = {i: name for name, i in idx.items()}
    umi = np.zeros(n_cells, dtype=np.int64)
    n_genes = np.zeros(n_cells, dtype=np.int32)
    mat = {name: np.zeros(n_cells, dtype=np.float32) for name in idx}
    opener = gzip.open if str(mtx_path).endswith(".gz") else open
    with opener(mtx_path, "rt") as f:
        header_done = False
        for line in f:
            if line.startswith("%"):
                continue
            if not header_done:
                header_done = True
                continue
            a, b, c = line.split()
            gi = int(a) - 1
            ci = int(b) - 1
            v = float(c)
            if 0 <= ci < n_cells:
                umi[ci] += v
                n_genes[ci] += 1
                if gi in want:
                    mat[want[gi]][ci] = v
    if min_umi is not None:
        keep = umi >= min_umi
    else:
        keep = umi > 0
    mat = {k: v[keep] for k, v in mat.items()}
    return mat, umi[keep], n_genes[keep], int(keep.sum()), int((~keep).sum()), keep


def _pos(mat, name, n):
    if name in mat:
        return mat[name] > 0
    return np.zeros(n, dtype=bool)


def classify(mat):
    """Epithelial/tumor = Epcam or (Cdh1 and Krt8) or SCLC markers. T/NK = Cd3/Cd8/Nkg7/Ncr1.

    Epcam/tumor markers win over ambient T UMIs. Do not require Ptprc==0.
    """
    n = len(next(iter(mat.values())))
    tnk = (
        _pos(mat, "Cd3e", n)
        | _pos(mat, "Cd3d", n)
        | _pos(mat, "Cd8a", n)
        | _pos(mat, "Nkg7", n)
        | _pos(mat, "Ncr1", n)
    )
    epi = (
        _pos(mat, "Epcam", n)
        | (_pos(mat, "Cdh1", n) & _pos(mat, "Krt8", n))
        | _pos(mat, "Ascl1", n)
        | _pos(mat, "Chga", n)
        | _pos(mat, "Insm1", n)
    )
    lab = np.array(["other"] * n, dtype=object)
    lab[tnk] = "tnk"
    lab[epi] = "epithelial"
    return lab


def log1p_mat(mat):
    return {k: np.log1p(v) for k, v in mat.items()}


def mean_of(mat, names, mask=None):
    arrs = [mat[n] for n in names if n in mat]
    if not arrs:
        return None
    stacked = np.vstack(arrs)
    if mask is None:
        return stacked.mean(axis=0)
    return stacked[:, mask].mean(axis=0)


def spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 8:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(x[m], y[m])
    return {"n": n, "rho": float(rho), "p": float(p)}


def welch_or_nan(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    t, p = stats.ttest_ind(a, b, equal_var=False)
    return float(t), float(p)


def mwu_or_nan(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    try:
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        return float(u), float(p)
    except ValueError:
        return np.nan, np.nan


def wilcoxon_signed(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    n = int(a.size)
    if n < 3:
        return {"n": n, "p": np.nan, "n_a_gt_b": int((a > b).sum()) if n else 0}
    try:
        p = float(stats.wilcoxon(a - b, alternative="two-sided").pvalue)
    except ValueError:
        p = np.nan
    return {"n": n, "p": p, "n_a_gt_b": int((a > b).sum())}


def high_low_cut(cldn4, immune, how="median"):
    """Boolean mask: Cldn4-high AND immune-low. how='median' or 'quartile'."""
    c = np.asarray(cldn4, float)
    im = np.asarray(immune, float)
    ok = np.isfinite(c) & np.isfinite(im)
    mask = np.zeros(c.size, dtype=bool)
    if ok.sum() < 8:
        return mask, {"n_ok": int(ok.sum()), "cut": how, "n_hi_lo": 0}
    if how == "quartile":
        c_cut = np.nanpercentile(c[ok], 75)
        i_cut = np.nanpercentile(im[ok], 25)
    else:
        c_cut = np.nanmedian(c[ok])
        i_cut = np.nanmedian(im[ok])
    mask[ok] = (c[ok] >= c_cut) & (im[ok] <= i_cut)
    return mask, {
        "n_ok": int(ok.sum()),
        "cut": how,
        "cldn4_cut": float(c_cut),
        "immune_cut": float(i_cut),
        "n_hi_lo": int(mask.sum()),
        "frac_hi_lo": float(mask.sum() / ok.sum()),
    }


def dump_json(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=str) + "\n")


def hex_neighbors(r, c):
    """Visium / hex grid neighbors (even-row offset)."""
    if int(r) % 2 == 0:
        return [(r - 1, c - 1), (r - 1, c), (r, c - 1), (r, c + 1), (r + 1, c - 1), (r + 1, c)]
    return [(r - 1, c), (r - 1, c + 1), (r, c - 1), (r, c + 1), (r + 1, c), (r + 1, c + 1)]
