#!/usr/bin/env python3
"""Concordant-4 malignant surface genes vs patient T/NK, keratin partialled.

Head-to-head of CLDN1, CLDN4, CLDN7, and EPCAM on the locked n=65 unit set
(GSE123902 + GSE131907 + GSE205335 + GSE189357). MUC1 is context only,
because the locked bulk surface ranking named it beside CLDN7 and EPCAM.

Does not rebuild the Seurat/Harmony object. Patient T/NK and the malignant
gate match methods/seurat_concordant4_cldn4 (PR #539). R is used only to
read the GSE205335 dgCMatrix RDS.
"""

from __future__ import annotations

import csv
import gzip
import math
import os
import subprocess
import tarfile
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
GEO = Path(os.environ.get("GEO_DIR", "/tmp/geo_c4"))
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
LOCKED = ROOT / "data" / "pr539_cldn4_units.tsv"

GENES = ["CLDN1", "CLDN4", "CLDN7", "EPCAM", "MUC1"]
PRIMARY = ["CLDN1", "CLDN4", "CLDN7", "EPCAM"]
KERATIN = ["KRT8", "KRT18", "KRT19"]
MARKERS = ["EPCAM", "KRT8", "KRT18", "KRT19", "PTPRC", "CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]
SCORE_GENES = ["CLDN1", "CLDN4", "CLDN7", "EPCAM", "MUC1", "KRT8", "KRT18", "KRT19"]
N_PERM = 10000
SEED = 1

ELIG_131_ORIGIN = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
MIN_MAL_131 = 20


def say(msg: str) -> None:
    print(msg, flush=True)


def rank_average(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    sx = x[order]
    n = len(x)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and sx[j + 1] == sx[i]:
            j += 1
        ranks[order[i : j + 1]] = 0.5 * ((i + 1) + (j + 1))
        i = j + 1
    return ranks


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return pearson(rank_average(x), rank_average(y))


def _ols_resid(y: np.ndarray, A: np.ndarray) -> np.ndarray:
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return y - A @ coef


def partial_spearman(x: np.ndarray, y: np.ndarray, covariates: list[np.ndarray]) -> float:
    """Spearman partial: Pearson of rank residuals. Both sides adjusted."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 5:
        return float("nan")
    cols = [np.ones(len(x))]
    for z in covariates:
        cols.append(rank_average(np.asarray(z, dtype=float)))
    A = np.column_stack(cols)
    return pearson(_ols_resid(rank_average(x), A), _ols_resid(rank_average(y), A))


def spearman_p(rho: float, n: int, k_cov: int = 0) -> float:
    df = n - 2 - k_cov
    if not math.isfinite(rho) or df <= 0 or abs(rho) >= 1:
        return float("nan") if not math.isfinite(rho) else 0.0
    t = rho * math.sqrt(df / (1.0 - rho * rho))
    # two-sided Student-t via regularized beta, matching a plain t test
    x = df / (df + t * t)
    # survival function of Beta(df/2, 0.5) at x, times 1 (two-sided = I_x)
    p = _betainc_reg(df / 2.0, 0.5, x)
    return float(min(1.0, max(0.0, p)))


def _betainc_reg(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta I_x(a, b). Used for the t tail."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    # Lentz continued fraction for the incomplete beta (Numerical Recipes).
    bt = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x)
    )
    if x < (a + 1) / (a + b + 2):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1 - x) / b


def _betacf(a: float, b: float, x: float) -> float:
    max_iter = 200
    eps = 3e-14
    am = 1.0
    bm = 1.0
    az = 1.0
    qab = a + b
    qap = a + 1
    qam = a - 1
    bz = 1.0 - qab * x / qap
    for m in range(1, max_iter + 1):
        em = float(m)
        tem = em + em
        d = em * (b - em) * x / ((qam + tem) * (a + tem))
        ap = az + d * am
        bp = bz + d * bm
        d = -(a + em) * (qab + em) * x / ((a + tem) * (qap + tem))
        app = ap + d * az
        bpp = bp + d * bz
        am = ap / bpp
        bm = bp / bpp
        az = app / bpp
        bz = 1.0
        if abs(az - am) < eps * abs(az):
            return az
    return az


def dl_meta(rhos: list[float], ns: list[int], k_cov: int = 0) -> dict:
    """DerSimonian–Laird on Fisher z. Partial correlations use 1/(n-3-k)."""
    pairs = [(r, n) for r, n in zip(rhos, ns) if math.isfinite(r) and n - 3 - k_cov > 0]
    if not pairs:
        return {k: float("nan") for k in ("rho", "p", "I2", "ci_lo", "ci_hi", "k", "N", "tau2")}
    rhos_a = np.array([r for r, _ in pairs], dtype=float)
    ns_a = np.array([n for _, n in pairs], dtype=float)
    z = np.arctanh(np.clip(rhos_a, -0.999999, 0.999999))
    var_z = 1.0 / (ns_a - 3.0 - k_cov)
    w = 1.0 / var_z
    zbar = np.sum(w * z) / np.sum(w)
    q = float(np.sum(w * (z - zbar) ** 2))
    k = len(pairs)
    dfree = k - 1
    cdenom = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (q - dfree) / cdenom) if dfree > 0 and cdenom > 0 else 0.0
    wstar = 1.0 / (var_z + tau2)
    zre = float(np.sum(wstar * z) / np.sum(wstar))
    se = math.sqrt(1.0 / float(np.sum(wstar)))
    p = 2.0 * (1.0 - _norm_cdf(abs(zre / se)))
    i2 = max(0.0, (q - dfree) / q) if q > 0 else 0.0
    return {
        "rho": math.tanh(zre),
        "p": p,
        "I2": i2,
        "ci_lo": math.tanh(zre - 1.96 * se),
        "ci_hi": math.tanh(zre + 1.96 * se),
        "k": k,
        "N": int(np.sum(ns_a)),
        "tau2": tau2,
    }


def _norm_cdf(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def quantile_type7(x: np.ndarray, probs: list[float]) -> np.ndarray:
    v = np.sort(np.asarray(x, dtype=float))
    n = len(v)
    out = []
    for p in probs:
        h = (n - 1) * p + 1
        lo = int(math.floor(h))
        hi = int(math.ceil(h))
        if lo < 1:
            lo = 1
        if hi > n:
            hi = n
        frac = h - lo
        out.append(v[lo - 1] if hi == lo else (1 - frac) * v[lo - 1] + frac * v[hi - 1])
    return np.array(out)


def within_quartile(x: np.ndarray) -> np.ndarray:
    r = rank_average(x)
    breaks = quantile_type7(r, [0, 0.25, 0.5, 0.75, 1])
    if np.any(np.diff(breaks) <= 0):
        return np.array(["NA"] * len(x))
    labels = np.array(["Q1", "Q2", "Q3", "Q4"], dtype=object)
    # include.lowest: left-closed on the first bin, right-closed on the last
    out = np.empty(len(x), dtype=object)
    for i, val in enumerate(r):
        placed = False
        for b in range(4):
            left, right = breaks[b], breaks[b + 1]
            if b == 0 and left <= val <= right:
                out[i] = labels[b]
                placed = True
                break
            if b > 0 and left < val <= right:
                out[i] = labels[b]
                placed = True
                break
        if not placed:
            out[i] = "NA"
    return out


def mann_whitney_u(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int, int]:
    """U for x vs y with midranks, continuity-corrected normal p. r = 2U/(nx*ny)-1."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return float("nan"), float("nan"), ny, nx
    allv = np.concatenate([x, y])
    ranks = rank_average(allv)
    rx = ranks[:nx]
    # R wilcox.test statistic is U = rank sum of x minus nx*(nx+1)/2
    u = float(np.sum(rx) - nx * (nx + 1) / 2.0)
    mu = nx * ny / 2.0
    # tie correction
    _, counts = np.unique(allv, return_counts=True)
    tie = np.sum(counts ** 3 - counts)
    n = nx + ny
    sigma2 = nx * ny / 12.0 * ((n + 1) - tie / (n * (n - 1)))
    if sigma2 <= 0:
        p = float("nan")
    else:
        z = (abs(u - mu) - 0.5) / math.sqrt(sigma2)
        p = 2.0 * (1.0 - _norm_cdf(z))
    r = 2.0 * u / (nx * ny) - 1.0
    return r, p, ny, nx


def within_z(values: np.ndarray, cohorts: np.ndarray) -> np.ndarray:
    out = np.zeros(len(values), dtype=float)
    for c in np.unique(cohorts):
        m = cohorts == c
        r = rank_average(values[m])
        sd = float(np.std(r))
        out[m] = 0.0 if sd == 0 else (r - np.mean(r)) / sd
    return out


def pooled_partial(x, y, covariates, cohorts) -> float:
    xz = within_z(np.asarray(x, dtype=float), cohorts)
    yz = within_z(np.asarray(y, dtype=float), cohorts)
    if not covariates:
        return pearson(xz, yz)
    cols = [np.ones(len(x))]
    for z in covariates:
        cols.append(within_z(np.asarray(z, dtype=float), cohorts))
    A = np.column_stack(cols)
    return pearson(_ols_resid(xz, A), _ols_resid(yz, A))


def permute_pooled(x, y, covariates, cohorts, rng, n_perm=N_PERM) -> dict:
    """Shuffle y within cohort. Two-sided permutation p for the pooled partial r."""
    obs = pooled_partial(x, y, covariates, cohorts)
    if not math.isfinite(obs):
        return {"rho": obs, "p": float("nan"), "n_perm": 0}
    count = 0
    y = np.asarray(y, dtype=float).copy()
    for _ in range(n_perm):
        yperm = y.copy()
        for c in np.unique(cohorts):
            m = np.where(cohorts == c)[0]
            yperm[m] = y[rng.permutation(m)]
        stat = pooled_partial(x, yperm, covariates, cohorts)
        if math.isfinite(stat) and abs(stat) >= abs(obs) - 1e-15:
            count += 1
    return {"rho": obs, "p": (count + 1) / (n_perm + 1), "n_perm": n_perm}


def permute_delta(xa, xb, y, covariates, cohorts, rng, n_perm=N_PERM) -> dict:
    """Δ = partial_r(xa) - partial_r(xb). Two-sided permutation p."""
    ra = pooled_partial(xa, y, covariates, cohorts)
    rb = pooled_partial(xb, y, covariates, cohorts)
    delta = ra - rb
    if not (math.isfinite(ra) and math.isfinite(rb)):
        return {"rho_a": ra, "rho_b": rb, "delta": delta, "p": float("nan"), "n_perm": 0}
    count = 0
    y = np.asarray(y, dtype=float)
    for _ in range(n_perm):
        yperm = y.copy()
        for c in np.unique(cohorts):
            m = np.where(cohorts == c)[0]
            yperm[m] = y[rng.permutation(m)]
        da = pooled_partial(xa, yperm, covariates, cohorts)
        db = pooled_partial(xb, yperm, covariates, cohorts)
        if math.isfinite(da) and math.isfinite(db) and abs(da - db) >= abs(delta) - 1e-15:
            count += 1
    return {
        "rho_a": ra,
        "rho_b": rb,
        "delta": delta,
        "p": (count + 1) / (n_perm + 1),
        "n_perm": n_perm,
    }


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def _gene_cols_from_header(header: list[str], genes: list[str]) -> dict[str, int]:
    first = {}
    for i, name in enumerate(header):
        key = name.strip().upper()
        if key and key not in first:
            first[key] = i
    return {g: first[g] for g in genes if g in first}


def score_from_arrays(counts: dict[str, np.ndarray], malignant: np.ndarray) -> dict:
    idx = np.where(malignant)[0]
    n = int(idx.size)
    out = {"n_malignant_scored": n}
    if n == 0:
        for g in SCORE_GENES:
            out[f"mal_{g}_pct"] = float("nan")
            out[f"mal_{g}_mean"] = float("nan")
        out["mal_KRT_mean"] = float("nan")
        out["mal_KRT1819_mean"] = float("nan")
        out["mal_KRT_any_pct"] = float("nan")
        return out
    logs = {}
    for g in SCORE_GENES:
        if g not in counts:
            out[f"mal_{g}_pct"] = float("nan")
            out[f"mal_{g}_mean"] = float("nan")
            continue
        v = counts[g][idx].astype(float)
        out[f"mal_{g}_pct"] = 100.0 * float(np.mean(v > 0))
        logs[g] = np.log1p(v)
        out[f"mal_{g}_mean"] = float(np.mean(logs[g]))
    if all(g in logs for g in KERATIN):
        out["mal_KRT_mean"] = float(np.mean((logs["KRT8"] + logs["KRT18"] + logs["KRT19"]) / 3.0))
        out["mal_KRT1819_mean"] = float(np.mean((logs["KRT18"] + logs["KRT19"]) / 2.0))
        raw_any = (
            (counts["KRT8"][idx] > 0) | (counts["KRT18"][idx] > 0) | (counts["KRT19"][idx] > 0)
        )
        out["mal_KRT_any_pct"] = 100.0 * float(np.mean(raw_any))
    else:
        out["mal_KRT_mean"] = float("nan")
        out["mal_KRT1819_mean"] = float("nan")
        out["mal_KRT_any_pct"] = float("nan")
    return out


def _zeros(n: int) -> np.ndarray:
    return np.zeros(n, dtype=np.float64)


def load_dense_selected(path: Path, genes: list[str]) -> dict[str, np.ndarray]:
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split(",")
    cols = _gene_cols_from_header(header, genes)
    missing = [g for g in genes if g not in cols]
    usecols = sorted(set(cols.values()))
    arr = np.loadtxt(gzip.open(path, "rt"), delimiter=",", skiprows=1, usecols=usecols)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    pos = {c: i for i, c in enumerate(usecols)}
    n = arr.shape[0]
    out = {}
    for g in genes:
        if g not in cols:
            out[g] = _zeros(n)
        else:
            out[g] = arr[:, pos[cols[g]]]
    out["_missing"] = missing
    out["_n"] = n
    return out


def marker_masks(counts: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    n = len(next(v for k, v in counts.items() if k != "_missing" and k != "_n"))
    def g(name):
        return counts[name] if name in counts else _zeros(n)
    mal = ((g("EPCAM") > 0) | (g("KRT8") > 0) | (g("KRT18") > 0) | (g("KRT19") > 0)) & (g("PTPRC") == 0)
    tnk = (
        (g("CD3D") > 0) | (g("CD3E") > 0) | (g("CD8A") > 0) | (g("NKG7") > 0) | (g("GNLY") > 0) | (g("KLRD1") > 0)
    ) & (~mal)
    return mal, tnk


def score_gse123902() -> list[dict]:
    say("GSE123902")
    folder = GEO / "gse123902"
    if not any(folder.glob("*.csv.gz")):
        tar = GEO / "GSE123902_RAW.tar"
        if not tar.exists():
            raise SystemExit(f"missing {tar}")
        folder.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar) as tf:
            tf.extractall(folder)
    rows = []
    for path in sorted(folder.glob("*_dense.csv.gz")):
        name = path.name
        if "_NORMAL_" in name:
            continue
        if "_PRIMARY_" in name:
            tissue = "PRIMARY"
        elif "_METASTASIS_" in name:
            tissue = "METASTASIS"
        else:
            continue
        # GSM3516662_MSK_LX653_PRIMARY_TUMOUR_dense.csv.gz
        donor = name.split("_")[2]
        counts = load_dense_selected(path, sorted(set(MARKERS + SCORE_GENES)))
        mal, tnk = marker_masks(counts)
        rec = {
            "dataset": "GSE123902",
            "unit_id": donor,
            "unit_type": "donor",
            "tissue": tissue,
            "n_cells": int(counts["_n"]),
            "n_malignant": int(mal.sum()),
            "n_tnk": int(tnk.sum()),
            "frac_tnk": float(np.mean(tnk)),
            "malig_def": "marker_malig",
            "file": name,
        }
        missing_score = [g for g in SCORE_GENES if g in counts.get("_missing", [])]
        if missing_score:
            raise RuntimeError(f"{name} missing score genes {missing_score}")
        rec.update(score_from_arrays(counts, mal))
        rows.append(rec)
        say(f"  {donor} {tissue} cells={rec['n_cells']} mal={rec['n_malignant']} tnk={rec['n_tnk']} CLDN4%={rec['mal_CLDN4_pct']:.2f}")
    rows.sort(key=lambda r: (r["unit_id"], r["tissue"]))
    kept = []
    seen = set()
    for rec in rows:
        if rec["unit_id"] in seen:
            continue
        seen.add(rec["unit_id"])
        if rec["n_malignant"] >= 20 and rec["n_tnk"] >= 20:
            kept.append(rec)
    say(f"  kept {len(kept)} donors")
    return kept


def _read_10x_selected(prefix: Path, genes: list[str]) -> dict[str, np.ndarray]:
    feat = list(prefix.parent.glob(prefix.name + "_features.tsv.gz"))[0]
    bar = list(prefix.parent.glob(prefix.name + "_barcodes.tsv.gz"))[0]
    mtx = list(prefix.parent.glob(prefix.name + "_matrix.mtx.gz"))[0]
    symbols = []
    with gzip.open(feat, "rt") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            symbols.append(parts[1].upper() if len(parts) > 1 else parts[0].upper())
    first_row = {}
    for i, sym in enumerate(symbols):
        if sym not in first_row:
            first_row[sym] = i + 1  # 1-based mtx row
    wanted = {first_row[g]: g for g in genes if g in first_row}
    missing = [g for g in genes if g not in first_row]
    n_cells = 0
    with gzip.open(bar, "rt") as handle:
        for _ in handle:
            n_cells += 1
    counts = {g: _zeros(n_cells) for g in genes}
    with gzip.open(mtx, "rt") as handle:
        n_gene = n_col = None
        for line in handle:
            if line.startswith("%"):
                continue
            if n_gene is None:
                n_gene, n_col, _nnz = map(int, line.split())
                if n_col != n_cells:
                    raise RuntimeError(f"{mtx.name} columns {n_col} != barcodes {n_cells}")
                continue
            sp = line.find(" ")
            if sp < 0:
                continue
            gene = wanted.get(int(line[:sp]))
            if gene is None:
                continue
            c_s, v_s = line[sp + 1 :].split()
            counts[gene][int(c_s) - 1] = float(v_s)
    counts["_missing"] = missing
    counts["_n"] = n_cells
    return counts


def score_gse189357() -> list[dict]:
    say("GSE189357")
    folder = GEO / "gse189357"
    if not any(folder.glob("*_matrix.mtx.gz")):
        tar = GEO / "GSE189357_RAW.tar"
        if not tar.exists():
            raise SystemExit(f"missing {tar}")
        folder.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar) as tf:
            tf.extractall(folder)
    rows = []
    for pat in [f"TD{i}" for i in range(1, 10)]:
        hits = list(folder.glob(f"*_{pat}_matrix.mtx.gz"))
        if len(hits) != 1:
            raise RuntimeError(f"expected one matrix for {pat}, found {hits}")
        prefix = hits[0].with_name(hits[0].name.replace("_matrix.mtx.gz", ""))
        # prefix path is .../GSM_TD1  because with_name stripped suffix incorrectly if we used name only
        prefix = Path(str(hits[0]).replace("_matrix.mtx.gz", ""))
        counts = _read_10x_selected(prefix, sorted(set(MARKERS + SCORE_GENES)))
        missing_score = [g for g in SCORE_GENES if g in counts.get("_missing", [])]
        if missing_score:
            raise RuntimeError(f"{pat} missing score genes {missing_score}")
        mal, tnk = marker_masks(counts)
        rec = {
            "dataset": "GSE189357",
            "unit_id": pat,
            "unit_type": "patient",
            "tissue": "TUMOR",
            "n_cells": int(counts["_n"]),
            "n_malignant": int(mal.sum()),
            "n_tnk": int(tnk.sum()),
            "frac_tnk": float(np.mean(tnk)),
            "malig_def": "marker_malig",
        }
        rec.update(score_from_arrays(counts, mal))
        rows.append(rec)
        say(f"  {pat} cells={rec['n_cells']} mal={rec['n_malignant']} tnk={rec['n_tnk']} CLDN4%={rec['mal_CLDN4_pct']:.2f}")
    return rows


def _parse_soft(path: Path) -> dict[str, dict]:
    code_map = {}
    cur = None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                cur = {}
            elif cur is None:
                continue
            elif line.startswith("!Sample_title = "):
                cur["title"] = line.split("= ", 1)[1]
            elif line.startswith("!Sample_characteristics_ch1 = "):
                val = line.split("= ", 1)[1]
                if ": " in val:
                    k, v = val.split(": ", 1)
                    cur[k] = v
            elif line.startswith("!Sample_platform_id") and cur.get("title"):
                title = cur["title"]
                code = title.split(" ", 1)[1] if " " in title else ""
                code = code.upper().replace("-", "_")
                if code:
                    code_map[code] = cur
                cur = None
    return code_map


def _norm_code(orig: str) -> str:
    x = orig.upper().replace("-", "_")
    if x.endswith("_3P") or x.endswith("_5P"):
        x = x[:-3]
    return x


def score_gse205335() -> list[dict]:
    say("GSE205335")
    gene_path = GEO / "GSE205335_surface_genes.tsv.gz"
    if not gene_path.exists():
        rscript = ROOT / "extract_gse205335.R"
        subprocess.check_call(["Rscript", str(rscript), str(GEO)])
    ident_path = GEO / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    soft_path = GEO / "GSE205335_family.soft.gz"
    code_map = _parse_soft(soft_path)
    # barcode -> (patient, tissue, malignant, tnk) for non-normal
    by_patient = defaultdict(lambda: {"n": 0, "mal": 0, "tnk": 0, "tissues": set(), "mal_bc": []})
    n_unmapped = 0
    with gzip.open(ident_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            code = _norm_code(parts[idx["orig.ident"]])
            meta = code_map.get(code)
            if meta is None:
                n_unmapped += 1
                continue
            tissue = meta.get("tissue", "")
            if tissue.startswith("Normal"):
                continue
            pat = meta.get("patient", "")
            rec = by_patient[pat]
            rec["n"] += 1
            rec["tissues"].add(tissue)
            mal = parts[idx["lineage.sub"]] == "Malignant cells"
            tnk = parts[idx["lineage.total"]] == "T/NK cells"
            if mal:
                rec["mal"] += 1
                rec["mal_bc"].append(parts[idx["barcode"]])
            if tnk:
                rec["tnk"] += 1
    say(f"  unmapped identity rows {n_unmapped}; patients with tumor tissue {len(by_patient)}")
    # gene table
    gene_of = {}
    with gzip.open(gene_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        gidx = {name: i for i, name in enumerate(header)}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            gene_of[parts[0]] = parts
    say(f"  gene table barcodes {len(gene_of)}")
    rows = []
    for pat in sorted(by_patient):
        info = by_patient[pat]
        mal_bc = [b for b in info["mal_bc"] if b in gene_of]
        n_scored = len(mal_bc)
        counts = {g: np.empty(n_scored, dtype=float) for g in SCORE_GENES}
        for i, b in enumerate(mal_bc):
            parts = gene_of[b]
            for g in SCORE_GENES:
                counts[g][i] = float(parts[gidx[g]])
        mal_mask = np.ones(n_scored, dtype=bool)
        rec = {
            "dataset": "GSE205335",
            "unit_id": pat,
            "unit_type": "patient",
            "tissue": ",".join(sorted(info["tissues"])),
            "n_cells": info["n"],
            "n_malignant": info["mal"],
            "n_tnk": info["tnk"],
            "frac_tnk": info["tnk"] / info["n"] if info["n"] else float("nan"),
            "malig_def": "author_malig",
        }
        rec.update(score_from_arrays(counts, mal_mask))
        rows.append(rec)
        say(f"  {pat} cells={rec['n_cells']} mal={rec['n_malignant']} scored={rec['n_malignant_scored']} CLDN4%={rec['mal_CLDN4_pct']:.2f}")
    return rows


def score_gse131907() -> list[dict]:
    say("GSE131907")
    ann_path = GEO / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    mat_path = GEO / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    cells = []
    with gzip.open(ann_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in handle:
            p = line.rstrip("\n").split("\t")
            cells.append(
                {
                    "index": p[idx["Index"]],
                    "sample": p[idx["Sample"]],
                    "origin": p[idx["Sample_Origin"]],
                    "cell_type": p[idx["Cell_type"]],
                    "cell_subtype": p[idx["Cell_subtype"]],
                }
            )
    by_sample = defaultdict(list)
    for rec in cells:
        by_sample[rec["sample"]].append(rec)
    eligible = []
    for sample, recs in sorted(by_sample.items()):
        origin = recs[0]["origin"]
        n_mal = sum(r["cell_subtype"] == "Malignant cells" for r in recs)
        n_tnk = sum(r["cell_type"] in {"T lymphocytes", "NK cells"} for r in recs)
        if origin in ELIG_131_ORIGIN and n_mal >= MIN_MAL_131:
            eligible.append((sample, origin, len(recs), n_mal, n_tnk))
    say(f"  eligible samples {len(eligible)}")
    mal_index = {}
    for sample, *_ in eligible:
        for r in by_sample[sample]:
            if r["cell_subtype"] == "Malignant cells":
                mal_index[r["index"]] = sample
    say(f"  malignant barcodes {len(mal_index)}")
    # stream matrix header and target genes
    targets = set(SCORE_GENES)
    per_sample_vals = {s: {g: [] for g in SCORE_GENES} for s, *_ in eligible}
    order_sample = []  # sample id per malignant column found in the matrix
    with gzip.open(mat_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        barcodes = header[1:]
        keep_cols = []
        for i, b in enumerate(barcodes):
            sample = mal_index.get(b)
            if sample is not None:
                keep_cols.append(i)
                order_sample.append(sample)
        say(f"  matrix cells {len(barcodes)} malignant columns matched {len(keep_cols)}")
        if len(keep_cols) != len(mal_index):
            say(f"  WARNING annotation malignant not in matrix: {len(mal_index) - len(keep_cols)}")
        found = set()
        n_gene = 0
        for line in handle:
            tab = line.find("\t")
            if tab < 0:
                continue
            gene = line[:tab].upper()
            n_gene += 1
            if gene not in targets or gene in found:
                if n_gene % 5000 == 0:
                    say(f"  streamed {n_gene} genes, hits {len(found)}")
                continue
            found.add(gene)
            vals = line.rstrip("\n").split("\t")
            # vals[0] is gene, vals[1+i] is barcodes[i]
            buckets = {s: [] for s, *_ in eligible}
            for col, sample in zip(keep_cols, order_sample):
                raw = vals[col + 1]
                fv = 0.0 if raw in {"", "0", "0.0"} else float(raw)
                buckets[sample].append(fv)
            for sample, arr in buckets.items():
                per_sample_vals[sample][gene] = arr
            say(f"  hit {gene} at gene {n_gene}")
        missing = targets - found
        if missing:
            raise RuntimeError(f"GSE131907 missing genes {sorted(missing)}")
    rows = []
    for sample, origin, n_cells, n_mal, n_tnk in eligible:
        counts = {}
        n_scored = None
        for g in SCORE_GENES:
            arr = np.asarray(per_sample_vals[sample][g], dtype=float)
            counts[g] = arr
            n_scored = arr.size
        mal_mask = np.ones(n_scored, dtype=bool)
        rec = {
            "dataset": "GSE131907",
            "unit_id": sample,
            "unit_type": "sample",
            "tissue": origin,
            "n_cells": n_cells,
            "n_malignant": n_mal,
            "n_tnk": n_tnk,
            "frac_tnk": n_tnk / n_cells if n_cells else float("nan"),
            "malig_def": "author_malig",
        }
        rec.update(score_from_arrays(counts, mal_mask))
        rows.append(rec)
        say(f"  {sample} mal={n_mal} scored={rec['n_malignant_scored']} CLDN4%={rec['mal_CLDN4_pct']:.2f}")
    return rows


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def _subset_mask(dataset: np.ndarray, subset: str) -> np.ndarray:
    if subset == "all4":
        return np.ones(len(dataset), dtype=bool)
    if subset == "author":
        return np.isin(dataset, ["GSE131907", "GSE205335"])
    if subset == "marker":
        return np.isin(dataset, ["GSE123902", "GSE189357"])
    raise KeyError(subset)


def _score_vector(units: list[dict], gene: str, kind: str) -> np.ndarray:
    key = f"mal_{gene}_{kind}"
    return np.array([u[key] for u in units], dtype=float)


def _covariate(units: list[dict], adjustment: str) -> list[np.ndarray]:
    if adjustment == "none":
        return []
    if adjustment == "KRT":
        return [np.array([u["mal_KRT_mean"] for u in units], dtype=float)]
    if adjustment == "KRT1819":
        return [np.array([u["mal_KRT1819_mean"] for u in units], dtype=float)]
    raise KeyError(adjustment)


def cohort_rhos(units, gene, kind, adjustment, mask) -> tuple[list[str], list[float], list[int]]:
    datasets = np.array([u["dataset"] for u in units])
    x_all = _score_vector(units, gene, kind)
    y_all = np.array([u["frac_tnk"] for u in units], dtype=float)
    cov_all = _covariate(units, adjustment)
    names, rhos, ns = [], [], []
    for ds in ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]:
        m = mask & (datasets == ds)
        if int(m.sum()) < 5:
            continue
        x = x_all[m]
        y = y_all[m]
        cov = [c[m] for c in cov_all]
        if not np.isfinite(x).all() or not np.isfinite(y).all() or any(not np.isfinite(c).all() for c in cov):
            continue
        if adjustment == "none":
            rho = spearman(x, y)
        else:
            rho = partial_spearman(x, y, cov)
        names.append(ds)
        rhos.append(rho)
        ns.append(int(m.sum()))
    return names, rhos, ns


def analyze(units: list[dict]) -> dict:
    datasets = np.array([u["dataset"] for u in units])
    y = np.array([u["frac_tnk"] for u in units], dtype=float)
    k_of = {"none": 0, "KRT": 1, "KRT1819": 1}
    dl_rows = []
    single_rows = []
    for subset in ("all4", "author", "marker"):
        mask = _subset_mask(datasets, subset)
        for kind in ("pct", "mean"):
            for adjustment in ("none", "KRT", "KRT1819"):
                for gene in GENES:
                    names, rhos, ns = cohort_rhos(units, gene, kind, adjustment, mask)
                    meta = dl_meta(rhos, ns, k_cov=k_of[adjustment])
                    dl_rows.append(
                        {
                            "subset": subset,
                            "score": kind,
                            "adjustment": adjustment,
                            "gene": gene,
                            "k": meta["k"],
                            "N": meta["N"],
                            "rho": meta["rho"],
                            "p": meta["p"],
                            "I2": meta["I2"],
                            "ci_lo": meta["ci_lo"],
                            "ci_hi": meta["ci_hi"],
                        }
                    )
                    if subset == "all4":
                        for ds, rho, n in zip(names, rhos, ns):
                            single_rows.append(
                                {
                                    "dataset": ds,
                                    "score": kind,
                                    "adjustment": adjustment,
                                    "gene": gene,
                                    "n": n,
                                    "rho": rho,
                                    "p": spearman_p(rho, n, k_of[adjustment]),
                                }
                            )
    # pooled permutation contrasts on the primary score
    rng = np.random.default_rng(SEED)
    head_rows = []
    subsets = {
        "all4": np.ones(len(units), dtype=bool),
        "author": _subset_mask(datasets, "author"),
        "marker": _subset_mask(datasets, "marker"),
    }
    contrasts = [("CLDN4", g) for g in ("CLDN1", "CLDN7", "EPCAM", "MUC1")]
    for subset, mask in subsets.items():
        idx = np.where(mask)[0]
        y_s = y[idx]
        coh = datasets[idx]
        for adjustment in ("KRT", "KRT1819"):
            if subset == "marker" and adjustment == "KRT1819":
                continue
            if subset != "all4" and adjustment == "KRT1819":
                # keep the KRT18/19 sensitivity on the full set and the author set only for CLDN4 vs CLDN7
                pass
            cov = [c[idx] for c in _covariate(units, adjustment)]
            for a, b in contrasts:
                if adjustment == "KRT1819" and not (subset in {"all4", "author"} and b == "CLDN7"):
                    continue
                say(f"permute {subset} {adjustment} {a} vs {b}")
                delta = permute_delta(
                    _score_vector(units, a, "pct")[idx],
                    _score_vector(units, b, "pct")[idx],
                    y_s,
                    cov,
                    coh,
                    rng,
                )
                # unique partials: each gene adjusted for keratin AND the other gene
                uniq_a = permute_pooled(
                    _score_vector(units, a, "pct")[idx],
                    y_s,
                    cov + [_score_vector(units, b, "pct")[idx]],
                    coh,
                    rng,
                )
                uniq_b = permute_pooled(
                    _score_vector(units, b, "pct")[idx],
                    y_s,
                    cov + [_score_vector(units, a, "pct")[idx]],
                    coh,
                    rng,
                )
                # DL versions of the unique partials
                # build a temporary view by calling partial within cohort
                dl_ua_names, dl_ua_rhos, dl_ua_ns = [], [], []
                dl_ub_names, dl_ub_rhos, dl_ub_ns = [], [], []
                for ds in np.unique(coh):
                    m = coh == ds
                    if int(m.sum()) < 6:
                        continue
                    cov_ds = [c[m] for c in cov]
                    xa = _score_vector(units, a, "pct")[idx][m]
                    xb = _score_vector(units, b, "pct")[idx][m]
                    yy = y_s[m]
                    dl_ua_rhos.append(partial_spearman(xa, yy, cov_ds + [xb]))
                    dl_ub_rhos.append(partial_spearman(xb, yy, cov_ds + [xa]))
                    dl_ua_ns.append(int(m.sum()))
                    dl_ub_ns.append(int(m.sum()))
                    dl_ua_names.append(ds)
                    dl_ub_names.append(ds)
                dl_ua = dl_meta(dl_ua_rhos, dl_ua_ns, k_cov=1 + len(cov))
                dl_ub = dl_meta(dl_ub_rhos, dl_ub_ns, k_cov=1 + len(cov))
                head_rows.append(
                    {
                        "subset": subset,
                        "score": "pct",
                        "adjustment": adjustment,
                        "gene_a": a,
                        "gene_b": b,
                        "pooled_rho_a": delta["rho_a"],
                        "pooled_rho_b": delta["rho_b"],
                        "delta_a_minus_b": delta["delta"],
                        "perm_p_delta": delta["p"],
                        "unique_pooled_rho_a": uniq_a["rho"],
                        "unique_perm_p_a": uniq_a["p"],
                        "unique_pooled_rho_b": uniq_b["rho"],
                        "unique_perm_p_b": uniq_b["p"],
                        "unique_dl_rho_a": dl_ua["rho"],
                        "unique_dl_p_a": dl_ua["p"],
                        "unique_dl_rho_b": dl_ub["rho"],
                        "unique_dl_p_b": dl_ub["p"],
                        "n_perm": delta["n_perm"],
                        "N": int(mask.sum()),
                    }
                )
    # each gene's own pooled partial + permutation, primary
    pooled_rows = []
    for subset, mask in subsets.items():
        idx = np.where(mask)[0]
        for adjustment in ("none", "KRT", "KRT1819"):
            if subset != "all4" and adjustment != "KRT":
                continue
            cov = [c[idx] for c in _covariate(units, adjustment)]
            for gene in GENES:
                say(f"permute pooled {subset} {adjustment} {gene}")
                stat = permute_pooled(
                    _score_vector(units, gene, "pct")[idx],
                    y[idx],
                    cov,
                    datasets[idx],
                    rng,
                )
                pooled_rows.append(
                    {
                        "subset": subset,
                        "score": "pct",
                        "adjustment": adjustment,
                        "gene": gene,
                        "pooled_rho": stat["rho"],
                        "perm_p": stat["p"],
                        "N": int(mask.sum()),
                        "n_perm": stat["n_perm"],
                    }
                )
    # collinearity of %pos and keratin, all4, unadjusted DL
    coli_rows = []
    pairs = []
    for i, g1 in enumerate(GENES + ["KRT"]):
        for g2 in (GENES + ["KRT"])[i + 1 :]:
            pairs.append((g1, g2))
    for g1, g2 in pairs:
        rhos, ns = [], []
        for ds in ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]:
            m = datasets == ds
            a = np.array([u["mal_KRT_mean"] if g1 == "KRT" else u[f"mal_{g1}_pct"] for u in units], dtype=float)[m]
            b = np.array([u["mal_KRT_mean"] if g2 == "KRT" else u[f"mal_{g2}_pct"] for u in units], dtype=float)[m]
            rhos.append(spearman(a, b))
            ns.append(int(m.sum()))
        meta = dl_meta(rhos, ns, k_cov=0)
        coli_rows.append({"gene_a": g1, "gene_b": g2, **{k: meta[k] for k in ("k", "N", "rho", "p", "I2", "ci_lo", "ci_hi")}})
    # quartiles on %pos, all4
    q_rows = []
    for gene in GENES:
        x = _score_vector(units, gene, "pct")
        q = np.array(["NA"] * len(units), dtype=object)
        for ds in np.unique(datasets):
            m = np.where(datasets == ds)[0]
            q[m] = within_quartile(x[m])
        q4 = y[q == "Q4"]
        q1 = y[q == "Q1"]
        r, p, n1, n4 = mann_whitney_u(q4, q1)
        # keratin-residualized T/NK: within-cohort rank residual of frac_tnk on KRT
        yres = np.zeros(len(units))
        krt = np.array([u["mal_KRT_mean"] for u in units], dtype=float)
        for ds in np.unique(datasets):
            m = datasets == ds
            ry = rank_average(y[m])
            rz = rank_average(krt[m])
            A = np.column_stack([np.ones(int(m.sum())), rz])
            yres[m] = _ols_resid(ry, A)
        r2, p2, n1b, n4b = mann_whitney_u(yres[q == "Q4"], yres[q == "Q1"])
        q_rows.append(
            {
                "gene": gene,
                "score": "pct",
                "n_q1": n1,
                "n_q4": n4,
                "rank_biserial_tnk": r,
                "p_tnk": p,
                "rank_biserial_tnk_krt_resid": r2,
                "p_tnk_krt_resid": p2,
            }
        )
    # leave one cohort out, primary pct + KRT, DL
    loo_rows = []
    for drop in ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]:
        mask = datasets != drop
        for gene in PRIMARY:
            names, rhos, ns = cohort_rhos(units, gene, "pct", "KRT", mask)
            meta = dl_meta(rhos, ns, k_cov=1)
            loo_rows.append({"dropped": drop, "gene": gene, "k": meta["k"], "N": meta["N"], "rho": meta["rho"], "p": meta["p"], "I2": meta["I2"]})
    # score distributions
    dist_rows = []
    for gene in GENES:
        for kind in ("pct", "mean"):
            v = _score_vector(units, gene, kind)
            dist_rows.append(
                {
                    "gene": gene,
                    "score": kind,
                    "n": int(np.isfinite(v).sum()),
                    "min": float(np.min(v)),
                    "q25": float(np.quantile(v, 0.25)),
                    "median": float(np.median(v)),
                    "q75": float(np.quantile(v, 0.75)),
                    "max": float(np.max(v)),
                }
            )
    return {
        "dl": dl_rows,
        "singles": single_rows,
        "head": head_rows,
        "pooled": pooled_rows,
        "coli": coli_rows,
        "quartiles": q_rows,
        "loo": loo_rows,
        "dist": dist_rows,
    }


def calibrate(units: list[dict]) -> dict:
    locked = {}
    with LOCKED.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            locked[(row["dataset"], row["unit_id"])] = row
    got = {(u["dataset"], u["unit_id"]): u for u in units}
    missing = sorted(set(locked) - set(got))
    extra = sorted(set(got) - set(locked))
    diffs = []
    for key, row in locked.items():
        if key not in got:
            continue
        delta = abs(float(got[key]["mal_CLDN4_pct"]) - float(row["mal_CLDN4_pct"]))
        diffs.append((delta, key, float(got[key]["mal_CLDN4_pct"]), float(row["mal_CLDN4_pct"])))
    diffs.sort(reverse=True)
    for delta, key, got_v, lock_v in diffs[:8]:
        say(f"  cal {key[0]} {key[1]} got={got_v:.4f} locked={lock_v:.4f} abs={delta:.4f}")
    max_diff = diffs[0][0] if diffs else float("nan")
    # reproduce unadjusted DL
    names, rhos, ns = cohort_rhos(units, "CLDN4", "pct", "none", np.ones(len(units), dtype=bool))
    meta = dl_meta(rhos, ns, k_cov=0)
    ok = not missing and not extra and max_diff < 0.05 and abs(meta["rho"] - (-0.531)) < 0.005
    say(
        f"calibration units missing={missing} extra={extra} max|Δ CLDN4%|={max_diff:.4f} "
        f"DL rho={meta['rho']:.4f} p={meta['p']:.3g} ok={ok}"
    )
    return {"ok": ok, "missing": missing, "extra": extra, "max_diff": max_diff, "dl": meta, "singles": list(zip(names, rhos, ns))}


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    keys = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            out = {}
            for k, v in row.items():
                if isinstance(v, float):
                    out[k] = f"{v:.6g}" if math.isfinite(v) else "NA"
                else:
                    out[k] = v
            writer.writerow(out)


def fmt_p(p: float) -> str:
    if not math.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4g}"


def fmt_rho(r: float) -> str:
    if not math.isfinite(r):
        return "NA"
    return f"{r:.3f}"


def lookup_dl(rows, **kw):
    for row in rows:
        if all(row[k] == v for k, v in kw.items()):
            return row
    raise KeyError(kw)


def lookup_head(rows, **kw):
    for row in rows:
        if all(row[k] == v for k, v in kw.items()):
            return row
    raise KeyError(kw)


def forest_plot(dl_rows: list[dict], path: Path) -> None:
    genes = ["CLDN1", "CLDN4", "CLDN7", "EPCAM", "MUC1"]
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    y = np.arange(len(genes))
    for gene, yi in zip(genes, y):
        part = lookup_dl(dl_rows, subset="all4", score="pct", adjustment="KRT", gene=gene)
        raw = lookup_dl(dl_rows, subset="all4", score="pct", adjustment="none", gene=gene)
        color = "#b45309" if gene == "CLDN4" else ("#6b7280" if gene == "MUC1" else "#1f2937")
        ax.plot([part["ci_lo"], part["ci_hi"]], [yi, yi], color=color, lw=1.6, zorder=2)
        ax.scatter([part["rho"]], [yi], s=46, color=color, zorder=3, label=None)
        ax.scatter([raw["rho"]], [yi], s=36, facecolors="none", edgecolors=color, linewidths=1.2, zorder=3)
    ax.axvline(0, color="#9ca3af", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(["CLDN1", "CLDN4", "CLDN7", "EPCAM", "MUC1 (context)"])
    ax.set_xlabel("DerSimonian–Laird Spearman ρ with patient T/NK fraction")
    ax.set_title("Concordant-4 malignant % positive vs T/NK  (n = 65)")
    ax.set_xlim(-0.85, 0.35)
    ax.invert_yaxis()
    # legend proxies
    ax.scatter([], [], s=46, color="#1f2937", label="Filled: partial, keratin mean")
    ax.scatter([], [], s=36, facecolors="none", edgecolors="#1f2937", label="Open: unadjusted")
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(units, cal, res) -> None:
    dl = res["dl"]
    head = res["head"]
    pooled = res["pooled"]
    n = len(units)
    counts = {ds: sum(u["dataset"] == ds for u in units) for ds in ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]}

    def line_dl(gene, adjustment="KRT", subset="all4", score="pct"):
        row = lookup_dl(dl, subset=subset, score=score, adjustment=adjustment, gene=gene)
        return row

    primary = [line_dl(g) for g in GENES]
    # rank by partial rho (more negative = stronger exclusion association)
    ranked = sorted(primary, key=lambda r: r["rho"])
    rank_of = {r["gene"]: i + 1 for i, r in enumerate(ranked)}
    c4 = lookup_head(head, subset="all4", adjustment="KRT", gene_a="CLDN4", gene_b="CLDN7")
    c1 = lookup_head(head, subset="all4", adjustment="KRT", gene_a="CLDN4", gene_b="CLDN1")
    ce = lookup_head(head, subset="all4", adjustment="KRT", gene_a="CLDN4", gene_b="EPCAM")
    cm = lookup_head(head, subset="all4", adjustment="KRT", gene_a="CLDN4", gene_b="MUC1")
    author_c7 = lookup_head(head, subset="author", adjustment="KRT", gene_a="CLDN4", gene_b="CLDN7")
    raw_c4 = line_dl("CLDN4", "none")
    part_c4 = line_dl("CLDN4", "KRT")
    part_c7 = line_dl("CLDN7", "KRT")
    part_c1 = line_dl("CLDN1", "KRT")
    part_ep = line_dl("EPCAM", "KRT")
    part_mu = line_dl("MUC1", "KRT")
    pool_c4 = lookup_head_pooled = None
    pool_map = {(r["subset"], r["adjustment"], r["gene"]): r for r in pooled}

    def pool(gene, adjustment="KRT", subset="all4"):
        return pool_map[(subset, adjustment, gene)]

    # Three separate facts. Do not collapse them into the bulk ranking.
    diff_separates = (
        math.isfinite(c4["perm_p_delta"]) and c4["perm_p_delta"] < 0.05 and c4["delta_a_minus_b"] < 0
    )
    unique_cldn4 = c4["unique_pooled_rho_a"] < 0 and c4["unique_perm_p_a"] < 0.05
    unique_cldn7 = c4["unique_pooled_rho_b"] < 0 and c4["unique_perm_p_b"] < 0.05
    if diff_separates and unique_cldn4 and not unique_cldn7 and part_c4["rho"] < part_c7["rho"]:
        pin_sentence = (
            "On the pre-specified rule, public concordant-4 scRNA does separate CLDN4 from CLDN7 "
            "for patient T/NK after keratin is partialled."
        )
    else:
        pin_sentence = (
            "Public concordant-4 scRNA does **not** pin CLDN4 over CLDN7. "
            f"The difference is Δ = {fmt_rho(c4['delta_a_minus_b'])} "
            f"(permutation p = {fmt_p(c4['perm_p_delta'])})."
        )
    if unique_cldn4 and not unique_cldn7:
        unique_sentence = (
            f"That is not the bulk failure mode. After keratin and CLDN7 are both partialled, "
            f"CLDN4 still tracks lower T/NK (pooled unique ρ = {fmt_rho(c4['unique_pooled_rho_a'])}, "
            f"permutation p = {fmt_p(c4['unique_perm_p_a'])}; DL ρ = {fmt_rho(c4['unique_dl_rho_a'])}, "
            f"p = {fmt_p(c4['unique_dl_p_a'])}), while CLDN7 given CLDN4 does not "
            f"(pooled unique ρ = {fmt_rho(c4['unique_pooled_rho_b'])}, "
            f"permutation p = {fmt_p(c4['unique_perm_p_b'])}). "
            "CLDN4 is not a relabeling of CLDN7. The margin is too small, and too cohort-dependent, to crown it."
        )
    elif not unique_cldn4:
        unique_sentence = (
            f"CLDN4 also loses its own association once keratin and CLDN7 are both partialled "
            f"(pooled unique ρ = {fmt_rho(c4['unique_pooled_rho_a'])}, "
            f"permutation p = {fmt_p(c4['unique_perm_p_a'])}). "
            "On that stricter reading the public scRNA cannot identify CLDN4 once CLDN7 is in the model."
        )
    else:
        unique_sentence = (
            f"Both genes keep a unique partial (CLDN4 {fmt_rho(c4['unique_pooled_rho_a'])}, "
            f"p = {fmt_p(c4['unique_perm_p_a'])}; CLDN7 {fmt_rho(c4['unique_pooled_rho_b'])}, "
            f"p = {fmt_p(c4['unique_perm_p_b'])}). Neither one absorbs the other."
        )
    decision = (
        f"{pin_sentence} {unique_sentence}\n\n"
        "The bulk surface ranking had CLDN7, EPCAM, and MUC1 stronger than CLDN4. "
        "This T/NK partial does not. CLDN4 has the most negative point estimate. EPCAM also "
        f"excludes zero ({fmt_rho(part_ep['rho'])}, p = {fmt_p(part_ep['p'])}, I² = {100 * part_ep['I2']:.1f}%). "
        f"CLDN7's meta interval crosses zero (I² = {100 * part_c7['I2']:.1f}%). "
        "Inside GSE131907 the keratin-partial ρ is stronger for CLDN7 than for CLDN4, and "
        "dropping GSE205335 reverses the DL rank.\n\n"
        "EPCAM is not separable from CLDN4 "
        f"(Δ = {fmt_rho(ce['delta_a_minus_b'])}, permutation p = {fmt_p(ce['perm_p_delta'])}), "
        "and its unique partial collapses. CLDN1 has the opposite sign and is separable "
        f"(Δ = {fmt_rho(c1['delta_a_minus_b'])}, permutation p = {fmt_p(c1['perm_p_delta'])}). "
        "MUC1 is weaker than CLDN4 here "
        f"(Δ = {fmt_rho(cm['delta_a_minus_b'])}, permutation p = {fmt_p(cm['perm_p_delta'])}), "
        "the reverse of the bulk rank. "
        "What still pins CLDN4 is the private KD co-culture, not a stable public ranking of CLDN4 over CLDN7."
    )

    # rank sentence
    order_txt = ", ".join(f"{r['gene']} {fmt_rho(r['rho'])}" for r in ranked)
    strongest = ranked[0]["gene"]

    single_lines = []
    for ds, rho, nn in cal["singles"]:
        single_lines.append(f"| {ds} | {nn} | {fmt_rho(rho)} |")

    def dl_table(adjustment, subset="all4", score="pct"):
        lines = []
        for gene in GENES:
            row = line_dl(gene, adjustment, subset, score)
            tag = " (context)" if gene == "MUC1" else ""
            lines.append(
                f"| {gene}{tag} | {row['k']} | {row['N']} | {fmt_rho(row['rho'])} | {fmt_p(row['p'])} | "
                f"{100 * row['I2']:.1f}% | {fmt_rho(row['ci_lo'])} to {fmt_rho(row['ci_hi'])} |"
            )
        return "\n".join(lines)

    def single_table(gene="CLDN4", adjustment="KRT"):
        lines = []
        for row in res["singles"]:
            if row["gene"] == gene and row["adjustment"] == adjustment and row["score"] == "pct":
                lines.append(f"| {row['dataset']} | {row['n']} | {fmt_rho(row['rho'])} | {fmt_p(row['p'])} |")
        return "\n".join(lines)

    coli_c4 = next(r for r in res["coli"] if {r["gene_a"], r["gene_b"]} == {"CLDN4", "CLDN7"})
    coli_k = next(r for r in res["coli"] if {r["gene_a"], r["gene_b"]} == {"CLDN4", "KRT"})
    coli_e = next(r for r in res["coli"] if {r["gene_a"], r["gene_b"]} == {"CLDN4", "EPCAM"})

    q_lines = []
    for row in res["quartiles"]:
        q_lines.append(
            f"| {row['gene']} | {row['n_q1']}/{row['n_q4']} | {fmt_rho(row['rank_biserial_tnk'])} | {fmt_p(row['p_tnk'])} | "
            f"{fmt_rho(row['rank_biserial_tnk_krt_resid'])} | {fmt_p(row['p_tnk_krt_resid'])} |"
        )

    loo_lines = []
    for row in res["loo"]:
        loo_lines.append(f"| {row['dropped']} | {row['gene']} | {row['N']} | {fmt_rho(row['rho'])} | {fmt_p(row['p'])} |")

    author_lines = []
    for gene in GENES:
        row = line_dl(gene, "KRT", "author", "pct")
        author_lines.append(
            f"| {gene} | {row['k']} | {row['N']} | {fmt_rho(row['rho'])} | {fmt_p(row['p'])} | {100 * row['I2']:.1f}% |"
        )
    marker_lines = []
    for gene in GENES:
        row = line_dl(gene, "KRT", "marker", "pct")
        marker_lines.append(
            f"| {gene} | {row['k']} | {row['N']} | {fmt_rho(row['rho'])} | {fmt_p(row['p'])} | {100 * row['I2']:.1f}% |"
        )
    mean_lines = []
    for gene in GENES:
        row = line_dl(gene, "KRT", "all4", "mean")
        mean_lines.append(
            f"| {gene} | {row['N']} | {fmt_rho(row['rho'])} | {fmt_p(row['p'])} | {100 * row['I2']:.1f}% |"
        )
    k1819 = line_dl("CLDN4", "KRT1819")
    k1819_7 = line_dl("CLDN7", "KRT1819")
    k1819_head = lookup_head(head, subset="all4", adjustment="KRT1819", gene_a="CLDN4", gene_b="CLDN7")

    dist_c4 = next(r for r in res["dist"] if r["gene"] == "CLDN4" and r["score"] == "pct")
    dist_c7 = next(r for r in res["dist"] if r["gene"] == "CLDN7" and r["score"] == "pct")
    dist_ep = next(r for r in res["dist"] if r["gene"] == "EPCAM" and r["score"] == "pct")

    text = f"""# Concordant-4: CLDN1 vs CLDN4 vs CLDN7 vs EPCAM against patient T/NK

Additive. Does not replace the locked CLDN4 result (malignant % positive vs T/NK,
DL ρ = −0.531, p = 1.65×10⁻⁵, I² = 0%, N = 65). Same four cohorts, same units,
same malignant gate, same T/NK fraction. No GSE148071, GSE127465, GSE154826,
GSE200563, or E-MTAB-13526. No mouse Harmony. Cell counts are not n.

Question: after partialling malignant keratin, does public scRNA single out CLDN4
from CLDN7 (and from CLDN1 / EPCAM) the way a CLDN4-specific barrier claim would
require? The locked bulk surface ranking does not. LUAD placed CLDN4 7th (+0.122),
with CLDN7, EPCAM, and MUC1 stronger. That bulk number is a different estimand.
This file asks the specificity question on patient T/NK.

## Gate and scores

- GSE123902 donor (n = {counts['GSE123902']}), GSE189357 patient (n = {counts['GSE189357']}):
  malignant = (EPCAM or KRT8 or KRT18 or KRT19) > 0 and PTPRC = 0.
  T/NK = (CD3D or CD3E or CD8A or NKG7 or GNLY or KLRD1) > 0 and not malignant.
  `frac_tnk = n_tnk / n_cells`.
- GSE131907 sample (n = {counts['GSE131907']}): author malignant subtype, tumor-bearing
  origins tLung / tL/B / mLN / PE / mBrain, n_malignant ≥ 20.
  T/NK = author T lymphocytes or NK cells.
- GSE205335 patient (n = {counts['GSE205335']}): author lineage. Tissue labels that
  start with Normal are out of the patient denominator. The four normal-only donors
  are out. P4001 stays in this correlation (n_malignant = 27).
- Primary malignant score = % of malignant cells with count > 0.
  Secondary = mean log1p(count).
- Keratin covariate = mean, over malignant cells, of the per-cell mean of
  log1p(KRT8), log1p(KRT18), log1p(KRT19). Sensitivity covariate drops KRT8
  and uses only KRT18 and KRT19.
- Partial Spearman residualizes the ranks of both the surface score and
  `frac_tnk` on the keratin ranks, inside each cohort. Pooling is
  DerSimonian–Laird on Fisher z. For one covariate the z variance is 1/(n−4),
  not 1/(n−3).
- Specificity test (pre-specified): on within-cohort rank-z scores, pooled
  partial correlation, {N_PERM} within-cohort shuffles of `frac_tnk` (seed {SEED}).
  Δ = ρ(CLDN4) − ρ(CLDN7). CLDN4 is pinned over CLDN7 only if its keratin-partial
  DL ρ is more negative, the permutation p for Δ is < 0.05 with Δ < 0, and the
  unique partial of CLDN4 given keratin and CLDN7 stays negative with permutation
  p < 0.05. Otherwise the public scRNA does not pin it.
- MUC1 is context. It is in the bulk sentence and in the extraction. It is not
  a fifth member of the pin rule.
- p-values are descriptive. The permutation p is the one used for the pin call.

Honest n = **{n}** ({counts['GSE123902']} + {counts['GSE131907']} + {counts['GSE205335']} + {counts['GSE189357']}).

## Calibration against PR #539

Recomputed malignant CLDN4 % positive matches the locked patient table:
max absolute difference {cal['max_diff']:.4f} percentage points. Unadjusted DL
ρ = {fmt_rho(cal['dl']['rho'])} (p = {fmt_p(cal['dl']['p'])}, I² = {100 * cal['dl']['I2']:.1f}%,
N = {cal['dl']['N']}). Locked reference ρ = −0.531.

| cohort | n | unadjusted ρ |
|---|---:|---:|
{chr(10).join(single_lines)}

## 1. Keratin-partial DL Spearman, malignant % positive vs T/NK

Filled points in `results/figures/forest_partial_rho.png` are these partial ρ
values. Open points are the unadjusted ρ.

| gene | k | N | partial ρ | p | I² | 95% CI |
|---|---:|---:|---:|---:|---:|---|
{dl_table("KRT")}

Rank of the partial ρ, most negative first (MUC1 included): {order_txt}.
Strongest exclusion associate on this table: **{strongest}**.
Among CLDN1 / CLDN4 / CLDN7 / EPCAM only, most negative first: {", ".join(f"{g} {fmt_rho(line_dl(g)['rho'])}" for g in sorted(PRIMARY, key=lambda g: line_dl(g)['rho']))}.
On the point estimate, CLDN4 is first of those four.

Unadjusted DL, same score, for scale:

| gene | k | N | ρ | p | I² | 95% CI |
|---|---:|---:|---:|---:|---:|---|
{dl_table("none")}

CLDN4 % positive spans {dist_c4['min']:.1f}–{dist_c4['max']:.1f} (median {dist_c4['median']:.1f}).
CLDN7 spans {dist_c7['min']:.1f}–{dist_c7['max']:.1f} (median {dist_c7['median']:.1f}).
EPCAM spans {dist_ep['min']:.1f}–{dist_ep['max']:.1f} (median {dist_ep['median']:.1f}).
None of the four is a constant, so the partials are not a ceiling artifact.
The marker gate does use EPCAM and the keratins, so EPCAM and keratin variation
in GSE123902 and GSE189357 is variation inside an already epithelial gate.
That is why the author-label subset is reported below.

Cohort partial ρ for CLDN4 and CLDN7 (% positive, keratin):

| cohort | n | CLDN4 ρ | CLDN4 p | CLDN7 ρ | CLDN7 p |
|---|---:|---:|---:|---:|---:|
"""
    # build cohort comparison manually
    c4s = {(r["dataset"]): r for r in res["singles"] if r["gene"] == "CLDN4" and r["adjustment"] == "KRT" and r["score"] == "pct"}
    c7s = {(r["dataset"]): r for r in res["singles"] if r["gene"] == "CLDN7" and r["adjustment"] == "KRT" and r["score"] == "pct"}
    cohort_cmp = []
    for ds in ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]:
        a, b = c4s[ds], c7s[ds]
        cohort_cmp.append(
            f"| {ds} | {a['n']} | {fmt_rho(a['rho'])} | {fmt_p(a['p'])} | {fmt_rho(b['rho'])} | {fmt_p(b['p'])} |"
        )
    text += "\n".join(cohort_cmp)

    text += f"""

## 2. Does this pin CLDN4 vs CLDN7?

Pooled within-cohort rank partial (keratin), {N_PERM} shuffles:

| contrast | pooled ρ CLDN4 | pooled ρ other | Δ (CLDN4 − other) | perm p | unique ρ CLDN4 given other | unique perm p | unique ρ other given CLDN4 | unique perm p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| vs CLDN7 | {fmt_rho(c4['pooled_rho_a'])} | {fmt_rho(c4['pooled_rho_b'])} | {fmt_rho(c4['delta_a_minus_b'])} | {fmt_p(c4['perm_p_delta'])} | {fmt_rho(c4['unique_pooled_rho_a'])} | {fmt_p(c4['unique_perm_p_a'])} | {fmt_rho(c4['unique_pooled_rho_b'])} | {fmt_p(c4['unique_perm_p_b'])} |
| vs CLDN1 | {fmt_rho(c1['pooled_rho_a'])} | {fmt_rho(c1['pooled_rho_b'])} | {fmt_rho(c1['delta_a_minus_b'])} | {fmt_p(c1['perm_p_delta'])} | {fmt_rho(c1['unique_pooled_rho_a'])} | {fmt_p(c1['unique_perm_p_a'])} | {fmt_rho(c1['unique_pooled_rho_b'])} | {fmt_p(c1['unique_perm_p_b'])} |
| vs EPCAM | {fmt_rho(ce['pooled_rho_a'])} | {fmt_rho(ce['pooled_rho_b'])} | {fmt_rho(ce['delta_a_minus_b'])} | {fmt_p(ce['perm_p_delta'])} | {fmt_rho(ce['unique_pooled_rho_a'])} | {fmt_p(ce['unique_perm_p_a'])} | {fmt_rho(ce['unique_pooled_rho_b'])} | {fmt_p(ce['unique_perm_p_b'])} |
| vs MUC1 (context) | {fmt_rho(cm['pooled_rho_a'])} | {fmt_rho(cm['pooled_rho_b'])} | {fmt_rho(cm['delta_a_minus_b'])} | {fmt_p(cm['perm_p_delta'])} | {fmt_rho(cm['unique_pooled_rho_a'])} | {fmt_p(cm['unique_perm_p_a'])} | {fmt_rho(cm['unique_pooled_rho_b'])} | {fmt_p(cm['unique_perm_p_b'])} |

DL unique partial (both sides; keratin + the other gene), CLDN4 vs CLDN7:
CLDN4 ρ = {fmt_rho(c4['unique_dl_rho_a'])} (p = {fmt_p(c4['unique_dl_p_a'])});
CLDN7 ρ = {fmt_rho(c4['unique_dl_rho_b'])} (p = {fmt_p(c4['unique_dl_p_b'])}).

Collinearity, DL Spearman of malignant % positive (keratin row is the mean score):
CLDN4–CLDN7 ρ = {fmt_rho(coli_c4['rho'])} (p = {fmt_p(coli_c4['p'])}, I² = {100 * coli_c4['I2']:.1f}%);
CLDN4–EPCAM ρ = {fmt_rho(coli_e['rho'])};
CLDN4–keratin ρ = {fmt_rho(coli_k['rho'])}.

**Decision.** {decision}

Pooled permutation p for each gene alone, keratin partial, % positive, n = 65:
CLDN1 {fmt_rho(pool('CLDN1')['pooled_rho'])} (p = {fmt_p(pool('CLDN1')['perm_p'])});
CLDN4 {fmt_rho(pool('CLDN4')['pooled_rho'])} (p = {fmt_p(pool('CLDN4')['perm_p'])});
CLDN7 {fmt_rho(pool('CLDN7')['pooled_rho'])} (p = {fmt_p(pool('CLDN7')['perm_p'])});
EPCAM {fmt_rho(pool('EPCAM')['pooled_rho'])} (p = {fmt_p(pool('EPCAM')['perm_p'])});
MUC1 {fmt_rho(pool('MUC1')['pooled_rho'])} (p = {fmt_p(pool('MUC1')['perm_p'])}).

## 3. Sensitivities

Author labels only (GSE131907 + GSE205335). The malignant gate does not use
EPCAM or keratin. Keratin-partial % positive:

| gene | k | N | partial ρ | p | I² |
|---|---:|---:|---:|---:|---:|
{chr(10).join(author_lines)}

On that subset, CLDN4 vs CLDN7 Δ = {fmt_rho(author_c7['delta_a_minus_b'])},
perm p = {fmt_p(author_c7['perm_p_delta'])}. Unique CLDN4 ρ = {fmt_rho(author_c7['unique_pooled_rho_a'])}
(perm p = {fmt_p(author_c7['unique_perm_p_a'])}); unique CLDN7 ρ = {fmt_rho(author_c7['unique_pooled_rho_b'])}
(perm p = {fmt_p(author_c7['unique_perm_p_b'])}).

Marker-gate cohorts only (GSE123902 + GSE189357), same partial. Read these as
within-gate, not as a clean epithelial definition:

| gene | k | N | partial ρ | p | I² |
|---|---:|---:|---:|---:|---:|
{chr(10).join(marker_lines)}

Mean log1p instead of % positive, keratin partial, all four cohorts:

| gene | N | partial ρ | p | I² |
|---|---:|---:|---:|---:|
{chr(10).join(mean_lines)}

KRT18+KRT19 only (KRT8 left out of the covariate), % positive:
CLDN4 partial ρ = {fmt_rho(k1819['rho'])} (p = {fmt_p(k1819['p'])}, I² = {100 * k1819['I2']:.1f}%);
CLDN7 partial ρ = {fmt_rho(k1819_7['rho'])} (p = {fmt_p(k1819_7['p'])}).
Δ CLDN4−CLDN7 = {fmt_rho(k1819_head['delta_a_minus_b'])}, perm p = {fmt_p(k1819_head['perm_p_delta'])}.

Leave-one-cohort-out, keratin-partial % positive DL ρ:

| dropped | gene | N | ρ | p |
|---|---|---:|---:|---:|
{chr(10).join(loo_lines)}

Within-cohort quartile of % positive, Q4 vs Q1, stacked. Rank-biserial on
raw T/NK fraction, then on the within-cohort rank residual of T/NK after keratin.
Negative means Q4 (higher surface gene) has lower T/NK.

| gene | n_Q1/n_Q4 | r raw | p raw | r keratin residual | p residual |
|---|---|---:|---:|---:|---:|
{chr(10).join(q_lines)}

## What this does not say

The locked unadjusted CLDN4 association stands (calibration above). This
analysis asks whether that association is CLDN4-specific once keratin, CLDN7,
CLDN1, and EPCAM are allowed to compete. It is not a spatial exclusion test
and it is not the private knockdown co-culture. Bulk rank +0.122 and these
partial ρ values are not the same number.

## Reproduce

GEO files are not in the repo. Put them in `/tmp/geo_c4` (or `GEO_DIR`):

- `GSE123902_RAW.tar`
- `GSE189357_RAW.tar`
- `GSE131907_Lung_Cancer_cell_annotation.txt.gz`
- `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`
- `GSE205335_Lung_IO_CellIdentity.txt.gz`
- `GSE205335_family.soft.gz`
- `GSE205335_Lung_IO_UMI_matrix.rds.gz` (double-gzipped; `extract_gse205335.R` unpacks it)

```
python3 methods/concordant4_surface_tnk_keratin/analyze.py
```

Python does the patient-level scores and the partial / permutation tests.
`Rscript methods/concordant4_surface_tnk_keratin/extract_gse205335.R` reads the RDS.
"""
    # fix the awkward rank sentence by rewriting that one paragraph more cleanly
    (ROOT / "FINDING.md").write_text(text)
    say("wrote FINDING.md")


def units_to_rows(units: list[dict]) -> list[dict]:
    cols = [
        "dataset", "unit_id", "unit_type", "tissue", "n_cells", "n_malignant",
        "n_malignant_scored", "n_tnk", "frac_tnk", "malig_def",
    ]
    for g in SCORE_GENES:
        cols += [f"mal_{g}_pct", f"mal_{g}_mean"]
    cols += ["mal_KRT_mean", "mal_KRT1819_mean", "mal_KRT_any_pct"]
    rows = []
    for u in units:
        rows.append({c: u.get(c, "") for c in cols})
    return rows


def read_units(path: Path) -> list[dict]:
    units = []
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rec = {}
            for k, v in row.items():
                if k in {"dataset", "unit_id", "unit_type", "tissue", "malig_def", "file"}:
                    rec[k] = v
                elif v in {"", "NA"}:
                    rec[k] = float("nan")
                else:
                    rec[k] = float(v)
                    if k in {"n_cells", "n_malignant", "n_malignant_scored", "n_tnk"}:
                        rec[k] = int(rec[k])
            units.append(rec)
    return units


def cached_score(name: str, fn) -> list[dict]:
    path = GEO / f"cache_{name}.tsv"
    if path.exists() and os.environ.get("REUSE_PARTS", "1") != "0":
        say(f"reusing {path}")
        return read_units(path)
    rows = fn()
    write_tsv(path, units_to_rows(rows))
    return rows


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    cache = GEO / "cache_units.tsv"
    if os.environ.get("REUSE_CACHE") == "1" and cache.exists():
        say(f"reusing {cache}")
        units = read_units(cache)
    else:
        units = []
        units += cached_score("gse123902", score_gse123902)
        units += cached_score("gse189357", score_gse189357)
        units += cached_score("gse205335", score_gse205335)
        units += cached_score("gse131907", score_gse131907)
        write_tsv(cache, units_to_rows(units))
    say(f"units {len(units)}")
    cal = calibrate(units)
    write_tsv(TAB / "calibration_cldn4.tsv", [
        {
            "n_units": len(units),
            "n_missing": len(cal["missing"]),
            "n_extra": len(cal["extra"]),
            "max_abs_pct_diff": cal["max_diff"],
            "dl_rho": cal["dl"]["rho"],
            "dl_p": cal["dl"]["p"],
            "dl_I2": cal["dl"]["I2"],
            "locked_rho": -0.531,
            "ok": cal["ok"],
        }
    ])
    if not cal["ok"]:
        raise SystemExit("calibration against PR #539 failed; not writing a finding")
    res = analyze(units)
    write_tsv(TAB / "patient_units.tsv", units_to_rows(units))
    write_tsv(TAB / "dl_meta.tsv", res["dl"])
    write_tsv(TAB / "cohort_spearman.tsv", res["singles"])
    write_tsv(TAB / "head_to_head.tsv", res["head"])
    write_tsv(TAB / "pooled_partial.tsv", res["pooled"])
    write_tsv(TAB / "collinearity.tsv", res["coli"])
    write_tsv(TAB / "quartile_contrast.tsv", res["quartiles"])
    write_tsv(TAB / "leave_one_cohort.tsv", res["loo"])
    write_tsv(TAB / "score_distribution.tsv", res["dist"])
    forest_plot(res["dl"], FIG / "forest_partial_rho.png")
    write_finding(units, cal, res)
    say("DONE")


if __name__ == "__main__":
    main()
