#!/usr/bin/env python3
"""Concordant-4: which TJ gene best tracks immune-cold (patient T/NK).

Panel: CLDN4 vs CLDN3, CLDN7, OCLN, F11R, CDH1 on the locked n=65 units
(GSE123902 + GSE131907 + GSE205335 + GSE189357).

Primary estimand matches PR #644 / #503: malignant % positive vs same-patient
T/NK fraction, with optional keratin (KRT8/18/19 mean) partial Spearman.
Does not invent missing genes. Calibrates CLDN4 %pos to the locked ρ = −0.531.
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
LOCKED = ROOT / "data" / "pr644_patient_units.tsv"

# Pre-specified TJ gene panel for the PPT funnel.
TJ_PANEL = ["CLDN4", "CLDN3", "CLDN7", "OCLN", "F11R", "CDH1"]
KERATIN = ["KRT8", "KRT18", "KRT19"]
MARKERS = [
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "PTPRC",
    "CD3D",
    "CD3E",
    "CD8A",
    "NKG7",
    "GNLY",
    "KLRD1",
]
SCORE_GENES = TJ_PANEL + KERATIN
N_PERM = 10000
SEED = 1
ELIG_131_ORIGIN = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
MIN_MAL_131 = 20
LOCKED_CLDN4_RHO = -0.531168
LOCKED_CLDN4_P = 1.64622e-05


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
        return float("nan")
    t = rho * math.sqrt(df / (1.0 - rho * rho))
    # two-sided via regularized incomplete beta
    x = df / (df + t * t)
    a = 0.5 * df
    b = 0.5
    return float(_betainc_reg(a, b, x))


def _betainc_reg(a: float, b: float, x: float) -> float:
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    if x < (a + 1) / (a + b + 2):
        return math.exp(math.log(x) * a + math.log(1 - x) * b - lbeta) * _betacf(a, b, x) / a
    return 1.0 - math.exp(math.log(x) * a + math.log(1 - x) * b - lbeta) * _betacf(b, a, 1 - x) / b


def _betacf(a: float, b: float, x: float) -> float:
    max_it = 200
    eps = 3e-7
    am, bm, az = 1.0, 1.0, 1.0
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    bz = 1.0 - qab * x / qap
    for m in range(1, max_it + 1):
        em = float(m)
        tem = em + em
        d = em * (b - em) * x / ((qam + tem) * (a + tem))
        ap = az + d * am
        bp = bz + d * bm
        d = -(a + em) * (qab + em) * x / ((a + tem) * (qap + tem))
        app = ap + d * az
        bpp = bp + d * bz
        am, bm, az, bz = ap / bpp, bp / bpp, app / bpp, 1.0
        if abs(az - am) < eps * abs(az):
            return az
    return az


def dl_meta(rhos: list[float], ns: list[int], k_cov: int = 0) -> dict:
    zs, ws = [], []
    for r, n in zip(rhos, ns):
        if not math.isfinite(r) or n < 5 + k_cov:
            continue
        z = math.atanh(max(min(r, 0.999999), -0.999999))
        w = float(n - 3 - k_cov)
        if w <= 0:
            continue
        zs.append(z)
        ws.append(w)
    k = len(zs)
    if k == 0:
        return {"k": 0, "N": 0, "rho": float("nan"), "p": float("nan"), "I2": float("nan"),
                "ci_lo": float("nan"), "ci_hi": float("nan")}
    W = sum(ws)
    zbar = sum(w * z for w, z in zip(ws, zs)) / W
    Q = sum(w * (z - zbar) ** 2 for w, z in zip(ws, zs))
    df = k - 1
    C = W - sum(w * w for w in ws) / W if k > 1 else 0.0
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    ws_star = [1.0 / (1.0 / w + tau2) for w in ws]
    Wstar = sum(ws_star)
    zstar = sum(w * z for w, z in zip(ws_star, zs)) / Wstar
    se = math.sqrt(1.0 / Wstar)
    rho = math.tanh(zstar)
    # normal approx two-sided
    from math import erfc
    p = float(erfc(abs(zstar / se) / math.sqrt(2.0)))
    I2 = max(0.0, (Q - df) / Q) if Q > 0 and df > 0 else 0.0
    return {
        "k": k,
        "N": int(sum(ns)),
        "rho": rho,
        "p": p,
        "I2": I2,
        "ci_lo": math.tanh(zstar - 1.96 * se),
        "ci_hi": math.tanh(zstar + 1.96 * se),
    }


def _zeros(n: int) -> np.ndarray:
    return np.zeros(n, dtype=float)


def score_from_arrays(counts: dict[str, np.ndarray], malignant: np.ndarray) -> dict:
    out = {}
    n_mal = int(malignant.sum())
    out["n_malignant_scored"] = n_mal
    if n_mal == 0:
        for g in SCORE_GENES:
            out[f"mal_{g}_pct"] = float("nan")
            out[f"mal_{g}_mean"] = float("nan")
        out["mal_KRT_mean"] = float("nan")
        return out
    logs = {}
    for g in SCORE_GENES:
        if g not in counts:
            out[f"mal_{g}_pct"] = float("nan")
            out[f"mal_{g}_mean"] = float("nan")
            continue
        v = counts[g][malignant]
        out[f"mal_{g}_pct"] = 100.0 * float(np.mean(v > 0))
        logs[g] = np.log1p(v)
        out[f"mal_{g}_mean"] = float(np.mean(logs[g]))
    if all(g in logs for g in KERATIN):
        out["mal_KRT_mean"] = float(np.mean((logs["KRT8"] + logs["KRT18"] + logs["KRT19"]) / 3.0))
    else:
        out["mal_KRT_mean"] = float("nan")
    return out


def load_dense_selected(path: Path, genes: list[str]) -> dict[str, np.ndarray]:
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split(",")
        gene_cols = {}
        for i, name in enumerate(header):
            sym = name.strip().strip('"').upper()
            if sym in genes and sym not in gene_cols:
                gene_cols[sym] = i
        missing = [g for g in genes if g not in gene_cols]
        cols_needed = sorted(set(gene_cols.values()))
        rows = []
        for line in handle:
            parts = line.rstrip("\n").split(",")
            rows.append(parts)
        n = len(rows)
        counts = {g: _zeros(n) for g in genes}
        for g, col in gene_cols.items():
            arr = counts[g]
            for i, parts in enumerate(rows):
                raw = parts[col].strip().strip('"')
                if raw and raw not in {"0", "0.0"}:
                    arr[i] = float(raw)
        counts["_missing"] = missing
        counts["_n"] = n
        return counts


def marker_masks(counts: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    n = len(next(v for k, v in counts.items() if k not in {"_missing", "_n"}))

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
        }
        missing_score = [g for g in SCORE_GENES if g in counts.get("_missing", [])]
        if missing_score:
            raise RuntimeError(f"{name} missing score genes {missing_score}")
        rec.update(score_from_arrays(counts, mal))
        rows.append(rec)
        say(f"  {donor} {tissue} mal={rec['n_malignant']} tnk={rec['n_tnk']} CLDN4%={rec['mal_CLDN4_pct']:.2f}")
    rows.sort(key=lambda r: (r["unit_id"], r["tissue"]))
    kept, seen = [], set()
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
            first_row[sym] = i + 1
    wanted = {first_row[g]: g for g in genes if g in first_row}
    missing = [g for g in genes if g not in first_row]
    n_cells = sum(1 for _ in gzip.open(bar, "rt"))
    counts = {g: _zeros(n_cells) for g in genes}
    with gzip.open(mtx, "rt") as handle:
        n_gene = n_col = None
        for line in handle:
            if line.startswith("%"):
                continue
            if n_gene is None:
                n_gene, n_col, _nnz = map(int, line.split())
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
        folder.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar) as tf:
            tf.extractall(folder)
    rows = []
    for pat in [f"TD{i}" for i in range(1, 10)]:
        hits = list(folder.glob(f"*_{pat}_matrix.mtx.gz"))
        if len(hits) != 1:
            raise RuntimeError(f"expected one matrix for {pat}, found {hits}")
        prefix = Path(str(hits[0]).replace("_matrix.mtx.gz", ""))
        counts = _read_10x_selected(prefix, sorted(set(MARKERS + SCORE_GENES)))
        missing_score = [g for g in SCORE_GENES if g in counts.get("_missing", [])]
        if missing_score:
            raise RuntimeError(f"{pat} missing {missing_score}")
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
        say(f"  {pat} mal={rec['n_malignant']} CLDN4%={rec['mal_CLDN4_pct']:.2f}")
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
    gene_path = GEO / "GSE205335_tj_lead_genes.tsv.gz"
    if not gene_path.exists():
        rscript = ROOT / "extract_gse205335.R"
        subprocess.check_call(["Rscript", str(rscript), str(GEO)])
    ident_path = GEO / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    soft_path = GEO / "GSE205335_family.soft.gz"
    code_map = _parse_soft(soft_path)
    by_patient = defaultdict(lambda: {"n": 0, "mal": 0, "tnk": 0, "tissues": set(), "mal_bc": []})
    with gzip.open(ident_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            code = _norm_code(parts[idx["orig.ident"]])
            meta = code_map.get(code)
            if meta is None:
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
    gene_of = {}
    with gzip.open(gene_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        gidx = {name: i for i, name in enumerate(header)}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            gene_of[parts[0]] = parts
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
        say(f"  {pat} mal={rec['n_malignant']} CLDN4%={rec['mal_CLDN4_pct']:.2f}")
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
    targets = set(SCORE_GENES)
    per_sample_vals = {s: {g: [] for g in SCORE_GENES} for s, *_ in eligible}
    order_sample = []
    with gzip.open(mat_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        barcodes = header[1:]
        keep_cols = []
        for i, b in enumerate(barcodes):
            sample = mal_index.get(b)
            if sample is not None:
                keep_cols.append(i)
                order_sample.append(sample)
        say(f"  malignant columns matched {len(keep_cols)}")
        found = set()
        n_gene = 0
        for line in handle:
            tab = line.find("\t")
            if tab < 0:
                continue
            gene = line[:tab].upper()
            n_gene += 1
            if gene not in targets or gene in found:
                continue
            found.add(gene)
            vals = line.rstrip("\n").split("\t")
            buckets = {s: [] for s, *_ in eligible}
            for col, sample in zip(keep_cols, order_sample):
                raw = vals[col + 1]
                fv = 0.0 if raw in {"", "0", "0.0"} else float(raw)
                buckets[sample].append(fv)
            for sample, arr in buckets.items():
                per_sample_vals[sample][gene] = arr
            say(f"  hit {gene}")
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
        say(f"  {sample} mal={n_mal} CLDN4%={rec['mal_CLDN4_pct']:.2f}")
    return rows


def _score_vector(units, gene, kind):
    return np.array([u[f"mal_{gene}_{kind}"] for u in units], dtype=float)


def cohort_rhos(units, gene, kind, adjustment):
    datasets = np.array([u["dataset"] for u in units])
    x_all = _score_vector(units, gene, kind)
    y_all = np.array([u["frac_tnk"] for u in units], dtype=float)
    if adjustment == "KRT":
        cov_all = [np.array([u["mal_KRT_mean"] for u in units], dtype=float)]
        k_cov = 1
    else:
        cov_all = []
        k_cov = 0
    names, rhos, ns = [], [], []
    for ds in ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]:
        m = datasets == ds
        if m.sum() < 5:
            continue
        x, y = x_all[m], y_all[m]
        cov = [c[m] for c in cov_all]
        rho = partial_spearman(x, y, cov) if cov else spearman(x, y)
        names.append(ds)
        rhos.append(rho)
        ns.append(int(m.sum()))
    return names, rhos, ns, k_cov


def permute_delta(xa, xb, y, covariates, cohorts, rng, n_perm=N_PERM):
    """Label-swap test of ρ(a,y|cov) − ρ(b,y|cov) within cohorts."""
    obs_a = pooled_partial(xa, y, covariates, cohorts)
    obs_b = pooled_partial(xb, y, covariates, cohorts)
    obs = obs_a - obs_b
    n_ge = 0
    for _ in range(n_perm):
        swap = rng.random(len(xa)) < 0.5
        a2 = np.where(swap, xb, xa)
        b2 = np.where(swap, xa, xb)
        d = pooled_partial(a2, y, covariates, cohorts) - pooled_partial(b2, y, covariates, cohorts)
        if abs(d) >= abs(obs):
            n_ge += 1
    return {"delta": obs, "rho_a": obs_a, "rho_b": obs_b, "p": (n_ge + 1) / (n_perm + 1)}


def pooled_partial(x, y, covariates, cohorts):
    """Within-cohort z-scored partial Spearman stacked, then Spearman."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    uniq = sorted(set(cohorts))
    xr, yr = [], []
    for c in uniq:
        m = cohorts == c
        if m.sum() < 5:
            continue
        cov = [z[m] for z in covariates]
        # residualize ranks then z-score within cohort
        rx = rank_average(x[m])
        ry = rank_average(y[m])
        if cov:
            A = np.column_stack([np.ones(m.sum())] + [rank_average(z) for z in cov])
            rx = _ols_resid(rx, A)
            ry = _ols_resid(ry, A)
        if np.std(rx) == 0 or np.std(ry) == 0:
            continue
        xr.append((rx - rx.mean()) / rx.std())
        yr.append((ry - ry.mean()) / ry.std())
    if not xr:
        return float("nan")
    return pearson(np.concatenate(xr), np.concatenate(yr))


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str] | None = None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore", delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def calibrate(units):
    locked = {}
    with LOCKED.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            locked[(row["dataset"], row["unit_id"])] = row
    diffs = []
    missing = 0
    for u in units:
        key = (u["dataset"], u["unit_id"])
        if key not in locked:
            missing += 1
            continue
        diffs.append(abs(float(u["mal_CLDN4_pct"]) - float(locked[key]["mal_CLDN4_pct"])))
    names, rhos, ns, _ = cohort_rhos(units, "CLDN4", "pct", "none")
    meta = dl_meta(rhos, ns, 0)
    return {
        "n_units": len(units),
        "locked_units": len(locked),
        "missing_vs_locked": missing,
        "max_abs_delta_cldn4_pct": max(diffs) if diffs else float("nan"),
        "meta_rho_cldn4_pct": meta["rho"],
        "meta_p_cldn4_pct": meta["p"],
        "locked_rho": LOCKED_CLDN4_RHO,
        "locked_p": LOCKED_CLDN4_P,
        "rho_delta_vs_locked": abs(meta["rho"] - LOCKED_CLDN4_RHO) if math.isfinite(meta["rho"]) else float("nan"),
    }


def main():
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    units = []
    units.extend(score_gse123902())
    units.extend(score_gse131907())
    units.extend(score_gse205335())
    units.extend(score_gse189357())
    say(f"total units {len(units)}")

    # write patient units
    field_order = [
        "dataset",
        "unit_id",
        "unit_type",
        "tissue",
        "n_cells",
        "n_malignant",
        "n_malignant_scored",
        "n_tnk",
        "frac_tnk",
        "malig_def",
        "mal_KRT_mean",
    ]
    for g in TJ_PANEL:
        field_order += [f"mal_{g}_pct", f"mal_{g}_mean"]
    write_tsv(TAB / "patient_units.tsv", units, field_order)

    cal = calibrate(units)
    write_tsv(TAB / "calibration_vs_pr644.tsv", [cal])
    say(f"calibration max|ΔCLDN4%|={cal['max_abs_delta_cldn4_pct']} rhoΔ={cal['rho_delta_vs_locked']}")

    meta_rows = []
    cohort_rows = []
    for kind in ("pct", "mean"):
        for adj in ("none", "KRT"):
            for gene in TJ_PANEL:
                names, rhos, ns, k_cov = cohort_rhos(units, gene, kind, adj)
                meta = dl_meta(rhos, ns, k_cov)
                meta_rows.append(
                    {
                        "score": kind,
                        "adjustment": adj,
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
                for name, rho, n in zip(names, rhos, ns):
                    cohort_rows.append(
                        {
                            "score": kind,
                            "adjustment": adj,
                            "gene": gene,
                            "dataset": name,
                            "n": n,
                            "rho": rho,
                            "p": spearman_p(rho, n, k_cov),
                        }
                    )
    write_tsv(TAB / "dl_meta.tsv", meta_rows)
    write_tsv(TAB / "cohort_spearman.tsv", cohort_rows)

    # Head-to-head: CLDN4 vs each other gene (pct, KRT)
    y = np.array([u["frac_tnk"] for u in units], dtype=float)
    krt = np.array([u["mal_KRT_mean"] for u in units], dtype=float)
    cohorts = np.array([u["dataset"] for u in units])
    xa = _score_vector(units, "CLDN4", "pct")
    rng = np.random.default_rng(SEED)
    h2h = []
    for gene in TJ_PANEL:
        if gene == "CLDN4":
            continue
        xb = _score_vector(units, gene, "pct")
        res = permute_delta(xa, xb, y, [krt], cohorts, rng)
        h2h.append(
            {
                "gene_a": "CLDN4",
                "gene_b": gene,
                "adjustment": "KRT",
                "score": "pct",
                "rho_a": res["rho_a"],
                "rho_b": res["rho_b"],
                "delta": res["delta"],
                "perm_p": res["p"],
                "n_perm": N_PERM,
            }
        )
        say(f"H2H CLDN4 vs {gene}: Δ={res['delta']:.3f} p={res['p']:.4f}")
    write_tsv(TAB / "head_to_head.tsv", h2h)

    # Ranking table for PPT
    primary = [r for r in meta_rows if r["score"] == "pct" and r["adjustment"] == "KRT"]
    primary_sorted = sorted(primary, key=lambda r: r["rho"] if math.isfinite(r["rho"]) else 999)
    rank_rows = []
    for i, r in enumerate(primary_sorted, 1):
        rank_rows.append({**r, "rank_most_negative": i, "is_lead_point_estimate": i == 1})
    write_tsv(TAB / "rank_keratin_partial_pct.tsv", rank_rows)

    # Forest plot
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    genes_ord = [r["gene"] for r in primary_sorted][::-1]
    ypos = np.arange(len(genes_ord))
    for i, gene in enumerate(genes_ord):
        r = next(x for x in primary if x["gene"] == gene)
        color = "#b45309" if gene == "CLDN4" else "#1f2937"
        ax.errorbar(
            r["rho"],
            i,
            xerr=[[r["rho"] - r["ci_lo"]], [r["ci_hi"] - r["rho"]]],
            fmt="o",
            color=color,
            capsize=3,
        )
        ax.text(0.98, i, f"ρ={r['rho']:.3f}  p={r['p']:.2g}  I²={100*r['I2']:.0f}%",
                transform=ax.get_yaxis_transform(), va="center", ha="right", fontsize=8, color=color)
    ax.axvline(0, color="#9ca3af", lw=1)
    ax.set_yticks(ypos)
    ax.set_yticklabels(genes_ord)
    ax.set_xlabel("Keratin-partial Spearman ρ (malignant %pos vs T/NK)")
    ax.set_title("Concordant-4 TJ panel vs immune-cold (n=65)")
    fig.tight_layout()
    fig.savefig(FIG / "c4_tj_forest.png", dpi=150)
    fig.savefig(FIG / "c4_tj_forest.pdf")
    plt.close(fig)
    say("wrote concordant-4 tables + forest")


if __name__ == "__main__":
    main()
