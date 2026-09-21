#!/usr/bin/env python3
"""Partial dependence of barrier-ligand outgoing on TACSTD2 vs CLDN4.

Concordant-4 malignant senders to T/NK. Patient is the unit.
Barrier ligands: F11R, NECTIN2, CDH1, LGALS9.

Crude contrast: Q4 vs Q1 of the gate gene inside malignant cells.
Does the TACSTD2 high−low gap shrink after the ligand is residualized on
CLDN4, or after Q4/Q1 cells are matched inside CLDN4 strata?

Co-primary scores, fixed before any cohort is ranked:
  expr_prop_pp, cellchat_hill, cpdb_means, liana_logfc.

A shrink is a drop of at least 20% of a positive crude mean, with a
two-sided paired Wilcoxon p < 0.05 on (crude − adjusted). The rule is
applied to every co-primary score. It is not used to pick a score.
"""

from __future__ import annotations

import math
import os
import tarfile
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import stats

HERE = Path(__file__).resolve().parents[1]
RAW = Path(os.environ.get("CONCORDANT4_RAW", "/tmp/concordant4_raw"))
CACHE = Path(os.environ.get("CONCORDANT4_CACHE", "/tmp/concordant4_cache"))
TAB = HERE / "results" / "tables"
FIG = HERE / "results" / "figures"

TRIM = 0.10
EXPR_PROP = 0.10
KH = 0.5
MIN_ARM = 10
MIN_TNK = 20
MIN_MAL = 40
MIN_PAIRS = 10
MIN_LOCAL = 5
MIN_LOCAL_STRATA = 2
SEED = 3979
N_FLIP = 10000
SHRINK_MIN = 0.20
LN2 = math.log(2.0)
N_STRATA = 4

ALIASES = {"PVRL2": "NECTIN2", "JAM1": "F11R"}
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19"]
TNK_MARKERS = ["CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]
COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
COHORT_COLOR = {
    "GSE123902": "#4C78A8",
    "GSE131907": "#F58518",
    "GSE205335": "#54A24B",
    "GSE189357": "#E45756",
}
EDGES = [
    ("F11R", "ITGAL_ITGB2"),
    ("F11R", "F11R"),
    ("NECTIN2", "TIGIT"),
    ("NECTIN2", "CD96"),
    ("CDH1", "ITGAE_ITGB7"),
    ("CDH1", "KLRG1"),
    ("LGALS9", "HAVCR2"),
    ("LGALS9", "CD44"),
    ("LGALS9", "PTPRC"),
]
LIGANDS = ["F11R", "NECTIN2", "CDH1", "LGALS9"]
SCORES = ["expr_prop_pp", "cellchat_hill", "cpdb_means", "liana_logfc"]
# Shrink is read off these two. The other two scores are reported with the same rule.
COPRIMARY = ["expr_prop_pp", "cellchat_hill"]
CONTRASTS = ["crude", "strata_resid", "linear_resid", "matched", "strata_local"]


def log(msg: str) -> None:
    print(msg, flush=True)


def subunits(token: str) -> list[str]:
    return token.split("_")


def wanted_genes() -> set[str]:
    genes = set(EPI + TNK_MARKERS + ["CLDN4", "TACSTD2", "PTPRC"])
    for lig, rec in EDGES:
        genes.add(lig)
        genes.update(subunits(rec))
    genes.update(ALIASES)
    return genes


def collapse_named(named: dict[str, np.ndarray], keep: set[str]) -> dict[str, np.ndarray]:
    buckets: dict[str, list[tuple[str, np.ndarray]]] = {}
    for gene, arr in named.items():
        g = gene.upper()
        canon = ALIASES.get(g, g)
        if canon not in keep and g not in keep:
            continue
        buckets.setdefault(canon, []).append((g, np.asarray(arr)))
    out = {}
    for canon, items in buckets.items():
        if canon not in keep:
            continue
        official = [a for g, a in items if g == canon]
        use = official if official else [a for _, a in items]
        out[canon] = use[0] if len(use) == 1 else np.sum(np.vstack(use), axis=0)
    return out


def matrix_from_named(named: dict[str, np.ndarray], n: int) -> tuple[np.ndarray, list[str]]:
    genes = sorted(named)
    if not genes:
        return np.zeros((n, 0), dtype=np.float32), []
    mat = np.column_stack([np.asarray(named[g], dtype=np.float32).reshape(-1) for g in genes])
    if mat.shape[0] != n:
        raise RuntimeError(f"rows {mat.shape[0]} != cells {n}")
    return mat, genes


def pos_any(mat: np.ndarray, genes: list[str], markers: list[str]) -> np.ndarray:
    idx = [genes.index(g) for g in markers if g in genes]
    if not idx:
        return np.zeros(mat.shape[0], dtype=bool)
    return np.any(mat[:, idx] > 0, axis=1)


def marker_compartments(mat: np.ndarray, genes: list[str]) -> np.ndarray:
    mal = pos_any(mat, genes, EPI) & ~pos_any(mat, genes, ["PTPRC"])
    tnk = pos_any(mat, genes, TNK_MARKERS) & ~mal
    out = np.full(mat.shape[0], "", dtype=object)
    out[tnk] = "TNK"
    out[mal] = "MAL"
    return out


def pack_unit(patient: str, mat: np.ndarray, genes: list[str], lib: np.ndarray, comp: np.ndarray) -> dict | None:
    keep = comp != ""
    if not np.any(keep):
        return None
    return {
        "patient": str(patient),
        "mat": np.asarray(mat[keep], dtype=np.float32),
        "genes": genes,
        "lib": np.asarray(lib[keep], dtype=np.float64),
        "comp": comp[keep],
    }


def read_dense_csv(path: Path) -> tuple[np.ndarray, list[str], np.ndarray]:
    df = pd.read_csv(path)
    genes = [str(c).upper() for c in df.columns[1:]]
    mat = df.iloc[:, 1:].to_numpy(dtype=np.float32)
    lib = mat.sum(axis=1).astype(np.float64)
    return mat, genes, lib


def load_gse123902(keep: set[str]) -> list[dict]:
    log("==== GSE123902 ====")
    units = pd.read_csv(HERE / "data" / "GSE123902_marker_units.tsv", sep="\t")
    tumor = units[units["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor["ord"] = np.where(tumor["tissue"].eq("PRIMARY"), 0, 1)
    tumor = tumor.sort_values(["patient", "ord"]).drop_duplicates("patient")
    csv_dir = RAW / "GSE123902" / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)
    if not list(csv_dir.glob("*_dense.csv.gz")):
        tar = RAW / "GSE123902" / "GSE123902_RAW.tar"
        log(f"  untar {tar.name}")
        with tarfile.open(tar) as tf:
            tf.extractall(csv_dir)
    out = []
    for row in tumor.itertuples(index=False):
        fp = csv_dir / row.file
        if not fp.exists():
            hits = sorted(p for p in csv_dir.glob(f"*{row.patient}*dense.csv.gz") if "NORMAL" not in p.name)
            if not hits:
                log(f"  missing {row.patient}")
                continue
            fp = hits[0]
        log(f"  {row.patient}")
        mat, genes, lib = read_dense_csv(fp)
        named = collapse_named({g: mat[:, i] for i, g in enumerate(genes)}, keep)
        mat, genes = matrix_from_named(named, mat.shape[0])
        comp = marker_compartments(mat, genes)
        packed = pack_unit(str(row.patient), mat, genes, lib, comp)
        if packed:
            packed["cohort"] = "GSE123902"
            out.append(packed)
    return out


def load_gse131907(keep: set[str]) -> list[dict]:
    log("==== GSE131907 ====")
    cache = CACHE / "gse131907"
    lib_path = cache / "libsize.f64"
    if not lib_path.exists() or not (cache / "genes_found.txt").exists():
        raise RuntimeError("GSE131907 slim matrix is missing. Run prepare_131907.py and extract_131907.")
    meta = pd.read_csv(cache / "cells.tsv", sep="\t")
    lib = np.fromfile(lib_path, dtype=np.float64)
    if lib.size != len(meta):
        raise RuntimeError(f"lib {lib.size} vs cells {len(meta)}")
    found = [ln.strip() for ln in (cache / "genes_found.txt").read_text().splitlines() if ln.strip()]
    named = {g: np.fromfile(cache / f"{g}.f32", dtype=np.float32) for g in found}
    named = collapse_named(named, keep)
    mat, genes = matrix_from_named(named, len(meta))
    comp = meta["comp"].astype(str).to_numpy()
    sample = meta["sample"].astype(str).to_numpy()
    out = []
    for s in sorted(set(sample)):
        which = np.where(sample == s)[0]
        packed = pack_unit(s, mat[which], genes, lib[which], comp[which])
        if packed:
            packed["cohort"] = "GSE131907"
            out.append(packed)
    log(f"  units {len(out)} cells {len(meta)}")
    return out


def load_gse205335(keep: set[str]) -> list[dict]:
    log("==== GSE205335 ====")
    export = CACHE / "gse205335"
    mtx = export / "matrix.mtx"
    if not mtx.exists():
        raise RuntimeError("GSE205335 matrix.mtx is missing. Run export_gse205335.R.")
    genes = [ln.strip() for ln in (export / "genes.txt").read_text().splitlines() if ln.strip()]
    cells = [ln.strip() for ln in (export / "cells.txt").read_text().splitlines() if ln.strip()]
    log(f"  MTX {len(genes)} x {len(cells)}")
    mm = spio.mmread(mtx).tocsr()
    if mm.shape == (len(genes), len(cells)):
        mm = mm.T.tocsr()
    elif mm.shape != (len(cells), len(genes)):
        raise RuntimeError(f"MTX shape {mm.shape}")
    lib_df = pd.read_csv(export / "libsize.tsv", sep="\t")
    lib_map = dict(zip(lib_df["cell"].astype(str), lib_df["libsize"].astype(float)))
    lib = np.array([lib_map[c] for c in cells], dtype=np.float64)
    named = collapse_named({g: np.asarray(mm[:, i].todense()).ravel() for i, g in enumerate(genes)}, keep)
    del mm
    mat, genes = matrix_from_named(named, len(cells))
    ident = pd.read_csv(RAW / "GSE205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    gsm = pd.read_csv(HERE / "data" / "GSE205335_gsm_map.tsv", sep="\t")
    ident = ident.merge(gsm[["orig.ident", "patient", "tissue"]], on="orig.ident", how="left")
    ident["barcode"] = ident["barcode"].astype(str)
    idmap = ident.drop_duplicates("barcode").set_index("barcode")
    common = [i for i, c in enumerate(cells) if c in idmap.index]
    if len(common) < 0.9 * len(cells):
        raise RuntimeError(f"identity overlap {len(common)}/{len(cells)}")
    idx = np.array(common)
    sub = idmap.loc[[cells[i] for i in idx]]
    normal = sub["tissue"].fillna("").astype(str).str.startswith("Normal").to_numpy()
    lineage_sub = sub["lineage.sub"].astype(str).to_numpy()
    lineage = sub["lineage.total"].astype(str).to_numpy()
    patient = sub["patient"].astype(str).to_numpy()
    mal = (lineage_sub == "Malignant cells") & ~normal
    tnk = (lineage == "T/NK cells") & ~normal & ~mal
    comp = np.full(idx.size, "", dtype=object)
    comp[tnk] = "TNK"
    comp[mal] = "MAL"
    mat, lib = mat[idx], lib[idx]
    locked = pd.read_csv(HERE / "data" / "GSE205335_patients.tsv", sep="\t")
    keep_pt = set(locked.loc[locked["n_malignant"] > 0, "patient"].astype(str))
    out = []
    for pt in sorted(keep_pt):
        which = np.where(patient == pt)[0]
        if which.size == 0:
            continue
        packed = pack_unit(pt, mat[which], genes, lib[which], comp[which])
        if packed:
            packed["cohort"] = "GSE205335"
            out.append(packed)
    log(f"  units {len(out)}")
    return out


def load_gse189357(keep: set[str]) -> list[dict]:
    log("==== GSE189357 ====")
    meta = pd.read_csv(HERE / "data" / "GSE189357_sample_metadata.tsv", sep="\t")
    ex = RAW / "GSE189357" / "raw"
    ex.mkdir(parents=True, exist_ok=True)
    needed = ex / f"{meta.gsm.iloc[0]}_{meta.patient.iloc[0]}_matrix.mtx.gz"
    if not needed.exists():
        tar = RAW / "GSE189357" / "GSE189357_RAW.tar"
        log(f"  untar {tar.name}")
        with tarfile.open(tar) as tf:
            tf.extractall(ex)
    out = []
    for row in meta.itertuples(index=False):
        prefix = ex / f"{row.gsm}_{row.patient}"
        mtx = Path(str(prefix) + "_matrix.mtx.gz")
        feat = Path(str(prefix) + "_features.tsv.gz")
        if not mtx.exists():
            hits = list(ex.glob(f"*{row.patient}_matrix.mtx.gz"))
            if not hits:
                log(f"  missing {row.patient}")
                continue
            mtx = hits[0]
            feat = Path(str(mtx).replace("_matrix.mtx.gz", "_features.tsv.gz"))
        log(f"  {row.patient}")
        mm = spio.mmread(mtx).tocsr()
        features = pd.read_csv(feat, sep="\t", header=None)
        symbols = features.iloc[:, 1 if features.shape[1] > 1 else 0].astype(str).str.upper()
        if mm.shape[0] != len(symbols):
            if mm.shape[1] == len(symbols):
                mm = mm.T.tocsr()
            else:
                raise RuntimeError(f"{row.patient} mtx {mm.shape} vs {len(symbols)}")
        n = mm.shape[1]
        lib = np.asarray(mm.sum(axis=0)).ravel().astype(np.float64)
        summed: dict[str, np.ndarray] = {}
        for i, g in enumerate(symbols):
            if g not in keep and g not in ALIASES:
                continue
            arr = np.asarray(mm[i].todense()).ravel()
            summed[g] = arr if g not in summed else summed[g] + arr
        del mm
        named = collapse_named(summed, keep)
        mat, genes = matrix_from_named(named, n)
        comp = marker_compartments(mat, genes)
        packed = pack_unit(str(row.patient), mat, genes, lib, comp)
        if packed:
            packed["cohort"] = "GSE189357"
            out.append(packed)
    return out


def qc_counts(units: list[dict]) -> pd.DataFrame:
    ref = pd.read_csv(HERE / "data" / "reference_cellchat_inventory.tsv", sep="\t")
    ref["patient"] = ref["patient"].astype(str)
    rows = []
    for u in units:
        comp = u["comp"]
        mal = comp == "MAL"
        rows.append(
            {
                "cohort": u["cohort"],
                "patient": u["patient"],
                "n_mal": int(mal.sum()),
                "n_tnk": int((comp == "TNK").sum()),
            }
        )
    inv = pd.DataFrame(rows)
    merged = inv.merge(ref, on=["cohort", "patient"], suffixes=("", "_ref"))
    bad = merged[(merged["n_mal"] != merged["n_mal_ref"]) | (merged["n_tnk"] != merged["n_tnk_ref"])]
    if len(merged) != 65 or not bad.empty:
        log(bad.to_string(index=False))
        raise RuntimeError(f"label QC failed: matched {len(merged)} / mismatches {len(bad)}")
    log("label QC matched 65 / 65")
    return inv


def hill(x: float) -> float:
    if x <= 0:
        return 0.0
    return x / (KH + x)


def gmean(vals: list[float]) -> float:
    arr = np.asarray(vals, dtype=float)
    if arr.size == 0 or np.any(arr <= 0):
        return 0.0
    return float(np.exp(np.mean(np.log(arr))))


def trimmed_mean(v: np.ndarray) -> float:
    v = np.sort(np.asarray(v, dtype=float))
    n = int(v.size)
    if n == 0:
        return 0.0
    k = int(n * TRIM)
    if k == 0 or n - 2 * k < 1:
        return float(v.mean())
    return float(v[k : n - k].mean())


def q_arms(x: np.ndarray, mode: str) -> tuple[np.ndarray, np.ndarray] | None:
    n = int(x.size)
    if n < 4:
        return None
    r = stats.rankdata(np.asarray(x, dtype=float), method="ordinal")
    if mode == "q4q1":
        if n < MIN_MAL:
            return None
        lo_cut = math.floor(n * 0.25)
        hi_cut = math.ceil(n * 0.75)
    elif mode == "decile":
        lo_cut = math.floor(n * 0.10)
        hi_cut = math.ceil(n * 0.90)
    else:
        raise ValueError(mode)
    if lo_cut < MIN_ARM or (n - hi_cut) < MIN_ARM:
        return None
    low = r <= lo_cut
    high = r > hi_cut
    if int(high.sum()) < MIN_ARM or int(low.sum()) < MIN_ARM:
        return None
    if float(x[high].mean()) <= float(x[low].mean()):
        return None
    return high, low


def value_strata(x: np.ndarray, n_target: int = N_STRATA) -> np.ndarray:
    """Quantile strata that do not split tied values."""
    n = int(x.size)
    order = np.argsort(np.asarray(x, dtype=float), kind="mergesort")
    xs = np.asarray(x, dtype=float)[order]
    bins_sorted = np.zeros(n, dtype=np.int32)
    b = 0
    start = 0
    for i in range(n):
        target = (b + 1) * n / n_target
        remaining = n_target - 1 - b
        at_boundary = (i + 1 == n) or (xs[i + 1] != xs[i])
        enough_left = (n - (i + 1)) >= remaining
        if b < n_target - 1 and (i + 1) >= target and at_boundary and enough_left:
            bins_sorted[start : i + 1] = b
            b += 1
            start = i + 1
    bins_sorted[start:] = b
    out = np.empty(n, dtype=np.int32)
    out[order] = bins_sorted
    return out


def demean(y: np.ndarray, strata: np.ndarray) -> np.ndarray:
    out = np.empty(y.shape[0], dtype=float)
    y = np.asarray(y, dtype=float)
    for s in np.unique(strata):
        m = strata == s
        out[m] = y[m] - float(y[m].mean())
    return out


def lin_resid(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    xc = x - x.mean()
    yc = y - y.mean()
    v = float(np.dot(xc, xc))
    if v <= 0:
        return yc
    beta = float(np.dot(xc, yc) / v)
    return yc - beta * xc


def between_r2(y: np.ndarray, strata: np.ndarray) -> float:
    y = np.asarray(y, dtype=float)
    ss = float(np.sum((y - y.mean()) ** 2))
    if ss <= 0:
        return float("nan")
    within = 0.0
    for s in np.unique(strata):
        m = strata == s
        part = y[m]
        within += float(np.sum((part - part.mean()) ** 2))
    return float(1.0 - within / ss)


def prepare_views(unit: dict) -> dict | None:
    comp = unit["comp"]
    mal = np.where(comp == "MAL")[0]
    tnk = np.where(comp == "TNK")[0]
    genes = unit["genes"]
    if mal.size < MIN_MAL or tnk.size < MIN_TNK:
        return None
    if "CLDN4" not in genes or "TACSTD2" not in genes:
        return None
    gindex = {g: i for i, g in enumerate(genes)}
    lib = np.maximum(unit["lib"], 1.0)
    cp = unit["mat"] / lib[:, None] * 1e4
    logv = np.log1p(cp)
    pos = (unit["mat"] > 0).astype(np.float64)
    return {
        "cohort": unit["cohort"],
        "patient": unit["patient"],
        "genes": genes,
        "gindex": gindex,
        "mal": mal,
        "tnk": tnk,
        "log_mal": logv[mal],
        "pos_mal": pos[mal],
        "log_tnk": logv[tnk],
        "pos_tnk": pos[tnk],
    }


def receptor_cache(view: dict) -> dict[str, tuple[float, float]]:
    """token -> (trimmed geometric-mean log, min subunit proportion) on T/NK."""
    gindex = view["gindex"]
    log_tnk = view["log_tnk"]
    pos_tnk = view["pos_tnk"]
    out = {}
    tokens = {rec for _, rec in EDGES}
    for token in tokens:
        parts = subunits(token)
        logs, props = [], []
        ok = True
        for g in parts:
            if g not in gindex:
                ok = False
                break
            i = gindex[g]
            logs.append(trimmed_mean(log_tnk[:, i]))
            props.append(float(pos_tnk[:, i].mean()))
        if not ok:
            out[token] = (0.0, 0.0)
        else:
            out[token] = (gmean(logs), min(props) if props else 0.0)
    return out


def score_arms(
    view: dict,
    recs: dict[str, tuple[float, float]],
    high: np.ndarray,
    low: np.ndarray,
    pos_use: np.ndarray,
    log_use: np.ndarray,
    detect_pos: np.ndarray,
) -> dict[str, float]:
    """Family means. high/low index the malignant axis."""
    gindex = view["gindex"]
    lig_edge: dict[str, dict[str, list[float]]] = {g: {s: [] for s in SCORES} for g in LIGANDS}
    for lig, rec in EDGES:
        if lig not in gindex:
            continue
        gi = gindex[lig]
        rh_t, rec_prop = recs[rec]
        ph = float(detect_pos[high, gi].mean()) if high.size else 0.0
        pl = float(detect_pos[low, gi].mean()) if low.size else 0.0
        detected = rec_prop >= EXPR_PROP and (ph >= EXPR_PROP or pl >= EXPR_PROP)
        if not detected:
            continue
        prop_d = 100.0 * (float(pos_use[high, gi].mean()) - float(pos_use[low, gi].mean()))
        lh = trimmed_mean(log_use[high, gi])
        ll = trimmed_mean(log_use[low, gi])
        hill_d = hill(lh) * hill(rh_t) - hill(ll) * hill(rh_t)
        mean_h = float(log_use[high, gi].mean())
        mean_l = float(log_use[low, gi].mean())
        cpdb_d = 0.5 * (mean_h - mean_l)
        logfc_d = (mean_h - mean_l) / LN2
        lig_edge[lig]["expr_prop_pp"].append(prop_d)
        lig_edge[lig]["cellchat_hill"].append(hill_d)
        lig_edge[lig]["cpdb_means"].append(cpdb_d)
        lig_edge[lig]["liana_logfc"].append(logfc_d)
    out: dict[str, float] = {}
    for score in SCORES:
        vals = []
        for lig in LIGANDS:
            edges = lig_edge[lig][score]
            vals.append(float(np.mean(edges)) if edges else 0.0)
            out[f"{score}__{lig}"] = vals[-1]
        out[score] = float(np.mean(vals))
    return out


def local_family(
    view: dict,
    recs: dict[str, tuple[float, float]],
    strata: np.ndarray,
    within: np.ndarray,
) -> tuple[dict[str, float] | None, int]:
    """Unweighted mean of within-stratum median splits on `within`."""
    acc = {s: [] for s in SCORES}
    lig_acc = {f"{s}__{g}": [] for s in SCORES for g in LIGANDS}
    n_used = 0
    for s in np.unique(strata):
        idx = np.where(strata == s)[0]
        if idx.size < 2 * MIN_LOCAL:
            continue
        r = stats.rankdata(within[idx], method="ordinal")
        hi = idx[r > (idx.size / 2.0)]
        lo = idx[r <= (idx.size / 2.0)]
        if hi.size < MIN_LOCAL or lo.size < MIN_LOCAL:
            continue
        if float(within[hi].mean()) <= float(within[lo].mean()):
            continue
        scored = score_arms(view, recs, hi, lo, view["pos_mal"], view["log_mal"], view["pos_mal"])
        n_used += 1
        for score in SCORES:
            acc[score].append(scored[score])
            for lig in LIGANDS:
                lig_acc[f"{score}__{lig}"].append(scored[f"{score}__{lig}"])
    if n_used < MIN_LOCAL_STRATA:
        return None, n_used
    out = {score: float(np.mean(acc[score])) for score in SCORES}
    for key, vals in lig_acc.items():
        out[key] = float(np.mean(vals)) if vals else 0.0
    return out, n_used


def match_q4q1(
    high: np.ndarray, low: np.ndarray, strata: np.ndarray, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray] | None:
    hi_parts, lo_parts = [], []
    for s in np.unique(strata):
        h = np.where(high & (strata == s))[0]
        l = np.where(low & (strata == s))[0]
        if h.size == 0 or l.size == 0:
            continue
        n = int(min(h.size, l.size))
        h_take = h if h.size == n else rng.choice(h, size=n, replace=False)
        l_take = l if l.size == n else rng.choice(l, size=n, replace=False)
        hi_parts.append(np.asarray(h_take))
        lo_parts.append(np.asarray(l_take))
    if not hi_parts:
        return None
    hi = np.concatenate(hi_parts)
    lo = np.concatenate(lo_parts)
    if hi.size < MIN_PAIRS:
        return None
    return hi, lo


def gap(x: np.ndarray, high: np.ndarray, low: np.ndarray) -> float:
    return float(x[high].mean() - x[low].mean())


def analyze_unit(view: dict, rng: np.random.Generator) -> tuple[list[dict], dict]:
    gindex = view["gindex"]
    log_mal = view["log_mal"]
    pos_mal = view["pos_mal"]
    tac = log_mal[:, gindex["TACSTD2"]]
    cldn = log_mal[:, gindex["CLDN4"]]
    recs = receptor_cache(view)
    cldn_strata = value_strata(cldn)
    tac_strata = value_strata(tac)
    rows = []
    info = {
        "cohort": view["cohort"],
        "patient": view["patient"],
        "n_mal": int(log_mal.shape[0]),
        "n_tnk": int(view["tnk"].size),
        "spearman_tac_cldn": float(stats.spearmanr(tac, cldn).statistic),
        "r2_tac_on_cldn_strata": between_r2(tac, cldn_strata),
        "r2_cldn_on_tac_strata": between_r2(cldn, tac_strata),
        "pct_tac": float(pos_mal[:, gindex["TACSTD2"]].mean()),
        "pct_cldn": float(pos_mal[:, gindex["CLDN4"]].mean()),
        "n_cldn_strata": int(np.unique(cldn_strata).size),
        "n_tac_strata": int(np.unique(tac_strata).size),
    }
    gates = {
        "TACSTD2": (tac, cldn, cldn_strata, tac),
        "CLDN4": (cldn, tac, tac_strata, cldn),
    }
    # gates: name -> (gate values, covariate values, covariate strata, within-gene for local split)
    for mode in ("q4q1", "decile"):
        for gate, (gvals, cvals, cov_strata, within) in gates.items():
            arms = q_arms(gvals, mode)
            if arms is None:
                continue
            high, low = arms
            crude = score_arms(view, recs, np.where(high)[0], np.where(low)[0], pos_mal, log_mal, pos_mal)
            # Residualize each ligand on the covariate. Detection stays on the raw arms.
            pos_s = pos_mal.copy()
            log_s = log_mal.copy()
            pos_l = pos_mal.copy()
            log_l = log_mal.copy()
            for lig in LIGANDS:
                if lig not in gindex:
                    continue
                gi = gindex[lig]
                pos_s[:, gi] = demean(pos_mal[:, gi], cov_strata)
                log_s[:, gi] = demean(log_mal[:, gi], cov_strata) + float(log_mal[:, gi].mean())
                pos_l[:, gi] = lin_resid(pos_mal[:, gi], cvals)
                log_l[:, gi] = lin_resid(log_mal[:, gi], cvals) + float(log_mal[:, gi].mean())
            hi_idx = np.where(high)[0]
            lo_idx = np.where(low)[0]
            strata_sc = score_arms(view, recs, hi_idx, lo_idx, pos_s, log_s, pos_mal)
            linear_sc = score_arms(view, recs, hi_idx, lo_idx, pos_l, log_l, pos_mal)
            matched = match_q4q1(high, low, cov_strata, rng)
            local, n_local = local_family(view, recs, cov_strata, within)
            pieces = {
                "crude": (crude, hi_idx, lo_idx, 0, float("nan")),
                "strata_resid": (strata_sc, hi_idx, lo_idx, 0, float("nan")),
                "linear_resid": (linear_sc, hi_idx, lo_idx, 0, float("nan")),
            }
            if matched is not None:
                mh, ml = matched
                matched_sc = score_arms(view, recs, mh, ml, pos_mal, log_mal, pos_mal)
                pieces["matched"] = (matched_sc, mh, ml, int(mh.size), float(mh.size / high.sum()))
            if local is not None:
                pieces["strata_local"] = (local, hi_idx, lo_idx, n_local, float("nan"))
            for contrast, (sc, hix, lix, extra, rate) in pieces.items():
                row = {
                    "cohort": view["cohort"],
                    "patient": view["patient"],
                    "split": mode,
                    "gate": gate,
                    "contrast": contrast,
                    "n_high": int(hix.size) if contrast != "strata_local" else int(high.sum()),
                    "n_low": int(lix.size) if contrast != "strata_local" else int(low.sum()),
                    "n_tnk": int(view["tnk"].size),
                    "gate_gap": gap(gvals, high, low) if contrast != "matched" else gap(gvals, _mask_from_idx(gvals.size, hix), _mask_from_idx(gvals.size, lix)),
                    "cov_gap": gap(cvals, high, low) if contrast != "matched" else gap(cvals, _mask_from_idx(cvals.size, hix), _mask_from_idx(cvals.size, lix)),
                    "n_extra": int(extra),
                    "match_rate": rate,
                }
                if contrast == "matched":
                    row["n_high"] = int(hix.size)
                    row["n_low"] = int(lix.size)
                for score in SCORES:
                    row[score] = sc[score]
                    for lig in LIGANDS:
                        row[f"{score}__{lig}"] = sc[f"{score}__{lig}"]
                rows.append(row)
    return rows, info


def _mask_from_idx(n: int, idx: np.ndarray) -> np.ndarray:
    m = np.zeros(n, dtype=bool)
    m[idx] = True
    return m


def wilcox_p(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 8 or np.allclose(x, 0):
        return float("nan")
    try:
        return float(stats.wilcoxon(x, zero_method="wilcox", alternative="two-sided", method="auto").pvalue)
    except ValueError:
        return float("nan")


def signflip_p(x: np.ndarray, rng: np.random.Generator) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 8:
        return float("nan")
    obs = float(x.mean())
    signs = rng.choice(np.array([-1.0, 1.0]), size=(N_FLIP, x.size))
    null = signs @ x / x.size
    return float((np.sum(null >= obs) + 1) / (N_FLIP + 1))


def cohort_i2(patient: pd.DataFrame, value: str) -> float:
    means, vars_ = [], []
    for cohort in COHORTS:
        part = patient.loc[patient["cohort"] == cohort, value].to_numpy(float)
        part = part[np.isfinite(part)]
        if part.size < 3:
            continue
        sd = float(part.std(ddof=1))
        se = sd / math.sqrt(part.size)
        se = max(se, 1e-8)
        means.append(float(part.mean()))
        vars_.append(se ** 2)
    y = np.asarray(means, dtype=float)
    v = np.asarray(vars_, dtype=float)
    if y.size < 2:
        return float("nan")
    w = 1.0 / v
    fe = float(np.sum(w * y) / np.sum(w))
    q = float(np.sum(w * (y - fe) ** 2))
    if q <= 0:
        return 0.0
    return float(max(0.0, (q - (y.size - 1)) / q) * 100)


def summarize(patient: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    keys = ["split", "gate", "contrast", "score"]
    long = patient.melt(
        id_vars=["cohort", "patient", "split", "gate", "contrast"],
        value_vars=SCORES,
        var_name="score",
        value_name="delta",
    )
    for key, g in long.groupby(keys, sort=False):
        x = g["delta"].to_numpy(float)
        rec = {
            "split": key[0],
            "gate": key[1],
            "contrast": key[2],
            "score": key[3],
            "n": int(np.isfinite(x).sum()),
            "mean": float(np.nanmean(x)),
            "median": float(np.nanmedian(x)),
            "frac_pos": float(np.mean(x[np.isfinite(x)] > 0)) if np.isfinite(x).any() else float("nan"),
            "p_wilcox": wilcox_p(x),
            "i2": cohort_i2(g, "delta"),
        }
        n_pos = 0
        for cohort in COHORTS:
            part = g.loc[g["cohort"] == cohort, "delta"].to_numpy(float)
            rec[f"n_{cohort}"] = int(np.isfinite(part).sum())
            rec[f"mean_{cohort}"] = float(np.nanmean(part)) if np.isfinite(part).any() else float("nan")
            if np.isfinite(part).sum() >= 3 and np.nanmean(part) > 0:
                n_pos += 1
        rec["n_cohorts_pos"] = n_pos
        do_flip = (
            key[0] == "q4q1"
            and key[3] in COPRIMARY
            and key[2] in ("crude", "strata_resid", "linear_resid")
        )
        rec["p_signflip"] = signflip_p(x, rng) if do_flip else float("nan")
        rows.append(rec)
    return pd.DataFrame(rows)


def shrink_table(patient: pd.DataFrame) -> pd.DataFrame:
    base = patient[patient["split"] == "q4q1"]
    rows = []
    for gate in ("TACSTD2", "CLDN4"):
        for score in SCORES:
            crude = base[(base["gate"] == gate) & (base["contrast"] == "crude")][
                ["cohort", "patient", score]
            ].rename(columns={score: "crude"})
            for contrast in ("strata_resid", "linear_resid", "matched", "strata_local"):
                adj = base[(base["gate"] == gate) & (base["contrast"] == contrast)][
                    ["patient", score, "match_rate", "cov_gap", "gate_gap"]
                ].rename(columns={score: "adjusted"})
                m = crude.merge(adj, on="patient", how="inner")
                m = m[np.isfinite(m["crude"]) & np.isfinite(m["adjusted"])]
                if m.empty:
                    continue
                drop = m["crude"].to_numpy(float) - m["adjusted"].to_numpy(float)
                mean_c = float(m["crude"].mean())
                mean_a = float(m["adjusted"].mean())
                frac = float(1.0 - mean_a / mean_c) if mean_c > 0 else float("nan")
                p = wilcox_p(drop)
                if not np.isfinite(mean_c) or mean_c <= 0:
                    call = "no positive crude gap"
                elif np.isfinite(frac) and frac >= SHRINK_MIN and np.isfinite(p) and p < 0.05 and (mean_c - mean_a) > 0:
                    call = "shrinks"
                elif np.isfinite(frac) and frac <= 0:
                    call = "does not shrink"
                else:
                    call = "does not materially shrink"
                rem_p = wilcox_p(m["adjusted"].to_numpy(float))
                rows.append(
                    {
                        "gate": gate,
                        "adjusted_for": "CLDN4" if gate == "TACSTD2" else "TACSTD2",
                        "contrast": contrast,
                        "score": score,
                        "n": int(len(m)),
                        "mean_crude": mean_c,
                        "mean_adjusted": mean_a,
                        "shrink_amount": mean_c - mean_a,
                        "shrink_frac": frac,
                        "frac_patients_drop": float(np.mean(drop > 0)),
                        "p_wilcox_drop": p,
                        "p_wilcox_adjusted": rem_p,
                        "mean_cov_gap_adjusted": float(m["cov_gap"].mean()),
                        "mean_gate_gap_adjusted": float(m["gate_gap"].mean()),
                        "mean_match_rate": float(np.nanmean(m["match_rate"])) if contrast == "matched" else float("nan"),
                        "call": call,
                        "remainder_positive": bool(mean_a > 0 and np.isfinite(rem_p) and rem_p < 0.05),
                    }
                )
    return pd.DataFrame(rows)


def ligand_summary(patient: pd.DataFrame) -> pd.DataFrame:
    rows = []
    sub = patient[(patient["split"] == "q4q1") & (patient["contrast"].isin(["crude", "strata_resid", "linear_resid", "matched"]))]
    for (gate, contrast, score), g in sub.groupby(["gate", "contrast", "score"], sort=False):
        for lig in LIGANDS:
            col = f"{score}__{lig}"
            x = g[col].to_numpy(float)
            rec = {
                "gate": gate,
                "contrast": contrast,
                "score": score,
                "ligand": lig,
                "n": int(np.isfinite(x).sum()),
                "mean": float(np.nanmean(x)),
                "frac_pos": float(np.mean(x[np.isfinite(x)] > 0)),
                "p_wilcox": wilcox_p(x),
            }
            n_pos = 0
            for cohort in COHORTS:
                part = g.loc[g["cohort"] == cohort, col].to_numpy(float)
                rec[cohort] = float(np.nanmean(part)) if np.isfinite(part).any() else float("nan")
                if np.isfinite(part).sum() >= 3 and float(np.nanmean(part)) > 0:
                    n_pos += 1
            rec["n_cohorts_pos"] = n_pos
            rows.append(rec)
    return pd.DataFrame(rows)


def draw(patient: pd.DataFrame, shrink: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    q = patient[(patient["split"] == "q4q1") & (patient["gate"] == "TACSTD2")]
    crude = q[q["contrast"] == "crude"][["patient", "cohort", "expr_prop_pp", "cellchat_hill"]].rename(
        columns={"expr_prop_pp": "crude_pp", "cellchat_hill": "crude_hill"}
    )
    resid = q[q["contrast"] == "strata_resid"][["patient", "expr_prop_pp", "cellchat_hill"]].rename(
        columns={"expr_prop_pp": "resid_pp", "cellchat_hill": "resid_hill"}
    )
    m = crude.merge(resid, on="patient")
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.2))
    for ax, xc, yc, title, xlab, ylab in (
        (axes[0], "crude_pp", "resid_pp", "Expression proportion", "Crude TACSTD2 Q4−Q1 (pp)", "After CLDN4 strata (pp)"),
        (axes[1], "crude_hill", "resid_hill", "CellChat Hill", "Crude TACSTD2 Q4−Q1", "After CLDN4 strata"),
    ):
        for cohort in COHORTS:
            part = m[m["cohort"] == cohort]
            ax.scatter(part[xc], part[yc], s=28, color=COHORT_COLOR[cohort], label=cohort, zorder=3)
        lims = [
            min(m[xc].min(), m[yc].min()),
            max(m[xc].max(), m[yc].max()),
        ]
        pad = 0.08 * (lims[1] - lims[0] + 1e-6)
        lo, hi = lims[0] - pad, lims[1] + pad
        ax.plot([lo, hi], [lo, hi], color="#888888", lw=0.8, zorder=1)
        ax.axhline(0, color="#bbbbbb", lw=0.6)
        ax.axvline(0, color="#bbbbbb", lw=0.6)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_title(title)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_aspect("equal", adjustable="box")
    axes[1].legend(frameon=False, fontsize=8, loc="upper left")
    fig.suptitle("TACSTD2 barrier outgoing, same patients", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG / "tacstd2_partial_patients.png", dpi=160)
    fig.savefig(FIG / "tacstd2_partial_patients.pdf")
    plt.close(fig)

    # Cohort means for the four headline contrasts.
    want = [
        ("TACSTD2", "crude", "TACSTD2 crude"),
        ("TACSTD2", "strata_resid", "TACSTD2 | CLDN4"),
        ("CLDN4", "crude", "CLDN4 crude"),
        ("CLDN4", "strata_resid", "CLDN4 | TACSTD2"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.0), sharey=False)
    for ax, score, title in (
        (axes[0], "expr_prop_pp", "Barrier ligands, percentage points"),
        (axes[1], "cellchat_hill", "Barrier ligands, CellChat Hill"),
    ):
        xs = np.arange(len(want))
        for cohort in COHORTS:
            ys = []
            for gate, contrast, _ in want:
                part = patient[
                    (patient["split"] == "q4q1")
                    & (patient["gate"] == gate)
                    & (patient["contrast"] == contrast)
                    & (patient["cohort"] == cohort)
                ]
                ys.append(float(part[score].mean()) if len(part) else np.nan)
            ax.plot(xs, ys, marker="o", color=COHORT_COLOR[cohort], label=cohort, lw=1.2)
        means = []
        for gate, contrast, _ in want:
            part = patient[
                (patient["split"] == "q4q1")
                & (patient["gate"] == gate)
                & (patient["contrast"] == contrast)
            ]
            means.append(float(part[score].mean()) if len(part) else np.nan)
        ax.scatter(xs, means, s=60, color="black", zorder=4, label="all patients")
        ax.axhline(0, color="#bbbbbb", lw=0.6)
        ax.set_xticks(xs)
        ax.set_xticklabels([w[2] for w in want], rotation=20, ha="right")
        ax.set_title(title)
        ax.set_ylabel("Family mean, high − low")
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "barrier_partial_cohorts.png", dpi=160)
    fig.savefig(FIG / "barrier_partial_cohorts.pdf")
    plt.close(fig)
    log(f"figures {FIG}")


def self_check() -> None:
    rng = np.random.default_rng(0)
    n_mal, n_tnk = 400, 80
    # Ligand is the CLDN4 stratum indicator. TACSTD2 copies CLDN4.
    cldn = np.repeat(np.arange(4, dtype=float), n_mal // 4)
    tac = cldn.copy()
    ligand = (cldn >= 2).astype(float)
    # Receptors on in every T/NK cell.
    genes = ["CLDN4", "TACSTD2", "F11R", "NECTIN2", "CDH1", "LGALS9", "ITGAL", "ITGB2", "TIGIT", "CD96", "ITGAE", "ITGB7", "KLRG1", "HAVCR2", "CD44", "PTPRC"]
    gindex = {g: i for i, g in enumerate(genes)}
    log_mal = np.zeros((n_mal, len(genes)))
    pos_mal = np.zeros((n_mal, len(genes)))
    log_mal[:, gindex["CLDN4"]] = cldn
    log_mal[:, gindex["TACSTD2"]] = tac
    for lig in LIGANDS:
        pos_mal[:, gindex[lig]] = ligand
        log_mal[:, gindex[lig]] = ligand * 2.0
    log_tnk = np.ones((n_tnk, len(genes)))
    pos_tnk = np.ones((n_tnk, len(genes)))
    view = {
        "cohort": "SYN",
        "patient": "S",
        "genes": genes,
        "gindex": gindex,
        "mal": np.arange(n_mal),
        "tnk": np.arange(n_tnk),
        "log_mal": log_mal,
        "pos_mal": pos_mal,
        "log_tnk": log_tnk,
        "pos_tnk": pos_tnk,
    }
    rows, info = analyze_unit(view, rng)
    df = pd.DataFrame(rows)
    crude = df[(df.split == "q4q1") & (df.gate == "TACSTD2") & (df.contrast == "crude")].iloc[0]
    resid = df[(df.split == "q4q1") & (df.gate == "TACSTD2") & (df.contrast == "strata_resid")].iloc[0]
    if not (crude["expr_prop_pp"] > 50 and abs(resid["expr_prop_pp"]) < 1e-6):
        raise RuntimeError(f"self-check collinear failed crude {crude['expr_prop_pp']} resid {resid['expr_prop_pp']}")
    # Independent ligand: positivity is TACSTD2 only, CLDN4 is a shuffled copy.
    tac2 = np.tile(np.linspace(0, 3, n_mal // 4), 4)
    rng2 = np.random.default_rng(1)
    cldn2 = rng2.permutation(tac2)
    log_mal2 = np.zeros_like(log_mal)
    pos_mal2 = np.zeros_like(pos_mal)
    log_mal2[:, gindex["TACSTD2"]] = tac2
    log_mal2[:, gindex["CLDN4"]] = cldn2
    ind = (tac2 > np.median(tac2)).astype(float)
    for lig in LIGANDS:
        pos_mal2[:, gindex[lig]] = ind
        log_mal2[:, gindex[lig]] = ind * 2.0
    view["log_mal"] = log_mal2
    view["pos_mal"] = pos_mal2
    rows2, _ = analyze_unit(view, np.random.default_rng(2))
    df2 = pd.DataFrame(rows2)
    c2 = df2[(df2.split == "q4q1") & (df2.gate == "TACSTD2") & (df2.contrast == "crude")].iloc[0]
    r2 = df2[(df2.split == "q4q1") & (df2.gate == "TACSTD2") & (df2.contrast == "strata_resid")].iloc[0]
    if not (c2["expr_prop_pp"] > 50 and r2["expr_prop_pp"] > 0.8 * c2["expr_prop_pp"]):
        raise RuntimeError(f"self-check independent failed crude {c2['expr_prop_pp']} resid {r2['expr_prop_pp']}")
    # Constant covariate: residual equals crude.
    log_mal3 = log_mal2.copy()
    log_mal3[:, gindex["CLDN4"]] = 1.0
    view["log_mal"] = log_mal3
    rows3, _ = analyze_unit(view, np.random.default_rng(3))
    df3 = pd.DataFrame(rows3)
    c3 = df3[(df3.split == "q4q1") & (df3.gate == "TACSTD2") & (df3.contrast == "crude")].iloc[0]
    r3 = df3[(df3.split == "q4q1") & (df3.gate == "TACSTD2") & (df3.contrast == "strata_resid")].iloc[0]
    if abs(c3["expr_prop_pp"] - r3["expr_prop_pp"]) > 1e-6:
        raise RuntimeError(f"self-check constant covariate failed {c3['expr_prop_pp']} vs {r3['expr_prop_pp']}")
    # Strata do not split ties.
    tied = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0])
    st = value_strata(tied, 4)
    if st[0] != st[1] or st[0] != st[2] or st[0] != st[3]:
        raise RuntimeError(f"tie strata split zeros: {st}")
    log(f"self-check ok  collinear resid {resid['expr_prop_pp']:.3g}  independent resid/crude {r2['expr_prop_pp']/c2['expr_prop_pp']:.3f}")


def versions() -> None:
    import scipy
    import sys

    text = "\n".join(
        [
            f"python {sys.version.split()[0]}",
            f"numpy {np.__version__}",
            f"scipy {scipy.__version__}",
            f"pandas {pd.__version__}",
            "R used only to export GSE205335 dgCMatrix",
        ]
    )
    (HERE / "results" / "versions.txt").write_text(text + "\n")
    log(text)


def main() -> None:
    self_check()
    keep = wanted_genes()
    units = []
    units.extend(load_gse123902(keep))
    units.extend(load_gse131907(keep))
    units.extend(load_gse205335(keep))
    units.extend(load_gse189357(keep))
    qc_counts(units)
    rng = np.random.default_rng(SEED)
    rows, infos = [], []
    for u in sorted(units, key=lambda d: (d["cohort"], d["patient"])):
        view = prepare_views(u)
        if view is None:
            log(f"skip {u['cohort']} {u['patient']}")
            continue
        got, info = analyze_unit(view, rng)
        rows.extend(got)
        infos.append(info)
        log(f"  scored {u['cohort']} {u['patient']} rows {len(got)}")
    patient = pd.DataFrame(rows)
    coexpr = pd.DataFrame(infos)
    if patient.empty:
        raise RuntimeError("no scored patients")
    flip_rng = np.random.default_rng(SEED)
    summary = summarize(patient, flip_rng)
    shrink = shrink_table(patient)
    ligands = ligand_summary(patient)
    TAB.mkdir(parents=True, exist_ok=True)
    patient.to_csv(TAB / "per_patient_contrasts.tsv", sep="\t", index=False)
    coexpr.to_csv(TAB / "coexpression.tsv", sep="\t", index=False)
    summary.to_csv(TAB / "summary_contrasts.tsv", sep="\t", index=False)
    shrink.to_csv(TAB / "shrink.tsv", sep="\t", index=False)
    ligands.to_csv(TAB / "summary_ligands.tsv", sep="\t", index=False)
    draw(patient, shrink)
    versions()
    # Locked CLDN4 Q4 expression-proportion QC against PR 713.
    hit = summary[
        (summary["split"] == "q4q1")
        & (summary["gate"] == "CLDN4")
        & (summary["contrast"] == "crude")
        & (summary["score"] == "expr_prop_pp")
    ]
    if len(hit) != 1:
        raise RuntimeError("missing CLDN4 crude expr_prop summary")
    mean_pp = float(hit["mean"].iloc[0])
    n = int(hit["n"].iloc[0])
    log(f"QC CLDN4 crude expr_prop family mean {mean_pp:.4f} n={n} (PR 713 reported +25.60, n=64)")
    head = shrink[(shrink["gate"] == "TACSTD2") & (shrink["score"].isin(COPRIMARY))]
    log(head.to_string(index=False))
    log("done")


if __name__ == "__main__":
    main()
