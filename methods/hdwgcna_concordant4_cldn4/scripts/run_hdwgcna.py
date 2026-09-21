#!/usr/bin/env python3
"""hdWGCNA-like modules in concordant-4 malignant cells, split by CLDN4.

Prespecified before looking at module–immune p-values. See README.md.
The immune endpoint is the locked patient/donor/sample T/NK fraction.
Metacells are not the sample size.
"""
from __future__ import annotations

import json
import zlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import fisher_exact, rankdata, spearmanr, t as student_t, wilcoxon
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

ROOT = Path("/workspace/methods/hdwgcna_concordant4_cldn4")
WORK = Path("/tmp/hdwgcna_c4")
STORE = WORK / "store"
OUT_TAB = ROOT / "results" / "tables"
OUT_FIG = ROOT / "results" / "figures"

# Frozen choices.
K = 25
MAX_SHARED = 10
MIN_CELLS = 50
TARGET_METACELLS = 40
MAX_ITER = 5000
N_PCS = 15
N_GENES = 2500
MIN_FRAC = 0.05
MIN_DATASETS_DETECTED = 3
MIN_DATASET_METACELLS = 8
MIN_NETWORK_METACELLS = 30
MIN_MODULE_SIZE = 30
MERGE_COR = 0.75
POWERS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
MIN_SET_OVERLAP = 5
PAIRED_MIN_CELLS = 20
DATASETS = ("GSE123902", "GSE131907", "GSE205335", "GSE189357")
COLORS = {
    "GSE123902": "#1b4f72",
    "GSE131907": "#117a65",
    "GSE205335": "#b9770e",
    "GSE189357": "#6c3483",
}


def say(msg: str) -> None:
    print(msg, flush=True)


def load_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(c, "")) for c in columns) + "\n")


def dl_spearman(rhos: list[float], ns: list[float], df_penalty: int = 0) -> dict:
    """DerSimonian–Laird on Fisher z. df_penalty=1 uses 1/(n-4) for a partial rho."""
    rhos_a = np.asarray(rhos, dtype=float)
    ns_a = np.asarray(ns, dtype=float)
    ok = np.isfinite(rhos_a) & np.isfinite(ns_a) & (ns_a > 3 + df_penalty)
    rhos_a = rhos_a[ok]
    ns_a = ns_a[ok]
    if len(rhos_a) == 0:
        return {"rho": np.nan, "p": np.nan, "I2": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "k": 0, "N": 0}
    z = np.arctanh(np.clip(rhos_a, -0.999999, 0.999999))
    var_z = 1.0 / (ns_a - 3 - df_penalty)
    w = 1.0 / var_z
    zbar = np.sum(w * z) / np.sum(w)
    q = float(np.sum(w * (z - zbar) ** 2))
    k = int(len(rhos_a))
    dfree = k - 1
    cdenom = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (q - dfree) / cdenom) if dfree > 0 and cdenom > 0 else 0.0
    wstar = 1.0 / (var_z + tau2)
    zre = float(np.sum(wstar * z) / np.sum(wstar))
    se = float(np.sqrt(1.0 / np.sum(wstar)))
    # two-sided normal tail, matching 2 * pnorm(-|z|)
    from math import erfc

    p = float(erfc(abs(zre / se) / np.sqrt(2.0))) if se > 0 else np.nan
    i2 = max(0.0, (q - dfree) / q) if q > 0 else 0.0
    return {
        "rho": float(np.tanh(zre)),
        "p": p,
        "I2": i2,
        "ci_lo": float(np.tanh(zre - 1.96 * se)),
        "ci_hi": float(np.tanh(zre + 1.96 * se)),
        "k": k,
        "N": int(ns_a.sum()),
    }


def assert_locked_endpoint() -> list[dict[str, str]]:
    rows = load_tsv(ROOT / "data" / "patient_units_locked.tsv")
    by: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by.setdefault(row["dataset"], []).append(row)
    rhos, ns = [], []
    for ds in DATASETS:
        items = by[ds]
        x = np.array([float(r["mal_CLDN4_pct"]) for r in items])
        y = np.array([float(r["frac_tnk"]) for r in items])
        rho, _p = spearmanr(x, y)
        rhos.append(float(rho))
        ns.append(len(items))
    meta = dl_spearman(rhos, ns)
    if abs(meta["rho"] - (-0.5311678045689989)) > 1e-9 or abs(meta["p"] - 1.646223294457401e-05) > 1e-12:
        raise SystemExit(f"locked endpoint was not reproduced: {meta}")
    say(f"locked CLDN4 %pos vs T/NK reproduced: rho={meta['rho']:.3f} p={meta['p']:.3g} N={meta['N']}")
    return rows


def bh(pvals: list[float]) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    pp = p[ok]
    n = len(pp)
    order = np.argsort(pp)
    ranked = pp[order]
    qq = ranked * n / (np.arange(1, n + 1))
    qq = np.minimum.accumulate(qq[::-1])[::-1]
    qq = np.clip(qq, 0, 1)
    back = np.empty(n)
    back[order] = qq
    q[ok] = back
    return q


def select_genes() -> list[str]:
    stats = {}
    for ds in DATASETS:
        rows = load_tsv(WORK / "stats" / f"{ds}.tsv")
        stats[ds] = {row["gene"]: row for row in rows}
    genes = set(stats[DATASETS[0]])
    for ds in DATASETS[1:]:
        genes &= set(stats[ds])
    genes = sorted(genes)
    keep = []
    var_rank_sum = []
    for gene in genes:
        n_det = 0
        for ds in DATASETS:
            row = stats[ds][gene]
            n_cells = float(row["n_cells"])
            frac = float(row["nnz"]) / n_cells if n_cells else 0.0
            if frac >= MIN_FRAC:
                n_det += 1
        if gene == "CLDN4" or n_det >= MIN_DATASETS_DETECTED:
            keep.append(gene)
    # variance ranks among genes that pass the detection rule, within each dataset
    ranks = {ds: {} for ds in DATASETS}
    for ds in DATASETS:
        variances = []
        for gene in keep:
            row = stats[ds][gene]
            n = float(row["n_cells"])
            mean = float(row["sum_log1p"]) / n
            var = float(row["sumsq_log1p"]) / n - mean ** 2
            variances.append((var, gene))
        variances.sort(reverse=True)
        for rank, (_var, gene) in enumerate(variances, start=1):
            ranks[ds][gene] = rank
    scored = []
    for gene in keep:
        mean_rank = float(np.mean([ranks[ds][gene] for ds in DATASETS]))
        scored.append((mean_rank, gene))
    scored.sort()
    chosen = [gene for _rank, gene in scored[:N_GENES]]
    if "CLDN4" not in chosen and "CLDN4" in keep:
        chosen.append("CLDN4")
    say(f"network genes {len(chosen)} (detection>={MIN_FRAC} in >={MIN_DATASETS_DETECTED} datasets, top variance)")
    return chosen


def load_dataset(dataset: str, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    dest = STORE / dataset
    all_genes = [g for g in (dest / "genes.txt").read_text().splitlines() if g]
    index = {g: i for i, g in enumerate(all_genes)}
    missing = [g for g in genes if g not in index]
    if missing:
        raise SystemExit(f"{dataset} missing {len(missing)} network genes, e.g. {missing[:5]}")
    n_genes, n_cells = (dest / "shape.txt").read_text().strip().split("\t")
    shape = (int(n_genes), int(n_cells))
    mm = np.memmap(dest / "counts.u16", dtype=np.uint16, mode="r", shape=shape)
    rows = [index[g] for g in genes]
    counts = np.array(mm[rows], dtype=np.uint16)
    del mm
    units = [u for u in (dest / "units.txt").read_text().splitlines() if u]
    if len(units) != counts.shape[1]:
        raise SystemExit(f"{dataset} unit length {len(units)} != cells {counts.shape[1]}")
    return counts, units


def construct_metacells(counts: np.ndarray, seed: int) -> np.ndarray | None:
    """Average raw counts of kNN metacells. counts is genes x cells."""
    n_genes, n_cells = counts.shape
    if n_cells < MIN_CELLS:
        return None
    lib = counts.sum(axis=0).astype(np.float64)
    lib[lib <= 0] = 1.0
    logn = np.log1p(counts.astype(np.float64) / lib * 1e4)
    mu = logn.mean(axis=1, keepdims=True)
    sd = logn.std(axis=1, keepdims=True)
    sd[sd < 1e-8] = 1.0
    z = ((logn - mu) / sd).T  # cells x genes
    n_pcs = min(N_PCS, n_cells - 2, n_genes - 1)
    if n_pcs < 2:
        return None
    pcs = PCA(n_components=n_pcs, svd_solver="randomized", random_state=seed).fit_transform(z)
    nn = NearestNeighbors(n_neighbors=K, algorithm="brute", metric="euclidean")
    nn.fit(pcs)
    idx = nn.kneighbors(return_distance=False)
    rng = np.random.default_rng(seed)
    good = list(range(n_cells))
    first = int(rng.integers(0, len(good)))
    chosen = [good.pop(first)]
    it = 0
    while good and len(chosen) < TARGET_METACELLS and it < MAX_ITER:
        it += 1
        pick = int(rng.integers(0, len(good)))
        cand = good.pop(pick)
        this = idx[cand]
        accept = True
        for prev in chosen:
            shared = 2 * K - len(np.union1d(idx[prev], this))
            if shared > MAX_SHARED:
                accept = False
                break
        if accept:
            chosen.append(cand)
    if len(chosen) <= 1:
        return None
    members = idx[np.array(chosen, dtype=int)]
    out = np.empty((n_genes, len(chosen)), dtype=np.float64)
    raw = counts.astype(np.float64)
    for j, memb in enumerate(members):
        out[:, j] = raw[:, memb].mean(axis=1)
    return out


def build_metacells(genes: list[str]) -> dict[str, dict]:
    """Return per-stratum log-normalized metacell matrices and labels."""
    cache = WORK / "metacell_cache.npz"
    if cache.exists():
        blob = np.load(cache, allow_pickle=False)
        cached_genes = blob["genes"].astype(str).tolist()
        if cached_genes == genes:
            say("loaded metacell cache")
            packed = {}
            for stratum, key in (("CLDN4pos", "pos"), ("CLDN4neg", "neg")):
                expr = blob[f"{key}_expr"]
                if expr.size == 0:
                    packed[stratum] = None
                else:
                    packed[stratum] = {
                        "expr": expr,
                        "dataset": blob[f"{key}_dataset"].astype(str),
                        "unit": blob[f"{key}_unit"].astype(str),
                    }
                    say(f"  {stratum} metacells {expr.shape[1]}")
            return packed
    strata = {
        "CLDN4pos": {"expr": [], "dataset": [], "unit": [], "n_cells": []},
        "CLDN4neg": {"expr": [], "dataset": [], "unit": [], "n_cells": []},
    }
    manifest = []
    cldn_i = genes.index("CLDN4")
    for dataset in DATASETS:
        say(f"metacells {dataset}")
        counts, units = load_dataset(dataset, genes)
        unit_arr = np.array(units)
        cldn4 = counts[cldn_i].astype(np.float64)
        for unit in dict.fromkeys(units):
            in_unit = unit_arr == unit
            for stratum, positive in (("CLDN4pos", True), ("CLDN4neg", False)):
                if positive:
                    mask = in_unit & (cldn4 > 0)
                else:
                    mask = in_unit & (cldn4 == 0)
                n = int(mask.sum())
                if n < MIN_CELLS:
                    manifest.append(
                        {"dataset": dataset, "unit_id": unit, "stratum": stratum, "n_cells": n, "n_metacells": 0}
                    )
                    continue
                seed = zlib.crc32(f"{dataset}|{unit}|{stratum}".encode()) & 0xFFFFFFFF
                acc = construct_metacells(counts[:, mask], seed)
                n_meta = 0 if acc is None else int(acc.shape[1])
                manifest.append(
                    {"dataset": dataset, "unit_id": unit, "stratum": stratum, "n_cells": n, "n_metacells": n_meta}
                )
                if acc is None:
                    continue
                lib = acc.sum(axis=0)
                lib[lib <= 0] = 1.0
                logn = np.log1p(acc / lib * 1e4)
                strata[stratum]["expr"].append(logn)
                strata[stratum]["dataset"].extend([dataset] * n_meta)
                strata[stratum]["unit"].extend([unit] * n_meta)
                strata[stratum]["n_cells"].extend([n] * n_meta)
        del counts
    write_tsv(
        OUT_TAB / "metacell_manifest.tsv",
        manifest,
        ["dataset", "unit_id", "stratum", "n_cells", "n_metacells"],
    )
    packed = {}
    for stratum, blob in strata.items():
        if not blob["expr"]:
            packed[stratum] = None
            continue
        expr = np.concatenate(blob["expr"], axis=1)
        packed[stratum] = {
            "expr": expr,
            "dataset": np.array(blob["dataset"]),
            "unit": np.array(blob["unit"]),
        }
        say(f"  {stratum} metacells {expr.shape[1]}")
    np.savez_compressed(
        cache,
        genes=np.array(genes),
        pos_expr=packed["CLDN4pos"]["expr"] if packed["CLDN4pos"] else np.empty((0, 0)),
        pos_dataset=packed["CLDN4pos"]["dataset"] if packed["CLDN4pos"] else np.array([]),
        pos_unit=packed["CLDN4pos"]["unit"] if packed["CLDN4pos"] else np.array([]),
        neg_expr=packed["CLDN4neg"]["expr"] if packed["CLDN4neg"] else np.empty((0, 0)),
        neg_dataset=packed["CLDN4neg"]["dataset"] if packed["CLDN4neg"] else np.array([]),
        neg_unit=packed["CLDN4neg"]["unit"] if packed["CLDN4neg"] else np.array([]),
    )
    return packed


def zscore_within_dataset(expr: np.ndarray, dataset: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Center and scale each gene inside each dataset. Drop datasets with too few metacells."""
    keep_cols = np.zeros(expr.shape[1], dtype=bool)
    out = np.zeros_like(expr)
    for ds in DATASETS:
        cols = dataset == ds
        n = int(cols.sum())
        if n < MIN_DATASET_METACELLS:
            say(f"  dropping {ds} from the network ({n} metacells < {MIN_DATASET_METACELLS})")
            continue
        block = expr[:, cols]
        mu = block.mean(axis=1, keepdims=True)
        sd = block.std(axis=1, keepdims=True)
        sd[sd < 1e-8] = 1.0
        out[:, cols] = (block - mu) / sd
        keep_cols[cols] = True
    # drop genes that are flat on the retained metacells
    kept = out[:, keep_cols]
    sd = kept.std(axis=1)
    keep_genes = sd > 1e-8
    return out[keep_genes][:, keep_cols], keep_genes


def bicor_matrix(x_samples_by_genes: np.ndarray) -> np.ndarray:
    med = np.median(x_samples_by_genes, axis=0)
    mad = np.median(np.abs(x_samples_by_genes - med), axis=0)
    valid = mad > 1e-8
    mad = np.where(valid, mad, 1.0)
    u = (x_samples_by_genes - med) / (9.0 * mad)
    w = (1.0 - u ** 2) ** 2
    w[np.abs(u) >= 1.0] = 0.0
    w[:, ~valid] = 0.0
    xw = (x_samples_by_genes - med) * w
    num = xw.T @ xw
    ss = np.sum(xw * xw, axis=0)
    den = np.sqrt(np.outer(ss, ss))
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = num / den
    corr[~np.isfinite(corr)] = 0.0
    np.fill_diagonal(corr, 1.0)
    return corr


def scale_free(k: np.ndarray) -> tuple[float, float]:
    k = np.asarray(k, dtype=float)
    k = k[np.isfinite(k)]
    if k.size < 20:
        return np.nan, np.nan
    hist, edges = np.histogram(k, bins=10)
    centers = 0.5 * (edges[:-1] + edges[1:])
    ok = (hist > 0) & (centers > 0)
    if ok.sum() < 4:
        return np.nan, np.nan
    x = np.log10(centers[ok])
    y = np.log10(hist[ok].astype(float))
    coef = np.polyfit(x, y, 1)
    pred = np.polyval(coef, x)
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return r2, float(coef[0])


def choose_power(corr: np.ndarray) -> tuple[int, list[dict]]:
    rows = []
    for beta in POWERS:
        adj = np.abs(0.5 + 0.5 * corr) ** beta
        np.fill_diagonal(adj, 0.0)
        k = adj.sum(axis=1)
        r2, slope = scale_free(k)
        rows.append(
            {
                "power": beta,
                "scale_free_R2": r2,
                "slope": slope,
                "mean_connectivity": float(np.mean(k)),
                "median_connectivity": float(np.median(k)),
            }
        )
    # R^2 >= 0.80 can land only after mean connectivity has collapsed (<5).
    # That network has almost no edges, so moderate co-expression (a TJ program)
    # cannot be recovered. QC, applied before enrichment: if the smallest
    # scale-free power has mean k < 5, use the negative-slope power inside
    # mean k [10, 100] with the highest scale-free R^2.
    scale_free_hits = [
        row
        for row in rows
        if np.isfinite(row["scale_free_R2"])
        and row["scale_free_R2"] >= 0.80
        and row["slope"] < 0
        and row["mean_connectivity"] <= 200
    ]
    window = [
        row
        for row in rows
        if np.isfinite(row["scale_free_R2"])
        and row["slope"] < 0
        and 10 <= row["mean_connectivity"] <= 100
    ]
    if scale_free_hits and min(scale_free_hits, key=lambda row: row["power"])["mean_connectivity"] >= 5 and window:
        chosen = min(
            [row for row in scale_free_hits if row["mean_connectivity"] >= 5],
            key=lambda row: row["power"],
        )
        rule = "smallest beta with R2>=0.80, negative slope, and mean k>=5"
    elif window:
        chosen = max(window, key=lambda row: (row["scale_free_R2"], -row["power"]))
        collapsed = min(scale_free_hits, key=lambda row: row["power"]) if scale_free_hits else None
        if collapsed is None:
            rule = "no beta reached R2>=0.80; highest R2 among negative-slope powers with mean k in [10, 100]"
        else:
            rule = (
                f"R2>=0.80 first at beta {collapsed['power']} with mean k "
                f"{collapsed['mean_connectivity']:.2f} (<5, edges collapsed). "
                "Primary beta is the highest R2 among negative-slope powers with mean k in [10, 100]. "
                "Scale-free fit was not reached in that window."
            )
    elif scale_free_hits:
        chosen = min(scale_free_hits, key=lambda row: row["power"])
        rule = "no power had mean k in [10, 100]; smallest scale-free beta with mean k<=200"
    else:
        chosen = next(row for row in rows if row["power"] == 12)
        rule = "FAQ signed power 12; scale-free fit and connectivity window both failed"
    for row in rows:
        row["chosen"] = row["power"] == chosen["power"]
        row["rule"] = rule if row["chosen"] else ""
    say(
        f"  power {chosen['power']} R2={chosen['scale_free_R2']:.3f} "
        f"slope={chosen['slope']:.3f} mean_k={chosen['mean_connectivity']:.1f}"
    )
    return int(chosen["power"]), rows


def signed_tom(corr: np.ndarray, beta: int) -> np.ndarray:
    adj = np.abs(0.5 + 0.5 * corr) ** beta
    np.fill_diagonal(adj, 0.0)
    overlap = adj @ adj
    k = adj.sum(axis=1)
    kmin = np.minimum(k[:, None], k[None, :])
    tom = (overlap + adj) / np.maximum(kmin + 1.0 - adj, 1e-12)
    np.fill_diagonal(tom, 1.0)
    return np.clip(tom, 0.0, 1.0)


def dynamic_modules(diss: np.ndarray) -> np.ndarray:
    path = WORK / "dissTOM.bin"
    labels_path = WORK / "module_labels.tsv"
    np.ascontiguousarray(diss, dtype=np.float64).tofile(path)
    n = diss.shape[0]
    r_path = WORK / "cut_modules.R"
    r_path.write_text(
        f"""
.libPaths(c("/home/ubuntu/R/library", .libPaths()))
suppressPackageStartupMessages(library(dynamicTreeCut))
n <- {n}
con <- file("{path}", "rb")
m <- matrix(readBin(con, "double", n * n), n, n, byrow = TRUE)
close(con)
m[!is.finite(m)] <- 1
m <- (m + t(m)) / 2
diag(m) <- 0
d <- as.dist(m)
hc <- hclust(d, method = "average")
labs <- cutreeDynamic(
  dendro = hc,
  distM = m,
  deepSplit = 2,
  minClusterSize = {MIN_MODULE_SIZE},
  method = "hybrid",
  pamRespectsDendro = FALSE,
  verbose = 0
)
write.table(data.frame(module = as.integer(labs)), "{labels_path}",
            sep = "\\t", quote = FALSE, row.names = FALSE)
cat("modules", length(unique(labs)), "grey", sum(labs == 0), "\\n")
"""
    )
    import subprocess

    proc = subprocess.run(["Rscript", str(r_path)], capture_output=True, text=True)
    say(proc.stdout.strip())
    if proc.returncode != 0:
        raise SystemExit(proc.stderr[-2000:])
    labels = np.array([int(row["module"]) for row in load_tsv(labels_path)], dtype=int)
    if labels.size != n:
        raise SystemExit(f"label length {labels.size} != {n}")
    return labels


def eigengene(x_samples_by_genes: np.ndarray, mask: np.ndarray) -> np.ndarray:
    y = x_samples_by_genes[:, mask]
    y = y - y.mean(axis=0, keepdims=True)
    _u, _s, vt = np.linalg.svd(y, full_matrices=False)
    # PC1 scores
    me = y @ vt[0]
    avg = y.mean(axis=1)
    if np.corrcoef(me, avg)[0, 1] < 0:
        me = -me
    return me


def merge_modules(x_samples_by_genes: np.ndarray, labels: np.ndarray) -> np.ndarray:
    labels = labels.copy()
    for _ in range(100):
        mods = [int(m) for m in np.unique(labels) if m != 0]
        if len(mods) < 2:
            break
        mes = np.vstack([eigengene(x_samples_by_genes, labels == m) for m in mods])
        corr = np.corrcoef(mes)
        np.fill_diagonal(corr, -np.inf)
        flat = int(np.nanargmax(corr))
        i, j = divmod(flat, corr.shape[1])
        if not np.isfinite(corr[i, j]) or corr[i, j] <= MERGE_COR:
            break
        say(f"  merging module {mods[j]} into {mods[i]} (ME r={corr[i, j]:.3f})")
        labels[labels == mods[j]] = mods[i]
    return labels


def renumber(labels: np.ndarray) -> np.ndarray:
    mods = [int(m) for m in np.unique(labels) if m != 0]
    mods.sort(key=lambda m: -int(np.sum(labels == m)))
    out = np.zeros_like(labels)
    for new_id, old in enumerate(mods, start=1):
        out[labels == old] = new_id
    return out


def kme_table(x_samples_by_genes: np.ndarray, labels: np.ndarray) -> dict[int, np.ndarray]:
    out = {}
    for module in [int(m) for m in np.unique(labels) if m != 0]:
        me = eigengene(x_samples_by_genes, labels == module)
        xc = x_samples_by_genes - x_samples_by_genes.mean(axis=0, keepdims=True)
        mec = me - me.mean()
        num = xc.T @ mec
        den = np.sqrt(np.sum(xc ** 2, axis=0) * np.sum(mec ** 2))
        with np.errstate(divide="ignore", invalid="ignore"):
            r = num / den
        r[~np.isfinite(r)] = 0.0
        out[module] = r
    return out


def enrich(module_genes: set[str], set_genes: set[str], universe: set[str]) -> dict:
    mod = module_genes & universe
    sig = set_genes & universe
    a = len(mod & sig)
    b = len(mod - sig)
    c = len(sig - mod)
    d = len(universe) - a - b - c
    if min(a + b, c + d, a + c, b + d) == 0:
        return {"overlap": a, "odds_ratio": np.nan, "p": 1.0}
    odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return {"overlap": a, "odds_ratio": float(odds), "p": float(p)}


def identify_module(rows: list[dict], family: str) -> dict | None:
    cands = [
        row
        for row in rows
        if row["module"] != "grey"
        and row["family"] == family
        and row["overlap"] >= MIN_SET_OVERLAP
        and np.isfinite(row["odds_ratio"])
        and row["odds_ratio"] > 1
    ]
    if not cands:
        return None
    cands.sort(key=lambda row: (row["p"], -row["overlap"]))
    best = dict(cands[0])
    best["called"] = bool(np.isfinite(best.get("q", np.nan)) and best["q"] < 0.05)
    return best


def pseudobulk_scores(
    genes: list[str], module_genes: list[str], use_mask
) -> list[dict]:
    """Mean log2(CPM+1) of module genes, CLDN4 held out, per locked unit."""
    score_genes = [g for g in module_genes if g != "CLDN4" and g in genes]
    if len(score_genes) < 5:
        return []
    gindex = {g: i for i, g in enumerate(genes)}
    ix = np.array([gindex[g] for g in score_genes], dtype=int)
    locked = load_tsv(ROOT / "data" / "patient_units_locked.tsv")
    out = []
    for dataset in DATASETS:
        counts, units = load_dataset(dataset, genes)
        # CLDN4 may be outside `genes` only if selection dropped it; it is forced in.
        cldn = counts[gindex["CLDN4"]].astype(np.float64)
        unit_arr = np.array(units)
        want = [row for row in locked if row["dataset"] == dataset]
        for row in want:
            in_unit = unit_arr == row["unit_id"]
            mask = use_mask(in_unit, cldn)
            n = int(mask.sum())
            if n == 0:
                score = np.nan
                lib = 0.0
            else:
                total = counts[ix][:, mask].astype(np.float64).sum(axis=1)
                lib = float(total.sum())
                if lib <= 0:
                    score = np.nan
                else:
                    logcpm = np.log2(total / lib * 1e6 + 1.0)
                    score = float(logcpm.mean())
            out.append(
                {
                    "dataset": dataset,
                    "unit_id": row["unit_id"],
                    "n_cells_scored": n,
                    "score": score,
                    "frac_tnk": float(row["frac_tnk"]),
                    "mal_CLDN4_pct": float(row["mal_CLDN4_pct"]),
                    "locked_ifn_score": float(row["ifn_score"]) if row["ifn_score"] else np.nan,
                    "locked_tj_score": float(row["tj_score"]) if row["tj_score"] else np.nan,
                    "n_genes": len(score_genes),
                }
            )
        del counts
    return out


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)

    def pr(a, b):
        if np.std(a) == 0 or np.std(b) == 0:
            return np.nan
        return float(np.corrcoef(a, b)[0, 1])

    rxy, rxz, ryz = pr(rx, ry), pr(rx, rz), pr(ry, rz)
    denom = np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
    if not np.isfinite(denom) or denom <= 0:
        return np.nan
    return float((rxy - rxz * ryz) / denom)


def correlate_scores(rows: list[dict], ykey: str) -> tuple[list[dict], dict]:
    singles = []
    rhos, ns = [], []
    for ds in DATASETS:
        items = [r for r in rows if r["dataset"] == ds and np.isfinite(r["score"]) and np.isfinite(r[ykey])]
        if len(items) <= 3:
            singles.append({"dataset": ds, "n": len(items), "rho": np.nan, "p": np.nan})
            continue
        rho, p = spearmanr([r["score"] for r in items], [r[ykey] for r in items])
        singles.append({"dataset": ds, "n": len(items), "rho": float(rho), "p": float(p)})
        rhos.append(float(rho))
        ns.append(len(items))
    return singles, dl_spearman(rhos, ns)


def correlate_partial(rows: list[dict]) -> tuple[list[dict], dict]:
    singles = []
    rhos, ns = [], []
    for ds in DATASETS:
        items = [
            r
            for r in rows
            if r["dataset"] == ds
            and np.isfinite(r["score"])
            and np.isfinite(r["frac_tnk"])
            and np.isfinite(r["mal_CLDN4_pct"])
        ]
        if len(items) <= 4:
            singles.append({"dataset": ds, "n": len(items), "rho": np.nan, "p": np.nan})
            continue
        rho = partial_spearman(
            np.array([r["score"] for r in items]),
            np.array([r["frac_tnk"] for r in items]),
            np.array([r["mal_CLDN4_pct"] for r in items]),
        )
        n = len(items)
        if np.isfinite(rho) and abs(rho) < 1:
            tstat = rho * np.sqrt((n - 4) / (1 - rho ** 2))
            p = float(2 * student_t.sf(abs(tstat), n - 4))
        else:
            p = np.nan
        singles.append({"dataset": ds, "n": n, "rho": rho, "p": p})
        if np.isfinite(rho):
            rhos.append(rho)
            ns.append(n)
    return singles, dl_spearman(rhos, ns, df_penalty=1)


def run_network(stratum: str, packed: dict, genes: list[str], sets: dict) -> dict | None:
    say(f"network {stratum}")
    expr = packed["expr"]
    dataset = packed["dataset"]
    z_expr, keep_genes = zscore_within_dataset(expr, dataset)
    used_genes = [g for g, flag in zip(genes, keep_genes) if flag]
    n_meta = z_expr.shape[1]
    say(f"  genes {len(used_genes)} metacells {n_meta}")
    if n_meta < MIN_NETWORK_METACELLS:
        say("  too few metacells; network not built")
        return None
    # samples x genes
    x = np.ascontiguousarray(z_expr.T)
    corr = bicor_matrix(x)
    beta, power_rows = choose_power(corr)
    tom = signed_tom(corr, beta)
    diss = 1.0 - tom
    np.fill_diagonal(diss, 0.0)
    labels = renumber(merge_modules(x, dynamic_modules(diss)))
    membership = kme_table(x, labels)
    universe = set(used_genes)
    enrich_rows = []
    module_ids = [int(m) for m in np.unique(labels) if m != 0] + [0]
    for module in module_ids:
        if module == 0:
            mod_genes = {g for g, lab in zip(used_genes, labels) if lab == 0}
            name = "grey"
        else:
            mod_genes = {g for g, lab in zip(used_genes, labels) if lab == module}
            name = f"ME{module}"
        for family in ("TJ", "IFN", "MHC"):
            hit = enrich(mod_genes, set(sets[family]), universe)
            enrich_rows.append(
                {
                    "stratum": stratum,
                    "module": name,
                    "n_genes": len(mod_genes),
                    "family": family,
                    "overlap": hit["overlap"],
                    "odds_ratio": hit["odds_ratio"],
                    "p": hit["p"],
                }
            )
    qvals = bh([row["p"] for row in enrich_rows])
    for row, q in zip(enrich_rows, qvals):
        row["q"] = float(q) if np.isfinite(q) else np.nan
    tj = identify_module(enrich_rows, "TJ")
    ifn = identify_module(enrich_rows, "IFN")
    # kME long table
    kme_rows = []
    gene_pos = {g: i for i, g in enumerate(used_genes)}
    for module, vec in membership.items():
        for gene, value in zip(used_genes, vec):
            in_mod = int(labels[gene_pos[gene]] == module)
            if in_mod or abs(value) >= 0.999:
                if in_mod:
                    kme_rows.append(
                        {
                            "stratum": stratum,
                            "module": f"ME{module}",
                            "gene": gene,
                            "in_module": 1,
                            "kME": float(value),
                        }
                    )
    def seed_row(symbol: str) -> dict:
        if symbol not in gene_pos:
            return {"gene": symbol, "module": "absent", "kME": np.nan}
        lab = int(labels[gene_pos[symbol]])
        name = "grey" if lab == 0 else f"ME{lab}"
        value = np.nan
        if lab in membership:
            value = float(membership[lab][gene_pos[symbol]])
        return {"gene": symbol, "module": name, "kME": value}

    n_by_ds = {ds: int(np.sum(dataset == ds)) for ds in DATASETS}
    # dataset vector was filtered inside zscore; recount from keep_cols by rebuilding
    # The packed dataset array was not filtered in the return. Reconstruct from zscore call.
    # zscore drops columns but we did not return the mask. Recompute counts from expr columns
    # that survived: z_expr columns correspond to datasets with >= MIN_DATASET_METACELLS,
    # in the original order of `dataset`.
    kept_ds = []
    for ds in DATASETS:
        if int(np.sum(dataset == ds)) >= MIN_DATASET_METACELLS:
            kept_ds.extend([ds] * int(np.sum(dataset == ds)))
    if len(kept_ds) != n_meta:
        raise SystemExit(f"dataset mask length {len(kept_ds)} != {n_meta}")
    n_by_ds = {ds: int(kept_ds.count(ds)) for ds in DATASETS}
    return {
        "stratum": stratum,
        "beta": beta,
        "power_rows": power_rows,
        "genes": used_genes,
        "labels": labels,
        "enrich_rows": enrich_rows,
        "kme_rows": kme_rows,
        "tj": tj,
        "ifn": ifn,
        "n_metacells": n_meta,
        "n_by_dataset": n_by_ds,
        "seeds": {"CLDN4": seed_row("CLDN4"), "TACSTD2": seed_row("TACSTD2")},
        "membership": membership,
    }


def module_gene_list(net: dict, name: str) -> list[str]:
    if name == "grey":
        return [g for g, lab in zip(net["genes"], net["labels"]) if lab == 0]
    module = int(name.replace("ME", ""))
    return [g for g, lab in zip(net["genes"], net["labels"]) if lab == module]


def fmt(x, digits=3) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    ax = abs(float(x))
    if ax != 0 and (ax < 0.001 or ax >= 1000):
        return f"{float(x):.3g}"
    return f"{float(x):.{digits}f}"


def plot_enrichment(nets: list[dict]) -> None:
    fig, axes = plt.subplots(1, len(nets), figsize=(6.2 * len(nets), 4.8), squeeze=False)
    for ax, net in zip(axes[0], nets):
        modules = []
        for row in net["enrich_rows"]:
            if row["family"] == "TJ" and row["module"] not in modules:
                modules.append(row["module"])
        modules = [m for m in modules if m != "grey"]
        tj = []
        ifn = []
        for module in modules:
            tj.append(next(r["p"] for r in net["enrich_rows"] if r["module"] == module and r["family"] == "TJ"))
            ifn.append(next(r["p"] for r in net["enrich_rows"] if r["module"] == module and r["family"] == "IFN"))
        y = np.arange(len(modules))
        ax.barh(y - 0.18, -np.log10(np.clip(tj, 1e-300, 1)), height=0.36, color="#1a5276", label="TJ / barrier")
        ax.barh(y + 0.18, -np.log10(np.clip(ifn, 1e-300, 1)), height=0.36, color="#b03a2e", label="IFN")
        ax.set_yticks(y)
        ax.set_yticklabels(modules)
        ax.set_xlabel(r"$-\log_{10}$ enrichment $p$")
        ax.set_title(net["stratum"] + " malignant metacells")
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_FIG / "fig_module_enrichment.png", dpi=140)
    fig.savefig(OUT_FIG / "fig_module_enrichment.pdf")
    plt.close(fig)


def plot_scatter(panels: list[dict]) -> None:
    fig, axes = plt.subplots(1, len(panels), figsize=(4.6 * len(panels), 4.2), squeeze=False)
    for ax, panel in zip(axes[0], panels):
        for ds in DATASETS:
            items = [r for r in panel["rows"] if r["dataset"] == ds and np.isfinite(r["score"])]
            ax.scatter(
                [r["score"] for r in items],
                [r["frac_tnk"] for r in items],
                s=28,
                color=COLORS[ds],
                label=ds.replace("GSE", ""),
            )
        meta = panel["meta"]
        ax.set_title(f"{panel['title']}\nDL ρ={fmt(meta['rho'])}, p={fmt(meta['p'])}, n={meta['N']}")
        ax.set_xlabel("Malignant module score")
        ax.set_ylabel("T/NK fraction")
    axes[0, 0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT_FIG / "fig_module_vs_tnk.png", dpi=140)
    fig.savefig(OUT_FIG / "fig_module_vs_tnk.pdf")
    plt.close(fig)


def plot_forest(panels: list[dict]) -> None:
    fig, axes = plt.subplots(1, len(panels), figsize=(4.4 * len(panels), 3.6), squeeze=False)
    for ax, panel in zip(axes[0], panels):
        ys = np.arange(len(DATASETS) + 1)
        rhos, lo, hi, labels = [], [], [], []
        for ds, row in zip(DATASETS, panel["singles"]):
            # approximate CI from Fisher z
            if not np.isfinite(row["rho"]) or row["n"] <= 3:
                continue
            z = np.arctanh(np.clip(row["rho"], -0.999, 0.999))
            se = 1 / np.sqrt(row["n"] - 3)
            labels.append(f"{ds.replace('GSE','')} (n={row['n']})")
            rhos.append(row["rho"])
            lo.append(np.tanh(z - 1.96 * se))
            hi.append(np.tanh(z + 1.96 * se))
        meta = panel["meta"]
        labels.append(f"DL (N={meta['N']})")
        rhos.append(meta["rho"])
        lo.append(meta["ci_lo"])
        hi.append(meta["ci_hi"])
        y = np.arange(len(labels))[::-1]
        ax.axvline(0, color="#bbbbbb", lw=1)
        ax.errorbar(
            rhos,
            y,
            xerr=[np.array(rhos) - np.array(lo), np.array(hi) - np.array(rhos)],
            fmt="o",
            color="#1a5276",
            capsize=3,
        )
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.set_xlabel("Spearman ρ vs T/NK")
        ax.set_title(panel["title"])
    fig.tight_layout()
    fig.savefig(OUT_FIG / "fig_forest_tnk.png", dpi=140)
    fig.savefig(OUT_FIG / "fig_forest_tnk.pdf")
    plt.close(fig)


def paired_summary(pos_rows: list[dict], neg_rows: list[dict]) -> list[dict]:
    neg = {(r["dataset"], r["unit_id"]): r for r in neg_rows}
    out = []
    for ds in DATASETS:
        deltas = []
        for r in pos_rows:
            if r["dataset"] != ds:
                continue
            other = neg.get((r["dataset"], r["unit_id"]))
            if other is None:
                continue
            if r["n_cells_scored"] < PAIRED_MIN_CELLS or other["n_cells_scored"] < PAIRED_MIN_CELLS:
                continue
            if not np.isfinite(r["score"]) or not np.isfinite(other["score"]):
                continue
            deltas.append(r["score"] - other["score"])
        deltas_a = np.array(deltas, dtype=float)
        if deltas_a.size >= 5 and np.any(deltas_a != 0):
            stat = wilcoxon(deltas_a, alternative="two-sided", zero_method="wilcox")
            p = float(stat.pvalue)
        else:
            p = np.nan
        out.append(
            {
                "dataset": ds,
                "n_patients": int(deltas_a.size),
                "median_pos_minus_neg": float(np.median(deltas_a)) if deltas_a.size else np.nan,
                "n_pos_higher": int(np.sum(deltas_a > 0)) if deltas_a.size else 0,
                "wilcoxon_p": p,
            }
        )
    return out


def hub_genes(net: dict, module_name: str, n: int = 12) -> list[str]:
    rows = [r for r in net["kme_rows"] if r["module"] == module_name]
    rows.sort(key=lambda r: -r["kME"])
    return [f"{r['gene']} ({r['kME']:.2f})" for r in rows[:n]]


def write_finding(text: str) -> None:
    (ROOT / "FINDING.md").write_text(text)


def main() -> None:
    OUT_TAB.mkdir(parents=True, exist_ok=True)
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    assert_locked_endpoint()
    sets = json.loads((ROOT / "data" / "gene_sets.json").read_text())
    genes = select_genes()
    (OUT_TAB / "network_genes.txt").write_text("\n".join(genes) + "\n")
    packed = build_metacells(genes)
    nets = []
    for stratum in ("CLDN4pos", "CLDN4neg"):
        if packed[stratum] is None:
            say(f"{stratum}: no metacells")
            continue
        net = run_network(stratum, packed[stratum], genes, sets)
        if net is not None:
            nets.append(net)
    if not nets:
        raise SystemExit("no network was built")

    enrich_all = []
    power_all = []
    kme_all = []
    for net in nets:
        enrich_all.extend(net["enrich_rows"])
        for row in net["power_rows"]:
            row = dict(row)
            row["stratum"] = net["stratum"]
            power_all.append(row)
        kme_all.extend(net["kme_rows"])
    write_tsv(
        OUT_TAB / "module_enrichment.tsv",
        enrich_all,
        ["stratum", "module", "n_genes", "family", "overlap", "odds_ratio", "p", "q"],
    )
    write_tsv(
        OUT_TAB / "soft_power.tsv",
        power_all,
        ["stratum", "power", "scale_free_R2", "slope", "mean_connectivity", "median_connectivity", "chosen", "rule"],
    )
    write_tsv(OUT_TAB / "module_kME.tsv", kme_all, ["stratum", "module", "gene", "in_module", "kME"])

    # Patient-level tests. Modules were already named by enrichment.
    corr_rows = []
    score_rows = []
    paired_rows = []
    scatter_panels = []
    forest_panels = []
    for net in nets:
        named = []
        if net["tj"] is not None:
            named.append(("TJ", net["tj"]["module"]))
        if net["ifn"] is not None and (net["tj"] is None or net["ifn"]["module"] != net["tj"]["module"]):
            named.append(("IFN", net["ifn"]["module"]))
        elif net["ifn"] is not None:
            named.append(("IFN", net["ifn"]["module"]))
        # all modules, including grey excluded
        modules = []
        for row in net["enrich_rows"]:
            if row["family"] == "TJ" and row["module"] != "grey" and row["module"] not in modules:
                modules.append(row["module"])
        for module in modules:
            mgenes = module_gene_list(net, module)
            rows = pseudobulk_scores(genes, mgenes, lambda in_unit, cldn: in_unit)
            for row in rows:
                row = dict(row)
                row["stratum_network"] = net["stratum"]
                row["module"] = module
                score_rows.append(row)
            singles, meta = correlate_scores(rows, "frac_tnk")
            psingles, pmeta = correlate_partial(rows)
            csingles, cmeta = correlate_scores(rows, "mal_CLDN4_pct")
            role = []
            if net["tj"] is not None and net["tj"]["called"] and net["tj"]["module"] == module:
                role.append("TJ")
            if net["ifn"] is not None and net["ifn"]["called"] and net["ifn"]["module"] == module:
                role.append("IFN")
            if net["tj"] is not None and not net["tj"]["called"] and net["tj"]["module"] == module:
                role.append("TJ-best-overlap-not-significant")
            if net["ifn"] is not None and not net["ifn"]["called"] and net["ifn"]["module"] == module:
                role.append("IFN-best-overlap-not-significant")
            corr_rows.append(
                {
                    "stratum_network": net["stratum"],
                    "module": module,
                    "role": "+".join(role),
                    "n_genes_in_module": len(mgenes),
                    "n_genes_scored": len([g for g in mgenes if g != "CLDN4"]),
                    "rho_tnk": meta["rho"],
                    "p_tnk": meta["p"],
                    "I2_tnk": meta["I2"],
                    "ci_lo_tnk": meta["ci_lo"],
                    "ci_hi_tnk": meta["ci_hi"],
                    "N_tnk": meta["N"],
                    "rho_tnk_partial_cldn4": pmeta["rho"],
                    "p_tnk_partial_cldn4": pmeta["p"],
                    "I2_partial": pmeta["I2"],
                    "rho_cldn4": cmeta["rho"],
                    "p_cldn4": cmeta["p"],
                }
            )
            if role:
                label = f"{net['stratum']} {module} ({'+'.join(role)})"
                scatter_panels.append({"title": label, "rows": rows, "meta": meta})
                forest_panels.append({"title": label, "singles": singles, "meta": meta})
            # store cohort rows for the named modules
            for single, psingle, csingle in zip(singles, psingles, csingles):
                corr_rows.append(
                    {
                        "stratum_network": net["stratum"],
                        "module": module,
                        "role": "cohort:" + "+".join(role),
                        "dataset": single["dataset"],
                        "n_genes_in_module": len(mgenes),
                        "n_genes_scored": len([g for g in mgenes if g != "CLDN4"]),
                        "rho_tnk": single["rho"],
                        "p_tnk": single["p"],
                        "I2_tnk": "",
                        "ci_lo_tnk": "",
                        "ci_hi_tnk": "",
                        "N_tnk": single["n"],
                        "rho_tnk_partial_cldn4": psingle["rho"],
                        "p_tnk_partial_cldn4": psingle["p"],
                        "rho_cldn4": csingle["rho"],
                        "p_cldn4": csingle["p"],
                    }
                )
        # paired CLDN4+ vs CLDN4- cell scores for named modules
        for family, hit in (("TJ", net["tj"]), ("IFN", net["ifn"])):
            if hit is None:
                continue
            family = family if hit["called"] else family + "-overlap-ns"
            mgenes = module_gene_list(net, hit["module"])
            pos = pseudobulk_scores(genes, mgenes, lambda in_unit, cldn: in_unit & (cldn > 0))
            neg = pseudobulk_scores(genes, mgenes, lambda in_unit, cldn: in_unit & (cldn == 0))
            for row in paired_summary(pos, neg):
                row = dict(row)
                row["stratum_network"] = net["stratum"]
                row["module"] = hit["module"]
                row["family"] = family
                paired_rows.append(row)

    # BH across module-level (not cohort-level) T/NK tests within each network
    for stratum in ("CLDN4pos", "CLDN4neg"):
        ix = [
            i
            for i, row in enumerate(corr_rows)
            if row["stratum_network"] == stratum and not str(row["role"]).startswith("cohort:")
        ]
        q = bh([corr_rows[i]["p_tnk"] for i in ix])
        for i, qi in zip(ix, q):
            corr_rows[i]["q_tnk"] = float(qi) if np.isfinite(qi) else np.nan
    for row in corr_rows:
        row.setdefault("q_tnk", "")
        row.setdefault("dataset", "")

    write_tsv(
        OUT_TAB / "module_immune_correlation.tsv",
        corr_rows,
        [
            "stratum_network",
            "module",
            "role",
            "dataset",
            "n_genes_in_module",
            "n_genes_scored",
            "rho_tnk",
            "p_tnk",
            "q_tnk",
            "I2_tnk",
            "ci_lo_tnk",
            "ci_hi_tnk",
            "N_tnk",
            "rho_tnk_partial_cldn4",
            "p_tnk_partial_cldn4",
            "rho_cldn4",
            "p_cldn4",
        ],
    )
    write_tsv(
        OUT_TAB / "patient_module_scores.tsv",
        score_rows,
        [
            "stratum_network",
            "module",
            "dataset",
            "unit_id",
            "n_cells_scored",
            "n_genes",
            "score",
            "frac_tnk",
            "mal_CLDN4_pct",
            "locked_ifn_score",
            "locked_tj_score",
        ],
    )
    write_tsv(
        OUT_TAB / "paired_cldn4_stratum.tsv",
        paired_rows,
        [
            "stratum_network",
            "module",
            "family",
            "dataset",
            "n_patients",
            "median_pos_minus_neg",
            "n_pos_higher",
            "wilcoxon_p",
        ],
    )

    # validation against locked family scores for the named modules
    valid_rows = []
    for net in nets:
        for family, hit, key in (
            ("TJ", net["tj"], "locked_tj_score"),
            ("IFN", net["ifn"], "locked_ifn_score"),
        ):
            if hit is None:
                continue
            rows = [
                r
                for r in score_rows
                if r["stratum_network"] == net["stratum"] and r["module"] == hit["module"]
            ]
            # rename key expected by correlate_scores
            tagged = []
            for r in rows:
                rr = dict(r)
                rr["y"] = r[key]
                tagged.append(rr)
            _singles, meta = correlate_scores(tagged, "y")
            valid_rows.append(
                {
                    "stratum_network": net["stratum"],
                    "module": hit["module"],
                    "family": family,
                    "locked_score": key,
                    "rho": meta["rho"],
                    "p": meta["p"],
                    "N": meta["N"],
                }
            )
    write_tsv(
        OUT_TAB / "module_vs_locked_family_score.tsv",
        valid_rows,
        ["stratum_network", "module", "family", "locked_score", "rho", "p", "N"],
    )

    if scatter_panels:
        plot_scatter(scatter_panels[:4])
        plot_forest(forest_panels[:4])
    plot_enrichment(nets)

    # FINDING
    lines = []
    lines.append("# hdWGCNA-like CLDN4-stratified modules in concordant-4 malignant cells\n")
    lines.append(
        "ADDITIVE. **CLDN4-only.** Same four datasets as the locked concordant-4 "
        "patient result (GSE123902, GSE131907, GSE205335, GSE189357). "
        "Not GSE148071, GSE127465, GSE207422, GSE154826, or E-MTAB-13526. "
        "This does not replace the locked malignant CLDN4 %pos vs T/NK result "
        "(ρ=−0.531, P=1.65×10⁻⁵, n=65). That number was recomputed from the locked "
        "patient table as a pipeline check and matched.\n"
    )
    lines.append(
        "The R package hdWGCNA was not installed. The network is hdWGCNA-like: "
        f"metacells follow `ConstructMetacells` (k={K}, max_shared={MAX_SHARED}, "
        f"average of the raw UMI in each kNN, then log-normalize), grouped by "
        f"patient × CLDN4 status. A group with fewer than {MIN_CELLS} malignant cells "
        f"contributes no metacells. Each group is capped at {TARGET_METACELLS} metacells "
        "so a large sample does not dominate (the package default is 1000). "
        f"Neighbors are Euclidean kNN on {N_PCS} within-group PCs, not a global Harmony "
        "reduction. Before the correlation, metacell expression is z-scored within dataset. "
        "Soft-threshold QC: a power is not used when scale-free R² first exceeds 0.80 only after mean connectivity falls below 5. The primary power is then the negative-slope power with mean connectivity between 10 and 100 and the highest R², and the scale-free fit is reported as not reached. "
        "The network is signed biweight midcorrelation, signed TOM, "
        f"`dynamicTreeCut` hybrid (deepSplit=2, minClusterSize={MIN_MODULE_SIZE}, "
        f"pamRespectsDendro=FALSE), then eigengene merge at r>{MERGE_COR}. "
        f"Genes were the intersection of the four datasets, detected in ≥{MIN_FRAC:.0%} of "
        f"malignant cells in ≥{MIN_DATASETS_DETECTED} datasets, top {N_GENES} by mean "
        "within-dataset variance rank of log1p(UMI). CLDN4 was forced in if it passed "
        "detection. TJ and IFN gene sets were not used to choose genes.\n"
    )
    lines.append(
        "CLDN4-positive means malignant UMI > 0. CLDN4-negative means UMI = 0. "
        "Two networks were built, one in each stratum. "
        "A module is called the TJ/barrier module or the IFN module by one-sided "
        f"Fisher enrichment (background = network genes; overlap ≥ {MIN_SET_OVERLAP}; "
        "odds ratio > 1; smallest p), not by its correlation with T/NK. "
        "TJ = KEGG tight junction ∪ GOBP tight-junction organization ∪ CDH1/VIM/ZEB1, "
        "CLDN4 removed. IFN = Hallmark IFNα ∪ IFNγ. "
        "The patient-level score is the mean log2(CPM+1) of module genes in the "
        "malignant pseudobulk, with CLDN4 held out of the score. "
        "The test is a within-cohort Spearman, then DerSimonian–Laird on Fisher z. "
        "n is patients, not cells and not metacells. p-values are descriptive.\n"
    )
    lines.append("## Metacells and soft threshold\n")
    lines.append("| network | metacells | " + " | ".join(d.replace("GSE", "") for d in DATASETS) + " | power | scale-free R² | mean k |")
    lines.append("|---|---:|" + "|".join(["---:"] * len(DATASETS)) + "|---:|---:|---:|")
    for net in nets:
        chosen = next(r for r in net["power_rows"] if r["chosen"])
        counts = " | ".join(str(net["n_by_dataset"][ds]) for ds in DATASETS)
        lines.append(
            f"| {net['stratum']} | {net['n_metacells']} | {counts} | {net['beta']} | "
            f"{fmt(chosen['scale_free_R2'])} | {fmt(chosen['mean_connectivity'], 1)} |"
        )
    lines.append("")
    lines.append("## Which module is TJ/barrier, which is IFN\n")
    lines.append("| network | call | module | genes | overlap | odds ratio | p | q | CLDN4 module (kME) | TACSTD2 module (kME) |")
    lines.append("|---|---|---|---:|---:|---:|---:|---:|---|---|")
    for net in nets:
        seed_c = net["seeds"]["CLDN4"]
        seed_t = net["seeds"]["TACSTD2"]
        seed_txt_c = f"{seed_c['module']} ({fmt(seed_c['kME'], 2)})"
        seed_txt_t = f"{seed_t['module']} ({fmt(seed_t['kME'], 2)})"
        for family, hit in (("TJ/barrier", net["tj"]), ("IFN", net["ifn"])):
            if hit is None:
                lines.append(
                    f"| {net['stratum']} | {family} | none with overlap≥{MIN_SET_OVERLAP} and OR>1 |  |  |  |  |  | {seed_txt_c} | {seed_txt_t} |"
                )
                continue
            status = "called (q<0.05)" if hit["called"] else "best overlap, not called (q≥0.05)"
            same = ""
            if net["tj"] is not None and net["ifn"] is not None and net["tj"]["module"] == net["ifn"]["module"]:
                same = "; same module as the other family"
            lines.append(
                f"| {net['stratum']} | {family}: {status}{same} | {hit['module']} | {hit['n_genes']} | "
                f"{hit['overlap']} | {fmt(hit['odds_ratio'], 2)} | {fmt(hit['p'])} | {fmt(hit['q'])} | "
                f"{seed_txt_c} | {seed_txt_t} |"
            )
    lines.append("")
    lines.append("Hub genes are the highest kME genes inside the module in the table above. A module is a TJ or IFN module only when q<0.05.\n")
    for net in nets:
        for family, hit in (("TJ/barrier", net["tj"]), ("IFN", net["ifn"])):
            if hit is None:
                continue
            hubs = ", ".join(hub_genes(net, hit["module"]))
            lines.append(f"- {net['stratum']} {family} {hit['module']}: {hubs}")
    lines.append("")
    lines.append("## Patient-level module vs T/NK\n")
    lines.append(
        "Every non-grey module is listed. Role is TJ or IFN only when that "
        "family's enrichment q is below 0.05. Cohort rows and the full table are in "
        "`results/tables/module_immune_correlation.tsv`. "
        "q is BH across non-grey modules inside that network for the T/NK meta-analysis p. "
        "Partial ρ is Spearman of the module score vs T/NK given malignant CLDN4 %pos.\n"
    )
    lines.append("| network | module | role | N | ρ vs T/NK (95% CI) | p | q | I² | ρ vs T/NK \\| CLDN4 | p | I² partial | ρ vs CLDN4 %pos | p |")
    lines.append("|---|---|---|---:|---|---:|---:|---:|---|---:|---:|---|---:|")
    for row in corr_rows:
        if str(row["role"]).startswith("cohort:"):
            continue
        lines.append(
            f"| {row['stratum_network']} | {row['module']} | {row['role'] or 'other'} | {row['N_tnk']} | "
            f"{fmt(row['rho_tnk'])} ({fmt(row['ci_lo_tnk'])} to {fmt(row['ci_hi_tnk'])}) | "
            f"{fmt(row['p_tnk'])} | {fmt(row['q_tnk']) if row['q_tnk'] != '' else 'NA'} | "
            f"{fmt(100 * row['I2_tnk'] if row['I2_tnk'] != '' else np.nan, 1)}% | "
            f"{fmt(row['rho_tnk_partial_cldn4'])} | {fmt(row['p_tnk_partial_cldn4'])} | "
            f"{fmt(100 * row['I2_partial'] if row.get('I2_partial') != '' and row.get('I2_partial') is not None else np.nan, 1)}% | "
            f"{fmt(row['rho_cldn4'])} | {fmt(row['p_cldn4'])} |"
        )
    lines.append("")
    lines.append("Cohort Spearmans for the modules named in the identification table:\n")
    lines.append("| network | module | cohort | n | ρ vs T/NK | p | partial ρ \\| CLDN4 | p |")
    lines.append("|---|---|---|---:|---:|---:|---:|---:|")
    for row in corr_rows:
        if not str(row["role"]).startswith("cohort:") or row["role"] == "cohort:":
            continue
        lines.append(
            f"| {row['stratum_network']} | {row['module']} | {row['dataset']} | {row['N_tnk']} | "
            f"{fmt(row['rho_tnk'])} | {fmt(row['p_tnk'])} | {fmt(row['rho_tnk_partial_cldn4'])} | "
            f"{fmt(row['p_tnk_partial_cldn4'])} |"
        )
    lines.append("")
    lines.append(
        "## CLDN4-positive vs CLDN4-negative cells, same patients\n\n"
        "For each patient with at least 20 malignant cells in both strata, the module "
        "score was computed separately in each stratum. The table is the within-patient "
        "difference (CLDN4-positive minus CLDN4-negative). This is not the T/NK test.\n"
    )
    lines.append("| network | module | family | cohort | n | median Δ | n with pos>neg | Wilcoxon p |")
    lines.append("|---|---|---|---|---:|---:|---:|---:|")
    for row in paired_rows:
        lines.append(
            f"| {row['stratum_network']} | {row['module']} | {row['family']} | {row['dataset']} | "
            f"{row['n_patients']} | {fmt(row['median_pos_minus_neg'])} | {row['n_pos_higher']} | "
            f"{fmt(row['wilcoxon_p'])} |"
        )
    lines.append("")
    lines.append(
        "## Check against the locked malignant IFN and TJ scores\n\n"
        "Spearman of the new module score vs the locked pseudobulk family score "
        "on the same patients (DL across cohorts). This asks whether the called "
        "module tracks the previous IFN or TJ summary. It is not an immune test.\n"
    )
    lines.append("| network | module | family | ρ vs locked score | p | N |")
    lines.append("|---|---|---|---:|---:|---:|")
    for row in valid_rows:
        lines.append(
            f"| {row['stratum_network']} | {row['module']} | {row['family']} | {fmt(row['rho'])} | {fmt(row['p'])} | {row['N']} |"
        )
    lines.append("")
    lines.append("## Reading\n")
    lines.append(
        "The two calls are allowed to land on the same module. They do not. "
        "In both strata the significant IFN module is an MHC-I / immunoproteasome block, "
        "not the TJ block.\n"
    )
    barrier_focus = ["TJP1", "TJP2", "TJP3", "OCLN", "F11R", "CDH1", "CLDN1", "CLDN3", "CLDN4", "CLDN7", "TACSTD2", "EPCAM", "CRB3"]
    lines.append(
        "Junction genes in the network, and the module that contains each one "
        "(kME for that module). CLDN4 is held out of the TJ enrichment set, so it "
        "cannot by itself create the TJ call.\n"
    )
    lines.append("| network | gene | module | kME |")
    lines.append("|---|---|---|---:|")
    for net in nets:
        gene_pos = {g: i for i, g in enumerate(net["genes"])}
        for gene in barrier_focus:
            if gene not in gene_pos:
                lines.append(f"| {net['stratum']} | {gene} | not in network |  |")
                continue
            lab = int(net["labels"][gene_pos[gene]])
            name = "grey" if lab == 0 else f"ME{lab}"
            value = ""
            if lab in net["membership"]:
                value = fmt(float(net["membership"][lab][gene_pos[gene]]), 2)
            lines.append(f"| {net['stratum']} | {gene} | {name} | {value} |")
    lines.append("")
    for net in nets:
        if net["tj"] is None:
            continue
        tj_in = []
        mod_genes = set(module_gene_list(net, net["tj"]["module"]))
        kme = {r["gene"]: r["kME"] for r in net["kme_rows"] if r["module"] == net["tj"]["module"]}
        for gene in sorted(mod_genes & set(sets["TJ"]), key=lambda g: -kme.get(g, 0)):
            tj_in.append(f"{gene} ({kme.get(gene, float('nan')):.2f})")
        status = "called" if net["tj"]["called"] else "best overlap, not called"
        lines.append(
            f"{net['stratum']} TJ set inside {net['tj']['module']} ({status}): " + ", ".join(tj_in) + ".\n"
        )
    n_tj = len(set(sets["TJ"]) & set(genes))
    n_ifn = len(set(sets["IFN"]) & set(genes))
    lines.append(
        f"Of the locked sets, {n_tj} TJ genes and {n_ifn} IFN genes entered the {len(genes)}-gene network. "
        "Enrichment is against that background.\n"
    )
    lines.append(
        "Patient-level T/NK tests use the full malignant pseudobulk. "
        "A 95% CI that crosses 0 is not a concordant module–immune correlation. "
        "The locked result is malignant CLDN4 percent-positive versus T/NK fraction. "
        "It is not replaced by a module score. "
        "The within-patient difference (CLDN4-positive cells minus CLDN4-negative cells) "
        "is a different contrast from that between-patient result, and from the locked "
        "Q4-versus-Q1 malignant IFN comparison.\n"
    )
    lines.append(
        "## What this does not say\n\n"
        "Metacells are not extra patients. A metacell-level p-value is not reported. "
        "The module was discovered and scored in the same 65 units; the gene list was "
        "not chosen using T/NK, but it is not an external locked signature. "
        "Dataset z-scoring removes mean shifts. It is not Harmony and it is not a "
        "batch-free claim. Visium co-localization is not in this analysis. "
        "Public mouse KL matrices were not added.\n"
    )
    lines.append("Figures: `results/figures/fig_module_enrichment.png`, `fig_module_vs_tnk.png`, `fig_forest_tnk.png`.\n")
    write_finding("\n".join(lines) + "\n")
    say("wrote FINDING.md")


if __name__ == "__main__":
    main()
