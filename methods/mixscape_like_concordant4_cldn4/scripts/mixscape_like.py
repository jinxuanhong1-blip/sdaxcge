#!/usr/bin/env python3
"""Mixscape-like local perturbation signatures for concordant-4 malignant cells.

This is not a CRISPR screen and it does not run Seurat RunMixscape.
There are no gRNAs, no non-targeting guides, and no escapee GMM.

Observational analog of CalcPerturbSig (Papalexi et al.; Seurat):
  * NT-like pool = within-unit CLDN4 value-quartile Q4 (CLDN4-high malignant).
  * KD-like / query cells = lower CLDN4 quartiles, especially Q1.
  * For each malignant cell, subtract the mean log-normalized expression of its
    k nearest NT-like neighbors in PCA space (neighbors found inside the unit).
  * The residual is the local perturbation signature.
  * That signature is scored on Hallmark IFN (alpha union gamma), with a
    detection-matched gene set and Hallmark spermatogenesis as specificity controls.

PCA features exclude CLDN4 and every gene that is later used as a readout
(IFN, MHC-I/APM, spermatogenesis, and the detection-matched controls) so the
neighbor graph is not built on the genes being scored.

Honest unit = patient / donor / sample. Cell counts are not n.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
import traceback
from array import array
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import mmread
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[1]
GEO = Path("/tmp/geo_mixscape")
K_TARGET = 20
N_PCS = 15
MIN_K = 10
MIN_NT = 15  # Q4 cells; k = min(K_TARGET, n_Q4 - 1) must be >= MIN_K
MIN_QUERY = 20
MIN_HVG = 100
N_HVG = 2000
N_PERM = 499
QC_MIN_GENES = 200
QC_MIN_UMI = 500
QC_MAX_MT = 20.0
SEED = 1


def say(msg: str) -> None:
    print(msg, flush=True)


def load_gene_sets(path: Path) -> dict[str, list[str]]:
    raw = json.loads(path.read_text())["sets"]
    ifn_a = [g.upper() for g in raw["HALLMARK_INTERFERON_ALPHA_RESPONSE"]]
    ifn_g = [g.upper() for g in raw["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]
    ifn = list(dict.fromkeys(ifn_a + ifn_g))
    ifn_set = set(ifn)
    sperm = [g.upper() for g in raw["HALLMARK_SPERMATOGENESIS"] if g.upper() not in ifn_set]
    mhc = [g.upper() for g in raw["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"] if g.upper() not in ifn_set]
    return {
        "IFN": ifn,
        "IFNA": ifn_a,
        "IFNG": ifn_g,
        "SPERM": sperm,
        "MHC": mhc,
    }


def load_locked(path: Path) -> dict[str, list[dict]]:
    by: dict[str, list[dict]] = {}
    with path.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            by.setdefault(row["dataset"], []).append(row)
    return by


def upper_index(genes: list[str]) -> tuple[list[str], dict[str, int]]:
    """Keep the first column of each upper-cased symbol."""
    keep_pos = []
    index: dict[str, int] = {}
    out = []
    for i, g in enumerate(genes):
        u = g.upper()
        if u in index:
            continue
        index[u] = len(out)
        out.append(u)
        keep_pos.append(i)
    return out, index, keep_pos


def value_quartiles(x: np.ndarray) -> tuple[np.ndarray, tuple[float, float, float]]:
    """Value quartiles. Ties stay in the lower bin, so empty bins are allowed.

    Q1 is x <= q25, Q2 is q25 < x <= q50, Q3 is q50 < x <= q75, Q4 is x > q75.
    When CLDN4 is zero-inflated the lower breaks collapse and Q2/Q3 can be empty.
    Cells are not randomly split across tied zeros.
    """
    b1, b2, b3 = [float(v) for v in np.quantile(x, [0.25, 0.5, 0.75])]
    lab = np.full(len(x), "Q4", dtype=object)
    lab[x <= b3] = "Q3"
    lab[x <= b2] = "Q2"
    lab[x <= b1] = "Q1"
    return lab, (b1, b2, b3)


def qc_mask(X: sparse.csr_matrix, genes: list[str]) -> np.ndarray:
    ncount = np.asarray(X.sum(axis=0)).ravel()
    nfeat = np.diff(X.tocsc().indptr).astype(np.float64) if False else np.asarray((X > 0).sum(axis=0)).ravel()
    mt = np.array([g.startswith("MT-") for g in genes])
    if mt.any():
        mt_count = np.asarray(X[mt].sum(axis=0)).ravel()
        mt_pct = 100.0 * mt_count / np.maximum(ncount, 1.0)
    else:
        mt_pct = np.zeros_like(ncount, dtype=float)
    return (nfeat >= QC_MIN_GENES) & (ncount >= QC_MIN_UMI) & (mt_pct < QC_MAX_MT)


def row_vector(X: sparse.csr_matrix, idx: int | None, n: int) -> np.ndarray:
    if idx is None:
        return np.zeros(n, dtype=np.float64)
    return np.asarray(X[idx].todense()).ravel().astype(np.float64)


def lognorm_rows(X: sparse.csr_matrix, row_idx: np.ndarray, lib: np.ndarray) -> np.ndarray:
    """Dense log1p(CP10k) for selected genes. Shape (n_genes, n_cells)."""
    if len(row_idx) == 0:
        return np.zeros((0, X.shape[1]), dtype=np.float32)
    sub = X[row_idx].tocsr()
    out = np.zeros((len(row_idx), X.shape[1]), dtype=np.float32)
    scale = (10000.0 / lib).astype(np.float64)
    indptr, indices, data = sub.indptr, sub.indices, sub.data
    for i in range(sub.shape[0]):
        s, e = indptr[i], indptr[i + 1]
        if s == e:
            continue
        cols = indices[s:e]
        out[i, cols] = np.log1p(data[s:e] * scale[cols]).astype(np.float32)
    return out


def gene_variances(X: sparse.csr_matrix, lib: np.ndarray) -> np.ndarray:
    n = X.shape[1]
    sum_log = np.zeros(X.shape[0], dtype=np.float64)
    sum_sq = np.zeros(X.shape[0], dtype=np.float64)
    scale = 10000.0 / lib
    indptr, indices, data = X.indptr, X.indices, X.data
    for i in range(X.shape[0]):
        s, e = indptr[i], indptr[i + 1]
        if s == e:
            continue
        cols = indices[s:e]
        ln = np.log1p(data[s:e] * scale[cols])
        sum_log[i] = ln.sum()
        sum_sq[i] = np.dot(ln, ln)
    mean = sum_log / n
    var = sum_sq / n - mean ** 2
    var[var < 0] = 0
    return var


def match_by_detection(
    det: np.ndarray,
    targets: list[int],
    pool: list[int],
) -> list[int]:
    """Greedy nearest detection-rate match without replacement. Deterministic."""
    if not targets or not pool:
        return []
    remaining = list(pool)
    det_rem = det[np.array(remaining)]
    chosen = []
    order = sorted(range(len(targets)), key=lambda i: det[targets[i]])
    for ti in order:
        if not remaining:
            break
        d = det[targets[ti]]
        j = int(np.argmin(np.abs(det_rem - d)))
        chosen.append(remaining[j])
        remaining.pop(j)
        det_rem = np.delete(det_rem, j)
    return chosen


def nearest_nt(emb: np.ndarray, nt_idx: np.ndarray, k: int) -> np.ndarray:
    """k nearest NT-like cells for every cell. Self is excluded."""
    nn = NearestNeighbors(n_neighbors=min(k + 1, len(nt_idx)), metric="euclidean")
    nn.fit(emb[nt_idx])
    ind = nn.kneighbors(emb, return_distance=False)
    mapped = nt_idx[ind]
    n = emb.shape[0]
    out = np.empty((n, k), dtype=np.int32)
    ar = np.arange(n)
    for i in range(n):
        row = mapped[i]
        row = row[row != ar[i]]
        if len(row) < k:
            raise RuntimeError(f"cell {i} has {len(row)} NT neighbors, need {k}")
        out[i] = row[:k]
    return out


def mean_gene_residual(G: np.ndarray, nn: np.ndarray, cells: np.ndarray | None = None) -> np.ndarray:
    """Mean over `cells` of (expression - neighbor mean). Returns one value per gene.

    G is genes x cells, float. Computed in gene chunks.
    """
    if cells is None:
        cells = np.arange(G.shape[1])
    acc = np.zeros(G.shape[0], dtype=np.float64)
    nn_c = nn[cells]
    step = 250
    for s in range(0, G.shape[0], step):
        e = min(G.shape[0], s + step)
        sub = G[s:e]
        neigh = sub[:, nn_c].mean(axis=2)
        acc[s:e] = (sub[:, cells] - neigh).mean(axis=1)
    return acc


def cell_mean_residual(G: np.ndarray, nn: np.ndarray) -> np.ndarray:
    """Per-cell mean residual across genes. G is genes x cells."""
    if G.shape[0] == 0:
        return np.full(G.shape[1], np.nan)
    acc = np.zeros(G.shape[1], dtype=np.float64)
    step = 200
    for s in range(0, G.shape[0], step):
        e = min(G.shape[0], s + step)
        sub = G[s:e]
        neigh = sub[:, nn].mean(axis=2)
        acc += (sub - neigh).sum(axis=0)
    return acc / G.shape[0]


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na == 0 or nb == 0:
        return float("nan")
    return float(np.dot(a, b) / (na * nb))


def cosine_perm(sig: np.ndarray, mask: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    ind = mask.astype(np.float64)
    obs = cosine(sig, ind)
    if not np.isfinite(obs):
        return obs, float("nan")
    null = np.empty(N_PERM, dtype=np.float64)
    for i in range(N_PERM):
        null[i] = cosine(rng.permutation(sig), ind)
    # one-sided: signature aligned with IFN more than a gene-label shuffle
    p = (1.0 + np.sum(null >= obs)) / (1.0 + N_PERM)
    return obs, float(p)


def present_idx(index: dict[str, int], names: list[str]) -> list[int]:
    return [index[g] for g in names if g in index]


def analyze_unit(
    X: sparse.csr_matrix,
    genes: list[str],
    sets: dict[str, list[str]],
    rng: np.random.Generator,
) -> tuple[dict, list[dict]]:
    """Return (summary, cell_rows). X is genes x QC-passed malignant cells."""
    n = X.shape[1]
    index = {g: i for i, g in enumerate(genes)}
    summary: dict = {"n_qc": n, "status": "ok"}
    if "CLDN4" not in index:
        summary["status"] = "no_cldn4"
        return summary, []
    cldn4 = row_vector(X, index["CLDN4"], n)
    if np.unique(cldn4).size < 2:
        summary["status"] = "cldn4_constant"
        return summary, []
    labels, breaks = value_quartiles(cldn4)
    counts = {q: int(np.sum(labels == q)) for q in ("Q1", "Q2", "Q3", "Q4")}
    summary.update(
        {
            "n_Q1": counts["Q1"],
            "n_Q2": counts["Q2"],
            "n_Q3": counts["Q3"],
            "n_Q4": counts["Q4"],
            "q25": breaks[0],
            "q50": breaks[1],
            "q75": breaks[2],
            "cldn4_pct": float(100.0 * np.mean(cldn4 > 0)),
            "cldn4_mean_log1p": float(np.mean(np.log1p(cldn4))),
            "frac_zero_cldn4": float(np.mean(cldn4 == 0)),
        }
    )
    n_nt = counts["Q4"]
    if n_nt < MIN_NT or (n_nt - 1) < MIN_K or n < (MIN_QUERY + MIN_NT):
        summary["status"] = "too_few_nt"
        return summary, []
    k = min(K_TARGET, n_nt - 1)
    summary["k"] = k

    lib = np.asarray(X.sum(axis=0)).ravel().astype(np.float64)
    lib[lib <= 0] = 1.0
    # Detection is per gene (fraction of cells), not per cell.
    det = np.asarray((X > 0).sum(axis=1)).ravel().astype(np.float64) / n
    if det.shape[0] != X.shape[0]:
        raise RuntimeError(f"detection vector {det.shape[0]} != n_genes {X.shape[0]}")

    ifn_i = present_idx(index, sets["IFN"])
    ifna_i = present_idx(index, sets["IFNA"])
    ifng_i = present_idx(index, sets["IFNG"])
    mhc_i = present_idx(index, sets["MHC"])
    sperm_i = present_idx(index, sets["SPERM"])
    summary["n_ifn"] = len(ifn_i)
    summary["n_ifna"] = len(ifna_i)
    summary["n_ifng"] = len(ifng_i)
    summary["n_mhc"] = len(mhc_i)
    summary["n_sperm"] = len(sperm_i)
    if len(ifn_i) < 5:
        summary["status"] = "too_few_ifn_genes"
        return summary, []

    blocked = set(ifn_i) | set(mhc_i) | set(sperm_i) | {index["CLDN4"]}
    pool = [
        i
        for i, g in enumerate(genes)
        if i not in blocked and not g.startswith("MT-") and det[i] >= 0.01
    ]
    matched_ifn = match_by_detection(det, ifn_i, pool)
    pool2 = [i for i in pool if i not in set(matched_ifn)]
    matched_sperm = match_by_detection(det, sperm_i, pool2) if sperm_i else []
    summary["n_matched_ifn"] = len(matched_ifn)
    summary["n_matched_sperm"] = len(matched_sperm)

    held = blocked | set(matched_ifn) | set(matched_sperm)
    var = gene_variances(X, lib)
    hvg_candidates = [
        i
        for i, g in enumerate(genes)
        if i not in held and not g.startswith("MT-") and det[i] >= 0.05 and var[i] > 0
    ]
    hvg_candidates.sort(key=lambda i: var[i], reverse=True)
    hvg = hvg_candidates[:N_HVG]
    summary["n_hvg"] = len(hvg)
    if len(hvg) < MIN_HVG:
        summary["status"] = "too_few_hvg"
        return summary, []

    # PCA on held-out-cleaned HVGs. Cells x genes, z-scored, clipped.
    H = lognorm_rows(X, np.array(hvg, dtype=int), lib).T.astype(np.float64)  # cells x genes
    sd = H.std(axis=0, ddof=1)
    keep_sd = sd > 1e-8
    H = H[:, keep_sd]
    sd = sd[keep_sd]
    H = (H - H.mean(axis=0)) / sd
    np.clip(H, -10, 10, out=H)
    n_pcs = int(min(N_PCS, H.shape[0] - 2, H.shape[1]))
    if n_pcs < 5:
        summary["status"] = "too_few_pcs"
        return summary, []
    summary["ndims"] = n_pcs
    pca = PCA(n_components=n_pcs, svd_solver="full", random_state=SEED)
    emb = pca.fit_transform(H)

    nt_idx = np.where(labels == "Q4")[0]
    nn = nearest_nt(emb, nt_idx, k)

    def score(idxs: list[int]) -> np.ndarray:
        if not idxs:
            return np.full(n, np.nan)
        G = lognorm_rows(X, np.array(idxs, dtype=int), lib)
        return cell_mean_residual(G, nn)

    ifn = score(ifn_i)
    ifna = score(ifna_i)
    ifng = score(ifng_i)
    mhc = score(mhc_i)
    sperm = score(sperm_i)
    matched = score(matched_ifn)
    matched_sp = score(matched_sperm) if matched_sperm else np.full(n, np.nan)

    # Unmatched (no neighbor subtraction) IFN mean, for a side-by-side estimand.
    Gifn = lognorm_rows(X, np.array(ifn_i, dtype=int), lib)
    unmatched = Gifn.mean(axis=0)

    # KD-like = lowest occupied quartile that is not Q4 and has >= MIN_QUERY cells,
    # preferring Q1.
    kd_q = None
    for q in ("Q1", "Q2", "Q3"):
        if counts[q] >= MIN_QUERY:
            kd_q = q
            break
    summary["kd_quartile"] = kd_q or ""
    if kd_q is None:
        summary["status"] = "too_few_query"
        return summary, []

    kd_cells = np.where(labels == kd_q)[0]
    non_nt = np.where(labels != "Q4")[0]
    summary["n_non_nt"] = int(len(non_nt))
    summary["n_kd"] = int(len(kd_cells))

    if np.unique(cldn4[non_nt]).size >= 2 and len(non_nt) >= MIN_QUERY:
        rho, rho_p = spearmanr(cldn4[non_nt], ifn[non_nt])
        summary["spearman_cldn4_ifn"] = float(rho)
        summary["spearman_p_cells"] = float(rho_p)
    else:
        summary["spearman_cldn4_ifn"] = float("nan")
        summary["spearman_p_cells"] = float("nan")

    # Gene-level KD signature on HVG ∪ IFN, cosine to the IFN indicator.
    universe = list(dict.fromkeys(list(hvg) + ifn_i))
    Gu = lognorm_rows(X, np.array(universe, dtype=int), lib)
    sig = mean_gene_residual(Gu, nn, kd_cells)
    ifn_mask = np.array([universe[i] in set(ifn_i) for i in range(len(universe))])
    cos, cos_p = cosine_perm(sig, ifn_mask, rng)
    summary["cosine_ifn"] = cos
    summary["cosine_ifn_p"] = cos_p
    summary["n_universe"] = len(universe)

    def qmean(v: np.ndarray, q: str) -> float:
        m = labels == q
        if not np.any(m) or not np.isfinite(v[m]).any():
            return float("nan")
        return float(np.nanmean(v[m]))

    for q in ("Q1", "Q2", "Q3", "Q4"):
        summary[f"ifn_{q}"] = qmean(ifn, q)
        summary[f"ifna_{q}"] = qmean(ifna, q)
        summary[f"ifng_{q}"] = qmean(ifng, q)
        summary[f"mhc_{q}"] = qmean(mhc, q)
        summary[f"sperm_{q}"] = qmean(sperm, q)
        summary[f"matched_{q}"] = qmean(matched, q)
        summary[f"delta_{q}"] = qmean(ifn - matched, q)
        summary[f"sperm_delta_{q}"] = qmean(sperm - matched_sp, q) if matched_sperm else float("nan")
        summary[f"unmatched_ifn_{q}"] = qmean(unmatched, q)

    q1 = labels == "Q1"
    q4 = labels == "Q4"
    if np.any(q1) and np.any(q4):
        summary["unmatched_ifn_q1_minus_q4"] = float(np.mean(unmatched[q1]) - np.mean(unmatched[q4]))
    else:
        summary["unmatched_ifn_q1_minus_q4"] = float("nan")

    cell_rows = []
    for i in range(n):
        cell_rows.append(
            {
                "quartile": labels[i],
                "cldn4_count": float(cldn4[i]),
                "ifn_prtb": float(ifn[i]),
                "ifna_prtb": float(ifna[i]),
                "ifng_prtb": float(ifng[i]),
                "mhc_prtb": float(mhc[i]),
                "sperm_prtb": float(sperm[i]),
                "matched_prtb": float(matched[i]),
                "delta_prtb": float(ifn[i] - matched[i]),
                "sperm_delta_prtb": float(sperm[i] - matched_sp[i]) if matched_sperm else "",
            }
        )
    return summary, cell_rows


def simulate() -> None:
    rng = np.random.default_rng(SEED)
    # More genes than cells, so a per-cell vector cannot be indexed as genes.
    n = 200
    n_genes = 350
    X = rng.poisson(1.0, size=(n_genes, n)).astype(np.float64)
    cldn4 = np.repeat([0, 2, 8, 20], 50).astype(np.float64)
    ifn_level = np.repeat([12, 6, 2, 0], 50).astype(np.float64)
    X[0] = cldn4
    for g in range(1, 31):
        X[g] = ifn_level
    genes = [f"G{i}" for i in range(n_genes)]
    genes[0] = "CLDN4"
    sets = {
        "IFN": [f"G{i}" for i in range(1, 31)],
        "IFNA": [f"G{i}" for i in range(1, 16)],
        "IFNG": [f"G{i}" for i in range(16, 31)],
        "MHC": ["G31", "G32"],
        "SPERM": [f"G{i}" for i in range(33, 48)],
    }
    summary, rows = analyze_unit(sparse.csr_matrix(X), genes, sets, np.random.default_rng(SEED))
    say(f"SIM status={summary.get('status')} spearman={summary.get('spearman_cldn4_ifn')} "
        f"IFN Q1={summary.get('ifn_Q1')} Q2={summary.get('ifn_Q2')} Q3={summary.get('ifn_Q3')} Q4={summary.get('ifn_Q4')} "
        f"deltaQ1={summary.get('delta_Q1')} spermQ1={summary.get('sperm_Q1')} cosine={summary.get('cosine_ifn')}")
    if summary.get("status") != "ok":
        raise SystemExit(f"simulation did not run: {summary.get('status')}")
    rho = summary["spearman_cldn4_ifn"]
    if not (rho < -0.3):
        raise SystemExit(f"simulation Spearman not negative: {rho}")
    if not (summary["ifn_Q1"] > summary["ifn_Q3"] > summary["ifn_Q4"]):
        raise SystemExit("simulation quartile order failed")
    if not (summary["delta_Q1"] > 0):
        raise SystemExit("simulation matched IFN delta not positive")
    if not (abs(summary["sperm_Q1"]) < abs(summary["ifn_Q1"]) * 0.5):
        raise SystemExit("simulation spermatogenesis control is not near zero relative to IFN")
    say(f"SIM ok ({len(rows)} cells)")


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def dedup_csr(X: sparse.spmatrix, genes: list[str]) -> tuple[sparse.csr_matrix, list[str]]:
    genes_u, index, keep = upper_index(genes)
    if len(keep) != len(genes):
        X = X.tocsr()[keep]
    else:
        X = X.tocsr()
    X.data = X.data.astype(np.float64)
    return X, genes_u


def marker_malignant(X: sparse.csr_matrix, index: dict[str, int]) -> np.ndarray:
    n = X.shape[1]
    def g(name: str) -> np.ndarray:
        return row_vector(X, index.get(name), n)
    mal = ((g("EPCAM") > 0) | (g("KRT8") > 0) | (g("KRT18") > 0) | (g("KRT19") > 0)) & (g("PTPRC") == 0)
    return mal


def prepare_malignant(X: sparse.csr_matrix, genes: list[str], mal: np.ndarray) -> tuple[sparse.csr_matrix, list[str], int]:
    n_raw = int(mal.sum())
    if n_raw == 0:
        return X[:, []], genes, 0
    sub = X[:, mal].tocsr()
    keep = qc_mask(sub, genes)
    return sub[:, keep].tocsr(), genes, n_raw


def select_gse123902(marker_tsv: Path, locked_ids: set[str]) -> list[dict]:
    rows = list(csv.DictReader(marker_tsv.open(), delimiter="\t"))
    rows = [r for r in rows if r["tissue"] in ("PRIMARY", "METASTASIS")]
    rows = sorted(rows, key=lambda r: (r["patient"], r["tissue"]))
    out = []
    seen = set()
    for r in rows:
        if r["patient"] in seen:
            continue
        seen.add(r["patient"])
        if int(float(r["n_malignant"])) >= 20 and int(float(r["n_tnk"])) >= 20:
            out.append(r)
    got = {r["patient"] for r in out}
    if got != locked_ids:
        raise SystemExit(f"GSE123902 donors {sorted(got)} != locked {sorted(locked_ids)}")
    return out


def run_gse123902(locked: list[dict], sets, rng, geo: Path, root: Path) -> list[tuple[dict, list[dict]]]:
    chosen = select_gse123902(root / "data" / "GSE123902_marker_units.tsv", {r["unit_id"] for r in locked})
    by_id = {r["patient"]: r for r in chosen}
    results = []
    d = geo / "gse123902"
    for unit in locked:
        uid = unit["unit_id"]
        meta = by_id[uid]
        path = d / meta["file"]
        say(f"GSE123902 {uid} {meta['file']}")
        df_genes = None
        import pandas as pd
        df = pd.read_csv(path, index_col=0)
        genes = [str(c) for c in df.columns]
        X = sparse.csr_matrix(df.to_numpy(dtype=np.float64).T)
        del df
        X, genes = dedup_csr(X, genes)
        index = {g: i for i, g in enumerate(genes)}
        mal = marker_malignant(X, index)
        sub, genes, n_raw = prepare_malignant(X, genes, mal)
        del X
        say(f"  malignant {n_raw} QC {sub.shape[1]} locked {unit['n_malignant']}")
        rec, cells = _finish_unit("GSE123902", unit, sub, genes, sets, rng, n_raw)
        results.append((rec, cells))
    return results


def _finish_unit(dataset, unit, sub, genes, sets, rng, n_raw) -> tuple[dict, list[dict]]:
    base = {
        "dataset": dataset,
        "unit_id": unit["unit_id"],
        "unit_type": unit["unit_type"],
        "tissue": unit["tissue"],
        "n_malignant_locked": unit["n_malignant"],
        "n_malignant_raw": n_raw,
    }
    try:
        if sub.shape[1] == 0:
            summary = {"status": "no_malignant_qc", "n_qc": 0}
            cells = []
        else:
            summary, cells = analyze_unit(sub, genes, sets, rng)
    except Exception:
        traceback.print_exc()
        summary = {"status": "error", "n_qc": int(sub.shape[1])}
        cells = []
    summary.update(base)
    for c in cells:
        c["dataset"] = dataset
        c["unit_id"] = unit["unit_id"]
    say(f"  status={summary.get('status')} k={summary.get('k')} "
        f"rho={summary.get('spearman_cldn4_ifn')} IFN Q1={summary.get('ifn_Q1')} Q4={summary.get('ifn_Q4')}")
    return summary, cells


def run_gse189357(locked, sets, rng, geo: Path) -> list[tuple[dict, list[dict]]]:
    d = geo / "gse189357"
    results = []
    for unit in locked:
        uid = unit["unit_id"]
        mtx = list(d.glob(f"*_{uid}_matrix.mtx.gz"))
        if len(mtx) != 1:
            raise SystemExit(f"GSE189357 {uid} matrix files: {mtx}")
        prefix = str(mtx[0]).replace("_matrix.mtx.gz", "")
        say(f"GSE189357 {uid}")
        feat_path = prefix + "_features.tsv.gz"
        genes = []
        with gzip.open(feat_path, "rt") as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                sym = p[1] if len(p) > 1 and p[1] else p[0]
                genes.append(sym)
        X = mmread(mtx[0]).tocsr()
        X, genes = dedup_csr(X, genes)
        index = {g: i for i, g in enumerate(genes)}
        mal = marker_malignant(X, index)
        sub, genes, n_raw = prepare_malignant(X, genes, mal)
        del X
        say(f"  malignant {n_raw} QC {sub.shape[1]} locked {unit['n_malignant']}")
        results.append(_finish_unit("GSE189357", unit, sub, genes, sets, rng, n_raw))
    return results


def parse_sample_map(path: Path) -> dict[str, dict]:
    """Map normalized sample code -> patient / tissue.

    Parsed from the public GEO family soft (sample title + characteristics).
    The soft file itself is not committed; this table is the join key.
    """
    out = {}
    with path.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            out[row["code"]] = {"patient": row["patient"], "tissue": row["tissue"], "title": row.get("title", "")}
    return out


def norm_code(x: str) -> str:
    x = x.upper().replace("-", "_")
    for suf in ("_5P", "_3P"):
        if x.endswith(suf):
            x = x[: -len(suf)]
    return x


def export_gse205335_via_r(geo: Path, keep_tsv: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    r = out_dir / "export.R"
    r.write_text(
        f"""
suppressPackageStartupMessages(library(Matrix))
src <- "{geo / "GSE205335_Lung_IO_UMI_matrix.rds.gz"}"
plain <- "{geo / "GSE205335_plain.rds"}"
if (!file.exists(plain)) {{
  message("double-gunzip")
  system(sprintf("gzip -dc %s | gzip -dc > %s", src, plain))
}}
message("readRDS")
mat <- readRDS(plain)
message("class ", paste(class(mat), collapse=","), " dim ", paste(dim(mat), collapse=" x "))
if (!inherits(mat, "dgCMatrix")) mat <- as(as(mat, "CsparseMatrix"), "dgCMatrix")
rn <- toupper(rownames(mat))
keepg <- !duplicated(rn)
mat <- mat[keepg, , drop = FALSE]
rownames(mat) <- rn[keepg]
message("CLDN4 ", "CLDN4" %in% rownames(mat), " col1 ", colnames(mat)[1])
id <- read.delim("{keep_tsv}", stringsAsFactors = FALSE)
message("keep rows ", nrow(id), " example ", id$barcode[1], " mat col1 ", colnames(mat)[1])
common <- intersect(colnames(mat), id$barcode)
message("intersect ", length(common))
if (length(common) < 1000) {{
  alt <- gsub("_", "-", colnames(mat))
  names(alt) <- colnames(mat)
  hit <- !is.na(match(alt, id$barcode))
  message("dash-intersect ", sum(hit))
  if (sum(hit) < 1000) stop("barcode intersect too small")
  common <- colnames(mat)[hit]
  id <- id[match(alt[common], id$barcode), , drop = FALSE]
  id$barcode <- common
}} else {{
  id <- id[match(common, id$barcode), , drop = FALSE]
}}
mat <- mat[, id$barcode, drop = FALSE]
dir.create("{out_dir}", showWarnings = FALSE, recursive = TRUE)
pats <- unique(id$patient)
for (p in pats) {{
  bc <- id$barcode[id$patient == p]
  sub <- mat[, bc, drop = FALSE]
  od <- file.path("{out_dir}", p)
  dir.create(od, showWarnings = FALSE)
  writeLines(rownames(sub), file.path(od, "genes.txt"))
  writeLines(colnames(sub), file.path(od, "barcodes.txt"))
  Matrix::writeMM(sub, file.path(od, "matrix.mtx"))
  message("wrote ", p, " ", ncol(sub))
}}
message("EXPORT_DONE")
"""
    )
    import subprocess
    say("R export GSE205335 (this loads the full UMI matrix)")
    proc = subprocess.run(["Rscript", str(r)], check=False)
    if proc.returncode != 0:
        raise SystemExit(f"GSE205335 export failed with code {proc.returncode}")


def write_gse205335_keep(geo: Path, sample_map: Path, locked_ids: set[str], dest: Path) -> None:
    import gzip as gz
    soft = parse_sample_map(sample_map)
    rows_out = []
    with gz.open(geo / "GSE205335_Lung_IO_CellIdentity.txt.gz", "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[idx["lineage.sub"]] != "Malignant cells":
                continue
            code = norm_code(p[idx["orig.ident"]])
            info = soft.get(code)
            if info is None:
                continue
            tissue = info["tissue"]
            if tissue.startswith("Normal"):
                continue
            patient = info["patient"]
            if patient not in locked_ids:
                continue
            rows_out.append((p[idx["barcode"]], patient, tissue))
    if len(rows_out) < 1000:
        raise SystemExit(f"GSE205335 malignant keep list only {len(rows_out)}")
    with dest.open("w") as f:
        f.write("barcode\tpatient\ttissue\n")
        for b, p, t in rows_out:
            f.write(f"{b}\t{p}\t{t}\n")
    say(f"GSE205335 keep malignant non-normal {len(rows_out)}")


def run_gse205335(locked, sets, rng, geo: Path, root: Path) -> list[tuple[dict, list[dict]]]:
    out_dir = geo / "gse205335_export"
    keep_tsv = geo / "gse205335_keep.tsv"
    if not (out_dir / locked[0]["unit_id"] / "matrix.mtx").exists():
        write_gse205335_keep(geo, root / "data" / "GSE205335_sample_map.tsv", {r["unit_id"] for r in locked}, keep_tsv)
        export_gse205335_via_r(geo, keep_tsv, out_dir)
    tissue = {}
    with keep_tsv.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            tissue.setdefault(row["patient"], set()).add(row["tissue"])
    results = []
    for unit in locked:
        uid = unit["unit_id"]
        od = out_dir / uid
        say(f"GSE205335 {uid}")
        if not (od / "matrix.mtx").exists():
            rec = {
                "dataset": "GSE205335",
                "unit_id": uid,
                "unit_type": unit["unit_type"],
                "tissue": unit["tissue"],
                "n_malignant_locked": unit["n_malignant"],
                "n_malignant_raw": 0,
                "status": "missing_export",
                "n_qc": 0,
            }
            results.append((rec, []))
            continue
        genes = (od / "genes.txt").read_text().splitlines()
        X = mmread(od / "matrix.mtx").tocsr()
        X, genes = dedup_csr(X, genes)
        # already malignant; still QC
        keep = qc_mask(X, genes)
        n_raw = X.shape[1]
        sub = X[:, keep].tocsr()
        del X
        unit = dict(unit)
        if uid in tissue:
            unit["tissue"] = ",".join(sorted(tissue[uid]))
        say(f"  malignant {n_raw} QC {sub.shape[1]} locked {unit['n_malignant']}")
        results.append(_finish_unit("GSE205335", unit, sub, genes, sets, rng, n_raw))
    return results


def stream_gse131907(geo: Path, locked_ids: set[str]) -> Path:
    """Stream the genes x cells UMI text into one npz per locked sample. Cached."""
    out = geo / "gse131907_mal"
    out.mkdir(parents=True, exist_ok=True)
    flag = out / "DONE"
    if flag.exists():
        say("GSE131907 malignant npz cache present")
        return out
    ann = geo / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    mat = geo / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    if not mat.exists():
        raise SystemExit(f"missing {mat}")
    keep_index = {}
    with gzip.open(ann, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[idx["Cell_subtype"]] != "Malignant cells":
                continue
            sample = p[idx["Sample"]]
            if sample not in locked_ids:
                continue
            keep_index[p[idx["Index"]]] = sample
    say(f"GSE131907 malignant cells to keep {len(keep_index)}")
    say("streaming UMI text (one pass)")
    sample_ids = sorted(locked_ids)
    sid_of = {s: i for i, s in enumerate(sample_ids)}
    genes: list[str] = []
    rows_a = array("I")
    cols_a = array("I")
    vals_a = array("I")
    sample_a = array("B")
    with gzip.open(mat, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        barcodes = header[1:]
        col_of = {b: i for i, b in enumerate(barcodes)}
        missing = 0
        sample_cols: dict[str, list[int]] = {s: [] for s in sample_ids}
        # parallel arrays over kept cells, global column -> (sample code, local column)
        kept_global: list[int] = []
        kept_sid: list[int] = []
        kept_local: list[int] = []
        for index, sample in keep_index.items():
            c = col_of.get(index)
            if c is None:
                missing += 1
                continue
            local = len(sample_cols[sample])
            sample_cols[sample].append(c)
            kept_global.append(c)
            kept_sid.append(sid_of[sample])
            kept_local.append(local)
        say(f"  matched columns {len(kept_global)} missing {missing}")
        for gi, line in enumerate(f):
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[0])
            for j, c in enumerate(kept_global):
                v = parts[c + 1]
                if v in ("0", "0.0", ""):
                    continue
                fv = float(v)
                if fv == 0.0:
                    continue
                rows_a.append(gi)
                cols_a.append(kept_local[j])
                sample_a.append(kept_sid[j])
                vals_a.append(int(fv) if fv.is_integer() and fv < 4294967295 else int(min(fv, 4294967294)))
            if gi and gi % 4000 == 0:
                say(f"  streamed {gi} genes nnz={len(vals_a)}")
    say(f"stream done genes={len(genes)} nnz={len(vals_a)}")
    (out / "genes.txt").write_text("\n".join(genes) + "\n")
    if len(vals_a):
        sample_np = np.frombuffer(sample_a, dtype=np.uint8)
        # frombuffer on array('B') may be a view; copy so we can free sample_a later
        sample_np = np.array(sample_np, copy=True)
    else:
        sample_np = np.zeros(0, dtype=np.uint8)
    for sample in sample_ids:
        ncells = len(sample_cols[sample])
        code = sid_of[sample]
        sel = np.flatnonzero(sample_np == code) if len(sample_np) else np.zeros(0, dtype=int)
        if len(sel) == 0:
            sparse.save_npz(out / f"{sample}.npz", sparse.csr_matrix((len(genes), max(ncells, 0))))
            say(f"  wrote {sample} cells={ncells} nnz=0")
            continue
        rr = np.fromiter((rows_a[i] for i in sel), dtype=np.int32, count=len(sel))
        cc = np.fromiter((cols_a[i] for i in sel), dtype=np.int32, count=len(sel))
        vv = np.fromiter((vals_a[i] for i in sel), dtype=np.float64, count=len(sel))
        M = sparse.csr_matrix((vv, (rr, cc)), shape=(len(genes), ncells))
        sparse.save_npz(out / f"{sample}.npz", M)
        say(f"  wrote {sample} cells={ncells} nnz={M.nnz}")
    flag.write_text("ok\n")
    return out


def run_gse131907(locked, sets, rng, geo: Path) -> list[tuple[dict, list[dict]]]:
    cache = stream_gse131907(geo, {r["unit_id"] for r in locked})
    genes = (cache / "genes.txt").read_text().splitlines()
    results = []
    for unit in locked:
        uid = unit["unit_id"]
        say(f"GSE131907 {uid}")
        X = sparse.load_npz(cache / f"{uid}.npz").tocsr()
        X, genes_u = dedup_csr(X, genes)
        keep = qc_mask(X, genes_u)
        n_raw = X.shape[1]
        sub = X[:, keep].tocsr()
        del X
        say(f"  malignant {n_raw} QC {sub.shape[1]} locked {unit['n_malignant']}")
        results.append(_finish_unit("GSE131907", unit, sub, genes_u, sets, rng, n_raw))
    return results


UNIT_FIELDS = [
    "dataset", "unit_id", "unit_type", "tissue", "status",
    "n_malignant_locked", "n_malignant_raw", "n_qc", "n_non_nt", "n_kd", "kd_quartile",
    "n_Q1", "n_Q2", "n_Q3", "n_Q4", "q25", "q50", "q75",
    "cldn4_pct", "cldn4_mean_log1p", "frac_zero_cldn4",
    "k", "ndims", "n_hvg", "n_universe",
    "n_ifn", "n_ifna", "n_ifng", "n_mhc", "n_sperm", "n_matched_ifn", "n_matched_sperm",
    "spearman_cldn4_ifn", "spearman_p_cells",
    "cosine_ifn", "cosine_ifn_p",
    "unmatched_ifn_q1_minus_q4",
]
for _q in ("Q1", "Q2", "Q3", "Q4"):
    UNIT_FIELDS += [
        f"ifn_{_q}", f"ifna_{_q}", f"ifng_{_q}", f"mhc_{_q}", f"sperm_{_q}",
        f"matched_{_q}", f"delta_{_q}", f"sperm_delta_{_q}", f"unmatched_ifn_{_q}",
    ]

CELL_FIELDS = [
    "dataset", "unit_id", "quartile", "cldn4_count",
    "ifn_prtb", "ifna_prtb", "ifng_prtb", "mhc_prtb", "sperm_prtb",
    "matched_prtb", "delta_prtb", "sperm_delta_prtb",
]


def _write_unit_cell(results: list[tuple[dict, list[dict]]], unit_path: Path, cell_path: Path) -> None:
    unit_path.parent.mkdir(parents=True, exist_ok=True)
    with unit_path.open("w") as f:
        w = csv.DictWriter(f, UNIT_FIELDS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for rec, _ in results:
            w.writerow({k: rec.get(k, "") for k in UNIT_FIELDS})
    with gzip.open(cell_path, "wt") as f:
        w = csv.DictWriter(f, CELL_FIELDS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for _, cells in results:
            for c in cells:
                w.writerow({k: c.get(k, "") for k in CELL_FIELDS})


def merge_parts(out_dir: Path) -> None:
    tab = out_dir / "tables"
    parts = sorted((tab / "parts").glob("*_unit.tsv"))
    with (tab / "unit_summary.tsv").open("w") as out:
        header_written = False
        for p in parts:
            with p.open() as f:
                lines = f.readlines()
            if not lines:
                continue
            if not header_written:
                out.writelines(lines)
                header_written = True
            else:
                out.writelines(lines[1:])
    with gzip.open(tab / "cell_scores.tsv.gz", "wt") as out:
        header_written = False
        for p in sorted((tab / "parts").glob("*_cells.tsv.gz")):
            with gzip.open(p, "rt") as f:
                lines = f.readlines()
            if not lines:
                continue
            if not header_written:
                out.writelines(lines)
                header_written = True
            else:
                out.writelines(lines[1:])
    say(f"merged {len(parts)} cohort parts into {tab / 'unit_summary.tsv'}")


def write_results(results: list[tuple[dict, list[dict]]], out_dir: Path) -> None:
    by: dict[str, list[tuple[dict, list[dict]]]] = {}
    for rec, cells in results:
        by.setdefault(rec["dataset"], []).append((rec, cells))
    part = out_dir / "tables" / "parts"
    for dataset, chunk in by.items():
        _write_unit_cell(chunk, part / f"{dataset}_unit.tsv", part / f"{dataset}_cells.tsv.gz")
    merge_parts(out_dir)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim-only", action="store_true")
    ap.add_argument("--datasets", default="GSE123902,GSE189357,GSE205335,GSE131907")
    ap.add_argument("--geo", default=str(GEO))
    args = ap.parse_args()
    simulate()
    if args.sim_only:
        return
    geo = Path(args.geo)
    sets = load_gene_sets(ROOT / "data" / "gene_sets.json")
    locked = load_locked(ROOT / "data" / "locked_patient_units.tsv")
    rng = np.random.default_rng(SEED)
    want = [d for d in args.datasets.split(",") if d]
    runners = {
        "GSE123902": lambda: run_gse123902(locked["GSE123902"], sets, rng, geo, ROOT),
        "GSE189357": lambda: run_gse189357(locked["GSE189357"], sets, rng, geo),
        "GSE205335": lambda: run_gse205335(locked["GSE205335"], sets, rng, geo, ROOT),
        "GSE131907": lambda: run_gse131907(locked["GSE131907"], sets, rng, geo),
    }
    results = []
    for d in want:
        say(f"==== {d} ====")
        results.extend(runners[d]())
        # checkpoint after each cohort
        write_results(results, ROOT / "results")
    say("ALL_DATASETS_DONE")


if __name__ == "__main__":
    main()
