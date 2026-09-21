#!/usr/bin/env python3
"""CLDN4-high vs CLDN4-low malignant embeddings, GSE131907 only.

Pre-specified before looking at test statistics. See METHODS.md.
"""

from __future__ import annotations

import gzip
import json
import os
import pickle
import time
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.stats import rankdata, wilcoxon
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from transformers import BertForMaskedLM

from models_embed import (
    embed_geneformer,
    embed_scgpt,
    load_scgpt,
    tokenize_geneformer,
    tokenize_scgpt,
)

ROOT = Path(__file__).resolve().parent
CACHE = Path(os.environ.get("GF_CACHE", "/tmp/gf_cache"))
TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "results" / "figures"

SEED = 1
PRE_PER_ARM = 70
CAP_PER_ARM = 50
MIN_ARM = 5
CALIPER = 0.25
MIN_PAIRS = 5
MIN_UMI = 500
MIN_GENES = 200
MAX_MITO = 0.20
ELIG_ORIGIN = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
MIN_MAL = 20
CLDN4_ENS = "ENSG00000189143"

# Seconds/cell above which a heavier model is not run on the full cohort.
HEAVY_SEC_PER_CELL = 1.5

PROGRAM_ORDER = [
    "IFN",
    "MHC-I/APM",
    "TJ",
    "chemokine",
    "EMT",
    "hypoxia",
    "G2M",
    "E2F",
    "TNFA",
    "inflammatory",
    "apical_junction",
    "barrier_keratin",
    "AT2",
    "AT1",
    "club",
    "basal",
    "ciliated",
]
CONTROLS = ["CLDN4", "TACSTD2"]

CHEMOKINE = [
    "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CXCL16", "CXCL8", "CX3CL1",
    "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL19", "CCL21", "CCL22",
    "XCL1", "XCL2", "IL15", "IL2", "IL7", "IL21",
    "CXCR3", "CXCR4", "CXCR5", "CCR5", "CCR7",
]
STATES = {
    "AT2": ["SFTPC", "SFTPB", "SFTPA1", "NAPSA", "LAMP3", "ABCA3"],
    "AT1": ["AGER", "PDPN", "CAV1"],
    "club": ["SCGB1A1", "SCGB3A2", "SCGB3A1"],
    "basal": ["KRT5", "KRT15", "TP63", "NGFR"],
    "ciliated": ["FOXJ1", "TPPP3", "PIFO"],
    "barrier_keratin": ["KRT8", "KRT18", "KRT19", "KRT7", "CDKN1A", "PLAUR"],
}


def fmt_p(p: float) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_n(x: float, digits: int = 3) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.{digits}g}"


def bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out.tolist()
    pv = p[ok]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    placed = np.empty(n)
    placed[order] = q
    out[np.flatnonzero(ok)] = placed
    return out.tolist()


def paired_means(values: np.ndarray, high: np.ndarray, sample: np.ndarray, min_arm: int = MIN_ARM):
    deltas = []
    samples = []
    low_m = []
    high_m = []
    origins_later = []
    for s in pd.unique(sample):
        m = sample == s
        h = m & high
        lo = m & ~high
        if int(h.sum()) < min_arm or int(lo.sum()) < min_arm:
            continue
        if not np.isfinite(values[m]).all():
            # Keep the sample only if both arms have finite values.
            hv = values[h]
            lv = values[lo]
            hv = hv[np.isfinite(hv)]
            lv = lv[np.isfinite(lv)]
            if hv.size < min_arm or lv.size < min_arm:
                continue
            hm = float(hv.mean())
            lm = float(lv.mean())
        else:
            hm = float(values[h].mean())
            lm = float(values[lo].mean())
        deltas.append(hm - lm)
        samples.append(s)
        low_m.append(lm)
        high_m.append(hm)
    deltas_a = np.asarray(deltas, dtype=float)
    p = np.nan
    if deltas_a.size >= 6 and np.any(deltas_a != 0):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p = float(wilcoxon(deltas_a, alternative="two-sided", zero_method="wilcox").pvalue)
    elif deltas_a.size >= 6:
        p = 1.0
    return {
        "n": int(deltas_a.size),
        "median": float(np.median(deltas_a)) if deltas_a.size else np.nan,
        "mean": float(deltas_a.mean()) if deltas_a.size else np.nan,
        "n_pos": int(np.sum(deltas_a > 0)) if deltas_a.size else 0,
        "n_neg": int(np.sum(deltas_a < 0)) if deltas_a.size else 0,
        "p": p,
        "samples": samples,
        "deltas": deltas_a,
        "low_means": np.asarray(low_m, dtype=float),
        "high_means": np.asarray(high_m, dtype=float),
    }


def loso_projection(emb: np.ndarray, high: np.ndarray, sample: np.ndarray, min_arm: int = MIN_ARM):
    """Project each cell onto the high-minus-low direction from the other samples."""
    usable = []
    cent_h = {}
    cent_l = {}
    for s in pd.unique(sample):
        m = sample == s
        h = m & high
        lo = m & ~high
        if int(h.sum()) < min_arm or int(lo.sum()) < min_arm:
            continue
        if not np.isfinite(emb[h]).all() or not np.isfinite(emb[lo]).all():
            continue
        usable.append(s)
        cent_h[s] = emb[h].mean(axis=0)
        cent_l[s] = emb[lo].mean(axis=0)
    proj = np.full(emb.shape[0], np.nan, dtype=np.float64)
    for s in usable:
        others = [o for o in usable if o != s]
        if len(others) < 5:
            continue
        delta = np.mean([cent_h[o] - cent_l[o] for o in others], axis=0)
        norm = np.linalg.norm(delta)
        if norm == 0 or not np.isfinite(norm):
            continue
        delta = delta / norm
        m = sample == s
        proj[m] = emb[m] @ delta
    return proj, usable


def residualize_within_sample(y: np.ndarray, x: np.ndarray, sample: np.ndarray) -> np.ndarray:
    resid = np.full(y.shape, np.nan, dtype=np.float64)
    for s in pd.unique(sample):
        m = (sample == s) & np.isfinite(y) & np.isfinite(x)
        if int(m.sum()) < 4:
            continue
        xx = x[m] - x[m].mean()
        yy = y[m]
        var = float(np.dot(xx, xx))
        slope = float(np.dot(xx, yy - yy.mean()) / var) if var > 0 else 0.0
        resid[m] = yy - (yy.mean() + slope * xx)
    return resid


def load_programs() -> dict[str, set[str]]:
    a8 = json.loads((ROOT / "data" / "a8_sets.json").read_text())
    sets = a8["sets"]
    krt = set(sets["KRT_EPITHELIAL"])
    tj = set(sets["KEGG_TIGHT_JUNCTION"]) | set(sets["GOBP_TIGHT_JUNCTION_ORGANIZATION"])
    tj |= {g for g in a8["focal_genes"] if g not in krt}
    programs = {
        "IFN": set(sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"])
        | set(sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"]),
        "MHC-I/APM": set(sets["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"]),
        "TJ": tj,
        "chemokine": set(CHEMOKINE),
        "EMT": set(sets["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]),
        "hypoxia": set(sets["HALLMARK_HYPOXIA"]),
        "G2M": set(sets["HALLMARK_G2M_CHECKPOINT"]),
        "E2F": set(sets["HALLMARK_E2F_TARGETS"]),
        "TNFA": set(sets["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]),
        "inflammatory": set(sets["HALLMARK_INFLAMMATORY_RESPONSE"]),
        "apical_junction": set(sets["HALLMARK_APICAL_JUNCTION"]),
        "barrier_keratin": set(STATES["barrier_keratin"]),
        "AT2": set(STATES["AT2"]),
        "AT1": set(STATES["AT1"]),
        "club": set(STATES["club"]),
        "basal": set(STATES["basal"]),
        "ciliated": set(STATES["ciliated"]),
        "CLDN4": {"CLDN4"},
        "TACSTD2": {"TACSTD2"},
    }
    for genes in programs.values():
        genes.discard("CLDN4")
    programs["CLDN4"].add("CLDN4")
    return programs


def read_annotation(path: Path) -> list[dict]:
    rows = []
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in handle:
            p = line.rstrip("\n").split("\t")
            if p[idx["Cell_subtype"]] != "Malignant cells":
                continue
            origin = p[idx["Sample_Origin"]]
            if origin not in ELIG_ORIGIN:
                continue
            rows.append(
                {
                    "barcode": p[idx["Index"]],
                    "sample": p[idx["Sample"]],
                    "origin": origin,
                }
            )
    return rows


def matrix_header(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    if header[0] != "Index":
        raise RuntimeError(f"unexpected matrix header start: {header[0]}")
    return header[1:]


def cldn4_for_columns(path: Path, columns: np.ndarray) -> np.ndarray:
    """One pass stopped at CLDN4. columns are 0-based barcode positions."""
    out = np.full(columns.size, np.nan, dtype=np.float32)
    want = {int(c): i for i, c in enumerate(columns)}
    with gzip.open(path, "rt") as handle:
        handle.readline()
        for line in handle:
            tab = line.find("\t")
            gene = line[:tab]
            if gene < "CLDN4":
                continue
            if gene != "CLDN4":
                break
            parts = line.rstrip("\n").split("\t")
            for col, i in want.items():
                out[i] = float(parts[col + 1])
            break
    if not np.isfinite(out).all():
        raise RuntimeError("CLDN4 row was not found for every malignant column")
    return out


def stream_cells(path: Path, columns: np.ndarray) -> tuple[list[str], sparse.csr_matrix, np.ndarray, np.ndarray]:
    """Return gene symbols, genes x cells CSR, nUMI, mito UMI."""
    n_cells = int(columns.size)
    umi = np.zeros(n_cells, dtype=np.float64)
    mito = np.zeros(n_cells, dtype=np.float64)
    n_genes = np.zeros(n_cells, dtype=np.int32)
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    genes: list[str] = []
    gene_row: dict[str, int] = {}
    keep = columns.astype(int)
    t0 = time.time()
    with gzip.open(path, "rt") as handle:
        handle.readline()
        for li, line in enumerate(handle):
            parts = line.rstrip("\n").split("\t")
            gene = parts[0]
            row_i = gene_row.get(gene)
            if row_i is None:
                row_i = len(genes)
                gene_row[gene] = row_i
                genes.append(gene)
            is_mito = gene.startswith("MT-")
            for j, col in enumerate(keep):
                v = parts[col + 1]
                if v == "0" or v == "0.0" or v == "":
                    continue
                fv = float(v)
                if fv == 0.0:
                    continue
                rows.append(row_i)
                cols.append(j)
                data.append(fv)
                umi[j] += fv
                n_genes[j] += 1
                if is_mito:
                    mito[j] += fv
            if li and li % 4000 == 0:
                print(f"  streamed {li} genes, nnz={len(data)}, {time.time() - t0:.0f}s", flush=True)
    mat = sparse.coo_matrix((data, (rows, cols)), shape=(len(genes), n_cells), dtype=np.float32).tocsr()
    mat.sum_duplicates()
    print(f"  stream done genes={len(genes)} nnz={mat.nnz} {time.time() - t0:.0f}s", flush=True)
    return genes, mat, umi, n_genes, mito


def select_cells(ann: list[dict], cldn4: np.ndarray) -> tuple[pd.DataFrame, np.ndarray]:
    df = pd.DataFrame(ann)
    df["cldn4_umi"] = cldn4
    df["mal_i"] = np.arange(len(df))
    rng = np.random.default_rng(SEED)
    chosen = []
    inventory = []
    for sample, sub in df.groupby("sample", sort=True):
        origin = sub["origin"].iloc[0]
        n_mal = int(len(sub))
        n_pos = int((sub["cldn4_umi"] > 0).sum())
        n_neg = n_mal - n_pos
        eligible = n_mal >= MIN_MAL and n_pos >= MIN_ARM and n_neg >= MIN_ARM
        take_idx = []
        if eligible:
            pos = sub.index[sub["cldn4_umi"] > 0].to_numpy()
            neg = sub.index[sub["cldn4_umi"] <= 0].to_numpy()
            n_pos_take = min(PRE_PER_ARM, pos.size)
            n_neg_take = min(PRE_PER_ARM, neg.size)
            take_idx = np.concatenate(
                [
                    rng.choice(pos, size=n_pos_take, replace=False),
                    rng.choice(neg, size=n_neg_take, replace=False),
                ]
            )
            chosen.append(take_idx)
        inventory.append(
            {
                "sample": sample,
                "origin": origin,
                "n_malignant": n_mal,
                "n_cldn4_pos": n_pos,
                "n_cldn4_neg": n_neg,
                "cldn4_pct": 100.0 * n_pos / n_mal if n_mal else np.nan,
                "contrast_eligible": bool(eligible),
                "n_preselected": int(len(take_idx)),
            }
        )
    pre_idx = np.concatenate(chosen) if chosen else np.array([], dtype=int)
    inv = pd.DataFrame(inventory)
    return inv, pre_idx


def qc_and_cap(meta: pd.DataFrame) -> pd.DataFrame:
    meta = meta.copy()
    meta["pct_mito"] = np.where(meta["n_umi"] > 0, meta["mito_umi"] / meta["n_umi"], np.nan)
    meta["qc_pass"] = (
        (meta["n_umi"] >= MIN_UMI)
        & (meta["n_genes"] >= MIN_GENES)
        & (meta["pct_mito"] < MAX_MITO)
    )
    rng = np.random.default_rng(SEED)
    keep_flags = np.zeros(len(meta), dtype=bool)
    for sample, sub in meta.groupby("sample", sort=True):
        ok = sub[sub["qc_pass"]]
        for arm, bit in ((True, ok["cldn4_umi"] > 0), (False, ok["cldn4_umi"] <= 0)):
            ix = ok.index[bit].to_numpy()
            if ix.size > CAP_PER_ARM:
                ix = rng.choice(ix, size=CAP_PER_ARM, replace=False)
            keep_flags[ix] = True
        _ = arm
    meta["embedded"] = keep_flags
    return meta


def collapse_ensembl(
    genes: list[str],
    mat_gc: sparse.csr_matrix,
    name_to_ens: dict,
    prefer: set[str],
) -> tuple[list[str], list[str], np.ndarray]:
    """mat_gc is genes x cells. Return symbols, ensembl ids, cells x genes dense."""
    buckets: dict[str, list[int]] = {}
    symbol_of: dict[str, str] = {}
    for i, g in enumerate(genes):
        ens = name_to_ens.get(g)
        if ens is None:
            continue
        buckets.setdefault(ens, []).append(i)
        current = symbol_of.get(ens)
        if current is None or (g in prefer and current not in prefer):
            symbol_of[ens] = g
    ens_ids = list(buckets.keys())
    symbols = [symbol_of[e] for e in ens_ids]
    mapper = np.full(len(genes), -1, dtype=np.int32)
    ens_index = {e: i for i, e in enumerate(ens_ids)}
    for ens, js in buckets.items():
        for j in js:
            mapper[j] = ens_index[ens]
    coo = mat_gc.tocoo()
    keep = mapper[coo.row] >= 0
    collapsed = sparse.coo_matrix(
        (coo.data[keep], (mapper[coo.row[keep]], coo.col[keep])),
        shape=(len(ens_ids), mat_gc.shape[1]),
        dtype=np.float32,
    ).tocsr()
    collapsed.sum_duplicates()
    dense = collapsed.T.toarray().astype(np.float32)
    return symbols, ens_ids, dense


def program_scores(lognorm: np.ndarray, symbols: list[str], programs: dict[str, set[str]]):
    index = {g: i for i, g in enumerate(symbols)}
    scores = {}
    n_used = {}
    for name, genes in programs.items():
        js = [index[g] for g in genes if g in index]
        n_used[name] = len(js)
        if name in CONTROLS:
            if not js:
                scores[name] = np.full(lognorm.shape[0], np.nan)
            else:
                scores[name] = lognorm[:, js[0]]
        elif len(js) < 3:
            scores[name] = np.full(lognorm.shape[0], np.nan)
        else:
            scores[name] = lognorm[:, js].mean(axis=1)
    return scores, n_used


def diffusion_map(x_pca: np.ndarray, n_comps: int = 11, n_neighbors: int = 15):
    n = x_pca.shape[0]
    k = min(n_neighbors + 1, n - 1)
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean")
    nn.fit(x_pca)
    dist, idx = nn.kneighbors(x_pca)
    sigma = dist[:, -1].astype(np.float64)
    sigma[sigma < 1e-8] = 1e-8
    rows = []
    cols = []
    vals = []
    for i in range(n):
        for t in range(1, k):
            j = int(idx[i, t])
            dij = float(dist[i, t])
            w = np.exp(-(dij * dij) / (sigma[i] * sigma[j]))
            rows.append(i)
            cols.append(j)
            vals.append(w)
    w = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=np.float64)
    w = 0.5 * (w + w.T)
    degree = np.asarray(w.sum(axis=1)).ravel()
    degree[degree <= 0] = 1.0
    d_inv_sqrt = sparse.diags(1.0 / np.sqrt(degree))
    sym = d_inv_sqrt @ w @ d_inv_sqrt
    evals, evecs = eigsh(sym, k=min(n_comps, n - 2), which="LM")
    order = np.argsort(evals)[::-1]
    evals = np.real(evals[order])
    evecs = np.real(evecs[:, order])
    right = np.asarray(d_inv_sqrt @ evecs)
    embedded = right * evals
    # Drop the trivial first component.
    return evals, embedded[:, 1:]


def mean_gene_rho(proj: np.ndarray, lognorm: np.ndarray, sample: np.ndarray, high: np.ndarray) -> np.ndarray:
    rhos = []
    for s in pd.unique(sample):
        m = (sample == s) & np.isfinite(proj)
        if int((m & high).sum()) < MIN_ARM or int((m & ~high).sum()) < MIN_ARM:
            continue
        x = proj[m]
        if np.nanstd(x) < 1e-8:
            continue
        rx = rankdata(x)
        rx = (rx - rx.mean()) / rx.std()
        y = lognorm[m]
        ry = rankdata(y, axis=0).astype(np.float64)
        std = ry.std(axis=0)
        std[std < 1e-8] = np.nan
        ry = (ry - ry.mean(axis=0)) / std
        rho = np.nanmean(rx[:, None] * ry, axis=0)
        rhos.append(rho)
    if not rhos:
        return np.full(lognorm.shape[1], np.nan)
    return np.nanmean(np.vstack(rhos), axis=0)


def token_cosine(gene_sum: np.ndarray, gene_n: np.ndarray, delta: np.ndarray, min_n: int = 30):
    ok = gene_n >= min_n
    vec = np.zeros_like(gene_sum)
    vec[ok] = gene_sum[ok] / gene_n[ok, None]
    norms = np.linalg.norm(vec, axis=1)
    dnorm = np.linalg.norm(delta)
    cos = np.full(gene_sum.shape[0], np.nan)
    good = ok & (norms > 0) & np.isfinite(dnorm) & (dnorm > 0)
    cos[good] = vec[good] @ (delta / dnorm) / norms[good]
    return cos


def full_delta(emb: np.ndarray, high: np.ndarray, sample: np.ndarray) -> np.ndarray:
    pieces = []
    for s in pd.unique(sample):
        m = sample == s
        h = m & high
        lo = m & ~high
        if int(h.sum()) < MIN_ARM or int(lo.sum()) < MIN_ARM:
            continue
        if not np.isfinite(emb[h]).all() or not np.isfinite(emb[lo]).all():
            continue
        pieces.append(emb[h].mean(0) - emb[lo].mean(0))
    if not pieces:
        return np.full(emb.shape[1], np.nan)
    return np.mean(pieces, axis=0)


def plot_paired(test: dict, origin_of: dict, path: Path, title: str, ylabel: str):
    if test["n"] == 0:
        return
    origins = [origin_of[s] for s in test["samples"]]
    order = np.argsort(np.asarray(origins))
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    for i in order:
        ax.plot(
            [0, 1],
            [test["low_means"][i], test["high_means"][i]],
            color="#c8c8c8",
            lw=1.0,
            zorder=1,
        )
    ax.scatter(np.zeros(test["n"]), test["low_means"], c="#4c78a8", s=28, zorder=2, label="CLDN4 low")
    ax.scatter(np.ones(test["n"]), test["high_means"], c="#e45756", s=28, zorder=2, label="CLDN4 high")
    ax.set_xticks([0, 1], ["CLDN4 low", "CLDN4 high"])
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_deltas(names: list[str], medians: list[float], qs: list[float], path: Path, title: str, xlabel: str):
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    y = np.arange(len(names))[::-1]
    colors = ["#e45756" if (np.isfinite(m) and m > 0) else "#4c78a8" for m in medians]
    ax.barh(y, medians, color=colors, height=0.72)
    ax.axvline(0, color="#333333", lw=0.8)
    labels = []
    for name, q in zip(names, qs):
        if np.isfinite(q):
            labels.append(f"{name}  q={fmt_p(q)}")
        else:
            labels.append(name)
    ax.set_yticks(y, labels)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_genes(table: pd.DataFrame, path: Path, title: str, value: str):
    sub = table.head(20).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.0, 5.8))
    ax.barh(sub["symbol"], sub[value], color="#54a24b")
    ax.axvline(0, color="#333333", lw=0.8)
    ax.set_xlabel(value)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_tsv(df: pd.DataFrame, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def self_check_models():
    print("self-check Geneformer V1", flush=True)
    model = BertForMaskedLM.from_pretrained(CACHE / "gf_v1", local_files_only=True)
    model.eval()
    layer = int(model.config.num_hidden_layers) - 1
    ids = torch.tensor([[2, 3, 4, 5], [6, 7, 8, 9]])
    with torch.inference_mode():
        out = model(input_ids=ids, attention_mask=torch.ones_like(ids), output_hidden_states=True)
    hs = out.hidden_states
    if len(hs) != model.config.num_hidden_layers + 1:
        raise RuntimeError("unexpected hidden_states length")
    if not torch.isfinite(hs[layer]).all():
        raise RuntimeError("non-finite Geneformer hidden state")
    print(f"  V1 hidden {model.config.hidden_size} layers {model.config.num_hidden_layers} pool index {layer}", flush=True)
    del model
    print("self-check scGPT eager forward", flush=True)
    vocab = json.loads((CACHE / "scgpt" / "vocab.json").read_text())
    pad_id = int(vocab["<pad>"])
    cls_id = int(vocab["<cls>"])
    sc_model = load_scgpt(str(CACHE / "scgpt" / "best_model.pt"), pad_id, n_token=int(max(vocab.values())) + 1)
    genes = torch.tensor([[cls_id, 10, 11, 12]])
    values = torch.tensor([[-2.0, 10.0, 20.0, 5.0]])
    pad = genes.eq(pad_id)
    with torch.inference_mode():
        h = sc_model.encode(genes, values, pad)
    if not torch.isfinite(h).all():
        raise RuntimeError("non-finite scGPT hidden state")
    print(f"  scGPT cls norm {float(h[0, 0].norm()):.3f}", flush=True)
    del sc_model
    return layer


def build_cell_table() -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    stamp = {
        "seed": SEED,
        "pre_per_arm": PRE_PER_ARM,
        "cap_per_arm": CAP_PER_ARM,
        "min_arm": MIN_ARM,
        "min_umi": MIN_UMI,
        "min_genes": MIN_GENES,
        "max_mito": MAX_MITO,
    }
    cache_npz = CACHE / "gse131907_cldn4_cells.npz"
    cache_stamp = CACHE / "gse131907_cldn4_cells_stamp.json"
    if cache_npz.exists() and cache_stamp.exists() and json.loads(cache_stamp.read_text()) == stamp:
        print("reusing streamed cell cache", flush=True)
        z = np.load(cache_npz, allow_pickle=True)
        meta = pd.read_csv(CACHE / "gse131907_cldn4_meta.tsv", sep="\t")
        counts = z["counts"]
        symbols = z["symbols"].tolist()
        return meta, counts, symbols

    ann_path = CACHE / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    mat_path = CACHE / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    print("reading annotation", flush=True)
    ann = read_annotation(ann_path)
    print(f"  author-malignant in eligible origins: {len(ann)}", flush=True)
    header = matrix_header(mat_path)
    col_of = {b: i for i, b in enumerate(header)}
    missing = [r["barcode"] for r in ann if r["barcode"] not in col_of]
    if missing:
        raise RuntimeError(f"{len(missing)} malignant barcodes missing from the matrix header")
    columns = np.asarray([col_of[r["barcode"]] for r in ann], dtype=np.int32)
    print("CLDN4 census", flush=True)
    cldn4 = cldn4_for_columns(mat_path, columns)
    inv, pre_idx = select_cells(ann, cldn4)
    print(f"  contrast-eligible samples {int(inv['contrast_eligible'].sum())} / {len(inv)}", flush=True)
    print(f"  preselected cells {pre_idx.size}", flush=True)
    pre_cols = columns[pre_idx]
    print("streaming preselected cells", flush=True)
    genes, mat_gc, umi, n_genes_arr, mito = stream_cells(mat_path, pre_cols)
    # mat_gc is genes x preselected cells, column order matches pre_idx.
    meta = pd.DataFrame([ann[i] for i in pre_idx]).reset_index(drop=True)
    meta["cldn4_umi_census"] = cldn4[pre_idx]
    meta["n_umi"] = umi
    meta["n_genes"] = n_genes_arr
    meta["mito_umi"] = mito
    # Prefer the matrix CLDN4 value when the symbol is present.
    if "CLDN4" in genes:
        cldn_row = mat_gc[genes.index("CLDN4")].toarray().ravel()
        meta["cldn4_umi"] = cldn_row
    else:
        meta["cldn4_umi"] = meta["cldn4_umi_census"]
    meta = qc_and_cap(meta)
    inv_path_rows = []
    for rec in inv.to_dict(orient="records"):
        sub = meta[meta["sample"] == rec["sample"]]
        emb = sub[sub["embedded"]]
        rec["n_qc_pass"] = int(sub["qc_pass"].sum())
        rec["n_embedded"] = int(sub["embedded"].sum())
        rec["n_embedded_high"] = int(((emb["cldn4_umi"] > 0)).sum())
        rec["n_embedded_low"] = int(((emb["cldn4_umi"] <= 0)).sum())
        rec["in_paired_test"] = bool(
            rec["n_embedded_high"] >= MIN_ARM and rec["n_embedded_low"] >= MIN_ARM
        )
        inv_path_rows.append(rec)
    inv2 = pd.DataFrame(inv_path_rows)
    write_tsv(inv2, TABLES / "sample_inventory.tsv")

    keep = meta["embedded"].to_numpy()
    meta_k = meta.loc[keep].reset_index(drop=True)
    mat_k = mat_gc[:, np.flatnonzero(keep)]
    print("mapping symbols to Ensembl", flush=True)
    with open(CACHE / "gene_name_id_dict_gc30M.pkl", "rb") as handle:
        name_to_ens = pickle.load(handle)
    prefer = set()
    for geneset in load_programs().values():
        prefer |= set(geneset)
    symbols, ens_ids, counts = collapse_ensembl(genes, mat_k, name_to_ens, prefer)
    # Recompute n_umi from the original matrix (all genes), already on meta.
    mapped_umi = counts.sum(axis=1)
    total_umi = meta_k["n_umi"].to_numpy()
    frac = float(mapped_umi.sum() / total_umi.sum()) if total_umi.sum() else np.nan
    print(f"  vocab-mapped UMI fraction {frac:.3f} genes {len(symbols)} cells {len(meta_k)}", flush=True)
    meta_k["mapped_umi_frac_global"] = frac
    # Attach ensembl list via a sidecar.
    np.savez_compressed(
        cache_npz,
        counts=counts.astype(np.float32),
        symbols=np.asarray(symbols),
        ensembl=np.asarray(ens_ids),
    )
    meta_k.to_csv(CACHE / "gse131907_cldn4_meta.tsv", sep="\t", index=False)
    cache_stamp.write_text(json.dumps(stamp))
    # ensembl is inside the npz; return symbols aligned to counts columns.
    meta_k.attrs = {}  # placate type checkers
    # Stash ensembl on a column-aligned file already in npz. Return via counts.
    meta_k["ensembl_ncols"] = len(ens_ids)
    return meta_k, counts, symbols


def load_ensembl(symbols: list[str]) -> list[str]:
    z = np.load(CACHE / "gse131907_cldn4_cells.npz", allow_pickle=True)
    ens = z["ensembl"].tolist()
    if ens and len(ens) != len(symbols):
        raise RuntimeError("ensembl/symbol length mismatch")
    return ens


def orient_components(comps: np.ndarray, high: np.ndarray, sample: np.ndarray) -> np.ndarray:
    out = comps.copy()
    for j in range(out.shape[1]):
        test = paired_means(out[:, j], high, sample)
        if np.isfinite(test["median"]) and test["median"] < 0:
            out[:, j] *= -1
    return out


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(4)
    t_all = time.time()
    attempts = []
    _ = self_check_models()
    attempts.append(
        {
            "model": "Geneformer-V1-10M",
            "status": "self_check_ok",
            "note": "BertForMaskedLM 6x256, 41MB safetensors, CPU",
        }
    )
    attempts.append(
        {
            "model": "scGPT-human",
            "status": "self_check_ok",
            "note": "eager attention in place of flash-attn; weights 196MB",
        }
    )

    meta, counts, symbols = build_cell_table()
    ensembl = load_ensembl(symbols)
    high = (meta["cldn4_umi"].to_numpy() > 0)
    sample = meta["sample"].to_numpy()
    origin = meta["origin"].to_numpy()
    n_umi = meta["n_umi"].to_numpy().astype(np.float64)
    lognorm = np.log1p(counts / n_umi[:, None] * 10000.0).astype(np.float32)
    programs = load_programs()
    scores, n_used = program_scores(lognorm, symbols, programs)
    origin_of = (
        pd.read_csv(TABLES / "sample_inventory.tsv", sep="\t")
        .set_index("sample")["origin"]
        .to_dict()
    )

    # --- expression programs (not an embedding) ---
    prog_rows = []
    for name in PROGRAM_ORDER + CONTROLS:
        test = paired_means(scores[name], high, sample)
        family = "control" if name in CONTROLS else "program"
        prog_rows.append(
            {
                "family": family,
                "program": name,
                "n_genes_in_matrix": n_used[name],
                "n_samples": test["n"],
                "median_delta_high_minus_low": test["median"],
                "mean_delta": test["mean"],
                "n_samples_delta_pos": test["n_pos"],
                "n_samples_delta_neg": test["n_neg"],
                "p_wilcoxon": test["p"],
            }
        )
    prog_df = pd.DataFrame(prog_rows)
    mask_prog = prog_df["family"].eq("program")
    qvals = bh(prog_df.loc[mask_prog, "p_wilcoxon"].tolist())
    prog_df["q_bh"] = np.nan
    prog_df.loc[mask_prog, "q_bh"] = qvals
    write_tsv(prog_df, TABLES / "program_paired_tests.tsv")
    plot_deltas(
        prog_df.loc[mask_prog, "program"].tolist(),
        prog_df.loc[mask_prog, "median_delta_high_minus_low"].tolist(),
        prog_df.loc[mask_prog, "q_bh"].tolist(),
        FIGS / "fig_program_paired.png",
        "Within-sample program score, CLDN4-high minus low",
        "Median paired difference of mean log1p(CP10k)",
    )

    # --- PCA + diffusion, CLDN4 column removed ---
    print("PCA and diffusion map", flush=True)
    cldn_col = symbols.index("CLDN4") if "CLDN4" in symbols else None
    detected = (counts > 0).sum(axis=0)
    hvg_ok = detected >= 10
    if cldn_col is not None:
        hvg_ok[cldn_col] = False
    var = lognorm[:, hvg_ok].var(axis=0)
    hvg_local = np.flatnonzero(hvg_ok)[np.argsort(var)[-2000:]]
    x = lognorm[:, hvg_local].astype(np.float64)
    x = (x - x.mean(axis=0)) / np.clip(x.std(axis=0), 1e-6, None)
    n_pcs = 30
    pca = PCA(n_components=n_pcs, svd_solver="full", random_state=SEED)
    x_pca = pca.fit_transform(x)
    evals, diff = diffusion_map(x_pca, n_comps=11, n_neighbors=15)
    x_pca_o = orient_components(x_pca[:, :10], high, sample)
    diff_o = orient_components(diff[:, :10], high, sample)
    geom_rows = []
    for kind, arr in (("PC", x_pca_o), ("DC", diff_o)):
        tests = []
        for j in range(arr.shape[1]):
            tests.append(paired_means(arr[:, j], high, sample))
        qs = bh([t["p"] for t in tests])
        for j, (t, q) in enumerate(zip(tests, qs)):
            geom_rows.append(
                {
                    "space": kind,
                    "component": f"{kind}{j + 1}",
                    "n_samples": t["n"],
                    "median_delta_high_minus_low": t["median"],
                    "n_samples_delta_pos": t["n_pos"],
                    "p_wilcoxon": t["p"],
                    "q_bh": q,
                    "variance_explained": float(pca.explained_variance_ratio_[j]) if kind == "PC" else np.nan,
                }
            )
    geom_df = pd.DataFrame(geom_rows)
    geom_df["diffusion_eigenvalue_trivial"] = float(evals[0])
    write_tsv(geom_df, TABLES / "pca_diffusion_paired_tests.tsv")
    plot_deltas(
        geom_df["component"].tolist(),
        geom_df["median_delta_high_minus_low"].tolist(),
        geom_df["q_bh"].tolist(),
        FIGS / "fig_pca_diffusion_paired.png",
        "PCA and diffusion components (CLDN4 gene held out of the features)",
        "Median paired difference after orienting median ≥ 0",
    )

    # --- Geneformer V1 ---
    print("Geneformer V1-10M embeddings", flush=True)
    with open(CACHE / "token_dictionary_gc30M.pkl", "rb") as handle:
        token_of = pickle.load(handle)
    with open(CACHE / "gene_median_dictionary_gc30M.pkl", "rb") as handle:
        median_of = pickle.load(handle)
    gf_emb = {}
    gf_raw = {}
    gf_cos_tables = []
    gf_cache = CACHE / "gf_v1_raw.npz"
    reuse_gf = False
    if gf_cache.exists():
        z = np.load(gf_cache)
        if z["holdout"].shape[0] == len(meta) and z["full"].shape[0] == len(meta):
            gf_raw["holdout"] = z["holdout"]
            gf_raw["full"] = z["full"]
            for label in ("holdout", "full"):
                proj, _ = loso_projection(gf_raw[label], high, sample)
                gf_emb[label] = proj
            reuse_gf = True
            attempts.append(
                {
                    "model": "Geneformer-V1-10M",
                    "status": "reused_cache",
                    "n_cells": int(len(meta)),
                    "note": str(gf_cache),
                }
            )
            print("reused Geneformer V1 embeddings", flush=True)
    if not reuse_gf:
        gf_model = BertForMaskedLM.from_pretrained(CACHE / "gf_v1", local_files_only=True)
        gf_model.eval()
        layer = int(gf_model.config.num_hidden_layers) - 1
        vocab_size = int(gf_model.config.vocab_size)
    id_to_ens = {int(v): k for k, v in token_of.items() if isinstance(k, str) and k.startswith("ENSG")}
    ens_to_symbol = {e: s for s, e in zip(symbols, ensembl)}
    t_gf = time.time()
    if reuse_gf:
        label_iter = []
    else:
        label_iter = (("holdout", CLDN4_ENS), ("full", None))
    for label, hold in label_iter:
        print(f" tokenize {label}", flush=True)
        seqs, info = tokenize_geneformer(
            counts,
            ensembl,
            token_of,
            median_of,
            n_umi,
            holdout_ensembl=hold,
            max_len=2048,
        )
        print(f"  {info}", flush=True)
        t1 = time.time()
        emb, gene_sum, gene_n = embed_geneformer(
            gf_model,
            seqs,
            layer_index=layer,
            special_tokens=False,
            batch_size=8,
            vocab_size=vocab_size,
            collect_tokens=True,
        )
        sec = (time.time() - t1) / max(len(seqs), 1)
        print(f"  {label} {sec:.3f} s/cell", flush=True)
        attempts.append(
            {
                "model": f"Geneformer-V1-10M-{label}",
                "status": "ran",
                "sec_per_cell": sec,
                "n_cells": int(len(seqs)),
                "note": json.dumps(info),
            }
        )
        proj, _ = loso_projection(emb, high, sample)
        gf_emb[label] = proj
        gf_raw[label] = emb
        delta = full_delta(emb, high, sample)
        cos = token_cosine(gene_sum, gene_n, delta, min_n=30)
        order = np.argsort(np.nan_to_num(cos, nan=-np.inf))[::-1]
        for rank, tok_i in enumerate(order[:25], start=1):
            ens = id_to_ens.get(int(tok_i), "")
            gf_cos_tables.append(
                {
                    "model": f"Geneformer-V1-10M-{label}",
                    "rank": rank,
                    "ensembl": ens,
                    "symbol": ens_to_symbol.get(ens, ""),
                    "cosine_to_high_minus_low": float(cos[tok_i]),
                    "n_token_observations": int(gene_n[tok_i]),
                }
            )
    if not reuse_gf:
        np.savez_compressed(gf_cache, holdout=gf_raw["holdout"], full=gf_raw["full"])
        del gf_model
    print(f"Geneformer V1 done in {time.time() - t_gf:.0f}s", flush=True)

    # depth residual of the primary (holdout) projection
    gf_resid = residualize_within_sample(gf_emb["holdout"], np.log1p(n_umi), sample)

    embed_rows = []
    projections = {
        "gf_v1_holdout_loso": gf_emb["holdout"],
        "gf_v1_full_loso": gf_emb["full"],
        "gf_v1_holdout_loso_depth_resid": gf_resid,
    }
    for name, vec in projections.items():
        test = paired_means(vec, high, sample)
        embed_rows.append(
            {
                "embedding": name,
                "role": "primary" if name == "gf_v1_holdout_loso" else "sensitivity",
                "n_samples": test["n"],
                "median_delta_high_minus_low": test["median"],
                "mean_delta": test["mean"],
                "n_samples_delta_pos": test["n_pos"],
                "n_samples_delta_neg": test["n_neg"],
                "p_wilcoxon": test["p"],
            }
        )
        if name == "gf_v1_holdout_loso":
            plot_paired(
                test,
                origin_of,
                FIGS / "fig_geneformer_v1_loso_paired.png",
                f"Geneformer V1-10M holdout LOSO  n={test['n']} samples  p={fmt_p(test['p'])}",
                "Leave-one-sample-out projection",
            )
            # origin split, descriptive
            origin_rows = []
            for og in sorted(set(origin_of.values())):
                samples_og = [s for s in test["samples"] if origin_of.get(s) == og]
                if not samples_og:
                    continue
                mask = np.isin(test["samples"], samples_og)
                d = test["deltas"][mask]
                p_og = np.nan
                if d.size >= 6 and np.any(d != 0):
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        p_og = float(wilcoxon(d, alternative="two-sided", zero_method="wilcox").pvalue)
                origin_rows.append(
                    {
                        "embedding": name,
                        "origin": og,
                        "n_samples": int(d.size),
                        "median_delta": float(np.median(d)),
                        "n_pos": int(np.sum(d > 0)),
                        "p_wilcoxon": p_og,
                    }
                )
            write_tsv(pd.DataFrame(origin_rows), TABLES / "geneformer_loso_by_origin.tsv")

    # program vs holdout axis
    axis_rows = []
    for name in PROGRAM_ORDER + CONTROLS:
        # per-sample Spearman of program score vs LOSO projection
        rhos = []
        for s in pd.unique(sample):
            m = (sample == s) & np.isfinite(gf_emb["holdout"]) & np.isfinite(scores[name])
            if int((m & high).sum()) < MIN_ARM or int((m & ~high).sum()) < MIN_ARM:
                continue
            if np.std(gf_emb["holdout"][m]) < 1e-8 or np.std(scores[name][m]) < 1e-8:
                continue
            rho = np.corrcoef(rankdata(gf_emb["holdout"][m]), rankdata(scores[name][m]))[0, 1]
            rhos.append(float(rho))
        rhos_a = np.asarray(rhos, dtype=float)
        p = np.nan
        if rhos_a.size >= 6 and np.any(rhos_a != 0):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                p = float(wilcoxon(rhos_a, alternative="two-sided", zero_method="wilcox").pvalue)
        elif rhos_a.size >= 6:
            p = 1.0
        axis_rows.append(
            {
                "program": name,
                "family": "control" if name in CONTROLS else "program",
                "n_samples": int(rhos_a.size),
                "median_spearman_vs_gf_holdout_loso": float(np.median(rhos_a)) if rhos_a.size else np.nan,
                "p_wilcoxon": p,
            }
        )
    axis_df = pd.DataFrame(axis_rows)
    q_axis = bh(axis_df.loc[axis_df["family"].eq("program"), "p_wilcoxon"].tolist())
    axis_df["q_bh"] = np.nan
    axis_df.loc[axis_df["family"].eq("program"), "q_bh"] = q_axis
    write_tsv(axis_df, TABLES / "program_vs_geneformer_axis.tsv")

    # nearest genes by within-sample Spearman vs holdout LOSO
    print("nearest genes along Geneformer axis", flush=True)
    detect_frac = (counts > 0).mean(axis=0)
    gene_ok = detect_frac >= 0.05
    rho_all = mean_gene_rho(gf_emb["holdout"], lognorm, sample, high)
    gene_tbl = pd.DataFrame(
        {
            "symbol": symbols,
            "ensembl": ensembl,
            "mean_within_sample_spearman": rho_all,
            "detection_frac": detect_frac,
        }
    )
    gene_tbl = gene_tbl[gene_ok].sort_values("mean_within_sample_spearman", ascending=False)
    write_tsv(gene_tbl.head(40), TABLES / "nearest_genes_expression_axis_top.tsv")
    write_tsv(gene_tbl.tail(40).sort_values("mean_within_sample_spearman"), TABLES / "nearest_genes_expression_axis_bottom.tsv")
    plot_genes(
        gene_tbl.head(20),
        FIGS / "fig_nearest_genes_axis.png",
        "Genes tracking the Geneformer holdout axis",
        "mean_within_sample_spearman",
    )
    cos_df = pd.DataFrame(gf_cos_tables)
    if len(cos_df):
        write_tsv(cos_df, TABLES / "nearest_genes_token_cosine.tsv")
    hold_cos = (
        cos_df[cos_df["model"].eq("Geneformer-V1-10M-holdout")].head(20)
        if len(cos_df)
        else cos_df
    )
    if len(hold_cos) and hold_cos["symbol"].astype(bool).any():
        plot_genes(
            hold_cos[hold_cos["symbol"].astype(bool)],
            FIGS / "fig_token_cosine_holdout.png",
            "Token states nearest the high−low vector (CLDN4 token removed)",
            "cosine_to_high_minus_low",
        )

    # --- scGPT ---
    scgpt_proj = {}
    sc_raw = {}
    try:
        print("scGPT-human embeddings", flush=True)
        vocab = json.loads((CACHE / "scgpt" / "vocab.json").read_text())
        pad_id = int(vocab["<pad>"])
        cls_id = int(vocab["<cls>"])
        n_token = int(max(vocab.values())) + 1
        gene_ids = np.asarray([int(vocab[g]) if g in vocab else -1 for g in symbols], dtype=np.int64)
        cldn_id = int(vocab.get("CLDN4", -1))
        sc_model = load_scgpt(str(CACHE / "scgpt" / "best_model.pt"), pad_id, n_token=n_token)
        # speed gate on 4 cells
        probe_counts = counts[:4]
        t_probe = time.time()
        probe_seq = tokenize_scgpt(
            probe_counts, gene_ids, holdout_gene_id=cldn_id, cls_id=cls_id, pad_id=pad_id, seed=SEED
        )
        _ = embed_scgpt(sc_model, probe_seq, pad_id, batch_size=2)
        sec_probe = (time.time() - t_probe) / 4
        print(f"  scGPT probe {sec_probe:.3f} s/cell", flush=True)
        passes = [("holdout", cldn_id)]
        if sec_probe <= HEAVY_SEC_PER_CELL:
            passes.append(("full", None))
        else:
            attempts.append(
                {
                    "model": "scGPT-human-full",
                    "status": "skipped_time",
                    "sec_per_cell": sec_probe,
                    "note": f"probe exceeded {HEAVY_SEC_PER_CELL} s/cell; holdout pass only",
                }
            )
        for label, hold in passes:
            t1 = time.time()
            seqs = tokenize_scgpt(
                counts,
                gene_ids,
                holdout_gene_id=hold,
                cls_id=cls_id,
                pad_id=pad_id,
                seed=SEED,
            )
            emb = embed_scgpt(sc_model, seqs, pad_id, batch_size=2)
            sec = (time.time() - t1) / max(len(seqs), 1)
            attempts.append(
                {
                    "model": f"scGPT-human-{label}",
                    "status": "ran",
                    "sec_per_cell": sec,
                    "n_cells": int(len(seqs)),
                    "note": "CLS L2-normalized; eager attention; CLDN4 held out"
                    if label == "holdout"
                    else "CLS L2-normalized; eager attention; CLDN4 token kept",
                }
            )
            proj, _ = loso_projection(emb, high, sample)
            scgpt_proj[label] = proj
            sc_raw[label] = emb
            projections[f"scgpt_{label}_loso"] = proj
        del sc_model
    except Exception as exc:
        attempts.append({"model": "scGPT-human", "status": "failed", "note": f"{type(exc).__name__}: {exc}"})
        print("scGPT failed", exc, flush=True)

    for name, vec in scgpt_proj.items():
        test = paired_means(vec, high, sample)
        embed_rows.append(
            {
                "embedding": f"scgpt_{name}_loso",
                "role": "primary" if name == "holdout" else "sensitivity",
                "n_samples": test["n"],
                "median_delta_high_minus_low": test["median"],
                "mean_delta": test["mean"],
                "n_samples_delta_pos": test["n_pos"],
                "n_samples_delta_neg": test["n_neg"],
                "p_wilcoxon": test["p"],
            }
        )
        if name == "holdout":
            plot_paired(
                test,
                origin_of,
                FIGS / "fig_scgpt_loso_paired.png",
                f"scGPT-human holdout LOSO  n={test['n']} samples  p={fmt_p(test['p'])}",
                "Leave-one-sample-out projection",
            )

    # --- Geneformer V2-104M cancer: benchmark, run only if fast ---
    v2_path = CACHE / "gf_v2_cancer" / "model.safetensors"
    v2_raw = None
    try:
        if not v2_path.exists() or v2_path.stat().st_size < 400_000_000:
            raise RuntimeError(f"V2 weights incomplete ({v2_path})")
        print("Geneformer V2-104M_CLcancer benchmark", flush=True)
        with open(CACHE / "token_dictionary_gc104M.pkl", "rb") as handle:
            token_v2 = pickle.load(handle)
        with open(CACHE / "gene_median_dictionary_gc104M.pkl", "rb") as handle:
            median_v2 = pickle.load(handle)
        v2 = BertForMaskedLM.from_pretrained(CACHE / "gf_v2_cancer", local_files_only=True)
        v2.eval()
        layer_v2 = int(v2.config.num_hidden_layers) - 1
        seqs_b, _ = tokenize_geneformer(
            counts[:8],
            ensembl,
            token_v2,
            median_v2,
            n_umi[:8],
            holdout_ensembl=CLDN4_ENS,
            max_len=4096,
            cls_id=int(token_v2["<cls>"]),
            eos_id=int(token_v2["<eos>"]),
        )
        t1 = time.time()
        _emb_b, _, _ = embed_geneformer(
            v2, seqs_b, layer_index=layer_v2, special_tokens=True, batch_size=1, collect_tokens=False
        )
        sec = (time.time() - t1) / max(len(seqs_b), 1)
        print(f"  V2 probe {sec:.3f} s/cell", flush=True)
        if sec > HEAVY_SEC_PER_CELL:
            attempts.append(
                {
                    "model": "Geneformer-V2-104M_CLcancer",
                    "status": "skipped_time",
                    "sec_per_cell": sec,
                    "n_cells": 8,
                    "note": f"benchmark only; full cohort not run because probe > {HEAVY_SEC_PER_CELL} s/cell on CPU",
                }
            )
            del v2
        else:
            t1 = time.time()
            seqs, info = tokenize_geneformer(
                counts,
                ensembl,
                token_v2,
                median_v2,
                n_umi,
                holdout_ensembl=CLDN4_ENS,
                max_len=4096,
                cls_id=int(token_v2["<cls>"]),
                eos_id=int(token_v2["<eos>"]),
            )
            emb, _, _ = embed_geneformer(
                v2, seqs, layer_index=layer_v2, special_tokens=True, batch_size=2, collect_tokens=False
            )
            sec_full = (time.time() - t1) / max(len(seqs), 1)
            attempts.append(
                {
                    "model": "Geneformer-V2-104M_CLcancer-holdout",
                    "status": "ran",
                    "sec_per_cell": sec_full,
                    "n_cells": int(len(seqs)),
                    "note": json.dumps(info),
                }
            )
            proj, _ = loso_projection(emb, high, sample)
            v2_raw = emb
            projections["gf_v2_cancer_holdout_loso"] = proj
            test = paired_means(proj, high, sample)
            embed_rows.append(
                {
                    "embedding": "gf_v2_cancer_holdout_loso",
                    "role": "sensitivity",
                    "n_samples": test["n"],
                    "median_delta_high_minus_low": test["median"],
                    "mean_delta": test["mean"],
                    "n_samples_delta_pos": test["n_pos"],
                    "n_samples_delta_neg": test["n_neg"],
                    "p_wilcoxon": test["p"],
                }
            )
            del v2
    except Exception as exc:
        attempts.append(
            {
                "model": "Geneformer-V2-104M_CLcancer",
                "status": "failed",
                "note": f"{type(exc).__name__}: {exc}",
            }
        )
        print("V2 failed", exc, flush=True)

    attempts.append(
        {
            "model": "Geneformer-V2-316M",
            "status": "skipped_weight",
            "note": "1.26GB root checkpoint not loaded. Heavier than V2-104M on a 16GB CPU host; not run.",
        }
    )

    # BH across foundation-model primary tests that ran (holdout LOSO only is primary;
    # report q across the holdout tests that exist: V1 and scGPT).
    embed_df = pd.DataFrame(embed_rows)
    primary_mask = embed_df["embedding"].isin(["gf_v1_holdout_loso", "scgpt_holdout_loso"])
    embed_df["q_bh_primary_holdouts"] = np.nan
    if primary_mask.any():
        embed_df.loc[primary_mask, "q_bh_primary_holdouts"] = bh(
            embed_df.loc[primary_mask, "p_wilcoxon"].tolist()
        )
    write_tsv(embed_df, TABLES / "embedding_loso_tests.tsv")
    write_tsv(pd.DataFrame(attempts), TABLES / "model_attempts.tsv")

    # cell table
    cell = meta[["barcode", "sample", "origin", "cldn4_umi", "n_umi", "n_genes", "pct_mito"]].copy()
    cell["cldn4_high"] = high
    for j in range(5):
        cell[f"PC{j + 1}"] = x_pca_o[:, j]
        cell[f"DC{j + 1}"] = diff_o[:, j]
    for name, vec in projections.items():
        cell[name] = vec
    write_tsv(cell, TABLES / "cell_projections.tsv")

    from depth_match import run_depth_match

    raw_embeddings = {
        "gf_v1_holdout": gf_raw.get("holdout"),
        "gf_v1_full": gf_raw.get("full"),
    }
    raw_embeddings.update({f"scgpt_{k}": v for k, v in sc_raw.items()})
    if v2_raw is not None:
        raw_embeddings["gf_v2_cancer_holdout"] = v2_raw
    run_depth_match(
        meta,
        counts,
        lognorm,
        symbols,
        ensembl,
        programs,
        raw_embeddings,
        tables=TABLES,
        figs=FIGS,
        caliper=CALIPER,
        min_pairs=MIN_PAIRS,
        helpers={
            "paired_means": paired_means,
            "bh": bh,
            "loso_projection": loso_projection,
            "mean_gene_rho": mean_gene_rho,
            "diffusion_map": diffusion_map,
            "orient_components": orient_components,
            "plot_deltas": plot_deltas,
            "plot_paired": plot_paired,
            "plot_genes": plot_genes,
            "write_tsv": write_tsv,
            "program_scores": program_scores,
            "PROGRAM_ORDER": PROGRAM_ORDER,
            "CONTROLS": CONTROLS,
            "SEED": SEED,
        },
    )

    write_finding(embed_df, prog_df, geom_df, axis_df, gene_tbl, pd.DataFrame(attempts))
    print(f"ALL DONE {time.time() - t_all:.0f}s", flush=True)


def write_finding(embed_df, prog_df, geom_df, axis_df, gene_tbl, attempts):
    inv = pd.read_csv(TABLES / "sample_inventory.tsv", sep="\t")
    n_mal_samples = int(len(inv))
    n_paired = int(inv["in_paired_test"].sum())
    n_cells = int(inv["n_embedded"].sum())
    n_high = int(inv["n_embedded_high"].sum())
    n_low = int(inv["n_embedded_low"].sum())

    def row(df, key, value):
        hit = df[df[key] == value]
        return None if hit.empty else hit.iloc[0]

    def line_embed(name: str) -> str:
        r = row(embed_df, "embedding", name)
        if r is None:
            return f"| {name} | NA | NA | NA | NA | NA |"
        return (
            f"| {name} | {int(r['n_samples'])} | {fmt_n(r['median_delta_high_minus_low'])} | "
            f"{int(r['n_samples_delta_pos'])}/{int(r['n_samples'])} | {fmt_p(r['p_wilcoxon'])} | "
            f"{fmt_p(r['q_bh_primary_holdouts']) if np.isfinite(r['q_bh_primary_holdouts']) else 'NA'} |"
        )

    gf = row(embed_df, "embedding", "gf_v1_holdout_loso")
    sanity = row(prog_df, "program", "CLDN4")
    sanity_note = "CLDN4 control paired test is in `program_paired_tests.tsv`."
    if sanity is not None:
        sanity_note = (
            f"CLDN4 itself (positive control) median paired Δ log1p(CP10k) = "
            f"{fmt_n(sanity['median_delta_high_minus_low'])}, "
            f"p={fmt_p(sanity['p_wilcoxon'])}, n={int(sanity['n_samples'])}."
        )

    prog_show = prog_df[prog_df["family"].eq("program")].copy()
    prog_lines = [
        "| program | n_genes | n_samples | median Δ | pos/n | p | q |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for _, r in prog_show.iterrows():
        prog_lines.append(
            f"| {r['program']} | {int(r['n_genes_in_matrix'])} | {int(r['n_samples'])} | "
            f"{fmt_n(r['median_delta_high_minus_low'])} | {int(r['n_samples_delta_pos'])}/{int(r['n_samples'])} | "
            f"{fmt_p(r['p_wilcoxon'])} | {fmt_p(r['q_bh'])} |"
        )
    sig = prog_show[prog_show["q_bh"].lt(0.10)]
    if sig.empty:
        prog_sentence = "No pre-specified program survives BH q<0.10 on the within-sample paired test."
    else:
        bits = [
            f"{r['program']} Δ={fmt_n(r['median_delta_high_minus_low'])} q={fmt_p(r['q_bh'])}"
            for _, r in sig.iterrows()
        ]
        prog_sentence = "Programs with BH q<0.10: " + "; ".join(bits) + "."

    geom_lines = [
        "| component | n | median Δ | pos/n | p | q |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for _, r in geom_df.iterrows():
        geom_lines.append(
            f"| {r['component']} | {int(r['n_samples'])} | {fmt_n(r['median_delta_high_minus_low'])} | "
            f"{int(r['n_samples_delta_pos'])}/{int(r['n_samples'])} | {fmt_p(r['p_wilcoxon'])} | {fmt_p(r['q_bh'])} |"
        )
    dc_sig = geom_df[geom_df["space"].eq("DC") & geom_df["q_bh"].lt(0.10)]
    pc_sig = geom_df[geom_df["space"].eq("PC") & geom_df["q_bh"].lt(0.10)]

    def sig_sentence(frame, label):
        if frame.empty:
            return f"None of the tested {label} survive BH q<0.10."
        bits = [f"{r['component']} Δ={fmt_n(r['median_delta_high_minus_low'])} q={fmt_p(r['q_bh'])}" for _, r in frame.iterrows()]
        return f"{label} with BH q<0.10: " + "; ".join(bits) + "."

    top_genes = ", ".join(gene_tbl.head(12)["symbol"].astype(str).tolist())
    bot = gene_tbl.tail(8).sort_values("mean_within_sample_spearman")
    bot_genes = ", ".join(bot["symbol"].astype(str).tolist())
    cldn_rho = gene_tbl.loc[gene_tbl["symbol"].eq("CLDN4"), "mean_within_sample_spearman"]
    cldn_rho_s = fmt_n(float(cldn_rho.iloc[0])) if len(cldn_rho) else "NA"

    axis_sig = axis_df[axis_df["family"].eq("program") & axis_df["q_bh"].lt(0.10)]
    if axis_sig.empty:
        axis_sentence = "No pre-specified program tracks the Geneformer holdout axis at BH q<0.10."
    else:
        bits = [
            f"{r['program']} median ρ={fmt_n(r['median_spearman_vs_gf_holdout_loso'])} q={fmt_p(r['q_bh'])}"
            for _, r in axis_sig.iterrows()
        ]
        axis_sentence = "Programs tracking the holdout axis at BH q<0.10: " + "; ".join(bits) + "."

    attempt_lines = [
        "| model | status | sec/cell | note |",
        "|---|---|---:|---|",
    ]
    for _, r in attempts.iterrows():
        sec = r["sec_per_cell"] if "sec_per_cell" in r and pd.notna(r["sec_per_cell"]) else np.nan
        note = str(r["note"]).replace("|", "/")
        attempt_lines.append(
            f"| {r['model']} | {r['status']} | {fmt_n(sec) if np.isfinite(sec) else 'NA'} | {note} |"
        )

    def _opt(name: str):
        path = TABLES / name
        if not path.exists():
            return None
        return pd.read_csv(path, sep="\t")

    def _md_table(df: pd.DataFrame | None, cols: list[str], formats: dict) -> str:
        if df is None or df.empty:
            return "_not run_"
        lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
        for _, r in df.iterrows():
            cells = []
            for c in cols:
                val = r[c]
                fmt = formats.get(c)
                if fmt == "p":
                    cells.append(fmt_p(val))
                elif fmt == "n":
                    cells.append(fmt_n(val) if pd.notna(val) else "NA")
                elif fmt == "i":
                    cells.append(str(int(val)) if pd.notna(val) else "NA")
                else:
                    cells.append(str(val))
            lines.append("| " + " | ".join(cells) + " |")
        return "\n".join(lines)

    depth_tests = _opt("depth_gap_tests.tsv")
    depth_arms = _opt("depth_gap_by_sample.tsv")
    matched_prog = _opt("program_paired_tests_matched.tsv")
    matched_emb = _opt("embedding_loso_tests_matched.tsv")
    matched_geom = _opt("pca_diffusion_paired_tests_matched.tsv")
    matched_genes = _opt("nearest_genes_expression_axis_top_matched.tsv")
    matched_axis = _opt("program_vs_geneformer_axis_matched.tsv")
    pairs = _opt("depth_caliper_pairs.tsv")

    if depth_arms is not None and len(depth_arms):
        depth_sentence = (
            f"Unmatched embedded cells: median nUMI high {int(round(float(depth_arms['median_umi_high'].median())))} "
            f"vs low {int(round(float(depth_arms['median_umi_low'].median())))}; "
            f"median genes detected {int(round(float(depth_arms['median_genes_high'].median())))} vs "
            f"{int(round(float(depth_arms['median_genes_low'].median())))} "
            f"({int(len(depth_arms))} samples)."
        )
    else:
        depth_sentence = "Depth table was not written."
    if pairs is not None and len(pairs):
        match_sentence = (
            f"Caliper match |Δ log1p(nUMI)| ≤ {CALIPER}, greedy within sample, ≥{MIN_PAIRS} pairs: "
            f"**{int(len(pairs))} pairs**, **{int(pairs['sample'].nunique())} samples**, "
            f"median |Δ log1p(nUMI)| = {fmt_n(float(pairs['abs_d_log1p_umi'].median()))}."
        )
    else:
        match_sentence = "Caliper match did not return pairs."

    def _sig_prog(df):
        if df is None or df.empty:
            return "Matched program table was not written."
        sig = df[df["family"].eq("program") & df["q_bh"].lt(0.10)]
        if sig.empty:
            return "After the depth caliper, no pre-specified program survives BH q<0.10."
        bits = [
            f"{r['program']} Δ={fmt_n(r['median_delta_high_minus_low'])} q={fmt_p(r['q_bh'])}"
            for _, r in sig.iterrows()
        ]
        return "Depth-matched programs with BH q<0.10: " + "; ".join(bits) + "."

    cos = _opt("nearest_genes_token_cosine.tsv")
    if cos is not None and len(cos) and "symbol" in cos.columns:
        full_cos = cos[cos["model"].astype(str).str.contains("full")]
        hit = full_cos[full_cos["symbol"].eq("CLDN4")]
        if hit.empty:
            token_note = (
                "CLDN4 is not among the top 25 full-token cosine genes, so that "
                "probe is not used as evidence that the marker token was recovered."
            )
        else:
            token_note = (
                f"CLDN4 token-cosine rank in the full model is {int(hit.iloc[0]['rank'])} "
                f"(cosine {fmt_n(float(hit.iloc[0]['cosine_to_high_minus_low']))})."
            )
    else:
        token_note = "Token-cosine table was not written."

    matched_gene_s = "NA"
    if matched_genes is not None and len(matched_genes):
        matched_gene_s = ", ".join(matched_genes.head(12)["symbol"].astype(str).tolist())
    matched_axis_s = "Matched axis table was not written."
    if matched_axis is not None and len(matched_axis):
        sig = matched_axis[matched_axis["family"].eq("program") & matched_axis["q_bh"].lt(0.10)]
        if sig.empty:
            matched_axis_s = "After the depth caliper, no pre-specified program tracks the Geneformer axis at BH q<0.10."
        else:
            bits = [
                f"{r['program']} median ρ={fmt_n(r['median_spearman_vs_gf_holdout_loso'])} q={fmt_p(r['q_bh'])}"
                for _, r in sig.iterrows()
            ]
            matched_axis_s = "Depth-matched programs tracking the holdout axis at BH q<0.10: " + "; ".join(bits) + "."

    depth_md = _md_table(
        depth_tests,
        ["measure", "n_samples", "median_delta_high_minus_low", "n_samples_delta_pos", "p_wilcoxon"],
        {"n_samples": "i", "median_delta_high_minus_low": "n", "n_samples_delta_pos": "i", "p_wilcoxon": "p"},
    )
    matched_prog_md = _md_table(
        None if matched_prog is None else matched_prog[matched_prog["family"].eq("program")],
        ["program", "n_samples", "median_delta_high_minus_low", "n_samples_delta_pos", "p_wilcoxon", "q_bh"],
        {
            "n_samples": "i",
            "median_delta_high_minus_low": "n",
            "n_samples_delta_pos": "i",
            "p_wilcoxon": "p",
            "q_bh": "p",
        },
    )
    matched_emb_md = _md_table(
        matched_emb,
        ["embedding", "n_samples", "n_pairs", "median_delta_high_minus_low", "n_samples_delta_pos", "p_wilcoxon", "q_bh_holdouts"],
        {
            "n_samples": "i",
            "n_pairs": "i",
            "median_delta_high_minus_low": "n",
            "n_samples_delta_pos": "i",
            "p_wilcoxon": "p",
            "q_bh_holdouts": "p",
        },
    )
    matched_geom_md = _md_table(
        matched_geom,
        ["component", "n_samples", "median_delta_high_minus_low", "n_samples_delta_pos", "p_wilcoxon", "q_bh"],
        {
            "n_samples": "i",
            "median_delta_high_minus_low": "n",
            "n_samples_delta_pos": "i",
            "p_wilcoxon": "p",
            "q_bh": "p",
        },
    )

    primary_p = fmt_p(gf["p_wilcoxon"]) if gf is not None else "NA"
    primary_d = fmt_n(gf["median_delta_high_minus_low"]) if gf is not None else "NA"
    primary_n = int(gf["n_samples"]) if gf is not None else 0
    primary_pos = int(gf["n_samples_delta_pos"]) if gf is not None else 0

    text = f"""# FINDING — Geneformer / scGPT embeddings of CLDN4-high vs low malignant cells

ADDITIVE. **CLDN4-only.** GSE131907 author-malignant cells only.
This does **not** re-estimate the concordant-4 patient Spearman
(n=65, ρ=−0.53). It does not use Visium, and it does not merge any
mouse matrix. Sample is the unit (Kim et al. sample id, not patient).
Cells are a **stratified** draw so both arms exist: up to {CAP_PER_ARM}
CLDN4 UMI>0 and {CAP_PER_ARM} CLDN4 UMI=0 per sample after QC. They are
not a prevalence-weighted sample of the malignant compartment.

Primary embedding test: Geneformer **V1-10M** (the 10M model; V2-104M and
V2-316M are the heavier checkpoints), **CLDN4 token removed** before
ranking, leave-one-sample-out projection onto the other samples'
high-minus-low direction. Positive median Δ means the held-out sample's
CLDN4-high cells sit on the high side of a direction they did not help fit.

## Honest n

- Author-malignant samples in tumor-bearing origins with ≥{MIN_MAL} malignant cells: **{n_mal_samples}**.
- Samples in the paired test (≥{MIN_ARM} embedded cells in each arm): **{n_paired}**.
- Embedded cells: **{n_cells}** ({n_high} high / {n_low} low). This is not the test n.
- tLung epithelial states labeled tS1/tS2 are not `Cell_subtype == Malignant cells` and are not in this object.
- Full-sample CLDN4 %pos (all author-malignant cells, not the cap) is in `results/tables/sample_inventory.tsv`.

## Library size comes first

CLDN4 UMI>0 vs UMI=0 inside a sample is not a depth-balanced contrast. {depth_sentence}

{depth_md}

The unmatched program table below moves together with that gap: large gene sets rise in the deeper arm. That is not evidence for a specific program. The caliper match is a sensitivity added after this gap was measured, and it is the contrast used for program claims.

{match_sentence}

## Foundation-model separation (unmatched draw)

Primary result, Geneformer V1-10M holdout LOSO: median Δ = {primary_d}, {primary_pos}/{primary_n} samples positive, Wilcoxon p = {primary_p}.

| embedding | n_samples | median Δ | pos/n | p | q (holdout family) |
|---|---:|---:|---:|---:|---:|
{line_embed("gf_v1_holdout_loso")}
{line_embed("gf_v1_full_loso")}
{line_embed("gf_v1_holdout_loso_depth_resid")}
{line_embed("scgpt_holdout_loso")}
{line_embed("scgpt_full_loso")}
{line_embed("gf_v2_cancer_holdout_loso")}

`gf_v1_holdout_loso_depth_resid` residualizes the projection on log1p(nUMI) inside each sample. q is BH across the holdout tests that were fit (V1 and scGPT when both ran). Full-token and V2 rows are sensitivities and are not in that q.

{sanity_note}

Figure: `results/figures/fig_geneformer_v1_loso_paired.png`.

## Depth-matched contrast

Same embeddings, restricted to the caliper pairs. PCA and the diffusion map are refit on those cells with CLDN4 still held out of the features. Geneformer / scGPT vectors are not refit; only the leave-one-sample-out contrast is recomputed on the matched cells.

{matched_emb_md}

{_sig_prog(matched_prog)}

This cell-level contrast is not the concordant-4 patient pseudobulk (IFN/MHC lower in CLDN4-high units). A positive cell-level Δ does not replace that result.

{matched_prog_md}

{matched_axis_s}

Depth-matched genes with the highest mean within-sample Spearman vs the Geneformer holdout axis: {matched_gene_s}.

{matched_geom_md}

Figures: `fig_geneformer_v1_loso_paired_matched.png`, `fig_program_paired_matched.png`, `fig_nearest_genes_axis_matched.png`, `fig_pca_diffusion_paired_matched.png`.

## Nearest gene programs (unmatched draw)

Expression programs are mean log1p(CP10k). TJ / apical junction / every multi-gene set has **CLDN4 removed**. IFN is Hallmark IFNα ∪ IFNγ. MHC-I/APM is the repo custom set. Barrier/keratin matches the prior epithelial list and does not contain CLDN4. BH is across the 17 programs, not across the CLDN4 or TACSTD2 controls.

{prog_sentence}

{chr(10).join(prog_lines)}

Along the Geneformer holdout axis (within-sample Spearman of the program score vs the LOSO projection, Wilcoxon across samples): {axis_sentence}

Genes with the highest mean within-sample Spearman vs that axis (detection ≥5%; CLDN4's own ρ = {cldn_rho_s}): {top_genes}.

Lowest: {bot_genes}.

Token-space probe (same layer as the cell vector; cosine of the mean contextual token state to the high−low vector; tokens seen in ≥30 cells): `results/tables/nearest_genes_token_cosine.tsv`. The holdout run cannot list CLDN4 because that token was removed. {token_note}

## PCA / diffusion map (CLDN4 held out of the features)

2000 HVGs by variance of log1p(CP10k), CLDN4 excluded, z-scored, 30 PCs (full SVD). Diffusion map: k=15 adaptive Gaussian on those PCs, symmetric normalization, 10 nontrivial components. Components are sign-oriented so the median paired Δ is ≥0; p-values are two-sided and unchanged by that flip. BH is within the 10 PCs and, separately, within the 10 DCs.

{sig_sentence(pc_sig, "PCs")} {sig_sentence(dc_sig, "diffusion components")}

{chr(10).join(geom_lines)}

Trivial diffusion eigenvalue (component dropped): {fmt_n(float(geom_df['diffusion_eigenvalue_trivial'].iloc[0]))}.

## What was installed and what was skipped

| model | decision |
|---|---|
| Geneformer-V1-10M | Ran. 6 layers, hidden 256, max 2048, no CLS. This is the lite checkpoint relative to V2. |
| Geneformer-V2-104M and V2-104M_CLcancer | 418MB, 12×768, max 4096, CLS+EOS. Cancer weights were benchmarked; the full cohort was run only if the probe was ≤{HEAVY_SEC_PER_CELL} s/cell. |
| Geneformer-V2-316M | Not loaded (1.26GB). Too heavy for this 16GB CPU host once V2-104M is the larger model under test. |
| scGPT whole-human (`wanglab/scGPT-human`) | Weights loaded. Official flash-attn kernel is not installed; the Wqkv projection runs as eager scaled-dot-product attention, post-norm, ReLU FFN, CLS pooling, L2-normalized. Binning is the official 51-bin within-cell quantile. Sequences longer than 1199 genes are a seeded subset, not a new sample of cells. |

{chr(10).join(attempt_lines)}

## Reproduce

```bash
python3 -m venv /tmp/venv
/tmp/venv/bin/pip install -r methods/geneformer_cldn4_gse131907/requirements.txt
# torch CPU wheel if the default torch build is unwanted
GF_CACHE=/tmp/gf_cache /tmp/venv/bin/python methods/geneformer_cldn4_gse131907/analyze.py
```

GEO inputs (not committed): `GSE131907_Lung_Cancer_cell_annotation.txt.gz`, `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`.
Model weights (not committed): `ctheodoris/Geneformer` `Geneformer-V1-10M` plus the gc30M dictionaries; optional V2-104M_CLcancer and gc104M dictionaries; `wanglab/scGPT-human`.
"""
    (ROOT / "FINDING.md").write_text(text)
    print("wrote FINDING.md", flush=True)


if __name__ == "__main__":
    main()
