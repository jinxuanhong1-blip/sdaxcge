#!/usr/bin/env python3
"""Concordant-4 Milo neighbourhood DA vs patient-level malignant CLDN4.

Graph: undirected kNN on the committed Harmony embedding (dataset-corrected).
Index sampling: Milo reduced-dimension median refinement.
DA: edgeR glmQLFit + glmQLFTest (the test inside miloR::testNhoods).

The tested covariate is patient-level malignant CLDN4 (high-rich vs low-rich,
and the continuous %pos). Each locked unit is one sample. A patient intercept
is fit only as an identifiability check: with one sample per patient it is
collinear with the group indicator and is not used as a result.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import coo_matrix
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CACHE = ROOT / "results" / "cache"
TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "figures"
RDS = DATA / "seurat_harmony_embeddings.rds"
PATIENTS = DATA / "patient_units.tsv"

K = 30
D = 30
PROP = 0.1
SEED = 1
MIN_PATIENTS = 5
MIN_MALIGNANT_CELLS = 10
SPATIAL_FDR_CUTS = (0.05, 0.10)

# GEO series-matrix patient ids for the 21 tumour-bearing GSE131907 samples
# kept in the concordant-4 object. Each sample is a different patient.
GSE131907_PATIENT = {
    "EBUS_06": "P1006",
    "EBUS_28": "P1028",
    "EBUS_49": "P1049",
    "BRONCHO_58": "P1058",
    "EBUS_10": "P1010",
    "BRONCHO_11": "P1011",
    "EBUS_12": "P1012",
    "EBUS_13": "P1013",
    "EBUS_15": "P1015",
    "EBUS_19": "P1019",
    "EBUS_51": "P1051",
    "NS_02": "P3002",
    "NS_03": "P3003",
    "NS_04": "P3004",
    "NS_06": "P3006",
    "NS_07": "P3007",
    "NS_12": "P3012",
    "NS_13": "P3013",
    "NS_16": "P3016",
    "NS_17": "P3017",
    "NS_19": "P3019",
}

CLASS_COLORS = {
    "malignant": "#b2182b",
    "T/NK": "#2166ac",
    "T": "#2166ac",
    "NK": "#92c5de",
    "myeloid": "#4d9221",
    "B": "#762a83",
    "other": "#737373",
    "mixed": "#bdbdbd",
}


def _json_default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, dict):
        return {str(k): _json_default(v) for k, v in obj.items()}
    raise TypeError(type(obj))


def spatial_fdr_kdistance(pvalues: np.ndarray, k_distance: np.ndarray) -> np.ndarray:
    """miloR::graphSpatialFDR weighting='k-distance' (Dann / cydar)."""
    p = np.asarray(pvalues, dtype=float)
    d = np.asarray(k_distance, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    ok = np.isfinite(p) & np.isfinite(d) & (d > 0)
    if ok.sum() == 0:
        return out
    pv = p[ok]
    w = 1.0 / d[ok]
    w[~np.isfinite(w)] = 1.0
    o = np.argsort(pv, kind="mergesort")
    pv_o = pv[o]
    w_o = w[o]
    raw = w_o.sum() * pv_o / np.cumsum(w_o)
    adj_o = np.minimum(1.0, np.minimum.accumulate(raw[::-1])[::-1])
    adj = np.empty_like(adj_o)
    adj[o] = adj_o
    out[ok] = adj
    return out


def export_rds() -> tuple[pd.DataFrame, np.ndarray]:
    CACHE.mkdir(parents=True, exist_ok=True)
    meta_path = CACHE / "meta.tsv"
    harm_path = CACHE / "harmony.tsv.gz"
    r_code = f"""
    x <- readRDS({json.dumps(str(RDS))})
    stopifnot(identical(rownames(x$meta), rownames(x$harmony)))
    write.table(x$meta, file={json.dumps(str(meta_path))}, sep="\\t", quote=FALSE, col.names=NA)
    hm <- x$harmony[, seq_len(min({D}, ncol(x$harmony))), drop=FALSE]
    gz <- gzfile({json.dumps(str(harm_path))}, "w")
    write.table(hm, file=gz, sep="\\t", quote=FALSE, col.names=NA)
    close(gz)
    cat("exported", nrow(x$meta), ncol(hm), "\\n")
    """
    subprocess.check_call(["Rscript", "-e", r_code])
    meta = pd.read_csv(meta_path, sep="\t", index_col=0)
    harmony = pd.read_csv(harm_path, sep="\t", index_col=0)
    harmony = harmony.loc[meta.index]
    return meta, harmony.to_numpy(dtype=np.float64)


def refined_indices(X: np.ndarray, knn_idx: np.ndarray, prop: float, seed: int) -> np.ndarray:
    """Milo reduced-dim refinement: median of kNN profiles, nearest cell."""
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    n_seed = int(math.floor(prop * n))
    seeds = rng.choice(n, size=n_seed, replace=False)
    med = np.median(X[knn_idx[seeds]], axis=1)
    nn = NearestNeighbors(n_neighbors=1, algorithm="auto", metric="euclidean")
    nn.fit(X)
    idx = nn.kneighbors(med, return_distance=False)[:, 0]
    _, first = np.unique(idx, return_index=True)
    return idx[np.sort(first)].astype(np.int32)


def undirected_members(knn_idx: np.ndarray, indices: np.ndarray) -> list[np.ndarray]:
    """1-hop neighbourhood on the undirected kNN graph, including the index."""
    n, k = knn_idx.shape
    src = np.repeat(np.arange(n, dtype=np.int32), k)
    dst = knn_idx.ravel()
    mask = src != dst
    src, dst = src[mask], dst[mask]
    rows = np.concatenate([src, dst])
    cols = np.concatenate([dst, src])
    data = np.ones(rows.shape[0], dtype=np.float32)
    A = coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    A.sum_duplicates()
    A.data[:] = 1.0
    members = []
    for i in indices:
        nbrs = A.indices[A.indptr[i] : A.indptr[i + 1]]
        cells = np.unique(np.concatenate(([np.int32(i)], nbrs)))
        members.append(cells)
    return members


def count_by_sample(members: list[np.ndarray], sample_codes: np.ndarray, n_samples: int) -> np.ndarray:
    n_cells = sample_codes.shape[0]
    indptr = [0]
    indices = []
    for cells in members:
        indices.extend(int(c) for c in cells)
        indptr.append(len(indices))
    data = np.ones(len(indices), dtype=np.float64)
    M = coo_matrix(
        (data, (np.repeat(np.arange(len(members)), np.diff(indptr)), indices)),
        shape=(len(members), n_cells),
    ).tocsr()
    S = coo_matrix(
        (np.ones(n_cells), (np.arange(n_cells), sample_codes)),
        shape=(n_cells, n_samples),
    ).tocsr()
    return np.asarray((M @ S).toarray(), dtype=np.float64)


def class_fractions(members: list[np.ndarray], class_codes: np.ndarray, n_classes: int) -> np.ndarray:
    n_cells = class_codes.shape[0]
    indptr = [0]
    indices = []
    for cells in members:
        indices.extend(int(c) for c in cells)
        indptr.append(len(indices))
    data = np.ones(len(indices), dtype=np.float64)
    M = coo_matrix(
        (data, (np.repeat(np.arange(len(members)), np.diff(indptr)), indices)),
        shape=(len(members), n_cells),
    ).tocsr()
    C = coo_matrix(
        (np.ones(n_cells), (np.arange(n_cells), class_codes)),
        shape=(n_cells, n_classes),
    ).tocsr()
    return np.asarray((M @ C).toarray(), dtype=np.float64)


def load_patients(meta: pd.DataFrame) -> pd.DataFrame:
    pat = pd.read_csv(PATIENTS, sep="\t")
    pat["unit_key"] = pat["dataset"] + "|" + pat["unit_id"].astype(str)
    meta = meta.copy()
    meta["unit_key"] = meta["dataset"].astype(str) + "|" + meta["unit_id"].astype(str)
    if not set(meta["unit_key"]).issubset(set(pat["unit_key"])):
        missing = sorted(set(meta["unit_key"]) - set(pat["unit_key"]))
        raise SystemExit(f"embedding units missing from patient table: {missing[:8]}")
    # Locked quartile labels are the high/low definition.
    pat["cldn4_high"] = pat["cldn4_quartile"].isin(["Q3", "Q4"]).astype(int)
    pat["cldn4_per10"] = pat["mal_CLDN4_pct"].astype(float) / 10.0
    pat["q4q1"] = np.where(pat["cldn4_quartile"].isin(["Q1", "Q4"]), 1, 0)
    pat["patient_id"] = [
        f"{r.dataset}|{GSE131907_PATIENT.get(str(r.unit_id), str(r.unit_id))}"
        for r in pat.itertuples(index=False)
    ]
    # Join check against the locked single-cohort Spearmans.
    locked = {
        "GSE123902": -0.659,
        "GSE131907": -0.522,
        "GSE205335": -0.435,
        "GSE189357": -0.600,
    }
    for ds, rho_locked in locked.items():
        d = pat.loc[pat["dataset"] == ds]
        rho = stats.spearmanr(d["mal_CLDN4_pct"], d["frac_tnk"]).statistic
        if abs(rho - rho_locked) > 0.02:
            raise SystemExit(f"patient table drifted for {ds}: rho={rho} locked={rho_locked}")
    return pat


def build_sample_meta(meta: pd.DataFrame, pat: pd.DataFrame) -> pd.DataFrame:
    vc = meta.groupby("unit_key").size().rename("n_cells_embedded")
    n_mal = (
        meta.loc[meta["cell_class"] == "malignant"].groupby("unit_key").size().rename("n_malignant_embedded")
    )
    n_tnk = (
        meta.loc[meta["cell_class"].isin(["T", "NK"])].groupby("unit_key").size().rename("n_tnk_embedded")
    )
    sm = pat.set_index("unit_key").join(vc).join(n_mal).join(n_tnk)
    sm["n_malignant_embedded"] = sm["n_malignant_embedded"].fillna(0).astype(int)
    sm["n_tnk_embedded"] = sm["n_tnk_embedded"].fillna(0).astype(int)
    sm["n_cells_embedded"] = sm["n_cells_embedded"].astype(int)
    sm["frac_tnk_embedded"] = sm["n_tnk_embedded"] / sm["n_cells_embedded"]
    sm["site_mLN"] = ((sm["dataset"] == "GSE131907") & (sm["tissue"] == "mLN")).astype(int)
    sm["site_tLB"] = ((sm["dataset"] == "GSE131907") & (sm["tissue"] == "tL/B")).astype(int)
    sm["keep_all"] = 1
    sm["keep_q4q1"] = sm["q4q1"].astype(int)
    sm["keep_131907"] = (sm["dataset"] == "GSE131907").astype(int)
    for ds in sorted(sm["dataset"].unique()):
        sm[f"keep_{ds}"] = (sm["dataset"] == ds).astype(int)
    return sm.sort_index()


def annotate_nhoods(
    members: list[np.ndarray],
    indices: np.ndarray,
    knn_dist: np.ndarray,
    meta: pd.DataFrame,
    sample_codes: np.ndarray,
    class_names: list[str],
    class_counts: np.ndarray,
) -> pd.DataFrame:
    rows = []
    classes = meta["cell_class"].to_numpy()
    for i, cells in enumerate(members):
        cc = class_counts[i]
        total = cc.sum()
        fr = cc / total
        order = np.argsort(-fr)
        top = class_names[int(order[0])]
        if fr[order[0]] < 0.5 or (len(order) > 1 and fr[order[0]] == fr[order[1]]):
            group = "mixed"
        else:
            group = "T/NK" if top in ("T", "NK") else top
        present = np.unique(sample_codes[cells])
        rows.append(
            {
                "nhood_id": f"nh{i}",
                "index_row": int(indices[i]),
                "index_cell": meta.index[int(indices[i])],
                "k_distance": float(knn_dist[int(indices[i]), -1]),
                "n_cells": int(total),
                "n_patients": int(present.size),
                "majority_class": top,
                "majority_frac": float(fr[order[0]]),
                "class_group": group,
                **{f"n_{c}": int(cc[j]) for j, c in enumerate(class_names)},
                **{f"frac_{c}": float(fr[j]) for j, c in enumerate(class_names)},
            }
        )
    return pd.DataFrame(rows)


def write_counts(path: Path, counts: np.ndarray, nhood_ids: list[str], samples: list[str]) -> None:
    df = pd.DataFrame(counts, index=nhood_ids, columns=samples)
    df.index.name = "nhood_id"
    df.to_csv(path, sep="\t")


def contingency(sm: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for ds, d in sm.groupby("dataset"):
        for tissue, t in d.groupby("tissue"):
            rows.append(
                {
                    "dataset": ds,
                    "tissue": tissue,
                    "n": int(len(t)),
                    "n_high": int(t["cldn4_high"].sum()),
                    "n_low": int((1 - t["cldn4_high"]).sum()),
                    "n_q4": int((t["cldn4_quartile"] == "Q4").sum()),
                    "n_q1": int((t["cldn4_quartile"] == "Q1").sum()),
                }
            )
    return pd.DataFrame(rows)


def identifiability(sm: pd.DataFrame) -> dict:
    n_per = sm.groupby("patient_id").size()
    multi = n_per[n_per > 1]
    # Within-patient CLDN4-group variation is undefined when each patient has one row.
    return {
        "n_units": int(len(sm)),
        "n_patients": int(n_per.size),
        "n_patients_with_two_or_more_samples": int(multi.size),
        "within_patient_design": "unidentified",
        "reason": (
            "Each concordant-4 unit is one dissociated sample from one patient. "
            "model.matrix(~ patient_id + cldn4_high) has more columns than rows "
            "and is rank-deficient. There is no second sampled condition inside a patient, "
            "so a within-patient malignant-neighbourhood contrast is not a Milo design."
        ),
    }


def summarize_da(da: pd.DataFrame, annot: pd.DataFrame, model: str) -> dict:
    m = da.merge(annot, on="nhood_id", how="left")
    out = {
        "model": model,
        "n_nhoods": int(len(m)),
        "min_p": float(m["PValue"].min()) if len(m) else None,
        "min_spatial_fdr": float(m["SpatialFDR"].min()) if len(m) else None,
        "min_bh_fdr": float(m["FDR"].min()) if len(m) else None,
    }
    for cut in SPATIAL_FDR_CUTS:
        hit = m["SpatialFDR"] < cut
        out[f"n_spatial_fdr_{cut}"] = int(hit.sum())
        out[f"n_spatial_fdr_{cut}_patients_ge_{MIN_PATIENTS}"] = int(
            (hit & (m["n_patients"] >= MIN_PATIENTS)).sum()
        )
    out["n_bh_fdr_0.10"] = int((m["FDR"] < 0.10).sum())
    # Descriptive direction among class groups. Not an extra test.
    for g, sub in m.groupby("class_group"):
        out[f"n_{g}"] = int(len(sub))
        out[f"median_logFC_{g}"] = float(sub["logFC"].median())
        out[f"n_spatial_fdr_0.10_{g}"] = int((sub["SpatialFDR"] < 0.10).sum())
    m.to_csv(TABLES / f"da_{model}.tsv", sep="\t", index=False)
    return out


def plot_volcano(da_path: Path, title: str, stem: str) -> None:
    m = pd.read_csv(da_path, sep="\t")
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    y = -np.log10(np.clip(m["SpatialFDR"].to_numpy(), 1e-300, 1))
    for g, sub_idx in m.groupby("class_group").groups.items():
        sub = m.loc[sub_idx]
        ax.scatter(
            sub["logFC"],
            -np.log10(np.clip(sub["SpatialFDR"], 1e-300, 1)),
            s=8,
            alpha=0.75,
            c=CLASS_COLORS.get(g, "#333333"),
            label=f"{g} ({len(sub)} nhoods)",
            linewidths=0,
        )
    ax.axhline(-math.log10(0.10), color="#666666", lw=0.6, ls="--")
    ax.axvline(0, color="#666666", lw=0.4)
    ax.set_xlabel("logFC (CLDN4-high-rich − low-rich)")
    ax.set_ylabel("−log10 SpatialFDR")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(FIGS / f"{stem}.png", dpi=160)
    fig.savefig(FIGS / f"{stem}.pdf")
    plt.close(fig)


def plot_class_logfc(da_path: Path, stem: str) -> None:
    m = pd.read_csv(da_path, sep="\t")
    order = ["T/NK", "myeloid", "B", "malignant", "other", "mixed"]
    present = [g for g in order if (m["class_group"] == g).any()]
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    data = [m.loc[m["class_group"] == g, "logFC"].to_numpy() for g in present]
    bp = ax.boxplot(data, tick_labels=present, showfliers=False, patch_artist=True)
    for patch, g in zip(bp["boxes"], present):
        patch.set_facecolor(CLASS_COLORS.get(g, "#cccccc"))
        patch.set_alpha(0.8)
    ax.axhline(0, color="#666666", lw=0.5)
    ax.set_ylabel("neighbourhood logFC")
    ax.set_title("Primary model logFC by neighbourhood class")
    fig.tight_layout()
    fig.savefig(FIGS / f"{stem}.png", dpi=160)
    fig.savefig(FIGS / f"{stem}.pdf")
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    print("export RDS")
    meta, X = export_rds()
    if X.shape[1] != D:
        raise SystemExit(f"expected {D} Harmony dimensions, got {X.shape[1]}")
    pat = load_patients(meta)
    meta = meta.copy()
    meta["unit_key"] = meta["dataset"].astype(str) + "|" + meta["unit_id"].astype(str)
    sm = build_sample_meta(meta, pat)
    # Restrict the graph to units in the locked table (all 65).
    keep_cells = meta["unit_key"].isin(sm.index).to_numpy()
    meta = meta.loc[keep_cells]
    X = X[keep_cells]
    sm = sm.loc[sm.index.isin(meta["unit_key"].unique())].sort_index()
    samples = list(sm.index)
    sample_index = {s: i for i, s in enumerate(samples)}
    sample_codes = meta["unit_key"].map(sample_index).to_numpy(dtype=np.int32)

    print(f"kNN k={K} d={D} cells={len(meta)} units={len(samples)}")
    nn = NearestNeighbors(n_neighbors=K + 1, algorithm="auto", metric="euclidean")
    nn.fit(X)
    dist, idx = nn.kneighbors(X)
    knn_idx = idx[:, 1:].astype(np.int32)
    knn_dist = dist[:, 1:]
    indices = refined_indices(X, knn_idx, PROP, SEED)
    print(f"refined indices {len(indices)} from floor({PROP}*{len(meta)})")
    members = undirected_members(knn_idx, indices)
    sizes = np.array([len(c) for c in members])
    print(f"nhood size median {np.median(sizes):.0f} min {sizes.min()} max {sizes.max()}")

    class_names = ["B", "T", "NK", "myeloid", "malignant", "other"]
    class_map = {c: i for i, c in enumerate(class_names)}
    unknown = set(meta["cell_class"]) - set(class_names)
    if unknown:
        raise SystemExit(f"unexpected cell classes: {unknown}")
    class_codes = meta["cell_class"].map(class_map).to_numpy(dtype=np.int32)
    counts = count_by_sample(members, sample_codes, len(samples))
    class_counts = class_fractions(members, class_codes, len(class_names))
    # Malignant-only abundance: zero the non-malignant contribution.
    mal_code = class_map["malignant"]
    mal_cell = (class_codes == mal_code).astype(np.int32)
    # Recount using only malignant cells by zeroing sample codes of others via a mask.
    mal_members = [cells[mal_cell[cells] == 1] for cells in members]
    counts_mal = count_by_sample(
        [c if len(c) else np.array([0], dtype=np.int32)[:0] for c in mal_members],
        sample_codes,
        len(samples),
    )

    annot = annotate_nhoods(members, indices, knn_dist, meta, sample_codes, class_names, class_counts)
    annot["frac_tnk"] = annot["frac_T"] + annot["frac_NK"]
    annot.to_csv(TABLES / "nhood_annotation.tsv", sep="\t", index=False)
    n_testable = int((annot["n_patients"] >= MIN_PATIENTS).sum())
    print(f"nhoods {len(annot)} with >={MIN_PATIENTS} patients: {n_testable}")
    print(annot["n_patients"].describe().to_string())
    print(annot["class_group"].value_counts().to_string())

    nhood_ids = annot["nhood_id"].tolist()
    write_counts(CACHE / "counts_all.tsv", counts, nhood_ids, samples)
    write_counts(CACHE / "counts_malignant.tsv", counts_mal, nhood_ids, samples)

    # Class totals in the embedded object (calibration, not neighbourhoods).
    class_tot = (
        pd.crosstab(meta["unit_key"], meta["cell_class"])
        .reindex(index=samples, columns=class_names)
        .fillna(0)
    )
    # Collapse T and NK for a T/NK row plus the other classes.
    class_tot_out = pd.DataFrame(
        {
            "T_NK": class_tot["T"] + class_tot["NK"],
            "myeloid": class_tot["myeloid"],
            "B": class_tot["B"],
            "malignant": class_tot["malignant"],
            "other": class_tot["other"],
        },
        index=samples,
    ).T
    class_tot_out.index.name = "nhood_id"
    class_tot_out.to_csv(CACHE / "counts_class.tsv", sep="\t")

    sm_out = sm.reset_index()
    sm_out.to_csv(TABLES / "sample_meta.tsv", sep="\t", index=False)
    contingency(sm).to_csv(TABLES / "group_by_tissue.tsv", sep="\t", index=False)
    ident = identifiability(sm)
    (TABLES / "within_patient_design.json").write_text(json.dumps(ident, indent=2) + "\n")

    # Embedded composition still tracks the locked patient score? Descriptive only.
    rho_emb, p_emb = stats.spearmanr(sm["mal_CLDN4_pct"], sm["frac_tnk_embedded"])
    rho_full, p_full = stats.spearmanr(sm["mal_CLDN4_pct"], sm["frac_tnk"])

    testable = annot.loc[annot["n_patients"] >= MIN_PATIENTS, "nhood_id"]
    # Subset count matrices to the pre-specified testable set for neighbourhood models.
    def subset_counts(src: Path, dest: Path) -> None:
        df = pd.read_csv(src, sep="\t", index_col=0)
        df.loc[testable].to_csv(dest, sep="\t")

    subset_counts(CACHE / "counts_all.tsv", CACHE / "counts_all_testable.tsv")
    mal_keep = annot.loc[annot["n_malignant"] >= MIN_MALIGNANT_CELLS, "nhood_id"]
    mal_df = pd.read_csv(CACHE / "counts_malignant.tsv", sep="\t", index_col=0)
    mal_df.loc[mal_keep].to_csv(CACHE / "counts_malignant_testable.tsv", sep="\t")

    meta_path = str(TABLES / "sample_meta.tsv")
    jobs = [
        {
            "name": "binary_tmm",
            "counts": str(CACHE / "counts_all_testable.tsv"),
            "meta": meta_path,
            "formula": "~ dataset + cldn4_high",
            "coef": "cldn4_high",
            "norm": "TMM",
            "lib_mode": "colsum",
            "lib_column": "",
            "sample_column": "unit_key",
            "keep_column": "keep_all",
        },
        {
            "name": "binary_logms",
            "counts": str(CACHE / "counts_all_testable.tsv"),
            "meta": meta_path,
            "formula": "~ dataset + cldn4_high",
            "coef": "cldn4_high",
            "norm": "logMS",
            "lib_mode": "column",
            "lib_column": "n_cells_embedded",
            "sample_column": "unit_key",
            "keep_column": "keep_all",
        },
        {
            "name": "continuous_tmm",
            "counts": str(CACHE / "counts_all_testable.tsv"),
            "meta": meta_path,
            "formula": "~ dataset + cldn4_per10",
            "coef": "cldn4_per10",
            "norm": "TMM",
            "lib_mode": "colsum",
            "lib_column": "",
            "sample_column": "unit_key",
            "keep_column": "keep_all",
        },
        {
            "name": "site_binary_tmm",
            "counts": str(CACHE / "counts_all_testable.tsv"),
            "meta": meta_path,
            "formula": "~ dataset + site_mLN + site_tLB + cldn4_high",
            "coef": "cldn4_high",
            "norm": "TMM",
            "lib_mode": "colsum",
            "lib_column": "",
            "sample_column": "unit_key",
            "keep_column": "keep_all",
        },
        {
            "name": "q4q1_tmm",
            "counts": str(CACHE / "counts_all_testable.tsv"),
            "meta": meta_path,
            "formula": "~ dataset + cldn4_high",
            "coef": "cldn4_high",
            "norm": "TMM",
            "lib_mode": "colsum",
            "lib_column": "",
            "sample_column": "unit_key",
            "keep_column": "keep_q4q1",
        },
        {
            "name": "gse131907_origin_tmm",
            "counts": str(CACHE / "counts_all_testable.tsv"),
            "meta": meta_path,
            "formula": "~ tissue + cldn4_high",
            "coef": "cldn4_high",
            "norm": "TMM",
            "lib_mode": "colsum",
            "lib_column": "",
            "sample_column": "unit_key",
            "keep_column": "keep_131907",
        },
        {
            "name": "malignant_state_logms",
            "counts": str(CACHE / "counts_malignant_testable.tsv"),
            "meta": meta_path,
            "formula": "~ dataset + cldn4_high",
            "coef": "cldn4_high",
            "norm": "logMS",
            "lib_mode": "column",
            "lib_column": "n_malignant_embedded",
            "sample_column": "unit_key",
            "keep_column": "keep_all",
        },
        {
            "name": "class_logms",
            "counts": str(CACHE / "counts_class.tsv"),
            "meta": meta_path,
            "formula": "~ dataset + cldn4_high",
            "coef": "cldn4_high",
            "norm": "logMS",
            "lib_mode": "column",
            "lib_column": "n_cells_embedded",
            "sample_column": "unit_key",
            "keep_column": "keep_all",
        },
        {
            "name": "patient_intercept_check",
            "counts": str(CACHE / "counts_all_testable.tsv"),
            "meta": meta_path,
            "formula": "~ patient_id + cldn4_high",
            "coef": "cldn4_high",
            "norm": "TMM",
            "lib_mode": "colsum",
            "lib_column": "",
            "sample_column": "unit_key",
            "keep_column": "keep_all",
        },
    ]
    for ds in ["GSE123902", "GSE131907", "GSE189357", "GSE205335"]:
        jobs.append(
            {
                "name": f"binary_tmm_{ds}",
                "counts": str(CACHE / "counts_all_testable.tsv"),
                "meta": meta_path,
                "formula": "~ cldn4_high",
                "coef": "cldn4_high",
                "norm": "TMM",
                "lib_mode": "colsum",
                "lib_column": "",
                "sample_column": "unit_key",
                "keep_column": f"keep_{ds}",
            }
        )

    spec_path = CACHE / "models.json"
    spec_path.write_text(json.dumps(jobs))
    edger_out = CACHE / "edger"
    edger_out.mkdir(parents=True, exist_ok=True)
    print("edgeR")
    subprocess.check_call(
        ["Rscript", str(ROOT / "scripts" / "edger_qlf.R"), str(spec_path), str(edger_out)]
    )

    summaries = []
    design_rows = []
    for job in jobs:
        design_path = edger_out / f"{job['name']}.design.tsv"
        if design_path.exists():
            design_rows.append(pd.read_csv(design_path, sep="\t"))
        da_path = edger_out / f"{job['name']}.tsv"
        if not da_path.exists():
            summaries.append({"model": job["name"], "status": "not_fit"})
            continue
        da = pd.read_csv(da_path, sep="\t")
        if job["name"] == "class_logms":
            da["SpatialFDR"] = da["FDR"]
            # class rows use class name as nhood_id; attach a dummy annotation
            class_annot = pd.DataFrame(
                {
                    "nhood_id": da["nhood_id"],
                    "class_group": da["nhood_id"].replace({"T_NK": "T/NK"}),
                    "n_patients": len(samples),
                    "k_distance": 1.0,
                }
            )
            # SpatialFDR is not defined for 5 disjoint classes; keep BH in SpatialFDR column
            # only so summarize_da can count them, and label the model in the table.
            merged = da.merge(class_annot, on="nhood_id")
            merged.to_csv(TABLES / "da_class_logms.tsv", sep="\t", index=False)
            summaries.append(
                {
                    "model": job["name"],
                    "n_nhoods": int(len(merged)),
                    "min_p": float(merged["PValue"].min()),
                    "min_bh_fdr": float(merged["FDR"].min()),
                    "note": "class totals, BH-FDR, not neighbourhood SpatialFDR",
                    **{
                        f"logFC_{r.nhood_id}": float(r.logFC)
                        for r in merged.itertuples(index=False)
                    },
                    **{
                        f"fdr_{r.nhood_id}": float(r.FDR)
                        for r in merged.itertuples(index=False)
                    },
                }
            )
            continue
        da["SpatialFDR"] = spatial_fdr_kdistance(
            da["PValue"].to_numpy(),
            annot.set_index("nhood_id").loc[da["nhood_id"], "k_distance"].to_numpy(),
        )
        summaries.append(summarize_da(da, annot, job["name"]))

    if design_rows:
        pd.concat(design_rows, ignore_index=True).to_csv(TABLES / "design_rank.tsv", sep="\t", index=False)

    # Group sizes
    n_high = int(sm["cldn4_high"].sum())
    n_low = int((sm["cldn4_high"] == 0).sum())
    n_q4 = int((sm["cldn4_quartile"] == "Q4").sum())
    n_q1 = int((sm["cldn4_quartile"] == "Q1").sum())
    mixed_interface = int(((annot["n_malignant"] >= 3) & ((annot["n_T"] + annot["n_NK"]) >= 3)).sum())
    summary = {
        "n_units": int(len(sm)),
        "n_high_q3q4": n_high,
        "n_low_q1q2": n_low,
        "n_q4": n_q4,
        "n_q1": n_q1,
        "n_cells": int(len(meta)),
        "k": K,
        "d": D,
        "prop": PROP,
        "seed": SEED,
        "n_nhoods_sampled": int(len(annot)),
        "n_nhoods_testable_min_patients": n_testable,
        "min_patients": MIN_PATIENTS,
        "nhood_size_median": float(np.median(sizes)),
        "n_interface_nhoods_mal_and_tnk_ge3": mixed_interface,
        "spearman_embedded_frac_tnk": {"rho": float(rho_emb), "p": float(p_emb)},
        "spearman_full_frac_tnk_pooled_not_dl": {"rho": float(rho_full), "p": float(p_full)},
        "within_patient": ident,
        "models": summaries,
        "class_group_counts": annot["class_group"].value_counts().to_dict(),
        "n_patients_per_nhood": {
            "min": int(annot["n_patients"].min()),
            "median": float(annot["n_patients"].median()),
            "max": int(annot["n_patients"].max()),
        },
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=_json_default) + "\n")
    if (TABLES / "da_binary_tmm.tsv").exists():
        plot_volcano(
            TABLES / "da_binary_tmm.tsv",
            "Concordant-4 Milo: CLDN4-high-rich vs low-rich",
            "volcano_binary_tmm",
        )
        plot_class_logfc(TABLES / "da_binary_tmm.tsv", "logfc_by_class_binary_tmm")
    print(json.dumps({k: summary[k] for k in summary if k != "models"}, indent=2))
    for s in summaries:
        brief = {
            k: s[k]
            for k in s
            if k
            in (
                "model",
                "status",
                "n_nhoods",
                "min_p",
                "min_spatial_fdr",
                "min_bh_fdr",
                "n_spatial_fdr_0.05",
                "n_spatial_fdr_0.1",
                "n_bh_fdr_0.10",
            )
            or str(k).startswith("n_spatial_fdr_0.1_")
            or str(k).startswith("median_logFC_")
            or str(k).startswith("logFC_")
            or str(k).startswith("fdr_")
        }
        print(json.dumps(brief))


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as exc:
        sys.exit(exc.returncode)
