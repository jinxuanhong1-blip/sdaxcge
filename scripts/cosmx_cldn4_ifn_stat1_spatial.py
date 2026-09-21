#!/usr/bin/env python3
"""CosMx He 2022: CLDN4 vs IFN/STAT1 among tumor cells (same cell vs neighbors).

Official public object: figshare 25976224, cosmx_human_nsclc_clustered.h5ad
(He et al. Nat Biotechnol 2022; CellCharter clustering). 8 sections, 5 patients.

Question, pre-specified before looking at correlations:
  Among cells labeled tumor, is CLDN4 anti-correlated with an IFN/STAT1
  module in the same cell, in neighboring tumor cells, or only because
  neighbor IFN tracks the same cell?

NHEJ ligation genes are not on the 960 panel. This script inventories them
and does not score an NHEJ module. cGAS/STING1/TBK1 are inventoried the
same way and are not scored.

Inference unit is the patient (n=5), with section (n=8) reported alongside.
Cell-level Spearman p-values are not used.
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
H5AD = ROOT / "data" / "cosmx" / "cosmx_human_nsclc_clustered.h5ad"
OUT = ROOT / "results" / "cosmx_cldn4_ifn_stat1"
GENE_DIR = ROOT / "resources"

# He et al. 2022 CosMx NSCLC prototype pixel size.
PX_TO_UM = 0.18
RADII_UM = (20, 50, 100)
PRIMARY_RADIUS_UM = 50
MIN_COUNTS = 50
MIN_FOV_TUMOR = 40
MIN_NEIGHBORS = 5
DETECT_MIN = 0.02
MIN_MODULE_GENES = 5
N_PERM = 299
RNG = np.random.default_rng(20260921)

# Direct IFN/STAT1 transcriptional readout. No HLA, no JAK/IFNGR, no CCL5.
# Fixed before correlation. Genes absent from the 960 panel are dropped.
IFN_STAT1_CORE = [
    "STAT1",
    "IFIT1",
    "IFITM1",
    "IFITM3",
    "MX1",
    "OAS1",
    "OAS2",
    "OAS3",
    "OASL",
    "IFI27",
    "CXCL9",
    "CXCL10",
    "IDO1",
    "CD274",
    "DDX58",
    "IFIH1",
    "NLRC5",
]
IFN_APM = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "CIITA"]
IFN_SIGNAL = ["JAK1", "JAK2", "IFNGR1", "IFNGR2", "CCL5"]
KRT_GENES = ["KRT8", "KRT18", "KRT19"]

# Ligation/core NHEJ. Not scored. Adjacent DDR genes are inventory-only.
NHEJ_CORE = [
    "PRKDC",
    "XRCC4",
    "XRCC5",
    "XRCC6",
    "LIG4",
    "NHEJ1",
    "DCLRE1C",
    "POLL",
    "POLM",
    "PAXX",
    "TP53BP1",
    "RIF1",
    "MAD2L2",
]
NHEJ_ADJACENT = ["DNTT", "ATM", "ATR", "CHEK1", "CHEK2", "XRCC1", "LIG3"]
STING_MACHINERY = ["CGAS", "MB21D1", "STING1", "TMEM173", "TBK1", "IRF3", "MAVS", "IFI16"]

IMMUNE_TYPES = {
    "neutrophil",
    "macrophage",
    "plasmablast",
    "T CD4 naive",
    "T CD4 memory",
    "T CD8 naive",
    "T CD8 memory",
    "B-cell",
    "mDC",
    "pDC",
    "mast",
    "monocyte",
    "NK",
    "Treg",
}


def _self_check() -> None:
    scores = np.array([1.0, 2.0, 3.0, 4.0])
    indices = np.array([1, 2, 0, 2, 0, 1], dtype=np.int32)
    indptr = np.array([0, 2, 4, 6], dtype=np.int64)
    means = neighbor_means(scores, indices, indptr)
    expect = np.array([2.5, 2.0, 1.5])
    if not np.allclose(means, expect):
        raise AssertionError((means, expect))


def decode(arr) -> list[str]:
    out = []
    for x in arr:
        out.append(x.decode() if isinstance(x, bytes) else str(x))
    return out


def read_cat(f: h5py.File, name: str) -> np.ndarray:
    cats = decode(f["obs"][name]["categories"][:])
    codes = f["obs"][name]["codes"][:]
    return np.asarray(cats, dtype=object)[codes]


def load_gene_set(path: Path) -> list[str]:
    genes = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(">") or line.startswith("HALLMARK"):
            continue
        genes.append(line.split()[0])
    return genes


def neighbor_means(scores: np.ndarray, indices: np.ndarray, indptr: np.ndarray) -> np.ndarray:
    sums = np.add.reduceat(scores[indices], indptr[:-1])
    return sums / np.diff(indptr)


def residualize(y: np.ndarray, cov: np.ndarray) -> np.ndarray:
    if cov.ndim == 1:
        cov = cov.reshape(-1, 1)
    design = np.column_stack([np.ones(len(y)), cov])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    return y - design @ beta


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 30:
        return np.nan
    xx, yy = x[m], y[m]
    if np.unique(xx).size < 2 or np.unique(yy).size < 2:
        return np.nan
    rho, _ = stats.spearmanr(xx, yy)
    return float(rho)


def residual_spearman(x: np.ndarray, y: np.ndarray, cov: np.ndarray | None) -> float:
    m = np.isfinite(x) & np.isfinite(y)
    if cov is not None:
        m &= np.isfinite(cov).all(axis=1)
    if m.sum() < 30:
        return np.nan
    xx, yy = x[m], y[m]
    if cov is not None:
        cc = cov[m]
        xx = residualize(xx, cc)
        yy = residualize(yy, cc)
    return spearman(xx, yy)


def partial_neighbor_rho(
    cldn4: np.ndarray,
    own: np.ndarray,
    neigh: np.ndarray,
    cov: np.ndarray | None,
) -> float:
    """Spearman of CLDN4 vs neighbor IFN after removing own IFN (and covariates)."""
    m = np.isfinite(cldn4) & np.isfinite(own) & np.isfinite(neigh)
    if cov is not None:
        m &= np.isfinite(cov).all(axis=1)
    if m.sum() < 30:
        return np.nan
    xx, own_m, nn = cldn4[m], own[m], neigh[m]
    if cov is not None:
        cc = cov[m]
        xx = residualize(xx, cc)
        own_m = residualize(own_m, cc)
        nn = residualize(nn, cc)
    nn = residualize(nn, own_m)
    return spearman(xx, nn)


def winsorize(col: np.ndarray, q: float = 0.995) -> np.ndarray:
    hi = np.quantile(col, q)
    return np.minimum(col, hi)


def zscore(col: np.ndarray) -> np.ndarray | None:
    sd = float(col.std())
    if sd < 1e-8:
        return None
    return (col - col.mean()) / sd


def module_from_lognorm(lognorm: np.ndarray, detect: np.ndarray) -> tuple[np.ndarray | None, np.ndarray]:
    """Mean of within-section z-scores. Drops genes below detection or with no variance."""
    use = []
    cols = []
    for j in range(lognorm.shape[1]):
        if detect[j] < DETECT_MIN:
            continue
        z = zscore(winsorize(lognorm[:, j]))
        if z is None:
            continue
        use.append(j)
        cols.append(z)
    if len(cols) < MIN_MODULE_GENES:
        return None, np.zeros(lognorm.shape[1], dtype=bool)
    mask = np.zeros(lognorm.shape[1], dtype=bool)
    mask[np.array(use)] = True
    return np.vstack(cols).mean(axis=0), mask


def median_split_delta(cldn4: np.ndarray, score: np.ndarray) -> tuple[float, str]:
    m = np.isfinite(cldn4) & np.isfinite(score)
    x, y = cldn4[m], score[m]
    if len(x) < 30:
        return np.nan, "too_few"
    med = float(np.median(x))
    if med <= 0:
        hi, lo = x > 0, x <= 0
        rule = "detected_vs_zero"
    else:
        hi, lo = x >= med, x < med
        rule = "median"
    if hi.sum() < 30 or lo.sum() < 30:
        return np.nan, rule
    return float(y[hi].mean() - y[lo].mean()), rule


def sign_test(rhos: list[float]) -> dict:
    vals = [r for r in rhos if np.isfinite(r)]
    n = len(vals)
    n_neg = int(sum(r < 0 for r in vals))
    if n == 0:
        return {"n": 0, "n_neg": 0, "mean": np.nan, "sign_p_onesided": np.nan, "sign_p_twosided": np.nan}
    greater = stats.binomtest(n_neg, n, 0.5, alternative="greater")
    two = stats.binomtest(n_neg, n, 0.5, alternative="two-sided")
    return {
        "n": n,
        "n_neg": n_neg,
        "mean": float(np.mean(vals)),
        "sign_p_onesided": float(greater.pvalue),
        "sign_p_twosided": float(two.pvalue),
    }


def wilcoxon_less(rhos: list[float]) -> float:
    vals = np.array([r for r in rhos if np.isfinite(r)], dtype=float)
    if len(vals) < 5 or np.allclose(vals, 0):
        return np.nan
    try:
        res = stats.wilcoxon(vals, alternative="less", zero_method="wilcox")
    except ValueError:
        return np.nan
    return float(res.pvalue)


def patient_mean_from_sections(patient: list[str], rho: list[float]) -> float:
    df = pd.DataFrame({"patient": patient, "rho": rho})
    df = df[np.isfinite(df["rho"])]
    if df.empty:
        return np.nan
    return float(df.groupby("patient")["rho"].mean().mean())


def ball_query(tree: cKDTree, xy: np.ndarray, radius_px: float):
    try:
        return tree.query_ball_point(xy, r=radius_px, workers=-1)
    except TypeError:
        return tree.query_ball_point(xy, r=radius_px)


def shuffle_within_fov(values: np.ndarray, fov: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    out = values.copy()
    for f in np.unique(fov):
        ix = np.flatnonzero(fov == f)
        if len(ix) < 2:
            continue
        out[ix] = values[ix][rng.permutation(len(ix))]
    return out


def load_matrix():
    if not H5AD.exists():
        raise SystemExit(
            f"Missing {H5AD}. Download figshare file 46841842 "
            "(doi:10.6084/m9.figshare.25976224) to that path."
        )
    f = h5py.File(H5AD, "r")
    genes = decode(f["var"]["_index"][:])
    gmap = {g: i for i, g in enumerate(genes)}
    cell_type = read_cat(f, "cell_type")
    patient = read_cat(f, "patient")
    sample = read_cat(f, "sample")
    fov = f["obs"]["fov"][:].astype(np.int32)
    n_counts_obs = f["obs"]["n_counts"][:].astype(np.float64)
    spatial = f["obsm"]["spatial"][:].astype(np.float64)

    hallmark_g = load_gene_set(GENE_DIR / "hallmark_interferon_gamma_response.txt")
    core = [g for g in IFN_STAT1_CORE if g in gmap]
    apm = [g for g in IFN_APM if g in gmap]
    signal = [g for g in IFN_SIGNAL if g in gmap]
    broad = list(dict.fromkeys(core + apm + signal))
    hallmark = [g for g in hallmark_g if g in gmap and g != "CLDN4"]
    krt = [g for g in KRT_GENES if g in gmap]
    needed = list(dict.fromkeys(["CLDN4", "STAT1"] + broad + hallmark + krt))
    missing_needed = [g for g in needed if g not in gmap]
    if missing_needed:
        raise SystemExit(f"panel missing required genes: {missing_needed}")
    if "CLDN4" not in gmap or "STAT1" not in core:
        raise SystemExit("CLDN4 or STAT1 absent; cannot run the IFN/STAT1 test")
    if len(krt) < 2:
        raise SystemExit(f"keratin controls missing: have {krt}")

    cols = [gmap[g] for g in needed]
    indptr = f["layers"]["counts"]["indptr"][:]
    indices = f["layers"]["counts"]["indices"][:]
    data = f["layers"]["counts"]["data"][:]
    mat = sparse.csr_matrix((data, indices, indptr), shape=(len(cell_type), len(genes)))
    # Row sums of the stored counts, to confirm obs n_counts is the library size.
    row_sum = np.asarray(mat.sum(axis=1)).ravel()
    counts = mat[:, cols].toarray().astype(np.float32)
    del mat, indptr, indices, data
    f.close()

    agree = np.corrcoef(row_sum[:5000], n_counts_obs[:5000])[0, 1]
    if not np.isfinite(agree) or agree < 0.99:
        raise SystemExit(f"obs n_counts does not match count-matrix row sums (r={agree})")

    panel = set(genes)
    inventory = {
        "n_genes_panel": len(genes),
        "nhej_core_present": [g for g in NHEJ_CORE if g in panel],
        "nhej_core_absent": [g for g in NHEJ_CORE if g not in panel],
        "nhej_adjacent_present": [g for g in NHEJ_ADJACENT if g in panel],
        "nhej_adjacent_absent": [g for g in NHEJ_ADJACENT if g not in panel],
        "sting_machinery_present": [g for g in STING_MACHINERY if g in panel],
        "sting_machinery_absent": [g for g in STING_MACHINERY if g not in panel],
        "ifn_stat1_core_present": core,
        "ifn_stat1_core_absent": [g for g in IFN_STAT1_CORE if g not in panel],
        "ifn_broad_present": broad,
        "hallmark_ifng_present": hallmark,
        "hallmark_ifng_n_gmt": len(hallmark_g),
        "keratin_present": krt,
        "score_nhej": False,
        "reason_no_nhej": (
            "Core NHEJ ligation genes (PRKDC, XRCC4, XRCC5, XRCC6, LIG4, NHEJ1, "
            "DCLRE1C) are absent from the 960 panel. ATM/ATR/CHEK are inventory-only "
            "and are not an NHEJ module."
        ),
    }
    meta = {
        "cell_type": cell_type,
        "patient": patient,
        "sample": sample,
        "fov": fov,
        "n_counts": n_counts_obs,
        "spatial": spatial,
        "genes": needed,
        "counts": counts,
        "inventory": inventory,
    }
    return meta


def prepare_sections(meta: dict) -> tuple[list[dict], dict]:
    cell_type = meta["cell_type"]
    tumor = np.array([str(x).startswith("tumor") for x in cell_type])
    immune = np.array([str(x) in IMMUNE_TYPES for x in cell_type])
    qc = meta["n_counts"] >= MIN_COUNTS
    genes = meta["genes"]
    gix = {g: i for i, g in enumerate(genes)}
    core_idx = np.array([gix[g] for g in meta["inventory"]["ifn_stat1_core_present"]])
    broad_idx = np.array([gix[g] for g in meta["inventory"]["ifn_broad_present"]])
    hall_idx = np.array([gix[g] for g in meta["inventory"]["hallmark_ifng_present"]])
    krt_idx = np.array([gix[g] for g in meta["inventory"]["keratin_present"]])
    cldn_i = gix["CLDN4"]
    stat1_i = gix["STAT1"]

    # Pixel scale check on one FOV of the first sample.
    samples = list(dict.fromkeys(meta["sample"].tolist()))
    s0 = samples[0]
    m0 = (meta["sample"] == s0) & (meta["fov"] == meta["fov"][meta["sample"] == s0][0])
    xy0 = meta["spatial"][m0]
    d0, _ = cKDTree(xy0).query(xy0, k=2)
    med_nn_um = float(np.median(d0[:, 1]) * PX_TO_UM)

    sections = []
    for sample in samples:
        in_sample = meta["sample"] == sample
        base = in_sample & qc
        fov_s = meta["fov"][base & tumor]
        # FOVs large enough, counted on QC tumor cells.
        fov_counts = pd.Series(fov_s).value_counts()
        big_fovs = set(fov_counts[fov_counts >= MIN_FOV_TUMOR].index.tolist())
        tumor_m = base & tumor & np.isin(meta["fov"], list(big_fovs))
        univ_m = base & np.isin(meta["fov"], list(big_fovs))
        if tumor_m.sum() < 100:
            print(f"{sample}: only {tumor_m.sum()} QC tumor cells, skipped", flush=True)
            continue
        xy_t = meta["spatial"][tumor_m]
        fov_t = meta["fov"][tumor_m]
        ncount_t = meta["n_counts"][tumor_m]
        raw_t = meta["counts"][tumor_m]
        log_t = np.log1p(raw_t * (1e4 / np.maximum(ncount_t, 1.0))[:, None]).astype(np.float32)

        detect = (raw_t > 0).mean(axis=0)
        core_score, core_use = module_from_lognorm(log_t[:, core_idx], detect[core_idx])
        broad_score, broad_use = module_from_lognorm(log_t[:, broad_idx], detect[broad_idx])
        hall_score, hall_use = module_from_lognorm(log_t[:, hall_idx], detect[hall_idx])
        if core_score is None:
            raise SystemExit(f"{sample}: IFN/STAT1 core module has <{MIN_MODULE_GENES} genes")
        krt_cols = []
        for j in krt_idx:
            z = zscore(winsorize(log_t[:, j]))
            if z is not None:
                krt_cols.append(z)
        if len(krt_cols) < 2:
            raise SystemExit(f"{sample}: keratin score failed")
        krt = np.vstack(krt_cols).mean(axis=0).astype(np.float64)
        cldn4 = log_t[:, cldn_i].astype(np.float64)
        stat1 = log_t[:, stat1_i].astype(np.float64)
        logn = np.log1p(ncount_t).astype(np.float64)

        # Immune fraction among QC cells in the same kept FOVs.
        xy_u = meta["spatial"][univ_m]
        immune_u = immune[univ_m]
        tumor_pos = np.flatnonzero(tumor[univ_m])
        # tumor_m is univ tumor, so tumor_pos order matches xy_t if both are stable.
        # Recompute tumor coordinates from the universe to keep index alignment explicit.
        if len(tumor_pos) != len(xy_t):
            raise SystemExit(f"{sample}: tumor index mismatch {len(tumor_pos)} vs {len(xy_t)}")
        tree_u = cKDTree(xy_u)
        balls = ball_query(tree_u, xy_u[tumor_pos], PRIMARY_RADIUS_UM / PX_TO_UM)
        imm_frac = np.empty(len(tumor_pos), dtype=np.float64)
        n_all_nb = np.empty(len(tumor_pos), dtype=np.int32)
        for i, ball in enumerate(balls):
            # ball indices are into univ; drop self (tumor_pos[i]).
            nb = [j for j in ball if j != tumor_pos[i]]
            n_all_nb[i] = len(nb)
            if len(nb) < MIN_NEIGHBORS:
                imm_frac[i] = np.nan
            else:
                nb_a = np.asarray(nb, dtype=np.int32)
                imm_frac[i] = float(immune_u[nb_a].mean())

        # Tumor-tumor neighbor graphs at each radius. Indices into this section's tumor array.
        tree_t = cKDTree(xy_t)
        graphs = {}
        for radius in RADII_UM:
            balls_t = ball_query(tree_t, xy_t, radius / PX_TO_UM)
            elig = []
            lists = []
            for i, ball in enumerate(balls_t):
                nb = [j for j in ball if j != i]
                if len(nb) >= MIN_NEIGHBORS:
                    elig.append(i)
                    lists.append(nb)
            if len(elig) < 100:
                graphs[radius] = None
                continue
            lengths = np.fromiter((len(x) for x in lists), dtype=np.int32, count=len(lists))
            indptr = np.zeros(len(lists) + 1, dtype=np.int64)
            np.cumsum(lengths, out=indptr[1:])
            indices = np.empty(int(indptr[-1]), dtype=np.int32)
            pos = 0
            for nb in lists:
                n = len(nb)
                indices[pos : pos + n] = nb
                pos += n
            graphs[radius] = {
                "elig": np.asarray(elig, dtype=np.int32),
                "indices": indices,
                "indptr": indptr,
            }
        patient = str(meta["patient"][tumor_m][0])
        # Guard: a section should be one patient.
        pats = set(meta["patient"][tumor_m].tolist())
        if len(pats) != 1:
            raise SystemExit(f"{sample}: multiple patients {pats}")
        core_genes = [meta["inventory"]["ifn_stat1_core_present"][j] for j in np.flatnonzero(core_use)]
        print(
            f"{sample} {patient}: tumor QC {tumor_m.sum()} "
            f"eligible50 {0 if graphs[50] is None else len(graphs[50]['elig'])} "
            f"core genes {len(core_genes)}",
            flush=True,
        )
        # Gene-level lognorm for the core+broad list (section-level Spearman later).
        gene_log = {}
        for name in list(dict.fromkeys(meta["inventory"]["ifn_broad_present"] + ["CLDN4"])):
            gene_log[name] = log_t[:, gix[name]].astype(np.float64)
        sections.append(
            {
                "sample": sample,
                "patient": patient,
                "histology": "LUSC" if "LUSC" in sample else "LUAD",
                "n_tumor_qc": int(tumor_m.sum()),
                "n_fov": int(len(big_fovs)),
                "fov": fov_t,
                "cldn4": cldn4,
                "stat1": stat1,
                "krt": krt,
                "logn": logn,
                "imm_frac": imm_frac,
                "core": core_score.astype(np.float64),
                "broad": None if broad_score is None else broad_score.astype(np.float64),
                "hallmark": None if hall_score is None else hall_score.astype(np.float64),
                "core_genes": core_genes,
                "n_broad_genes": int(broad_use.sum()) if broad_score is not None else 0,
                "n_hall_genes": int(hall_use.sum()) if hall_score is not None else 0,
                "detect_cldn4": float(detect[cldn_i]),
                "detect_stat1": float(detect[stat1_i]),
                "graphs": graphs,
                "gene_log": gene_log,
                "detect": {name: float(detect[gix[name]]) for name in gene_log},
            }
        )
    audit = {"median_nn_um_one_fov": med_nn_um, "n_tumor_labeled": int(tumor.sum())}
    return sections, audit


def metrics_for_score(sec: dict, score_name: str, radius: int, neigh: np.ndarray | None) -> dict:
    score = sec[score_name]
    if score is None:
        return {}
    krt = sec["krt"].reshape(-1, 1)
    out = {
        "same": spearman(sec["cldn4"], score),
        "same_krt": residual_spearman(sec["cldn4"], score, krt),
        "same_krt_lib": residual_spearman(
            sec["cldn4"], score, np.column_stack([sec["krt"], sec["logn"]])
        ),
        "same_imm": residual_spearman(sec["cldn4"], score, sec["imm_frac"].reshape(-1, 1)),
    }
    delta, rule = median_split_delta(sec["cldn4"], score)
    out["delta_high_minus_low"] = delta
    out["split_rule"] = rule
    if neigh is None:
        out["neigh"] = np.nan
        out["partial"] = np.nan
        out["neigh_krt"] = np.nan
        out["partial_krt"] = np.nan
        out["n_elig"] = 0
        return out
    elig = sec["graphs"][radius]["elig"]
    cldn = sec["cldn4"][elig]
    own = score[elig]
    krt_e = sec["krt"][elig].reshape(-1, 1)
    out["n_elig"] = int(len(elig))
    out["neigh"] = spearman(cldn, neigh)
    out["partial"] = partial_neighbor_rho(cldn, own, neigh, None)
    out["neigh_krt"] = residual_spearman(cldn, neigh, krt_e)
    out["partial_krt"] = partial_neighbor_rho(cldn, own, neigh, krt_e)
    out["same_on_elig"] = spearman(cldn, own)
    out["same_on_elig_krt"] = residual_spearman(cldn, own, krt_e)
    return out


def observed_tables(sections: list[dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    gene_rows = []
    for sec in sections:
        g50 = sec["graphs"][PRIMARY_RADIUS_UM]
        neigh50 = None if g50 is None else neighbor_means(sec["core"], g50["indices"], g50["indptr"])
        core_m = metrics_for_score(sec, "core", PRIMARY_RADIUS_UM, neigh50)
        # Other radii, core score only, unadjusted neighbor + partial.
        extra = {}
        for radius in RADII_UM:
            if radius == PRIMARY_RADIUS_UM:
                continue
            g = sec["graphs"][radius]
            if g is None:
                extra[f"neigh_{radius}"] = np.nan
                extra[f"partial_{radius}"] = np.nan
                continue
            nm = neighbor_means(sec["core"], g["indices"], g["indptr"])
            elig = g["elig"]
            extra[f"neigh_{radius}"] = spearman(sec["cldn4"][elig], nm)
            extra[f"partial_{radius}"] = partial_neighbor_rho(
                sec["cldn4"][elig], sec["core"][elig], nm, None
            )
        broad_m = {}
        hall_m = {}
        if sec["broad"] is not None and g50 is not None:
            nb = neighbor_means(sec["broad"], g50["indices"], g50["indptr"])
            broad_m = metrics_for_score(sec, "broad", PRIMARY_RADIUS_UM, nb)
        if sec["hallmark"] is not None and g50 is not None:
            nh = neighbor_means(sec["hallmark"], g50["indices"], g50["indptr"])
            hall_m = metrics_for_score(sec, "hallmark", PRIMARY_RADIUS_UM, nh)
        row = {
            "sample": sec["sample"],
            "patient": sec["patient"],
            "histology": sec["histology"],
            "n_tumor_qc": sec["n_tumor_qc"],
            "n_fov": sec["n_fov"],
            "n_core_genes": len(sec["core_genes"]),
            "core_genes": ",".join(sec["core_genes"]),
            "detect_cldn4": sec["detect_cldn4"],
            "detect_stat1": sec["detect_stat1"],
            "rho_stat1_same": spearman(sec["cldn4"], sec["stat1"]),
            "rho_stat1_same_krt": residual_spearman(
                sec["cldn4"], sec["stat1"], sec["krt"].reshape(-1, 1)
            ),
        }
        for key, val in core_m.items():
            row[f"core_{key}"] = val
        for key, val in extra.items():
            row[f"core_{key}"] = val
        for key in ("same", "same_krt", "neigh", "partial", "partial_krt", "delta_high_minus_low"):
            row[f"broad_{key}"] = broad_m.get(key, np.nan)
            row[f"hallmark_{key}"] = hall_m.get(key, np.nan)
        rows.append(row)
        for gene, vec in sec["gene_log"].items():
            if gene == "CLDN4":
                continue
            gene_rows.append(
                {
                    "sample": sec["sample"],
                    "patient": sec["patient"],
                    "gene": gene,
                    "detect": sec["detect"][gene],
                    "rho_same": spearman(sec["cldn4"], vec),
                    "rho_same_krt": residual_spearman(
                        sec["cldn4"], vec, sec["krt"].reshape(-1, 1)
                    ),
                    "in_core": gene in sec["core_genes"] or gene in IFN_STAT1_CORE,
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(gene_rows)


def summarize_patients(section_df: pd.DataFrame, value_col: str) -> dict:
    sub = section_df[["patient", "sample", value_col]].dropna()
    if sub.empty:
        return sign_test([])
    per_patient = sub.groupby("patient")[value_col].mean()
    per_section = sub[value_col].tolist()
    out = {
        "patient": sign_test(per_patient.tolist()),
        "section": sign_test(per_section),
        "patient_values": {k: float(v) for k, v in per_patient.items()},
        "wilcoxon_section_p_less": wilcoxon_less(per_section),
        "wilcoxon_patient_p_less": wilcoxon_less(per_patient.tolist()),
    }
    return out


def permute_spatial(sections: list[dict]) -> dict:
    """Within-FOV shuffle of the IFN/STAT1 score. CLDN4 stays put.

    Null for neighbor association. Same-cell pairing is broken by the shuffle
    of the score used to build neighbor means; own-IFN in the partial test
    stays the original score, so the partial null is spatial.
    """
    keys = ["neigh", "partial", "neigh_krt", "partial_krt"]
    # Observed patient-mean of each key, and nulls.
    observed_sec = {k: [] for k in keys}
    patients = []
    for sec in sections:
        g = sec["graphs"][PRIMARY_RADIUS_UM]
        nm = neighbor_means(sec["core"], g["indices"], g["indptr"])
        met = metrics_for_score(sec, "core", PRIMARY_RADIUS_UM, nm)
        patients.append(sec["patient"])
        for k in keys:
            observed_sec[k].append(met[k])
    observed = {k: patient_mean_from_sections(patients, observed_sec[k]) for k in keys}
    nulls = {k: np.empty(N_PERM, dtype=np.float64) for k in keys}
    for b in range(N_PERM):
        sec_vals = {k: [] for k in keys}
        pats = []
        for sec in sections:
            g = sec["graphs"][PRIMARY_RADIUS_UM]
            shuffled = shuffle_within_fov(sec["core"], sec["fov"], RNG)
            nm = neighbor_means(shuffled, g["indices"], g["indptr"])
            elig = g["elig"]
            cldn = sec["cldn4"][elig]
            own = sec["core"][elig]
            krt = sec["krt"][elig].reshape(-1, 1)
            sec_vals["neigh"].append(spearman(cldn, nm))
            sec_vals["partial"].append(partial_neighbor_rho(cldn, own, nm, None))
            sec_vals["neigh_krt"].append(residual_spearman(cldn, nm, krt))
            sec_vals["partial_krt"].append(partial_neighbor_rho(cldn, own, nm, krt))
            pats.append(sec["patient"])
        for k in keys:
            nulls[k][b] = patient_mean_from_sections(pats, sec_vals[k])
        if (b + 1) % 50 == 0 or b == 0:
            print(f"permutation {b + 1}/{N_PERM}", flush=True)
    out = {"n_perm": N_PERM, "null": "within_FOV_shuffle_of_IFNscore", "stats": {}}
    for k in keys:
        obs = observed[k]
        # One-sided: more negative than the spatial null.
        p = (1 + int(np.sum(nulls[k] <= obs))) / (N_PERM + 1)
        out["stats"][k] = {
            "observed_patient_mean_rho": obs,
            "null_mean": float(np.mean(nulls[k])),
            "null_p05": float(np.quantile(nulls[k], 0.05)),
            "perm_p_onesided_less": float(p),
        }
    return out


def plot_forest(section_df: pd.DataFrame, path: Path) -> None:
    # The shipped figure includes the count-quintile and library-size
    # controls from the companion scripts. Do not replace it with the
    # unadjusted-only forest on a re-run.
    if path.exists() and (OUT / "tables" / "neighbor_library.tsv").exists():
        print(f"keeping depth-controlled figure {path}", flush=True)
        return
    specs = [
        ("core_same", "Same cell", "#1f4e79"),
        ("core_neigh", "Tumor neighbors, 50 µm", "#c46b1a"),
        ("core_partial", "Neighbors | own IFN", "#5c4b8a"),
        ("core_same_krt", "Same cell | KRT8/18/19", "#4c78a8"),
        ("core_partial_krt", "Neighbors | own IFN + KRT", "#b279a2"),
    ]
    patients = ["Lung5", "Lung6", "Lung9", "Lung12", "Lung13"]
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    n_stats = len(specs)
    y_base = np.arange(len(patients))
    plotted = []
    for s, (col, label, color) in enumerate(specs):
        offset = (s - (n_stats - 1) / 2) * 0.13
        xs, ys = [], []
        for i, pat in enumerate(patients):
            sub = section_df.loc[section_df["patient"] == pat, col].dropna()
            if sub.empty:
                continue
            jitter = np.linspace(-0.03, 0.03, len(sub)) if len(sub) > 1 else np.array([0.0])
            ax.scatter(
                sub.to_numpy(),
                np.full(len(sub), i) + offset + jitter,
                s=18,
                color=color,
                alpha=0.45,
                zorder=2,
            )
            plotted.extend(sub.tolist())
            xs.append(float(sub.mean()))
            ys.append(i + offset)
        ax.scatter(xs, ys, s=46, color=color, marker="D", zorder=3, label=label)
        plotted.extend(xs)
    ax.axvline(0, color="#444444", lw=0.8)
    ax.set_yticks(y_base)
    ax.set_yticklabels(patients)
    ax.set_xlabel("Spearman ρ  (CLDN4 vs IFN/STAT1 module)")
    ax.set_title("Tumor cells, CosMx NSCLC (He 2022)")
    ax.legend(frameon=False, fontsize=8, loc="best")
    if plotted:
        lo, hi = min(plotted), max(plotted)
        pad = max(0.02, 0.08 * (hi - lo))
        ax.set_xlim(lo - pad, hi + pad)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_genes(gene_df: pd.DataFrame, core_genes: list[str], path: Path) -> None:
    sub = gene_df[gene_df["gene"].isin(core_genes)].copy()
    # Patient-mean rho.
    pat = sub.groupby(["patient", "gene"])["rho_same"].mean().reset_index()
    order = (
        pat.groupby("gene")["rho_same"].mean().sort_values().index.tolist()
    )
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    y = np.arange(len(order))
    for i, gene in enumerate(order):
        vals = pat.loc[pat["gene"] == gene, "rho_same"].to_numpy()
        ax.scatter(vals, np.full(len(vals), i), s=22, color="#9aa0a6", zorder=2)
        ax.scatter([vals.mean()], [i], s=36, color="#1f4e79", zorder=3)
    ax.axvline(0, color="#444444", lw=0.8)
    ax.set_yticks(y)
    labels = [g if g != "STAT1" else "STAT1" for g in order]
    ax.set_yticklabels(labels)
    for tick, gene in zip(ax.get_yticklabels(), order):
        if gene == "STAT1":
            tick.set_fontweight("bold")
    ax.set_xlabel("Same-cell Spearman ρ with CLDN4 (patient mean as diamond)")
    ax.set_title("IFN/STAT1 core genes in tumor cells")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fmt_rho(x: float) -> str:
    if not np.isfinite(x):
        return "NA"
    return f"{x:.3f}"


def fmt_p(x: float) -> str:
    if not np.isfinite(x):
        return "NA"
    if x < 0.001:
        return f"{x:.2e}"
    return f"{x:.3f}"


def write_results(section_df, gene_df, summaries, perm, inventory, audit) -> None:
    dest = OUT / "RESULTS.md"
    if dest.exists() and "library-size effect" in dest.read_text():
        print("keeping depth-controlled RESULTS.md", flush=True)
        return
    def line(col: str) -> str:
        s = summaries[col]
        pv = s["patient_values"]
        bits = ", ".join(
            f"{k} {fmt_rho(pv[k])}"
            for k in ["Lung5", "Lung6", "Lung9", "Lung12", "Lung13"]
            if k in pv
        )
        return (
            f"patient mean ρ {fmt_rho(s['patient']['mean'])} "
            f"({s['patient']['n_neg']}/{s['patient']['n']} negative, "
            f"one-sided sign P={fmt_p(s['patient']['sign_p_onesided'])}; "
            f"section {s['section']['n_neg']}/{s['section']['n']}, "
            f"sign P={fmt_p(s['section']['sign_p_onesided'])}; "
            f"section Wilcoxon P={fmt_p(s['wilcoxon_section_p_less'])}). "
            f"Patients: {bits}."
        )

    core_genes = inventory["ifn_stat1_core_present"]
    same = summaries["core_same"]["patient"]
    partial = summaries["core_partial"]["patient"]
    same_k = summaries["core_same_krt"]["patient"]
    part_k = summaries["core_partial_krt"]["patient"]
    perm_partial = perm["stats"]["partial"]["perm_p_onesided_less"]
    perm_part_k = perm["stats"]["partial_krt"]["perm_p_onesided_less"]
    perm_neigh = perm["stats"]["neigh"]["perm_p_onesided_less"]

    same_cell = same["n_neg"] == same["n"] and same["mean"] < 0
    neighbor_extra = (
        partial["n_neg"] == partial["n"] and partial["mean"] < 0 and perm_partial <= 0.05
    )
    if same_cell and neighbor_extra:
        where = "the same tumor cell and, beyond that, neighboring tumor cells"
    elif same_cell and not neighbor_extra:
        where = "the same tumor cell, not neighboring tumor cells beyond that cell"
    elif neighbor_extra and not same_cell:
        where = "neighboring tumor cells rather than the same cell"
    else:
        where = None
    if where:
        lead = (
            "Among tumor cells, CLDN4 is anti-correlated with the IFN/STAT1 module in "
            + where
            + "."
        )
    elif same["mean"] < 0 or partial["mean"] < 0:
        lead = (
            "CLDN4 and the IFN/STAT1 module do not show a 5/5 patient-level "
            "anti-correlation in the same tumor cell or in neighboring tumor cells "
            "beyond the same cell."
        )
    else:
        lead = (
            "Among tumor cells, CLDN4 is not anti-correlated with the IFN/STAT1 "
            "module in the same cell or in neighboring tumor cells."
        )
    if same_k["n_neg"] == same_k["n"] and same_k["mean"] < 0:
        krt_sentence = (
            "The same-cell anti-correlation remains after residualizing CLDN4 and "
            "the module on KRT8/KRT18/KRT19."
        )
    elif same["mean"] < 0 and same_k["mean"] > same["mean"]:
        krt_sentence = (
            "Residualizing CLDN4 and the module on KRT8/KRT18/KRT19 moves the "
            f"same-cell patient-mean ρ from {fmt_rho(same['mean'])} to {fmt_rho(same_k['mean'])}."
        )
    else:
        krt_sentence = (
            "Keratin residualization (KRT8/KRT18/KRT19) leaves the same-cell "
            f"patient-mean ρ at {fmt_rho(same_k['mean'])}."
        )
    if part_k["n_neg"] == part_k["n"] and part_k["mean"] < 0 and perm_part_k <= 0.05:
        krt_sentence += (
            " Neighbor IFN beyond own IFN also stays anti-correlated with CLDN4 after that residualization."
        )
    else:
        krt_sentence += (
            " Neighbor-tumor IFN beyond the cell's own IFN is not an additional "
            f"anti-correlation (patient-mean ρ {fmt_rho(partial['mean'])}, "
            f"keratin-adjusted {fmt_rho(part_k['mean'])})."
        )

    text = f"""# CosMx He 2022: CLDN4 vs IFN/STAT1 among tumor cells

Public CellCharter object of the He et al. 2022 CosMx SMI NSCLC 960-plex
(`cosmx_human_nsclc_clustered.h5ad`, figshare 25976224). 8 sections, 5 patients
(Lung5 ×3, Lung9 ×2, Lung6, Lung12, Lung13). Tumor cells are `cell_type`
labels that start with `tumor` (patient tumor clusters from the released
object). Normal `epithelial` cells are not in the test. CLDN4 only. No TACSTD2
gate. No private 8-KL.

This is the IFN side of an NHEJ–STING–IFN chain. It is not an NHEJ test.

## Panel limit (NHEJ and cGAS–STING)

Core NHEJ genes on the 960 panel: {inventory['nhej_core_present'] or 'none'}.
Absent: {', '.join(inventory['nhej_core_absent'])}.
Adjacent DDR genes present and not scored: {', '.join(inventory['nhej_adjacent_present']) or 'none'}.
cGAS/STING machinery present: {inventory['sting_machinery_present'] or 'none'}.
Absent: {', '.join(inventory['sting_machinery_absent'])}.
No NHEJ score and no cGAS–STING score are computed.

## Pre-specified IFN/STAT1 module

Core module (mean of within-section z-scores, genes detected in ≥2% of that
section's QC tumor cells): {', '.join(core_genes)}.
Queried but absent from the panel: {', '.join(inventory['ifn_stat1_core_absent'])}.
A broader score adds HLA/TAP/CIITA, JAK1/2, IFNGR1/2, and CCL5.
Hallmark interferon-gamma response contributes {len(inventory['hallmark_ifng_present'])} of {inventory['hallmark_ifng_n_gmt']} genes and is secondary, because that set includes immune-lineage markers.

Expression is log1p(CP10k) from the counts layer, including zeros. Cells with
n_counts < {MIN_COUNTS} are dropped. FOVs with < {MIN_FOV_TUMOR} QC tumor cells
are dropped. Neighbor means use other tumor cells inside the radius with at
least {MIN_NEIGHBORS} tumor neighbors. Primary radius is {PRIMARY_RADIUS_UM} µm
(20 and 100 µm are sensitivity). Pixel size {PX_TO_UM} µm; median nearest-neighbor
spacing in one checked FOV is {audit['median_nn_um_one_fov']:.1f} µm.

Inference is the patient (n=5). A one-sided sign test at n=5 reaches P=0.031
only at 5/5. Section Wilcoxon is reported for n=8. Spatial null: {N_PERM}
within-FOV shuffles of the IFN/STAT1 score, statistic = unweighted mean of
patient-mean ρ. CLDN4 is not shuffled. The partial test residualizes neighbor
IFN on the cell's own IFN before the Spearman with CLDN4.

## Result

{lead} {krt_sentence}

Same cell, core module: {line('core_same')}

Same cell, STAT1 alone: {line('rho_stat1_same')}

Same cell after KRT8/18/19 residualization: {line('core_same_krt')}

Same cell after KRT and log library size: {line('core_same_krt_lib')}

Tumor neighbors at 50 µm: {line('core_neigh')}
Spatial permutation P (one-sided, patient-mean ρ): {fmt_p(perm_neigh)}.
Null mean ρ {fmt_rho(perm['stats']['neigh']['null_mean'])}.

Neighbors after removing own IFN: {line('core_partial')}
Spatial permutation P: {fmt_p(perm_partial)}.
Null mean ρ {fmt_rho(perm['stats']['partial']['null_mean'])}.

Neighbors after removing own IFN and KRT8/18/19: {line('core_partial_krt')}
Spatial permutation P: {fmt_p(perm['stats']['partial_krt']['perm_p_onesided_less'])}.

Median-split Δ (CLDN4-high minus CLDN4-low, core module z): {line('core_delta_high_minus_low')}

20 µm neighbors: {line('core_neigh_20')}

100 µm neighbors: {line('core_neigh_100')}

20 µm partial (neighbors | own): {line('core_partial_20')}

100 µm partial (neighbors | own): {line('core_partial_100')}

Same cell residualized on local immune fraction (50 µm): {line('core_same_imm')}

Broad module, same cell: {line('broad_same')}

Broad module, neighbors | own: {line('broad_partial')}

Hallmark IFN-γ ∩ panel, same cell: {line('hallmark_same')}

Hallmark IFN-γ ∩ panel, neighbors | own: {line('hallmark_partial')}

## What this does not say

The locked CosMx exclusion result (fewer cytotoxic cells around CLDN4-high
tumor; 50/100 µm; 8/8 and 5/5) is not re-estimated here and is not replaced.
A correlation between CLDN4 and IFN in immune neighbors would restate that
exclusion, so immune-neighbor IFN is not the primary test. The primary
neighbor test is other tumor cells. Visium same-spot correlation is not used.
NHEJ is not supported or refuted: the panel does not contain it.

## Files

- `tables/section_correlations.tsv`
- `tables/gene_spearman.tsv`
- `tables/panel_inventory.json`
- `tables/summary.json`
- `figures/patient_forest.png`
- `figures/gene_spearman.png`
"""
    (OUT / "RESULTS.md").write_text(text)


def main() -> None:
    _self_check()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "tables").mkdir(exist_ok=True)
    (OUT / "figures").mkdir(exist_ok=True)
    print("loading counts", flush=True)
    meta = load_matrix()
    inv = meta["inventory"]
    (OUT / "tables" / "panel_inventory.json").write_text(json.dumps(inv, indent=2))
    print("NHEJ core present", inv["nhej_core_present"], flush=True)
    print("STING present", inv["sting_machinery_present"], flush=True)
    print("core IFN genes", inv["ifn_stat1_core_present"], flush=True)
    print("building sections", flush=True)
    sections, audit = prepare_sections(meta)
    del meta
    print("audit", audit, flush=True)
    print("correlations", flush=True)
    section_df, gene_df = observed_tables(sections)
    section_df.to_csv(OUT / "tables" / "section_correlations.tsv", sep="\t", index=False)
    gene_df.to_csv(OUT / "tables" / "gene_spearman.tsv", sep="\t", index=False)
    value_cols = [
        c
        for c in section_df.columns
        if c.startswith("core_") or c.startswith("broad_") or c.startswith("hallmark_") or c.startswith("rho_stat1")
    ]
    # Keep numeric rho-like columns only.
    skip_suffix = ("split_rule", "genes", "n_elig")
    value_cols = [c for c in value_cols if not c.endswith(skip_suffix) and pd.api.types.is_numeric_dtype(section_df[c])]
    summaries = {c: summarize_patients(section_df, c) for c in value_cols}
    print("permuting", flush=True)
    perm = permute_spatial(sections)
    summary = {"audit": audit, "summaries": summaries, "permutation": perm, "inventory": inv}
    (OUT / "tables" / "summary.json").write_text(json.dumps(summary, indent=2))
    plot_forest(section_df, OUT / "figures" / "patient_forest.png")
    plot_genes(gene_df, inv["ifn_stat1_core_present"], OUT / "figures" / "gene_spearman.png")
    write_results(section_df, gene_df, summaries, perm, inv, audit)
    print("wrote", OUT / "RESULTS.md", flush=True)


if __name__ == "__main__":
    main()
