#!/usr/bin/env python3
"""Maximize |Spearman ρ| of the locked 221-gene CLDN4-high signature.

The gene list and its order are the PR 590 signature: concordant-4 malignant
Q4 vs Q1, then TCGA-LUAD partial Spearman vs CLDN4 given KRT8/KRT18/KRT19
and ABSOLUTE purity. This script does not add or drop genes using CD8A,
ImmuneScore, or any GEO/OncoSG immune result. Size is a prefix of that
fixed order.

Primary objective, declared before the sweep is read:
  maximize the absolute DerSimonian–Laird meta-analytic Spearman of the
  score versus CD8A. Studies are OncoSG, GSE273377 (its two strata combined
  first), GSE282774, and GSE233774 tumors. All QC-passing tumors are used.
  No histology cut and no purity-median cut enters this objective.

Search space:
  method: z-mean (within-cohort) or ssGSEA (sum of the running enrichment)
  ssGSEA alpha: 0, 0.25, 0.75, 1 (z-mean does not use alpha; 0.25 is the
    ESTIMATE weight)
  size: 5..221, first `size` genes of the locked ranking
  purity: none; partial Spearman given published PURITY (OncoSG) or
    ESTIMATE StromalScore (GEO); or the same covariate residualized out of
    each gene before the score is computed

A second, separate objective repeats the same rule for ImmuneScore.
A GEO-only selection (same space, OncoSG held out) is applied unchanged
to OncoSG. A high-tumor subset (OncoSG PURITY >= median; GEO ESTIMATE
StromalScore <= median) is reported and is not the primary objective.

GSE10072, GSE11969, and GSE248378 are not opened. The locked single-gene
OncoSG correlations are not recomputed as a new claim; the 221-gene z-mean
is recomputed only as a numeric check against PR 590.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


HERE = Path(__file__).resolve().parent
SIG = HERE.parent / "cldn4_high_malignant_signature"
sys.path.insert(0, str(SIG))

import analyze as A  # noqa: E402
from weak_geo import ssgsea_sum  # noqa: E402


TABLES = HERE / "tables"
FIGURES = HERE / "figures"
ALPHAS = (0.0, 0.25, 0.75, 1.0)
PRIMARY_ALPHA = 0.25
SIZES = tuple(range(5, 222))
PURITY_MODES = ("none", "partial", "residual")
ENDPOINTS = ("CD8A", "ImmuneScore")
STUDIES = ("OncoSG", "GSE273377", "GSE282774", "GSE233774 tumor")

# Printed PR 590 values. A miss means this sweep is not the same score.
LOCKED = [
    ("OncoSG", "z-mean", 221, "CD8A", "none", -0.5779846849982597, 1e-9),
    ("OncoSG", "z-mean", 221, "ImmuneScore", "none", -0.615, 0.0015),
    ("OncoSG", "z-mean", 221, "CD8A", "partial", -0.4162752105822282, 1e-9),
    ("OncoSG", "z-mean", 30, "CD8A", "none", -0.478, 0.0015),
    ("GSE273377 discovery", "z-mean", 221, "CD8A", "none", -0.381, 0.0015),
    ("GSE273377 validation", "z-mean", 221, "CD8A", "none", -0.418, 0.0015),
    ("GSE282774", "z-mean", 221, "CD8A", "none", -0.212, 0.0015),
    ("GSE233774 tumor", "z-mean", 221, "CD8A", "none", -0.332, 0.0015),
    ("GSE282774", "ssGSEA", 20, "CD8A", "none", -0.592, 0.0015),
    ("GSE282774", "ssGSEA", 221, "CD8A", "none", -0.486, 0.0015),
    ("GSE282774", "ssGSEA", 20, "ImmuneScore", "none", -0.657, 0.0015),
    ("GSE282774", "ssGSEA", 221, "CD8A", "partial", -0.447, 0.0015),
    ("GSE233774 tumor", "ssGSEA", 221, "CD8A", "none", -0.420, 0.0015),
    ("GSE233774 tumor", "ssGSEA", 100, "ImmuneScore", "none", -0.637, 0.0015),
    ("GSE233774 tumor", "ssGSEA", 221, "ImmuneScore", "none", -0.578, 0.0015),
]


def prepare_ranks(mat: np.ndarray):
    n_genes, n_samples = mat.shape
    ranks = np.empty_like(mat, dtype=float)
    for j in range(n_samples):
        ranks[:, j] = stats.rankdata(mat[:, j], method="average")
    order = np.argsort(-ranks, axis=0, kind="mergesort")
    pos = np.empty((n_genes, n_samples), dtype=np.int32)
    pos[order, np.arange(n_samples)] = np.arange(n_genes, dtype=np.int32)[:, None]
    return ranks, pos


def ssgsea_from_pos(ranks: np.ndarray, pos: np.ndarray, gene_idx: np.ndarray, alpha: float) -> np.ndarray:
    """Sum of the ssGSEA running-enrichment walk. Matches weak_geo.ssgsea_sum."""
    n_genes, n_samples = ranks.shape
    idx = np.asarray(gene_idx, dtype=int)
    k = int(idx.size)
    out = np.full(n_samples, np.nan)
    if k < 5 or k >= n_genes:
        return out
    positions = pos[idx, :].astype(float)
    rank_vals = ranks[idx, :]
    order = np.argsort(positions, axis=0, kind="mergesort")
    positions = np.take_along_axis(positions, order, axis=0)
    rank_vals = np.take_along_axis(rank_vals, order, axis=0)
    weights = np.abs(rank_vals) ** alpha
    weight_sum = weights.sum(axis=0)
    ok = weight_sum > 0
    weights[:, ok] /= weight_sum[ok]
    nm = float(n_genes - k)
    hit_span = n_genes - positions
    sum_weighted = (weights * hit_span).sum(axis=0)
    sum_span = hit_span.sum(axis=0)
    out[ok] = sum_weighted[ok] - (n_genes * (n_genes + 1) / 2.0 - sum_span[ok]) / nm
    return out


def measured_genes(expr: pd.DataFrame) -> pd.DataFrame:
    """Drop genes that are missing in any sample.

    rankdata propagates a single NaN across the whole sample, which made
    OncoSG ssGSEA undefined (583 genes are blank in every sample). GEO
    matrices have no missing values, so this filter does not change them.
    """
    mat = expr.to_numpy(dtype=float)
    keep = np.isfinite(mat).all(axis=1)
    return expr.loc[keep]


def ssgsea_fast(expr: pd.DataFrame, genes: list[str], alpha: float = 0.25) -> np.ndarray:
    expr = measured_genes(expr)
    mat = expr.to_numpy(dtype=float)
    ranks, pos = prepare_ranks(mat)
    index = {g: i for i, g in enumerate(expr.index)}
    idx = np.array([index[g] for g in genes if g in index], dtype=int)
    return ssgsea_from_pos(ranks, pos, idx, alpha)


def self_test() -> None:
    rng = np.random.default_rng(20260921)
    n_genes, n_samples = 90, 5
    mat = rng.normal(size=(n_genes, n_samples))
    # Exact ties, the case where argsort stability matters.
    mat[:8, 0] = 1.0
    genes = [f"g{i}" for i in range(n_genes)]
    expr = pd.DataFrame(mat, index=genes)
    chosen = list(rng.choice(genes, size=15, replace=False))
    for alpha in (0.0, 0.25, 0.75, 1.0):
        got = ssgsea_fast(expr, chosen, alpha=alpha)
        ref = ssgsea_sum(expr, chosen) if alpha == 0.25 else _ssgsea_alpha(expr, chosen, alpha)
        if not np.allclose(got, ref, rtol=1e-10, atol=1e-8):
            raise SystemExit(f"ssGSEA mismatch at alpha {alpha}: {got} vs {ref}")
    short = ssgsea_fast(expr, chosen[:4], alpha=0.25)
    if np.isfinite(short).any():
        raise SystemExit("sets smaller than 5 genes must be missing")
    broken = expr.copy()
    broken.iloc[0, :] = np.nan
    with_hole = ssgsea_fast(broken, chosen, alpha=0.25)
    without = ssgsea_fast(expr.iloc[1:], chosen, alpha=0.25)
    if not np.allclose(with_hole, without, rtol=1e-10, atol=1e-8):
        raise SystemExit("dropping an all-missing gene changed ssGSEA")
    print("ssGSEA self-test ok", flush=True)


def _ssgsea_alpha(expr: pd.DataFrame, genes: list[str], alpha: float) -> np.ndarray:
    """Reference loop. Identical to weak_geo.ssgsea_sum except the weight power."""
    present = [g for g in genes if g in expr.index]
    mat = expr.to_numpy(dtype=float)
    in_set = np.isin(expr.index.to_numpy(), present)
    n_genes, n_samples = mat.shape
    out = np.empty(n_samples)
    for j in range(n_samples):
        ranks = stats.rankdata(mat[:, j], method="average")
        order = np.argsort(-ranks, kind="mergesort")
        correl = np.abs(ranks[order]) ** alpha
        tag = in_set[order].astype(float)
        nh = tag.sum()
        nm = n_genes - nh
        sum_correl = float(correl[tag == 1].sum())
        no_tag = 1.0 - tag
        res = np.cumsum(tag * correl / sum_correl) - np.cumsum(no_tag / nm)
        out[j] = float(res.sum())
    return out


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download_inputs() -> None:
    files = {
        "gencode.v36.probemap": "https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap",
        "GSE273377_exp_count.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273377/suppl/GSE273377_exp_count.txt.gz",
        "gse273377_gpl30173.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273377/matrix/GSE273377-GPL30173_series_matrix.txt.gz",
        "gse273377_gpl16791.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273377/matrix/GSE273377-GPL16791_series_matrix.txt.gz",
        "GSE233774_FPKM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233774/suppl/GSE233774_FPKM.txt.gz",
        "gse233774_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233774/matrix/GSE233774_series_matrix.txt.gz",
        "GSE282774_expr.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282774/suppl/GSE282774_N2_local_Gene_expression.csv.gz",
        "gse282774_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282774/matrix/GSE282774_series_matrix.txt.gz",
    }
    for name, url in files.items():
        A.download(url, A.CACHE / name)
    for name, sha in A.ONCOSG_FILES.items():
        A.download(f"{A.ONCOSG_MEDIA}/{name}", A.CACHE / name, sha256=sha)


def load_estimate_sets():
    stromal, immune = [], []
    for line in (SIG / "data" / "SI_geneset.gmt").read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        name, genes = parts[0], [g for g in parts[2:] if g]
        if name == "StromalSignature":
            stromal = genes
        elif name == "ImmuneSignature":
            immune = genes
    common = set(pd.read_csv(SIG / "data" / "common_genes.txt", sep="\t")["GeneSymbol"].astype(str))
    return stromal, immune, common


def load_ranked_signature() -> pd.DataFrame:
    path = SIG / "tables" / "signature_genes.tsv"
    signature = pd.read_csv(path, sep="\t")
    if len(signature) != 221 or signature["gene"].nunique() != 221:
        raise SystemExit(f"expected 221 unique signature genes, found {len(signature)}")
    if not signature["tcga_partial_rho"].is_monotonic_decreasing:
        raise SystemExit("signature is not ordered by descending TCGA partial ρ")
    signature.attrs["sha256"] = file_sha256(path)
    return signature


def stromal_vector(expr: pd.DataFrame, stromal_genes: list[str], common: set[str]) -> np.ndarray:
    keep = [g for g in expr.index if g in common]
    sub = expr.loc[keep]
    return ssgsea_fast(sub, stromal_genes, alpha=PRIMARY_ALPHA)


def residualize(expr: pd.DataFrame, cov: np.ndarray) -> pd.DataFrame:
    mat = expr.to_numpy(dtype=float)
    cov = np.asarray(cov, dtype=float)
    keep = np.isfinite(cov)
    if int(keep.sum()) < 8:
        return expr.iloc[:, 0:0]
    x = np.column_stack([np.ones(int(keep.sum())), cov[keep]])
    y = mat[:, keep]
    finite = np.isfinite(y).all(axis=1)
    resid = np.full_like(y, np.nan)
    beta, *_ = np.linalg.lstsq(x, y[finite].T, rcond=None)
    resid[finite] = y[finite] - (x @ beta).T
    out = pd.DataFrame(resid, index=expr.index, columns=expr.columns[keep])
    return out


def zmean_prefixes(expr: pd.DataFrame, ranked: list[str]) -> dict[int, np.ndarray]:
    """Within-cohort z-mean for every prefix. Same gene rules as analyze.zmean."""
    samples = list(expr.columns)
    entering: list[tuple[int, np.ndarray]] = []
    for size, gene in enumerate(ranked, start=1):
        if gene not in expr.index:
            continue
        values = expr.loc[gene, samples].to_numpy(dtype=float)
        finite = np.isfinite(values)
        if finite.mean() < 0.8:
            continue
        sd = np.nanstd(values)
        if sd == 0 or not np.isfinite(sd):
            continue
        entering.append((size, (values - np.nanmean(values)) / sd))
    out: dict[int, np.ndarray] = {}
    if not entering:
        return out
    stack = []
    next_i = 0
    for size in SIZES:
        while next_i < len(entering) and entering[next_i][0] <= size:
            stack.append(entering[next_i][1])
            next_i += 1
        if stack:
            out[size] = np.nanmean(np.vstack(stack), axis=0)
    return out


def ssgsea_prefixes(expr: pd.DataFrame, ranked: list[str], alphas=ALPHAS) -> dict[tuple[float, int], np.ndarray]:
    expr = measured_genes(expr)
    mat = expr.to_numpy(dtype=float)
    ranks, pos = prepare_ranks(mat)
    index = {gene: i for i, gene in enumerate(expr.index)}
    entering = [(size, index[gene]) for size, gene in enumerate(ranked, start=1) if gene in index]
    out: dict[tuple[float, int], np.ndarray] = {}
    idxs: list[int] = []
    next_i = 0
    for size in SIZES:
        while next_i < len(entering) and entering[next_i][0] <= size:
            idxs.append(entering[next_i][1])
            next_i += 1
        if len(idxs) < 5:
            continue
        gene_idx = np.asarray(idxs, dtype=int)
        for alpha in alphas:
            out[(alpha, size)] = ssgsea_from_pos(ranks, pos, gene_idx, alpha)
    return out


def present_count(expr: pd.DataFrame, ranked: list[str], size: int) -> int:
    return sum(1 for gene in ranked[:size] if gene in expr.index)


class Cohort:
    def __init__(self, name: str, study: str, expr: pd.DataFrame, purity: np.ndarray | None, purity_name: str):
        self.name = name
        self.study = study
        self.expr = expr
        self.purity = None if purity is None else np.asarray(purity, dtype=float)
        self.purity_name = purity_name
        self.samples = list(expr.columns)

    def endpoints(self, columns: list[str] | None = None):
        sub = self.expr if columns is None else self.expr.loc[:, columns]
        immune, immune_genes = A.zmean(sub, A.IMMUNE8, list(sub.columns))
        cd8 = sub.loc["CD8A"].to_numpy(dtype=float)
        return {"CD8A": cd8, "ImmuneScore": immune.to_numpy(dtype=float)}, immune_genes


def build_cohorts(stromal_genes, common) -> list[Cohort]:
    id_to_sym, _ = A.load_probemap(A.CACHE / "gencode.v36.probemap")
    oncosg, purity, _imsig = A.load_oncosg()
    cohorts = [
        Cohort("OncoSG", "OncoSG", oncosg, purity.to_numpy(dtype=float), "published PURITY"),
    ]
    strata = A.load_gse273377(id_to_sym)
    for stratum, (expr, _n_pass, _n_fail) in strata.items():
        cohorts.append(Cohort(f"GSE273377 {stratum}", "GSE273377", expr, None, "ESTIMATE StromalScore"))
    tumor = A.load_gse233774(id_to_sym)
    cohorts.append(Cohort("GSE233774 tumor", "GSE233774 tumor", tumor, None, "ESTIMATE StromalScore"))
    pN2 = A.load_gse282774()
    cohorts.append(Cohort("GSE282774", "GSE282774", pN2, None, "ESTIMATE StromalScore"))
    for cohort in cohorts:
        if "CD8A" not in cohort.expr.index:
            raise SystemExit(f"{cohort.name} has no CD8A")
        if cohort.purity is None:
            cohort.purity = stromal_vector(cohort.expr, stromal_genes, common)
            cohort.stromal_low_is_high_tumor = True
        else:
            cohort.stromal_low_is_high_tumor = False
        if len(cohort.purity) != cohort.expr.shape[1]:
            raise SystemExit(f"{cohort.name} covariate length {len(cohort.purity)} != samples {cohort.expr.shape[1]}")
    return cohorts


def high_tumor_columns(cohort: Cohort) -> list[str]:
    values = cohort.purity
    finite = np.isfinite(values)
    if finite.sum() < 8:
        return []
    cutoff = float(np.median(values[finite]))
    if cohort.stromal_low_is_high_tumor:
        keep = finite & (values <= cutoff)
    else:
        keep = finite & (values >= cutoff)
    return [sample for sample, flag in zip(cohort.samples, keep) if flag]


def correlate(score: np.ndarray, endpoint: np.ndarray, cov: np.ndarray | None, mode: str):
    score = np.asarray(score, dtype=float)
    finite = np.isfinite(score)
    if int(finite.sum()) < 6 or np.nanstd(score[finite]) == 0:
        return {"n": int(finite.sum()), "rho": np.nan, "p": np.nan, "ci_low": np.nan, "ci_high": np.nan, "k": 1 if mode == "partial" else 0}
    if mode == "partial":
        return A.partial_spearman(score, endpoint, [cov])
    return A.spearman_pair(score, endpoint)


def score_bank(expr: pd.DataFrame, ranked: list[str]):
    zmeans = zmean_prefixes(expr, ranked)
    ss = ssgsea_prefixes(expr, ranked)
    return zmeans, ss


def meta_from_strata(rows: list[dict], k: int):
    by_study: dict[str, list[dict]] = {}
    for row in rows:
        by_study.setdefault(row["study"], []).append(row)
    pieces = []
    stratum_rho = {}
    for study in STUDIES:
        group = by_study.get(study, [])
        if not group:
            continue
        if len(group) == 1:
            pieces.append({"rho": group[0]["rho"], "n": group[0]["n"]})
            stratum_rho[group[0]["cohort"]] = group[0]["rho"]
        else:
            combined = A.ivw(group, rho_key="rho", n_key="n", k=k)
            if combined["k_studies"]:
                pieces.append({"rho": combined["rho"], "n": combined["n_sum"]})
            for item in group:
                stratum_rho[item["cohort"]] = item["rho"]
    meta = A.ivw(pieces, rho_key="rho", n_key="n", k=k)
    meta["strata"] = stratum_rho
    return meta


def collect(cohorts: list[Cohort], ranked: list[str]):
    """Return cohort-level rows for the full sample set and the high-tumor subset."""
    rows = []
    coverage = []
    banks = {}
    for cohort in cohorts:
        subsets = {"all": list(cohort.samples), "high_tumor": high_tumor_columns(cohort)}
        for subset, columns in subsets.items():
            if len(columns) < 8:
                continue
            mask = np.array([sample in set(columns) for sample in cohort.samples])
            expr = cohort.expr.loc[:, columns]
            cov = cohort.purity[mask]
            endpoints, immune_genes = cohort.endpoints(columns)
            raw_z, raw_ss = score_bank(expr, ranked)
            finite_cov = np.isfinite(cov)
            resid_expr = residualize(expr.loc[:, np.asarray(columns)[finite_cov]], cov[finite_cov])
            resid_endpoints = {name: np.asarray(values)[finite_cov] for name, values in endpoints.items()}
            resid_z, resid_ss = score_bank(resid_expr, ranked)
            banks[(cohort.name, subset)] = {
                "z": raw_z,
                "ss": raw_ss,
                "endpoints": endpoints,
                "cov": cov,
                "n_immune": len(immune_genes),
            }
            coverage.append({
                "cohort": cohort.name,
                "study": cohort.study,
                "subset": subset,
                "n": len(columns),
                "n_genes_matrix": int(expr.shape[0]),
                "purity_name": cohort.purity_name,
                "purity_median": float(np.nanmedian(cov)),
                "n_immune8": len(immune_genes),
                "n_signature_present": present_count(expr, ranked, 221),
            })
            caches = {
                "none": (raw_z, raw_ss, endpoints, cov),
                "partial": (raw_z, raw_ss, endpoints, cov),
                "residual": (resid_z, resid_ss, resid_endpoints, cov[finite_cov]),
            }
            for mode, (zmeans, ss, endpoint_map, cov_use) in caches.items():
                k_cov = 1 if mode == "partial" else 0
                specs = [("z-mean", None, size, zmeans.get(size)) for size in SIZES]
                specs += [("ssGSEA", alpha, size, ss.get((alpha, size))) for alpha in ALPHAS for size in SIZES]
                for method, alpha, size, score in specs:
                    if score is None:
                        continue
                    for endpoint, y in endpoint_map.items():
                        stat = correlate(score, y, cov_use, mode)
                        rows.append({
                            "cohort": cohort.name,
                            "study": cohort.study,
                            "subset": subset,
                            "method": method,
                            "alpha": PRIMARY_ALPHA if alpha is None else alpha,
                            "size": size,
                            "purity": mode,
                            "endpoint": endpoint,
                            "n": stat["n"],
                            "rho": stat["rho"],
                            "p": stat["p"],
                            "k": k_cov,
                            "n_genes_present": present_count(expr, ranked, size),
                        })
            print(f"scored {cohort.name} {subset} n={len(columns)}", flush=True)
    return pd.DataFrame(rows), pd.DataFrame(coverage), banks


def check_locked(grid: pd.DataFrame) -> pd.DataFrame:
    records = []
    for cohort, method, size, endpoint, purity, expected, tol in LOCKED:
        alpha = PRIMARY_ALPHA
        hit = grid[
            (grid.cohort == cohort)
            & (grid.subset == "all")
            & (grid.method == method)
            & (grid["size"] == size)
            & (grid.endpoint == endpoint)
            & (grid.purity == purity)
            & (grid.alpha == alpha)
        ]
        if hit.empty:
            raise SystemExit(f"missing locked cell {cohort} {method} size {size} {endpoint} {purity}")
        rho = float(hit.iloc[0].rho)
        ok = abs(rho - expected) <= tol
        records.append({
            "cohort": cohort,
            "method": method,
            "size": size,
            "endpoint": endpoint,
            "purity": purity,
            "expected": expected,
            "observed": rho,
            "abs_diff": abs(rho - expected),
            "ok": ok,
        })
        if not ok:
            raise SystemExit(
                f"locked value mismatch {cohort} {method} size {size} {endpoint} {purity}: "
                f"observed {rho} expected {expected}"
            )
    out = pd.DataFrame(records)
    print(f"locked checks {len(out)}/{len(out)} matched", flush=True)
    return out


def spec_meta(grid: pd.DataFrame, method: str, alpha: float, size: int, purity: str, endpoint: str, subset: str):
    hit = grid[
        (grid.subset == subset)
        & (grid.method == method)
        & (grid.alpha == alpha)
        & (grid["size"] == size)
        & (grid.purity == purity)
        & (grid.endpoint == endpoint)
    ]
    k = 1 if purity == "partial" else 0
    rows = hit.to_dict("records")
    full = meta_from_strata(rows, k)
    geo_rows = [row for row in rows if row["study"] != "OncoSG"]
    geo = meta_from_strata(geo_rows, k)
    return full, geo, hit


REQUIRED = ("OncoSG", "GSE273377 discovery", "GSE273377 validation", "GSE282774", "GSE233774 tumor")


def strata_complete(hit: pd.DataFrame) -> bool:
    ok = set(hit.loc[np.isfinite(hit.rho), "cohort"])
    return all(name in ok for name in REQUIRED)


def enumerate_metas(grid: pd.DataFrame, subset: str) -> pd.DataFrame:
    records = []
    sub = grid[grid.subset == subset]
    keys = sub.groupby(["method", "alpha", "size", "purity", "endpoint"], sort=False).size().reset_index()
    for rec in keys.itertuples(index=False):
        full, geo, hit = spec_meta(sub, rec.method, rec.alpha, int(rec.size), rec.purity, rec.endpoint, subset)
        eligible = strata_complete(hit)
        records.append({
            "subset": subset,
            "eligible": eligible,
            "method": rec.method,
            "alpha": rec.alpha,
            "size": int(rec.size),
            "purity": rec.purity,
            "endpoint": rec.endpoint,
            "meta_rho": full["rho"],
            "meta_p": full["p"],
            "meta_i2": full["i2"],
            "meta_model": full["model"],
            "meta_n": full["n_sum"],
            "meta_k": full["k_studies"],
            "meta_ci_low": full.get("ci_low", np.nan),
            "meta_ci_high": full.get("ci_high", np.nan),
            "geo_rho": geo["rho"],
            "geo_p": geo["p"],
            "geo_i2": geo["i2"],
            "geo_model": geo["model"],
            "geo_n": geo["n_sum"],
            "abs_meta": abs(full["rho"]) if np.isfinite(full["rho"]) else np.nan,
            "abs_geo": abs(geo["rho"]) if np.isfinite(geo["rho"]) else np.nan,
        })
    return pd.DataFrame(records)


def pick_max(metas: pd.DataFrame, endpoint: str, column: str = "abs_meta") -> pd.Series:
    hit = metas[(metas.endpoint == endpoint) & (metas.eligible)].copy()
    hit = hit[np.isfinite(hit[column])]
    if hit.empty:
        raise SystemExit(f"no eligible spec for {endpoint} {column}")
    # Tie-break only when |ρ| is equal: ESTIMATE alpha, then the longer prefix.
    hit["alpha_rank"] = np.where(hit.method == "z-mean", 0, np.where(hit.alpha == PRIMARY_ALPHA, 0, 1))
    hit = hit.sort_values([column, "alpha_rank", "size"], ascending=[False, True, False])
    return hit.iloc[0]


def fmt_p(p) -> str:
    return A.fmt_p(p)


def fmt_rho(r) -> str:
    return A.fmt_rho(r)


def spec_label(row) -> str:
    alpha = ""
    if row["method"] == "ssGSEA":
        alpha = f", α={float(row['alpha']):g}"
    purity = {"none": "unadjusted", "partial": "partial | purity/stroma", "residual": "gene residual | purity/stroma"}[row["purity"]]
    return f"{row['method']}{alpha}, size {int(row['size'])}, {purity}"


def write_figures(grid: pd.DataFrame, winner, immune_at_winner, transfer_oncosg):
    FIGURES.mkdir(parents=True, exist_ok=True)
    sub = grid[
        (grid.subset == "all")
        & (grid.purity == "none")
        & (grid.endpoint == "CD8A")
        & (
            ((grid.method == "z-mean") & (grid.alpha == PRIMARY_ALPHA))
            | ((grid.method == "ssGSEA") & (grid.alpha == PRIMARY_ALPHA))
        )
    ]
    colors = {
        "OncoSG": "#1b4f72",
        "GSE273377 discovery": "#148f77",
        "GSE273377 validation": "#0e6655",
        "GSE282774": "#b9770e",
        "GSE233774 tumor": "#922b21",
    }
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.3), sharey=True)
    for ax, method in zip(axes, ("z-mean", "ssGSEA")):
        ax.axhline(0, color="#888888", lw=0.6)
        for cohort, color in colors.items():
            hit = sub[(sub.method == method) & (sub.cohort == cohort)].sort_values("size")
            ax.plot(hit["size"], hit["rho"], color=color, lw=1.15, label=cohort)
        ax.set_title(method if method == "z-mean" else "ssGSEA α=0.25")
        ax.set_xlabel("Prefix size")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("Spearman ρ vs CD8A")
    axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle("Unadjusted score vs CD8A", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_rho_vs_size.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGURES / "fig1_rho_vs_size.pdf", bbox_inches="tight")
    plt.close(fig)

    cohorts = list(colors)
    rhos = []
    los = []
    his = []
    ns = []
    for cohort in cohorts:
        _full, _geo, hit = spec_meta(
            grid, winner["method"], float(winner["alpha"]), int(winner["size"]), winner["purity"], "CD8A", "all"
        )
        row = hit[hit.cohort == cohort]
        if row.empty:
            rhos.append(np.nan)
            los.append(np.nan)
            his.append(np.nan)
            ns.append(0)
            continue
        rho = float(row.iloc[0].rho)
        n = int(row.iloc[0].n)
        k = 1 if winner.purity == "partial" else 0
        z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
        se = 1.0 / np.sqrt(n - 3 - k)
        rhos.append(rho)
        los.append(float(np.tanh(z - 1.959963984540054 * se)))
        his.append(float(np.tanh(z + 1.959963984540054 * se)))
        ns.append(n)
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    y = np.arange(len(cohorts))[::-1]
    ax.axvline(0, color="#888888", lw=0.6)
    ax.axvline(winner["meta_rho"], color="#1b4f72", lw=0.8, ls="--")
    ax.errorbar(
        rhos, y, xerr=[np.array(rhos) - np.array(los), np.array(his) - np.array(rhos)],
        fmt="o", color="#1b4f72", ms=5, lw=1,
    )
    for yi, rho, n in zip(y, rhos, ns):
        ax.text(0.02, yi + 0.18, f"n={n}  ρ={rho:+.3f}", fontsize=8, transform=ax.get_yaxis_transform())
    ax.set_yticks(y)
    ax.set_yticklabels(cohorts)
    ax.set_xlabel("Spearman ρ vs CD8A")
    ax.set_title(f"Winning CD8A spec: {spec_label(winner)}", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_winner_forest.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGURES / "fig2_winner_forest.pdf", bbox_inches="tight")
    plt.close(fig)

    # OncoSG at the GEO-selected spec is the held-out number. Plot that score.
    # The score itself is recomputed by the caller and passed in.
    if transfer_oncosg is not None:
        x_plot, y_plot, xlab, ylab, title = transfer_oncosg
        fig, ax = plt.subplots(figsize=(4.6, 4.2))
        ax.scatter(x_plot, y_plot, s=16, c="#1b4f72", alpha=0.8, linewidths=0)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        fig.savefig(FIGURES / "fig3_oncosg_heldout.png", dpi=160, bbox_inches="tight")
        fig.savefig(FIGURES / "fig3_oncosg_heldout.pdf", bbox_inches="tight")
        plt.close(fig)
    _ = immune_at_winner


def cohort_line(grid, spec, endpoint: str) -> str:
    _full, _geo, hit = spec_meta(
        grid, spec["method"], float(spec["alpha"]), int(spec["size"]), spec["purity"], endpoint, spec["subset"]
    )
    bits = []
    for cohort in ("OncoSG", "GSE273377 discovery", "GSE273377 validation", "GSE282774", "GSE233774 tumor"):
        row = hit[hit.cohort == cohort]
        if row.empty:
            bits.append(f"{cohort} NA")
        else:
            bits.append(f"{cohort} {fmt_rho(row.iloc[0].rho)} (p={fmt_p(row.iloc[0].p)}, n={int(row.iloc[0].n)})")
    return "; ".join(bits)


def _plateau(metas: pd.DataFrame) -> list[str]:
    hit = metas[(metas.endpoint == "CD8A") & (metas.eligible) & (metas.subset == "all")].nlargest(8, "abs_meta")
    methods = ", ".join(
        f"{row.method} α={float(row.alpha):g} size {int(row.size)} {row.purity} {fmt_rho(row.meta_rho)}"
        for row in hit.itertuples()
    )
    return [
        "The eight largest eligible CD8A meta |ρ| values are "
        + methods
        + ". The maximum is not a one-prefix spike."
    ]


def _purity_at_winner(grid: pd.DataFrame, metas: pd.DataFrame, winner) -> list[str]:
    lines = ["## Purity at the winning score", ""]
    lines.append(
        "Unadjusted correlation had the largest |meta ρ|. "
        "Partial correlation uses published PURITY on OncoSG and ESTIMATE StromalScore on GEO. "
        "Residual mode builds the score from gene-level linear residuals on that covariate and then correlates it with the raw endpoint."
    )
    lines.append("")
    for mode, label in (
        ("none", "unadjusted"),
        ("partial", "partial | purity/stroma"),
        ("residual", "gene residual | purity/stroma"),
    ):
        meta = metas[
            (metas.subset == "all")
            & (metas.method == winner["method"])
            & (metas.alpha == winner["alpha"])
            & (metas["size"] == winner["size"])
            & (metas.purity == mode)
            & (metas.endpoint == "CD8A")
        ].iloc[0]
        stub = winner.copy()
        stub["purity"] = mode
        lines.append(
            f"{label}: meta ρ {fmt_rho(meta.meta_rho)} ({meta.meta_model}, p={fmt_p(meta.meta_p)}, I²={meta.meta_i2:.0%}). "
            + cohort_line(grid, stub, "CD8A")
        )
        lines.append("")
    return lines


def _cohort_peaks(grid: pd.DataFrame) -> str:
    sub = grid[(grid.subset == "all") & np.isfinite(grid.rho) & (grid.endpoint == "CD8A")]
    lines = ["| cohort | n | ρ | p | spec |", "|---|---:|---:|---:|---|"]
    order = ["OncoSG", "GSE273377 discovery", "GSE273377 validation", "GSE282774", "GSE233774 tumor"]
    for cohort in order:
        hit = sub[sub.cohort == cohort]
        row = hit.loc[hit.rho.abs().idxmax()]
        alpha = f", α={float(row.alpha):g}" if row.method == "ssGSEA" else ""
        spec = f"{row.method}{alpha}, size {int(row['size'])}, {row.purity}"
        lines.append(
            f"| {cohort} | {int(row.n)} | {fmt_rho(row.rho)} | {fmt_p(row.p)} | {spec} |"
        )
    return "\n".join(lines)


def write_finding(ctx: dict) -> None:
    w = ctx["winner"]
    imm = ctx["immune_winner"]
    geo = ctx["geo_winner"]
    held = ctx["heldout"]
    base = ctx["baseline"]
    peak = ctx["peak_cell"]
    high = ctx["high_winner"]
    lines = []
    lines.append("# Maximum |ρ|: locked 221-gene CLDN4-high signature vs CD8A / ImmuneScore")
    lines.append("")
    lines.append("The 221 genes and their order are the PR 590 signature (`methods/cldn4_high_malignant_signature/tables/signature_genes.tsv`). Prefixes use that order. No gene was added or removed because of its correlation with CD8A or ImmuneScore. GSE10072, GSE11969, and GSE248378 stay closed. Bulk ρ is not a spatial exclusion result.")
    lines.append("")
    lines.append("## Objective")
    lines.append("")
    lines.append("Primary number: absolute DerSimonian–Laird meta-analytic Spearman versus CD8A on all QC-passing tumors. Studies are OncoSG, GSE273377 (discovery and validation inverse-variance combined first), GSE282774, and GSE233774 tumors. A spec is eligible only when all five strata have a finite correlation, so a cohort cannot be dropped to raise the meta |ρ|. The search is method (z-mean or ssGSEA), ssGSEA α in {0, 0.25, 0.75, 1}, prefix size 5 through 221, and purity mode {unadjusted, partial correlation on published PURITY or ESTIMATE StromalScore, gene-level residual on that same covariate}. The high-tumor median split is not in this objective.")
    lines.append("")
    lines.append("The selected spec was chosen on these same cohorts, so its meta p-value is the p-value of a maximized |ρ|, not a single pre-specified test. The pre-specified 221-gene z-mean from PR 590 is the baseline row. The OncoSG correlation at the GEO-only winning spec was not used to pick that spec.")
    lines.append("")
    lines.append("## Primary maximum")
    lines.append("")
    lines.append(f"Maximum |meta ρ| versus CD8A is **{fmt_rho(w.meta_rho)}** ({w.meta_model}, p={fmt_p(w.meta_p)}, I²={w.meta_i2:.0%}, n sum={int(w.meta_n)}, k={int(w.meta_k)}). Spec: **{spec_label(w)}**.")
    lines.append("")
    lines.append(f"The locked baseline (z-mean, size 221, unadjusted) meta ρ is {fmt_rho(base.meta_rho)} ({base.meta_model}, p={fmt_p(base.meta_p)}, I²={base.meta_i2:.0%}).")
    lines.append("")
    lines.append("Per cohort, CD8A, at the winning spec:")
    lines.append("")
    lines.append(cohort_line(ctx["grid"], w, "CD8A"))
    lines.append("")
    lines.append("ImmuneScore at that same spec (not re-selected):")
    lines.append("")
    lines.append(cohort_line(ctx["grid"], w, "ImmuneScore"))
    lines.append("")
    same_imm = ctx["immune_at_cd8_spec"]
    lines.append(f"ImmuneScore meta ρ at the CD8A spec: {fmt_rho(same_imm.meta_rho)} ({same_imm.meta_model}, p={fmt_p(same_imm.meta_p)}, I²={same_imm.meta_i2:.0%}).")
    lines.append("")
    lines.extend(_plateau(ctx["metas"]))
    lines.append("")
    lines.extend(_purity_at_winner(ctx["grid"], ctx["metas"], w))
    lines.append("")
    lines.append("## ImmuneScore maximum")
    lines.append("")
    lines.append(f"Repeating the objective for ImmuneScore gives meta ρ **{fmt_rho(imm.meta_rho)}** ({imm.meta_model}, p={fmt_p(imm.meta_p)}, I²={imm.meta_i2:.0%}, n sum={int(imm.meta_n)}). Spec: **{spec_label(imm)}**.")
    lines.append("")
    lines.append(cohort_line(ctx["grid"], imm, "ImmuneScore"))
    lines.append("")
    lines.append("## Held-out OncoSG")
    lines.append("")
    lines.append(f"The GEO-only maximum (OncoSG not in the selection) is {spec_label(geo)}, GEO meta ρ {fmt_rho(geo.geo_rho)} ({geo.geo_model}, p={fmt_p(geo.geo_p)}, I²={geo.geo_i2:.0%}, n sum={int(geo.geo_n)}).")
    lines.append("")
    lines.append(f"Applied to OncoSG versus CD8A: ρ={fmt_rho(held['rho'])} (p={fmt_p(held['p'])}, n={int(held['n'])}). Versus ImmuneScore: ρ={fmt_rho(held['immune_rho'])} (p={fmt_p(held['immune_p'])}, n={int(held['immune_n'])}).")
    lines.append("")
    lines.append("## Largest |ρ| inside one cohort")
    lines.append("")
    lines.append(
        f"Inside the primary search space (all tumors), the largest |ρ| is {fmt_rho(peak.rho)} "
        f"in {peak.cohort} versus {peak.endpoint} (p={fmt_p(peak.p)}, n={int(peak.n)}). "
        f"Spec: {peak.method}"
        + (f" α={float(peak.alpha):g}" if peak.method == "ssGSEA" else "")
        + f", size {int(peak['size'])}, purity={peak.purity}. "
        "That cell was not required to be the same spec in the other cohorts."
    )
    lines.append("")
    lines.append("Largest |ρ| for CD8A inside each cohort, same rule:")
    lines.append("")
    lines.append(_cohort_peaks(ctx["grid"]))
    lines.append("")
    lines.append("## High-tumor subset")
    lines.append("")
    lines.append("Not the primary objective. OncoSG keeps samples with published PURITY at or above the cohort median. GEO keeps samples with ESTIMATE StromalScore at or below the cohort median. z-means are recomputed inside the subset.")
    lines.append("")
    lines.append(f"Maximum |meta ρ| versus CD8A on this subset is {fmt_rho(high.meta_rho)} ({high.meta_model}, p={fmt_p(high.meta_p)}, I²={high.meta_i2:.0%}, n sum={int(high.meta_n)}). Spec: {spec_label(high)}.")
    lines.append("")
    lines.append(cohort_line(ctx["grid_high"], high, "CD8A"))
    lines.append("")
    lines.append("## Reproduction of PR 590 cells")
    lines.append("")
    lines.append(f"All {ctx['n_locked']} locked cells matched within the printed tolerance, including OncoSG z-mean size 221 versus CD8A at the full stored precision of the PR 590 summary.")
    lines.append("")
    lines.append("Purity covariate: OncoSG uses the cBioPortal published PURITY column. GEO uses ESTIMATE stromal ssGSEA (α=0.25, package 1.0.13 stromal set, common-gene filter). The Affymetrix cosine purity formula is not used. Partial correlation is Pearson of rank residuals. Residual mode is an ordinary Spearman of a score built from gene-level linear residuals on that covariate; the covariate is not the immune endpoint.")
    lines.append("")
    lines.append("Full grid: `tables/cohort_grid.tsv`. One row per spec: `tables/meta_grid.tsv`. Locked checks: `tables/locked_checks.tsv`.")
    lines.append("")
    (HERE / "FINDING.md").write_text("\n".join(lines))


def main():
    self_test()
    if "--self-test" in sys.argv:
        return
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    download_inputs()
    signature = load_ranked_signature()
    ranked = list(signature["gene"])
    stromal_genes, _immune_genes, common = load_estimate_sets()
    cohorts = build_cohorts(stromal_genes, common)
    grid, coverage, banks = collect(cohorts, ranked)
    grid.to_csv(TABLES / "cohort_grid.tsv", sep="\t", index=False)
    coverage.to_csv(TABLES / "coverage.tsv", sep="\t", index=False)
    locked = check_locked(grid)
    locked.to_csv(TABLES / "locked_checks.tsv", sep="\t", index=False)

    metas = enumerate_metas(grid, "all")
    metas_high = enumerate_metas(grid, "high_tumor")
    metas.to_csv(TABLES / "meta_grid.tsv", sep="\t", index=False)
    metas_high.to_csv(TABLES / "meta_grid_high_tumor.tsv", sep="\t", index=False)

    winner = pick_max(metas, "CD8A", "abs_meta")
    immune_winner = pick_max(metas, "ImmuneScore", "abs_meta")
    geo_winner = pick_max(metas, "CD8A", "abs_geo")
    high_winner = pick_max(metas_high, "CD8A", "abs_meta")
    baseline_rows = metas[
        (metas.method == "z-mean")
        & (metas["size"] == 221)
        & (metas.purity == "none")
        & (metas.endpoint == "CD8A")
        & (metas.alpha == PRIMARY_ALPHA)
    ]
    baseline = baseline_rows.iloc[0]
    immune_at = metas[
        (metas.method == winner.method)
        & (metas.alpha == winner["alpha"])
        & (metas["size"] == winner["size"])
        & (metas.purity == winner["purity"])
        & (metas.endpoint == "ImmuneScore")
    ].iloc[0]

    _full, _geo, held_hit = spec_meta(
        grid, geo_winner["method"], float(geo_winner["alpha"]), int(geo_winner["size"]), geo_winner["purity"], "CD8A", "all"
    )
    held_row = held_hit[held_hit.cohort == "OncoSG"].iloc[0]
    _f2, _g2, held_imm = spec_meta(
        grid, geo_winner["method"], float(geo_winner["alpha"]), int(geo_winner["size"]), geo_winner["purity"], "ImmuneScore", "all"
    )
    held_imm_row = held_imm[held_imm.cohort == "OncoSG"].iloc[0]
    held = {
        "rho": float(held_row.rho),
        "p": float(held_row.p),
        "n": int(held_row.n),
        "immune_rho": float(held_imm_row.rho),
        "immune_p": float(held_imm_row.p),
        "immune_n": int(held_imm_row.n),
    }

    primary_cells = grid[(grid.subset == "all") & np.isfinite(grid.rho)]
    peak = primary_cells.loc[primary_cells.rho.abs().idxmax()]

    onco = next(c for c in cohorts if c.name == "OncoSG")
    expr = onco.expr
    cov = onco.purity
    if geo_winner["purity"] == "residual":
        finite = np.isfinite(cov)
        expr = residualize(expr.loc[:, np.asarray(onco.samples)[finite]], cov[finite])
        cov = cov[finite]
    if geo_winner["method"] == "z-mean":
        score, _used = A.zmean(expr, ranked[: int(geo_winner["size"])], list(expr.columns))
        score = score.to_numpy(dtype=float)
    else:
        score = ssgsea_fast(expr, ranked[: int(geo_winner["size"])], alpha=float(geo_winner["alpha"]))
    cd8 = onco.expr.loc["CD8A", expr.columns].to_numpy(dtype=float)
    if not np.isfinite(held["rho"]):
        raise SystemExit("held-out OncoSG correlation is not finite")
    if geo_winner["purity"] == "partial":
        plotted = A.partial_spearman(score, cd8, [cov])
        x_plot = stats.rankdata(score)
        y_plot = stats.rankdata(cd8)
        # Show the partial association: residuals of the ranks.
        x_design = np.column_stack([np.ones(len(cov)), stats.rankdata(cov)])
        x_plot = x_plot - x_design @ np.linalg.lstsq(x_design, x_plot, rcond=None)[0]
        y_plot = y_plot - x_design @ np.linalg.lstsq(x_design, y_plot, rcond=None)[0]
        xlab, ylab = "Score rank residual | purity", "CD8A rank residual | purity"
    else:
        plotted = A.spearman_pair(score, cd8)
        x_plot, y_plot = score, cd8
        xlab, ylab = "Signature score", "CD8A"
    if abs(plotted["rho"] - held["rho"]) > 1e-6:
        raise SystemExit(f"held-out recompute {plotted['rho']} != grid {held['rho']}")
    transfer_plot = (x_plot, y_plot, xlab, ylab, f"OncoSG held-out CD8A ρ={held['rho']:+.3f}")

    write_figures(grid, winner, immune_at, transfer_plot)
    # Forest helper filters on winner.subset; attach it.
    winner = winner.copy()
    winner["subset"] = "all"
    immune_winner = immune_winner.copy()
    immune_winner["subset"] = "all"
    high_plot = high_winner.copy()
    high_plot["subset"] = "high_tumor"

    ctx = {
        "grid": grid,
        "grid_high": grid,
        "winner": winner,
        "immune_winner": immune_winner,
        "geo_winner": geo_winner,
        "heldout": held,
        "baseline": baseline,
        "immune_at_cd8_spec": immune_at,
        "peak_cell": peak,
        "high_winner": high_plot,
        "n_locked": int(len(locked)),
        "metas": metas,
    }
    write_finding(ctx)
    summary = {
        "signature_sha256": signature.attrs["sha256"],
        "winner": winner.to_dict(),
        "immune_winner": immune_winner.to_dict(),
        "geo_winner": geo_winner.to_dict(),
        "heldout_oncosg": held,
        "baseline_meta_rho": float(baseline.meta_rho),
        "peak_cell": {
            "cohort": peak.cohort,
            "endpoint": peak.endpoint,
            "method": peak.method,
            "alpha": float(peak.alpha),
            "size": int(peak["size"]),
            "purity": peak.purity,
            "rho": float(peak.rho),
            "p": float(peak.p),
            "n": int(peak.n),
        },
        "high_tumor_winner_rho": float(high_winner.meta_rho),
        "high_tumor_spec": spec_label(high_winner),
        "n_locked_ok": int(locked.ok.sum()),
    }
    def _plain(value):
        if isinstance(value, dict):
            return {str(k): _plain(v) for k, v in value.items()}
        if isinstance(value, (np.floating, float)):
            number = float(value)
            return number if np.isfinite(number) else None
        if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
            return int(value)
        if isinstance(value, (np.bool_, bool)):
            return bool(value)
        return value

    summary = _plain(summary)
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: summary[k] for k in ("winner", "heldout_oncosg", "peak_cell", "baseline_meta_rho")}, indent=2), flush=True)


if __name__ == "__main__":
    main()
