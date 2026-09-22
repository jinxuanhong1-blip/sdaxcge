#!/usr/bin/env python3
"""Maximize the TCGA TACSTD2 versus CD8/CD3 association after keratin.

The grid below is the specification family. It is fixed in this file before
the cohort matrices are scored.

Question
--------
In eight TCGA primary-tumor cohorts (LUAD, LUSC, BRCA, CESC, KIRC, STAD,
BLCA, PAAD), how negative is TACSTD2 versus a CD8 or CD3 score once keratin
expression is in the model?

Immune scores
-------------
Single genes CD3D, CD3E, CD3G, CD8A, CD8B.
Arithmetic mean, within-sample z-mean, rank-mean, and first principal
component of CD3 (CD3D/E/G), CD8 (CD8A/B), and CD3+CD8.
No cytotoxic-only, interferon, or extended T-cell signature is in the grid.

Keratin covariates
------------------
Basal: KRT5, KRT6A, KRT6B, KRT14, the broader basal set, and that basal set
plus TP63.
Simple epithelial: KRT8, KRT18, KRT19, with KRT7 and with KRT7/KRT20.
Summaries: z-mean and PC1 of the basal four, of KRT8/18/19, and of the full
keratin list.
Both programs together: the basal four plus KRT8/18/19.
Optional extra column: Yoshihara ESTIMATE purity on the basal four, on
KRT8/18/19, or beside the full keratin z-mean.
`stratum_only` has no covariate. It is scored only inside a keratin median
or tertile split, and it is excluded from the keratin-adjusted maximum.

Strata, cut inside each cohort
------------------------------
All patients.
Basal score = mean of KRT5, KRT6A, KRT6B: at/above median, below median,
top tertile, bottom tertile.
Simple score = mean of KRT8, KRT18, KRT19: the same four splits.
ESTIMATE purity: at/above median, below median.

Methods and tail cuts
---------------------
Partial Spearman, partial Pearson, and Pearson after a 1%/99% winsorization.
Cliff's delta compares the outer 25%, 20%, and 10% of TACSTD2. The immune
score is the keratin residual. The TACSTD2 cut is either raw expression or
that same residual. Arms smaller than 12 patients are untested and are not
counted as negative.

Eligibility and objective
-------------------------
Keratin-adjusted: at least one covariate (k >= 1), random-effects pooled
effect negative, at least 6 of 8 cohorts negative, at least 6 cohorts in
the pool.
Continuous objective: maximize |pooled Spearman/Pearson ρ|.
Cliff objective: maximize |pooled Cliff's δ|.
Tie-break: more negative cohorts, then more cohorts in the pool, then
smaller meta p.
Correlations are pooled by DerSimonian-Laird on the Fisher z scale.
Deltas are pooled by DerSimonian-Laird on the delta scale with the Cliff
1993 variance.

This maximum is the extreme of the search. Selecting it inflates |effect|
relative to any one pre-specified test. Benjamini-Hochberg q-values use the
meta p of every specification in that family, including specifications that
are not eligible.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import sys
import time
import urllib.request
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats import (  # noqa: E402
    bh_fdr,
    cliffs_batch,
    dl_random_effects,
    fisher_ci,
    fisher_random_effects,
    matrix_partial_corr,
    residualize,
    select_abs,
    select_positive,
)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "results", "tcga_tacstd2_keratin_cd8")
TABLES = os.path.join(OUT, "tables")
FIGS = os.path.join(OUT, "figures")
CACHE = os.environ.get("TCGA_TACSTD2_CACHE", "/tmp/tcga_tacstd2")
ART = os.environ.get("TCGA_TACSTD2_ART", "/opt/cursor/artifacts/tcga_tacstd2")
GDC = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"
PROBEMAP_URL = f"{GDC}/gencode.v36.annotation.gtf.gene.probemap"
MDACC_BASE = "https://bioinformatics.mdanderson.org/estimate/tables"
ESTIMATE_A = 0.6049872018
ESTIMATE_B = 0.0001467884

COHORTS = ["LUAD", "LUSC", "BRCA", "CESC", "KIRC", "STAD", "BLCA", "PAAD"]
MIN_N = 40
MIN_ARM = 12
CUTS = (0.25, 0.20, 0.10)
METHODS = ("spearman", "pearson", "pearson_winsor")
MIN_NEG = 6

BASAL4 = ["KRT5", "KRT6A", "KRT6B", "KRT14"]
SIMPLE3 = ["KRT8", "KRT18", "KRT19"]
BASAL_BROAD = ["KRT5", "KRT6A", "KRT6B", "KRT6C", "KRT14", "KRT16", "KRT17"]
ALL_KRT = [
    "KRT5", "KRT6A", "KRT6B", "KRT6C", "KRT7", "KRT8", "KRT13", "KRT14",
    "KRT16", "KRT17", "KRT18", "KRT19", "KRT20",
]
# (name, genes, kind, needs_estimate). kind: columns, z, pc1, z_plus, none.
COVARIATE_SETS = [
    ("KRT5_KRT6A", ["KRT5", "KRT6A"], "columns", False),
    ("KRT5_6AB", ["KRT5", "KRT6A", "KRT6B"], "columns", False),
    ("KRT5_6_14", BASAL4, "columns", False),
    ("basal_broad", BASAL_BROAD, "columns", False),
    ("KRT5_6_14_TP63", BASAL4 + ["TP63"], "columns", False),
    ("basal_z", BASAL4, "z", False),
    ("basal_PC1", BASAL4, "pc1", False),
    ("KRT8_18_19", SIMPLE3, "columns", False),
    ("KRT8_18_19_7", SIMPLE3 + ["KRT7"], "columns", False),
    ("KRT7_8_18_19_20", ["KRT7", "KRT8", "KRT18", "KRT19", "KRT20"], "columns", False),
    ("simple_z", SIMPLE3, "z", False),
    ("simple_PC1", SIMPLE3, "pc1", False),
    ("both_programs", BASAL4 + SIMPLE3, "columns", False),
    ("keratin_z", ALL_KRT, "z", False),
    ("keratin_PC1", ALL_KRT, "pc1", False),
    ("KRT5_6_14_ESTIMATE", BASAL4 + ["ESTIMATE"], "columns", True),
    ("KRT8_18_19_ESTIMATE", SIMPLE3 + ["ESTIMATE"], "columns", True),
    ("keratin_z_ESTIMATE", ALL_KRT + ["ESTIMATE"], "z_plus", True),
    ("stratum_only", [], "none", False),
]
GENE_SETS = {
    "CD3": ["CD3D", "CD3E", "CD3G"],
    "CD8": ["CD8A", "CD8B"],
    "CD3CD8": ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B"],
}
SINGLES = ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B"]
PANEL_OUTCOMES = ["CD8A", "CD8_mean", "CD3D", "CD3_mean", "CD3CD8_z"]
PANEL_COVARS = ["KRT5_6_14", "KRT8_18_19", "both_programs", "keratin_z", "simple_z", "basal_z"]
KRT_STRATA_PREFIX = ("basal_", "simple_")

MDACC_ESTIMATE = {
    "LUAD": "lung_adenocarcinoma_RNAseqV2.txt",
    "LUSC": "lung_squamous_cell_carcinoma_RNAseqV2.txt",
    "BRCA": "breast_cancer_RNAseqV2.txt",
    "CESC": "cervical_carcinoma_RNAseqV2.txt",
    "KIRC": "kidney_renal_clear_cell_carcinoma_RNAseqV2.txt",
    "STAD": "stomach_adenocarcinoma_RNAseqV2.txt",
    "BLCA": "bladder_urothelial_carcinoma_RNAseqV2.txt",
    "PAAD": "pancreatic_ductal_adenocarcinoma_RNAseqV2.txt",
}

# Prior partial Spearman, all patients, KRT5+KRT6A+KRT6B+KRT14, TACSTD2 only.
# Six-decimal values from the KRT5/6 CD8/CD3 grid.
ANCHOR_RHO = {
    "CD8_mean": {
        "BLCA": -0.218597, "BRCA": -0.147181, "CESC": 0.080593, "KIRC": -0.074849,
        "LUAD": -0.142577, "LUSC": -0.141121, "PAAD": -0.260992, "STAD": -0.039386,
    },
    "CD3_mean": {
        "BLCA": -0.159745, "BRCA": -0.149325, "CESC": 0.079720, "KIRC": -0.043927,
        "LUAD": -0.104046, "LUSC": -0.094443, "PAAD": -0.250227, "STAD": -0.106670,
    },
}

DATA = {}


def symbols_needed():
    genes = {"TACSTD2", "TP63"}
    genes.update(SINGLES)
    genes.update(ALL_KRT)
    for items in GENE_SETS.values():
        genes.update(items)
    for _, items, _, _ in COVARIATE_SETS:
        genes.update(g for g in items if g != "ESTIMATE")
    return sorted(genes)


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url, dest, retries=4):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "tcga-tacstd2-keratin/1.0"})
            with urllib.request.urlopen(req, timeout=300) as resp, open(tmp, "wb") as out:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
            os.replace(tmp, dest)
            return
        except Exception as exc:  # noqa: BLE001
            last = exc
            if os.path.exists(tmp):
                os.remove(tmp)
            print(f"[retry {attempt}] {url} :: {exc}", flush=True)
            time.sleep(min(8, 2 ** (attempt - 1)))
    raise RuntimeError(f"failed to download {url}: {last}")


def load_probemap(symbols):
    path = os.path.join(CACHE, "gencode.v36.annotation.gtf.gene.probemap")
    download(PROBEMAP_URL, path)
    table = pd.read_csv(path, sep="\t")
    mapping = {}
    for symbol in symbols:
        sub = table[table["gene"] == symbol].copy()
        if sub.empty:
            raise SystemExit(f"probemap missing {symbol}")
        sub = sub[sub["chrom"].astype(str).str.match(r"^chr([0-9]+|X|Y)$")]
        if sub.empty:
            raise SystemExit(f"no primary-chrom id for {symbol}")
        sub["span"] = sub["chromEnd"] - sub["chromStart"]
        sub = sub.sort_values(["span", "id"], ascending=[False, True])
        mapping[symbol] = str(sub.iloc[0]["id"])
    if len(set(mapping.values())) != len(mapping):
        raise SystemExit(f"Ensembl id collision: {mapping}")
    return mapping


def patient_id(barcode):
    parts = str(barcode).replace(".", "-").split("-")
    if len(parts) < 4 or not parts[3].startswith("01"):
        return None
    return "-".join(parts[:3])


def extract_cohort(cohort, id_to_gene):
    cache_path = os.path.join(CACHE, f"{cohort}.primary01.tsv.gz")
    url = f"{GDC}/TCGA-{cohort}.star_tpm.tsv.gz"
    if os.path.exists(cache_path) and os.path.getsize(cache_path) > 0:
        table = pd.read_csv(cache_path, sep="\t")
        return {
            "cohort": cohort,
            "url": url,
            "cache": cache_path,
            "sha256": sha256_file(cache_path),
            "n_patients": int(table.shape[0]),
            "from_cache": True,
        }
    want = set(id_to_gene)
    print(f"[get] {cohort} {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "tcga-tacstd2-keratin/1.0"})
    collected = {}
    with urllib.request.urlopen(req, timeout=300) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            header = gz.readline().decode("utf-8").rstrip("\n").split("\t")
            samples = header[1:]
            keep_idx = []
            patients = []
            for i, sample in enumerate(samples):
                pid = patient_id(sample)
                if pid is None:
                    continue
                keep_idx.append(i)
                patients.append(pid)
            if len(keep_idx) < 15:
                raise RuntimeError(f"{cohort}: only {len(keep_idx)} primary-tumor columns")
            keep_idx_arr = np.asarray(keep_idx, dtype=int)
            for raw in gz:
                gid, rest = raw.decode("utf-8").split("\t", 1)
                if gid not in want:
                    continue
                values = np.fromstring(rest, sep="\t", dtype=float)
                if values.size != len(samples):
                    parts = rest.rstrip("\n").split("\t")
                    values = np.array([float(x) if x else np.nan for x in parts], dtype=float)
                collected[id_to_gene[gid]] = values[keep_idx_arr]
    missing = sorted(set(id_to_gene.values()) - set(collected))
    if missing:
        raise RuntimeError(f"{cohort} missing genes: {missing}")
    frame = pd.DataFrame(collected)
    frame.insert(0, "patient", patients)
    collapsed = frame.groupby("patient", as_index=False).mean(numeric_only=True)
    os.makedirs(CACHE, exist_ok=True)
    collapsed.to_csv(cache_path, sep="\t", index=False, compression="gzip")
    return {
        "cohort": cohort,
        "url": url,
        "cache": cache_path,
        "sha256": sha256_file(cache_path),
        "n_patients": int(collapsed.shape[0]),
        "from_cache": False,
    }


def estimate_purity_from_score(score):
    purity = np.cos(ESTIMATE_A + ESTIMATE_B * pd.to_numeric(score, errors="coerce"))
    return pd.Series(np.clip(purity, 0.0, 1.0), index=getattr(score, "index", None))


def load_purity():
    rows = []
    metas = []
    for cohort, fname in MDACC_ESTIMATE.items():
        path = os.path.join(CACHE, "estimate", fname)
        download(f"{MDACC_BASE}/{fname}", path)
        table = pd.read_csv(path, sep="\t")
        table["patient"] = table["ID"].map(lambda s: patient_id(str(s)))
        table = table.dropna(subset=["patient"])
        table["ESTIMATE"] = estimate_purity_from_score(table["ESTIMATE_score"])
        table["cohort"] = cohort
        rows.append(table[["cohort", "patient", "ESTIMATE"]])
        metas.append({"cohort": cohort, "file": fname, "sha256": sha256_file(path), "n": int(len(table))})
    estimate = pd.concat(rows, ignore_index=True)
    estimate = estimate.groupby(["cohort", "patient"], as_index=False)["ESTIMATE"].mean()
    return estimate, metas


def k_of(kind, genes):
    if kind == "none":
        return 0
    if kind in ("z", "pc1"):
        return 1
    if kind == "z_plus":
        return 2
    return len(genes)


def zmean(mat):
    sd = mat.std(axis=0, ddof=1)
    keep = sd > 1e-8
    if not np.any(keep):
        return np.full(mat.shape[0], np.nan)
    z = (mat[:, keep] - mat[:, keep].mean(axis=0)) / sd[keep]
    return z.mean(axis=1)


def rankmean(mat):
    ranks = np.column_stack([stats.rankdata(mat[:, j]) for j in range(mat.shape[1])])
    return ranks.mean(axis=1)


def first_pc(mat, orient):
    sd = mat.std(axis=0, ddof=1)
    keep = sd > 1e-8
    if keep.sum() < 2:
        return np.full(mat.shape[0], np.nan)
    z = (mat[:, keep] - mat[:, keep].mean(axis=0)) / sd[keep]
    _, _, vt = np.linalg.svd(z, full_matrices=False)
    scores = z @ vt[0]
    o = np.asarray(orient, dtype=float)
    if np.isfinite(o).all() and np.std(o) > 1e-12 and np.std(scores) > 1e-12:
        if np.corrcoef(scores, o)[0, 1] < 0:
            scores = -scores
    return scores


def scale_matrix(mat, method):
    mat = np.asarray(mat, dtype=float)
    single = mat.ndim == 1
    if single:
        mat = mat[:, None]
    if method == "spearman":
        out = np.column_stack([stats.rankdata(mat[:, j]).astype(float) for j in range(mat.shape[1])])
    elif method == "pearson_winsor":
        lo = np.quantile(mat, 0.01, axis=0)
        hi = np.quantile(mat, 0.99, axis=0)
        out = np.clip(mat, lo, hi)
    elif method == "pearson":
        out = np.array(mat, dtype=float, copy=True)
    else:
        raise ValueError(method)
    return out[:, 0] if single else out


def build_outcomes(columns):
    names = []
    cols = []
    n = len(columns["TACSTD2"])
    for gene in SINGLES:
        names.append(gene)
        cols.append(columns[gene])
    for set_name, genes in GENE_SETS.items():
        mat = np.column_stack([columns[g] for g in genes])
        orient_gene = next(g for g in ("CD3D", "CD3E", "CD3G", "CD8A", "CD8B") if g in genes)
        names.append(f"{set_name}_mean")
        cols.append(mat.mean(axis=1))
        names.append(f"{set_name}_z")
        cols.append(zmean(mat))
        names.append(f"{set_name}_rank")
        cols.append(rankmean(mat))
        names.append(f"{set_name}_pc1")
        cols.append(first_pc(mat, columns[orient_gene]))
    out = np.column_stack(cols)
    if out.shape != (n, len(names)):
        raise RuntimeError("outcome matrix shape mismatch")
    return names, out


def covariate_vectors(columns, name, genes, kind, method):
    if kind == "none":
        return []
    if kind == "pc1":
        mat = np.column_stack([columns[g] for g in genes])
        if name == "basal_PC1":
            orient = columns["KRT5"]
        elif name == "simple_PC1":
            orient = columns["KRT8"]
        else:
            orient = 0.5 * (columns["KRT5"] + columns["KRT8"])
        score = first_pc(mat, orient)
        return [scale_matrix(score, method)]
    if kind == "z":
        mat = np.column_stack([columns[g] for g in genes])
        return [scale_matrix(zmean(mat), method)]
    if kind == "z_plus":
        krt = [g for g in genes if g != "ESTIMATE"]
        mat = np.column_stack([columns[g] for g in krt])
        score = scale_matrix(zmean(mat), method)
        est = scale_matrix(columns["ESTIMATE"], method)
        return [score, est]
    mat = np.column_stack([columns[g] for g in genes])
    scaled = scale_matrix(mat, method)
    return [scaled[:, j] for j in range(scaled.shape[1])]


def stratum_masks(df):
    basal = df[["KRT5", "KRT6A", "KRT6B"]].mean(axis=1).to_numpy(dtype=float)
    simple = df[["KRT8", "KRT18", "KRT19"]].mean(axis=1).to_numpy(dtype=float)
    masks = {"all": np.ones(len(df), dtype=bool)}

    def add_splits(prefix, score):
        median = np.median(score)
        q33, q66 = np.quantile(score, [1.0 / 3.0, 2.0 / 3.0])
        masks[f"{prefix}_high"] = score >= median
        masks[f"{prefix}_low"] = score < median
        masks[f"{prefix}_top_tertile"] = score >= q66
        masks[f"{prefix}_bot_tertile"] = score <= q33

    add_splits("basal", basal)
    add_splits("simple", simple)
    purity = df["ESTIMATE"].to_numpy(dtype=float)
    ok = np.isfinite(purity)
    if int(ok.sum()) >= 80:
        pmed = np.median(purity[ok])
        masks["pur_high"] = ok & (purity >= pmed)
        masks["pur_low"] = ok & (purity < pmed)
    return masks


def groups(x, cut):
    q_lo = np.quantile(x, cut)
    q_hi = np.quantile(x, 1.0 - cut)
    low = x <= q_lo
    high = x >= q_hi
    both = low & high
    return low & ~both, high & ~both


def evaluate_block(cohort, stratum, block, cov_items, bucket):
    columns = {
        col: block[col].to_numpy(dtype=float)
        for col in block.columns
        if col != "patient"
    }
    n = len(block)
    if n < MIN_N:
        return
    outcome_names, outcome_raw = build_outcomes(columns)
    pred_raw = columns["TACSTD2"]
    n_out = len(outcome_names)
    for method in METHODS:
        pred_s = scale_matrix(pred_raw, method)
        out_s = scale_matrix(outcome_raw, method)
        for cov_name, genes, kind, _needs in cov_items:
            if kind == "none" and not stratum.startswith(KRT_STRATA_PREFIX):
                continue
            covs = covariate_vectors(columns, cov_name, genes, kind, method)
            k = k_of(kind, genes)
            if n < k + 5:
                continue
            rho, pval = matrix_partial_corr(pred_s[:, None], out_s, covs, n - 2 - k)
            for i, outcome in enumerate(outcome_names):
                bucket.append((
                    "continuous", cohort, stratum, cov_name, method, outcome,
                    "", "", float(rho[0, i]), float(pval[0, i]), np.nan,
                    int(k), int(n), 0, 0,
                ))
            pred_resid = residualize(pred_s, covs)
            out_resid = residualize(out_s, covs)
            for cut in CUTS:
                cut_label = f"{cut:.2f}"
                for mode, xvec in (("expression", pred_raw), ("residual", pred_resid)):
                    low, high = groups(np.asarray(xvec, dtype=float), cut)
                    n_low = int(low.sum())
                    n_high = int(high.sum())
                    if n_low < MIN_ARM or n_high < MIN_ARM:
                        delta = np.full(n_out, np.nan)
                        p_nan = np.full(n_out, np.nan)
                        v_nan = np.full(n_out, np.nan)
                    else:
                        delta, p_nan, v_nan = cliffs_batch(out_resid[high], out_resid[low])
                    for i, outcome in enumerate(outcome_names):
                        bucket.append((
                            "quantile", cohort, stratum, cov_name, method, outcome,
                            cut_label, mode, float(delta[i]), float(p_nan[i]), float(v_nan[i]),
                            int(k), int(n), int(n_high), int(n_low),
                        ))


def evaluate_cohort(cohort):
    frame = DATA[cohort]
    bucket = []
    masks = stratum_masks(frame)
    plain = [item for item in COVARIATE_SETS if not item[3]]
    with_purity = [item for item in COVARIATE_SETS if item[3]]
    for stratum, mask in masks.items():
        sub = frame.loc[mask]
        if len(sub) < MIN_N:
            print(f"[sweep] {cohort} {stratum} skip n={len(sub)}", flush=True)
            continue
        evaluate_block(cohort, stratum, sub, plain, bucket)
        pure = sub[np.isfinite(sub["ESTIMATE"].to_numpy(dtype=float))]
        if len(pure) >= MIN_N:
            evaluate_block(cohort, stratum, pure, with_purity, bucket)
    print(f"[sweep] {cohort} n={len(frame)} rows={len(bucket)}", flush=True)
    return bucket


COLUMNS = [
    "family", "cohort", "stratum", "covariates", "method", "outcome",
    "cut", "predictor_cut", "effect", "p", "var", "k", "n", "n_high", "n_low",
]


def _meta_row(key, sub, family):
    rec = dict(zip(
        ["stratum", "covariates", "method", "outcome"]
        + (["cut", "predictor_cut"] if family == "quantile" else []),
        key if isinstance(key, tuple) else (key,),
    ))
    if sub["cohort"].duplicated().any():
        raise RuntimeError(f"duplicate cohort in {rec}")
    k = int(sub["k"].iloc[0])
    rec["k"] = k
    if family == "continuous":
        meta = fisher_random_effects(sub["effect"].to_numpy(), sub["n"].to_numpy(), k)
    else:
        meta = dl_random_effects(sub["effect"].to_numpy(), sub["var"].to_numpy())
    finite = sub["effect"].to_numpy(dtype=float)
    rec["effect"] = meta["effect"]
    rec["fe"] = meta["fe"]
    rec["se"] = meta["se"]
    rec["p"] = meta["p"]
    rec["I2"] = meta["I2"]
    rec["ci_low"] = meta["ci_low"]
    rec["ci_high"] = meta["ci_high"]
    rec["tau2"] = meta["tau2"]
    rec["neg"] = int(np.sum(finite < 0))
    rec["npos"] = int(np.sum(finite > 0))
    rec["cohorts"] = int(meta["n_cohorts"])
    rec["eligible"] = bool(
        np.isfinite(meta["effect"])
        and meta["effect"] < 0
        and rec["neg"] >= MIN_NEG
        and rec["cohorts"] >= MIN_NEG
        and k >= 1
    )
    return rec


def summarize(detail, family):
    frame = detail[detail["family"] == family]
    keys = ["stratum", "covariates", "method", "outcome"]
    if family == "quantile":
        keys += ["cut", "predictor_cut"]
    rows = [_meta_row(key, sub, family) for key, sub in frame.groupby(keys, sort=False)]
    panel = pd.DataFrame(rows)
    panel["q"] = bh_fdr(panel["p"].to_numpy())
    return panel


def pick(panel, mask):
    sub = panel.loc[mask]
    if sub.empty:
        return None
    local = select_abs(
        sub["effect"].to_numpy(),
        sub["neg"].to_numpy(),
        sub["cohorts"].to_numpy(),
        sub["p"].to_numpy(),
    )
    if local < 0:
        return None
    return sub.iloc[int(local)]


def pick_positive(panel, mask):
    sub = panel.loc[mask]
    if sub.empty:
        return None
    local = select_positive(sub["effect"].to_numpy(), sub["cohorts"].to_numpy(), sub["p"].to_numpy())
    if local < 0:
        return None
    return sub.iloc[int(local)]


def find_spec(panel, **kw):
    hit = panel
    for key, val in kw.items():
        hit = hit[hit[key] == val]
    if len(hit) != 1:
        raise RuntimeError(f"spec {kw} matched {len(hit)} rows")
    return hit.iloc[0]


def fmt(x, digits=3):
    if x is None or not np.isfinite(x):
        return "NA"
    if abs(x) >= 0.001 or x == 0:
        return f"{x:.{digits}f}"
    return f"{x:.2e}"


def fmt_p(x):
    if x is None or not np.isfinite(x):
        return "NA"
    if abs(x) < 1e-4:
        return f"{x:.2e}"
    return f"{x:.4f}"


def md_table(records, columns):
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for rec in records:
        cells = []
        for col in columns:
            val = rec[col]
            if isinstance(val, (float, np.floating)):
                cells.append(fmt_p(val) if col in ("p", "q", "worst_p", "q_worst") else fmt(val))
            else:
                cells.append(str(val))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def spec_label(spec):
    if spec is None:
        return "none"
    base = f"{spec['outcome']} | {spec['covariates']} | {spec['stratum']} | {spec['method']}"
    if "cut" in spec and spec["cut"] != "" and pd.notna(spec["cut"]):
        base += f" | cut {spec['cut']} | {spec['predictor_cut']}"
    return base


def dist_of(panel, mask):
    values = np.abs(panel.loc[mask, "effect"].to_numpy(dtype=float))
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {"n": 0, "median": np.nan, "p95": np.nan, "max": np.nan}
    return {
        "n": int(values.size),
        "median": float(np.median(values)),
        "p95": float(np.quantile(values, 0.95)),
        "max": float(np.max(values)),
    }


def check_anchor(detail):
    rows = []
    diffs = []
    for outcome, expected in ANCHOR_RHO.items():
        for cohort, rho_hat in expected.items():
            hit = detail[
                (detail["family"] == "continuous")
                & (detail["cohort"] == cohort)
                & (detail["stratum"] == "all")
                & (detail["covariates"] == "KRT5_6_14")
                & (detail["method"] == "spearman")
                & (detail["outcome"] == outcome)
            ]
            if hit.empty:
                raise RuntimeError(f"anchor row missing {outcome} {cohort}")
            got = float(hit["effect"].iloc[0])
            diff = abs(got - rho_hat)
            diffs.append(diff)
            rows.append({
                "outcome": outcome,
                "cohort": cohort,
                "expected_rho": rho_hat,
                "got_rho": got,
                "abs_diff": diff,
                "n": int(hit["n"].iloc[0]),
            })
    return pd.DataFrame(rows), float(np.max(diffs))


def cohort_effects(detail, spec, family):
    hit = detail[
        (detail["family"] == family)
        & (detail["stratum"] == spec["stratum"])
        & (detail["covariates"] == spec["covariates"])
        & (detail["method"] == spec["method"])
        & (detail["outcome"] == spec["outcome"])
    ]
    if family == "quantile":
        hit = hit[(hit["cut"] == spec["cut"]) & (hit["predictor_cut"] == spec["predictor_cut"])]
    rows = []
    k = int(spec["k"])
    indexed = {r.cohort: r for r in hit.itertuples(index=False)}
    for cohort in COHORTS:
        if cohort not in indexed:
            rows.append({
                "cohort": cohort, "n": np.nan, "n_low": np.nan, "n_high": np.nan,
                "effect": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p": np.nan,
            })
            continue
        r = indexed[cohort]
        effect = float(r.effect)
        if family == "continuous":
            lo, hi = fisher_ci(effect, int(r.n), k)
            rows.append({
                "cohort": cohort, "n": int(r.n), "n_low": "", "n_high": "",
                "effect": effect, "ci_low": lo, "ci_high": hi, "p": float(r.p),
            })
        else:
            se = np.sqrt(float(r.var)) if np.isfinite(r.var) and r.var > 0 else np.nan
            zcrit = stats.norm.ppf(0.975)
            rows.append({
                "cohort": cohort,
                "n": int(r.n),
                "n_low": int(r.n_low),
                "n_high": int(r.n_high),
                "effect": effect,
                "ci_low": float(effect - zcrit * se) if np.isfinite(se) else np.nan,
                "ci_high": float(effect + zcrit * se) if np.isfinite(se) else np.nan,
                "p": float(r.p),
            })
    return rows


def sign_note(rows):
    nonneg = [r["cohort"] for r in rows if np.isfinite(r["effect"]) and r["effect"] >= 0]
    missing = [r["cohort"] for r in rows if not np.isfinite(r["effect"])]
    parts = []
    if nonneg:
        parts.append("Non-negative cohorts: " + ", ".join(nonneg) + ".")
    else:
        parts.append("Every tested cohort is negative.")
    if missing:
        parts.append("Untested cohorts, not counted as negative: " + ", ".join(missing) + ".")
    return " ".join(parts)


def effect_record(spec):
    return {
        "pooled": spec["effect"],
        "ci_low": spec["ci_low"],
        "ci_high": spec["ci_high"],
        "p": spec["p"],
        "q": spec["q"],
        "I2": spec["I2"],
        "negative_cohorts": int(spec["neg"]),
        "cohorts_in_pool": int(spec["cohorts"]),
        "fixed_effect": spec["fe"],
        "k": int(spec["k"]),
    }


def top_rows(panel, n=12):
    elig = panel.loc[panel["eligible"]].copy()
    elig["abs_effect"] = np.abs(elig["effect"])
    elig = elig.sort_values(["abs_effect", "neg", "cohorts", "p"], ascending=[False, False, False, True]).head(n)
    rows = []
    for rec in elig.to_dict(orient="records"):
        row = {
            "outcome": rec["outcome"],
            "covariates": rec["covariates"],
            "stratum": rec["stratum"],
            "method": rec["method"],
            "abs_effect": rec["abs_effect"],
            "effect": rec["effect"],
            "neg": int(rec["neg"]),
            "cohorts": int(rec["cohorts"]),
            "p": rec["p"],
            "q": rec["q"],
        }
        if "cut" in rec:
            row["cut"] = rec["cut"]
            row["predictor_cut"] = rec["predictor_cut"]
        rows.append(row)
    return rows


def panel_rows(panel, family):
    rows = []
    covs = PANEL_COVARS if family == "continuous" else [c for c in PANEL_COVARS if c in ("KRT5_6_14", "KRT8_18_19", "both_programs", "keratin_z")]
    for outcome in PANEL_OUTCOMES:
        for cov in covs:
            kw = {
                "outcome": outcome,
                "covariates": cov,
                "stratum": "all",
                "method": "spearman",
            }
            if family == "quantile":
                kw["cut"] = "0.25"
                kw["predictor_cut"] = "residual"
            spec = find_spec(panel, **kw)
            rows.append({
                "outcome": outcome,
                "covariates": cov,
                "pooled": spec["effect"],
                "ci_low": spec["ci_low"],
                "ci_high": spec["ci_high"],
                "p": spec["p"],
                "I2": spec["I2"],
                "neg": int(spec["neg"]),
                "cohorts": int(spec["cohorts"]),
            })
    return rows


def save_fig(fig, name):
    os.makedirs(FIGS, exist_ok=True)
    os.makedirs(ART, exist_ok=True)
    path = os.path.join(FIGS, name)
    fig.savefig(path, dpi=150)
    fig.savefig(os.path.join(ART, name), dpi=150)
    plt.close(fig)
    return path


def forest_plot(cohort_rows, spec, family, path_name, title):
    labels = []
    effects = []
    lo = []
    hi = []
    for row in cohort_rows:
        arm = ""
        if family == "quantile" and row["n_low"] != "" and np.isfinite(row["n"]):
            arm = f" ({int(row['n_low'])}/{int(row['n_high'])})"
        if np.isfinite(row["effect"]):
            n_label = int(row["n"]) if np.isfinite(row["n"]) else "NA"
            labels.append(f"{row['cohort']} n={n_label}{arm}")
            effects.append(row["effect"])
            lo.append(row["ci_low"])
            hi.append(row["ci_high"])
        else:
            labels.append(f"{row['cohort']} untested")
            effects.append(np.nan)
            lo.append(np.nan)
            hi.append(np.nan)
    labels.append("pooled")
    effects.append(float(spec["effect"]))
    lo.append(float(spec["ci_low"]))
    hi.append(float(spec["ci_high"]))
    y = np.arange(len(labels))[::-1]
    effects = np.asarray(effects, dtype=float)
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    ax.axvline(0, color="#888888", lw=0.8)
    finite = np.isfinite(effects) & np.isfinite(lo) & np.isfinite(hi)
    ax.errorbar(
        effects[finite], y[finite],
        xerr=[effects[finite] - lo[finite], hi[finite] - effects[finite]],
        fmt="o", color="#1f4e79", ecolor="#1f4e79", capsize=2, markersize=5,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    xlabel = "partial correlation" if family == "continuous" else "Cliff's delta"
    ax.set_xlabel(xlabel)
    wrapped = "\n".join(__import__("textwrap").wrap(title, 78))
    ax.set_title(wrapped, fontsize=9)
    fig.tight_layout()
    save_fig(fig, path_name)


def hist_plot(panel, winner, path_name, title, xlabel):
    values = np.abs(panel.loc[panel["eligible"], "effect"].to_numpy(dtype=float))
    values = values[np.isfinite(values)]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    if values.size:
        ax.hist(values, bins=40, color="#5d6d7e", edgecolor="white")
    if winner is not None and np.isfinite(winner["effect"]):
        ax.axvline(abs(float(winner["effect"])), color="#922b21", lw=1.6, label=f"maximum {abs(float(winner['effect'])):.3f}")
        ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("eligible specifications")
    ax.set_title(title)
    fig.tight_layout()
    save_fig(fig, path_name)


def rec_to_json(spec):
    if spec is None:
        return None
    out = {}
    for key, val in spec.to_dict().items():
        if isinstance(val, (np.floating, float)):
            out[key] = None if not np.isfinite(val) else float(val)
        elif isinstance(val, (np.integer, int)):
            out[key] = int(val)
        elif isinstance(val, (np.bool_, bool)):
            out[key] = bool(val)
        else:
            out[key] = val
    out["label"] = spec_label(spec)
    return out


def write_spec_block(lines, spec, cohort_rows, family):
    lines.append(f"`{spec_label(spec)}`.")
    lines.append("")
    lines.append(md_table([effect_record(spec)], [
        "pooled", "ci_low", "ci_high", "p", "q", "I2", "negative_cohorts", "cohorts_in_pool", "fixed_effect", "k",
    ]))
    lines.append("")
    lines.append(sign_note(cohort_rows))
    lines.append("")
    cols = ["cohort", "n", "effect", "ci_low", "ci_high", "p"]
    if family == "quantile":
        cols = ["cohort", "n", "n_low", "n_high", "effect", "ci_low", "ci_high", "p"]
    lines.append(md_table(cohort_rows, cols))
    lines.append("")


def write_results(path, context):
    c = context
    lines = []
    lines.append("# TCGA TACSTD2 versus CD8/CD3 after keratin")
    lines.append("")
    lines.append("Primary tumors only (sample type 01), one row per patient, Xena GDC STAR log2(TPM+1). Eight cohorts: LUAD, LUSC, BRCA, CESC, KIRC, STAD, BLCA, PAAD. This file is written by `scripts/tcga_tacstd2_keratin_cd8/sweep.py` from the tables next to it.")
    lines.append("")
    lines.append("## Objective")
    lines.append("")
    lines.append("The grid is declared in `sweep.py` before the matrices are scored. Immune scores are CD3D, CD3E, CD3G, CD8A, CD8B, and the mean, z-mean, rank-mean, and first principal component of CD3, of CD8, and of CD3+CD8. Covariates are basal keratins (KRT5/6/14, the broader basal set, and that set plus TP63), simple keratins (KRT8/18/19, plus KRT7 and KRT20), z-means and principal components of those sets, the basal-plus-simple set, and the same keratin scores with ESTIMATE purity. Strata are all patients, basal and simple keratin median and tertile splits, and an ESTIMATE median split.")
    lines.append("")
    lines.append("A keratin-adjusted specification has at least one covariate. It is eligible when the random-effects pooled association is negative, at least 6 of the 8 cohorts are negative, and at least 6 cohorts enter the pool. Untested cohorts are not counted as negative. The continuous objective maximizes |pooled ρ|. The Cliff objective maximizes |pooled δ|. Negative δ means the upper TACSTD2 group has the lower immune score. Correlations use a Fisher-z DerSimonian-Laird pool. Deltas use a DerSimonian-Laird pool on the delta scale with the Cliff 1993 variance.")
    lines.append("")
    lines.append("The maximum is the extreme of this search. Selecting it pushes |effect| upward relative to a single pre-specified test. `q` is Benjamini-Hochberg on the meta p across every specification in that family, eligible or not.")
    lines.append("")
    lines.append("`stratum_only` drops the covariate inside a keratin split. Those rows are stored and are excluded from the keratin-adjusted maximum, because that maximum is a partial association.")
    lines.append("")
    lines.append("## Data")
    lines.append("")
    lines.append(md_table(c["counts"], ["cohort", "n_patients", "n_with_ESTIMATE"]))
    lines.append("")
    lines.append(
        f"Continuous specifications: {c['n_cont']}. Quantile specifications: {c['n_quant']}. "
        f"Keratin-adjusted eligible continuous: {c['dist_cont']['n']}. "
        f"Keratin-adjusted eligible quantile: {c['dist_quant']['n']}."
    )
    lines.append("")
    lines.append("Partial Spearman of CD8_mean and CD3_mean on KRT5+KRT6A+KRT6B+KRT14, all patients, was recomputed and compared with the previous KRT5/6 grid at 6 decimal places. Maximum absolute difference across 16 cohort-level TACSTD2 correlations: " + f"{c['anchor_diff']:.3e}.")
    lines.append("")
    lines.append("## Pre-specified all-patient Spearman panel")
    lines.append("")
    lines.append("These rows are fixed in the script. They are the keratin-adjusted CD8/CD3 associations on all patients, before any search for a larger effect.")
    lines.append("")
    lines.append(md_table(c["cont_panel"], ["outcome", "covariates", "pooled", "ci_low", "ci_high", "p", "I2", "neg", "cohorts"]))
    lines.append("")
    lines.append("## Continuous maximum")
    lines.append("")
    d_all_s = c["dist_all_spearman"]
    lines.append(
        f"Eligible all-patient Spearman specifications: {d_all_s['n']}. "
        f"Median |ρ| {fmt(d_all_s['median'])}, 95th percentile {fmt(d_all_s['p95'])}, maximum {fmt(d_all_s['max'])}."
    )
    lines.append("")
    lines.append("All-patient Spearman maximum:")
    lines.append("")
    write_spec_block(lines, c["rho_all_s"], c["rho_all_s_cohort"], "continuous")
    if spec_label(c["rho_all"]) != spec_label(c["rho_all_s"]):
        d_all = c["dist_all"]
        lines.append(
            f"Eligible all-patient specifications across Spearman, Pearson, and winsorized Pearson: {d_all['n']}. "
            f"Median |ρ| {fmt(d_all['median'])}, 95th percentile {fmt(d_all['p95'])}, maximum {fmt(d_all['max'])}."
        )
        lines.append("")
        lines.append("All-patient maximum:")
        lines.append("")
        write_spec_block(lines, c["rho_all"], c["rho_all_cohort"], "continuous")
    d = c["dist_cont"]
    lines.append(
        f"Eligible keratin-adjusted specifications in the full continuous grid: {d['n']}. "
        f"Median |ρ| {fmt(d['median'])}, 95th percentile {fmt(d['p95'])}, maximum {fmt(d['max'])}."
    )
    lines.append("")
    lines.append("Full-grid maximum:")
    lines.append("")
    write_spec_block(lines, c["rho_max"], c["rho_max_cohort"], "continuous")
    if c["rho_max8"] is not None and spec_label(c["rho_max8"]) != spec_label(c["rho_max"]):
        lines.append("Largest |ρ| among eligible specifications that keep all 8 cohorts in the pool:")
        lines.append("")
        write_spec_block(lines, c["rho_max8"], c["rho_max8_cohort"], "continuous")
    if c["rho_stable"] is not None and spec_label(c["rho_stable"]) != spec_label(c["rho_max"]):
        lines.append("Largest |ρ| among eligible specifications whose pooled 95% interval lies entirely below 0:")
        lines.append("")
        write_spec_block(lines, c["rho_stable"], c["rho_stable_cohort"], "continuous")
    lines.append("Largest eligible continuous specifications:")
    lines.append("")
    lines.append(md_table(c["cont_top"], ["outcome", "covariates", "stratum", "method", "abs_effect", "effect", "neg", "cohorts", "p", "q"]))
    lines.append("")
    opp = c["rho_pos"]
    if opp is not None:
        lines.append(
            f"Largest positive pooled ρ with at least 6 cohorts in the pool: {fmt(opp['effect'])} "
            f"at `{spec_label(opp)}` (p={fmt_p(opp['p'])}, {int(opp['npos'])} cohorts positive, {int(opp['cohorts'])} cohorts in the pool)."
        )
        lines.append("")
    lines.append("## Cliff's delta")
    lines.append("")
    lines.append("Quartiles are the outer 25%. The immune score is the keratin residual. `residual` cuts TACSTD2 on that residual; `expression` cuts TACSTD2 on log2(TPM+1). Arms below 12 patients are untested.")
    lines.append("")
    lines.append("Pre-specified all-patient quartile contrasts, Spearman residuals, cut 0.25:")
    lines.append("")
    lines.append(md_table(c["quant_panel"], ["outcome", "covariates", "pooled", "ci_low", "ci_high", "p", "I2", "neg", "cohorts"]))
    lines.append("")
    dq = c["dist_quartile"]
    lines.append(
        f"Eligible all-patient Spearman residual quartile contrasts: {dq['n']}. "
        f"Median |δ| {fmt(dq['median'])}, 95th percentile {fmt(dq['p95'])}, maximum {fmt(dq['max'])}."
    )
    lines.append("")
    lines.append("All-patient Spearman residual quartile maximum:")
    lines.append("")
    write_spec_block(lines, c["delta_q"], c["delta_q_cohort"], "quantile")
    dd = c["dist_quant"]
    lines.append(
        f"Eligible keratin-adjusted quantile specifications: {dd['n']}. "
        f"Median |δ| {fmt(dd['median'])}, 95th percentile {fmt(dd['p95'])}, maximum {fmt(dd['max'])}."
    )
    lines.append("")
    lines.append("Full-grid Cliff maximum:")
    lines.append("")
    write_spec_block(lines, c["delta_max"], c["delta_max_cohort"], "quantile")
    if c["delta_max8"] is not None and spec_label(c["delta_max8"]) != spec_label(c["delta_max"]):
        lines.append("Largest |δ| among eligible specifications that keep all 8 cohorts in the pool:")
        lines.append("")
        write_spec_block(lines, c["delta_max8"], c["delta_max8_cohort"], "quantile")
    lines.append("Largest eligible Cliff specifications:")
    lines.append("")
    lines.append(md_table(c["quant_top"], [
        "outcome", "covariates", "stratum", "method", "cut", "predictor_cut",
        "abs_effect", "effect", "neg", "cohorts", "p", "q",
    ]))
    lines.append("")
    lines.append("## Reading")
    lines.append("")
    lines.append(c["reading"])
    lines.append("")
    lines.append("These are bulk RNA associations in TCGA primary tumors after the named keratin covariates. The maximum is the extreme of the pre-declared grid.")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 scripts/tcga_tacstd2_keratin_cd8/test_stats.py")
    lines.append("python3 scripts/tcga_tacstd2_keratin_cd8/sweep.py")
    lines.append("```")
    lines.append("")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


def reading_paragraph(context):
    a = context["anchor_cd8"]
    s = context["simple_cd8"]
    r = context["rho_max"]
    d = context["delta_max"]
    rs = context["rho_all_s"]
    dq = context["delta_q"]
    parts = [
        f"On all patients, partial Spearman of CD8_mean on KRT5+KRT6A+KRT6B+KRT14 is pooled ρ {fmt(a['effect'])} ({int(a['neg'])}/8 negative, p={fmt_p(a['p'])}).",
        f"The same CD8_mean score on KRT8+KRT18+KRT19 is pooled ρ {fmt(s['effect'])} ({int(s['neg'])}/8 negative, p={fmt_p(s['p'])}).",
        f"The largest eligible all-patient Spearman |ρ| is {fmt(abs(rs['effect']))} at `{spec_label(rs)}` (pooled ρ {fmt(rs['effect'])}, 95% CI {fmt(rs['ci_low'])} to {fmt(rs['ci_high'])}, p={fmt_p(rs['p'])}, q={fmt_p(rs['q'])}, {int(rs['neg'])}/8 negative).",
        f"The keratin-adjusted continuous maximum is |ρ| {fmt(abs(r['effect']))} at `{spec_label(r)}` (pooled ρ {fmt(r['effect'])}, 95% CI {fmt(r['ci_low'])} to {fmt(r['ci_high'])}, p={fmt_p(r['p'])}, q={fmt_p(r['q'])}, {int(r['neg'])}/8 negative, I²={fmt(r['I2'])}).",
        (
            f"The same outcome, covariates, and stratum without winsorization are Spearman pooled ρ {fmt(context['rho_twin_spearman']['effect'])} "
            f"and Pearson pooled ρ {fmt(context['rho_twin_pearson']['effect'])}, both with {int(context['rho_twin_spearman']['neg'])}/8 and {int(context['rho_twin_pearson']['neg'])}/8 cohorts negative."
        ),
        f"The all-patient Spearman residual quartile maximum is |δ| {fmt(abs(dq['effect']))} at `{spec_label(dq)}` (pooled δ {fmt(dq['effect'])}, {int(dq['neg'])}/8 negative, p={fmt_p(dq['p'])}).",
        f"The keratin-adjusted Cliff maximum is |δ| {fmt(abs(d['effect']))} at `{spec_label(d)}` (pooled δ {fmt(d['effect'])}, 95% CI {fmt(d['ci_low'])} to {fmt(d['ci_high'])}, p={fmt_p(d['p'])}, q={fmt_p(d['q'])}, {int(d['neg'])}/8 negative).",
    ]
    return " ".join(parts)


def row_from_spec(role, spec, family):
    row = {
        "role": role,
        "family": family,
        "label": spec_label(spec),
        "outcome": spec["outcome"],
        "covariates": spec["covariates"],
        "stratum": spec["stratum"],
        "method": spec["method"],
        "cut": spec["cut"] if "cut" in spec.index else "",
        "predictor_cut": spec["predictor_cut"] if "predictor_cut" in spec.index else "",
        "effect": spec["effect"],
        "ci_low": spec["ci_low"],
        "ci_high": spec["ci_high"],
        "p": spec["p"],
        "q": spec["q"],
        "I2": spec["I2"],
        "neg": int(spec["neg"]),
        "cohorts": int(spec["cohorts"]),
        "k": int(spec["k"]),
        "eligible": bool(spec["eligible"]),
    }
    return row


def main():
    os.makedirs(TABLES, exist_ok=True)
    os.makedirs(FIGS, exist_ok=True)
    os.makedirs(ART, exist_ok=True)
    symbols = symbols_needed()
    mapping = load_probemap(symbols)
    id_to_gene = {ens: gene for gene, ens in mapping.items()}
    print(f"[map] {len(mapping)} genes", flush=True)
    metas = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(extract_cohort, cohort, id_to_gene): cohort for cohort in COHORTS}
        for fut in as_completed(futures):
            meta = fut.result()
            metas.append(meta)
            print(f"[cohort] {meta['cohort']} n={meta['n_patients']} cache={meta['from_cache']}", flush=True)
    purity, purity_meta = load_purity()
    frames = {}
    counts = []
    for cohort in COHORTS:
        frame = pd.read_csv(os.path.join(CACHE, f"{cohort}.primary01.tsv.gz"), sep="\t")
        sub_p = purity.loc[purity["cohort"] == cohort, ["patient", "ESTIMATE"]]
        frame = frame.merge(sub_p, on="patient", how="left")
        frames[cohort] = frame
        counts.append({
            "cohort": cohort,
            "n_patients": int(frame.shape[0]),
            "n_with_ESTIMATE": int(np.isfinite(frame["ESTIMATE"]).sum()),
        })
        # Spot-check the scale: log2(TPM+1) sits well below raw TPM.
        tac = frame["TACSTD2"].to_numpy(dtype=float)
        if np.nanmax(tac) > 40:
            raise SystemExit(f"{cohort} TACSTD2 max {np.nanmax(tac)} is outside log2(TPM+1) range")
    global DATA
    DATA = frames
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(max_workers=4, mp_context=ctx) as pool:
        parts = list(pool.map(evaluate_cohort, COHORTS))
    detail = pd.DataFrame([row for part in parts for row in part], columns=COLUMNS)
    print(f"[detail] rows={len(detail)}", flush=True)
    anchor, anchor_diff = check_anchor(detail)
    anchor.to_csv(os.path.join(TABLES, "anchor_reproduction.tsv"), sep="\t", index=False)
    print(f"[anchor] max abs diff {anchor_diff:.3e}", flush=True)
    if anchor_diff > 1e-4:
        raise SystemExit(f"anchor reproduction failed: max abs diff {anchor_diff}")
    cont = summarize(detail, "continuous")
    quant = summarize(detail, "quantile")
    cont.to_csv(os.path.join(TABLES, "continuous_panel.tsv"), sep="\t", index=False)
    quant.to_csv(os.path.join(TABLES, "quantile_panel.tsv"), sep="\t", index=False)
    detail.to_csv(os.path.join(TABLES, "cohort_effects.tsv.gz"), sep="\t", index=False, compression="gzip")

    rho_all_s = pick(cont, (cont["eligible"]) & (cont["stratum"] == "all") & (cont["method"] == "spearman"))
    rho_all = pick(cont, (cont["eligible"]) & (cont["stratum"] == "all"))
    rho_max = pick(cont, cont["eligible"])
    rho_max8 = pick(cont, (cont["eligible"]) & (cont["cohorts"] == 8))
    rho_stable = pick(cont, (cont["eligible"]) & (cont["ci_high"] < 0))
    rho_pos = pick_positive(cont, cont["k"] >= 1)
    delta_q = pick(
        quant,
        (quant["eligible"])
        & (quant["stratum"] == "all")
        & (quant["method"] == "spearman")
        & (quant["cut"] == "0.25")
        & (quant["predictor_cut"] == "residual"),
    )
    delta_max = pick(quant, quant["eligible"])
    delta_max8 = pick(quant, (quant["eligible"]) & (quant["cohorts"] == 8))
    for name, spec in [
        ("rho_all_s", rho_all_s), ("rho_all", rho_all), ("rho_max", rho_max),
        ("delta_q", delta_q), ("delta_max", delta_max),
    ]:
        if spec is None:
            raise SystemExit(f"no eligible specification for {name}")

    winner_specs = [
        ("rho_all_spearman", rho_all_s, "continuous"),
        ("rho_all", rho_all, "continuous"),
        ("rho_grid_max", rho_max, "continuous"),
        ("rho_grid_max_8cohorts", rho_max8, "continuous"),
        ("rho_ci_entirely_negative", rho_stable, "continuous"),
        ("delta_quartile_all_spearman_residual", delta_q, "quantile"),
        ("delta_grid_max", delta_max, "quantile"),
        ("delta_grid_max_8cohorts", delta_max8, "quantile"),
        ("rho_largest_positive", rho_pos, "continuous"),
    ]
    winners = [row_from_spec(role, spec, family) for role, spec, family in winner_specs if spec is not None]
    pd.DataFrame(winners).to_csv(os.path.join(TABLES, "winners.tsv"), sep="\t", index=False)

    context = {
        "counts": counts,
        "n_cont": int((detail["family"] == "continuous").sum() and cont.shape[0]),
        "n_quant": int(quant.shape[0]),
        "anchor_diff": anchor_diff,
        "dist_cont": dist_of(cont, cont["eligible"]),
        "dist_all": dist_of(cont, (cont["eligible"]) & (cont["stratum"] == "all")),
        "dist_all_spearman": dist_of(cont, (cont["eligible"]) & (cont["stratum"] == "all") & (cont["method"] == "spearman")),
        "dist_quant": dist_of(quant, quant["eligible"]),
        "dist_quartile": dist_of(
            quant,
            (quant["eligible"])
            & (quant["stratum"] == "all")
            & (quant["method"] == "spearman")
            & (quant["cut"] == "0.25")
            & (quant["predictor_cut"] == "residual"),
        ),
        "cont_panel": panel_rows(cont, "continuous"),
        "quant_panel": panel_rows(quant, "quantile"),
        "cont_top": top_rows(cont),
        "quant_top": top_rows(quant),
        "rho_all_s": rho_all_s,
        "rho_all": rho_all,
        "rho_max": rho_max,
        "rho_max8": rho_max8,
        "rho_stable": rho_stable,
        "rho_pos": rho_pos,
        "delta_q": delta_q,
        "delta_max": delta_max,
        "delta_max8": delta_max8,
        "anchor_cd8": find_spec(cont, outcome="CD8_mean", covariates="KRT5_6_14", stratum="all", method="spearman"),
        "simple_cd8": find_spec(cont, outcome="CD8_mean", covariates="KRT8_18_19", stratum="all", method="spearman"),
    }
    context["n_cont"] = int(cont.shape[0])
    context["rho_all_s_cohort"] = cohort_effects(detail, rho_all_s, "continuous")
    context["rho_all_cohort"] = cohort_effects(detail, rho_all, "continuous")
    context["rho_max_cohort"] = cohort_effects(detail, rho_max, "continuous")
    context["rho_max8_cohort"] = cohort_effects(detail, rho_max8, "continuous")
    context["rho_stable_cohort"] = cohort_effects(detail, rho_stable, "continuous") if rho_stable is not None else []
    context["delta_q_cohort"] = cohort_effects(detail, delta_q, "quantile")
    context["delta_max_cohort"] = cohort_effects(detail, delta_max, "quantile")
    context["delta_max8_cohort"] = cohort_effects(detail, delta_max8, "quantile")
    context["rho_twin_spearman"] = find_spec(
        cont,
        outcome=rho_max["outcome"],
        covariates=rho_max["covariates"],
        stratum=rho_max["stratum"],
        method="spearman",
    )
    context["rho_twin_pearson"] = find_spec(
        cont,
        outcome=rho_max["outcome"],
        covariates=rho_max["covariates"],
        stratum=rho_max["stratum"],
        method="pearson",
    )
    context["reading"] = reading_paragraph(context)

    forest_plot(context["rho_all_s_cohort"], rho_all_s, "continuous", "rho_all_spearman_forest.png", spec_label(rho_all_s))
    forest_plot(context["rho_max_cohort"], rho_max, "continuous", "rho_grid_max_forest.png", spec_label(rho_max))
    forest_plot(context["delta_q_cohort"], delta_q, "quantile", "delta_quartile_forest.png", spec_label(delta_q))
    forest_plot(context["delta_max_cohort"], delta_max, "quantile", "delta_grid_max_forest.png", spec_label(delta_max))
    hist_plot(cont, rho_max, "rho_eligible_hist.png", "Eligible keratin-adjusted partial correlations", "|pooled ρ| of TACSTD2")
    hist_plot(quant, delta_max, "delta_eligible_hist.png", "Eligible keratin-adjusted Cliff contrasts", "|pooled Cliff's δ| of TACSTD2")

    write_results(os.path.join(OUT, "RESULTS.md"), context)
    summary = {
        "anchor_max_abs_diff": anchor_diff,
        "n_continuous": int(cont.shape[0]),
        "n_quantile": int(quant.shape[0]),
        "n_eligible_continuous": context["dist_cont"]["n"],
        "n_eligible_quantile": context["dist_quant"]["n"],
        "dist_continuous": context["dist_cont"],
        "dist_all_spearman": context["dist_all_spearman"],
        "dist_quantile": context["dist_quant"],
        "rho_all_spearman": rec_to_json(rho_all_s),
        "rho_grid_max": rec_to_json(rho_max),
        "rho_grid_max_8": rec_to_json(rho_max8),
        "rho_stable": rec_to_json(rho_stable),
        "rho_positive": rec_to_json(rho_pos),
        "delta_quartile": rec_to_json(delta_q),
        "delta_grid_max": rec_to_json(delta_max),
        "delta_grid_max_8": rec_to_json(delta_max8),
        "anchor_cd8_mean": rec_to_json(context["anchor_cd8"]),
        "simple_cd8_mean": rec_to_json(context["simple_cd8"]),
        "rho_grid_max_cohorts": context["rho_max_cohort"],
        "delta_grid_max_cohorts": context["delta_max_cohort"],
        "reading": context["reading"],
    }
    provenance = {
        "source": "UCSC Xena GDC hub STAR log2(TPM+1), primary tumor sample type 01",
        "probemap": PROBEMAP_URL,
        "probemap_sha256": sha256_file(os.path.join(CACHE, "gencode.v36.annotation.gtf.gene.probemap")),
        "genes": mapping,
        "cohorts": metas,
        "estimate": purity_meta,
        "estimate_transform": "clip(cos(0.6049872018 + 0.0001467884 * ESTIMATE_score), 0, 1)",
        "objective": "maximize |pooled rho| and |Cliff delta| for TACSTD2 vs CD8/CD3 after keratin, pooled effect negative, >=6/8 cohorts negative",
        "counts": counts,
    }
    with open(os.path.join(TABLES, "summary.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    with open(os.path.join(TABLES, "provenance.json"), "w", encoding="utf-8") as handle:
        json.dump(provenance, handle, indent=2)
    pd.DataFrame(counts).to_csv(os.path.join(TABLES, "cohort_counts.tsv"), sep="\t", index=False)
    print("[done]", context["reading"], flush=True)


if __name__ == "__main__":
    main()
