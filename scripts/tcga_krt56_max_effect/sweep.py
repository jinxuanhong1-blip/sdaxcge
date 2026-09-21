#!/usr/bin/env python3
"""Maximize the TCGA CD3/CD8 association of TACSTD2, CLDN4, and CLDN7 after KRT5/6.

The objective is fixed here, before the grid is evaluated.

Question
--------
Across the same eight primary-tumor cohorts used previously (LUAD, LUSC,
BRCA, CESC, KIRC, STAD, BLCA, PAAD), which pre-declared CD3 or CD8 score,
KRT5/6 covariate set, keratin or purity stratum, and quartile cut makes
TACSTD2, CLDN4, and CLDN7 move together toward a lower T-cell score?

Eligibility, both for partial correlations and for Cliff's delta
-----------------------------------------------------------------
- The covariate set contains KRT5 and a KRT6 gene, or it is the first
  principal component of KRT5+KRT6A+KRT6B+KRT14, or (only inside a KRT5/6
  stratum) there is no further covariate.
- Each of the three predictors is negative in at least 7 of the 8 cohorts.
- The random-effects pooled effect for each predictor is finite and negative.
- At least 7 cohorts enter the pool.

Objective
---------
Maximize the minimum of the three absolute pooled effects. Tie-break: larger
mean absolute pooled effect, then smaller worst of the three meta p-values.
Correlations are pooled by DerSimonian-Laird on the Fisher z scale.
Cliff's deltas are pooled by DerSimonian-Laird on the delta scale, using the
Cliff 1993 sampling variance.

This maximum is a specification maximum. It is not a confirmatory estimate,
and selecting it inflates |effect| relative to any single pre-specified test.
Benjamini-Hochberg q-values are computed across every specification in that
family, including specifications that are not eligible. Those q-values do not
remove the winner's curse on the point estimate.

The locked CLDN4-versus-KRT8 surface-gene screen is not rerun.
Cytotoxic-only scores are not in this grid. The question is CD3/CD8.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

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
    fisher_ci,
    fisher_random_effects,
    matrix_partial_corr,
    residualize,
    select_max,
    dl_random_effects,
)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "results", "tcga_krt56_max_effect")
TABLES = os.path.join(OUT, "tables")
FIGS = os.path.join(OUT, "figures")
CACHE = os.environ.get("TCGA_KRT56_CACHE", "/tmp/tcga_krt56")
GDC = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"
PROBEMAP_URL = f"{GDC}/gencode.v36.annotation.gtf.gene.probemap"
MDACC_BASE = "https://bioinformatics.mdanderson.org/estimate/tables"
ESTIMATE_A = 0.6049872018
ESTIMATE_B = 0.0001467884

COHORTS = ["LUAD", "LUSC", "BRCA", "CESC", "KIRC", "STAD", "BLCA", "PAAD"]
PREDICTORS = ["TACSTD2", "CLDN4", "CLDN7"]
MIN_N = 40
MIN_ARM = 12
CUTS = (0.25, 0.20, 0.10)
METHODS = ("spearman", "pearson", "pearson_winsor")
KRT_STRATA = {"krt_high", "krt_low", "krt_top_tertile", "krt_bot_tertile"}

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

# Published multi-gene CD3/CD8 scores. Single genes are entered separately.
# Danaher et al., J Immunother Cancer 2017: T cells and CD8 T cells.
# Signature-H: Mangino et al. style 15-gene T-cell list, Int J Mol Sci 2019
# (CD2, CD247, CD28, CD3D, CD3G, CD6, GPR171, GZMK, ICOS, ITK, KLRB1, PYHIN1,
# TIGIT, TRAT1, TRBC1). Tsig is the prior CD3-centered T score.
GENE_SETS = {
    "CD3": ["CD3D", "CD3E", "CD3G"],
    "CD8": ["CD8A", "CD8B"],
    "CD3CD8": ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B"],
    "DanaherT": ["CD6", "CD3D", "CD3E", "SH2D1A", "TRAT1", "CD3G"],
    "SigH": [
        "CD2", "CD247", "CD28", "CD3D", "CD3G", "CD6", "GPR171", "GZMK",
        "ICOS", "ITK", "KLRB1", "PYHIN1", "TIGIT", "TRAT1", "TRBC1",
    ],
    "Tsig": ["CD3D", "CD3E", "CD3G", "CD2", "CD247", "LCK"],
}
SINGLES = ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B"]
ANCHOR = ["KRT5", "KRT6A", "KRT6B", "KRT14"]

# (name, genes, needs ESTIMATE). squamous_PC1 and stratum_only are special-cased.
COVARIATE_SETS = [
    ("KRT5_KRT6A", ["KRT5", "KRT6A"], False),
    ("KRT5_KRT6AB", ["KRT5", "KRT6A", "KRT6B"], False),
    ("KRT5_6_14", ANCHOR, False),
    ("basal_broad", ["KRT5", "KRT6A", "KRT6B", "KRT6C", "KRT14", "KRT16", "KRT17"], False),
    ("KRT5_6_14_TP63", ANCHOR + ["TP63"], False),
    ("KRT5_6_14_KRT13", ANCHOR + ["KRT13"], False),
    ("KRT5_6_14_MKI67", ANCHOR + ["MKI67"], False),
    ("KRT5_6_14_EPCAM", ANCHOR + ["EPCAM"], False),
    ("KRT5_6_14_CD68", ANCHOR + ["CD68"], False),
    ("KRT5_6_14_S100A8_A9", ANCHOR + ["S100A8", "S100A9"], False),
    ("KRT5_6_14_stroma", ANCHOR + ["FAP", "COL1A1"], False),
    ("KRT5_6_14_simpleKRT", ANCHOR + ["KRT8", "KRT18", "KRT19"], False),
    ("KRT5_6_14_ESTIMATE", ANCHOR + ["ESTIMATE"], True),
    ("KRT5_KRT6A_ESTIMATE", ["KRT5", "KRT6A", "ESTIMATE"], True),
    ("squamous_PC1", ANCHOR, False),
    ("stratum_only", [], False),
]

# Prior shared-panel partial Spearman, KRT5+KRT6A+KRT6B+KRT14, no purity.
# Rounded to 6 decimals from that grid. Used only as a pipeline check.
ANCHOR_RHO = {
    ("CD8_mean", "TACSTD2"): {
        "BLCA": -0.218597, "BRCA": -0.147181, "CESC": 0.080593, "KIRC": -0.074849,
        "LUAD": -0.142577, "LUSC": -0.141121, "PAAD": -0.260992, "STAD": -0.039386,
    },
    ("CD8_mean", "CLDN4"): {
        "BLCA": -0.134772, "BRCA": -0.012995, "CESC": -0.087846, "KIRC": -0.095033,
        "LUAD": -0.116250, "LUSC": -0.118110, "PAAD": -0.347770, "STAD": -0.085426,
    },
    ("CD8_mean", "CLDN7"): {
        "BLCA": -0.136085, "BRCA": -0.164944, "CESC": -0.000992, "KIRC": -0.027037,
        "LUAD": -0.281196, "LUSC": -0.107626, "PAAD": -0.374256, "STAD": -0.064884,
    },
    ("CD3_mean", "TACSTD2"): {
        "BLCA": -0.159745, "BRCA": -0.149325, "CESC": 0.079720, "KIRC": -0.043927,
        "LUAD": -0.104046, "LUSC": -0.094443, "PAAD": -0.250227, "STAD": -0.106670,
    },
    ("CD3_mean", "CLDN4"): {
        "BLCA": -0.111376, "BRCA": -0.044279, "CESC": -0.131997, "KIRC": -0.138020,
        "LUAD": -0.146944, "LUSC": -0.130279, "PAAD": -0.299999, "STAD": -0.153603,
    },
    ("CD3_mean", "CLDN7"): {
        "BLCA": -0.175195, "BRCA": -0.187964, "CESC": -0.016413, "KIRC": -0.000723,
        "LUAD": -0.341412, "LUSC": -0.141246, "PAAD": -0.348581, "STAD": -0.109433,
    },
}


def symbols_needed():
    genes = set(PREDICTORS)
    genes.update(SINGLES)
    for items in GENE_SETS.values():
        genes.update(items)
    for _, items, _ in COVARIATE_SETS:
        genes.update(items)
    genes.discard("ESTIMATE")
    return sorted(genes)


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url, dest, retries=3):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "tcga-krt56-max/1.0"})
            with urllib.request.urlopen(req, timeout=180) as resp, open(tmp, "wb") as out:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
            os.replace(tmp, dest)
            return
        except Exception as exc:  # noqa: BLE001 - retry network errors
            last = exc
            if os.path.exists(tmp):
                os.remove(tmp)
            print(f"[retry {attempt}] {url} :: {exc}", flush=True)
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
    req = urllib.request.Request(url, headers={"User-Agent": "tcga-krt56-max/1.0"})
    collected = {}
    with urllib.request.urlopen(req, timeout=180) as resp:
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
    return pd.Series(np.clip(purity, 0.0, 1.0), index=score.index)


def load_purity():
    rows = []
    for cohort, fname in MDACC_ESTIMATE.items():
        path = os.path.join(CACHE, "estimate", fname)
        download(f"{MDACC_BASE}/{fname}", path)
        table = pd.read_csv(path, sep="\t")
        table["patient"] = table["ID"].map(lambda s: patient_id(str(s)))
        table = table.dropna(subset=["patient"])
        table["ESTIMATE"] = estimate_purity_from_score(table["ESTIMATE_score"])
        table["cohort"] = cohort
        rows.append(table[["cohort", "patient", "ESTIMATE"]])
    estimate = pd.concat(rows, ignore_index=True)
    return estimate.groupby(["cohort", "patient"], as_index=False)["ESTIMATE"].mean()


def k_of(name, genes):
    if name == "stratum_only":
        return 0
    if name == "squamous_PC1":
        return 1
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


def build_outcomes(columns):
    """Return outcome names and an (n, p) matrix on the raw log2(TPM+1) scale."""
    names = []
    cols = []
    n = len(next(iter(columns.values())))
    for gene in SINGLES:
        names.append(gene)
        cols.append(columns[gene])
    for set_name, genes in GENE_SETS.items():
        mat = np.column_stack([columns[g] for g in genes])
        orient_gene = next((g for g in ("CD3D", "CD3E", "CD3G", "CD8A", "CD8B") if g in genes), genes[0])
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


def scale_matrix(mat, method):
    if method == "spearman":
        out = np.empty(mat.shape, dtype=float)
        for j in range(mat.shape[1]):
            out[:, j] = stats.rankdata(mat[:, j]).astype(float)
        return out
    if method == "pearson_winsor":
        lo = np.quantile(mat, 0.01, axis=0)
        hi = np.quantile(mat, 0.99, axis=0)
        return np.clip(mat, lo, hi)
    if method == "pearson":
        return np.array(mat, dtype=float, copy=True)
    raise ValueError(method)


def covariate_columns(columns, name, genes, method):
    if name == "stratum_only":
        return []
    if name == "squamous_PC1":
        mat = np.column_stack([columns[g] for g in genes])
        score = first_pc(mat, columns["KRT5"])
        scaled = scale_matrix(score[:, None], method)
        return [scaled[:, 0]]
    mat = np.column_stack([columns[g] for g in genes])
    scaled = scale_matrix(mat, method)
    return [scaled[:, j] for j in range(scaled.shape[1])]


def stratum_masks(df):
    score = df[["KRT5", "KRT6A", "KRT6B"]].mean(axis=1).to_numpy(dtype=float)
    median = np.median(score)
    q33, q66 = np.quantile(score, [1.0 / 3.0, 2.0 / 3.0])
    masks = {
        "all": np.ones(len(df), dtype=bool),
        "krt_high": score >= median,
        "krt_low": score < median,
        "krt_top_tertile": score >= q66,
        "krt_bot_tertile": score <= q33,
    }
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


def append_continuous(bucket, cohort, stratum, cov_name, method, outcome_names, rho, pval, n, k):
    for j, predictor in enumerate(PREDICTORS):
        for i, outcome in enumerate(outcome_names):
            bucket.append((
                "continuous", cohort, stratum, cov_name, method, outcome,
                np.nan, "", predictor, float(rho[j, i]), float(pval[j, i]),
                float(k), int(n), 0, 0,
            ))


def append_quantile(bucket, cohort, stratum, cov_name, method, outcome_names, cut, mode,
                    predictor, delta, pval, var, n, n_high, n_low):
    for i, outcome in enumerate(outcome_names):
        bucket.append((
            "quantile", cohort, stratum, cov_name, method, outcome,
            float(cut), mode, predictor, float(delta[i]), float(pval[i]),
            float(var[i]), int(n), int(n_high), int(n_low),
        ))


def evaluate_block(cohort, stratum, block, cov_items, bucket):
    columns = {col: block[col].to_numpy(dtype=float) for col in block.columns}
    n = len(block)
    if n < MIN_N:
        return
    outcome_names, outcome_raw = build_outcomes(columns)
    pred_raw = np.column_stack([columns[g] for g in PREDICTORS])
    for method in METHODS:
        pred_s = scale_matrix(pred_raw, method)
        out_s = scale_matrix(outcome_raw, method)
        for cov_name, genes, _needs in cov_items:
            if cov_name == "stratum_only" and stratum not in KRT_STRATA:
                continue
            covs = covariate_columns(columns, cov_name, genes, method)
            k = k_of(cov_name, genes)
            if n < k + 5:
                continue
            rho, pval = matrix_partial_corr(pred_s, out_s, covs, n - 2 - k)
            append_continuous(bucket, cohort, stratum, cov_name, method, outcome_names, rho, pval, n, k)
            pred_resid = residualize(pred_s, covs)
            out_resid = residualize(out_s, covs)
            for cut in CUTS:
                for mode, xmat in (("expression", pred_raw), ("residual", pred_resid)):
                    for j, predictor in enumerate(PREDICTORS):
                        low, high = groups(xmat[:, j], cut)
                        n_low = int(low.sum())
                        n_high = int(high.sum())
                        if n_low < MIN_ARM or n_high < MIN_ARM:
                            delta = np.full(len(outcome_names), np.nan)
                            p_nan = np.full(len(outcome_names), np.nan)
                            v_nan = np.full(len(outcome_names), np.nan)
                            append_quantile(
                                bucket, cohort, stratum, cov_name, method, outcome_names,
                                cut, mode, predictor, delta, p_nan, v_nan, n, n_high, n_low,
                            )
                            continue
                        delta, p_d, var = cliffs_batch(out_resid[high], out_resid[low])
                        append_quantile(
                            bucket, cohort, stratum, cov_name, method, outcome_names,
                            cut, mode, predictor, delta, p_d, var, n, n_high, n_low,
                        )


def run_cohort(cohort, frame):
    bucket = []
    masks = stratum_masks(frame)
    print(f"[sweep] {cohort} n={len(frame)} strata={','.join(masks)}", flush=True)
    for stratum, mask in masks.items():
        sub = frame.loc[mask]
        if len(sub) < MIN_N:
            print(f"[sweep] {cohort} {stratum} skip n={len(sub)}", flush=True)
            continue
        plain = [item for item in COVARIATE_SETS if not item[2]]
        with_purity = [item for item in COVARIATE_SETS if item[2]]
        evaluate_block(cohort, stratum, sub, plain, bucket)
        pure = sub[np.isfinite(sub["ESTIMATE"].to_numpy(dtype=float))]
        if len(pure) >= MIN_N:
            evaluate_block(cohort, stratum, pure, with_purity, bucket)
        print(f"[sweep] {cohort} {stratum} n={len(sub)} rows={len(bucket)}", flush=True)
    return bucket


COLUMNS = [
    "family", "cohort", "stratum", "covariates", "method", "outcome",
    "cut", "predictor_cut", "predictor", "effect", "p", "aux", "n", "n_high", "n_low",
]


def summarize_continuous(frame):
    keys = ["stratum", "covariates", "method", "outcome"]
    rows = []
    for key, sub in frame.groupby(keys, sort=False):
        rec = dict(zip(keys, key if isinstance(key, tuple) else (key,)))
        k = int(sub["aux"].iloc[0])
        effects = []
        ps = []
        negs = []
        for predictor in PREDICTORS:
            hit = sub[sub["predictor"] == predictor]
            meta = fisher_random_effects(hit["effect"].to_numpy(), hit["n"].to_numpy(), k)
            rec[f"effect_{predictor}"] = meta["effect"]
            rec[f"fe_{predictor}"] = meta["fe"]
            rec[f"p_{predictor}"] = meta["p"]
            rec[f"I2_{predictor}"] = meta["I2"]
            rec[f"ci_low_{predictor}"] = meta["ci_low"]
            rec[f"ci_high_{predictor}"] = meta["ci_high"]
            rec[f"neg_{predictor}"] = int(np.sum(hit["effect"].to_numpy() < 0))
            rec[f"cohorts_{predictor}"] = int(meta["n_cohorts"])
            effects.append(meta["effect"])
            ps.append(meta["p"])
            negs.append(rec[f"neg_{predictor}"])
        effects = np.asarray(effects, dtype=float)
        rec["min_abs"] = float(np.min(np.abs(effects))) if np.all(np.isfinite(effects)) else np.nan
        rec["mean_abs"] = float(np.mean(np.abs(effects))) if np.all(np.isfinite(effects)) else np.nan
        rec["worst_p"] = float(np.nanmax(ps)) if np.all(np.isfinite(ps)) else np.nan
        rec["k"] = k
        rec["eligible"] = bool(
            np.all(np.isfinite(effects))
            and np.all(effects < 0)
            and all(n >= 7 for n in negs)
            and all(rec[f"cohorts_{p}"] >= 7 for p in PREDICTORS)
        )
        rows.append(rec)
    panel = pd.DataFrame(rows)
    panel["q_worst"] = bh_fdr(panel["worst_p"].to_numpy())
    return panel


def summarize_quantile(frame):
    keys = ["stratum", "covariates", "method", "outcome", "cut", "predictor_cut"]
    rows = []
    for key, sub in frame.groupby(keys, sort=False):
        rec = dict(zip(keys, key if isinstance(key, tuple) else (key,)))
        effects = []
        ps = []
        negs = []
        for predictor in PREDICTORS:
            hit = sub[sub["predictor"] == predictor]
            meta = dl_random_effects(hit["effect"].to_numpy(), hit["aux"].to_numpy())
            finite = hit["effect"].to_numpy()
            rec[f"effect_{predictor}"] = meta["effect"]
            rec[f"fe_{predictor}"] = meta["fe"]
            rec[f"p_{predictor}"] = meta["p"]
            rec[f"I2_{predictor}"] = meta["I2"]
            rec[f"ci_low_{predictor}"] = meta["ci_low"]
            rec[f"ci_high_{predictor}"] = meta["ci_high"]
            rec[f"neg_{predictor}"] = int(np.sum(finite < 0))
            rec[f"cohorts_{predictor}"] = int(meta["n_cohorts"])
            rec[f"median_{predictor}"] = float(np.nanmedian(finite))
            effects.append(meta["effect"])
            ps.append(meta["p"])
            negs.append(rec[f"neg_{predictor}"])
        effects = np.asarray(effects, dtype=float)
        rec["min_abs"] = float(np.min(np.abs(effects))) if np.all(np.isfinite(effects)) else np.nan
        rec["mean_abs"] = float(np.mean(np.abs(effects))) if np.all(np.isfinite(effects)) else np.nan
        rec["worst_p"] = float(np.nanmax(ps)) if np.all(np.isfinite(ps)) else np.nan
        rec["eligible"] = bool(
            np.all(np.isfinite(effects))
            and np.all(effects < 0)
            and all(n >= 7 for n in negs)
            and all(rec[f"cohorts_{p}"] >= 7 for p in PREDICTORS)
        )
        rows.append(rec)
    panel = pd.DataFrame(rows)
    panel["q_worst"] = bh_fdr(panel["worst_p"].to_numpy())
    return panel


def pick_row(panel, mask):
    sub = panel.loc[mask].copy()
    if sub.empty:
        return None
    effects = sub[[f"effect_{p}" for p in PREDICTORS]].to_numpy()
    negs = sub[[f"neg_{p}" for p in PREDICTORS]].to_numpy()
    idx = select_max(effects, negs, sub["worst_p"].to_numpy(), sub["mean_abs"].to_numpy())
    if idx < 0:
        return None
    return sub.iloc[idx]


def fmt(x, digits=3):
    if x is None or not np.isfinite(x):
        return "NA"
    if abs(x) >= 0.001 or x == 0:
        return f"{x:.{digits}f}"
    return f"{x:.2e}"


def fmt_p(x):
    if x is None or not np.isfinite(x):
        return "NA"
    if x < 1e-4:
        return f"{x:.2e}"
    return f"{x:.4f}"


def md_table(records, columns):
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for rec in records:
        cells = []
        for col in columns:
            val = rec[col]
            if isinstance(val, (float, np.floating)):
                as_p = col == "p" or col.startswith("p_") or col.startswith("q") or col == "q_worst" or col.endswith("_p")
                cells.append(fmt_p(float(val)) if as_p else fmt(float(val)))
            else:
                cells.append(str(val))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def winner_dict(row, family):
    if row is None:
        return None
    out = row.to_dict()
    out["family"] = family
    return out


def distribution_stats(panel):
    elig = panel.loc[panel["eligible"], "min_abs"].to_numpy(dtype=float)
    elig = elig[np.isfinite(elig)]
    if elig.size == 0:
        return {"n_eligible": 0}
    return {
        "n_eligible": int(elig.size),
        "median_min_abs": float(np.median(elig)),
        "p90_min_abs": float(np.quantile(elig, 0.90)),
        "p95_min_abs": float(np.quantile(elig, 0.95)),
        "max_min_abs": float(np.max(elig)),
    }


def cohort_slice(detail, spec, family):
    mask = (
        (detail["family"] == family)
        & (detail["stratum"] == spec["stratum"])
        & (detail["covariates"] == spec["covariates"])
        & (detail["method"] == spec["method"])
        & (detail["outcome"] == spec["outcome"])
    )
    if family == "quantile":
        mask = mask & np.isclose(detail["cut"].to_numpy(dtype=float), float(spec["cut"]))
        mask = mask & (detail["predictor_cut"] == spec["predictor_cut"])
    return detail.loc[mask].copy()


def add_ci(slice_df, spec, family):
    rows = []
    for rec in slice_df.itertuples(index=False):
        if family == "continuous":
            lo, hi = fisher_ci(rec.effect, int(rec.n), int(spec["k"]))
        else:
            se = np.sqrt(rec.aux) if np.isfinite(rec.aux) and rec.aux > 0 else np.nan
            if np.isfinite(se) and np.isfinite(rec.effect):
                lo, hi = float(rec.effect - 1.96 * se), float(rec.effect + 1.96 * se)
            else:
                lo, hi = np.nan, np.nan
        rows.append({
            "cohort": rec.cohort,
            "predictor": rec.predictor,
            "effect": rec.effect,
            "p": rec.p,
            "n": rec.n,
            "n_high": rec.n_high,
            "n_low": rec.n_low,
            "ci_low": lo,
            "ci_high": hi,
        })
    return pd.DataFrame(rows)


def write_forest(per_cohort, spec, path, title, xlabel):
    plot = per_cohort.copy()
    plot["cohort"] = pd.Categorical(plot["cohort"], COHORTS, ordered=True)
    plot["predictor"] = pd.Categorical(plot["predictor"], PREDICTORS, ordered=True)
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.8), sharex=True)
    colors = {"TACSTD2": "#1b4f72", "CLDN4": "#b9770e", "CLDN7": "#196f3d"}
    for ax, predictor in zip(axes, PREDICTORS):
        sub = plot[plot["predictor"] == predictor].sort_values("cohort", ascending=False)
        ax.axvline(0, color="#444444", lw=0.8)
        for i, row in enumerate(sub.itertuples(index=False)):
            ax.plot([row.ci_low, row.ci_high], [i, i], color=colors[predictor], lw=1.4)
            ax.plot(row.effect, i, "o", color=colors[predictor], ms=5)
        # Pooled row.
        i = len(sub)
        ax.plot([spec[f"ci_low_{predictor}"], spec[f"ci_high_{predictor}"]], [i, i], color="#111111", lw=1.8)
        ax.plot(spec[f"effect_{predictor}"], i, "D", color="#111111", ms=5)
        labels = [f"{r.cohort}  n={int(r.n)}" for r in sub.itertuples(index=False)] + ["pooled"]
        ax.set_yticks(np.arange(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_title(predictor, color=colors[predictor])
        neg = int(spec[f"neg_{predictor}"])
        ax.set_xlabel(xlabel)
        ax.text(
            0.02, 0.02,
            f"pool {fmt(spec[f'effect_{predictor}'])}\n{neg}/8 neg",
            transform=ax.transAxes, fontsize=8, va="bottom",
        )
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def write_distribution(panel, winner, anchor_min, path, title):
    elig = panel.loc[panel["eligible"], "min_abs"].dropna()
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    if len(elig):
        ax.hist(elig, bins=40, color="#5d6d7e", edgecolor="white")
    if anchor_min is not None and np.isfinite(anchor_min):
        ax.axvline(anchor_min, color="#1b4f72", lw=1.4, label=f"prior CD8 anchor {anchor_min:.3f}")
    if winner is not None:
        ax.axvline(float(winner["min_abs"]), color="#922b21", lw=1.6, label=f"grid max {float(winner['min_abs']):.3f}")
    ax.set_xlabel("minimum |pooled effect| across TACSTD2, CLDN4, CLDN7")
    ax.set_ylabel("eligible specifications")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def spec_label(spec, family):
    if spec is None:
        return "none"
    base = f"{spec['outcome']} | {spec['covariates']} | {spec['stratum']} | {spec['method']}"
    if family == "quantile":
        base += f" | cut {float(spec['cut']):.2f} | {spec['predictor_cut']}"
    return base


def effect_table(spec):
    rows = []
    for predictor in PREDICTORS:
        rows.append({
            "predictor": predictor,
            "pooled": spec[f"effect_{predictor}"],
            "ci_low": spec[f"ci_low_{predictor}"],
            "ci_high": spec[f"ci_high_{predictor}"],
            "p": spec[f"p_{predictor}"],
            "I2": spec[f"I2_{predictor}"],
            "negative_cohorts": spec[f"neg_{predictor}"],
            "cohorts_in_pool": spec[f"cohorts_{predictor}"],
            "fixed_effect": spec[f"fe_{predictor}"],
        })
    return rows


def check_anchor(detail):
    diffs = []
    rows = []
    for (outcome, predictor), expected in ANCHOR_RHO.items():
        for cohort, rho_hat in expected.items():
            hit = detail[
                (detail["family"] == "continuous")
                & (detail["cohort"] == cohort)
                & (detail["stratum"] == "all")
                & (detail["covariates"] == "KRT5_6_14")
                & (detail["method"] == "spearman")
                & (detail["outcome"] == outcome)
                & (detail["predictor"] == predictor)
            ]
            if hit.empty:
                raise RuntimeError(f"anchor row missing {outcome} {predictor} {cohort}")
            got = float(hit["effect"].iloc[0])
            diffs.append(abs(got - rho_hat))
            rows.append({
                "outcome": outcome,
                "predictor": predictor,
                "cohort": cohort,
                "expected_rho": rho_hat,
                "got_rho": got,
                "abs_diff": abs(got - rho_hat),
                "n": int(hit["n"].iloc[0]),
            })
    return pd.DataFrame(rows), float(np.max(diffs))


def top_table(panel, n=12):
    elig = panel.loc[panel["eligible"]].sort_values(
        ["min_abs", "mean_abs", "worst_p"], ascending=[False, False, True]
    ).head(n)
    rows = []
    for rec in elig.to_dict(orient="records"):
        row = {
            "outcome": rec["outcome"],
            "covariates": rec["covariates"],
            "stratum": rec["stratum"],
            "method": rec["method"],
            "min_abs": rec["min_abs"],
            "mean_abs": rec["mean_abs"],
            "TACSTD2": rec["effect_TACSTD2"],
            "CLDN4": rec["effect_CLDN4"],
            "CLDN7": rec["effect_CLDN7"],
            "neg_T2": rec["neg_TACSTD2"],
            "neg_C4": rec["neg_CLDN4"],
            "neg_C7": rec["neg_CLDN7"],
            "worst_p": rec["worst_p"],
            "q_worst": rec["q_worst"],
        }
        if "cut" in rec:
            row["cut"] = rec["cut"]
            row["predictor_cut"] = rec["predictor_cut"]
        rows.append(row)
    return rows


def write_results(path, context):
    c_all = context["cont_all"]
    c_max = context["cont_max"]
    q_q = context["quant_quartile"]
    q_max = context["quant_max"]
    lines = []
    lines.append("# TCGA TACSTD2 / CLDN4 / CLDN7 versus CD3 / CD8 after KRT5/6")
    lines.append("")
    lines.append("Primary tumors only (sample type 01), one row per patient, Xena GDC STAR log2(TPM+1). Eight cohorts: LUAD, LUSC, BRCA, CESC, KIRC, STAD, BLCA, PAAD. This file is written by `scripts/tcga_krt56_max_effect/sweep.py` from the tables next to it.")
    lines.append("")
    lines.append("## Objective")
    lines.append("")
    lines.append("The grid below was declared in the script before it was evaluated. A specification is eligible only when TACSTD2, CLDN4, and CLDN7 each have a negative random-effects pooled effect and each is negative in at least 7 of the 8 cohorts. The reported maximum maximizes the minimum of those three absolute pooled effects. Correlations use a Fisher-z DerSimonian-Laird pool. Quartile and tail contrasts use Cliff's delta (negative means the upper predictor group has the lower immune score), pooled on the delta scale with the Cliff 1993 variance.")
    lines.append("")
    lines.append("The maximum is the extreme of this search. Selecting it pushes |effect| upward relative to a single pre-specified test. `q_worst` is Benjamini-Hochberg on the worst of the three meta p-values across every specification in that family, eligible or not. It does not undo that selection.")
    lines.append("")
    lines.append("Every covariate set contains KRT5 and a KRT6 gene, except `squamous_PC1` (first principal component of KRT5, KRT6A, KRT6B, and KRT14, oriented with KRT5) and `stratum_only` (allowed only inside a KRT5/6 median or tertile split, where the split itself is the KRT5/6 control). Cytotoxic-only gene sets are not in the grid.")
    lines.append("")
    lines.append("## Data")
    lines.append("")
    lines.append(md_table(context["counts"], ["cohort", "n_patients", "n_with_ESTIMATE"]))
    lines.append("")
    lines.append(f"Continuous specifications: {context['n_cont']}. Quantile specifications: {context['n_quant']}. Eligible continuous: {context['dist_cont']['n_eligible']}. Eligible quantile: {context['dist_quant']['n_eligible']}.")
    lines.append("")
    lines.append("Immune scores: CD3D, CD3E, CD3G, CD8A, CD8B; arithmetic mean, within-sample z-mean, rank-mean, and first principal component of CD3 (CD3D/E/G), CD8 (CD8A/B), CD3+CD8, the Danaher T-cell set (CD6, CD3D, CD3E, SH2D1A, TRAT1, CD3G), signature-H (15 T-cell genes), and Tsig (CD3D/E/G, CD2, CD247, LCK). CD8_mean is the two-gene Danaher CD8 score.")
    lines.append("")
    lines.append("Strata, defined inside each cohort: all patients; KRT5/6 score (mean of KRT5, KRT6A, KRT6B) at or above the median, below the median, top tertile, bottom tertile; ESTIMATE purity at or above the median and below the median. ESTIMATE is the Yoshihara cosine transform of the MD Anderson RNAseqV2 ESTIMATE score.")
    lines.append("")
    lines.append("## Pipeline check against the previous KRT5/6 panel")
    lines.append("")
    lines.append(f"Partial Spearman of CD3_mean and CD8_mean on KRT5+KRT6A+KRT6B+KRT14, all patients, was recomputed and compared with the previous grid at 6 decimal places. Maximum absolute difference across 48 cohort-level correlations: {context['anchor_diff']:.3e}.")
    lines.append("")
    lines.append("That previous CD8_mean specification, which is inside this grid, has minimum |pooled ρ| = " + f"{context['anchor_min_abs']:.3f}" + " (TACSTD2, CLDN4, CLDN7 pooled ρ " + ", ".join(fmt(context['anchor_pooled'][p]) for p in PREDICTORS) + "). TACSTD2 is the gene that misses a negative sign, in CESC. CD3_mean on the same covariates has a smaller minimum |ρ| because TACSTD2 is weaker.")
    lines.append("")
    lines.append("## Continuous maximum")
    lines.append("")
    lines.append(f"Among eligible continuous specifications restricted to all patients, the minimum |pooled ρ| has median {fmt(context['dist_cont_all'].get('median_min_abs', np.nan))}, 95th percentile {fmt(context['dist_cont_all'].get('p95_min_abs', np.nan))}, and maximum {fmt(context['dist_cont_all'].get('max_min_abs', np.nan))} ({context['dist_cont_all'].get('n_eligible', 0)} eligible).")
    lines.append("")
    lines.append(f"All-patient maximum: `{spec_label(c_all, 'continuous')}`.")
    lines.append("")
    if c_all is not None:
        lines.append(md_table(effect_table(c_all), ["predictor", "pooled", "ci_low", "ci_high", "p", "I2", "negative_cohorts", "cohorts_in_pool", "fixed_effect"]))
        lines.append("")
        lines.append(f"Worst meta p = {fmt_p(c_all['worst_p'])}. Sweep q = {fmt_p(c_all['q_worst'])}.")
        lines.append("")
        lines.append(md_table(context["cont_all_cohort"], ["cohort", "n", "TACSTD2", "CLDN4", "CLDN7"]))
        lines.append("")
    lines.append(f"Across the full continuous grid, including KRT5/6 and purity strata, {context['dist_cont']['n_eligible']} specifications are eligible. Median minimum |ρ| = {fmt(context['dist_cont'].get('median_min_abs', np.nan))}, 95th percentile = {fmt(context['dist_cont'].get('p95_min_abs', np.nan))}, maximum = {fmt(context['dist_cont'].get('max_min_abs', np.nan))}.")
    lines.append("")
    lines.append(f"Grid maximum: `{spec_label(c_max, 'continuous')}`.")
    lines.append("")
    if c_max is not None:
        lines.append(md_table(effect_table(c_max), ["predictor", "pooled", "ci_low", "ci_high", "p", "I2", "negative_cohorts", "cohorts_in_pool", "fixed_effect"]))
        lines.append("")
        lines.append(f"Worst meta p = {fmt_p(c_max['worst_p'])}. Sweep q = {fmt_p(c_max['q_worst'])}.")
        lines.append("")
        lines.append(md_table(context["cont_max_cohort"], ["cohort", "n", "TACSTD2", "CLDN4", "CLDN7"]))
        lines.append("")
    lines.append("Largest eligible continuous specifications, ordered by the minimum |pooled ρ|:")
    lines.append("")
    if context["cont_top"]:
        cols = ["outcome", "covariates", "stratum", "method", "min_abs", "TACSTD2", "CLDN4", "CLDN7", "neg_T2", "neg_C4", "neg_C7", "worst_p", "q_worst"]
        lines.append(md_table(context["cont_top"], cols))
        lines.append("")
    lines.append("## Cliff's delta maximum")
    lines.append("")
    lines.append("Quartiles are the outer 25% (cut 0.25). The grid also contains outer 20% and outer 10% cuts, because those tail contrasts are where |Cliff's delta| can grow when the association is monotone. Groups are formed either on raw log2(TPM+1) (`expression`) or on the same residual used for the immune score (`residual`). The immune score is always the residual after the named covariates. Arms smaller than 12 patients are not tested.")
    lines.append("")
    lines.append(f"Eligible quartile contrasts in all patients: {context['dist_q_all'].get('n_eligible', 0)}. Median minimum |pooled δ| = {fmt(context['dist_q_all'].get('median_min_abs', np.nan))}, 95th percentile = {fmt(context['dist_q_all'].get('p95_min_abs', np.nan))}, maximum = {fmt(context['dist_q_all'].get('max_min_abs', np.nan))}.")
    lines.append("")
    lines.append(f"All-patient quartile maximum: `{spec_label(q_q, 'quantile')}`.")
    lines.append("")
    if q_q is not None:
        lines.append(md_table(effect_table(q_q), ["predictor", "pooled", "ci_low", "ci_high", "p", "I2", "negative_cohorts", "cohorts_in_pool", "fixed_effect"]))
        lines.append("")
        lines.append(f"Worst meta p = {fmt_p(q_q['worst_p'])}. Sweep q = {fmt_p(q_q['q_worst'])}.")
        lines.append("")
        lines.append(md_table(context["quant_q_cohort"], ["cohort", "n", "n_low", "n_high", "TACSTD2", "CLDN4", "CLDN7"]))
        lines.append("")
    lines.append(f"Across every cut and stratum, {context['dist_quant']['n_eligible']} quantile specifications are eligible. Median minimum |pooled δ| = {fmt(context['dist_quant'].get('median_min_abs', np.nan))}, 95th percentile = {fmt(context['dist_quant'].get('p95_min_abs', np.nan))}, maximum = {fmt(context['dist_quant'].get('max_min_abs', np.nan))}.")
    lines.append("")
    lines.append(f"Grid maximum: `{spec_label(q_max, 'quantile')}`.")
    lines.append("")
    if q_max is not None:
        lines.append(md_table(effect_table(q_max), ["predictor", "pooled", "ci_low", "ci_high", "p", "I2", "negative_cohorts", "cohorts_in_pool", "fixed_effect"]))
        lines.append("")
        lines.append(f"Worst meta p = {fmt_p(q_max['worst_p'])}. Sweep q = {fmt_p(q_max['q_worst'])}.")
        lines.append("")
        lines.append(md_table(context["quant_max_cohort"], ["cohort", "n", "n_low", "n_high", "TACSTD2", "CLDN4", "CLDN7"]))
        lines.append("")
    lines.append("Largest eligible Cliff's delta specifications:")
    lines.append("")
    if context["quant_top"]:
        cols = ["outcome", "covariates", "stratum", "method", "cut", "predictor_cut", "min_abs", "TACSTD2", "CLDN4", "CLDN7", "neg_T2", "neg_C4", "neg_C7", "worst_p", "q_worst"]
        lines.append(md_table(context["quant_top"], cols))
        lines.append("")
    lines.append("## Reading")
    lines.append("")
    lines.append(context["reading"])
    lines.append("")
    lines.append("These are bulk RNA associations after a keratin adjustment. They are not a spatial exclusion measurement, not a knockdown result, and not an estimate of what a single pre-specified test would have returned.")
    lines.append("")
    text = "\n".join(lines)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def wide_cohort(per_cohort):
    rows = []
    for cohort in COHORTS:
        sub = per_cohort[per_cohort["cohort"] == cohort]
        if sub.empty:
            continue
        rec = {
            "cohort": cohort,
            "n": int(sub["n"].iloc[0]),
            "n_low": int(sub["n_low"].iloc[0]) if int(sub["n_low"].iloc[0]) else "",
            "n_high": int(sub["n_high"].iloc[0]) if int(sub["n_high"].iloc[0]) else "",
        }
        for predictor in PREDICTORS:
            hit = sub[sub["predictor"] == predictor]
            rec[predictor] = float(hit["effect"].iloc[0]) if len(hit) else np.nan
        rows.append(rec)
    # Drop empty arm columns for continuous tables in the caller if needed.
    return rows


def reading_paragraph(context):
    c_all = context["cont_all"]
    c_max = context["cont_max"]
    q_q = context["quant_quartile"]
    q_max = context["quant_max"]
    bits = []
    if c_all is not None:
        bits.append(
            "The strongest all-patient partial correlation that keeps all three predictors negative in at least 7 cohorts is "
            f"{spec_label(c_all, 'continuous')}, with pooled ρ "
            + ", ".join(f"{p} {fmt(c_all[f'effect_{p}'])}" for p in PREDICTORS)
            + f" (minimum |ρ| {fmt(c_all['min_abs'])})."
        )
    if c_max is not None and (c_all is None or spec_label(c_max, "continuous") != spec_label(c_all, "continuous")):
        bits.append(
            "Allowing the pre-declared strata raises the minimum |pooled ρ| to "
            f"{fmt(c_max['min_abs'])} at {spec_label(c_max, 'continuous')}."
        )
    elif c_max is not None:
        bits.append("No keratin or purity stratum beat that all-patient continuous specification on the declared objective.")
    if q_q is not None:
        bits.append(
            "The strongest all-patient quartile contrast on the same rule has pooled Cliff's delta "
            + ", ".join(f"{p} {fmt(q_q[f'effect_{p}'])}" for p in PREDICTORS)
            + f" (minimum |δ| {fmt(q_q['min_abs'])}) at {spec_label(q_q, 'quantile')}."
        )
    if q_max is not None:
        bits.append(
            "The largest Cliff's delta in the full cut-and-stratum grid has minimum |δ| "
            f"{fmt(q_max['min_abs'])} at {spec_label(q_max, 'quantile')}."
        )
        if int(q_max.get("cohorts_TACSTD2", 8)) < 8:
            bits.append(
                "That Cliff maximum does not include every cohort: an arm below the pre-set floor of 12 patients is left untested and is not counted as negative."
            )
    runner = context.get("quant_top") or []
    if len(runner) > 1 and q_max is not None:
        second = runner[1]
        if int(second["neg_T2"]) >= 8 and int(second["neg_C4"]) >= 8 and int(second["neg_C7"]) >= 8:
            bits.append(
                "The next Cliff specification keeps a negative sign in all 8 cohorts: "
                f"{second['outcome']} | {second['covariates']} | {second['stratum']} | {second['method']} | "
                f"cut {float(second['cut']):.2f} | {second['predictor_cut']}, minimum |δ| {fmt(second['min_abs'])} "
                f"(TACSTD2 {fmt(second['TACSTD2'])}, CLDN4 {fmt(second['CLDN4'])}, CLDN7 {fmt(second['CLDN7'])})."
            )
    anchor = context["anchor_min_abs"]
    if c_all is not None and np.isfinite(anchor):
        bits.append(
            f"On all patients, the balanced continuous effect moves only from the previous CD8_mean minimum |ρ| of {anchor:.3f} to {float(c_all['min_abs']):.3f}. "
            "The three pooled correlations remain small."
        )
    if c_max is not None and np.isfinite(anchor):
        bits.append(
            f"The continuous grid maximum minimum |ρ| is {float(c_max['min_abs']):.3f}, and it is a KRT5/6-stratum result rather than an all-patient result."
        )
    return " ".join(bits)


def main():
    os.makedirs(TABLES, exist_ok=True)
    os.makedirs(FIGS, exist_ok=True)
    symbols = symbols_needed()
    symbol_to_id = load_probemap(symbols)
    id_to_gene = {ens: gene for gene, ens in symbol_to_id.items()}
    provenance = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(extract_cohort, cohort, id_to_gene) for cohort in COHORTS]
        for fut in as_completed(futures):
            info = fut.result()
            provenance.append(info)
            print(f"[expr] {info['cohort']} n={info['n_patients']} cache={info['from_cache']}", flush=True)
    purity = load_purity()
    frames = {}
    counts = []
    for cohort in COHORTS:
        frame = pd.read_csv(os.path.join(CACHE, f"{cohort}.primary01.tsv.gz"), sep="\t").set_index("patient")
        pur = purity.loc[purity["cohort"] == cohort].set_index("patient")["ESTIMATE"]
        frame = frame.join(pur, how="left")
        frames[cohort] = frame
        counts.append({
            "cohort": cohort,
            "n_patients": int(len(frame)),
            "n_with_ESTIMATE": int(np.isfinite(frame["ESTIMATE"]).sum()),
        })
        print(f"[n] {cohort} {counts[-1]}", flush=True)

    records = []
    for cohort in COHORTS:
        records.extend(run_cohort(cohort, frames[cohort]))
    detail = pd.DataFrame.from_records(records, columns=COLUMNS)
    print(f"[rows] {len(detail)}", flush=True)
    anchor_df, anchor_diff = check_anchor(detail)
    anchor_df.to_csv(os.path.join(TABLES, "anchor_reproduction.tsv"), sep="\t", index=False)
    print(f"[anchor] max abs diff {anchor_diff:.3e}", flush=True)
    if anchor_diff > 5e-4:
        raise SystemExit(f"anchor reproduction failed, max abs diff {anchor_diff}")

    continuous = detail[detail["family"] == "continuous"].copy()
    quantile = detail[detail["family"] == "quantile"].copy()
    cont_panel = summarize_continuous(continuous)
    quant_panel = summarize_quantile(quantile)
    cont_panel.to_csv(os.path.join(TABLES, "continuous_panels.tsv"), sep="\t", index=False)
    quant_panel.to_csv(os.path.join(TABLES, "quantile_panels.tsv"), sep="\t", index=False)
    print(f"[panels] continuous {len(cont_panel)} eligible {int(cont_panel['eligible'].sum())}", flush=True)
    print(f"[panels] quantile {len(quant_panel)} eligible {int(quant_panel['eligible'].sum())}", flush=True)

    cont_all = pick_row(cont_panel, (cont_panel["eligible"]) & (cont_panel["stratum"] == "all"))
    cont_max = pick_row(cont_panel, cont_panel["eligible"])
    q_mask = quant_panel["eligible"] & (quant_panel["stratum"] == "all") & np.isclose(quant_panel["cut"], 0.25)
    quant_quartile = pick_row(quant_panel, q_mask)
    quant_max = pick_row(quant_panel, quant_panel["eligible"])

    def dist_for(panel, mask):
        return distribution_stats(panel.loc[mask].copy() if mask is not None else panel)

    # distribution_stats expects a panel with an eligible column and uses it.
    # For a subset, mark rows outside the subset ineligible by passing a filtered copy
    # that still contains only the subset, with eligible already set.
    dist_cont = distribution_stats(cont_panel)
    dist_quant = distribution_stats(quant_panel)
    dist_cont_all = distribution_stats(cont_panel.loc[cont_panel["stratum"] == "all"].copy())
    dist_q_all = distribution_stats(
        quant_panel.loc[(quant_panel["stratum"] == "all") & np.isclose(quant_panel["cut"], 0.25)].copy()
    )

    winners = []
    slices = {}
    for key, spec, family in (
        ("cont_all", cont_all, "continuous"),
        ("cont_max", cont_max, "continuous"),
        ("quant_quartile", quant_quartile, "quantile"),
        ("quant_max", quant_max, "quantile"),
    ):
        if spec is None:
            slices[key] = pd.DataFrame()
            continue
        raw = cohort_slice(detail, spec, family)
        raw.to_csv(os.path.join(TABLES, f"{key}_by_cohort.tsv"), sep="\t", index=False)
        plotted = add_ci(raw, spec, family)
        slices[key] = plotted
        winners.append(winner_dict(spec, family) | {"role": key})
        xlabel = "partial ρ" if family == "continuous" else "Cliff's delta"
        write_forest(
            plotted, spec,
            os.path.join(FIGS, f"{key}_forest.png"),
            spec_label(spec, family),
            xlabel,
        )
    pd.DataFrame(winners).to_csv(os.path.join(TABLES, "winners.tsv"), sep="\t", index=False)

    anchor_hit = cont_panel[
        (cont_panel["stratum"] == "all")
        & (cont_panel["covariates"] == "KRT5_6_14")
        & (cont_panel["method"] == "spearman")
        & (cont_panel["outcome"] == "CD8_mean")
    ]
    if anchor_hit.empty:
        raise SystemExit("anchor panel row missing")
    anchor_row = anchor_hit.iloc[0]
    anchor_min_abs = float(anchor_row["min_abs"])
    anchor_pooled = {p: float(anchor_row[f"effect_{p}"]) for p in PREDICTORS}
    if cont_all is not None and float(cont_all["min_abs"]) + 1e-9 < anchor_min_abs:
        raise SystemExit("all-patient maximum is below the CD8 anchor; selection bug")

    write_distribution(
        cont_panel, cont_max, anchor_min_abs,
        os.path.join(FIGS, "continuous_min_abs_distribution.png"),
        "Eligible continuous specifications",
    )
    write_distribution(
        quant_panel, quant_max, None,
        os.path.join(FIGS, "quantile_min_abs_distribution.png"),
        "Eligible Cliff's delta specifications",
    )

    def cohort_rows(plotted, with_arms):
        rows = wide_cohort(plotted)
        if not with_arms:
            for rec in rows:
                rec.pop("n_low", None)
                rec.pop("n_high", None)
        return rows

    context = {
        "counts": counts,
        "n_cont": int(len(cont_panel)),
        "n_quant": int(len(quant_panel)),
        "dist_cont": dist_cont,
        "dist_quant": dist_quant,
        "dist_cont_all": dist_cont_all,
        "dist_q_all": dist_q_all,
        "cont_all": None if cont_all is None else cont_all.to_dict(),
        "cont_max": None if cont_max is None else cont_max.to_dict(),
        "quant_quartile": None if quant_quartile is None else quant_quartile.to_dict(),
        "quant_max": None if quant_max is None else quant_max.to_dict(),
        "cont_all_cohort": [] if cont_all is None else cohort_rows(slices["cont_all"], False),
        "cont_max_cohort": [] if cont_max is None else cohort_rows(slices["cont_max"], False),
        "quant_q_cohort": [] if quant_quartile is None else cohort_rows(slices["quant_quartile"], True),
        "quant_max_cohort": [] if quant_max is None else cohort_rows(slices["quant_max"], True),
        "cont_top": top_table(cont_panel),
        "quant_top": top_table(quant_panel),
        "anchor_diff": anchor_diff,
        "anchor_min_abs": anchor_min_abs,
        "anchor_pooled": anchor_pooled,
    }
    context["reading"] = reading_paragraph(context)
    write_results(os.path.join(OUT, "RESULTS.md"), context)
    meta = {
        "cohorts": provenance,
        "sample_counts": counts,
        "gene_ids": symbol_to_id,
        "anchor_max_abs_diff": anchor_diff,
        "anchor_min_abs_cd8": anchor_min_abs,
        "anchor_pooled_cd8": anchor_pooled,
        "n_continuous_specs": int(len(cont_panel)),
        "n_quantile_specs": int(len(quant_panel)),
        "continuous_distribution": dist_cont,
        "quantile_distribution": dist_quant,
        "winners": {
            "cont_all": spec_label(context["cont_all"], "continuous"),
            "cont_max": spec_label(context["cont_max"], "continuous"),
            "quant_quartile": spec_label(context["quant_quartile"], "quantile"),
            "quant_max": spec_label(context["quant_max"], "quantile"),
        },
    }
    with open(os.path.join(TABLES, "provenance.json"), "w", encoding="utf-8") as handle:
        json.dump(meta, handle, indent=2, default=float)
    art = "/opt/cursor/artifacts/tcga_krt56"
    os.makedirs(art, exist_ok=True)
    for name in os.listdir(FIGS):
        shutil.copy2(os.path.join(FIGS, name), os.path.join(art, name))
    shutil.copy2(os.path.join(OUT, "RESULTS.md"), os.path.join(art, "RESULTS.md"))
    print("[done]", meta["winners"], flush=True)
    print(context["reading"], flush=True)


if __name__ == "__main__":
    main()
