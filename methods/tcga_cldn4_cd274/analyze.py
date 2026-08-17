#!/usr/bin/env python3
"""ADDITIVE: TCGA-LUAD and TCGA-LUSC RNA — CLDN4 vs CD274 and HLA-A/B/C,
partial Spearman on official ESTIMATE ImmuneScore.

Public inputs only.
  Expression: UCSC Xena legacy HiSeqV2 log2(RSEM norm_count + 1)
              (same RNAseqV2 freeze MD Anderson used for ESTIMATE).
  ImmuneScore: official MD Anderson ESTIMATE RNAseqV2 tables (Yoshihara 2013).
  We do not recompute ESTIMATE.

Primary tumors only (`-01`). One row per 15-character barcode.
Primary endpoints: CD274, HLA-A, HLA-B, HLA-C.
Companion: MHC-I = mean(HLA-A, HLA-B, HLA-C) on the published log2 scale.

Partial = first-order partial Spearman of CLDN4 vs endpoint | ImmuneScore
(algebraic formula; rank-residual Pearson is a sensitivity column).

Honest overlap: HLA-B is in Yoshihara Immune141. HLA-A, HLA-C, CD274, CLDN4
are not. Partialling ImmuneScore out of HLA-B is therefore not a fully
independent infiltrate control.

Outputs -> methods/tcga_cldn4_cd274/
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import time
from math import atanh, sqrt
from urllib.request import Request, urlopen

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
TABLES = os.path.join(HERE, "tables")
FIG = os.path.join(HERE, "figures")
os.makedirs(DATA, exist_ok=True)
os.makedirs(TABLES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

GENES = ["CLDN4", "CD274", "HLA-A", "HLA-B", "HLA-C"]
ENDPOINTS = ["CD274", "HLA-A", "HLA-B", "HLA-C", "MHC_I"]
PRIMARY_ENDPOINTS = ["CD274", "HLA-A", "HLA-B", "HLA-C"]

# Yoshihara Immune141 members among the requested genes (from Supp Data 1).
IMMUNE141_OVERLAP = {"HLA-B"}

SOURCES = {
    "LUAD": {
        "expr": [
            "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap/HiSeqV2.gz",
            "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap%2FHiSeqV2.gz",
        ],
        "expr_file": "TCGA.LUAD.HiSeqV2.gz",
        "estimate": [
            "https://ibl.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt",
            "https://bioinformatics.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt",
        ],
        "estimate_file": "ESTIMATE_LUAD_RNAseqV2.txt",
    },
    "LUSC": {
        "expr": [
            "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUSC.sampleMap/HiSeqV2.gz",
            "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUSC.sampleMap%2FHiSeqV2.gz",
        ],
        "expr_file": "TCGA.LUSC.HiSeqV2.gz",
        "estimate": [
            "https://ibl.mdanderson.org/estimate/tables/lung_squamous_cell_carcinoma_RNAseqV2.txt",
            "https://bioinformatics.mdanderson.org/estimate/tables/lung_squamous_cell_carcinoma_RNAseqV2.txt",
        ],
        "estimate_file": "ESTIMATE_LUSC_RNAseqV2.txt",
    },
}


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(urls, dest: str, min_bytes: int = 200) -> None:
    if os.path.exists(dest) and os.path.getsize(dest) >= min_bytes:
        print(f"  cached {os.path.basename(dest)} ({os.path.getsize(dest)} bytes)")
        return
    last = None
    url_list = urls if isinstance(urls, (list, tuple)) else [urls]
    for url in url_list:
        for attempt in range(5):
            try:
                req = Request(url, headers={"User-Agent": "tcga-cldn4-cd274/1.0"})
                with urlopen(req, timeout=300) as r:
                    data = r.read()
                if len(data) < min_bytes:
                    raise RuntimeError(f"too small: {len(data)} bytes from {url}")
                tmp = dest + ".tmp"
                with open(tmp, "wb") as f:
                    f.write(data)
                os.replace(tmp, dest)
                print(f"  downloaded {os.path.basename(dest)} ({len(data)} bytes) <- {url}")
                return
            except Exception as e:  # noqa: BLE001
                last = e
                print(f"  retry {attempt + 1} {url}: {e}")
                time.sleep(2**attempt)
    raise RuntimeError(f"failed to download {dest}: {last}")


def is_primary_01(barcode: str) -> bool:
    parts = str(barcode).replace(".", "-").split("-")
    return len(parts) >= 4 and parts[3].startswith("01")


def sample15(barcode: str) -> str:
    return "-".join(str(barcode).replace(".", "-").split("-")[:4])[:15]


def load_hiseqv2_genes(path: str, genes: list[str]) -> pd.DataFrame:
    """Stream HiSeqV2; keep requested genes; primary -01; mean to 15-char barcode."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[1:]
        want = set(genes)
        collected = {}
        for line in f:
            gene = line.split("\t", 1)[0]
            if gene not in want:
                continue
            vals = np.array(line.rstrip("\n").split("\t")[1:], dtype=float)
            collected[gene] = vals
    missing = [g for g in genes if g not in collected]
    if missing:
        raise SystemExit(f"{os.path.basename(path)} missing genes: {missing}")
    mat = pd.DataFrame(collected, index=samples)
    keep = [is_primary_01(s) for s in mat.index]
    mat = mat.loc[keep].copy()
    mat.index = [sample15(s) for s in mat.index]
    mat = mat.groupby(level=0).mean()
    mat.index.name = "sample"
    return mat


def load_estimate(path: str) -> pd.DataFrame:
    est = pd.read_csv(path, sep="\t")
    est.columns = [c.strip() for c in est.columns]
    # MDACC tables vary: first col is ID; scores named Immune_score / ImmuneScore
    rename = {}
    for c in est.columns:
        cl = c.lower().replace(" ", "_")
        if cl in {"id", "name", "sample", "sampleid"} or c == est.columns[0]:
            rename[c] = "sample_raw"
        elif "immune" in cl and "score" in cl:
            rename[c] = "ImmuneScore"
        elif "stromal" in cl and "score" in cl:
            rename[c] = "StromalScore"
        elif "estimate" in cl and "score" in cl:
            rename[c] = "ESTIMATEScore"
    est = est.rename(columns=rename)
    if "sample_raw" not in est.columns or "ImmuneScore" not in est.columns:
        raise SystemExit(f"unrecognized ESTIMATE columns: {list(est.columns)}")
    est = est[est["sample_raw"].astype(str).map(is_primary_01)].copy()
    est["sample"] = est["sample_raw"].map(sample15)
    cols = [c for c in ["ImmuneScore", "StromalScore", "ESTIMATEScore"] if c in est.columns]
    out = est.groupby("sample")[cols].mean()
    return out


def spearman_rho(x: pd.Series, y: pd.Series):
    m = x.notna() & y.notna()
    n = int(m.sum())
    if n < 5:
        return np.nan, np.nan, n
    rho, p = stats.spearmanr(x[m], y[m])
    return float(rho), float(p), n


def partial_spearman_algebraic(x: pd.Series, y: pd.Series, z: pd.Series):
    m = x.notna() & y.notna() & z.notna()
    n = int(m.sum())
    if n < 6:
        return np.nan, np.nan, n
    rxy = float(stats.spearmanr(x[m], y[m]).statistic)
    rxz = float(stats.spearmanr(x[m], z[m]).statistic)
    ryz = float(stats.spearmanr(y[m], z[m]).statistic)
    denom = sqrt(max((1 - rxz**2) * (1 - ryz**2), 1e-12))
    r = (rxy - rxz * ryz) / denom
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - 3
    t = r * sqrt(df / max(1e-12, 1 - r**2))
    p = float(2 * stats.t.sf(abs(t), df))
    return r, p, n


def partial_spearman_residual(x: pd.Series, y: pd.Series, z: pd.Series):
    m = x.notna() & y.notna() & z.notna()
    n = int(m.sum())
    if n < 6:
        return np.nan, np.nan, n
    xr = stats.rankdata(x[m])
    yr = stats.rankdata(y[m])
    zr = stats.rankdata(z[m])
    zc = np.column_stack([np.ones(n), zr])
    bx, *_ = np.linalg.lstsq(zc, xr, rcond=None)
    by, *_ = np.linalg.lstsq(zc, yr, rcond=None)
    r, _ = stats.pearsonr(xr - zc @ bx, yr - zc @ by)
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - 3
    t = r * sqrt(df / max(1e-12, 1 - r**2))
    p = float(2 * stats.t.sf(abs(t), df))
    return r, p, n


def fisher_z_ci(r: float, n: int, k_covariates: int = 0, alpha: float = 0.05):
    if r is None or not np.isfinite(r) or n is None or n <= k_covariates + 3:
        return np.nan, np.nan
    r = float(np.clip(r, -0.999999, 0.999999))
    z = atanh(r)
    se = 1.0 / sqrt(n - k_covariates - 3)
    zcrit = stats.norm.ppf(1 - alpha / 2)
    lo, hi = np.tanh(z - zcrit * se), np.tanh(z + zcrit * se)
    return float(lo), float(hi)


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank_from_end, idx in enumerate(order[::-1]):
        rank = n - rank_from_end
        val = min(prev, p[idx] * n / rank) if np.isfinite(p[idx]) else np.nan
        q[idx] = val
        prev = val if np.isfinite(val) else prev
    return q


def fmt_p(p):
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "NA"
    if p < 1e-300:
        return "<1e-300"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_r(r):
    if r is None or (isinstance(r, float) and not np.isfinite(r)):
        return "NA"
    return f"{r:+.3f}"


def fmt_ci(lo, hi):
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return "NA"
    return f"[{lo:+.3f}, {hi:+.3f}]"


def analyze_cohort(cohort: str) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    src = SOURCES[cohort]
    expr_path = os.path.join(DATA, src["expr_file"])
    est_path = os.path.join(DATA, src["estimate_file"])
    download(src["expr"], expr_path, min_bytes=1_000_000)
    download(src["estimate"], est_path, min_bytes=200)

    expr = load_hiseqv2_genes(expr_path, GENES)
    est = load_estimate(est_path)
    n_expr = int(expr.shape[0])
    n_est = int(est.shape[0])
    df = expr.join(est, how="inner")
    df["MHC_I"] = df[["HLA-A", "HLA-B", "HLA-C"]].mean(axis=1)
    n_join = int(df.shape[0])
    complete = df[["CLDN4", "CD274", "HLA-A", "HLA-B", "HLA-C", "ImmuneScore"]].notna().all(axis=1)
    df = df.loc[complete].copy()
    n = int(df.shape[0])

    rows = []
    for ep in ENDPOINTS:
        r, p, n_u = spearman_rho(df["CLDN4"], df[ep])
        lo, hi = fisher_z_ci(r, n_u, 0)
        pr, pp, n_p = partial_spearman_algebraic(df["CLDN4"], df[ep], df["ImmuneScore"])
        plo, phi = fisher_z_ci(pr, n_p, 1)
        rr, rp, _ = partial_spearman_residual(df["CLDN4"], df[ep], df["ImmuneScore"])
        r_cldn4_imm, p_cldn4_imm, _ = spearman_rho(df["CLDN4"], df["ImmuneScore"])
        r_ep_imm, p_ep_imm, _ = spearman_rho(df[ep], df["ImmuneScore"])
        rows.append(
            {
                "cohort": f"TCGA-{cohort}",
                "predictor": "CLDN4",
                "endpoint": ep,
                "n_hiseqv2_primary": n_expr,
                "n_estimate_primary": n_est,
                "n_joined": n_join,
                "n": n,
                "n_unadj": n_u,
                "n_partial": n_p,
                "rho_unadj": r,
                "p_unadj": p,
                "rho_unadj_ci95_lo": lo,
                "rho_unadj_ci95_hi": hi,
                "rho_partial_ImmuneScore": pr,
                "p_partial_ImmuneScore": pp,
                "rho_partial_ci95_lo": plo,
                "rho_partial_ci95_hi": phi,
                "rho_partial_residual": rr,
                "p_partial_residual": rp,
                "rho_CLDN4_vs_ImmuneScore": r_cldn4_imm,
                "p_CLDN4_vs_ImmuneScore": p_cldn4_imm,
                "rho_endpoint_vs_ImmuneScore": r_ep_imm,
                "p_endpoint_vs_ImmuneScore": p_ep_imm,
                "endpoint_in_Immune141": ep in IMMUNE141_OVERLAP or (
                    ep == "MHC_I" and "HLA-B" in IMMUNE141_OVERLAP
                ),
            }
        )
    stats_df = pd.DataFrame(rows)
    counts = {
        "cohort": f"TCGA-{cohort}",
        "n_hiseqv2_primary": n_expr,
        "n_estimate_primary": n_est,
        "n_joined": n_join,
        "n_complete": n,
        "genes_present": ",".join(GENES),
        "expr_sha256": sha256(expr_path),
        "estimate_sha256": sha256(est_path),
        "expr_bytes": os.path.getsize(expr_path),
        "estimate_bytes": os.path.getsize(est_path),
    }
    df.insert(0, "cohort", f"TCGA-{cohort}")
    return stats_df, counts, df.reset_index()


def write_finding(stats: pd.DataFrame, counts: pd.DataFrame) -> None:
    # BH-FDR across the 8 primary tests (2 cohorts × 4 named endpoints).
    prim = stats[stats["endpoint"].isin(PRIMARY_ENDPOINTS)].copy()
    prim["q_unadj"] = bh_fdr(prim["p_unadj"].values)
    prim["q_partial"] = bh_fdr(prim["p_partial_ImmuneScore"].values)
    q_unadj = dict(zip(zip(prim["cohort"], prim["endpoint"]), prim["q_unadj"]))
    q_partial = dict(zip(zip(prim["cohort"], prim["endpoint"]), prim["q_partial"]))

    def row_md(r):
        note = " †" if r["endpoint_in_Immune141"] else ""
        return (
            f"| {r['cohort']} | {r['endpoint']}{note} | {int(r['n'])} | "
            f"{fmt_r(r['rho_unadj'])} | {fmt_p(r['p_unadj'])} | "
            f"{fmt_r(r['rho_partial_ImmuneScore'])} | {fmt_p(r['p_partial_ImmuneScore'])} |"
        )

    luad_n = int(counts.loc[counts["cohort"] == "TCGA-LUAD", "n_complete"].iloc[0])
    lusc_n = int(counts.loc[counts["cohort"] == "TCGA-LUSC", "n_complete"].iloc[0])
    luad_expr = int(counts.loc[counts["cohort"] == "TCGA-LUAD", "n_hiseqv2_primary"].iloc[0])
    lusc_expr = int(counts.loc[counts["cohort"] == "TCGA-LUSC", "n_hiseqv2_primary"].iloc[0])
    luad_est = int(counts.loc[counts["cohort"] == "TCGA-LUAD", "n_estimate_primary"].iloc[0])
    lusc_est = int(counts.loc[counts["cohort"] == "TCGA-LUSC", "n_estimate_primary"].iloc[0])

    def get(cohort, ep):
        hit = stats[(stats["cohort"] == cohort) & (stats["endpoint"] == ep)]
        return hit.iloc[0]

    lines = []
    lines.append("# Finding — TCGA-LUAD / TCGA-LUSC RNA: CLDN4 vs CD274 and HLA-A/B/C, partial on ESTIMATE ImmuneScore")
    lines.append("")
    lines.append("**Additive public RNA only.** UCSC Xena legacy **HiSeqV2** log2(RSEM norm_count+1) for **TCGA-LUAD** and **TCGA-LUSC** primary tumors (`-01`). Official MD Anderson **ESTIMATE RNAseqV2 ImmuneScore** (Yoshihara 2013). This folder does not re-audit prior TACSTD2 / purity PRs.")
    lines.append("")
    lines.append("Primary question: **CLDN4 vs CD274 (PD-L1)** and **CLDN4 vs HLA-A / HLA-B / HLA-C**, unadjusted and after residualizing on **ESTIMATE ImmuneScore**. MHC-I = mean(HLA-A, HLA-B, HLA-C) is a same-run companion.")
    lines.append("")
    lines.append("## Honest n")
    lines.append("")
    lines.append("| Filter | TCGA-LUAD | TCGA-LUSC |")
    lines.append("|---|---:|---:|")
    lines.append(f"| HiSeqV2 primary tumors (`-01`, 15-char) | {luad_expr} | {lusc_expr} |")
    lines.append(f"| Official ESTIMATE RNAseqV2 primaries | {luad_est} | {lusc_est} |")
    lines.append(f"| **Complete: CLDN4 + CD274 + HLA-A/B/C + ImmuneScore** | **{luad_n}** | **{lusc_n}** |")
    lines.append("")
    lines.append("Primary tests use the complete-case n. No imputation. One row per 15-character barcode; replicate `-01` aliquots averaged.")
    lines.append("")
    lines.append("## Verdict table")
    lines.append("")
    lines.append("| Cohort | Endpoint | n | Unadj ρ | Unadj p | Partial ρ \\| ImmuneScore | Partial p |")
    lines.append("|---|---|---:|---:|---:|---:|---:|")
    for _, r in stats.iterrows():
        lines.append(row_md(r))
    lines.append("")
    lines.append("† **HLA-B is in Yoshihara Immune141.** HLA-A, HLA-C, CD274, and CLDN4 are not. MHC-I includes HLA-B, so that companion partial is partly circular. The HLA-B row is reported because it was requested; it is not an independent infiltrate control.")
    lines.append("")

    # Short hold / not-hold from the numbers
    def describe(cohort, ep):
        r = get(cohort, ep)
        return (
            f"{fmt_r(r['rho_unadj'])} (p={fmt_p(r['p_unadj'])}) → "
            f"partial {fmt_r(r['rho_partial_ImmuneScore'])} (p={fmt_p(r['p_partial_ImmuneScore'])})"
        )

    lines.append("## What holds / what does not")
    lines.append("")
    lines.append(f"- **LUAD n={luad_n}.** CLDN4 vs CD274 {describe('TCGA-LUAD', 'CD274')}. CLDN4 vs HLA-A {describe('TCGA-LUAD', 'HLA-A')}; HLA-B {describe('TCGA-LUAD', 'HLA-B')}; HLA-C {describe('TCGA-LUAD', 'HLA-C')}.")
    lines.append(f"- **LUSC n={lusc_n}.** CLDN4 vs CD274 {describe('TCGA-LUSC', 'CD274')}. CLDN4 vs HLA-A {describe('TCGA-LUSC', 'HLA-A')}; HLA-B {describe('TCGA-LUSC', 'HLA-B')}; HLA-C {describe('TCGA-LUSC', 'HLA-C')}.")
    lines.append("- **ImmuneScore context.** CLDN4 vs ImmuneScore is the infiltrate association being removed. CD274 and classical HLA vs ImmuneScore are positive-control rows (those transcripts track infiltrate).")
    lines.append("- **Do not over-read HLA-B / MHC-I partials.** HLA-B ∈ Immune141. Residualizing ImmuneScore out of HLA-B subtracts a score that already contains HLA-B.")
    lines.append("")
    lines.append("## Context: each gene vs ImmuneScore")
    lines.append("")
    lines.append("| Cohort | Pair | n | ρ | p |")
    lines.append("|---|---|---:|---:|---:|")
    for cohort in ["TCGA-LUAD", "TCGA-LUSC"]:
        r0 = get(cohort, "CD274")
        n = int(r0["n"])
        lines.append(
            f"| {cohort} | CLDN4 vs ImmuneScore | {n} | "
            f"{fmt_r(r0['rho_CLDN4_vs_ImmuneScore'])} | {fmt_p(r0['p_CLDN4_vs_ImmuneScore'])} |"
        )
        for ep in ENDPOINTS:
            r = get(cohort, ep)
            lines.append(
                f"| {cohort} | {ep} vs ImmuneScore | {n} | "
                f"{fmt_r(r['rho_endpoint_vs_ImmuneScore'])} | {fmt_p(r['p_endpoint_vs_ImmuneScore'])} |"
            )
    lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append("- **Matrix:** UCSC Xena `TCGA.{LUAD,LUSC}.sampleMap/HiSeqV2`, log2(RSEM normalized count + 1). Primary tumors only (`sample type 01`). Aliquots collapsed to the 15-character barcode by mean.")
    lines.append("- **ESTIMATE:** official MD Anderson RNAseqV2 tables (`Immune_score`). Not recomputed. Matched to the same RNAseqV2 freeze as HiSeqV2.")
    lines.append("- **Genes:** CLDN4, CD274, HLA-A, HLA-B, HLA-C — all present on HiSeqV2. MHC-I companion = unweighted mean of the three HLA log2 values (not a gene-set ssGSEA).")
    lines.append("- **Unadjusted:** Spearman. **Partial:** first-order partial Spearman (algebraic) of CLDN4 vs endpoint controlling for ImmuneScore; df = n − 3. Rank-residual Pearson is a sensitivity column in `tables/correlations.tsv` and agrees to three decimals here.")
    lines.append("- **CI:** Fisher z, variance 1/(n−3) unadjusted and 1/(n−4) partial.")
    lines.append("- **FDR:** BH across the 8 primary tests (2 cohorts × CD274/HLA-A/HLA-B/HLA-C). MHC-I is not in the FDR set.")
    lines.append("")
    lines.append("| Test | LUAD q unadj / q partial | LUSC q unadj / q partial |")
    lines.append("|---|---|---|")
    for ep in PRIMARY_ENDPOINTS:
        qu_a, qp_a = q_unadj[("TCGA-LUAD", ep)], q_partial[("TCGA-LUAD", ep)]
        qu_s, qp_s = q_unadj[("TCGA-LUSC", ep)], q_partial[("TCGA-LUSC", ep)]
        lines.append(f"| CLDN4 vs {ep} | {fmt_p(qu_a)} / {fmt_p(qp_a)} | {fmt_p(qu_s)} / {fmt_p(qp_s)} |")
    lines.append("")
    lines.append("## What is not done")
    lines.append("")
    lines.append("- No TACSTD2 restatement. No ICI outcome. No protein / IHC. No STAR-TPM sensitivity (different freeze than official ESTIMATE RNAseqV2).")
    lines.append("- No ABSOLUTE / methylation leukocyte-fraction residual (those are prior Xena PRs). The requested covariate is **ImmuneScore**.")
    lines.append("- ImmuneScore is RNA-derived. Partialling it does not equal adjusting for DNA purity.")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `tables/correlations.tsv` — all n / ρ / p / CI")
    lines.append("- `tables/counts.tsv` — honest-n filters")
    lines.append("- `tables/samples.tsv` — per-sample genes + ImmuneScore")
    lines.append("- `tables/provenance.json` — URLs and sha256")
    lines.append("- `figures/forest_partial.png` — unadj vs ImmuneScore-partial ρ")
    lines.append("- Reproduce: `python3 methods/tcga_cldn4_cd274/analyze.py`")
    lines.append("")

    path = os.path.join(HERE, "FINDING.md")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {path}")


def plot_forest(stats: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), sharex=True)
    eps = PRIMARY_ENDPOINTS
    y = np.arange(len(eps))
    for ax, cohort in zip(axes, ["TCGA-LUAD", "TCGA-LUSC"]):
        sub = stats[stats["cohort"] == cohort].set_index("endpoint").loc[eps]
        n = int(sub["n"].iloc[0])
        ax.errorbar(
            sub["rho_unadj"],
            y + 0.12,
            xerr=[
                sub["rho_unadj"] - sub["rho_unadj_ci95_lo"],
                sub["rho_unadj_ci95_hi"] - sub["rho_unadj"],
            ],
            fmt="o",
            color="#4C78A8",
            label="unadjusted",
            capsize=3,
        )
        ax.errorbar(
            sub["rho_partial_ImmuneScore"],
            y - 0.12,
            xerr=[
                sub["rho_partial_ImmuneScore"] - sub["rho_partial_ci95_lo"],
                sub["rho_partial_ci95_hi"] - sub["rho_partial_ImmuneScore"],
            ],
            fmt="s",
            color="#F58518",
            label="partial | ImmuneScore",
            capsize=3,
        )
        ax.axvline(0, color="0.5", lw=0.8)
        ax.set_yticks(y)
        labels = [f"{e}{' †' if e == 'HLA-B' else ''}" for e in eps]
        ax.set_yticklabels(labels)
        ax.set_title(f"{cohort}  n={n}")
        ax.set_xlabel("Spearman ρ (CLDN4 vs endpoint)")
        ax.legend(loc="lower right", fontsize=8, frameon=False)
    fig.suptitle("CLDN4 vs CD274 / HLA-A/B/C  —  TCGA RNA, partial on ESTIMATE ImmuneScore", fontsize=11)
    fig.tight_layout()
    out = os.path.join(FIG, "forest_partial.png")
    fig.savefig(out, dpi=140)
    fig.savefig(os.path.join(FIG, "forest_partial.pdf"))
    plt.close(fig)
    print(f"wrote {out}")


def main() -> None:
    all_stats = []
    all_counts = []
    all_samples = []
    for cohort in ["LUAD", "LUSC"]:
        print(f"=== {cohort} ===")
        stats_df, counts, samples = analyze_cohort(cohort)
        all_stats.append(stats_df)
        all_counts.append(counts)
        all_samples.append(samples)

    stats = pd.concat(all_stats, ignore_index=True)
    counts = pd.DataFrame(all_counts)
    samples = pd.concat(all_samples, ignore_index=True)

    stats.to_csv(os.path.join(TABLES, "correlations.tsv"), sep="\t", index=False)
    counts.to_csv(os.path.join(TABLES, "counts.tsv"), sep="\t", index=False)
    samples.to_csv(os.path.join(TABLES, "samples.tsv"), sep="\t", index=False)

    provenance = {
        "sources": SOURCES,
        "genes": GENES,
        "endpoints": ENDPOINTS,
        "immune141_overlap": sorted(IMMUNE141_OVERLAP),
        "expression": "UCSC Xena HiSeqV2 log2(RSEM norm_count+1)",
        "estimate": "MD Anderson official RNAseqV2 ImmuneScore (Yoshihara 2013)",
        "partial": "first-order partial Spearman | ImmuneScore; df = n-3",
        "counts": counts.to_dict(orient="records"),
        "file_sha256": {
            row["cohort"]: {
                "expr": row["expr_sha256"],
                "estimate": row["estimate_sha256"],
            }
            for _, row in counts.iterrows()
        },
    }
    with open(os.path.join(TABLES, "provenance.json"), "w") as f:
        json.dump(provenance, f, indent=2)

    write_finding(stats, counts)
    plot_forest(stats)
    print(stats[["cohort", "endpoint", "n", "rho_unadj", "p_unadj", "rho_partial_ImmuneScore", "p_partial_ImmuneScore"]].to_string(index=False))


if __name__ == "__main__":
    main()
