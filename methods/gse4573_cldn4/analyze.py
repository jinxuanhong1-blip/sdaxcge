#!/usr/bin/env python3
"""GSE4573 LUSC Affy U133A: CLDN4 vs CD8A / immune genes.

TACSTD2 TLS/CD8 extras from PR 229 are treated as given and are not the
claim of this slice. This script scores CLDN4 on the public series matrix
and writes honest-n tables.

Downloads stay under $GSE4573_CLDN4_DATA (default /tmp/gse4573_cldn4)
and are not committed.
"""

from __future__ import annotations

import gzip
import io
import json
import os
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
DATA = Path(os.environ.get("GSE4573_CLDN4_DATA", "/tmp/gse4573_cldn4"))
SEED = 20260816
N_BOOT = 2000

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE4nnn/GSE4573/"
    "matrix/GSE4573_series_matrix.txt.gz"
)
GPL_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz"

# Pre-specified single-gene immune panel. CD8A is the primary partner.
IMMUNE_GENES = [
    "CD8A",
    "CD8B",
    "CD2",
    "CD3D",
    "CD3E",
    "GZMA",
    "GZMB",
    "GZMK",
    "PRF1",
    "IFNG",
    "STAT1",
    "CXCL9",
    "CXCL10",
    "IDO1",
    "HLA-DRA",
    "CXCL13",
    "CCL19",
    "CCL21",
    "MS4A1",
    "CD19",
    "CD79A",
    "CD79B",
]

# Same locked lists as PR 229 (scripts/opus_tls/opus_tls_lib.py).
SIGNATURES = {
    "TLS_Cabrita": ["CCL19", "CCL21", "CXCL13", "CCL17", "CCR7", "CXCR5", "SELL", "LAMP3"],
    "TLS_12chemokine": [
        "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18",
        "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13",
    ],
    "TLS_imprint": ["CD79B", "CD1D", "CCR6", "LAT", "SKAP1", "CETP", "EIF1AY", "RBP5", "PTGDS"],
    "B_cell": [
        "MS4A1", "CD19", "CD79A", "CD79B", "BANK1", "BLK", "CD22", "CR2",
        "FCRL2", "FCRL5", "PAX5", "TCL1A", "TNFRSF13C", "VPREB3", "CD37", "IRF8",
    ],
    "Plasma_cell": [
        "MZB1", "JCHAIN", "DERL3", "TNFRSF17", "POU2AF1", "XBP1", "PRDM1",
        "SDC1", "IGHG1", "IGHG3", "IGKC", "IGHA1", "SSR4", "FKBP11",
    ],
    "Tfh": ["CXCL13", "CD200", "BTLA", "ICOS", "PDCD1", "CD40LG", "MAF", "IL21", "TOX2"],
    "T_cell_CD8": ["CD8A", "CD8B", "GZMK", "GZMA", "PRF1", "CD2", "CD3D", "CD3E"],
    "IFNg_Ayers": ["IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"],
}

NAMED_PROBES = {
    "CLDN4": "201428_at",
    "CD8A": "205758_at",
    "TACSTD2": "202286_s_at",
}

# PR 229 HOLDS rule (given): n>=40, partial Spearman rho<0, p<0.05.
HOLDS_N = 40


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "gse4573-cldn4/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as out:
        out.write(r.read())
    return dest


def parse_soft_table(text: str) -> pd.DataFrame:
    lines = text.splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_begin") or l.startswith("!series_matrix_table_begin")) + 1
    end = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_end") or l.startswith("!series_matrix_table_end"))
    return pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", dtype=str)


def parse_series_meta(path: Path) -> tuple[pd.DataFrame, dict]:
    series = {}
    sample_fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Series_"):
                k = line.split("\t", 1)[0][8:]
                v = line.split("\t", 1)[1].strip().strip('"') if "\t" in line else ""
                series.setdefault(k, []).append(v)
            elif line.startswith("!Sample_"):
                key = line.split("\t", 1)[0][8:]
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                sample_fields[key] = vals
    meta = pd.DataFrame(sample_fields)
    if "geo_accession" not in meta.columns:
        raise KeyError(f"no geo_accession in sample fields: {list(meta.columns)}")
    meta.index = meta["geo_accession"].astype(str)
    meta.index.name = "gsm"
    return meta, {k: " | ".join(v) for k, v in series.items()}


def first_symbol(s: str) -> str:
    if not isinstance(s, str) or not s or s == "nan":
        return ""
    return s.split(" /// ")[0].strip()


def collapse_maxmean(probe_expr: pd.DataFrame, id2gene: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = probe_expr.index.map(lambda x: id2gene.get(str(x), ""))
    df = probe_expr.copy()
    df["_gene"] = genes
    df = df[df["_gene"].astype(str).str.len() > 0]
    means = df.drop(columns="_gene").astype(float).mean(axis=1)
    df["_mean"] = means
    df = df.sort_values("_mean", ascending=False)
    kept = df.loc[~df["_gene"].duplicated(keep="first")]
    gene_expr = kept.drop(columns=["_gene", "_mean"]).astype(float)
    gene_expr.index = kept["_gene"].values
    gene_expr = gene_expr.sort_index()
    tmp = df.reset_index()
    probe_col = tmp.columns[0]
    probe_audit = (
        tmp.rename(columns={probe_col: "probe", "_gene": "gene", "_mean": "mean"})
        [["gene", "probe", "mean"]]
        .sort_values(["gene", "mean"], ascending=[True, False])
    )
    return gene_expr, probe_audit


def spearman_ci(x, y, n_boot: int = N_BOOT, seed: int = SEED):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = int(x.size)
    if n < 8:
        return dict(n=n, rho=np.nan, p=np.nan, ci_low=np.nan, ci_high=np.nan)
    rho, p = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        if np.unique(x[idx]).size < 3 or np.unique(y[idx]).size < 3:
            boots[i] = np.nan
            continue
        boots[i] = stats.spearmanr(x[idx], y[idx])[0]
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return dict(n=n, rho=float(rho), p=float(p), ci_low=float(lo), ci_high=float(hi))


def partial_spearman(x, y, covar):
    x = pd.Series(np.asarray(x, float))
    y = pd.Series(np.asarray(y, float))
    c = pd.Series(np.asarray(covar, float), name="cov")
    df = pd.concat([x.rename("x"), y.rename("y"), c], axis=1).dropna()
    n = int(df.shape[0])
    if n < 12:
        return dict(n=n, rho_adj=np.nan, p_adj=np.nan)
    xr = stats.rankdata(df["x"].values)
    yr = stats.rankdata(df["y"].values)
    cr = np.column_stack([np.ones(n), stats.rankdata(df["cov"].values)])
    bx, *_ = np.linalg.lstsq(cr, xr, rcond=None)
    by, *_ = np.linalg.lstsq(cr, yr, rcond=None)
    rx = xr - cr @ bx
    ry = yr - cr @ by
    if np.std(rx) == 0 or np.std(ry) == 0:
        return dict(n=n, rho_adj=np.nan, p_adj=np.nan)
    r = float(np.corrcoef(rx, ry)[0, 1])
    dof = n - 3
    t = r * np.sqrt(dof / max(1e-12, 1 - r ** 2))
    p = float(2 * stats.t.sf(abs(t), dof))
    return dict(n=n, rho_adj=r, p_adj=p)


def zscore_rows(expr: pd.DataFrame) -> pd.DataFrame:
    mu = expr.mean(axis=1)
    sd = expr.std(axis=1, ddof=1).replace(0, np.nan)
    return expr.sub(mu, axis=0).div(sd, axis=0)


def signature_scores(expr: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    z = zscore_rows(expr)
    out, used = {}, {}
    for name, genes in SIGNATURES.items():
        present = [g for g in genes if g in z.index]
        if len(present) < 3:
            continue
        out[name] = z.loc[present].mean(axis=0)
        used[name] = present
    return pd.DataFrame(out), used


def verdict(rho_adj, p_adj, n) -> str:
    if n < HOLDS_N:
        return "UNDERPOWERED"
    if np.isfinite(rho_adj) and rho_adj < 0 and p_adj < 0.05:
        return "HOLDS"
    if np.isfinite(rho_adj) and rho_adj > 0 and p_adj < 0.05:
        return "OPPOSITE"
    return "NO_EVIDENCE"


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    matrix_path = dl(DATA / "GSE4573_series_matrix.txt.gz", GEO_MATRIX)
    annot_path = dl(DATA / "GPL96.annot.gz", GPL_ANNOT)

    with gzip.open(matrix_path, "rt", errors="replace") as fh:
        raw = fh.read()
    meta, series = parse_series_meta(matrix_path)
    lines = raw.split("\n")
    start = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin")) + 1
    end = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end"))
    probe_expr = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", index_col=0)
    probe_expr.index = probe_expr.index.astype(str).str.strip('"')
    probe_expr = probe_expr.astype(float)
    probe_expr = probe_expr.loc[:, [c for c in probe_expr.columns if c in meta.index]]
    meta = meta.loc[probe_expr.columns]

    with gzip.open(annot_path, "rt", errors="replace") as fh:
        annot = parse_soft_table(fh.read())
    annot["symbol"] = annot["Gene symbol"].map(first_symbol)
    id2gene = annot.drop_duplicates("ID").set_index("ID")["symbol"]

    gene_expr, probe_audit = collapse_maxmean(probe_expr, id2gene)

    # Probe confirmation for named targets.
    confirm_rows = []
    for gene, named in NAMED_PROBES.items():
        hits = annot[annot["symbol"] == gene][["ID", "Gene symbol", "Gene title", "Gene ID"]]
        for _, r in hits.iterrows():
            confirm_rows.append({
                "gene": gene,
                "named_probe": named,
                "probe": r["ID"],
                "gpl96_symbol": r["Gene symbol"],
                "entrez": r["Gene ID"],
                "title": r["Gene title"],
                "in_matrix": r["ID"] in probe_expr.index,
                "is_named": r["ID"] == named,
                "chosen_for_collapse": (
                    gene in gene_expr.index
                    and probe_audit.loc[probe_audit["gene"] == gene, "probe"].iloc[0] == r["ID"]
                ),
            })
    # CD8B is the expected missing / dirty U133A case.
    c8b = annot[annot["Gene symbol"].fillna("").str.contains("CD8B", regex=False)]
    for _, r in c8b.iterrows():
        confirm_rows.append({
            "gene": "CD8B",
            "named_probe": "",
            "probe": r["ID"],
            "gpl96_symbol": r["Gene symbol"],
            "entrez": r["Gene ID"],
            "title": r["Gene title"],
            "in_matrix": r["ID"] in probe_expr.index,
            "is_named": False,
            "chosen_for_collapse": False,
        })
    confirm = pd.DataFrame(confirm_rows)
    confirm.to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)

    scores, used = signature_scores(gene_expr)
    epithelial = scores["Epithelial"] if "Epithelial" in scores.columns else None

    # Coverage
    cov_rows = []
    for name, genes in SIGNATURES.items():
        present = [g for g in genes if g in gene_expr.index]
        cov_rows.append({
            "set": name,
            "n_listed": len(genes),
            "n_present": len(present),
            "missing": ",".join([g for g in genes if g not in gene_expr.index]),
            "present": ",".join(present),
        })
    for g in ["CLDN4", "TACSTD2"] + IMMUNE_GENES:
        cov_rows.append({
            "set": f"gene:{g}",
            "n_listed": 1,
            "n_present": int(g in gene_expr.index),
            "missing": "" if g in gene_expr.index else g,
            "present": g if g in gene_expr.index else "",
        })
    coverage = pd.DataFrame(cov_rows)
    coverage.to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)

    # Sample annotation / label inventory (nothing invented).
    sample = meta[["geo_accession", "title", "source_name_ch1", "characteristics_ch1", "description", "data_processing"]].copy()
    sample["cldn4_mas5"] = gene_expr.loc["CLDN4"].reindex(sample.index).values if "CLDN4" in gene_expr.index else np.nan
    sample["cd8a_mas5"] = gene_expr.loc["CD8A"].reindex(sample.index).values if "CD8A" in gene_expr.index else np.nan
    sample["tacstd2_mas5"] = gene_expr.loc["TACSTD2"].reindex(sample.index).values if "TACSTD2" in gene_expr.index else np.nan
    sample.to_csv(TABLES / "sample_annotation.tsv", sep="\t")

    labels = pd.DataFrame([
        {"field": "arrays_in_matrix", "public": "yes", "n": int(probe_expr.shape[1]), "note": "22283 probes x 130 GSM"},
        {"field": "unique_GSM", "public": "yes", "n": int(meta.index.nunique()), "note": "all unique"},
        {"field": "unique_title", "public": "yes", "n": int(meta["title"].nunique()), "note": "LS-* titles; all unique"},
        {"field": "GEO_overall_design_patients", "public": "yes (text only)", "n": 129, "note": "Series says '130 samples from 129 patients'; duplicate pair is not identified in sample metadata"},
        {"field": "histology", "public": "series-level only", "n": 130, "note": "all described as lung squamous carcinoma; no LUAD rows"},
        {"field": "stage", "public": "no", "n": 0, "note": "not a GEO sample characteristic"},
        {"field": "OS / DSS time or event", "public": "no", "n": 0, "note": "prognosis paper; labels not deposited on GEO. No suppl folder."},
        {"field": "ICI / treatment", "public": "no", "n": 0, "note": "resected LUSC atlas, not an ICI series"},
        {"field": "tumor_percent / purity", "public": "no", "n": 0, "note": "epithelial RNA score used as the only public purity proxy"},
        {"field": "CLDN4 finite", "public": "yes", "n": int(sample["cldn4_mas5"].notna().sum()), "note": "named probe 201428_at"},
        {"field": "CD8A finite", "public": "yes", "n": int(sample["cd8a_mas5"].notna().sum()), "note": "named probe 205758_at"},
        {"field": "pairwise CLDN4+CD8A", "public": "yes", "n": int((sample["cldn4_mas5"].notna() & sample["cd8a_mas5"].notna()).sum()), "note": "primary n"},
    ])
    labels.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    n_array = int(probe_expr.shape[1])
    assert n_array == 130
    assert "CLDN4" in gene_expr.index and "CD8A" in gene_expr.index

    cldn4 = gene_expr.loc["CLDN4"]
    gene_rows = []
    partners = ["CD8A"] + [g for g in IMMUNE_GENES if g != "CD8A"]
    if "TACSTD2" in gene_expr.index:
        partners = partners + ["TACSTD2"]
    for g in partners:
        if g not in gene_expr.index:
            gene_rows.append({
                "gene_x": "CLDN4", "gene_y": g, "n": 0, "n_arrays": n_array,
                "present": False, "rho": np.nan, "p": np.nan, "ci_low": np.nan, "ci_high": np.nan,
                "rho_adj_epithelial": np.nan, "p_adj_epithelial": np.nan, "verdict": "ABSENT",
            })
            continue
        y = gene_expr.loc[g]
        crude = spearman_ci(cldn4.values, y.values)
        adj = partial_spearman(cldn4.values, y.values, epithelial.values) if epithelial is not None else dict(n=crude["n"], rho_adj=np.nan, p_adj=np.nan)
        rec = {
            "gene_x": "CLDN4",
            "gene_y": g,
            "n": crude["n"],
            "n_arrays": n_array,
            "present": True,
            "rho": crude["rho"],
            "p": crude["p"],
            "ci_low": crude["ci_low"],
            "ci_high": crude["ci_high"],
            "rho_adj_epithelial": adj["rho_adj"],
            "p_adj_epithelial": adj["p_adj"],
            "verdict": verdict(adj["rho_adj"], adj["p_adj"], adj["n"]),
        }
        gene_rows.append(rec)
    genes_df = pd.DataFrame(gene_rows)
    genes_df["q_bh_immune"] = np.nan
    immune_mask = genes_df["gene_y"].isin(IMMUNE_GENES) & genes_df["present"]
    genes_df.loc[immune_mask, "q_bh_immune"] = bh_fdr(genes_df.loc[immune_mask, "p"].values)
    genes_df.to_csv(TABLES / "spearman_cldn4_vs_genes.tsv", sep="\t", index=False)

    # Companion: TACSTD2 vs CD8A (not the claim; contrast only).
    tac_rows = []
    if "TACSTD2" in gene_expr.index:
        crude = spearman_ci(gene_expr.loc["TACSTD2"].values, gene_expr.loc["CD8A"].values)
        adj = partial_spearman(gene_expr.loc["TACSTD2"].values, gene_expr.loc["CD8A"].values, epithelial.values)
        tac_rows.append({
            "gene_x": "TACSTD2", "gene_y": "CD8A", **crude,
            "rho_adj_epithelial": adj["rho_adj"], "p_adj_epithelial": adj["p_adj"],
            "verdict": verdict(adj["rho_adj"], adj["p_adj"], adj["n"]),
            "note": "companion only; TACSTD2 TLS/CD8 extras are given in PR 229",
        })
    tac_df = pd.DataFrame(tac_rows)
    tac_df.to_csv(TABLES / "spearman_tacstd2_companion.tsv", sep="\t", index=False)

    sig_rows = []
    for sig in SIGNATURES:
        if sig not in scores.columns:
            sig_rows.append({
                "gene": "CLDN4", "signature": sig, "n": n_array,
                "n_genes_in_sig": 0, "genes_used": "",
                "rho": np.nan, "p": np.nan, "ci_low": np.nan, "ci_high": np.nan,
                "rho_adj_epithelial": np.nan, "p_adj_epithelial": np.nan,
                "verdict": "ABSENT",
            })
            continue
        y = scores[sig]
        crude = spearman_ci(cldn4.reindex(y.index).values, y.values)
        if sig == "Epithelial":
            adj = {"rho_adj": np.nan, "p_adj": np.nan, "n": crude["n"]}
            v = "COMPANION"
        else:
            adj = partial_spearman(cldn4.reindex(y.index).values, y.values, epithelial.reindex(y.index).values)
            v = verdict(adj["rho_adj"], adj["p_adj"], adj["n"])
        sig_rows.append({
            "gene": "CLDN4",
            "signature": sig,
            "n": crude["n"],
            "n_genes_in_sig": len(used.get(sig, [])),
            "genes_used": ",".join(used.get(sig, [])),
            "rho": crude["rho"],
            "p": crude["p"],
            "ci_low": crude["ci_low"],
            "ci_high": crude["ci_high"],
            "rho_adj_epithelial": adj["rho_adj"],
            "p_adj_epithelial": adj["p_adj"],
            "verdict": v,
        })
    sig_df = pd.DataFrame(sig_rows)
    immune_sig = sig_df["signature"] != "Epithelial"
    sig_df["q_bh"] = np.nan
    sig_df.loc[immune_sig, "q_bh"] = bh_fdr(sig_df.loc[immune_sig, "p"].values)
    sig_df.to_csv(TABLES / "spearman_cldn4_vs_signatures.tsv", sep="\t", index=False)

    # CLDN4 quartile contrast on CD8A (descriptive).
    q = pd.qcut(cldn4, 4, labels=["Q1", "Q2", "Q3", "Q4"])
    a = gene_expr.loc["CD8A"][q == "Q4"]
    b = gene_expr.loc["CD8A"][q == "Q1"]
    U, p_mw = stats.mannwhitneyu(a.values, b.values, alternative="two-sided")
    rbc = 2 * U / (a.size * b.size) - 1
    hl = pd.DataFrame([{
        "endpoint": "CD8A",
        "n_Q4": int(a.size),
        "n_Q1": int(b.size),
        "median_Q4": float(np.median(a)),
        "median_Q1": float(np.median(b)),
        "U": float(U),
        "p": float(p_mw),
        "rank_biserial_Q4_minus_Q1": float(rbc),
    }])
    hl.to_csv(TABLES / "highlow_cldn4_cd8a.tsv", sep="\t", index=False)

    # Figures
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    x = np.log2(cldn4.values + 1)
    y = np.log2(gene_expr.loc["CD8A"].values + 1)
    ax.scatter(x, y, s=18, alpha=0.7, c="#2c5f8a", edgecolors="none")
    row = genes_df[genes_df.gene_y == "CD8A"].iloc[0]
    ax.set_xlabel("CLDN4  log2(MAS5+1)")
    ax.set_ylabel("CD8A  log2(MAS5+1)")
    ax.set_title(
        f"GSE4573 LUSC  n={int(row.n)} arrays\n"
        f"CLDN4 vs CD8A  ρ={row.rho:.3f}  p={fmt_p(row.p)}"
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a.png", dpi=160)
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a.pdf")
    plt.close(fig)

    plot = genes_df[genes_df.gene_y.isin(IMMUNE_GENES) & genes_df.present].copy()
    plot = plot.sort_values("rho")
    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    y_pos = np.arange(len(plot))
    ax.axvline(0, color="0.5", lw=0.8)
    ax.errorbar(
        plot["rho"], y_pos,
        xerr=[plot["rho"] - plot["ci_low"], plot["ci_high"] - plot["rho"]],
        fmt="o", color="#2c5f8a", ecolor="#2c5f8a", capsize=2, ms=5,
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot["gene_y"])
    ax.set_xlabel("Spearman ρ vs CLDN4 (bootstrap 95% CI)")
    ax.set_title(f"GSE4573 LUSC  n={n_array} arrays  CLDN4 vs immune genes")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cldn4_immune_forest.png", dpi=160)
    fig.savefig(FIGURES / "fig2_cldn4_immune_forest.pdf")
    plt.close(fig)

    if "T_cell_CD8" in scores.columns:
        fig, ax = plt.subplots(figsize=(5.2, 4.4))
        ax.scatter(np.log2(cldn4.reindex(scores.index).values + 1), scores["T_cell_CD8"].values,
                   s=18, alpha=0.7, c="#6b3a2a", edgecolors="none")
        srow = sig_df[sig_df.signature == "T_cell_CD8"].iloc[0]
        ax.set_xlabel("CLDN4  log2(MAS5+1)")
        ax.set_ylabel("T_cell_CD8  mean-z")
        ax.set_title(
            f"GSE4573 LUSC  n={int(srow.n)}\n"
            f"CLDN4 vs CD8 signature  ρ={srow.rho:.3f}  p={fmt_p(srow.p)}"
        )
        fig.tight_layout()
        fig.savefig(FIGURES / "fig3_cldn4_vs_cd8sig.png", dpi=160)
        fig.savefig(FIGURES / "fig3_cldn4_vs_cd8sig.pdf")
        plt.close(fig)

    c8 = genes_df[genes_df.gene_y == "CD8A"].iloc[0]
    epi = sig_df[sig_df.signature == "Epithelial"]
    epi_rec = epi.iloc[0].to_dict() if len(epi) else {}
    summary = {
        "dataset": "GSE4573",
        "pmid": "16885343",
        "platform": "GPL96 Affymetrix HG-U133A",
        "histology": "LUSC (series-level; every sample described as lung squamous carcinoma)",
        "n_arrays": n_array,
        "n_probes": int(probe_expr.shape[0]),
        "n_genes_collapsed": int(gene_expr.shape[0]),
        "n_patients_geo_text": 129,
        "n_primary_cldn4_cd8a": int(c8["n"]),
        "processing": "MAS5 as deposited (linear signal). Spearman is rank-based.",
        "cldn4_probe": "201428_at",
        "cd8a_probe": "205758_at",
        "cd8b": "not a clean U133A symbol (LOC100996919///CD8B); excluded by first-symbol collapse",
        "cldn4_vs_cd8a": {
            "n": int(c8["n"]),
            "rho": c8["rho"],
            "p": c8["p"],
            "ci": [c8["ci_low"], c8["ci_high"]],
            "rho_adj_epithelial": c8["rho_adj_epithelial"],
            "p_adj_epithelial": c8["p_adj_epithelial"],
            "verdict": c8["verdict"],
        },
        "cldn4_vs_epithelial": {
            "n": int(epi_rec.get("n", 0) or 0),
            "rho": epi_rec.get("rho"),
            "p": epi_rec.get("p"),
            "ci": [epi_rec.get("ci_low"), epi_rec.get("ci_high")],
        },
        "pr229_tacstd2_given": {
            "note": "TLS/CD8 TACSTD2 extras from PR 229; not re-claimed here",
            "n": 130,
            "T_cell_CD8_rho": -0.24674033572502613,
            "T_cell_CD8_p": 0.004654889763503568,
            "TLS_Cabrita_verdict": "NO_EVIDENCE",
            "T_cell_CD8_verdict": "HOLDS",
        },
        "missing_public_labels": ["OS", "stage", "ICI", "tumor_percent", "patient-duplicate ID"],
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print("n_arrays", n_array, "genes", gene_expr.shape[0])
    print(genes_df[genes_df.gene_y.isin(["CD8A", "IFNG", "CXCL13", "MS4A1", "TACSTD2"])].to_string(index=False))
    print(sig_df[["signature", "n", "n_genes_in_sig", "rho", "p", "rho_adj_epithelial", "p_adj_epithelial", "verdict"]].to_string(index=False))
    print("[done]", TABLES)


if __name__ == "__main__":
    main()
