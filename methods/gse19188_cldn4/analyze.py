#!/usr/bin/env python3
"""GSE19188 NSCLC Affy U133 Plus 2.0: CLDN4 vs CD8A / CD274.

Hou et al., PLoS One 2010, PMID 20421987. Public series matrix is RMA
then log2(ratio to geometric mean). Spearman is rank-based.

Primary n is tumor arrays with finite CLDN4 + partner. Adjacent-normal
arrays are inventoried and not mixed into the NSCLC pairwise tests.

Downloads stay under $GSE19188_CLDN4_DATA (default /tmp/gse19188_cldn4)
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
DATA = Path(os.environ.get("GSE19188_CLDN4_DATA", "/tmp/gse19188_cldn4"))
SEED = 20260817
N_BOOT = 2000

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19188/"
    "matrix/GSE19188_series_matrix.txt.gz"
)
GPL_SOFT = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    "?acc=GPL570&targ=self&form=text&view=data"
)

# Primary partners requested for this slice.
PRIMARY_PARTNERS = ["CD8A", "CD274"]

IMMUNE_GENES = [
    "CD8A",
    "CD8B",
    "CD274",
    "PDCD1",
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
    "MS4A1",
]

SIGNATURES = {
    "T_cell_CD8": ["CD8A", "CD8B", "GZMK", "GZMA", "PRF1", "CD2", "CD3D", "CD3E"],
    "IFNg_Ayers": ["IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"],
}

# Named Plus2 probes (confirmed against GPL570 Gene Symbol).
NAMED_PROBES = {
    "CLDN4": "201428_at",
    "CD8A": "205758_at",
    "CD274": "223834_at",
    "TACSTD2": "202286_s_at",
}

HOLDS_N = 40


def dl(dest: Path, url: str, ua: str = "gse19188-cldn4/1.0") -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": ua})
    with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as out:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    return dest


def parse_soft_table(text: str) -> pd.DataFrame:
    lines = text.splitlines()
    start = next(
        i
        for i, l in enumerate(lines)
        if l.startswith("!platform_table_begin") or l.startswith("!series_matrix_table_begin")
    ) + 1
    end = next(
        i
        for i, l in enumerate(lines)
        if l.startswith("!platform_table_end") or l.startswith("!series_matrix_table_end")
    )
    return pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", dtype=str)


def parse_series_meta(path: Path) -> tuple[pd.DataFrame, dict]:
    series: dict[str, list[str]] = {}
    sample_fields: dict[str, list[list[str]]] = {}
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
                sample_fields.setdefault(key, []).append(vals)
    n = len(sample_fields["geo_accession"][0])
    rows = []
    for i in range(n):
        rec = {
            "gsm": sample_fields["geo_accession"][0][i],
            "title": sample_fields["title"][0][i],
            "description": sample_fields.get("description", [[""] * n])[0][i],
        }
        for block in sample_fields.get("characteristics_ch1", []):
            v = block[i]
            if ": " in v:
                k, val = v.split(": ", 1)
                rec[k] = val
        rows.append(rec)
    meta = pd.DataFrame(rows).set_index("gsm")
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


def fmt_rho(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:+.3f}"


def pair_row(cldn4, y, epithelial, gene_y: str, subset: str, n_arrays: int) -> dict:
    crude = spearman_ci(cldn4.values, y.values)
    if epithelial is not None:
        adj = partial_spearman(cldn4.values, y.values, epithelial.values)
    else:
        adj = dict(n=crude["n"], rho_adj=np.nan, p_adj=np.nan)
    return {
        "subset": subset,
        "gene_x": "CLDN4",
        "gene_y": gene_y,
        "n": crude["n"],
        "n_arrays": n_arrays,
        "present": True,
        "rho": crude["rho"],
        "p": crude["p"],
        "ci_low": crude["ci_low"],
        "ci_high": crude["ci_high"],
        "rho_adj_epithelial": adj["rho_adj"],
        "p_adj_epithelial": adj["p_adj"],
        "verdict": verdict(adj["rho_adj"], adj["p_adj"], adj["n"]),
    }


def highlow(cldn4: pd.Series, y: pd.Series, endpoint: str) -> dict:
    q = pd.qcut(cldn4, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    a = y[q == "Q4"]
    b = y[q == "Q1"]
    U, p_mw = stats.mannwhitneyu(a.values, b.values, alternative="two-sided")
    rbc = 2 * U / (a.size * b.size) - 1
    return {
        "endpoint": endpoint,
        "n_Q4": int(a.size),
        "n_Q1": int(b.size),
        "median_Q4": float(np.median(a)),
        "median_Q1": float(np.median(b)),
        "U": float(U),
        "p": float(p_mw),
        "rank_biserial_Q4_minus_Q1": float(rbc),
    }


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    matrix_path = dl(DATA / "GSE19188_series_matrix.txt.gz", GEO_MATRIX)
    annot_path = dl(DATA / "GPL570.soft.txt", GPL_SOFT, ua="Mozilla/5.0 gse19188-cldn4/1.0")

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
    meta = meta.loc[probe_expr.columns].copy()

    annot = parse_soft_table(annot_path.read_text(errors="replace"))
    annot["symbol"] = annot["Gene Symbol"].map(first_symbol)
    id2gene = annot.drop_duplicates("ID").set_index("ID")["symbol"]

    gene_expr, probe_audit = collapse_maxmean(probe_expr, id2gene)

    # Confirm named probes.
    probe_rows = []
    for gene, pid in NAMED_PROBES.items():
        if pid not in probe_expr.index:
            probe_rows.append({"gene": gene, "named_probe": pid, "present": False, "gpl570_symbol": "", "used": "ABSENT"})
            continue
        sym = str(id2gene.get(pid, ""))
        used = "named" if gene in gene_expr.index else "named_not_collapsed"
        probe_rows.append({
            "gene": gene,
            "named_probe": pid,
            "present": True,
            "gpl570_symbol": sym,
            "used": used,
            "n_probes_for_symbol": int((probe_audit["gene"] == gene).sum()) if gene in set(probe_audit["gene"]) else 0,
            "collapsed_probe": (
                probe_audit.loc[probe_audit["gene"] == gene, "probe"].iloc[0]
                if gene in set(probe_audit["gene"]) else ""
            ),
        })
    # Extra CD274 probe on Plus2.
    extra_cd274 = "227458_at"
    if extra_cd274 in probe_expr.index:
        probe_rows.append({
            "gene": "CD274",
            "named_probe": extra_cd274,
            "present": True,
            "gpl570_symbol": str(id2gene.get(extra_cd274, "")),
            "used": "companion_probe",
            "n_probes_for_symbol": int((probe_audit["gene"] == "CD274").sum()),
            "collapsed_probe": probe_audit.loc[probe_audit["gene"] == "CD274", "probe"].iloc[0],
        })
    pd.DataFrame(probe_rows).to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)

    # Force named probes for the four headline genes (same rule as GSE4573).
    for gene, pid in NAMED_PROBES.items():
        if pid in probe_expr.index:
            gene_expr.loc[gene] = probe_expr.loc[pid].astype(float)

    meta["tissue_type"] = meta["tissue type"].astype(str).str.strip().str.lower()
    meta["histology"] = meta["cell type"].astype(str).str.strip()
    meta["is_tumor"] = meta["tissue_type"].eq("tumor")
    meta["is_healthy"] = meta["tissue_type"].eq("healthy")
    os_raw = meta["overall survival"].astype(str)
    meta["os_months"] = pd.to_numeric(os_raw.replace({"Not available": np.nan}), errors="coerce")
    meta["os_status"] = meta["status"].astype(str)
    meta["os_event"] = meta["os_status"].map({"deceased": 1, "alive": 0})

    tumor = meta.index[meta["is_tumor"]]
    healthy = meta.index[meta["is_healthy"]]
    n_array = int(probe_expr.shape[1])
    n_tumor = int(len(tumor))
    n_healthy = int(len(healthy))
    assert n_array == 156
    assert n_tumor == 91
    assert n_healthy == 65
    for g in ["CLDN4", "CD8A", "CD274"]:
        if g not in gene_expr.index:
            raise SystemExit(f"required gene missing after collapse: {g}")

    tumor_expr = gene_expr.loc[:, tumor]
    scores, used = signature_scores(tumor_expr)
    epithelial = scores["Epithelial"] if "Epithelial" in scores.columns else None

    sample = meta.copy()
    sample["cldn4"] = gene_expr.loc["CLDN4"]
    sample["cd8a"] = gene_expr.loc["CD8A"]
    sample["cd274"] = gene_expr.loc["CD274"]
    if "TACSTD2" in gene_expr.index:
        sample["tacstd2"] = gene_expr.loc["TACSTD2"]
    if epithelial is not None:
        sample["epithelial_mean_z"] = np.nan
        sample.loc[epithelial.index, "epithelial_mean_z"] = epithelial
    sample.to_csv(TABLES / "sample_annotation.tsv", sep="\t")

    hist_counts = (
        meta.loc[tumor, "histology"].value_counts().rename_axis("histology").reset_index(name="n")
    )
    hist_counts.to_csv(TABLES / "histology_counts.tsv", sep="\t", index=False)

    n_os_tumor = int(meta.loc[tumor, "os_months"].notna().sum())
    n_os_event = int(meta.loc[tumor, "os_event"].notna().sum())
    labels = pd.DataFrame([
        {"field": "arrays in series matrix", "public": "yes", "n": n_array, "note": "54,675 probes; RMA then log2(ratio to geometric mean)"},
        {"field": "unique GSM", "public": "yes", "n": int(meta.index.nunique()), "note": "all unique"},
        {"field": "patients (GEO design text)", "public": "text only", "n": 91, "note": "overall_design: 91 tumor + 65 adjacent normal from a cohort of 91 patients"},
        {"field": "tumor arrays (tissue type: tumor)", "public": "yes", "n": n_tumor, "note": "PRIMARY n; ADC 45 / SCC 27 / LCC 19"},
        {"field": "adjacent-normal arrays (tissue type: healthy)", "public": "yes", "n": n_healthy, "note": "not mixed into the NSCLC pairwise tests"},
        {"field": "title suffix T/N vs tissue type", "public": "yes", "n": n_array, "note": "title is not reliable (e.g. GSM475659 title 2344T is healthy); use tissue type"},
        {"field": "histology (cell type) on tumors", "public": "yes", "n": n_tumor, "note": "ADC/SCC/LCC; no mixed/other tumor rows"},
        {"field": "OS months on tumors", "public": "partial", "n": n_os_tumor, "note": "numeric overall survival; 9 tumors Not available"},
        {"field": "OS event on tumors", "public": "partial", "n": n_os_event, "note": "alive/deceased; not this claim"},
        {"field": "stage", "public": "no", "n": 0, "note": "paper is early-stage NSCLC; stage is not a GEO characteristic"},
        {"field": "ICI / treatment", "public": "no", "n": 0, "note": "resected atlas, not an ICI series"},
        {"field": "tumor % / ABSOLUTE purity", "public": "no", "n": 0, "note": "only public proxy is an RNA epithelial score"},
        {"field": "CLDN4 finite on tumors (201428_at)", "public": "yes", "n": int(sample.loc[tumor, "cldn4"].notna().sum()), "note": "named Plus2 probe; GPL570 = CLDN4 / Entrez 1364"},
        {"field": "CD8A finite on tumors (205758_at)", "public": "yes", "n": int(sample.loc[tumor, "cd8a"].notna().sum()), "note": "named Plus2 probe; GPL570 = CD8A / Entrez 925"},
        {"field": "CD274 finite on tumors (223834_at)", "public": "yes", "n": int(sample.loc[tumor, "cd274"].notna().sum()), "note": "named Plus2 probe; GPL570 = CD274 / Entrez 29126"},
        {"field": "pairwise CLDN4+CD8A tumors", "public": "yes", "n": int((sample.loc[tumor, "cldn4"].notna() & sample.loc[tumor, "cd8a"].notna()).sum()), "note": "primary n for CLDN4 vs CD8A"},
        {"field": "pairwise CLDN4+CD274 tumors", "public": "yes", "n": int((sample.loc[tumor, "cldn4"].notna() & sample.loc[tumor, "cd274"].notna()).sum()), "note": "primary n for CLDN4 vs CD274"},
    ])
    labels.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    cldn4 = tumor_expr.loc["CLDN4"]
    gene_rows = []
    partners = PRIMARY_PARTNERS + [g for g in IMMUNE_GENES if g not in PRIMARY_PARTNERS]
    if "TACSTD2" in tumor_expr.index:
        partners = partners + ["TACSTD2"]
    for g in partners:
        if g not in tumor_expr.index:
            gene_rows.append({
                "subset": "tumor",
                "gene_x": "CLDN4",
                "gene_y": g,
                "n": 0,
                "n_arrays": n_tumor,
                "present": False,
                "rho": np.nan,
                "p": np.nan,
                "ci_low": np.nan,
                "ci_high": np.nan,
                "rho_adj_epithelial": np.nan,
                "p_adj_epithelial": np.nan,
                "verdict": "ABSENT",
            })
            continue
        gene_rows.append(pair_row(cldn4, tumor_expr.loc[g], epithelial, g, "tumor", n_tumor))
    genes_df = pd.DataFrame(gene_rows)
    genes_df["q_bh_immune"] = np.nan
    immune_mask = genes_df["gene_y"].isin(IMMUNE_GENES) & genes_df["present"]
    genes_df.loc[immune_mask, "q_bh_immune"] = bh_fdr(genes_df.loc[immune_mask, "p"].values)
    genes_df.to_csv(TABLES / "spearman_cldn4_vs_genes.tsv", sep="\t", index=False)

    # Histology-stratified primary pairs (honest n; SCC/LCC are <40).
    histo_rows = []
    for hist in ["ADC", "SCC", "LCC"]:
        idx = meta.index[meta["is_tumor"] & meta["histology"].eq(hist)]
        if len(idx) < 8:
            continue
        sub_cldn4 = gene_expr.loc["CLDN4", idx]
        sub_epi = epithelial.reindex(idx) if epithelial is not None else None
        for g in PRIMARY_PARTNERS:
            histo_rows.append(pair_row(sub_cldn4, gene_expr.loc[g, idx], sub_epi, g, f"tumor_{hist}", int(len(idx))))
    histo_df = pd.DataFrame(histo_rows)
    histo_df.to_csv(TABLES / "spearman_cldn4_by_histology.tsv", sep="\t", index=False)

    # Companion: TACSTD2 vs CD8A / CD274 on tumors (not the claim).
    tac_rows = []
    if "TACSTD2" in tumor_expr.index:
        for g in PRIMARY_PARTNERS:
            crude = spearman_ci(tumor_expr.loc["TACSTD2"].values, tumor_expr.loc[g].values)
            adj = partial_spearman(tumor_expr.loc["TACSTD2"].values, tumor_expr.loc[g].values, epithelial.values)
            tac_rows.append({
                "gene_x": "TACSTD2",
                "gene_y": g,
                **crude,
                "rho_adj_epithelial": adj["rho_adj"],
                "p_adj_epithelial": adj["p_adj"],
                "verdict": verdict(adj["rho_adj"], adj["p_adj"], adj["n"]),
                "note": "companion only; this slice is CLDN4",
            })
    tac_df = pd.DataFrame(tac_rows)
    tac_df.to_csv(TABLES / "spearman_tacstd2_companion.tsv", sep="\t", index=False)

    # Alternate CD274 probe 227458_at (not used for the headline).
    alt_rows = []
    if extra_cd274 in probe_expr.index:
        y = probe_expr.loc[extra_cd274, tumor]
        crude = spearman_ci(cldn4.values, y.values)
        adj = partial_spearman(cldn4.values, y.values, epithelial.values)
        alt_rows.append({
            "gene_x": "CLDN4",
            "gene_y": "CD274",
            "probe_y": extra_cd274,
            **crude,
            "rho_adj_epithelial": adj["rho_adj"],
            "p_adj_epithelial": adj["p_adj"],
            "verdict": verdict(adj["rho_adj"], adj["p_adj"], adj["n"]),
            "note": "companion probe; headline uses named 223834_at",
        })
    pd.DataFrame(alt_rows).to_csv(TABLES / "spearman_cd274_alt_probe.tsv", sep="\t", index=False)

    sig_rows = []
    for sig in SIGNATURES:
        if sig not in scores.columns:
            sig_rows.append({
                "gene": "CLDN4",
                "signature": sig,
                "n": n_tumor,
                "n_genes_in_sig": 0,
                "genes_used": "",
                "rho": np.nan,
                "p": np.nan,
                "ci_low": np.nan,
                "ci_high": np.nan,
                "rho_adj_epithelial": np.nan,
                "p_adj_epithelial": np.nan,
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
    sig_df.to_csv(TABLES / "spearman_cldn4_vs_signatures.tsv", sep="\t", index=False)

    hl = pd.DataFrame([
        highlow(cldn4, tumor_expr.loc["CD8A"], "CD8A"),
        highlow(cldn4, tumor_expr.loc["CD274"], "CD274"),
    ])
    hl.to_csv(TABLES / "highlow_cldn4.tsv", sep="\t", index=False)

    # Figures
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.3))
    for ax, gene, color in zip(axes, PRIMARY_PARTNERS, ["#2c5f8a", "#6b3a2a"]):
        row = genes_df[genes_df.gene_y == gene].iloc[0]
        ax.scatter(cldn4.values, tumor_expr.loc[gene].values, s=18, alpha=0.75, c=color, edgecolors="none")
        ax.set_xlabel("CLDN4  201428_at  (deposited log2-ratio)")
        ax.set_ylabel(f"{gene}  (deposited log2-ratio)")
        ax.set_title(
            f"GSE19188 NSCLC tumor  n={int(row.n)}\n"
            f"CLDN4 vs {gene}  ρ={row.rho:.3f}  p={fmt_p(row.p)}"
        )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_cd274.png", dpi=160)
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_cd274.pdf")
    plt.close(fig)

    plot = genes_df[genes_df.gene_y.isin(IMMUNE_GENES) & genes_df.present].copy()
    plot = plot.sort_values("rho")
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    y_pos = np.arange(len(plot))
    ax.axvline(0, color="0.5", lw=0.8)
    ax.errorbar(
        plot["rho"],
        y_pos,
        xerr=[plot["rho"] - plot["ci_low"], plot["ci_high"] - plot["rho"]],
        fmt="o",
        color="#2c5f8a",
        ecolor="#2c5f8a",
        capsize=2,
        ms=5,
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot["gene_y"])
    ax.set_xlabel("Spearman ρ vs CLDN4 (bootstrap 95% CI)")
    ax.set_title(f"GSE19188 NSCLC tumor  n={n_tumor}  CLDN4 vs immune genes")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cldn4_immune_forest.png", dpi=160)
    fig.savefig(FIGURES / "fig2_cldn4_immune_forest.pdf")
    plt.close(fig)

    if "T_cell_CD8" in scores.columns:
        fig, ax = plt.subplots(figsize=(5.2, 4.4))
        ax.scatter(
            cldn4.reindex(scores.index).values,
            scores["T_cell_CD8"].values,
            s=18,
            alpha=0.7,
            c="#6b3a2a",
            edgecolors="none",
        )
        srow = sig_df[sig_df.signature == "T_cell_CD8"].iloc[0]
        ax.set_xlabel("CLDN4  (deposited log2-ratio)")
        ax.set_ylabel("T_cell_CD8  mean-z")
        ax.set_title(
            f"GSE19188 NSCLC tumor  n={int(srow.n)}\n"
            f"CLDN4 vs CD8 signature  ρ={srow.rho:.3f}  p={fmt_p(srow.p)}"
        )
        fig.tight_layout()
        fig.savefig(FIGURES / "fig3_cldn4_vs_cd8sig.png", dpi=160)
        fig.savefig(FIGURES / "fig3_cldn4_vs_cd8sig.pdf")
        plt.close(fig)

    c8 = genes_df[genes_df.gene_y == "CD8A"].iloc[0]
    pdl1 = genes_df[genes_df.gene_y == "CD274"].iloc[0]
    epi = sig_df[sig_df.signature == "Epithelial"]
    epi_rec = epi.iloc[0].to_dict() if len(epi) else {}
    summary = {
        "dataset": "GSE19188",
        "pmid": "20421987",
        "platform": "GPL570 Affymetrix HG-U133 Plus 2.0",
        "processing": "RMA; probe-set intensities <30 reset to 30; log2(ratio to geometric mean). Spearman is rank-based.",
        "histology": "NSCLC tumors: ADC 45, SCC 27, LCC 19; plus 65 adjacent-normal arrays not used in the primary pairs",
        "n_arrays": n_array,
        "n_tumor": n_tumor,
        "n_healthy": n_healthy,
        "n_patients_geo_text": 91,
        "n_probes": int(probe_expr.shape[0]),
        "n_genes_collapsed": int(gene_expr.shape[0]),
        "n_primary_cldn4_cd8a": int(c8["n"]),
        "n_primary_cldn4_cd274": int(pdl1["n"]),
        "cldn4_probe": NAMED_PROBES["CLDN4"],
        "cd8a_probe": NAMED_PROBES["CD8A"],
        "cd274_probe": NAMED_PROBES["CD274"],
        "cd274_alt_probe": extra_cd274,
        "cldn4_vs_cd8a": {
            "n": int(c8["n"]),
            "rho": c8["rho"],
            "p": c8["p"],
            "ci": [c8["ci_low"], c8["ci_high"]],
            "rho_adj_epithelial": c8["rho_adj_epithelial"],
            "p_adj_epithelial": c8["p_adj_epithelial"],
            "verdict": c8["verdict"],
        },
        "cldn4_vs_cd274": {
            "n": int(pdl1["n"]),
            "rho": pdl1["rho"],
            "p": pdl1["p"],
            "ci": [pdl1["ci_low"], pdl1["ci_high"]],
            "rho_adj_epithelial": pdl1["rho_adj_epithelial"],
            "p_adj_epithelial": pdl1["p_adj_epithelial"],
            "verdict": pdl1["verdict"],
        },
        "cldn4_vs_epithelial": {
            "n": int(epi_rec.get("n", 0) or 0),
            "rho": epi_rec.get("rho"),
            "p": epi_rec.get("p"),
            "ci": [epi_rec.get("ci_low"), epi_rec.get("ci_high")],
        },
        "missing_public_labels": ["stage", "ICI", "tumor_percent"],
        "os_on_geo_tumors": {"n_time": n_os_tumor, "n_event": n_os_event, "used_in_this_claim": False},
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print("n_arrays", n_array, "tumor", n_tumor, "healthy", n_healthy, "genes", gene_expr.shape[0])
    print(genes_df[genes_df.gene_y.isin(["CD8A", "CD274", "TACSTD2"])].to_string(index=False))
    print(histo_df.to_string(index=False))
    print(sig_df[["signature", "n", "n_genes_in_sig", "rho", "p", "rho_adj_epithelial", "p_adj_epithelial", "verdict"]].to_string(index=False))
    print("[done]", TABLES)


if __name__ == "__main__":
    main()
