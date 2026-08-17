#!/usr/bin/env python3
"""GSE10245 NSCLC Affy U133 Plus 2.0: CLDN4 vs CD8A / CD274.

Additive public series-matrix slice. Honest n from the deposited matrix
(58 arrays: 40 adenocarcinoma + 18 squamous). Lung Cancer Explorer's n=48
is a different reprocessed subset and is not used.

Downloads stay under $GSE10245_CLDN4_DATA (default /tmp/gse10245_cldn4)
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
DATA = Path(os.environ.get("GSE10245_CLDN4_DATA", "/tmp/gse10245_cldn4"))
SEED = 20260817
N_BOOT = 2000
HOLDS_N = 40

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10245/"
    "matrix/GSE10245_series_matrix.txt.gz"
)
GPL_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"

# Named probes (official GEO GPL570 annot). Collapse still uses max-mean.
NAMED_PROBES = {
    "CLDN4": "201428_at",  # 1569421_at is background-level on this matrix
    "CD8A": "205758_at",
    "CD274": "227458_at",  # higher-mean of two clean CD274 probes
    "TACSTD2": "202286_s_at",
}
CD274_PROBES = ["223834_at", "227458_at"]
CLDN4_PROBES = ["201428_at", "1569421_at"]

EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]
PRIMARY = ["CD8A", "CD274"]


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "gse10245-cldn4/1.0"})
    with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as out:
        out.write(r.read())
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


def residualize(y, Z):
    x = np.column_stack([np.ones(len(y)), Z])
    b, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ b


def partial_spearman(x, y, *covars):
    cols = [np.asarray(x, float), np.asarray(y, float)]
    for c in covars:
        cols.append(np.asarray(c, float))
    arr = np.column_stack(cols)
    ok = np.isfinite(arr).all(axis=1)
    arr = arr[ok]
    n = int(arr.shape[0])
    k = arr.shape[1] - 2
    if n < 8 + k:
        return dict(n=n, rho_adj=np.nan, p_adj=np.nan)
    xr = stats.rankdata(arr[:, 0])
    yr = stats.rankdata(arr[:, 1])
    Z = np.column_stack([stats.rankdata(arr[:, 2 + i]) for i in range(k)])
    rx = residualize(xr, Z)
    ry = residualize(yr, Z)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return dict(n=n, rho_adj=np.nan, p_adj=np.nan)
    r, _ = stats.pearsonr(rx, ry)
    r = float(np.clip(r, -0.999999, 0.999999))
    dof = n - 2 - k
    t = r * np.sqrt(dof / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), dof))
    return dict(n=n, rho_adj=r, p_adj=p)


def verdict(rho_adj, p_adj, n) -> str:
    if n < HOLDS_N:
        return "UNDERPOWERED"
    if np.isfinite(rho_adj) and rho_adj < 0 and p_adj < 0.05:
        return "HOLDS"
    if np.isfinite(rho_adj) and rho_adj > 0 and p_adj < 0.05:
        return "OPPOSITE"
    return "NO_EVIDENCE"


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_rho(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:+.3f}"


def fmt_ci(lo, hi) -> str:
    if lo is None or hi is None or not (np.isfinite(lo) and np.isfinite(hi)):
        return "NA"
    return f"{lo:+.3f} to {hi:+.3f}"


def zmean(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns), present
    sub = expr.loc[present].astype(float)
    mu = sub.mean(axis=1)
    sd = sub.std(axis=1, ddof=1).replace(0, np.nan)
    z = sub.sub(mu, axis=0).div(sd, axis=0)
    return z.mean(axis=0), present


def histolabel(s: str) -> str:
    s = (s or "").lower()
    if "adenocarcinoma" in s:
        return "ADC"
    if "squamous" in s:
        return "SCC"
    return "other"


def pair_row(cldn4, y, epithelial, histo_bin, subset, partner, source):
    crude = spearman_ci(cldn4.values, y.values)
    adj_epi = partial_spearman(cldn4.values, y.values, epithelial.values)
    rec = {
        "subset": subset,
        "partner": partner,
        "source": source,
        "n": crude["n"],
        "rho": crude["rho"],
        "p": crude["p"],
        "ci_low": crude["ci_low"],
        "ci_high": crude["ci_high"],
        "rho_adj_epithelial": adj_epi["rho_adj"],
        "p_adj_epithelial": adj_epi["p_adj"],
        "rho_adj_histology": np.nan,
        "p_adj_histology": np.nan,
        "verdict_epithelial": verdict(adj_epi["rho_adj"], adj_epi["p_adj"], adj_epi["n"]),
        "verdict_histology": "NA",
    }
    if histo_bin is not None and np.unique(histo_bin).size > 1:
        adj_h = partial_spearman(cldn4.values, y.values, histo_bin)
        rec["rho_adj_histology"] = adj_h["rho_adj"]
        rec["p_adj_histology"] = adj_h["p_adj"]
        rec["verdict_histology"] = verdict(adj_h["rho_adj"], adj_h["p_adj"], adj_h["n"])
    return rec


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    matrix_path = dl(DATA / "GSE10245_series_matrix.txt.gz", GEO_MATRIX)
    annot_path = dl(DATA / "GPL570.annot.gz", GPL_ANNOT)

    with gzip.open(matrix_path, "rt", errors="replace") as fh:
        raw = fh.read()
    meta, series = parse_series_meta(matrix_path)
    lines = raw.split("\n")
    start = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin")) + 1
    end = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end"))
    probe_expr = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", index_col=0)
    probe_expr.index = probe_expr.index.astype(str).str.strip('"')
    probe_expr.columns = [c.strip('"') for c in probe_expr.columns]
    probe_expr = probe_expr.astype(float)
    probe_expr = probe_expr.loc[:, [c for c in probe_expr.columns if c in meta.index]]
    meta = meta.loc[probe_expr.columns].copy()
    meta["histology"] = meta["characteristics_ch1"].map(histolabel)
    if (meta["histology"] == "other").any():
        raise ValueError(f"unexpected histology labels: {meta['histology'].value_counts().to_dict()}")

    with gzip.open(annot_path, "rt", errors="replace") as fh:
        annot = parse_soft_table(fh.read())
    annot["symbol"] = annot["Gene symbol"].map(first_symbol)
    id2gene = annot.drop_duplicates("ID").set_index("ID")["symbol"]

    gene_expr, probe_audit = collapse_maxmean(probe_expr, id2gene)
    epithelial, epi_used = zmean(gene_expr, EPITHELIAL)

    n_array = int(probe_expr.shape[1])
    n_adc = int((meta["histology"] == "ADC").sum())
    n_scc = int((meta["histology"] == "SCC").sum())
    n_nan = int(probe_expr.isna().sum().sum())
    assert n_array == 58
    assert n_adc == 40 and n_scc == 18
    assert n_nan == 0
    for g in ["CLDN4", "CD8A", "CD274"]:
        if g not in gene_expr.index:
            raise KeyError(f"{g} missing after collapse")

    # Probe confirmation
    confirm_rows = []
    wanted = {
        "CLDN4": CLDN4_PROBES,
        "CD8A": ["205758_at"],
        "CD274": CD274_PROBES,
        "TACSTD2": ["202285_s_at", "202286_s_at", "202287_s_at", "227128_s_at"],
    }
    for gene, named_list in wanted.items():
        hits = annot[annot["symbol"] == gene][["ID", "Gene symbol", "Gene title", "Gene ID"]]
        chosen = ""
        if gene in gene_expr.index:
            chosen = str(probe_audit.loc[probe_audit["gene"] == gene, "probe"].iloc[0])
        for _, r in hits.iterrows():
            confirm_rows.append({
                "gene": gene,
                "probe": r["ID"],
                "gpl570_symbol": r["Gene symbol"],
                "entrez": r["Gene ID"],
                "title": r["Gene title"],
                "in_matrix": r["ID"] in probe_expr.index,
                "is_named_primary": r["ID"] == NAMED_PROBES.get(gene, ""),
                "listed_named": r["ID"] in named_list,
                "chosen_for_collapse": r["ID"] == chosen,
                "mean_gcrma": float(probe_expr.loc[r["ID"]].mean()) if r["ID"] in probe_expr.index else np.nan,
            })
    confirm = pd.DataFrame(confirm_rows)
    confirm.to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)

    # Sample table
    sample = meta[["geo_accession", "title", "source_name_ch1", "characteristics_ch1", "histology"]].copy()
    sample["cldn4"] = gene_expr.loc["CLDN4"].reindex(sample.index).values
    sample["cd8a"] = gene_expr.loc["CD8A"].reindex(sample.index).values
    sample["cd274"] = gene_expr.loc["CD274"].reindex(sample.index).values
    sample["tacstd2"] = gene_expr.loc["TACSTD2"].reindex(sample.index).values if "TACSTD2" in gene_expr.index else np.nan
    sample["epithelial_z"] = epithelial.reindex(sample.index).values
    sample["cldn4_201428_at"] = probe_expr.loc["201428_at"].reindex(sample.index).values
    sample["cd274_223834_at"] = probe_expr.loc["223834_at"].reindex(sample.index).values
    sample["cd274_227458_at"] = probe_expr.loc["227458_at"].reindex(sample.index).values
    sample.to_csv(TABLES / "sample_annotation.tsv", sep="\t")

    # Honest n inventory — only what is on GEO
    titles_unique = int(meta["title"].nunique())
    gsm_unique = int(meta.index.nunique())
    labels = pd.DataFrame([
        {"item": "arrays_in_series_matrix", "public": "yes", "n": n_array,
         "note": "54675 probes × 58 GSM; 0 missing gcRMA values"},
        {"item": "unique_GSM", "public": "yes", "n": gsm_unique, "note": "all unique"},
        {"item": "unique_title", "public": "yes", "n": titles_unique,
         "note": "NSCLC_AC_* / NSCLC_SCC_*; all unique"},
        {"item": "patients", "public": "yes (1 array = 1 title)", "n": n_array,
         "note": "GEO overall_design states 40 AC + 18 SCC; no duplicate-patient sentence"},
        {"item": "adenocarcinoma", "public": "yes", "n": n_adc,
         "note": "characteristics_ch1 = disease state: adenocarcinoma"},
        {"item": "squamous_cell_carcinoma", "public": "yes", "n": n_scc,
         "note": "characteristics_ch1 = disease state: squamous cell carcinoma"},
        {"item": "normals / paired adjacent", "public": "no", "n": 0,
         "note": "source_name = human non-small cell lung cancer tumor tissue for all 58"},
        {"item": "OS / DSS time or event", "public": "no", "n": 0,
         "note": "junction / histology paper; survival not deposited"},
        {"item": "ICI / treatment response", "public": "no", "n": 0,
         "note": "resected atlas, not an ICI series"},
        {"item": "numeric tumor % / purity", "public": "no", "n": 0,
         "note": "protocol text says tumor/stromal content was estimated then macro-dissected; values not deposited"},
        {"item": "Lung_Cancer_Explorer_n48", "public": "external subset", "n": 48,
         "note": "LCE lists 48 tumors for this accession; that is not the GEO series matrix. Not used."},
        {"item": "CLDN4 finite (collapsed)", "public": "yes", "n": int(sample["cldn4"].notna().sum()),
         "note": f"max-mean chose {confirm.loc[(confirm.gene=='CLDN4') & confirm.chosen_for_collapse, 'probe'].iloc[0]}"},
        {"item": "CD8A finite", "public": "yes", "n": int(sample["cd8a"].notna().sum()),
         "note": "named U133 probe 205758_at; only CD8A probe on GPL570"},
        {"item": "CD274 finite (collapsed)", "public": "yes", "n": int(sample["cd274"].notna().sum()),
         "note": f"max-mean chose {confirm.loc[(confirm.gene=='CD274') & confirm.chosen_for_collapse, 'probe'].iloc[0]} of 223834_at / 227458_at"},
        {"item": "primary pairwise n (all NSCLC)", "public": "yes", "n": n_array,
         "note": "complete-case CLDN4 + CD8A + CD274 = 58 / 58"},
        {"item": "primary pairwise n (ADC)", "public": "yes", "n": n_adc, "note": "HOLDS-rule eligible (n>=40)"},
        {"item": "primary pairwise n (SCC)", "public": "yes", "n": n_scc,
         "note": "UNDERPOWERED for HOLDS rule (n>=40)"},
    ])
    labels.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    # Coverage
    cov_rows = []
    for g in ["CLDN4", "CD8A", "CD274", "TACSTD2"] + EPITHELIAL:
        cov_rows.append({
            "gene": g,
            "present": int(g in gene_expr.index),
            "collapse_probe": (
                str(probe_audit.loc[probe_audit["gene"] == g, "probe"].iloc[0])
                if g in gene_expr.index else ""
            ),
        })
    pd.DataFrame(cov_rows).to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)

    # Correlations
    rows = []
    subsets = {
        "all_NSCLC": meta.index,
        "ADC": meta.index[meta["histology"] == "ADC"],
        "SCC": meta.index[meta["histology"] == "SCC"],
    }
    histo_all = (meta["histology"] == "ADC").astype(float)

    for subset, idx in subsets.items():
        cldn4 = gene_expr.loc["CLDN4", idx]
        epi = epithelial.reindex(idx)
        hbin = histo_all.reindex(idx).values if subset == "all_NSCLC" else None
        for partner in PRIMARY:
            y = gene_expr.loc[partner, idx]
            rows.append(pair_row(cldn4, y, epi, hbin, subset, partner, "collapsed_maxmean"))
        # named-probe sensitivity
        rows.append(pair_row(
            probe_expr.loc["201428_at", idx],
            probe_expr.loc["205758_at", idx],
            epi, hbin, subset, "CD8A", "named_201428_at_vs_205758_at",
        ))
        for p274 in CD274_PROBES:
            rows.append(pair_row(
                probe_expr.loc["201428_at", idx],
                probe_expr.loc[p274, idx],
                epi, hbin, subset, "CD274", f"named_201428_at_vs_{p274}",
            ))
        if "TACSTD2" in gene_expr.index:
            rows.append(pair_row(
                gene_expr.loc["TACSTD2", idx],
                gene_expr.loc["CD8A", idx],
                epi, hbin, subset, "CD8A", "TACSTD2_companion_vs_CD8A",
            ))
            rows.append(pair_row(
                gene_expr.loc["TACSTD2", idx],
                gene_expr.loc["CD274", idx],
                epi, hbin, subset, "CD274", "TACSTD2_companion_vs_CD274",
            ))
            rows.append(pair_row(
                gene_expr.loc["CLDN4", idx],
                gene_expr.loc["TACSTD2", idx],
                epi, hbin, subset, "TACSTD2", "CLDN4_vs_TACSTD2_companion",
            ))

    # CLDN4 vs epithelial / histology (context, not a claim)
    for subset, idx in subsets.items():
        cldn4 = gene_expr.loc["CLDN4", idx]
        epi = epithelial.reindex(idx)
        crude = spearman_ci(cldn4.values, epi.values)
        rows.append({
            "subset": subset, "partner": "epithelial_z", "source": "context",
            "n": crude["n"], "rho": crude["rho"], "p": crude["p"],
            "ci_low": crude["ci_low"], "ci_high": crude["ci_high"],
            "rho_adj_epithelial": np.nan, "p_adj_epithelial": np.nan,
            "rho_adj_histology": np.nan, "p_adj_histology": np.nan,
            "verdict_epithelial": "COMPANION", "verdict_histology": "NA",
        })
    crude_h = spearman_ci(gene_expr.loc["CLDN4"].values, histo_all.values)
    rows.append({
        "subset": "all_NSCLC", "partner": "histology_ADC", "source": "context",
        "n": crude_h["n"], "rho": crude_h["rho"], "p": crude_h["p"],
        "ci_low": crude_h["ci_low"], "ci_high": crude_h["ci_high"],
        "rho_adj_epithelial": np.nan, "p_adj_epithelial": np.nan,
        "rho_adj_histology": np.nan, "p_adj_histology": np.nan,
        "verdict_epithelial": "COMPANION", "verdict_histology": "NA",
    })

    corr = pd.DataFrame(rows)
    corr.to_csv(TABLES / "spearman.tsv", sep="\t", index=False)

    # Q4 vs Q1 on CD8A / CD274 (all NSCLC and ADC only)
    hl_rows = []
    for subset, idx in subsets.items():
        cldn4 = gene_expr.loc["CLDN4", idx]
        if cldn4.size < 16:
            continue
        q = pd.qcut(cldn4, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
        if q.nunique() < 4:
            continue
        for partner in PRIMARY:
            a = gene_expr.loc[partner, idx][q == "Q4"]
            b = gene_expr.loc[partner, idx][q == "Q1"]
            U, p_mw = stats.mannwhitneyu(a.values, b.values, alternative="two-sided")
            rbc = 2 * U / (a.size * b.size) - 1
            hl_rows.append({
                "subset": subset,
                "endpoint": partner,
                "n_Q4": int(a.size),
                "n_Q1": int(b.size),
                "median_Q4": float(np.median(a)),
                "median_Q1": float(np.median(b)),
                "U": float(U),
                "p": float(p_mw),
                "rank_biserial_Q4_minus_Q1": float(rbc),
            })
    hl = pd.DataFrame(hl_rows)
    hl.to_csv(TABLES / "highlow_cldn4.tsv", sep="\t", index=False)

    # Histology contrast on CLDN4 / CD8A / CD274
    hist_rows = []
    for gene in ["CLDN4", "CD8A", "CD274", "TACSTD2"]:
        if gene not in gene_expr.index:
            continue
        a = gene_expr.loc[gene, meta.index[meta["histology"] == "ADC"]]
        b = gene_expr.loc[gene, meta.index[meta["histology"] == "SCC"]]
        U, p_mw = stats.mannwhitneyu(a.values, b.values, alternative="two-sided")
        rbc = 2 * U / (a.size * b.size) - 1
        hist_rows.append({
            "gene": gene,
            "n_ADC": int(a.size),
            "n_SCC": int(b.size),
            "median_ADC": float(np.median(a)),
            "median_SCC": float(np.median(b)),
            "U": float(U),
            "p": float(p_mw),
            "rank_biserial_ADC_minus_SCC": float(rbc),
        })
    hist_df = pd.DataFrame(hist_rows)
    hist_df.to_csv(TABLES / "histology_contrast.tsv", sep="\t", index=False)

    # Figures
    colors = {"ADC": "#2c5f8a", "SCC": "#b35c1e"}
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
    for ax, partner, ylab in zip(axes, ["CD8A", "CD274"], ["CD8A  gcRMA", "CD274  gcRMA (max-mean)"]):
        for h, c in colors.items():
            m = meta["histology"] == h
            ax.scatter(
                gene_expr.loc["CLDN4", m],
                gene_expr.loc[partner, m],
                s=22, alpha=0.8, c=c, edgecolors="none", label=f"{h} n={int(m.sum())}",
            )
        row = corr[(corr.subset == "all_NSCLC") & (corr.partner == partner) & (corr.source == "collapsed_maxmean")].iloc[0]
        ax.set_xlabel("CLDN4  gcRMA")
        ax.set_ylabel(ylab)
        ax.set_title(
            f"all NSCLC n={int(row.n)}\n"
            f"ρ={row.rho:+.3f}  p={fmt_p(row.p)}"
        )
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_cd274.png", dpi=160)
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_cd274.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    plot = corr[(corr.source == "collapsed_maxmean") & (corr.partner.isin(PRIMARY))].copy()
    y_pos = np.arange(len(plot))
    ax.axvline(0, color="0.5", lw=0.8)
    ax.errorbar(
        plot["rho"], y_pos,
        xerr=[plot["rho"] - plot["ci_low"], plot["ci_high"] - plot["rho"]],
        fmt="o", color="#2c5f8a", ecolor="#2c5f8a", capsize=2, ms=5,
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{r.subset}  CLDN4 vs {r.partner}  n={int(r.n)}" for r in plot.itertuples()])
    ax.set_xlabel("Spearman ρ (bootstrap 95% CI)")
    ax.set_title("GSE10245  CLDN4 vs CD8A / CD274")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_forest.png", dpi=160)
    fig.savefig(FIGURES / "fig2_forest.pdf")
    plt.close(fig)

    def pick(subset, partner, source="collapsed_maxmean"):
        hit = corr[(corr.subset == subset) & (corr.partner == partner) & (corr.source == source)]
        return hit.iloc[0].to_dict() if len(hit) else {}

    summary = {
        "dataset": "GSE10245",
        "pmid": "18486272",
        "platform": "GPL570 Affymetrix HG-U133 Plus 2.0",
        "processing": "gcRMA as deposited (Bioconductor). Spearman is rank-based.",
        "n_arrays": n_array,
        "n_ADC": n_adc,
        "n_SCC": n_scc,
        "n_probes": int(probe_expr.shape[0]),
        "n_genes_collapsed": int(gene_expr.shape[0]),
        "n_missing_values": n_nan,
        "epithelial_genes_used": epi_used,
        "cldn4_collapse_probe": str(probe_audit.loc[probe_audit["gene"] == "CLDN4", "probe"].iloc[0]),
        "cd8a_collapse_probe": str(probe_audit.loc[probe_audit["gene"] == "CD8A", "probe"].iloc[0]),
        "cd274_collapse_probe": str(probe_audit.loc[probe_audit["gene"] == "CD274", "probe"].iloc[0]),
        "holds_rule": "n>=40, residual ρ<0, residual p<0.05 (same as PR 229 / GSE4573)",
        "primary": {
            "all_NSCLC_CD8A": pick("all_NSCLC", "CD8A"),
            "all_NSCLC_CD274": pick("all_NSCLC", "CD274"),
            "ADC_CD8A": pick("ADC", "CD8A"),
            "ADC_CD274": pick("ADC", "CD274"),
            "SCC_CD8A": pick("SCC", "CD8A"),
            "SCC_CD274": pick("SCC", "CD274"),
        },
        "missing_public_labels": ["OS", "ICI", "numeric_tumor_percent", "stage"],
        "do_not_use_n": {"Lung_Cancer_Explorer": 48, "reason": "external reprocessed subset, not the GEO series matrix"},
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")

    print("n_arrays", n_array, "ADC", n_adc, "SCC", n_scc, "genes", gene_expr.shape[0])
    print(corr[corr.source == "collapsed_maxmean"][
        ["subset", "partner", "n", "rho", "p", "rho_adj_epithelial", "p_adj_epithelial",
         "rho_adj_histology", "p_adj_histology", "verdict_epithelial", "verdict_histology"]
    ].to_string(index=False))
    print(hist_df.to_string(index=False))
    print("[done]", TABLES)


if __name__ == "__main__":
    main()
