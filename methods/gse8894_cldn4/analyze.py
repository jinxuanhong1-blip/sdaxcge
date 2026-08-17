#!/usr/bin/env python3
"""GSE8894 NSCLC Affy U133 Plus 2.0: CLDN4 vs CD8A / CD274.

Additive public bulk slice. Unit is the array. No ICI arm. Downloads stay
under $GSE8894_CLDN4_DATA (default /tmp/gse8894_cldn4) and are not committed.
"""

from __future__ import annotations

import gzip
import io
import json
import os
import urllib.request
from collections import defaultdict
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
DATA = Path(os.environ.get("GSE8894_CLDN4_DATA", "/tmp/gse8894_cldn4"))
SEED = 20260817
N_BOOT = 2000

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE8nnn/GSE8894/"
    "matrix/GSE8894_series_matrix.txt.gz"
)
GPL_SOFT = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    "?acc=GPL570&targ=self&form=text&view=data"
)

# Primary partners requested. TACSTD2 is companion only.
PRIMARY = ["CD8A", "CD274"]
COMPANION = ["TACSTD2"]
EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]

NAMED_PROBES = {
    "CLDN4": "201428_at",
    "CD8A": "205758_at",
    "CD274": "223834_at",
    "TACSTD2": "202286_s_at",
}

HOLDS_N = 40


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "gse8894-cldn4/1.0"})
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
    sample_fields: dict[str, list[list[str]]] = defaultdict(list)
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
                sample_fields[key].append(vals)

    n = len(sample_fields["geo_accession"][0])
    rows = []
    for i in range(n):
        rec = {
            "gsm": sample_fields["geo_accession"][0][i],
            "title": sample_fields["title"][0][i],
            "source_name_ch1": sample_fields["source_name_ch1"][0][i].strip(),
            "description": sample_fields.get("description", [[""]])[0][i],
            "data_processing": sample_fields["data_processing"][0][i].strip(),
        }
        for row in sample_fields["characteristics_ch1"]:
            v = row[i]
            if not v or ":" not in v:
                continue
            k, val = v.split(":", 1)
            rec[k.strip().lower()] = val.strip()
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
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def histolabel(s: str) -> str:
    s = (s or "").strip().lower()
    if s == "adenocarcinoma":
        return "ADC"
    if s == "squamous cell carcinoma":
        return "SCC"
    return s or "NA"


def pair_row(cldn4, y, epithelial, gene_y: str, stratum: str) -> dict:
    crude = spearman_ci(cldn4.values, y.values)
    adj = (
        partial_spearman(cldn4.values, y.values, epithelial.values)
        if epithelial is not None
        else dict(n=crude["n"], rho_adj=np.nan, p_adj=np.nan)
    )
    return {
        "stratum": stratum,
        "gene_x": "CLDN4",
        "gene_y": gene_y,
        "n": crude["n"],
        "rho": crude["rho"],
        "p": crude["p"],
        "ci_low": crude["ci_low"],
        "ci_high": crude["ci_high"],
        "rho_adj_epithelial": adj["rho_adj"],
        "p_adj_epithelial": adj["p_adj"],
        "verdict": verdict(adj["rho_adj"], adj["p_adj"], adj["n"]),
    }


def highlow(cldn4: pd.Series, y: pd.Series, endpoint: str, stratum: str) -> dict:
    q = pd.qcut(cldn4, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    a = y[q == "Q4"]
    b = y[q == "Q1"]
    if a.size < 5 or b.size < 5:
        return {
            "stratum": stratum,
            "endpoint": endpoint,
            "n_Q4": int(a.size),
            "n_Q1": int(b.size),
            "median_Q4": np.nan,
            "median_Q1": np.nan,
            "U": np.nan,
            "p": np.nan,
            "rank_biserial_Q4_minus_Q1": np.nan,
        }
    U, p_mw = stats.mannwhitneyu(a.values, b.values, alternative="two-sided")
    rbc = 2 * U / (a.size * b.size) - 1
    return {
        "stratum": stratum,
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

    matrix_path = dl(DATA / "GSE8894_series_matrix.txt.gz", GEO_MATRIX)
    annot_path = dl(DATA / "GPL570_query.txt", GPL_SOFT)

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
    meta["histology"] = meta["cell type"].map(histolabel)
    meta["recurrence"] = meta["status (1=recurrence, 0=non-recurrence)"].map(
        {"0": "non-recurrence", "1": "recurrence"}
    )
    meta["rfs_months"] = pd.to_numeric(
        meta.get("recurrence free survival time (month)"), errors="coerce"
    )
    meta["age"] = pd.to_numeric(meta.get("age"), errors="coerce")

    with open(annot_path, "r", errors="replace") as fh:
        annot = parse_soft_table(fh.read())
    annot["symbol"] = annot["Gene Symbol"].map(first_symbol)
    id2gene = annot.drop_duplicates("ID").set_index("ID")["symbol"]

    gene_expr, probe_audit = collapse_maxmean(probe_expr, id2gene)

    confirm_rows = []
    for gene, named in NAMED_PROBES.items():
        hits = annot[annot["symbol"] == gene][["ID", "Gene Symbol", "Gene Title", "ENTREZ_GENE_ID"]]
        for _, r in hits.iterrows():
            confirm_rows.append(
                {
                    "gene": gene,
                    "named_probe": named,
                    "probe": r["ID"],
                    "gpl570_symbol": r["Gene Symbol"],
                    "entrez": r["ENTREZ_GENE_ID"],
                    "title": r["Gene Title"],
                    "in_matrix": r["ID"] in probe_expr.index,
                    "is_named": r["ID"] == named,
                    "chosen_for_collapse": (
                        gene in gene_expr.index
                        and probe_audit.loc[probe_audit["gene"] == gene, "probe"].iloc[0] == r["ID"]
                    ),
                }
            )
    confirm = pd.DataFrame(confirm_rows)
    confirm.to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)

    epi_present = [g for g in EPITHELIAL if g in gene_expr.index]
    z = gene_expr.loc[epi_present]
    mu = z.mean(axis=1)
    sd = z.std(axis=1, ddof=1).replace(0, np.nan)
    epithelial = z.sub(mu, axis=0).div(sd, axis=0).mean(axis=0)

    n_array = int(probe_expr.shape[1])
    n_adc = int((meta["histology"] == "ADC").sum())
    n_scc = int((meta["histology"] == "SCC").sum())
    n_rec = int((meta["recurrence"] == "recurrence").sum())
    n_nrec = int((meta["recurrence"] == "non-recurrence").sum())
    n_age = int(meta["age"].notna().sum())
    n_rfs = int(meta["rfs_months"].notna().sum())
    n_missing = int(probe_expr.isna().sum().sum())

    for g in ["CLDN4", "CD8A", "CD274"]:
        if g not in gene_expr.index:
            raise SystemExit(f"required gene absent after collapse: {g}")

    sample = meta.copy()
    for g, col in [("CLDN4", "cldn4"), ("CD8A", "cd8a"), ("CD274", "cd274"), ("TACSTD2", "tacstd2")]:
        sample[col] = gene_expr.loc[g].reindex(sample.index).values if g in gene_expr.index else np.nan
    sample["epithelial_mean_z"] = epithelial.reindex(sample.index).values
    sample.to_csv(TABLES / "sample_annotation.tsv", sep="\t")

    labels = pd.DataFrame(
        [
            {
                "field": "arrays_in_matrix",
                "public": "yes",
                "n": n_array,
                "note": f"{probe_expr.shape[0]} probes x {n_array} GSM; missing values={n_missing}",
            },
            {
                "field": "unique_GSM",
                "public": "yes",
                "n": int(meta.index.nunique()),
                "note": "all unique",
            },
            {
                "field": "unique_title",
                "public": "yes",
                "n": int(meta["title"].nunique()),
                "note": "LUNG CANCER SML* titles; all unique",
            },
            {
                "field": "tissue",
                "public": "yes",
                "n": int((meta["source_name_ch1"].str.strip() == "frozen tissue of primary lung tumor").sum()),
                "note": "every array is frozen primary lung tumor; 0 adjacent/normal",
            },
            {
                "field": "histology_ADC_GEO",
                "public": "yes",
                "n": n_adc,
                "note": "GEO `cell type: Adenocarcinoma`. Paper often cited as 69/69; GEO is 63 ADC",
            },
            {
                "field": "histology_SCC_GEO",
                "public": "yes",
                "n": n_scc,
                "note": "GEO `cell type: squamous cell carcinoma`. GEO is 75 SCC, not 69",
            },
            {
                "field": "recurrence_event",
                "public": "yes",
                "n": n_rec,
                "note": "status 1=recurrence; 69/69 split. Not the primary claim",
            },
            {
                "field": "non_recurrence",
                "public": "yes",
                "n": n_nrec,
                "note": "status 0=non-recurrence",
            },
            {
                "field": "RFS_months",
                "public": "yes",
                "n": n_rfs,
                "note": "recurrence free survival time (month) on every array",
            },
            {
                "field": "age",
                "public": "yes (incomplete)",
                "n": n_age,
                "note": "missing on GSM225774 and GSM225864",
            },
            {
                "field": "gender",
                "public": "yes",
                "n": n_array,
                "note": f"Male {(meta['gender']=='Male').sum()}; Female {(meta['gender']=='Female').sum()}",
            },
            {
                "field": "stage",
                "public": "no",
                "n": 0,
                "note": "not a GEO sample characteristic",
            },
            {
                "field": "ICI / treatment",
                "public": "no",
                "n": 0,
                "note": "resected 2008 atlas, not an ICI series",
            },
            {
                "field": "tumor_percent / purity",
                "public": "no",
                "n": 0,
                "note": "epithelial RNA mean-z is the only public purity proxy used here",
            },
            {
                "field": "CLDN4 finite",
                "public": "yes",
                "n": int(sample["cldn4"].notna().sum()),
                "note": "named probe 201428_at; max-mean collapse",
            },
            {
                "field": "CD8A finite",
                "public": "yes",
                "n": int(sample["cd8a"].notna().sum()),
                "note": "named probe 205758_at",
            },
            {
                "field": "CD274 finite",
                "public": "yes",
                "n": int(sample["cd274"].notna().sum()),
                "note": "named probe 223834_at; second probe 227458_at also on GPL570",
            },
            {
                "field": "pairwise_CLDN4_CD8A",
                "public": "yes",
                "n": int((sample["cldn4"].notna() & sample["cd8a"].notna()).sum()),
                "note": "primary n for CLDN4 vs CD8A",
            },
            {
                "field": "pairwise_CLDN4_CD274",
                "public": "yes",
                "n": int((sample["cldn4"].notna() & sample["cd274"].notna()).sum()),
                "note": "primary n for CLDN4 vs CD274",
            },
        ]
    )
    labels.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    cldn4 = gene_expr.loc["CLDN4"]
    rows = []
    hl_rows = []
    strata = {
        "all": meta.index,
        "ADC": meta.index[meta["histology"] == "ADC"],
        "SCC": meta.index[meta["histology"] == "SCC"],
    }
    partners = PRIMARY + [g for g in COMPANION if g in gene_expr.index] + ["epithelial_mean_z"]
    for stratum, idx in strata.items():
        epi_s = epithelial.reindex(idx)
        c = cldn4.reindex(idx)
        for g in partners:
            if g == "epithelial_mean_z":
                crude = spearman_ci(c.values, epi_s.values)
                rows.append(
                    {
                        "stratum": stratum,
                        "gene_x": "CLDN4",
                        "gene_y": "epithelial_mean_z",
                        "n": crude["n"],
                        "rho": crude["rho"],
                        "p": crude["p"],
                        "ci_low": crude["ci_low"],
                        "ci_high": crude["ci_high"],
                        "rho_adj_epithelial": np.nan,
                        "p_adj_epithelial": np.nan,
                        "verdict": "COMPANION",
                    }
                )
                continue
            if g not in gene_expr.index:
                rows.append(
                    {
                        "stratum": stratum,
                        "gene_x": "CLDN4",
                        "gene_y": g,
                        "n": 0,
                        "rho": np.nan,
                        "p": np.nan,
                        "ci_low": np.nan,
                        "ci_high": np.nan,
                        "rho_adj_epithelial": np.nan,
                        "p_adj_epithelial": np.nan,
                        "verdict": "ABSENT",
                    }
                )
                continue
            y = gene_expr.loc[g].reindex(idx)
            rec = pair_row(c, y, epi_s, g, stratum)
            if g == "TACSTD2":
                rec["verdict"] = "COMPANION"
            rows.append(rec)
            if g in PRIMARY:
                hl_rows.append(highlow(c, y, g, stratum))

        if "TACSTD2" in gene_expr.index:
            crude = spearman_ci(
                gene_expr.loc["TACSTD2"].reindex(idx).values,
                gene_expr.loc["CD8A"].reindex(idx).values,
            )
            adj = partial_spearman(
                gene_expr.loc["TACSTD2"].reindex(idx).values,
                gene_expr.loc["CD8A"].reindex(idx).values,
                epi_s.values,
            )
            rows.append(
                {
                    "stratum": stratum,
                    "gene_x": "TACSTD2",
                    "gene_y": "CD8A",
                    "n": crude["n"],
                    "rho": crude["rho"],
                    "p": crude["p"],
                    "ci_low": crude["ci_low"],
                    "ci_high": crude["ci_high"],
                    "rho_adj_epithelial": adj["rho_adj"],
                    "p_adj_epithelial": adj["p_adj"],
                    "verdict": "COMPANION",
                }
            )
            crude = spearman_ci(
                gene_expr.loc["TACSTD2"].reindex(idx).values,
                gene_expr.loc["CD274"].reindex(idx).values,
            )
            adj = partial_spearman(
                gene_expr.loc["TACSTD2"].reindex(idx).values,
                gene_expr.loc["CD274"].reindex(idx).values,
                epi_s.values,
            )
            rows.append(
                {
                    "stratum": stratum,
                    "gene_x": "TACSTD2",
                    "gene_y": "CD274",
                    "n": crude["n"],
                    "rho": crude["rho"],
                    "p": crude["p"],
                    "ci_low": crude["ci_low"],
                    "ci_high": crude["ci_high"],
                    "rho_adj_epithelial": adj["rho_adj"],
                    "p_adj_epithelial": adj["p_adj"],
                    "verdict": "COMPANION",
                }
            )

    # Histology residual on the pooled set (ADC=1, SCC=0). This is not a purity
    # control; it asks whether the mixed-NSCLC rho is only ADC-vs-SCC mixing.
    hist_code = (meta["histology"] == "ADC").astype(float)
    for g in PRIMARY:
        crude = spearman_ci(cldn4.values, gene_expr.loc[g].values)
        adj_h = partial_spearman(cldn4.values, gene_expr.loc[g].values, hist_code.values)
        rows.append(
            {
                "stratum": "all_partial_histology",
                "gene_x": "CLDN4",
                "gene_y": g,
                "n": crude["n"],
                "rho": crude["rho"],
                "p": crude["p"],
                "ci_low": crude["ci_low"],
                "ci_high": crude["ci_high"],
                "rho_adj_epithelial": adj_h["rho_adj"],
                "p_adj_epithelial": adj_h["p_adj"],
                "verdict": verdict(adj_h["rho_adj"], adj_h["p_adj"], adj_h["n"]),
            }
        )

    # Named-probe CD274 sensitivity (223834_at). Collapse chose 227458_at.
    if "223834_at" in probe_expr.index:
        y_named = probe_expr.loc["223834_at"].astype(float)
        rec = pair_row(cldn4, y_named.reindex(cldn4.index), epithelial, "CD274_223834_at", "all")
        rec["verdict"] = "SENSITIVITY"
        rows.append(rec)

    genes_df = pd.DataFrame(rows)
    genes_df.to_csv(TABLES / "spearman_cldn4_vs_cd8a_cd274.tsv", sep="\t", index=False)
    hl = pd.DataFrame(hl_rows)
    hl.to_csv(TABLES / "highlow_cldn4.tsv", sep="\t", index=False)

    cov = pd.DataFrame(
        [
            {
                "set": f"gene:{g}",
                "n_listed": 1,
                "n_present": int(g in gene_expr.index),
                "missing": "" if g in gene_expr.index else g,
                "present": g if g in gene_expr.index else "",
                "collapse_probe": (
                    probe_audit.loc[probe_audit["gene"] == g, "probe"].iloc[0]
                    if g in gene_expr.index
                    else ""
                ),
            }
            for g in ["CLDN4", "CD8A", "CD274", "TACSTD2"] + EPITHELIAL
        ]
    )
    cov.to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)

    # Figures
    def scatter(x, y, xlabel, ylabel, title, path_stem, color):
        fig, ax = plt.subplots(figsize=(5.2, 4.4))
        ax.scatter(x, y, s=18, alpha=0.7, c=color, edgecolors="none")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        fig.tight_layout()
        fig.savefig(FIGURES / f"{path_stem}.png", dpi=160)
        fig.savefig(FIGURES / f"{path_stem}.pdf")
        plt.close(fig)

    r_cd8 = genes_df[(genes_df.stratum == "all") & (genes_df.gene_x == "CLDN4") & (genes_df.gene_y == "CD8A")].iloc[0]
    r_pdl1 = genes_df[(genes_df.stratum == "all") & (genes_df.gene_x == "CLDN4") & (genes_df.gene_y == "CD274")].iloc[0]
    scatter(
        cldn4.values,
        gene_expr.loc["CD8A"].values,
        "CLDN4  GCRMA",
        "CD8A  GCRMA",
        f"GSE8894 NSCLC  n={int(r_cd8.n)} arrays\nCLDN4 vs CD8A  ρ={r_cd8.rho:.3f}  p={fmt_p(r_cd8.p)}",
        "fig1_cldn4_vs_cd8a",
        "#2c5f8a",
    )
    scatter(
        cldn4.values,
        gene_expr.loc["CD274"].values,
        "CLDN4  GCRMA",
        "CD274  GCRMA",
        f"GSE8894 NSCLC  n={int(r_pdl1.n)} arrays\nCLDN4 vs CD274  ρ={r_pdl1.rho:.3f}  p={fmt_p(r_pdl1.p)}",
        "fig2_cldn4_vs_cd274",
        "#6b3a2a",
    )

    plot = genes_df[
        (genes_df.gene_x == "CLDN4")
        & (genes_df.gene_y.isin(PRIMARY))
        & (genes_df.stratum.isin(["all", "ADC", "SCC"]))
    ].copy()
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
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
    ax.set_yticklabels([f"{r.stratum}  CLDN4–{r.gene_y}  n={int(r.n)}" for r in plot.itertuples()])
    ax.set_xlabel("Spearman ρ vs CLDN4 (bootstrap 95% CI)")
    ax.set_title("GSE8894  CLDN4 vs CD8A / CD274")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_cldn4_forest.png", dpi=160)
    fig.savefig(FIGURES / "fig3_cldn4_forest.pdf")
    plt.close(fig)

    one = []
    for stratum in ("all", "ADC", "SCC"):
        rec = {"dataset": "GSE8894", "stratum": stratum}
        rec["n"] = int((meta["histology"] == stratum).sum()) if stratum != "all" else n_array
        for g in PRIMARY:
            row = genes_df[
                (genes_df.stratum == stratum) & (genes_df.gene_x == "CLDN4") & (genes_df.gene_y == g)
            ].iloc[0]
            rec[f"CLDN4_{g}_rho"] = row["rho"]
            rec[f"CLDN4_{g}_p"] = row["p"]
            rec[f"CLDN4_{g}_rho_adj"] = row["rho_adj_epithelial"]
            rec[f"CLDN4_{g}_p_adj"] = row["p_adj_epithelial"]
            rec[f"CLDN4_{g}_verdict"] = row["verdict"]
        one.append(rec)
    one_df = pd.DataFrame(one)
    one_df.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    chosen = {}
    for g in ["CLDN4", "CD8A", "CD274", "TACSTD2"]:
        if g in gene_expr.index:
            chosen[g] = str(probe_audit.loc[probe_audit["gene"] == g, "probe"].iloc[0])

    summary = {
        "dataset": "GSE8894",
        "pmid": "19010856",
        "platform": "GPL570 Affymetrix HG-U133 Plus 2.0",
        "processing": "GCRMA as deposited. Spearman is rank-based.",
        "histology_geo": {"ADC": n_adc, "SCC": n_scc},
        "paper_often_cited": "138 NSCLC; many secondary sources write 69 ADC + 69 SCC",
        "honest_n_geo": {
            "arrays": n_array,
            "ADC": n_adc,
            "SCC": n_scc,
            "note": "Use GEO cell-type counts (63 ADC / 75 SCC), not 69/69",
        },
        "n_probes": int(probe_expr.shape[0]),
        "n_genes_collapsed": int(gene_expr.shape[0]),
        "n_missing_values": n_missing,
        "epithelial_genes_used": epi_present,
        "collapse_probes": chosen,
        "named_probes": NAMED_PROBES,
        "primary": {
            "CLDN4_vs_CD8A": genes_df[
                (genes_df.stratum == "all") & (genes_df.gene_x == "CLDN4") & (genes_df.gene_y == "CD8A")
            ]
            .iloc[0]
            .to_dict(),
            "CLDN4_vs_CD274": genes_df[
                (genes_df.stratum == "all") & (genes_df.gene_x == "CLDN4") & (genes_df.gene_y == "CD274")
            ]
            .iloc[0]
            .to_dict(),
        },
        "missing_public_labels": ["stage", "ICI", "tumor_percent"],
        "holds_rule": "n>=40, partial Spearman rho<0, p<0.05 (same as PR 229 / GSE4573)",
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    print("n_arrays", n_array, "ADC", n_adc, "SCC", n_scc, "genes", gene_expr.shape[0])
    print(genes_df.to_string(index=False))
    print(hl.to_string(index=False))
    print("[done]", TABLES)


if __name__ == "__main__":
    main()
