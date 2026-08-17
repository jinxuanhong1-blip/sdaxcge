#!/usr/bin/env python3
"""GSE3141 NSCLC Affy U133 Plus 2.0: CLDN4 vs CD8A / CD274.

Additive public bulk slice. Unit is the array. No ICI arm. No slide
was re-scored. Downloads stay under $GSE3141_CLDN4_DATA (default
/tmp/gse3141_cldn4) and are not committed.
"""

from __future__ import annotations

import gzip
import io
import json
import os
import re
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
DATA = Path(os.environ.get("GSE3141_CLDN4_DATA", "/tmp/gse3141_cldn4"))
SEED = 20260817
N_BOOT = 2000

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE3nnn/GSE3141/"
    "matrix/GSE3141_series_matrix.txt.gz"
)
GPL_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"

# Primary partners requested. TACSTD2 is companion only.
PRIMARY = ["CD8A", "CD274"]
COMPANION = ["TACSTD2"]
EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]

NAMED_PROBES = {
    "CLDN4": "201428_at",
    "CD8A": "205758_at",
    "CD274": "227458_at",  # higher-mean GPL570 probe; 223834_at is the other
    "TACSTD2": "202286_s_at",
}
CD274_PROBES = ["223834_at", "227458_at"]
CLDN4_PROBES = ["201428_at", "1569421_at"]

# Same HOLDS rule as PR 229 / GSE4573 CLDN4 slice: n>=40, adj ρ<0, p<0.05.
HOLDS_N = 40


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "gse3141-cldn4/1.0"})
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


def parse_characteristics(raw: str) -> dict:
    cell = re.search(r"Cell type:\s*([^;]+)", raw)
    surv = re.search(r"Surv\(months\):\s*([^;]+)", raw)
    status = re.search(r"STATUS\(0=alive, 1=dead\):\s*([^;]+)", raw)
    cell_code = cell.group(1).strip() if cell else ""
    hist = {"A": "LUAD", "S": "LUSC"}.get(cell_code, cell_code or "unknown")
    surv_raw = surv.group(1).strip() if surv else ""
    surv_num = np.nan
    if re.match(r"^-?\d", surv_raw):
        surv_num = float(surv_raw)
    return {
        "cell_type_code": cell_code,
        "histology": hist,
        "surv_months_raw": surv_raw,
        "surv_months": surv_num,
        "surv_numeric": bool(np.isfinite(surv_num)),
        "os_event": int(status.group(1).strip()) if status else np.nan,
    }


def pair_row(x, y, covar, gene_x, gene_y, subset, n_arrays, extra=None):
    crude = spearman_ci(x, y)
    adj = partial_spearman(x, y, covar)
    rec = {
        "gene_x": gene_x,
        "gene_y": gene_y,
        "subset": subset,
        "n": crude["n"],
        "n_arrays": n_arrays,
        "rho": crude["rho"],
        "p": crude["p"],
        "ci_low": crude["ci_low"],
        "ci_high": crude["ci_high"],
        "rho_adj_epithelial": adj["rho_adj"],
        "p_adj_epithelial": adj["p_adj"],
        "verdict": verdict(adj["rho_adj"], adj["p_adj"], adj["n"]),
    }
    if extra:
        rec.update(extra)
    return rec


def highlow(x, y, endpoint: str) -> dict:
    q = pd.qcut(pd.Series(x), 4, labels=["Q1", "Q2", "Q3", "Q4"])
    a = pd.Series(y)[q == "Q4"]
    b = pd.Series(y)[q == "Q1"]
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

    matrix_path = dl(DATA / "GSE3141_series_matrix.txt.gz", GEO_MATRIX)
    annot_path = dl(DATA / "GPL570.annot.gz", GPL_ANNOT)

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

    parsed = meta["characteristics_ch1"].map(parse_characteristics).apply(pd.Series)
    meta = pd.concat([meta, parsed], axis=1)

    with gzip.open(annot_path, "rt", errors="replace") as fh:
        annot = parse_soft_table(fh.read())
    annot["symbol"] = annot["Gene symbol"].map(first_symbol)
    id2gene = annot.drop_duplicates("ID").set_index("ID")["symbol"]

    n_missing = int((~np.isfinite(probe_expr.to_numpy())).sum())
    gene_expr, probe_audit = collapse_maxmean(probe_expr, id2gene)

    confirm_rows = []
    for gene, named in {**NAMED_PROBES, "CD274_alt": "223834_at", "CLDN4_alt": "1569421_at"}.items():
        gene_key = gene.replace("_alt", "")
        hits = annot[annot["symbol"] == gene_key][["ID", "Gene symbol", "Gene title", "Gene ID"]]
        for _, r in hits.iterrows():
            chosen = (
                gene_key in gene_expr.index
                and probe_audit.loc[probe_audit["gene"] == gene_key, "probe"].iloc[0] == r["ID"]
            )
            confirm_rows.append({
                "gene": gene_key,
                "named_probe": named if gene in NAMED_PROBES else "",
                "probe": r["ID"],
                "gpl570_symbol": r["Gene symbol"],
                "entrez": r["Gene ID"],
                "title": r["Gene title"],
                "in_matrix": r["ID"] in probe_expr.index,
                "n_finite": int(np.isfinite(probe_expr.loc[r["ID"]].to_numpy()).sum())
                if r["ID"] in probe_expr.index
                else 0,
                "mean_mas5": float(probe_expr.loc[r["ID"]].mean()) if r["ID"] in probe_expr.index else np.nan,
                "is_named_primary": r["ID"] == NAMED_PROBES.get(gene_key, ""),
                "chosen_for_collapse": chosen,
            })
    confirm = pd.DataFrame(confirm_rows).drop_duplicates(["gene", "probe"])
    confirm.to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)

    epi_genes = [g for g in EPITHELIAL if g in gene_expr.index]
    z = gene_expr.loc[epi_genes]
    mu = z.mean(axis=1)
    sd = z.std(axis=1, ddof=1).replace(0, np.nan)
    epithelial = z.sub(mu, axis=0).div(sd, axis=0).mean(axis=0)

    sample = meta[
        [
            "geo_accession",
            "title",
            "source_name_ch1",
            "characteristics_ch1",
            "cell_type_code",
            "histology",
            "surv_months_raw",
            "surv_months",
            "surv_numeric",
            "os_event",
            "data_processing",
        ]
    ].copy()
    for g in ["CLDN4", "CD8A", "CD274", "TACSTD2"]:
        sample[f"{g.lower()}_mas5"] = (
            gene_expr.loc[g].reindex(sample.index).values if g in gene_expr.index else np.nan
        )
    for p in CD274_PROBES + [NAMED_PROBES["CLDN4"], NAMED_PROBES["CD8A"]]:
        sample[p] = probe_expr.loc[p].reindex(sample.index).values if p in probe_expr.index else np.nan
    sample["epithelial_mean_z"] = epithelial.reindex(sample.index).values
    sample.to_csv(TABLES / "sample_annotation.tsv", sep="\t")

    n_array = int(probe_expr.shape[1])
    n_luad = int((meta["histology"] == "LUAD").sum())
    n_lusc = int((meta["histology"] == "LUSC").sum())
    n_surv_num = int(meta["surv_numeric"].sum())
    n_cldn4 = int(sample["cldn4_mas5"].notna().sum())
    n_cd8a = int(sample["cd8a_mas5"].notna().sum())
    n_cd274 = int(sample["cd274_mas5"].notna().sum())
    n_pair_cd8 = int((sample["cldn4_mas5"].notna() & sample["cd8a_mas5"].notna()).sum())
    n_pair_pdl1 = int((sample["cldn4_mas5"].notna() & sample["cd274_mas5"].notna()).sum())

    labels = pd.DataFrame([
        {"field": "arrays_in_matrix", "public": "yes", "n": n_array,
         "note": f"{probe_expr.shape[0]} probes x {n_array} GSM; missing MAS5 values={n_missing}"},
        {"field": "unique_GSM", "public": "yes", "n": int(meta.index.nunique()), "note": "all unique"},
        {"field": "unique_title", "public": "yes", "n": int(meta["title"].nunique()), "note": "Duke case IDs; all unique"},
        {"field": "source_name", "public": "yes", "n": n_array,
         "note": "every row is 'frozen tissue of primary lung tumor'; no adjacent-normal rows"},
        {"field": "histology_LUAD_A", "public": "yes", "n": n_luad, "note": "Cell type: A"},
        {"field": "histology_LUSC_S", "public": "yes", "n": n_lusc, "note": "Cell type: S"},
        {"field": "histology_other", "public": "yes", "n": n_array - n_luad - n_lusc, "note": "none"},
        {"field": "OS_months_numeric", "public": "yes", "n": n_surv_num,
         "note": "110 numeric; GSM70230 Surv(months)='>72' (alive LUSC) is not a number"},
        {"field": "OS_event", "public": "yes", "n": int(meta["os_event"].notna().sum()),
         "note": "STATUS 0=alive / 1=dead on all 111 arrays"},
        {"field": "stage", "public": "no", "n": 0, "note": "not a GEO sample characteristic"},
        {"field": "ICI / treatment", "public": "no", "n": 0,
         "note": "Bild 2006 oncogenic-pathway NSCLC atlas; not an ICI series"},
        {"field": "tumor_percent / purity", "public": "no", "n": 0,
         "note": "epithelial RNA mean-z is the only public purity proxy used here"},
        {"field": "CLDN4 finite (collapse)", "public": "yes", "n": n_cldn4, "note": "max-mean of 201428_at / 1569421_at"},
        {"field": "CD8A finite", "public": "yes", "n": n_cd8a, "note": "named probe 205758_at"},
        {"field": "CD274 finite (collapse)", "public": "yes", "n": n_cd274, "note": "max-mean of 223834_at / 227458_at"},
        {"field": "pairwise CLDN4+CD8A", "public": "yes", "n": n_pair_cd8, "note": "primary n for CD8A"},
        {"field": "pairwise CLDN4+CD274", "public": "yes", "n": n_pair_pdl1, "note": "primary n for CD274"},
    ])
    labels.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    assert n_array == 111
    for g in ["CLDN4", "CD8A", "CD274"]:
        if g not in gene_expr.index:
            raise SystemExit(f"required gene absent after collapse: {g}")

    cldn4 = gene_expr.loc["CLDN4"]
    rows = []
    subsets = {
        "all_tumors": meta.index,
        "LUAD": meta.index[meta["histology"] == "LUAD"],
        "LUSC": meta.index[meta["histology"] == "LUSC"],
    }
    for subset_name, idx in subsets.items():
        epi = epithelial.reindex(idx)
        for g in PRIMARY + COMPANION:
            extra = {"role": "primary" if g in PRIMARY else "companion"}
            rows.append(
                pair_row(
                    cldn4.reindex(idx).values,
                    gene_expr.loc[g].reindex(idx).values,
                    epi.values,
                    "CLDN4",
                    g,
                    subset_name,
                    n_arrays=int(len(idx)),
                    extra=extra,
                )
            )
        rows.append(
            pair_row(
                cldn4.reindex(idx).values,
                epi.values,
                epi.values,
                "CLDN4",
                "epithelial_mean_z",
                subset_name,
                n_arrays=int(len(idx)),
                extra={"role": "covariate", "verdict": "COMPANION"},
            )
        )
        # overwrite verdict for epithelial (partial vs self is meaningless)
        rows[-1]["rho_adj_epithelial"] = np.nan
        rows[-1]["p_adj_epithelial"] = np.nan
        rows[-1]["verdict"] = "COMPANION"

    # Named-probe sensitivity on all tumors (CD274 has two GPL570 probes).
    for p in CD274_PROBES:
        rows.append(
            pair_row(
                probe_expr.loc[NAMED_PROBES["CLDN4"]].values,
                probe_expr.loc[p].values,
                epithelial.values,
                "CLDN4_201428_at",
                f"CD274_{p}",
                "all_tumors_named_probe",
                n_arrays=n_array,
                extra={"role": "probe_sensitivity"},
            )
        )
    rows.append(
        pair_row(
            probe_expr.loc[NAMED_PROBES["CLDN4"]].values,
            probe_expr.loc[NAMED_PROBES["CD8A"]].values,
            epithelial.values,
            "CLDN4_201428_at",
            "CD8A_205758_at",
            "all_tumors_named_probe",
            n_arrays=n_array,
            extra={"role": "probe_sensitivity"},
        )
    )
    if "TACSTD2" in gene_expr.index:
        rows.append(
            pair_row(
                gene_expr.loc["TACSTD2"].values,
                gene_expr.loc["CD8A"].values,
                epithelial.values,
                "TACSTD2",
                "CD8A",
                "all_tumors",
                n_arrays=n_array,
                extra={"role": "companion"},
            )
        )
        rows.append(
            pair_row(
                gene_expr.loc["TACSTD2"].values,
                gene_expr.loc["CD274"].values,
                epithelial.values,
                "TACSTD2",
                "CD274",
                "all_tumors",
                n_arrays=n_array,
                extra={"role": "companion"},
            )
        )

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(TABLES / "spearman_cldn4_vs_cd8a_cd274.tsv", sep="\t", index=False)

    hl = pd.DataFrame([
        highlow(cldn4.values, gene_expr.loc["CD8A"].values, "CD8A"),
        highlow(cldn4.values, gene_expr.loc["CD274"].values, "CD274"),
    ])
    hl.to_csv(TABLES / "highlow_cldn4.tsv", sep="\t", index=False)

    # Coverage
    cov = pd.DataFrame([
        {"set": f"gene:{g}", "n_listed": 1, "n_present": int(g in gene_expr.index),
         "probe_used": (probe_audit.loc[probe_audit["gene"] == g, "probe"].iloc[0]
                        if g in gene_expr.index else ""),
         "missing": "" if g in gene_expr.index else g}
        for g in ["CLDN4", "CD8A", "CD274", "TACSTD2"] + EPITHELIAL
    ])
    cov.to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)

    # Figures
    def scatter(ax, x, y, title, xlab, ylab, color):
        ax.scatter(np.log2(np.asarray(x) + 1), np.log2(np.asarray(y) + 1),
                   s=18, alpha=0.75, c=color, edgecolors="none")
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(title)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.3))
    r8 = stats_df[(stats_df.gene_y == "CD8A") & (stats_df.subset == "all_tumors") & (stats_df.gene_x == "CLDN4")].iloc[0]
    rp = stats_df[(stats_df.gene_y == "CD274") & (stats_df.subset == "all_tumors") & (stats_df.gene_x == "CLDN4")].iloc[0]
    scatter(
        axes[0], cldn4.values, gene_expr.loc["CD8A"].values,
        f"GSE3141 NSCLC  n={int(r8.n)}\nCLDN4 vs CD8A  ρ={r8.rho:.3f}  p={fmt_p(r8.p)}",
        "CLDN4  log2(MAS5+1)", "CD8A  log2(MAS5+1)", "#2c5f8a",
    )
    scatter(
        axes[1], cldn4.values, gene_expr.loc["CD274"].values,
        f"GSE3141 NSCLC  n={int(rp.n)}\nCLDN4 vs CD274  ρ={rp.rho:.3f}  p={fmt_p(rp.p)}",
        "CLDN4  log2(MAS5+1)", "CD274  log2(MAS5+1)", "#6b3a2a",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_cd274.png", dpi=160)
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_cd274.pdf")
    plt.close(fig)

    forest = stats_df[
        (stats_df.gene_x == "CLDN4")
        & (stats_df.gene_y.isin(PRIMARY))
        & (stats_df.subset.isin(["all_tumors", "LUAD", "LUSC"]))
    ].copy()
    forest["label"] = forest["subset"] + "  " + forest["gene_y"]
    forest = forest.sort_values(["gene_y", "subset"])
    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    y_pos = np.arange(len(forest))
    ax.axvline(0, color="0.5", lw=0.8)
    ax.errorbar(
        forest["rho"], y_pos,
        xerr=[forest["rho"] - forest["ci_low"], forest["ci_high"] - forest["rho"]],
        fmt="o", color="#2c5f8a", ecolor="#2c5f8a", capsize=2, ms=5,
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(forest["label"])
    ax.set_xlabel("Spearman ρ vs CLDN4 (bootstrap 95% CI)")
    ax.set_title("GSE3141  CLDN4 vs CD8A / CD274  (honest n on each row)")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cldn4_forest.png", dpi=160)
    fig.savefig(FIGURES / "fig2_cldn4_forest.pdf")
    plt.close(fig)

    def pick(gene_y, subset, gene_x="CLDN4"):
        hit = stats_df[(stats_df.gene_x == gene_x) & (stats_df.gene_y == gene_y) & (stats_df.subset == subset)]
        return hit.iloc[0].to_dict() if len(hit) else {}

    one_row = pd.DataFrame([{
        "dataset": "GSE3141",
        "platform": "GPL570 U133 Plus 2.0 MAS5",
        "histology": f"NSCLC mixed (LUAD n={n_luad}; LUSC n={n_lusc})",
        "n_arrays": n_array,
        "n_CLDN4_CD8A": n_pair_cd8,
        "n_CLDN4_CD274": n_pair_pdl1,
        "CLDN4_probe": probe_audit.loc[probe_audit["gene"] == "CLDN4", "probe"].iloc[0],
        "CD8A_probe": probe_audit.loc[probe_audit["gene"] == "CD8A", "probe"].iloc[0],
        "CD274_probe": probe_audit.loc[probe_audit["gene"] == "CD274", "probe"].iloc[0],
        "CLDN4_CD8A_rho": pick("CD8A", "all_tumors")["rho"],
        "CLDN4_CD8A_p": pick("CD8A", "all_tumors")["p"],
        "CLDN4_CD8A_rho_adj": pick("CD8A", "all_tumors")["rho_adj_epithelial"],
        "CLDN4_CD8A_p_adj": pick("CD8A", "all_tumors")["p_adj_epithelial"],
        "CLDN4_CD8A_verdict": pick("CD8A", "all_tumors")["verdict"],
        "CLDN4_CD274_rho": pick("CD274", "all_tumors")["rho"],
        "CLDN4_CD274_p": pick("CD274", "all_tumors")["p"],
        "CLDN4_CD274_rho_adj": pick("CD274", "all_tumors")["rho_adj_epithelial"],
        "CLDN4_CD274_p_adj": pick("CD274", "all_tumors")["p_adj_epithelial"],
        "CLDN4_CD274_verdict": pick("CD274", "all_tumors")["verdict"],
        "ICI": "not deposited",
    }])
    one_row.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    summary = {
        "dataset": "GSE3141",
        "pmid": "16273092",
        "title": series.get("title", "Lung Cancer Dataset"),
        "platform": "GPL570 Affymetrix HG-U133 Plus 2.0",
        "processing": "MAS5 as deposited (linear signal). Spearman is rank-based.",
        "n_arrays": n_array,
        "n_probes": int(probe_expr.shape[0]),
        "n_missing_mas5": n_missing,
        "n_genes_collapsed": int(gene_expr.shape[0]),
        "n_LUAD": n_luad,
        "n_LUSC": n_lusc,
        "n_OS_numeric": n_surv_num,
        "n_OS_event_labeled": int(meta["os_event"].notna().sum()),
        "surv_non_numeric_gsm": meta.loc[~meta["surv_numeric"], "geo_accession"].tolist(),
        "cldn4_probe_collapse": probe_audit.loc[probe_audit["gene"] == "CLDN4", "probe"].iloc[0],
        "cd8a_probe_collapse": probe_audit.loc[probe_audit["gene"] == "CD8A", "probe"].iloc[0],
        "cd274_probe_collapse": probe_audit.loc[probe_audit["gene"] == "CD274", "probe"].iloc[0],
        "epithelial_genes_used": epi_genes,
        "primary": {
            "CLDN4_vs_CD8A_all": pick("CD8A", "all_tumors"),
            "CLDN4_vs_CD274_all": pick("CD274", "all_tumors"),
            "CLDN4_vs_CD8A_LUAD": pick("CD8A", "LUAD"),
            "CLDN4_vs_CD274_LUAD": pick("CD274", "LUAD"),
            "CLDN4_vs_CD8A_LUSC": pick("CD8A", "LUSC"),
            "CLDN4_vs_CD274_LUSC": pick("CD274", "LUSC"),
        },
        "highlow": hl.to_dict(orient="records"),
        "missing_public_labels": ["stage", "ICI", "tumor_percent", "one OS time encoded as >72"],
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    print("n_arrays", n_array, "LUAD", n_luad, "LUSC", n_lusc, "missing", n_missing)
    print(stats_df[stats_df.role.isin(["primary", "companion"])][
        ["gene_x", "gene_y", "subset", "n", "rho", "p", "rho_adj_epithelial", "p_adj_epithelial", "verdict"]
    ].to_string(index=False))
    print(hl.to_string(index=False))
    print("[done]", TABLES)


if __name__ == "__main__":
    main()
