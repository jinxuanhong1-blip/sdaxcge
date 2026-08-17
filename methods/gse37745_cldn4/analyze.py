#!/usr/bin/env python3
"""GSE37745 NSCLC Affy U133 Plus 2.0: CLDN4 vs CD8A / CD274.

Additive public bulk slice on the Botling Uppsala resected NSCLC series
(Botling et al., Clin Cancer Res 2013, PMID 23032747; GEO GSE37745).
Primary pairs are CLDN4–CD8A and CLDN4–CD274. Unit is the array.
One array per titled patient. No ICI arm.

Downloads stay under $GSE37745_CLDN4_DATA (default /tmp/gse37745_cldn4)
and are not committed.
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
DATA = Path(os.environ.get("GSE37745_CLDN4_DATA", "/tmp/gse37745_cldn4"))
SEED = 20260817
N_BOOT = 2000

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE37nnn/GSE37745/"
    "matrix/GSE37745_series_matrix.txt.gz"
)
GPL_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"

# Primary partners. TACSTD2 is companion only.
PRIMARY = ["CD8A", "CD274"]
COMPANION = ["TACSTD2"]

NAMED_PROBES = {
    "CLDN4": "201428_at",
    "CD8A": "205758_at",
    "CD274": "223834_at",
    "TACSTD2": "202286_s_at",
}

EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]

# Same HOLDS rule as the GSE4573 / PR 229 slices: n>=40, residual ρ<0, p<0.05.
HOLDS_N = 40


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "gse37745-cldn4/1.0"})
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
    char_rows: list[list[str]] = []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Series_"):
                k = line.split("\t", 1)[0][8:]
                v = line.split("\t", 1)[1].strip().strip('"') if "\t" in line else ""
                series.setdefault(k, []).append(v)
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                char_rows.append(vals)
            elif line.startswith("!Sample_"):
                key = line.split("\t", 1)[0][8:]
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                sample_fields[key] = vals
    meta = pd.DataFrame(sample_fields)
    if "geo_accession" not in meta.columns:
        raise KeyError(f"no geo_accession in sample fields: {list(meta.columns)}")
    meta.index = meta["geo_accession"].astype(str)
    meta.index.name = "gsm"
    for row in char_rows:
        keys = {v.split(": ", 1)[0] if ": " in v else "raw" for v in row}
        if len(keys) != 1:
            continue
        key = next(iter(keys))
        col = key.replace(" ", "_")
        meta[col] = [v.split(": ", 1)[1] if ": " in v else v for v in row]
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


def pair_row(cldn4, y, epithelial, gene_y: str, subset: str, n_arrays: int) -> dict:
    crude = spearman_ci(cldn4.values, y.values)
    adj = (
        partial_spearman(cldn4.values, y.values, epithelial.values)
        if epithelial is not None
        else dict(n=crude["n"], rho_adj=np.nan, p_adj=np.nan)
    )
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


def highlow(cldn4: pd.Series, y: pd.Series, endpoint: str, subset: str) -> dict:
    q = pd.qcut(cldn4, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    a = y[q == "Q4"]
    b = y[q == "Q1"]
    if a.size < 5 or b.size < 5:
        return {
            "subset": subset,
            "endpoint": endpoint,
            "n_Q4": int(a.size),
            "n_Q1": int(b.size),
            "median_Q4": float(np.median(a)) if a.size else np.nan,
            "median_Q1": float(np.median(b)) if b.size else np.nan,
            "U": np.nan,
            "p": np.nan,
            "rank_biserial_Q4_minus_Q1": np.nan,
        }
    U, p_mw = stats.mannwhitneyu(a.values, b.values, alternative="two-sided")
    rbc = 2 * U / (a.size * b.size) - 1
    return {
        "subset": subset,
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

    matrix_path = dl(DATA / "GSE37745_series_matrix.txt.gz", GEO_MATRIX)
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
    meta = meta.loc[probe_expr.columns].copy()

    meta["patient_id"] = meta["title"].str.extract(r"Patient\s+(\d+)", expand=False)
    meta["histology"] = meta["histology"].astype(str).str.strip().str.lower()
    hist_map = {"adeno": "LUAD", "squamous": "LUSC", "large": "LCC"}
    meta["histology_std"] = meta["histology"].map(hist_map).fillna(meta["histology"])

    with gzip.open(annot_path, "rt", errors="replace") as fh:
        annot = parse_soft_table(fh.read())
    annot["symbol"] = annot["Gene symbol"].map(first_symbol)
    id2gene = annot.drop_duplicates("ID").set_index("ID")["symbol"]

    gene_expr, probe_audit = collapse_maxmean(probe_expr, id2gene)

    confirm_rows = []
    for gene, named in NAMED_PROBES.items():
        hits = annot[annot["symbol"] == gene][["ID", "Gene symbol", "Gene title", "Gene ID"]]
        for _, r in hits.iterrows():
            chosen = (
                gene in gene_expr.index
                and probe_audit.loc[probe_audit["gene"] == gene, "probe"].iloc[0] == r["ID"]
            )
            confirm_rows.append({
                "gene": gene,
                "named_probe": named,
                "probe": r["ID"],
                "gpl570_symbol": r["Gene symbol"],
                "entrez": r["Gene ID"],
                "title": r["Gene title"],
                "in_matrix": r["ID"] in probe_expr.index,
                "is_named": r["ID"] == named,
                "chosen_for_collapse": chosen,
            })
    confirm = pd.DataFrame(confirm_rows)
    confirm.to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)

    # Probe-level sensitivity: named probe vs max-mean collapse (CD274 has two Plus-2 probes).
    probe_sens = []
    for gene, named in NAMED_PROBES.items():
        probes = confirm.loc[confirm.gene == gene, "probe"].tolist()
        for pr in probes:
            if pr not in probe_expr.index:
                continue
            probe_sens.append({
                "gene": gene,
                "probe": pr,
                "is_named": pr == named,
                "is_collapse": bool(
                    gene in gene_expr.index
                    and probe_audit.loc[probe_audit["gene"] == gene, "probe"].iloc[0] == pr
                ),
                "mean_rma": float(probe_expr.loc[pr].mean()),
            })
    pd.DataFrame(probe_sens).to_csv(TABLES / "probe_means.tsv", sep="\t", index=False)

    epi_present = [g for g in EPITHELIAL if g in gene_expr.index]
    if len(epi_present) < 3:
        raise RuntimeError(f"epithelial genes missing: present={epi_present}")
    z = gene_expr.loc[epi_present]
    mu = z.mean(axis=1)
    sd = z.std(axis=1, ddof=1).replace(0, np.nan)
    epithelial = z.sub(mu, axis=0).div(sd, axis=0).mean(axis=0)

    n_array = int(probe_expr.shape[1])
    n_na = int(probe_expr.isna().sum().sum())
    required = ["CLDN4", "CD8A", "CD274"]
    missing_genes = [g for g in required if g not in gene_expr.index]
    if missing_genes:
        raise RuntimeError(f"required genes absent after collapse: {missing_genes}")

    sample = meta[[
        "geo_accession", "title", "patient_id", "histology", "histology_std",
        "tumor_stage", "dead", "days_to_determined_death_status",
        "gender", "age", "adjuvant_treatment", "recurrence",
        "days_to_recurrence_/_to_last_visit",
        "performance_status_corresponding_to_who_criteria",
        "source_name_ch1", "data_processing",
    ]].copy()
    sample["cldn4_rma"] = gene_expr.loc["CLDN4"].reindex(sample.index).values
    sample["cd8a_rma"] = gene_expr.loc["CD8A"].reindex(sample.index).values
    sample["cd274_rma"] = gene_expr.loc["CD274"].reindex(sample.index).values
    sample["tacstd2_rma"] = (
        gene_expr.loc["TACSTD2"].reindex(sample.index).values
        if "TACSTD2" in gene_expr.index else np.nan
    )
    sample["epithelial_meanz"] = epithelial.reindex(sample.index).values
    sample.to_csv(TABLES / "sample_annotation.tsv", sep="\t")

    hist_counts = sample["histology"].value_counts().to_dict()
    stage_counts = sample["tumor_stage"].value_counts().to_dict()
    dead_counts = sample["dead"].value_counts().to_dict()
    rec_known = int((sample["recurrence"].isin(["yes", "no"])).sum())
    adj_known = int((sample["adjuvant_treatment"].isin(["yes", "no"])).sum())

    labels = pd.DataFrame([
        {"field": "arrays_in_matrix", "public": "yes", "n": n_array,
         "note": f"{probe_expr.shape[0]} probes x {n_array} GSM; missing values={n_na}"},
        {"field": "unique_GSM", "public": "yes", "n": int(meta.index.nunique()),
         "note": "all unique"},
        {"field": "unique_title", "public": "yes", "n": int(meta["title"].nunique()),
         "note": "Patient N, sex, histology; all unique"},
        {"field": "unique_patient_id", "public": "yes", "n": int(meta["patient_id"].nunique()),
         "note": "parsed from title; 1 array / patient; no duplicate pair"},
        {"field": "GEO_overall_design_patients", "public": "yes (text)", "n": 196,
         "note": "Series: '196 consecutive NSCLC patients, operated between 1995 and 2005'"},
        {"field": "histology_adeno", "public": "yes", "n": int(hist_counts.get("adeno", 0)),
         "note": "GEO characteristic histology: adeno"},
        {"field": "histology_squamous", "public": "yes", "n": int(hist_counts.get("squamous", 0)),
         "note": "GEO characteristic histology: squamous"},
        {"field": "histology_large", "public": "yes", "n": int(hist_counts.get("large", 0)),
         "note": "GEO characteristic histology: large; kept in NSCLC pool, reported separately"},
        {"field": "tumor_stage", "public": "yes", "n": n_array,
         "note": "counts=" + ",".join(f"{k}:{v}" for k, v in sorted(stage_counts.items()))},
        {"field": "OS_dead_and_days", "public": "yes", "n": n_array,
         "note": f"dead={dead_counts}; days finite for all {n_array}. Not the claim of this slice."},
        {"field": "recurrence_known", "public": "partial", "n": rec_known,
         "note": "100 of 196 are 'not known'; not used"},
        {"field": "adjuvant_known", "public": "partial", "n": adj_known,
         "note": "96 of 196 are 'not known'; surgical 1995–2005 cohort"},
        {"field": "ICI / PD-1 treatment", "public": "no", "n": 0,
         "note": "resected 1995–2005 Uppsala atlas; no ICI arm"},
        {"field": "tumor_percent / purity", "public": "no", "n": 0,
         "note": "epithelial RNA mean-z is the only public purity proxy"},
        {"field": "CLDN4 finite", "public": "yes",
         "n": int(sample["cldn4_rma"].notna().sum()),
         "note": "named probe 201428_at; RMA as deposited"},
        {"field": "CD8A finite", "public": "yes",
         "n": int(sample["cd8a_rma"].notna().sum()),
         "note": "named probe 205758_at"},
        {"field": "CD274 finite", "public": "yes",
         "n": int(sample["cd274_rma"].notna().sum()),
         "note": "named probe 223834_at (Plus 2.0; 227458_at also present)"},
        {"field": "pairwise CLDN4+CD8A", "public": "yes",
         "n": int((sample["cldn4_rma"].notna() & sample["cd8a_rma"].notna()).sum()),
         "note": "primary n, all NSCLC"},
        {"field": "pairwise CLDN4+CD274", "public": "yes",
         "n": int((sample["cldn4_rma"].notna() & sample["cd274_rma"].notna()).sum()),
         "note": "primary n, all NSCLC"},
    ])
    labels.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    if n_array != 196:
        raise AssertionError(f"expected 196 arrays, got {n_array}")
    if int(meta["patient_id"].nunique()) != 196:
        raise AssertionError("patient IDs are not 1:1 with arrays")

    subsets = {
        "NSCLC_all": sample.index,
        "LUAD": sample.index[sample["histology"] == "adeno"],
        "LUSC": sample.index[sample["histology"] == "squamous"],
        "LCC": sample.index[sample["histology"] == "large"],
    }

    gene_rows = []
    hl_rows = []
    partners = PRIMARY + [g for g in COMPANION if g in gene_expr.index]
    partners = partners + ["epithelial_meanz"]
    for subset, idx in subsets.items():
        cldn4 = gene_expr.loc["CLDN4", idx]
        epi = epithelial.reindex(idx)
        n_sub = int(len(idx))
        for g in partners:
            if g == "epithelial_meanz":
                y = epi
                rec = pair_row(cldn4, y, None, g, subset, n_sub)
                rec["rho_adj_epithelial"] = np.nan
                rec["p_adj_epithelial"] = np.nan
                rec["verdict"] = "COMPANION"
                gene_rows.append(rec)
                continue
            if g not in gene_expr.index:
                gene_rows.append({
                    "subset": subset, "gene_x": "CLDN4", "gene_y": g,
                    "n": 0, "n_arrays": n_sub, "present": False,
                    "rho": np.nan, "p": np.nan, "ci_low": np.nan, "ci_high": np.nan,
                    "rho_adj_epithelial": np.nan, "p_adj_epithelial": np.nan,
                    "verdict": "ABSENT",
                })
                continue
            y = gene_expr.loc[g, idx]
            gene_rows.append(pair_row(cldn4, y, epi, g, subset, n_sub))
            if g in PRIMARY:
                hl_rows.append(highlow(cldn4, y, g, subset))
    genes_df = pd.DataFrame(gene_rows)
    genes_df.to_csv(TABLES / "spearman_cldn4_vs_cd8a_cd274.tsv", sep="\t", index=False)
    pd.DataFrame(hl_rows).to_csv(TABLES / "highlow_cldn4_cd8a_cd274.tsv", sep="\t", index=False)

    tac_rows = []
    if "TACSTD2" in gene_expr.index:
        for subset, idx in subsets.items():
            x = gene_expr.loc["TACSTD2", idx]
            epi = epithelial.reindex(idx)
            for g in PRIMARY:
                y = gene_expr.loc[g, idx]
                crude = spearman_ci(x.values, y.values)
                adj = partial_spearman(x.values, y.values, epi.values)
                tac_rows.append({
                    "subset": subset,
                    "gene_x": "TACSTD2",
                    "gene_y": g,
                    **crude,
                    "rho_adj_epithelial": adj["rho_adj"],
                    "p_adj_epithelial": adj["p_adj"],
                    "verdict": verdict(adj["rho_adj"], adj["p_adj"], adj["n"]),
                    "note": "companion only; not the claim",
                })
    tac_df = pd.DataFrame(tac_rows)
    tac_df.to_csv(TABLES / "spearman_tacstd2_companion.tsv", sep="\t", index=False)

    # Named-probe vs collapse for the two primary pairs, by histology.
    named_rows = []
    for subset, idx in subsets.items():
        cldn4_named = probe_expr.loc[NAMED_PROBES["CLDN4"], idx]
        epi_sub = epithelial.reindex(idx)
        for gene in PRIMARY:
            named = NAMED_PROBES[gene]
            y_named = probe_expr.loc[named, idx]
            y_coll = gene_expr.loc[gene, idx]
            variants = [
                (f"{gene}_named_{named}", y_named),
                (f"{gene}_collapse", y_coll),
            ]
            extra = [p for p in confirm.loc[confirm.gene == gene, "probe"] if p != named]
            for pr in extra:
                variants.append((f"{gene}_probe_{pr}", probe_expr.loc[pr, idx]))
            seen = set()
            for label, y in variants:
                if label in seen:
                    continue
                seen.add(label)
                crude = spearman_ci(cldn4_named.values, y.values)
                adj = partial_spearman(cldn4_named.values, y.values, epi_sub.values)
                named_rows.append({
                    "subset": subset,
                    "cldn4": NAMED_PROBES["CLDN4"],
                    "partner": label,
                    **crude,
                    "rho_adj_epithelial": adj["rho_adj"],
                    "p_adj_epithelial": adj["p_adj"],
                    "verdict": verdict(adj["rho_adj"], adj["p_adj"], adj["n"]),
                })
    pd.DataFrame(named_rows).to_csv(TABLES / "spearman_named_probe_sensitivity.tsv", sep="\t", index=False)

    # Figures: all-NSCLC primary pairs + histology forest.
    cldn4_all = gene_expr.loc["CLDN4"]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
    for ax, gene, color in zip(axes, PRIMARY, ["#2c5f8a", "#6b3a2a"]):
        row = genes_df[(genes_df.subset == "NSCLC_all") & (genes_df.gene_y == gene)].iloc[0]
        ax.scatter(cldn4_all.values, gene_expr.loc[gene].values, s=16, alpha=0.7, c=color, edgecolors="none")
        ax.set_xlabel("CLDN4  RMA")
        ax.set_ylabel(f"{gene}  RMA")
        ax.set_title(
            f"GSE37745 NSCLC  n={int(row.n)}\n"
            f"CLDN4 vs {gene}  ρ={row.rho:.3f}  p={fmt_p(row.p)}"
        )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_cd274.png", dpi=160)
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_cd274.pdf")
    plt.close(fig)

    plot = genes_df[genes_df.gene_y.isin(PRIMARY)].copy()
    plot["label"] = plot["subset"] + "  " + plot["gene_y"]
    order = []
    for sub in ["NSCLC_all", "LUAD", "LUSC", "LCC"]:
        for g in PRIMARY:
            order.append(f"{sub}  {g}")
    plot = plot.set_index("label").loc[order].reset_index()
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    y_pos = np.arange(len(plot))
    ax.axvline(0, color="0.5", lw=0.8)
    ax.errorbar(
        plot["rho"], y_pos,
        xerr=[plot["rho"] - plot["ci_low"], plot["ci_high"] - plot["rho"]],
        fmt="o", color="#2c5f8a", ecolor="#2c5f8a", capsize=2, ms=5,
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot["label"])
    ax.set_xlabel("Spearman ρ vs CLDN4 (bootstrap 95% CI)")
    ax.set_title("GSE37745  CLDN4 vs CD8A / CD274  by histology")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cldn4_forest_by_histology.png", dpi=160)
    fig.savefig(FIGURES / "fig2_cldn4_forest_by_histology.pdf")
    plt.close(fig)

    def pick(subset, gene):
        return genes_df[(genes_df.subset == subset) & (genes_df.gene_y == gene)].iloc[0]

    summary = {
        "dataset": "GSE37745",
        "pmid": "23032747",
        "platform": "GPL570 Affymetrix HG-U133 Plus 2.0",
        "processing": "Bioconductor affy RMA as deposited in the series matrix. Spearman is rank-based.",
        "histology": "NSCLC mixed: adeno 106, squamous 66, large 24",
        "n_arrays": n_array,
        "n_patients": int(meta["patient_id"].nunique()),
        "n_probes": int(probe_expr.shape[0]),
        "n_genes_collapsed": int(gene_expr.shape[0]),
        "n_missing_matrix_values": n_na,
        "cldn4_probe": NAMED_PROBES["CLDN4"],
        "cd8a_probe": NAMED_PROBES["CD8A"],
        "cd274_probe": NAMED_PROBES["CD274"],
        "cd274_other_probes": confirm.loc[confirm.gene == "CD274", "probe"].tolist(),
        "epithelial_genes_used": epi_present,
        "primary": {
            "NSCLC_all": {
                "CD8A": pick("NSCLC_all", "CD8A").to_dict(),
                "CD274": pick("NSCLC_all", "CD274").to_dict(),
            },
            "LUAD": {
                "CD8A": pick("LUAD", "CD8A").to_dict(),
                "CD274": pick("LUAD", "CD274").to_dict(),
            },
            "LUSC": {
                "CD8A": pick("LUSC", "CD8A").to_dict(),
                "CD274": pick("LUSC", "CD274").to_dict(),
            },
            "LCC": {
                "CD8A": pick("LCC", "CD8A").to_dict(),
                "CD274": pick("LCC", "CD274").to_dict(),
            },
        },
        "missing_public_labels": ["ICI", "tumor_percent", "pathologist CD8/PD-L1 IHC"],
        "os_deposited": True,
        "os_used_in_this_slice": False,
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    print("n_arrays", n_array, "patients", meta["patient_id"].nunique(), "genes", gene_expr.shape[0])
    print(genes_df[genes_df.gene_y.isin(PRIMARY)].to_string(index=False))
    print("[done]", TABLES)


if __name__ == "__main__":
    main()
