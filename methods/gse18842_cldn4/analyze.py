#!/usr/bin/env python3
"""GSE18842 NSCLC Affy U133 Plus 2.0: CLDN4 vs CD8A / CD274 / ImmuneScore.

Additive CLDN4-only public bulk slice. No dual-high. Unit is the array.
Tumor arrays only; adjacent-normal / control arrays are dropped.
GEO does not label LUAD vs LUSC on this series, so no histology split
is invented. No ICI arm. No slide was re-scored.

Downloads stay under $GSE18842_CLDN4_DATA (default /tmp/gse18842_cldn4)
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
DATA = Path(os.environ.get("GSE18842_CLDN4_DATA", "/tmp/gse18842_cldn4"))
SIG = HERE / "estimate_yoshihara_2013.tsv"
SEED = 20260817
N_BOOT = 2000

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18842/"
    "matrix/GSE18842_series_matrix.txt.gz"
)
GPL_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"

PRIMARY = ["CD8A", "CD274"]
EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]

NAMED_PROBES = {
    "CLDN4": "201428_at",
    "CD8A": "205758_at",
    "CD274": "227458_at",
}
CD274_PROBES = ["223834_at", "227458_at"]
CLDN4_PROBES = ["201428_at", "1569421_at"]

# Same HOLDS rule as PR 229 / GSE4573 CLDN4 slice: n>=40, adj ρ<0, p<0.05.
HOLDS_N = 40

EST_A = 0.6049872018
EST_B = 0.0001467884


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 10_000:
        print(f"[cached] {dest} ({dest.stat().st_size} bytes)", flush=True)
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "gse18842-cldn4/1.0"})
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "wb") as out:
        out.write(r.read())
    tmp.replace(dest)
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
    char_n = 0
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
                if key == "characteristics_ch1":
                    char_n += 1
                    key = f"characteristics_ch1_{char_n}"
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


def fmt_rho(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:+.3f}"


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


def ssgsea(expr: pd.DataFrame, genes: list[str], tau: float = 0.25) -> pd.Series:
    present = [g for g in genes if g in expr.index]
    if len(present) < 10:
        return pd.Series(np.nan, index=expr.columns)
    n_genes = expr.shape[0]
    ranked = expr.rank(axis=0, method="average", ascending=True) * (10000.0 / n_genes)
    gene_set = set(present)
    scores = {}
    for sample in expr.columns:
        m = ranked[sample]
        order = m.sort_values(ascending=False).index
        m_ord = m.loc[order].to_numpy(float)
        hits = np.fromiter((g in gene_set for g in order), dtype=bool, count=len(order))
        w = np.abs(m_ord) ** tau
        w_hit = np.where(hits, w, 0.0)
        nhit = float(w_hit.sum())
        nmiss = float((~hits).sum())
        if nhit <= 0 or nmiss <= 0:
            scores[sample] = np.nan
            continue
        p_hit = np.cumsum(w_hit) / nhit
        p_miss = np.cumsum((~hits).astype(float)) / nmiss
        scores[sample] = float(np.sum(p_hit - p_miss))
    return pd.Series(scores)


def estimate_scores(gene_expr: pd.DataFrame, stromal: list[str], immune: list[str]) -> pd.DataFrame:
    strom = ssgsea(gene_expr, stromal)
    imm = ssgsea(gene_expr, immune)
    est = strom + imm
    pur = np.cos(EST_A + EST_B * est.to_numpy(float))
    return pd.DataFrame(
        {
            "ESTIMATE_StromalScore": strom,
            "ESTIMATE_ImmuneScore": imm,
            "ESTIMATE_Score": est,
            "ESTIMATE_TumorPurity": pur,
        },
        index=gene_expr.columns,
    )


def parse_pair_id(title: str) -> tuple[str, str]:
    t = str(title).strip()
    m = re.match(r"(Tumor|Control)\s+(.+)$", t)
    if not m:
        return "", ""
    return m.group(1), m.group(2)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    matrix_path = dl(DATA / "GSE18842_series_matrix.txt.gz", GEO_MATRIX)
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

    # Sample type: GEO characteristic "sample type: tumor|control"
    sample_type = meta.get("characteristics_ch1_2", pd.Series("", index=meta.index)).astype(str)
    sample_type = sample_type.str.replace(r"^sample type:\s*", "", regex=True).str.strip().str.lower()
    if sample_type.eq("").all():
        src = meta["source_name_ch1"].astype(str).str.lower()
        sample_type = np.where(src.str.contains("tumor"), "tumor",
                               np.where(src.str.contains("control"), "control", "unknown"))
        sample_type = pd.Series(sample_type, index=meta.index)
    meta["sample_type"] = sample_type
    meta["is_tumor"] = meta["sample_type"].eq("tumor")
    meta["is_control"] = meta["sample_type"].eq("control")

    pair_parsed = meta["title"].map(parse_pair_id)
    meta["title_class"] = [a for a, _ in pair_parsed]
    meta["pair_id"] = [b for _, b in pair_parsed]
    tumor_ids = set(meta.loc[meta["is_tumor"], "pair_id"])
    ctrl_ids = set(meta.loc[meta["is_control"], "pair_id"])
    meta["paired_by_title"] = meta["pair_id"].isin(tumor_ids & ctrl_ids) & meta["pair_id"].ne("")

    # No LUAD/LUSC characteristic is deposited.
    hist_fields = [c for c in meta.columns if re.search(r"histolog|cell type|disease|subtype", c, re.I)]
    meta["histology_geo"] = "not_deposited"

    with gzip.open(annot_path, "rt", errors="replace") as fh:
        annot = parse_soft_table(fh.read())
    symbol_col = "Gene symbol" if "Gene symbol" in annot.columns else (
        next(c for c in annot.columns if "symbol" in c.lower())
    )
    annot["symbol"] = annot[symbol_col].map(first_symbol)
    id2gene = annot.drop_duplicates("ID").set_index("ID")["symbol"]

    gene_expr, probe_audit = collapse_maxmean(probe_expr, id2gene)

    probe_rows = []
    for gene, pid in NAMED_PROBES.items():
        if pid not in probe_expr.index:
            probe_rows.append({"gene": gene, "named_probe": pid, "present": False, "gpl570_symbol": "", "used": "ABSENT"})
            continue
        probe_rows.append({
            "gene": gene,
            "named_probe": pid,
            "present": True,
            "gpl570_symbol": str(id2gene.get(pid, "")),
            "used": "named",
            "n_probes_for_symbol": int((probe_audit["gene"] == gene).sum()) if gene in set(probe_audit["gene"]) else 0,
            "collapsed_probe": (
                probe_audit.loc[probe_audit["gene"] == gene, "probe"].iloc[0]
                if gene in set(probe_audit["gene"]) else ""
            ),
        })
    for p in CD274_PROBES:
        if p == NAMED_PROBES["CD274"]:
            continue
        if p in probe_expr.index:
            probe_rows.append({
                "gene": "CD274",
                "named_probe": p,
                "present": True,
                "gpl570_symbol": str(id2gene.get(p, "")),
                "used": "companion_probe",
                "n_probes_for_symbol": int((probe_audit["gene"] == "CD274").sum()),
                "collapsed_probe": probe_audit.loc[probe_audit["gene"] == "CD274", "probe"].iloc[0],
            })
    pd.DataFrame(probe_rows).to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)

    # Force named probes for the three headline genes.
    for gene, pid in NAMED_PROBES.items():
        if pid in probe_expr.index:
            gene_expr.loc[gene] = probe_expr.loc[pid].astype(float)

    n_array = int(probe_expr.shape[1])
    n_missing = int(probe_expr.isna().sum().sum())
    tumor = meta.index[meta["is_tumor"]]
    control = meta.index[meta["is_control"]]
    n_tumor = int(len(tumor))
    n_control = int(len(control))
    n_paired = int(meta.loc[meta["is_tumor"], "paired_by_title"].sum())
    n_tumor_unpaired = n_tumor - n_paired
    n_control_unpaired = int((~meta.loc[meta["is_control"], "paired_by_title"]).sum())

    if n_array != 91:
        print(f"[warn] expected 91 arrays, got {n_array}", flush=True)
    if n_tumor != 46 or n_control != 45:
        print(f"[warn] expected 46 tumor / 45 control, got {n_tumor}/{n_control}", flush=True)
    for g in ["CLDN4", "CD8A", "CD274"]:
        if g not in gene_expr.index:
            raise SystemExit(f"required gene absent after collapse: {g}")

    tumor_expr = gene_expr.loc[:, tumor]
    epi_genes = [g for g in EPITHELIAL if g in gene_expr.index]
    if not epi_genes:
        raise SystemExit("no epithelial genes present")
    epi_z = tumor_expr.loc[epi_genes].apply(lambda s: (s - s.mean()) / (s.std(ddof=0) or 1.0), axis=1)
    epithelial = epi_z.mean(axis=0)

    sig = pd.read_csv(SIG, sep="\t")
    stromal_genes = sig.loc[sig["set"].str.contains("Stromal", case=False), "hugo"].astype(str).tolist()
    immune_genes = sig.loc[sig["set"].str.contains("Immune", case=False), "hugo"].astype(str).tolist()
    n_strom_present = sum(g in gene_expr.index for g in stromal_genes)
    n_imm_present = sum(g in gene_expr.index for g in immune_genes)
    est = estimate_scores(tumor_expr, stromal_genes, immune_genes)
    immune_score = est["ESTIMATE_ImmuneScore"]
    purity = est["ESTIMATE_TumorPurity"]
    n_purity_01 = int(((purity >= 0) & (purity <= 1)).sum())

    sample = meta.copy()
    sample["cldn4"] = gene_expr.loc["CLDN4"]
    sample["cd8a"] = gene_expr.loc["CD8A"]
    sample["cd274"] = gene_expr.loc["CD274"]
    sample["epithelial_mean_z"] = np.nan
    sample.loc[epithelial.index, "epithelial_mean_z"] = epithelial
    for col in est.columns:
        sample[col] = np.nan
        sample.loc[est.index, col] = est[col]
    sample.to_csv(TABLES / "sample_annotation.tsv", sep="\t")

    n_pair_cd8 = int(sample.loc[tumor, ["cldn4", "cd8a"]].notna().all(axis=1).sum())
    n_pair_pdl1 = int(sample.loc[tumor, ["cldn4", "cd274"]].notna().all(axis=1).sum())
    n_pair_imm = int(sample.loc[tumor, ["cldn4"]].notna().all(axis=1).sum() and immune_score.notna().sum())

    labels = pd.DataFrame([
        {"field": "arrays in series matrix", "public": "yes", "n": n_array,
         "note": "54,675 probes × 91 GSM; RMA as deposited"},
        {"field": "unique GSM / unique titles", "public": "yes", "n": int(meta.index.nunique()),
         "note": "Tumor N / Control N / NP titles; all unique"},
        {"field": "tumor arrays (sample type: tumor)", "public": "yes", "n": n_tumor,
         "note": "PRIMARY n; source_name Human Lung Tumor"},
        {"field": "control / adjacent-normal arrays", "public": "yes", "n": n_control,
         "note": "sample type: control; DROPPED from pairwise tests"},
        {"field": "title-paired tumor+control", "public": "yes", "n": n_paired,
         "note": "44 numeric pair IDs shared by Tumor and Control titles"},
        {"field": "unpaired tumors (title NP)", "public": "yes", "n": n_tumor_unpaired,
         "note": "Tumor NP_1, Tumor NP_2; GEO overall_design: all paired except three"},
        {"field": "unpaired controls (title NP)", "public": "yes", "n": n_control_unpaired,
         "note": "Control NP"},
        {"field": "LUAD / LUSC / histology characteristic", "public": "no", "n": 0,
         "note": "not a GEO characteristic; paper text mentions ADC and SCC but per-array labels are not deposited. No split invented."},
        {"field": "stage", "public": "no", "n": 0, "note": "paper discusses stage; not a GEO characteristic"},
        {"field": "OS / DFS time or event", "public": "no", "n": 0, "note": "not deposited"},
        {"field": "ICI / treatment", "public": "no", "n": 0, "note": "resected atlas, not an ICI series"},
        {"field": "tumor % / ABSOLUTE purity", "public": "no", "n": 0,
         "note": "public proxies used here: epithelial mean-z and ESTIMATE"},
        {"field": "CLDN4 finite on tumors (201428_at)", "public": "yes",
         "n": int(sample.loc[tumor, "cldn4"].notna().sum()),
         "note": "named Plus2 probe; official GPL570 = CLDN4"},
        {"field": "CD8A finite on tumors (205758_at)", "public": "yes",
         "n": int(sample.loc[tumor, "cd8a"].notna().sum()),
         "note": "named Plus2 probe; official GPL570 = CD8A"},
        {"field": "CD274 finite on tumors (227458_at)", "public": "yes",
         "n": int(sample.loc[tumor, "cd274"].notna().sum()),
         "note": "named Plus2 probe (higher-mean of two clean CD274 probes)"},
        {"field": "primary pairwise n (CLDN4 + CD8A tumors)", "public": "yes", "n": n_pair_cd8,
         "note": "this is the n used below"},
        {"field": "primary pairwise n (CLDN4 + CD274 tumors)", "public": "yes", "n": n_pair_pdl1,
         "note": "this is the n used below"},
        {"field": "ESTIMATE ImmuneScore on tumors", "public": "computed", "n": int(immune_score.notna().sum()),
         "note": f"Yoshihara Immune141 ssGSEA ({n_imm_present}/141); Stromal141 {n_strom_present}/141"},
        {"field": "ESTIMATE TumorPurity in [0,1]", "public": "computed", "n": n_purity_01,
         "note": "cosine transform; values outside [0,1] are not used as a 0–1 fraction"},
        {"field": "histology fields on series matrix", "public": "no", "n": len(hist_fields),
         "note": "none" if not hist_fields else ",".join(hist_fields)},
    ])
    labels.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    cldn4 = gene_expr.loc["CLDN4", tumor]
    rows = []
    idx = tumor
    epi = epithelial.reindex(idx)
    for g in PRIMARY:
        rows.append(
            pair_row(
                cldn4.values,
                gene_expr.loc[g, idx].values,
                epi.values,
                "CLDN4",
                g,
                "tumor",
                n_arrays=n_tumor,
                extra={"role": "primary"},
            )
        )
    rows.append(
        pair_row(
            cldn4.values,
            immune_score.reindex(idx).values,
            epi.values,
            "CLDN4",
            "ImmuneScore",
            "tumor",
            n_arrays=n_tumor,
            extra={"role": "primary"},
        )
    )
    rows.append(
        pair_row(
            cldn4.values,
            epi.values,
            epi.values,
            "CLDN4",
            "epithelial_mean_z",
            "tumor",
            n_arrays=n_tumor,
            extra={"role": "covariate"},
        )
    )
    rows[-1]["rho_adj_epithelial"] = np.nan
    rows[-1]["p_adj_epithelial"] = np.nan
    rows[-1]["verdict"] = "COMPANION"

    rows.append(
        pair_row(
            cldn4.values,
            purity.reindex(idx).values,
            epi.values,
            "CLDN4",
            "TumorPurity",
            "tumor",
            n_arrays=n_tumor,
            extra={"role": "covariate"},
        )
    )
    rows.append(
        pair_row(
            gene_expr.loc["CD8A", idx].values,
            immune_score.reindex(idx).values,
            epi.values,
            "CD8A",
            "ImmuneScore",
            "tumor",
            n_arrays=n_tumor,
            extra={"role": "positive_control"},
        )
    )
    rows.append(
        pair_row(
            gene_expr.loc["CD274", idx].values,
            gene_expr.loc["CD8A", idx].values,
            epi.values,
            "CD274",
            "CD8A",
            "tumor",
            n_arrays=n_tumor,
            extra={"role": "positive_control"},
        )
    )

    # Named-probe sensitivity (CD274 has two GPL570 probes).
    for p in CD274_PROBES:
        rows.append(
            pair_row(
                probe_expr.loc[NAMED_PROBES["CLDN4"], idx].values,
                probe_expr.loc[p, idx].values,
                epi.values,
                "CLDN4_201428_at",
                f"CD274_{p}",
                "tumor_named_probe",
                n_arrays=n_tumor,
                extra={"role": "probe_sensitivity"},
            )
        )
    rows.append(
        pair_row(
            probe_expr.loc[NAMED_PROBES["CLDN4"], idx].values,
            probe_expr.loc[NAMED_PROBES["CD8A"], idx].values,
            epi.values,
            "CLDN4_201428_at",
            "CD8A_205758_at",
            "tumor_named_probe",
            n_arrays=n_tumor,
            extra={"role": "probe_sensitivity"},
        )
    )

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(TABLES / "spearman_cldn4_vs_cd8a_cd274.tsv", sep="\t", index=False)

    hl = pd.DataFrame([
        highlow(cldn4.values, gene_expr.loc["CD8A", idx].values, "CD8A"),
        highlow(cldn4.values, gene_expr.loc["CD274", idx].values, "CD274"),
        highlow(cldn4.values, immune_score.reindex(idx).values, "ImmuneScore"),
    ])
    hl.to_csv(TABLES / "highlow_cldn4.tsv", sep="\t", index=False)

    cov = pd.DataFrame([
        {"set": f"gene:{g}", "n_listed": 1, "n_present": int(g in gene_expr.index),
         "probe_used": (probe_audit.loc[probe_audit["gene"] == g, "probe"].iloc[0]
                        if g in gene_expr.index else ""),
         "missing": "" if g in gene_expr.index else g}
        for g in ["CLDN4", "CD8A", "CD274"] + EPITHELIAL
    ] + [
        {"set": "ESTIMATE_Stromal141", "n_listed": len(stromal_genes), "n_present": n_strom_present,
         "probe_used": "ssGSEA first-symbol collapse", "missing": ""},
        {"set": "ESTIMATE_Immune141", "n_listed": len(immune_genes), "n_present": n_imm_present,
         "probe_used": "ssGSEA first-symbol collapse", "missing": ""},
        {"set": "epithelial_mean_z", "n_listed": len(EPITHELIAL), "n_present": len(epi_genes),
         "probe_used": ",".join(epi_genes), "missing": ",".join(g for g in EPITHELIAL if g not in epi_genes)},
    ])
    cov.to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)

    def scatter(ax, x, y, title, xlab, ylab, color):
        ax.scatter(np.asarray(x, float), np.asarray(y, float),
                   s=22, alpha=0.8, c=color, edgecolors="none")
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(title)

    def pick(gene_y, subset="tumor", gene_x="CLDN4"):
        hit = stats_df[(stats_df.gene_x == gene_x) & (stats_df.gene_y == gene_y) & (stats_df.subset == subset)]
        return hit.iloc[0].to_dict() if len(hit) else {}

    r8 = pick("CD8A")
    rp = pick("CD274")
    ri = pick("ImmuneScore")

    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.2))
    scatter(
        axes[0], cldn4.values, gene_expr.loc["CD8A", idx].values,
        f"GSE18842 tumor  n={int(r8['n'])}\nCLDN4 vs CD8A  ρ={r8['rho']:.3f}  p={fmt_p(r8['p'])}",
        "CLDN4  RMA (named 201428_at)", "CD8A  RMA (205758_at)", "#2c5f8a",
    )
    scatter(
        axes[1], cldn4.values, gene_expr.loc["CD274", idx].values,
        f"GSE18842 tumor  n={int(rp['n'])}\nCLDN4 vs CD274  ρ={rp['rho']:.3f}  p={fmt_p(rp['p'])}",
        "CLDN4  RMA (named 201428_at)", "CD274  RMA (227458_at)", "#6b3a2a",
    )
    scatter(
        axes[2], cldn4.values, immune_score.reindex(idx).values,
        f"GSE18842 tumor  n={int(ri['n'])}\nCLDN4 vs ImmuneScore  ρ={ri['rho']:.3f}  p={fmt_p(ri['p'])}",
        "CLDN4  RMA (named 201428_at)", "ESTIMATE ImmuneScore", "#2f6b4f",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_cd274_immunescore.png", dpi=160)
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_cd274_immunescore.pdf")
    plt.close(fig)

    forest = stats_df[
        (stats_df.gene_x == "CLDN4")
        & (stats_df.gene_y.isin(PRIMARY + ["ImmuneScore"]))
        & (stats_df.subset == "tumor")
    ].copy()
    forest["label"] = forest["gene_y"]
    # keep CD8A, CD274, ImmuneScore order
    order = {k: i for i, k in enumerate(["CD8A", "CD274", "ImmuneScore"])}
    forest = forest.assign(_ord=forest["gene_y"].map(order)).sort_values("_ord")
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    y_pos = np.arange(len(forest))
    ax.axvline(0, color="0.5", lw=0.8)
    ax.errorbar(
        forest["rho"], y_pos,
        xerr=[forest["rho"] - forest["ci_low"], forest["ci_high"] - forest["rho"]],
        fmt="o", color="#2c5f8a", ecolor="#2c5f8a", capsize=2, ms=6,
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(forest["label"])
    ax.set_xlabel("Spearman ρ vs CLDN4 (bootstrap 95% CI)")
    ax.set_title(f"GSE18842 tumor  CLDN4 vs CD8A / CD274 / ImmuneScore  n={n_tumor}")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cldn4_forest.png", dpi=160)
    fig.savefig(FIGURES / "fig2_cldn4_forest.pdf")
    plt.close(fig)

    one_row = pd.DataFrame([{
        "dataset": "GSE18842",
        "platform": "GPL570 U133 Plus 2.0 RMA",
        "histology": "NSCLC mixed; LUAD/LUSC not deposited on GEO",
        "n_arrays_matrix": n_array,
        "n_tumor": n_tumor,
        "n_control_dropped": n_control,
        "n_CLDN4_CD8A": n_pair_cd8,
        "n_CLDN4_CD274": n_pair_pdl1,
        "n_CLDN4_ImmuneScore": int(immune_score.notna().sum()),
        "CLDN4_probe": NAMED_PROBES["CLDN4"],
        "CD8A_probe": NAMED_PROBES["CD8A"],
        "CD274_probe": NAMED_PROBES["CD274"],
        "CLDN4_CD8A_rho": r8["rho"],
        "CLDN4_CD8A_p": r8["p"],
        "CLDN4_CD8A_rho_adj": r8["rho_adj_epithelial"],
        "CLDN4_CD8A_p_adj": r8["p_adj_epithelial"],
        "CLDN4_CD8A_verdict": r8["verdict"],
        "CLDN4_CD274_rho": rp["rho"],
        "CLDN4_CD274_p": rp["p"],
        "CLDN4_CD274_rho_adj": rp["rho_adj_epithelial"],
        "CLDN4_CD274_p_adj": rp["p_adj_epithelial"],
        "CLDN4_CD274_verdict": rp["verdict"],
        "CLDN4_ImmuneScore_rho": ri["rho"],
        "CLDN4_ImmuneScore_p": ri["p"],
        "CLDN4_ImmuneScore_rho_adj": ri["rho_adj_epithelial"],
        "CLDN4_ImmuneScore_p_adj": ri["p_adj_epithelial"],
        "CLDN4_ImmuneScore_verdict": ri["verdict"],
        "epithelial_genes": ",".join(epi_genes),
        "estimate_immune_coverage": f"{n_imm_present}/141",
        "ICI": "not deposited",
        "LUAD_LUSC_split": "not labeled on GEO; not invented",
    }])
    one_row.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    summary = {
        "dataset": "GSE18842",
        "pmid": "20878980",
        "title": series.get("title", ""),
        "overall_design": series.get("overall_design", ""),
        "platform": "GPL570 Affymetrix HG-U133 Plus 2.0",
        "processing": "RMA as deposited. Spearman is rank-based.",
        "n_arrays": n_array,
        "n_probes": int(probe_expr.shape[0]),
        "n_missing_rma": n_missing,
        "n_genes_collapsed": int(gene_expr.shape[0]),
        "n_tumor": n_tumor,
        "n_control_dropped": n_control,
        "n_title_paired_tumors": n_paired,
        "n_unpaired_tumors": n_tumor_unpaired,
        "histology_geo": "not_deposited",
        "cldn4_probe": NAMED_PROBES["CLDN4"],
        "cd8a_probe": NAMED_PROBES["CD8A"],
        "cd274_probe": NAMED_PROBES["CD274"],
        "epithelial_genes_used": epi_genes,
        "estimate_stromal_present": n_strom_present,
        "estimate_immune_present": n_imm_present,
        "n_purity_in_0_1": n_purity_01,
        "primary": {
            "CLDN4_vs_CD8A_tumor": r8,
            "CLDN4_vs_CD274_tumor": rp,
            "CLDN4_vs_ImmuneScore_tumor": ri,
            "CLDN4_vs_epithelial": pick("epithelial_mean_z"),
            "CD8A_vs_ImmuneScore": pick("ImmuneScore", gene_x="CD8A"),
            "CD274_vs_CD8A": pick("CD8A", gene_x="CD274"),
        },
        "highlow": hl.to_dict(orient="records"),
        "missing_public_labels": ["LUAD/LUSC", "stage", "OS", "ICI", "tumor_percent"],
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    print("n_arrays", n_array, "tumor", n_tumor, "control", n_control, "missing", n_missing)
    print("paired tumors", n_paired, "unpaired tumors", n_tumor_unpaired)
    print("epithelial", epi_genes, "Immune141", n_imm_present, "Stromal141", n_strom_present)
    print(stats_df[stats_df.role.isin(["primary", "positive_control", "covariate"])][
        ["gene_x", "gene_y", "subset", "n", "rho", "p", "rho_adj_epithelial", "p_adj_epithelial", "verdict"]
    ].to_string(index=False))
    print(hl.to_string(index=False))
    print("[done]", TABLES)


if __name__ == "__main__":
    main()
