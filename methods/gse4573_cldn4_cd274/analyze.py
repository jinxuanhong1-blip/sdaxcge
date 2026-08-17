#!/usr/bin/env python3
"""GSE4573 LUSC Affy U133A: CLDN4 vs CD274 and HLA-A/B/C.

Additive CD274 / classical MHC-I cut. CLDN4 vs CD8A is already NS in
PR 313 and is not the claim here. CD8A is used only as a residual
covariate (same rule as the ICI CD274 page).

Downloads stay under $GSE4573_CLDN4_CD274_DATA (default
/tmp/gse4573_cldn4_cd274) and are not committed.
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
DATA = Path(os.environ.get("GSE4573_CLDN4_CD274_DATA", "/tmp/gse4573_cldn4_cd274"))
SEED = 20260817
N_BOOT = 2000

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE4nnn/GSE4573/"
    "matrix/GSE4573_series_matrix.txt.gz"
)
GPL_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz"

# Locked claim genes. CD8A is covariate only (already NS vs CLDN4, PR 313).
CLAIM_GENES = ["CD274", "HLA-A", "HLA-B", "HLA-C"]
# Exact first-symbol / Entrez only. Do not substring-match "PDL1" (hits SPDL1).
ALIASES = {
    "CLDN4": ["CLDN4"],
    "CD274": ["CD274", "PDCD1LG1"],
    "HLA-A": ["HLA-A"],
    "HLA-B": ["HLA-B"],
    "HLA-C": ["HLA-C"],
    "CD8A": ["CD8A"],
    "TACSTD2": ["TACSTD2"],
}
ENTREZ = {
    "CLDN4": "1364",
    "CD274": "29126",
    "HLA-A": "3105",
    "HLA-B": "3106",
    "HLA-C": "3107",
    "CD8A": "925",
    "TACSTD2": "4070",
    "PDCD1LG2": "80380",  # PD-L2; on U133A but is not CD274
}
EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]
NAMED_PROBES = {
    "CLDN4": "201428_at",
    "CD8A": "205758_at",
    "TACSTD2": "202286_s_at",
}
# Common U133A HLA probes (official NetAffx / GPL96). Confirmed at runtime.
EXPECTED_HLA_PROBES = {
    "HLA-A": ["215313_x_at", "213932_x_at", "217456_x_at"],
    "HLA-B": ["209140_x_at", "208729_x_at", "211911_x_at", "204806_x_at"],
    "HLA-C": ["208812_x_at", "211799_x_at", "214459_x_at", "216526_x_at"],
}


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "gse4573-cldn4-cd274/1.0"})
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


def resolve_alias(index: pd.Index, aliases: list[str]) -> str | None:
    upper = {str(i).upper(): str(i) for i in index}
    for a in aliases:
        if a.upper() in upper:
            return upper[a.upper()]
    return None


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


def verdict(rho, p, n, present: bool) -> str:
    if not present:
        return "ABSENT"
    if n < 40:
        return "UNDERPOWERED"
    if not np.isfinite(rho) or not np.isfinite(p):
        return "NO_EVIDENCE"
    if p < 0.05 and rho < 0:
        return "NEGATIVE"
    if p < 0.05 and rho > 0:
        return "POSITIVE"
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
    sign = "+" if r > 0 else "−" if r < 0 else ""
    return f"{sign}{abs(r):.3f}"


def mwu_q4q1(x: pd.Series, y: pd.Series) -> dict:
    q = pd.qcut(x, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    a = y[q == "Q4"].dropna()
    b = y[q == "Q1"].dropna()
    if a.size < 4 or b.size < 4:
        return dict(
            n_Q4=int(a.size), n_Q1=int(b.size),
            median_Q4=np.nan, median_Q1=np.nan,
            U=np.nan, p=np.nan, rank_biserial_Q4_minus_Q1=np.nan,
        )
    U, p_mw = stats.mannwhitneyu(a.values, b.values, alternative="two-sided")
    rbc = 2 * U / (a.size * b.size) - 1
    return dict(
        n_Q4=int(a.size),
        n_Q1=int(b.size),
        median_Q4=float(np.median(a)),
        median_Q1=float(np.median(b)),
        U=float(U),
        p=float(p_mw),
        rank_biserial_Q4_minus_Q1=float(rbc),
    )


def search_annot_exact(annot: pd.DataFrame, aliases: list[str], entrez: str | None = None) -> pd.DataFrame:
    """Exact first-symbol or exact Entrez. No substring aliases (PDL1⊂SPDL1)."""
    cols = [c for c in ["ID", "Gene symbol", "Gene title", "Gene ID"] if c in annot.columns]
    first = annot["Gene symbol"].map(first_symbol)
    mask = first.isin(aliases)
    if entrez:
        mask = mask | (annot["Gene ID"].fillna("").astype(str) == str(entrez))
    return annot.loc[mask, cols].copy()


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
    n_array = int(probe_expr.shape[1])
    assert n_array == 130, n_array
    assert "CLDN4" in gene_expr.index

    # Probe confirmation: exact first-symbol + Entrez (no substring aliases).
    confirm_rows = []
    for gene, aliases in ALIASES.items():
        hits = search_annot_exact(annot, aliases, ENTREZ.get(gene))
        if hits.empty:
            confirm_rows.append({
                "gene": gene,
                "named_probe": NAMED_PROBES.get(gene, ""),
                "probe": "",
                "gpl96_symbol": "",
                "entrez": "",
                "title": "",
                "in_matrix": False,
                "is_named": False,
                "chosen_for_collapse": False,
                "note": "no GPL96 hit for locked aliases",
            })
            continue
        chosen = None
        if gene in gene_expr.index:
            chosen_hits = probe_audit.loc[probe_audit["gene"] == gene, "probe"]
            if len(chosen_hits):
                chosen = str(chosen_hits.iloc[0])
        for _, r in hits.iterrows():
            confirm_rows.append({
                "gene": gene,
                "named_probe": NAMED_PROBES.get(gene, ""),
                "probe": r["ID"],
                "gpl96_symbol": r.get("Gene symbol", ""),
                "entrez": r.get("Gene ID", ""),
                "title": r.get("Gene title", ""),
                "in_matrix": r["ID"] in probe_expr.index,
                "is_named": r["ID"] == NAMED_PROBES.get(gene, ""),
                "chosen_for_collapse": chosen is not None and r["ID"] == chosen,
                "note": "",
            })
    # PD-L2 is on U133A; it is not CD274 and is not substituted.
    pdl2 = search_annot_exact(annot, ["PDCD1LG2"], ENTREZ["PDCD1LG2"])
    if pdl2.empty:
        confirm_rows.append({
            "gene": "PDCD1LG2",
            "named_probe": "",
            "probe": "",
            "gpl96_symbol": "",
            "entrez": ENTREZ["PDCD1LG2"],
            "title": "",
            "in_matrix": False,
            "is_named": False,
            "chosen_for_collapse": False,
            "note": "PD-L2; related ligand, not CD274; not used as a substitute",
        })
    else:
        for _, r in pdl2.iterrows():
            confirm_rows.append({
                "gene": "PDCD1LG2",
                "named_probe": "",
                "probe": r["ID"],
                "gpl96_symbol": r.get("Gene symbol", ""),
                "entrez": r.get("Gene ID", ""),
                "title": r.get("Gene title", ""),
                "in_matrix": r["ID"] in probe_expr.index,
                "is_named": False,
                "chosen_for_collapse": False,
                "note": "PD-L2; related ligand, not CD274; not used as a substitute",
            })
    confirm = pd.DataFrame(confirm_rows)
    confirm.to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)

    resolved = {g: resolve_alias(gene_expr.index, ALIASES[g]) for g in ALIASES}
    cldn4 = gene_expr.loc[resolved["CLDN4"]]
    cd8a = gene_expr.loc[resolved["CD8A"]] if resolved["CD8A"] else None

    epi_present = [g for g in EPITHELIAL if g in gene_expr.index]
    epithelial = zscore_rows(gene_expr.loc[epi_present]).mean(axis=0) if len(epi_present) >= 3 else None

    mhci_genes = [resolved[g] for g in ["HLA-A", "HLA-B", "HLA-C"] if resolved[g]]
    mhci = zscore_rows(gene_expr.loc[mhci_genes]).mean(axis=0) if len(mhci_genes) >= 2 else None

    # Coverage
    cov_rows = []
    for g, aliases in ALIASES.items():
        hit = resolved[g]
        cov_rows.append({
            "set": f"gene:{g}",
            "n_listed": 1,
            "n_present": int(hit is not None),
            "resolved_symbol": hit or "",
            "aliases_tried": ",".join(aliases),
            "missing": "" if hit else g,
        })
    cov_rows.append({
        "set": "Epithelial",
        "n_listed": len(EPITHELIAL),
        "n_present": len(epi_present),
        "resolved_symbol": ",".join(epi_present),
        "aliases_tried": "",
        "missing": ",".join([g for g in EPITHELIAL if g not in gene_expr.index]),
    })
    cov_rows.append({
        "set": "MHC-I_HLAABC",
        "n_listed": 3,
        "n_present": len(mhci_genes),
        "resolved_symbol": ",".join(mhci_genes),
        "aliases_tried": "HLA-A,HLA-B,HLA-C",
        "missing": ",".join([g for g in ["HLA-A", "HLA-B", "HLA-C"] if not resolved[g]]),
    })
    pd.DataFrame(cov_rows).to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)

    # Sample annotation (nothing invented).
    keep_cols = [c for c in ["geo_accession", "title", "source_name_ch1", "characteristics_ch1", "description", "data_processing"] if c in meta.columns]
    sample = meta[keep_cols].copy()
    sample["cldn4_mas5"] = cldn4.reindex(sample.index).values
    if cd8a is not None:
        sample["cd8a_mas5"] = cd8a.reindex(sample.index).values
    if resolved["CD274"]:
        sample["cd274_mas5"] = gene_expr.loc[resolved["CD274"]].reindex(sample.index).values
    else:
        sample["cd274_mas5"] = np.nan
    for g in ["HLA-A", "HLA-B", "HLA-C"]:
        col = g.lower().replace("-", "") + "_mas5"
        sample[col] = gene_expr.loc[resolved[g]].reindex(sample.index).values if resolved[g] else np.nan
    if mhci is not None:
        sample["mhci_meanz"] = mhci.reindex(sample.index).values
    if epithelial is not None:
        sample["epithelial_meanz"] = epithelial.reindex(sample.index).values
    sample.to_csv(TABLES / "sample_annotation.tsv", sep="\t")

    n_cd274 = int(sample["cd274_mas5"].notna().sum())
    n_hla = int(sample[["hlaa_mas5", "hlab_mas5", "hlac_mas5"]].notna().all(axis=1).sum())
    labels = pd.DataFrame([
        {"field": "arrays_in_matrix", "public": "yes", "n": n_array, "note": "22283 probes x 130 GSM; 0 missing MAS5 values"},
        {"field": "unique_GSM", "public": "yes", "n": int(meta.index.nunique()), "note": "all unique"},
        {"field": "unique_title", "public": "yes", "n": int(meta["title"].nunique()), "note": "LS-* titles; all unique"},
        {"field": "GEO_overall_design_patients", "public": "yes (text only)", "n": 129, "note": "Series says '130 samples from 129 patients'; duplicate pair is not identified in sample metadata; not dropped"},
        {"field": "histology", "public": "series-level only", "n": 130, "note": "all described as lung squamous carcinoma; no LUAD rows"},
        {"field": "stage", "public": "no", "n": 0, "note": "not a GEO sample characteristic"},
        {"field": "OS / DSS time or event", "public": "no", "n": 0, "note": "prognosis paper; labels not deposited on GEO. No suppl folder."},
        {"field": "ICI / treatment", "public": "no", "n": 0, "note": "resected LUSC atlas, not an ICI series"},
        {"field": "tumor_percent / purity", "public": "no", "n": 0, "note": "epithelial RNA score is the only public purity proxy"},
        {"field": "CLDN4 finite", "public": "yes", "n": int(sample["cldn4_mas5"].notna().sum()), "note": "named probe 201428_at"},
        {"field": "CD274 finite", "public": "yes" if n_cd274 else "no", "n": n_cd274, "note": "Entrez 29126 / symbols CD274,PDCD1LG1: no GPL96 probe. Plus-2 probes 223834_at / 227458_at are not on U133A. PDCD1LG2 (PD-L2, 220049_s_at) is present and is not substituted."},
        {"field": "HLA-A/B/C finite (all three)", "public": "yes" if n_hla else "no", "n": n_hla, "note": "classical MHC-I; MHC-I mean-z uses these three only (not B2M/TAP)"},
        {"field": "pairwise CLDN4+CD274", "public": "yes" if n_cd274 else "no", "n": n_cd274, "note": "primary CD274 n; 0 means ABSENT on U133A"},
        {"field": "pairwise CLDN4+HLA-A/B/C", "public": "yes", "n": n_hla, "note": "primary MHC-I n"},
        {"field": "CLDN4 vs CD8A (PR 313, not this claim)", "public": "yes", "n": 130, "note": "already NS (ρ=−0.054, p=0.54); do not headline"},
    ])
    labels.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    # Spearman table: claim genes + MHC-I mean-z. CD8A recorded as covariate only.
    partners: list[tuple[str, pd.Series | None, bool]] = []
    for g in CLAIM_GENES:
        hit = resolved[g]
        partners.append((g, gene_expr.loc[hit] if hit else None, hit is not None))
    partners.append(("MHC-I_HLAABC", mhci, mhci is not None))

    gene_rows = []
    for name, y, present in partners:
        if not present or y is None:
            gene_rows.append({
                "gene_x": "CLDN4", "gene_y": name, "n": 0, "n_arrays": n_array,
                "present": False, "rho": np.nan, "p": np.nan, "ci_low": np.nan, "ci_high": np.nan,
                "rho_adj_epithelial": np.nan, "p_adj_epithelial": np.nan,
                "rho_adj_cd8a": np.nan, "p_adj_cd8a": np.nan,
                "verdict_crude": "ABSENT", "verdict_adj_epithelial": "ABSENT",
                "verdict_adj_cd8a": "ABSENT",
            })
            continue
        y = y.reindex(cldn4.index)
        crude = spearman_ci(cldn4.values, y.values)
        adj_e = partial_spearman(cldn4.values, y.values, epithelial.values) if epithelial is not None else dict(n=crude["n"], rho_adj=np.nan, p_adj=np.nan)
        adj_c = partial_spearman(cldn4.values, y.values, cd8a.values) if cd8a is not None else dict(n=crude["n"], rho_adj=np.nan, p_adj=np.nan)
        gene_rows.append({
            "gene_x": "CLDN4",
            "gene_y": name,
            "n": crude["n"],
            "n_arrays": n_array,
            "present": True,
            "rho": crude["rho"],
            "p": crude["p"],
            "ci_low": crude["ci_low"],
            "ci_high": crude["ci_high"],
            "rho_adj_epithelial": adj_e["rho_adj"],
            "p_adj_epithelial": adj_e["p_adj"],
            "rho_adj_cd8a": adj_c["rho_adj"],
            "p_adj_cd8a": adj_c["p_adj"],
            "verdict_crude": verdict(crude["rho"], crude["p"], crude["n"], True),
            "verdict_adj_epithelial": verdict(adj_e["rho_adj"], adj_e["p_adj"], adj_e["n"], True),
            "verdict_adj_cd8a": verdict(adj_c["rho_adj"], adj_c["p_adj"], adj_c["n"], True),
        })
    genes_df = pd.DataFrame(gene_rows)
    genes_df.to_csv(TABLES / "spearman_cldn4_vs_cd274_hla.tsv", sep="\t", index=False)

    # New CD274 cut: CLDN4 Q4 vs Q1 on CD274 / HLA / MHC-I.
    hl_rows = []
    for name, y, present in partners:
        rec = {"endpoint": name, "present": present}
        if not present or y is None:
            rec.update(mwu_q4q1(cldn4, pd.Series(np.nan, index=cldn4.index)))
            rec["verdict"] = "ABSENT"
        else:
            rec.update(mwu_q4q1(cldn4, y.reindex(cldn4.index)))
            rec["verdict"] = verdict(rec["rank_biserial_Q4_minus_Q1"], rec["p"], rec["n_Q4"] + rec["n_Q1"], True)
        hl_rows.append(rec)
    hl_df = pd.DataFrame(hl_rows)
    hl_df.to_csv(TABLES / "highlow_cldn4_cd274_hla.tsv", sep="\t", index=False)

    # One-row claim table (numeric).
    def pick(name: str, col: str):
        hit = genes_df[genes_df.gene_y == name]
        if hit.empty:
            return np.nan
        return hit.iloc[0][col]

    one = pd.DataFrame([{
        "dataset": "GSE4573",
        "histology": "LUSC",
        "platform": "GPL96 U133A MAS5",
        "n_arrays": n_array,
        "n_patients_geo_text": 129,
        "n_CLDN4_CD274": n_cd274,
        "n_CLDN4_HLAABC": n_hla,
        "CD274_present": bool(resolved["CD274"]),
        "CD274_rho": pick("CD274", "rho"),
        "CD274_p": pick("CD274", "p"),
        "CD274_rho_adj_epithelial": pick("CD274", "rho_adj_epithelial"),
        "CD274_p_adj_epithelial": pick("CD274", "p_adj_epithelial"),
        "CD274_rho_adj_cd8a": pick("CD274", "rho_adj_cd8a"),
        "CD274_p_adj_cd8a": pick("CD274", "p_adj_cd8a"),
        "HLA_A_rho": pick("HLA-A", "rho"),
        "HLA_A_p": pick("HLA-A", "p"),
        "HLA_B_rho": pick("HLA-B", "rho"),
        "HLA_B_p": pick("HLA-B", "p"),
        "HLA_C_rho": pick("HLA-C", "rho"),
        "HLA_C_p": pick("HLA-C", "p"),
        "MHCI_rho": pick("MHC-I_HLAABC", "rho"),
        "MHCI_p": pick("MHC-I_HLAABC", "p"),
        "MHCI_rho_adj_epithelial": pick("MHC-I_HLAABC", "rho_adj_epithelial"),
        "MHCI_p_adj_epithelial": pick("MHC-I_HLAABC", "p_adj_epithelial"),
        "MHCI_rho_adj_cd8a": pick("MHC-I_HLAABC", "rho_adj_cd8a"),
        "MHCI_p_adj_cd8a": pick("MHC-I_HLAABC", "p_adj_cd8a"),
        "CD274_Q4Q1_n": f"{int(hl_df.loc[hl_df.endpoint=='CD274','n_Q4'].iloc[0])} vs {int(hl_df.loc[hl_df.endpoint=='CD274','n_Q1'].iloc[0])}" if resolved["CD274"] else "ABSENT",
        "CD274_Q4Q1_p": hl_df.loc[hl_df.endpoint == "CD274", "p"].iloc[0],
        "note_cd8": "CLDN4 vs CD8A already NS in PR 313 (n=130, ρ=−0.054, p=0.54); not headlined",
    }])
    one.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    # Figures: CD274 (or HLA-A if CD274 absent) + MHC-I forest.
    def scatter(y, ylabel, title, dest_stem, color):
        fig, ax = plt.subplots(figsize=(5.2, 4.4))
        ax.scatter(np.log2(cldn4.values + 1), y, s=18, alpha=0.7, c=color, edgecolors="none")
        ax.set_xlabel("CLDN4  log2(MAS5+1)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        fig.tight_layout()
        fig.savefig(FIGURES / f"{dest_stem}.png", dpi=160)
        fig.savefig(FIGURES / f"{dest_stem}.pdf")
        plt.close(fig)

    if resolved["CD274"]:
        row = genes_df[genes_df.gene_y == "CD274"].iloc[0]
        scatter(
            np.log2(gene_expr.loc[resolved["CD274"]].reindex(cldn4.index).values + 1),
            "CD274  log2(MAS5+1)",
            f"GSE4573 LUSC  n={int(row.n)} arrays\nCLDN4 vs CD274  ρ={row.rho:.3f}  p={fmt_p(row.p)}",
            "fig1_cldn4_vs_cd274",
            "#1f4e79",
        )
    else:
        fig, ax = plt.subplots(figsize=(5.2, 3.2))
        ax.axis("off")
        ax.text(
            0.5, 0.5,
            "GSE4573 / GPL96 U133A\nCD274 (PD-L1) is ABSENT\nno official first-symbol or alias probe",
            ha="center", va="center", fontsize=12,
        )
        fig.tight_layout()
        fig.savefig(FIGURES / "fig1_cldn4_vs_cd274.png", dpi=160)
        fig.savefig(FIGURES / "fig1_cldn4_vs_cd274.pdf")
        plt.close(fig)

    plot = genes_df[genes_df.present].copy()
    if len(plot):
        plot = plot.sort_values("rho")
        fig, ax = plt.subplots(figsize=(6.2, 3.8))
        y_pos = np.arange(len(plot))
        ax.axvline(0, color="0.5", lw=0.8)
        ax.errorbar(
            plot["rho"], y_pos,
            xerr=[plot["rho"] - plot["ci_low"], plot["ci_high"] - plot["rho"]],
            fmt="o", color="#1f4e79", ecolor="#1f4e79", capsize=2, ms=6,
        )
        ax.set_yticks(y_pos)
        ax.set_yticklabels(plot["gene_y"])
        ax.set_xlabel("Spearman ρ vs CLDN4 (bootstrap 95% CI)")
        ax.set_title(f"GSE4573 LUSC  n={n_array} arrays  CLDN4 vs CD274 / HLA")
        fig.tight_layout()
        fig.savefig(FIGURES / "fig2_cldn4_cd274_hla_forest.png", dpi=160)
        fig.savefig(FIGURES / "fig2_cldn4_cd274_hla_forest.pdf")
        plt.close(fig)

    if mhci is not None:
        row = genes_df[genes_df.gene_y == "MHC-I_HLAABC"].iloc[0]
        scatter(
            mhci.reindex(cldn4.index).values,
            "MHC-I (HLA-A/B/C)  mean-z",
            f"GSE4573 LUSC  n={int(row.n)} arrays\nCLDN4 vs MHC-I  ρ={row.rho:.3f}  p={fmt_p(row.p)}",
            "fig3_cldn4_vs_mhci",
            "#6b3a2a",
        )

    summary = {
        "dataset": "GSE4573",
        "pmid": "16885343",
        "platform": "GPL96 Affymetrix HG-U133A",
        "histology": "LUSC (series-level; every sample described as lung squamous carcinoma)",
        "n_arrays": n_array,
        "n_probes": int(probe_expr.shape[0]),
        "n_genes_collapsed": int(gene_expr.shape[0]),
        "n_patients_geo_text": 129,
        "n_primary_cldn4_cd274": n_cd274,
        "n_primary_cldn4_hlaabc": n_hla,
        "processing": "MAS5 as deposited (linear signal). Spearman is rank-based.",
        "cldn4_probe": "201428_at",
        "cd274_present": bool(resolved["CD274"]),
        "cd274_resolved": resolved["CD274"],
        "hla_resolved": {g: resolved[g] for g in ["HLA-A", "HLA-B", "HLA-C"]},
        "mhci_definition": "mean of gene-wise z of HLA-A/B/C only (not B2M/TAP)",
        "epithelial_genes": epi_present,
        "cd8_not_headlined": {
            "reason": "PR 313 already NS",
            "n": 130,
            "rho": -0.05423268682178352,
            "p": 0.5399908033255645,
        },
        "pairs": genes_df.to_dict(orient="records"),
        "q4q1": hl_df.to_dict(orient="records"),
        "series_title": series.get("title", ""),
        "missing_public_labels": ["OS", "stage", "ICI", "tumor_percent", "patient-duplicate ID"],
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print("n_arrays", n_array, "CD274", resolved["CD274"], "HLA", {g: resolved[g] for g in ["HLA-A", "HLA-B", "HLA-C"]})
    print(genes_df.to_string(index=False))
    print(hl_df.to_string(index=False))
    print("[done]", TABLES)


if __name__ == "__main__":
    main()
