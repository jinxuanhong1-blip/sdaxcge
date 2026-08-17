#!/usr/bin/env python3
"""GSE101929 NSCLC Affy U133 Plus 2.0: CLDN4 vs CD8A / CD274 / ImmuneScore.

Additive CLDN4-only public series-matrix slice. No dual-high. No TACSTD2
companion. Unit is the array. Headline n is tumors, not the 66-array
series and not the ancestry / paper-histology splits.

Mitchell / Ryan NCI-Maryland AA vs EA NSCLC (PMID 29196495; GEO GSE101929,
GPL570). Series matrix is GCOS default as deposited. Paper Table 1
histology is aggregate only and is not mapped onto GSM IDs.

Downloads stay under $GSE101929_CLDN4_DATA (default /tmp/gse101929_cldn4)
and are not committed.
"""

from __future__ import annotations

import gzip
import io
import json
import os
import tarfile
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
DATA = Path(os.environ.get("GSE101929_CLDN4_DATA", "/tmp/gse101929_cldn4"))
SEED = 20260817
N_BOOT = 2000
HOLDS_N = 40

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE101nnn/GSE101929/"
    "matrix/GSE101929_series_matrix.txt.gz"
)
GPL_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"
ESTIMATE_TAR = "https://download.r-forge.r-project.org/src/contrib/estimate_1.0.13.tar.gz"

NAMED_PROBES = {
    "CLDN4": "201428_at",
    "CD8A": "205758_at",
    "CD274": "227458_at",
}
CD274_PROBES = ["223834_at", "227458_at"]
EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]
PARTNERS = ["CD8A", "CD274"]


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 10000:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "gse101929-cldn4/1.0"})
    with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as out:
        out.write(r.read())
    return dest


def parse_series_meta(path: Path) -> tuple[pd.DataFrame, dict]:
    series: dict[str, list[str]] = {}
    char_blocks: list[list[str]] = []
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
                if key == "characteristics_ch1":
                    char_blocks.append(vals)
                else:
                    sample_fields[key] = vals
    meta = pd.DataFrame(sample_fields)
    if "geo_accession" not in meta.columns:
        raise KeyError(f"no geo_accession in sample fields: {list(meta.columns)}")
    meta.index = meta["geo_accession"].astype(str)
    meta.index.name = "gsm"
    for block in char_blocks:
        keys = []
        values = []
        for v in block:
            if ": " in v:
                k, val = v.split(": ", 1)
            elif ":" in v:
                k, val = v.split(":", 1)
                val = val.strip()
            else:
                k, val = "unknown", v
            keys.append(k.strip())
            values.append(val.strip())
        key = keys[0] if len(set(keys)) == 1 else "mixed"
        col = key.replace(" ", "_").replace("-", "_")
        if col in meta.columns:
            col = col + "_2"
        meta[col] = values
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


def zscore_rows(expr: pd.DataFrame) -> pd.DataFrame:
    mu = expr.mean(axis=1)
    sd = expr.std(axis=1, ddof=1).replace(0, np.nan)
    return expr.sub(mu, axis=0).div(sd, axis=0)


def score_mean_z(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns), present
    z = zscore_rows(expr.loc[present].astype(float))
    return z.mean(axis=0), present


def load_estimate_immune_genes(tar_path: Path) -> list[str]:
    with tarfile.open(tar_path, "r:gz") as tf:
        member = next(m for m in tf.getmembers() if m.name.endswith("SI_geneset.gmt"))
        raw = tf.extractfile(member).read().decode("utf-8")
    immune = None
    for line in raw.splitlines():
        parts = line.strip().split("\t")
        if parts and parts[0] == "ImmuneSignature":
            immune = [g for g in parts[2:] if g]
    if not immune:
        raise RuntimeError("ImmuneSignature not found in SI_geneset.gmt")
    if len(immune) != 141:
        raise RuntimeError(f"expected 141 ImmuneSignature genes, got {len(immune)}")
    return immune


def pair_row(cldn4, y, epithelial, ancestry_bin, subset, partner, source):
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
        "rho_adj_ancestry": np.nan,
        "p_adj_ancestry": np.nan,
        "verdict_epithelial": verdict(adj_epi["rho_adj"], adj_epi["p_adj"], adj_epi["n"]),
        "verdict_ancestry": "NA",
    }
    if ancestry_bin is not None and np.unique(np.asarray(ancestry_bin, float)[np.isfinite(ancestry_bin)]).size > 1:
        adj_a = partial_spearman(cldn4.values, y.values, ancestry_bin)
        rec["rho_adj_ancestry"] = adj_a["rho_adj"]
        rec["p_adj_ancestry"] = adj_a["p_adj"]
        rec["verdict_ancestry"] = verdict(adj_a["rho_adj"], adj_a["p_adj"], adj_a["n"])
    return rec


def highlow(x: pd.Series, y: pd.Series, endpoint: str, subset: str) -> dict:
    if x.size < 16:
        return {
            "subset": subset,
            "endpoint": endpoint,
            "n_Q4": 0,
            "n_Q1": 0,
            "median_Q4": np.nan,
            "median_Q1": np.nan,
            "U": np.nan,
            "p": np.nan,
            "rank_biserial_Q4_minus_Q1": np.nan,
            "note": "n<16; quartiles not scored",
        }
    q = pd.qcut(x, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    a = y[q == "Q4"]
    b = y[q == "Q1"]
    if a.size < 5 or b.size < 5 or q.nunique() < 4:
        return {
            "subset": subset,
            "endpoint": endpoint,
            "n_Q4": int(a.size),
            "n_Q1": int(b.size),
            "median_Q4": np.nan,
            "median_Q1": np.nan,
            "U": np.nan,
            "p": np.nan,
            "rank_biserial_Q4_minus_Q1": np.nan,
            "note": "quartile split collapsed",
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
        "note": "",
    }


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    matrix_path = dl(DATA / "GSE101929_series_matrix.txt.gz", GEO_MATRIX)
    annot_path = dl(DATA / "GPL570.annot.gz", GPL_ANNOT)
    est_tar = dl(DATA / "estimate_1.0.13.tar.gz", ESTIMATE_TAR)
    immune_genes = load_estimate_immune_genes(est_tar)
    (DATA / "SI_geneset.ImmuneSignature.txt").write_text("\n".join(immune_genes) + "\n")

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

    if "tumor_normal_status" not in meta.columns or "race" not in meta.columns:
        raise KeyError(f"missing tumor/race columns: {list(meta.columns)}")
    meta["tissue"] = meta["tumor_normal_status"].astype(str).str.upper().str.strip()
    meta["ancestry"] = meta["race"].astype(str).str.upper().str.strip()
    if not set(meta["tissue"]).issubset({"T", "N"}):
        raise ValueError(f"unexpected tissue labels: {sorted(set(meta['tissue']))}")
    if not set(meta["ancestry"]).issubset({"AA", "EA"}):
        raise ValueError(f"unexpected ancestry labels: {sorted(set(meta['ancestry']))}")
    meta["patient"] = meta["individual"].astype(str).str.replace(r"^patient\s+", "", regex=True)
    meta["is_tumor"] = meta["tissue"] == "T"

    with gzip.open(annot_path, "rt", errors="replace") as fh:
        annot_text = fh.read()
    annot = None
    for sep_try in ("\t",):
        try:
            lines_a = annot_text.splitlines()
            a_start = next(i for i, l in enumerate(lines_a) if l.startswith("!platform_table_begin")) + 1
            a_end = next(i for i, l in enumerate(lines_a) if l.startswith("!platform_table_end"))
            annot = pd.read_csv(io.StringIO("\n".join(lines_a[a_start:a_end])), sep=sep_try, dtype=str)
            break
        except Exception:
            annot = None
    if annot is None:
        raise RuntimeError("failed to parse GPL570 annot")
    id_col = "ID" if "ID" in annot.columns else annot.columns[0]
    if "Gene symbol" in annot.columns:
        sym_col = "Gene symbol"
    elif "Gene Symbol" in annot.columns:
        sym_col = "Gene Symbol"
    else:
        raise KeyError(f"no Gene symbol column in GPL570 annot: {list(annot.columns)[:12]}")
    id2gene = pd.Series(
        annot[sym_col].map(first_symbol).values,
        index=annot[id_col].astype(str),
    )
    gene_expr, probe_audit = collapse_maxmean(probe_expr, id2gene)

    n_array = int(probe_expr.shape[1])
    n_nan = int(np.isnan(probe_expr.to_numpy()).sum())
    tumors = meta.index[meta["is_tumor"]].tolist()
    normals = meta.index[~meta["is_tumor"]].tolist()
    n_tumor = len(tumors)
    n_normal = len(normals)
    n_aa_t = int(((meta["is_tumor"]) & (meta["ancestry"] == "AA")).sum())
    n_ea_t = int(((meta["is_tumor"]) & (meta["ancestry"] == "EA")).sum())
    n_aa_n = int((~meta["is_tumor"] & (meta["ancestry"] == "AA")).sum())
    n_ea_n = int((~meta["is_tumor"] & (meta["ancestry"] == "EA")).sum())
    n_patients = int(meta["patient"].nunique())
    n_aa_pat = int(meta.loc[meta["ancestry"] == "AA", "patient"].nunique())
    n_ea_pat = int(meta.loc[meta["ancestry"] == "EA", "patient"].nunique())
    n_matched = int(
        meta.groupby("patient")["tissue"].apply(lambda s: set(s) == {"T", "N"}).sum()
    )

    for g, probe in NAMED_PROBES.items():
        if probe not in probe_expr.index:
            raise KeyError(f"named probe {probe} for {g} missing")
        if g not in gene_expr.index:
            raise KeyError(f"{g} missing after collapse")

    # Scores on tumors (headline) and on ancestry extras. z is within-subset.
    tumor_expr = gene_expr.loc[:, tumors]
    epithelial, epi_used = score_mean_z(tumor_expr, EPITHELIAL)
    immune, imm_used = score_mean_z(tumor_expr, immune_genes)

    # Coverage
    cov_rows = []
    for g in ["CLDN4", "CD8A", "CD274"] + EPITHELIAL:
        cov_rows.append({
            "gene": g,
            "class": "epithelial" if g in EPITHELIAL else "primary",
            "present": int(g in gene_expr.index),
            "collapse_probe": (
                str(probe_audit.loc[probe_audit["gene"] == g, "probe"].iloc[0])
                if g in gene_expr.index else ""
            ),
            "named_probe": NAMED_PROBES.get(g, ""),
        })
    imm_cov = []
    for g in immune_genes:
        present = g in gene_expr.index
        imm_cov.append({
            "gene": g,
            "present_after_collapse": int(present),
            "collapse_probe": (
                str(probe_audit.loc[probe_audit["gene"] == g, "probe"].iloc[0]) if present else ""
            ),
        })
    pd.DataFrame(cov_rows).to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(imm_cov).to_csv(TABLES / "immunescore_gene_coverage.tsv", sep="\t", index=False)

    probe_rows = []
    for g, named in NAMED_PROBES.items():
        hits = probe_audit[probe_audit["gene"] == g]
        collapse_probe = str(hits["probe"].iloc[0]) if len(hits) else ""
        probe_rows.append({
            "gene": g,
            "named_probe": named,
            "collapse_probe": collapse_probe,
            "named_equals_collapse": int(named == collapse_probe),
            "n_mapped_probes": int(len(hits)),
            "named_mean": float(probe_expr.loc[named, tumors].mean()) if named in probe_expr.index else np.nan,
            "collapse_mean": float(hits["mean"].iloc[0]) if len(hits) else np.nan,
        })
    for p274 in CD274_PROBES:
        probe_rows.append({
            "gene": "CD274",
            "named_probe": p274,
            "collapse_probe": str(probe_audit.loc[probe_audit["gene"] == "CD274", "probe"].iloc[0]),
            "named_equals_collapse": int(p274 == str(probe_audit.loc[probe_audit["gene"] == "CD274", "probe"].iloc[0])),
            "n_mapped_probes": 2,
            "named_mean": float(probe_expr.loc[p274, tumors].mean()),
            "collapse_mean": float(probe_audit.loc[probe_audit["gene"] == "CD274", "mean"].iloc[0]),
        })
    pd.DataFrame(probe_rows).to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)

    labels = pd.DataFrame([
        {"item": "arrays in series matrix", "public": "yes", "n": n_array,
         "note": "54,675 probes × 66 GSM; do not use as the pairwise n"},
        {"item": "unique GSM / unique titles", "public": "yes", "n": n_array,
         "note": "patient ID + T/N suffix; all unique"},
        {"item": "unique patients (GEO individual)", "public": "yes", "n": n_patients,
         "note": f"paper mRNA cohort 22 AA + 19 EA = 41; observed {n_aa_pat} AA + {n_ea_pat} EA"},
        {"item": "tumor arrays (tumor_normal status: T)", "public": "yes", "n": n_tumor,
         "note": "headline pairwise n; 16 AA + 16 EA"},
        {"item": "adjacent-normal arrays (N)", "public": "yes", "n": n_normal,
         "note": "inventoried; not mixed into the NSCLC pairwise tests"},
        {"item": "matched T+N patients", "public": "yes", "n": n_matched,
         "note": "paper text is 11 AA + 14 EA matched pairs; unmatched T or N remain"},
        {"item": "AA tumor arrays", "public": "yes", "n": n_aa_t,
         "note": "ancestry extra only; UNDERPOWERED for HOLDS (n<40)"},
        {"item": "EA tumor arrays", "public": "yes", "n": n_ea_t,
         "note": "ancestry extra only; UNDERPOWERED for HOLDS (n<40)"},
        {"item": "AA normal arrays", "public": "yes", "n": n_aa_n, "note": "inventory"},
        {"item": "EA normal arrays", "public": "yes", "n": n_ea_n, "note": "inventory"},
        {"item": "histology per array", "public": "no", "n": 0,
         "note": "not a GEO characteristic; paper Table 1 is aggregate only (ADC 27 / SCC 1 / other 4 of 32)"},
        {"item": "Stage", "public": "yes", "n": n_tumor,
         "note": "I/II/III on every array; not a covariate for this claim"},
        {"item": "OS days / months / years + lung-cancer death", "public": "yes", "n": n_array,
         "note": "deposited; not this claim"},
        {"item": "ICI / treatment", "public": "no", "n": 0,
         "note": "resected NCI-Maryland atlas, not an ICI series"},
        {"item": "numeric tumor % / ABSOLUTE purity", "public": "no", "n": 0,
         "note": "macro-dissected; only public proxy is RNA epithelial mean-z"},
        {"item": "CLDN4 finite on tumors (max-mean)", "public": "yes", "n": n_tumor,
         "note": f"collapse {probe_audit.loc[probe_audit['gene']=='CLDN4','probe'].iloc[0]}"},
        {"item": "CD8A finite on tumors", "public": "yes", "n": n_tumor,
         "note": "named 205758_at"},
        {"item": "CD274 finite on tumors", "public": "yes", "n": n_tumor,
         "note": f"collapse {probe_audit.loc[probe_audit['gene']=='CD274','probe'].iloc[0]}"},
        {"item": "ESTIMATE ImmuneSignature genes present", "public": "yes",
         "n": int(sum(g in gene_expr.index for g in immune_genes)),
         "note": f"{int(sum(g in gene_expr.index for g in immune_genes))}/141 after collapse"},
        {"item": "primary pairwise n (tumors, CLDN4 + partners)", "public": "yes", "n": n_tumor,
         "note": "this is the n used below; UNDERPOWERED for HOLDS (n<40)"},
        {"item": "missing values in series matrix", "public": "yes", "n": n_nan, "note": "probe × sample"},
    ])
    labels.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    # Per-array annotation (scores on tumors only; normals get NA scores)
    ann = meta.copy()
    for g in ["CLDN4", "CD8A", "CD274"]:
        ann[g] = gene_expr.loc[g]
    ann["epithelial_z"] = np.nan
    ann.loc[tumors, "epithelial_z"] = epithelial.reindex(tumors).values
    ann["ImmuneScore"] = np.nan
    ann.loc[tumors, "ImmuneScore"] = immune.reindex(tumors).values
    keep_cols = [
        "title", "patient", "tissue", "ancestry", "gender", "age", "Stage",
        "smoking_pack_years", "survival_(days)", "survival_(months)",
        "death_due_to_lung_cancer_(all_years)",
        "CLDN4", "CD8A", "CD274", "epithelial_z", "ImmuneScore",
    ]
    keep_cols = [c for c in keep_cols if c in ann.columns]
    ann[keep_cols].to_csv(TABLES / "sample_annotation.tsv", sep="\t")

    subsets = {
        "tumors": tumors,
        "AA_tumors": meta.index[(meta["is_tumor"]) & (meta["ancestry"] == "AA")].tolist(),
        "EA_tumors": meta.index[(meta["is_tumor"]) & (meta["ancestry"] == "EA")].tolist(),
    }
    ancestry_all = (meta.loc[tumors, "ancestry"] == "AA").astype(float)

    rows = []
    hl_rows = []
    for subset, idx in subsets.items():
        sub = gene_expr.loc[:, idx]
        cldn4 = sub.loc["CLDN4"]
        epi, _ = score_mean_z(sub, EPITHELIAL)
        imm, _ = score_mean_z(sub, immune_genes)
        abin = ancestry_all.reindex(idx).values if subset == "tumors" else None
        for partner in PARTNERS:
            rows.append(pair_row(cldn4, sub.loc[partner], epi, abin, subset, partner, "collapsed_maxmean"))
        rows.append(pair_row(cldn4, imm, epi, abin, subset, "ImmuneScore", "estimate_immunesignature_meanz"))
        # named-probe sensitivity on tumors / extras
        rows.append(pair_row(
            probe_expr.loc["201428_at", idx],
            probe_expr.loc["205758_at", idx],
            epi, abin, subset, "CD8A", "named_201428_at_vs_205758_at",
        ))
        for p274 in CD274_PROBES:
            rows.append(pair_row(
                probe_expr.loc["201428_at", idx],
                probe_expr.loc[p274, idx],
                epi, abin, subset, "CD274", f"named_201428_at_vs_{p274}",
            ))
        crude_epi = spearman_ci(cldn4.values, epi.values)
        rows.append({
            "subset": subset, "partner": "epithelial_z", "source": "context",
            "n": crude_epi["n"], "rho": crude_epi["rho"], "p": crude_epi["p"],
            "ci_low": crude_epi["ci_low"], "ci_high": crude_epi["ci_high"],
            "rho_adj_epithelial": np.nan, "p_adj_epithelial": np.nan,
            "rho_adj_ancestry": np.nan, "p_adj_ancestry": np.nan,
            "verdict_epithelial": "COMPANION", "verdict_ancestry": "NA",
        })
        for endpoint, series_y in [("CD8A", sub.loc["CD8A"]), ("CD274", sub.loc["CD274"]), ("ImmuneScore", imm)]:
            hl_rows.append(highlow(cldn4, series_y, endpoint, subset))

    # ancestry contrast on tumors (context)
    crude_a = spearman_ci(gene_expr.loc["CLDN4", tumors].values, ancestry_all.values)
    rows.append({
        "subset": "tumors", "partner": "ancestry_AA", "source": "context",
        "n": crude_a["n"], "rho": crude_a["rho"], "p": crude_a["p"],
        "ci_low": crude_a["ci_low"], "ci_high": crude_a["ci_high"],
        "rho_adj_epithelial": np.nan, "p_adj_epithelial": np.nan,
        "rho_adj_ancestry": np.nan, "p_adj_ancestry": np.nan,
        "verdict_epithelial": "COMPANION", "verdict_ancestry": "NA",
    })

    corr = pd.DataFrame(rows)
    corr.to_csv(TABLES / "spearman.tsv", sep="\t", index=False)
    pd.DataFrame(hl_rows).to_csv(TABLES / "highlow_cldn4.tsv", sep="\t", index=False)

    # Ancestry contrast MWU on tumors
    anc_rows = []
    for gene in ["CLDN4", "CD8A", "CD274"]:
        a = gene_expr.loc[gene, subsets["AA_tumors"]]
        b = gene_expr.loc[gene, subsets["EA_tumors"]]
        U, p_mw = stats.mannwhitneyu(a.values, b.values, alternative="two-sided")
        rbc = 2 * U / (a.size * b.size) - 1
        anc_rows.append({
            "gene": gene,
            "n_AA": int(a.size),
            "n_EA": int(b.size),
            "median_AA": float(np.median(a)),
            "median_EA": float(np.median(b)),
            "U": float(U),
            "p": float(p_mw),
            "rank_biserial_AA_minus_EA": float(rbc),
        })
    a = immune.reindex(subsets["AA_tumors"])
    b = immune.reindex(subsets["EA_tumors"])
    U, p_mw = stats.mannwhitneyu(a.values, b.values, alternative="two-sided")
    anc_rows.append({
        "gene": "ImmuneScore",
        "n_AA": int(a.size),
        "n_EA": int(b.size),
        "median_AA": float(np.median(a)),
        "median_EA": float(np.median(b)),
        "U": float(U),
        "p": float(p_mw),
        "rank_biserial_AA_minus_EA": float(2 * U / (a.size * b.size) - 1),
    })
    pd.DataFrame(anc_rows).to_csv(TABLES / "ancestry_contrast.tsv", sep="\t", index=False)

    def pick(subset, partner, source=None):
        hit = corr[(corr.subset == subset) & (corr.partner == partner)]
        if source is not None:
            hit = hit[hit.source == source]
        return hit.iloc[0].to_dict() if len(hit) else {}

    t_cd8 = pick("tumors", "CD8A", "collapsed_maxmean")
    t_pdl1 = pick("tumors", "CD274", "collapsed_maxmean")
    t_imm = pick("tumors", "ImmuneScore", "estimate_immunesignature_meanz")
    t_epi = pick("tumors", "epithelial_z", "context")

    one = pd.DataFrame([{
        "dataset": "GSE101929 Mitchell NCI-Maryland",
        "histology": "NSCLC mixed; per-array histology not deposited",
        "platform": "GPL570 U133 Plus 2.0 GCOS default as deposited",
        "n_arrays_series": n_array,
        "n_tumors": n_tumor,
        "n_AA_tumors": n_aa_t,
        "n_EA_tumors": n_ea_t,
        "CLDN4": str(probe_audit.loc[probe_audit["gene"] == "CLDN4", "probe"].iloc[0]),
        "CD8A": "205758_at",
        "CD274": str(probe_audit.loc[probe_audit["gene"] == "CD274", "probe"].iloc[0]),
        "ImmuneScore": f"ESTIMATE ImmuneSignature mean-z, {len(imm_used)}/{len(immune_genes)}",
        "CLDN4_CD8A_rho": t_cd8.get("rho"),
        "CLDN4_CD8A_p": t_cd8.get("p"),
        "CLDN4_CD8A_adj_epi_rho": t_cd8.get("rho_adj_epithelial"),
        "CLDN4_CD8A_adj_epi_p": t_cd8.get("p_adj_epithelial"),
        "CLDN4_CD274_rho": t_pdl1.get("rho"),
        "CLDN4_CD274_p": t_pdl1.get("p"),
        "CLDN4_CD274_adj_epi_rho": t_pdl1.get("rho_adj_epithelial"),
        "CLDN4_CD274_adj_epi_p": t_pdl1.get("p_adj_epithelial"),
        "CLDN4_ImmuneScore_rho": t_imm.get("rho"),
        "CLDN4_ImmuneScore_p": t_imm.get("p"),
        "CLDN4_ImmuneScore_adj_epi_rho": t_imm.get("rho_adj_epithelial"),
        "CLDN4_ImmuneScore_adj_epi_p": t_imm.get("p_adj_epithelial"),
        "verdict": t_cd8.get("verdict_epithelial"),
        "OS_ICI": "OS deposited; ICI not deposited",
    }])
    one.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    # Figures
    colors = {"AA": "#5c3d8a", "EA": "#2c5f8a"}
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))
    ys = [
        (gene_expr.loc["CD8A", tumors], "CD8A  GCOS", t_cd8),
        (gene_expr.loc["CD274", tumors], "CD274  GCOS (max-mean)", t_pdl1),
        (immune.reindex(tumors), "ImmuneScore  mean-z", t_imm),
    ]
    x = gene_expr.loc["CLDN4", tumors]
    anc = meta.loc[tumors, "ancestry"]
    for ax, (y, ylab, row) in zip(axes, ys):
        for a, c in colors.items():
            m = anc == a
            ax.scatter(x[m], y[m], s=26, alpha=0.85, c=c, edgecolors="none",
                       label=f"{a} n={int(m.sum())}")
        ax.set_xlabel("CLDN4  GCOS")
        ax.set_ylabel(ylab)
        ax.set_title(f"tumors n={int(row['n'])}\nρ={row['rho']:+.3f}  p={fmt_p(row['p'])}")
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_cd274_immunescore.png", dpi=160)
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_cd274_immunescore.pdf")
    plt.close(fig)

    plot = corr[
        corr.partner.isin(["CD8A", "CD274", "ImmuneScore"])
        & corr.source.isin(["collapsed_maxmean", "estimate_immunesignature_meanz"])
    ].copy()
    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    y_pos = np.arange(len(plot))
    ax.axvline(0, color="0.5", lw=0.8)
    ax.errorbar(
        plot["rho"], y_pos,
        xerr=[plot["rho"] - plot["ci_low"], plot["ci_high"] - plot["rho"]],
        fmt="o", color="#2c5f8a", ecolor="#2c5f8a", capsize=2, ms=5,
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{r.subset}  vs {r.partner}  n={int(r.n)}" for r in plot.itertuples()])
    ax.set_xlabel("Spearman ρ (bootstrap 95% CI)")
    ax.set_title("GSE101929  CLDN4 vs CD8A / CD274 / ImmuneScore")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_forest.png", dpi=160)
    fig.savefig(FIGURES / "fig2_forest.pdf")
    plt.close(fig)

    summary = {
        "dataset": "GSE101929",
        "pmid": "29196495",
        "platform": "GPL570 Affymetrix HG-U133 Plus 2.0",
        "processing": "GCOS default as deposited in the series matrix. Paper used Partek RMA; that matrix is not the GEO table. Spearman is rank-based.",
        "n_arrays": n_array,
        "n_tumors": n_tumor,
        "n_normals": n_normal,
        "n_patients": n_patients,
        "n_AA_tumors": n_aa_t,
        "n_EA_tumors": n_ea_t,
        "n_probes": int(probe_expr.shape[0]),
        "n_genes_collapsed": int(gene_expr.shape[0]),
        "n_missing_values": n_nan,
        "epithelial_genes_used": epi_used,
        "immune_genes_used": len(imm_used),
        "immune_genes_listed": len(immune_genes),
        "cldn4_collapse_probe": str(probe_audit.loc[probe_audit["gene"] == "CLDN4", "probe"].iloc[0]),
        "cd8a_collapse_probe": str(probe_audit.loc[probe_audit["gene"] == "CD8A", "probe"].iloc[0]),
        "cd274_collapse_probe": str(probe_audit.loc[probe_audit["gene"] == "CD274", "probe"].iloc[0]),
        "holds_rule": "n>=40, residual ρ<0, residual p<0.05 (same as PR 229 / GSE4573)",
        "no_dual_high": True,
        "histology_per_array": False,
        "paper_table1_histology_aggregate": {"ADC": 27, "SCC": 1, "other": 4, "n": 32},
        "primary": {
            "tumors_CD8A": t_cd8,
            "tumors_CD274": t_pdl1,
            "tumors_ImmuneScore": t_imm,
            "tumors_epithelial_z": t_epi,
        },
        "missing_public_labels": ["per_array_histology", "ICI", "numeric_tumor_percent"],
        "do_not_use_n": {
            "series_arrays_including_normal": 66,
            "paper_mRNA_patients": 41,
            "reason": "pairwise immune tests use tumor arrays only (n=32)",
        },
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")

    show = corr[corr.source.isin(["collapsed_maxmean", "estimate_immunesignature_meanz", "context"])]
    print("n_arrays", n_array, "tumors", n_tumor, "AA_T", n_aa_t, "EA_T", n_ea_t)
    print(show[["subset", "partner", "source", "n", "rho", "p",
                "rho_adj_epithelial", "p_adj_epithelial",
                "rho_adj_ancestry", "p_adj_ancestry", "verdict_epithelial"]].to_string(index=False))
    print("[done]", TABLES)


if __name__ == "__main__":
    main()
