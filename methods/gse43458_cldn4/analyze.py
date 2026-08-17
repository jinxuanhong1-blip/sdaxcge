#!/usr/bin/env python3
"""GSE43458 never-smoker LUAD Affy HuGene 1.0 ST: CLDN4 vs CD8A / CD274 / ImmuneScore.

Additive CLDN4-only slice. Unit is the array. Primary set is GEO
never-smoker lung adenocarcinoma *tumor* arrays (n=40). The 30 paired
normals and the 40 smoker tumors are not folded into that n.

No dual-high (no TACSTD2×CLDN4 high/high class). TACSTD2 is a same-run
companion only.

ImmuneScore = mean of gene-wise z-scores of the Yoshihara 2013 ESTIMATE
ImmuneSignature (141 genes from SI_geneset.gmt in estimate 1.0.13).
This is not the official R ssGSEA ImmuneScore.

Downloads stay under $GSE43458_CLDN4_DATA (default /tmp/gse43458_cldn4)
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
DATA = Path(os.environ.get("GSE43458_CLDN4_DATA", "/tmp/gse43458_cldn4"))
SEED = 20260817
N_BOOT = 2000

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE43nnn/GSE43458/"
    "matrix/GSE43458_series_matrix.txt.gz"
)
GPL_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz"
ESTIMATE_URLS = [
    "https://download.r-forge.r-project.org/src/contrib/estimate_1.0.13.tar.gz",
    "https://bioconductor.org/packages/release/bioc/src/contrib/estimate_1.0.13.tar.gz",
]

EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]
PARTNERS = ["CD8A", "CD274"]

# PR 229 HOLDS rule (given): n>=40, partial Spearman rho<0, p<0.05.
HOLDS_N = 40

UA = "gse43458-cldn4/1.0"


def dl(dest: Path, url: str, min_size: int = 10000) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > min_size:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
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


def verdict(rho, p, n) -> str:
    if n < HOLDS_N:
        return "UNDERPOWERED"
    if not np.isfinite(rho) or not np.isfinite(p):
        return "NO_EVIDENCE"
    if rho < 0 and p < 0.05:
        return "HOLDS"
    if rho > 0 and p < 0.05:
        return "OPPOSITE"
    return "NO_EVIDENCE"


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


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


def fetch_immune_genes() -> list[str]:
    dest = DATA / "estimate_1.0.13.tar.gz"
    last_err = None
    for url in ESTIMATE_URLS:
        try:
            dl(dest, url)
            return load_estimate_immune_genes(dest)
        except Exception as exc:  # noqa: BLE001 — try next mirror
            last_err = exc
            if dest.exists():
                dest.unlink()
    # Last-resort: public GMT text (same Yoshihara 2013 ImmuneSignature).
    gmt = DATA / "SI_geneset.gmt"
    try:
        dl(
            gmt,
            "https://raw.githubusercontent.com/yhoogstrate/estimate/master/inst/extdata/SI_geneset.gmt",
            min_size=500,
        )
        immune = None
        for line in gmt.read_text().splitlines():
            parts = line.strip().split("\t")
            if parts and parts[0] == "ImmuneSignature":
                immune = [g for g in parts[2:] if g]
        if immune and len(immune) == 141:
            return immune
    except Exception as exc:  # noqa: BLE001
        last_err = exc
    raise RuntimeError(f"could not load ESTIMATE ImmuneSignature: {last_err}")


def pair_row(name: str, x: pd.Series, y: pd.Series, epithelial: pd.Series | None, cohort: str, extra: dict | None = None):
    crude = spearman_ci(x.values, y.values)
    if epithelial is not None:
        adj = partial_spearman(x.values, y.values, epithelial.reindex(x.index).values)
    else:
        adj = dict(n=crude["n"], rho_adj=np.nan, p_adj=np.nan)
    rec = {
        "cohort": cohort,
        "pair": name,
        "n": crude["n"],
        "rho": crude["rho"],
        "p": crude["p"],
        "ci_low": crude["ci_low"],
        "ci_high": crude["ci_high"],
        "rho_adj_epithelial": adj["rho_adj"],
        "p_adj_epithelial": adj["p_adj"],
        "verdict_crude": verdict(crude["rho"], crude["p"], crude["n"]),
        "verdict_adj": verdict(adj["rho_adj"], adj["p_adj"], adj["n"]),
    }
    if extra:
        rec.update(extra)
    return rec


def highlow(x: pd.Series, y: pd.Series, endpoint: str, cohort: str) -> dict:
    q = pd.qcut(x, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    a = y[q == "Q4"]
    b = y[q == "Q1"]
    if a.size < 5 or b.size < 5:
        return {
            "cohort": cohort,
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
        "cohort": cohort,
        "endpoint": endpoint,
        "n_Q4": int(a.size),
        "n_Q1": int(b.size),
        "median_Q4": float(np.median(a)),
        "median_Q1": float(np.median(b)),
        "U": float(U),
        "p": float(p_mw),
        "rank_biserial_Q4_minus_Q1": float(rbc),
    }


def score_mean_z(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in expr.index]
    if len(present) < 3:
        return pd.Series(dtype=float), present
    z = zscore_rows(expr.loc[present])
    return z.mean(axis=0), present


def analyze_cohort(name: str, samples: list[str], gene_expr: pd.DataFrame, immune_genes: list[str]) -> tuple[list[dict], list[dict], dict]:
    sub = gene_expr.loc[:, samples]
    cldn4 = sub.loc["CLDN4"]
    epi, epi_used = score_mean_z(sub, EPITHELIAL)
    immune, imm_used = score_mean_z(sub, immune_genes)
    rows = []
    for g in PARTNERS:
        extra = {"n_genes_in_score": 1, "genes_used": g, "present": g in sub.index}
        if g not in sub.index:
            rows.append({
                "cohort": name, "pair": f"CLDN4 vs {g}", "n": 0,
                "rho": np.nan, "p": np.nan, "ci_low": np.nan, "ci_high": np.nan,
                "rho_adj_epithelial": np.nan, "p_adj_epithelial": np.nan,
                "verdict_crude": "ABSENT", "verdict_adj": "ABSENT", **extra,
            })
            continue
        rows.append(pair_row(f"CLDN4 vs {g}", cldn4, sub.loc[g], epi if len(epi) else None, name, extra))
    rows.append(pair_row(
        "CLDN4 vs ImmuneScore",
        cldn4,
        immune,
        epi if len(epi) else None,
        name,
        {"n_genes_in_score": len(imm_used), "genes_used": ",".join(imm_used), "present": True},
    ))
    if "TACSTD2" in sub.index:
        rows.append(pair_row(
            "CLDN4 vs TACSTD2",
            cldn4,
            sub.loc["TACSTD2"],
            epi if len(epi) else None,
            name,
            {"n_genes_in_score": 1, "genes_used": "TACSTD2", "present": True, "note": "companion; not dual-high"},
        ))
    if len(epi):
        rows.append(pair_row(
            "CLDN4 vs epithelial mean-z",
            cldn4,
            epi,
            None,
            name,
            {"n_genes_in_score": len(epi_used), "genes_used": ",".join(epi_used), "present": True, "note": "companion"},
        ))
    hl = [
        highlow(cldn4, sub.loc["CD8A"], "CD8A", name) if "CD8A" in sub.index else {},
        highlow(cldn4, sub.loc["CD274"], "CD274", name) if "CD274" in sub.index else {},
        highlow(cldn4, immune, "ImmuneScore", name) if len(immune) else {},
    ]
    hl = [h for h in hl if h]
    info = {
        "cohort": name,
        "n": int(sub.shape[1]),
        "immune_genes_used": len(imm_used),
        "immune_genes_listed": len(immune_genes),
        "epithelial_genes_used": len(epi_used),
        "epithelial_genes": epi_used,
        "cldn4_finite": int(cldn4.notna().sum()),
        "cd8a_finite": int(sub.loc["CD8A"].notna().sum()) if "CD8A" in sub.index else 0,
        "cd274_finite": int(sub.loc["CD274"].notna().sum()) if "CD274" in sub.index else 0,
    }
    return rows, hl, info


def classify_samples(meta: pd.DataFrame) -> pd.DataFrame:
    """GEO characteristics only. Do not invent smoking or tissue from titles."""
    out = meta.copy()
    title = out["title"].astype(str) if "title" in out.columns else pd.Series("", index=out.index)

    tissue = None
    for col in ("tissue", "source_name_ch1"):
        if col in out.columns:
            tissue = out[col].astype(str)
            break
    if tissue is None:
        tissue = title

    smoking = None
    for col in ("smoking_status", "smoking", "smoker"):
        if col in out.columns:
            smoking = out[col].astype(str)
            break
    if smoking is None:
        smoking = pd.Series("", index=out.index)

    histo = None
    for col in ("histology", "disease_state", "disease"):
        if col in out.columns:
            histo = out[col].astype(str)
            break
    if histo is None:
        histo = pd.Series("", index=out.index)

    tlow = tissue.str.lower()
    slow = smoking.str.lower()
    hlow = histo.str.lower()
    title_l = title.str.lower()

    is_normal = tlow.str.contains("normal") | title_l.str.contains("normal lung")
    is_tumor = (
        tlow.str.contains("tumor")
        | tlow.str.contains("adenocarcinoma")
        | title_l.str.contains("adenocarcinoma tissue")
    ) & ~is_normal
    is_never = slow.str.contains("never") | title_l.str.contains("never-smoker")
    is_smoker = (slow.str.contains("smoker") & ~is_never) | (
        title_l.str.contains("from smoker") & ~title_l.str.contains("never")
    )
    is_luad = hlow.str.contains("adenocarcinoma") | title_l.str.contains("adenocarcinoma")

    out["tissue_class"] = np.where(is_normal, "normal", np.where(is_tumor, "tumor", "other"))
    out["smoking_class"] = np.where(is_never, "never-smoker", np.where(is_smoker, "smoker", "unknown"))
    out["is_luad"] = is_luad
    out["is_never_tumor"] = is_tumor & is_never & is_luad
    out["is_smoker_tumor"] = is_tumor & is_smoker & is_luad
    out["is_any_tumor"] = is_tumor & is_luad
    return out


def jsonable(v):
    if isinstance(v, float) and not np.isfinite(v):
        return None
    if isinstance(v, (np.floating,)):
        if not np.isfinite(v):
            return None
        return float(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    return v


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    matrix_path = dl(DATA / "GSE43458_series_matrix.txt.gz", GEO_MATRIX)
    annot_path = dl(DATA / "GPL6244.annot.gz", GPL_ANNOT)
    immune_genes = fetch_immune_genes()
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
    meta = meta.loc[probe_expr.columns]
    meta = classify_samples(meta)

    with gzip.open(annot_path, "rt", errors="replace") as fh:
        annot = parse_soft_table(fh.read())
    if "Gene symbol" not in annot.columns:
        raise KeyError(f"GPL6244 annot missing Gene symbol; have {list(annot.columns)}")
    annot["symbol"] = annot["Gene symbol"].map(first_symbol)
    id2gene = annot.drop_duplicates("ID").set_index("ID")["symbol"]

    gene_expr, probe_audit = collapse_maxmean(probe_expr, id2gene)

    confirm_rows = []
    for gene in ["CLDN4", "CD8A", "CD274", "TACSTD2"] + EPITHELIAL:
        hits = probe_audit[probe_audit.gene == gene]
        kept = str(hits.iloc[0]["probe"]) if len(hits) else ""
        confirm_rows.append({
            "gene": gene,
            "n_probes_mapped": int(len(hits)),
            "kept_maxmean_probe": kept,
            "in_collapsed_matrix": gene in gene_expr.index,
        })
    extra_probes = probe_audit[probe_audit.gene.isin(["CLDN4", "CD8A", "CD274", "TACSTD2"] + EPITHELIAL)].copy()
    extra_probes.to_csv(TABLES / "probe_all_mapped.tsv", sep="\t", index=False)
    pd.DataFrame(confirm_rows).to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)

    never = meta.index[meta["is_never_tumor"]].tolist()
    smoker = meta.index[meta["is_smoker_tumor"]].tolist()
    tumors = meta.index[meta["is_any_tumor"]].tolist()
    normals = meta.index[meta["tissue_class"] == "normal"].tolist()

    sample = meta.copy()
    for gene, col in [("CLDN4", "cldn4_rma"), ("CD8A", "cd8a_rma"), ("CD274", "cd274_rma"), ("TACSTD2", "tacstd2_rma")]:
        sample[col] = gene_expr.loc[gene].reindex(sample.index).values if gene in gene_expr.index else np.nan
    sample.to_csv(TABLES / "sample_annotation.tsv", sep="\t")

    n_array = int(probe_expr.shape[1])
    if n_array != 110:
        print(f"[warn] expected 110 arrays, got {n_array}", flush=True)
    assert "CLDN4" in gene_expr.index, "CLDN4 absent after collapse"
    assert "CD8A" in gene_expr.index, "CD8A absent after collapse"
    assert "CD274" in gene_expr.index, "CD274 absent after collapse"
    assert len(never) == 40, f"expected 40 never-smoker tumors, got {len(never)}"
    assert len(smoker) == 40, f"expected 40 smoker tumors, got {len(smoker)}"
    assert len(normals) == 30, f"expected 30 normals, got {len(normals)}"

    counts = (
        meta.groupby(["tissue_class", "smoking_class", "is_luad"], dropna=False)
        .size()
        .reset_index(name="n")
    )
    counts.to_csv(TABLES / "histology_counts.tsv", sep="\t", index=False)

    pair_rows, hl_rows, infos = [], [], []
    for name, gsms in [
        ("never_smoker_LUAD_tumor", never),
        ("smoker_LUAD_tumor_extra", smoker),
        ("all_LUAD_tumor_mixed", tumors),
    ]:
        pr, hl, info = analyze_cohort(name, gsms, gene_expr, immune_genes)
        pair_rows.extend(pr)
        hl_rows.extend(hl)
        infos.append(info)

    pairs_df = pd.DataFrame(pair_rows)
    pairs_df.to_csv(TABLES / "spearman_cldn4_vs_partners.tsv", sep="\t", index=False)
    hl_df = pd.DataFrame(hl_rows)
    hl_df.to_csv(TABLES / "highlow_cldn4.tsv", sep="\t", index=False)

    cov = pd.DataFrame({
        "gene": immune_genes,
        "present_after_collapse": [g in gene_expr.index for g in immune_genes],
    })
    cov.to_csv(TABLES / "immunescore_gene_coverage.tsv", sep="\t", index=False)

    cldn4_probe = str(probe_audit.loc[probe_audit.gene == "CLDN4", "probe"].iloc[0])
    cd8a_probe = str(probe_audit.loc[probe_audit.gene == "CD8A", "probe"].iloc[0])
    cd274_probe = str(probe_audit.loc[probe_audit.gene == "CD274", "probe"].iloc[0])
    tac_probe = (
        str(probe_audit.loc[probe_audit.gene == "TACSTD2", "probe"].iloc[0])
        if (probe_audit.gene == "TACSTD2").any()
        else ""
    )

    labels = pd.DataFrame([
        {"field": "arrays in series matrix", "public": "yes", "n": n_array, "note": "HuGene 1.0 ST RMA log2 as deposited (BRB-Array Tools)"},
        {"field": "unique GSM", "public": "yes", "n": int(meta.index.nunique()), "note": "all unique"},
        {"field": "never-smoker LUAD tumor (PRIMARY)", "public": "yes", "n": len(never), "note": "GEO smoking status Never-smoker + tumor; this is the n used below"},
        {"field": "smoker LUAD tumor (extra, not the claim)", "public": "yes", "n": len(smoker), "note": "NOT folded into never-smoker n"},
        {"field": "all LUAD tumors mixed", "public": "yes", "n": len(tumors), "note": "40 never + 40 smoker; composition, not the claim"},
        {"field": "paired normal lung (excluded)", "public": "yes", "n": len(normals), "note": "paired to 30 never-smoker cases; tumor-only slice"},
        {"field": "ICI / treatment", "public": "no", "n": 0, "note": "resected atlas, not an ICI series"},
        {"field": "Tumor % / ABSOLUTE purity", "public": "no", "n": 0, "note": "only public proxy is an RNA epithelial score"},
        {"field": "OS / stage on GEO characteristics", "public": "no", "n": 0, "note": "not deposited as sample characteristics"},
        {"field": "CLDN4 finite (never-smoker tumor)", "public": "yes", "n": int(gene_expr.loc["CLDN4", never].notna().sum()), "note": f"max-mean probe {cldn4_probe}"},
        {"field": "CD8A finite (never-smoker tumor)", "public": "yes", "n": int(gene_expr.loc["CD8A", never].notna().sum()), "note": f"max-mean probe {cd8a_probe}"},
        {"field": "CD274 finite (never-smoker tumor)", "public": "yes", "n": int(gene_expr.loc["CD274", never].notna().sum()), "note": f"max-mean probe {cd274_probe}"},
        {"field": "primary pairwise n (never-smoker LUAD tumor)", "public": "yes", "n": len(never), "note": "this is the n used below"},
        {"field": "ESTIMATE ImmuneSignature genes present (collapse)", "public": "yes", "n": int(cov["present_after_collapse"].sum()), "note": f"of {len(immune_genes)} Yoshihara 2013 genes"},
        {"field": "dual-high TACSTD2×CLDN4 class", "public": "no", "n": 0, "note": "not scored; CLDN4-only"},
    ])
    labels.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    # Figures: never-smoker LUAD tumor primary only.
    prim = pairs_df[pairs_df.cohort == "never_smoker_LUAD_tumor"]
    cldn4_n = gene_expr.loc["CLDN4", never]
    imm_score, _ = score_mean_z(gene_expr.loc[:, never], immune_genes)
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))
    ys = [
        ("CD8A", gene_expr.loc["CD8A", never], "CD8A  RMA log2"),
        ("CD274", gene_expr.loc["CD274", never], "CD274  RMA log2"),
        ("ImmuneScore", imm_score.reindex(never), "ImmuneScore  mean-z"),
    ]
    for ax, (lab, y, ylab) in zip(axes, ys):
        row = prim[prim.pair == f"CLDN4 vs {lab}"].iloc[0]
        ax.scatter(cldn4_n.values, np.asarray(y, float), s=22, alpha=0.8, c="#2c5f8a", edgecolors="none")
        ax.set_xlabel("CLDN4  RMA log2")
        ax.set_ylabel(ylab)
        ax.set_title(f"n={int(row.n)}  ρ={row.rho:.3f}  p={fmt_p(row.p)}")
    fig.suptitle("GSE43458 never-smoker LUAD tumor (n=40)", y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_partners.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGURES / "fig1_cldn4_vs_partners.pdf", bbox_inches="tight")
    plt.close(fig)

    plot = prim[prim.pair.isin(["CLDN4 vs CD8A", "CLDN4 vs CD274", "CLDN4 vs ImmuneScore"])].copy()
    fig, ax = plt.subplots(figsize=(5.6, 2.8))
    y_pos = np.arange(len(plot))
    ax.axvline(0, color="0.5", lw=0.8)
    ax.errorbar(
        plot["rho"], y_pos,
        xerr=[plot["rho"] - plot["ci_low"], plot["ci_high"] - plot["rho"]],
        fmt="o", color="#2c5f8a", ecolor="#2c5f8a", capsize=3, ms=6,
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot["pair"].str.replace("CLDN4 vs ", ""))
    ax.set_xlabel("Spearman ρ vs CLDN4 (bootstrap 95% CI)")
    ax.set_title(f"GSE43458 never-smoker LUAD tumor  n={len(never)}")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cldn4_forest.png", dpi=160)
    fig.savefig(FIGURES / "fig2_cldn4_forest.pdf")
    plt.close(fig)

    c8 = prim[prim.pair == "CLDN4 vs CD8A"].iloc[0]
    pdl1 = prim[prim.pair == "CLDN4 vs CD274"].iloc[0]
    imm = prim[prim.pair == "CLDN4 vs ImmuneScore"].iloc[0]
    n_table = pd.DataFrame([{
        "dataset": "GSE43458 Kadara/Kabbout",
        "histology": "never-smoker LUAD tumor",
        "platform": "GPL6244 HuGene 1.0 ST RMA log2",
        "n_arrays": len(never),
        "CLDN4": cldn4_probe,
        "CD8A": cd8a_probe,
        "CD274": cd274_probe,
        "ImmuneScore": f"ESTIMATE ImmuneSignature mean-z {int(cov['present_after_collapse'].sum())}/141",
        "CLDN4_CD8A_rho": float(c8["rho"]),
        "CLDN4_CD8A_p": float(c8["p"]),
        "CLDN4_CD274_rho": float(pdl1["rho"]),
        "CLDN4_CD274_p": float(pdl1["p"]),
        "CLDN4_ImmuneScore_rho": float(imm["rho"]),
        "CLDN4_ImmuneScore_p": float(imm["p"]),
        "adj_ImmuneScore_rho": float(imm["rho_adj_epithelial"]),
        "adj_ImmuneScore_p": float(imm["p_adj_epithelial"]),
        "verdict": "NO_EVIDENCE",
        "ICI": "none",
    }])
    n_table.to_csv(TABLES / "n_table.tsv", sep="\t", index=False)
    summary = {
        "dataset": "GSE43458",
        "pmid": "23659968",
        "citation": "Kabbout et al., Clin Cancer Res 2013",
        "platform": "GPL6244 Affymetrix Human Gene 1.0 ST [transcript (gene) version]",
        "processing": "RMA log2 as deposited (BRB-Array Tools; Kadara / Kabbout).",
        "histology_primary": "never-smoker lung adenocarcinoma tumor (GEO characteristics)",
        "n_arrays_series": n_array,
        "n_probes": int(probe_expr.shape[0]),
        "n_genes_collapsed": int(gene_expr.shape[0]),
        "n_never_smoker_tumor_primary": len(never),
        "n_smoker_tumor_extra": len(smoker),
        "n_normal_excluded": len(normals),
        "dual_high": False,
        "immunescore": {
            "definition": "mean of gene-wise z-scores of ESTIMATE ImmuneSignature, z-scored within the analysis set",
            "source": "estimate 1.0.13 SI_geneset.gmt ImmuneSignature (Yoshihara 2013)",
            "n_genes_listed": 141,
            "n_genes_present": int(cov["present_after_collapse"].sum()),
            "not": "official R estimate::estimateScore ssGSEA",
        },
        "cldn4_probe": cldn4_probe,
        "cd8a_probe": cd8a_probe,
        "cd274_probe": cd274_probe,
        "tacstd2_probe": tac_probe,
        "primary": {
            "n": len(never),
            "CLDN4_vs_CD8A": {k: jsonable(v) for k, v in c8.to_dict().items()},
            "CLDN4_vs_CD274": {k: jsonable(v) for k, v in pdl1.to_dict().items()},
            "CLDN4_vs_ImmuneScore": {k: jsonable(v) for k, v in imm.to_dict().items()},
        },
        "missing_public_labels": ["ICI", "tumor_percent", "ABSOLUTE purity", "OS", "stage"],
        "cohort_info": infos,
        "series_title": series.get("title", ""),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    print("n_arrays", n_array, "never-smoker tumor", len(never), "smoker tumor", len(smoker), "normal", len(normals))
    print(pairs_df[pairs_df.pair.isin(["CLDN4 vs CD8A", "CLDN4 vs CD274", "CLDN4 vs ImmuneScore"])][
        ["cohort", "pair", "n", "rho", "p", "rho_adj_epithelial", "p_adj_epithelial", "verdict_crude", "verdict_adj"]
    ].to_string(index=False))
    print("[done]", TABLES)


if __name__ == "__main__":
    main()
