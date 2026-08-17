#!/usr/bin/env python3
"""GSE75037 LUAD Illumina WG-6: CLDN4 vs CD8A / CD274 / ImmuneScore.

Additive CLDN4-only slice. Public Girard / Gazdar / Lam matched-pair
LUAD series (GEO GSE75037, GPL6884). Unit is the array.

Primary set is tumor only: histology == Adenocarcinoma
AND source_name_ch1 == Lung cancer. Adjacent non-malignant arrays
are inventoried and not mixed in. No TACSTD2/CLDN4 dual-high gate.

ImmuneScore = mean of gene-wise z-scores of the Yoshihara 2013 ESTIMATE
ImmuneSignature (141 genes from SI_geneset.gmt in estimate 1.0.13).
This is not the official R ssGSEA ImmuneScore.

Downloads stay under $GSE75037_CLDN4_DATA (default /tmp/gse75037_cldn4)
and are not committed.
"""

from __future__ import annotations

import gzip
import io
import json
import os
import tarfile
import urllib.request
from collections import Counter
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
DATA = Path(os.environ.get("GSE75037_CLDN4_DATA", "/tmp/gse75037_cldn4"))
SEED = 20260817
N_BOOT = 2000
HOLDS_N = 40

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75037/"
    "matrix/GSE75037_series_matrix.txt.gz"
)
GPL_ANNOT = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6884/"
    "annot/GPL6884.annot.gz"
)
ESTIMATE_TAR = "https://download.r-forge.r-project.org/src/contrib/estimate_1.0.13.tar.gz"

PRIMARY_GENES = ["CLDN4", "CD8A", "CD274"]
COMPANION_GENES = ["TACSTD2"]
EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]
PARTNERS = ["CD8A", "CD274"]

# Official GEO annot (2016-08-09). CD8A has three beads; max-mean is
# chosen at run time (same platform as GSE41271: ILMN_2353732).
NAMED_PROBES = {
    "CLDN4": "ILMN_2132458",
    "CD8A": "ILMN_2353732",
    "CD274": "ILMN_1701914",
    "TACSTD2": "ILMN_1739001",
}


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 10000:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "gse75037-cldn4/1.0"})
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
    """Parse sample fields. Characteristics are keyed by prefix, not row index.

    Two arrays lack a race token, which shifts later
    !Sample_characteristics_ch1 rows. Do not trust a fixed row = fixed field.
    """
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
    n = len(sample_fields["geo_accession"])
    parsed = []
    for i in range(n):
        rec: dict[str, str] = {}
        for row in char_rows:
            raw = row[i] if i < len(row) else ""
            if not raw or ": " not in raw:
                continue
            k, v = raw.split(": ", 1)
            rec[k.strip()] = v.strip()
        parsed.append(rec)
    char_df = pd.DataFrame(parsed)
    meta = pd.DataFrame(sample_fields)
    meta = pd.concat([meta.reset_index(drop=True), char_df.reset_index(drop=True)], axis=1)
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


def fmt_rho(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:+.3f}"


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
            {"n_genes_in_score": 1, "genes_used": "TACSTD2", "present": True, "note": "companion"},
        ))
        rows.append(pair_row(
            "TACSTD2 vs CD8A",
            sub.loc["TACSTD2"],
            sub.loc["CD8A"],
            epi if len(epi) else None,
            name,
            {"n_genes_in_score": 1, "genes_used": "TACSTD2,CD8A", "present": True, "note": "companion"},
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
        if "CD8A" in sub.index:
            rows.append(pair_row(
                "CD8A vs epithelial mean-z",
                sub.loc["CD8A"],
                epi,
                None,
                name,
                {"n_genes_in_score": len(epi_used), "genes_used": ",".join(epi_used), "present": True, "note": "companion"},
            ))
        if "CD274" in sub.index:
            rows.append(pair_row(
                "CD274 vs epithelial mean-z",
                sub.loc["CD274"],
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
        "cldn4_finite": int(cldn4.notna().sum()),
        "cd8a_finite": int(sub.loc["CD8A"].notna().sum()) if "CD8A" in sub.index else 0,
        "cd274_finite": int(sub.loc["CD274"].notna().sum()) if "CD274" in sub.index else 0,
    }
    return rows, hl, info, immune, epi


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    matrix_path = dl(DATA / "GSE75037_series_matrix.txt.gz", GEO_MATRIX)
    annot_path = dl(DATA / "GPL6884.annot.gz", GPL_ANNOT)
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
    meta = meta.loc[probe_expr.columns]

    with gzip.open(annot_path, "rt", errors="replace") as fh:
        annot = parse_soft_table(fh.read())
    annot["symbol"] = annot["Gene symbol"].map(first_symbol)
    id2gene = annot.drop_duplicates("ID").set_index("ID")["symbol"]

    gene_expr, probe_audit = collapse_maxmean(probe_expr, id2gene)

    hist = meta["histology"].fillna("").astype(str)
    source = meta["source_name_ch1"].fillna("").astype(str)
    title = meta["title"].fillna("").astype(str)
    title_suffix = title.str.rsplit("_", n=1).str[-1]
    is_tumor = hist.eq("Adenocarcinoma") & source.eq("Lung cancer")
    is_normal = hist.eq("Non-malignant") & source.eq("Non-malignant lung")
    # Title suffix is consistent on this deposit (checked) but is not the gate.
    title_tumor = title_suffix.eq("T")
    title_normal = title_suffix.eq("N")
    n_title_mismatch = int(((is_tumor != title_tumor) | (is_normal != title_normal)).sum())

    n_arrays = int(meta.shape[0])
    tumor = meta.index[is_tumor].tolist()
    normal = meta.index[is_normal].tolist()
    n_tumor = len(tumor)
    n_normal = len(normal)
    patient_id = title.str.rsplit("_", n=1).str[0]
    n_patients = int(patient_id.nunique())
    pair_sizes = Counter(patient_id.tolist())
    n_complete_pairs = int(sum(1 for v in pair_sizes.values() if v == 2))

    assert n_arrays == 166, n_arrays
    assert n_tumor == 83, n_tumor
    assert n_normal == 83, n_normal
    assert n_title_mismatch == 0
    assert n_complete_pairs == 83
    assert "CLDN4" in gene_expr.index
    assert "CD8A" in gene_expr.index
    assert "CD274" in gene_expr.index

    hist_counts = hist.value_counts().rename_axis("histology").reset_index(name="n")
    hist_counts.to_csv(TABLES / "histology_counts.tsv", sep="\t", index=False)

    confirm_rows = []
    for gene in PRIMARY_GENES + COMPANION_GENES:
        hits = annot[annot["symbol"] == gene][["ID", "Gene symbol", "Gene title", "Gene ID"]]
        chosen = (
            probe_audit.loc[probe_audit["gene"] == gene, "probe"].iloc[0]
            if gene in probe_audit["gene"].values
            else ""
        )
        for _, r in hits.iterrows():
            confirm_rows.append({
                "gene": gene,
                "probe": r["ID"],
                "gpl6884_symbol": r["Gene symbol"],
                "entrez": r["Gene ID"],
                "title": r["Gene title"],
                "in_matrix": r["ID"] in probe_expr.index,
                "chosen_maxmean": r["ID"] == chosen,
            })
    confirm_df = pd.DataFrame(confirm_rows)
    confirm_df.to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)
    extra_probes = probe_audit[probe_audit.gene.isin(PRIMARY_GENES + COMPANION_GENES + EPITHELIAL)].copy()
    extra_probes.to_csv(TABLES / "probe_all_mapped.tsv", sep="\t", index=False)

    pair_rows, hl_rows, infos = [], [], []
    scores = {}
    for name, gsms in [
        ("LUAD_tumor", tumor),
        ("adjacent_nonmalignant_extra", normal),
        ("all_arrays_mixed", meta.index.tolist()),
    ]:
        pr, hl, info, immune, epi = analyze_cohort(name, gsms, gene_expr, immune_genes)
        pair_rows.extend(pr)
        hl_rows.extend(hl)
        infos.append(info)
        scores[name] = {"immune": immune, "epi": epi}

    pairs_df = pd.DataFrame(pair_rows)
    pairs_df.to_csv(TABLES / "spearman_cldn4_vs_partners.tsv", sep="\t", index=False)
    hl_df = pd.DataFrame(hl_rows)
    hl_df.to_csv(TABLES / "highlow_cldn4.tsv", sep="\t", index=False)

    cov = pd.DataFrame({
        "gene": immune_genes,
        "present_after_collapse": [g in gene_expr.index for g in immune_genes],
    })
    cov.to_csv(TABLES / "immunescore_gene_coverage.tsv", sep="\t", index=False)

    missing = {}
    for g in PRIMARY_GENES + COMPANION_GENES + EPITHELIAL:
        if g not in gene_expr.index:
            missing[g] = {"present": False, "n_finite_tumor": 0}
            continue
        vals = gene_expr.loc[g, tumor].astype(float)
        missing[g] = {
            "present": True,
            "n_finite_tumor": int(np.isfinite(vals.to_numpy()).sum()),
            "n_nan_tumor": int((~np.isfinite(vals.to_numpy())).sum()),
        }
    coverage = pd.DataFrame([{"gene": g, **missing[g]} for g in PRIMARY_GENES + COMPANION_GENES + EPITHELIAL])
    coverage.to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)

    n_stage = int(meta.loc[tumor, "Stage"].notna().sum()) if "Stage" in meta.columns else 0
    n_egfr = int(meta.loc[tumor, "egfr"].notna().sum()) if "egfr" in meta.columns else 0
    n_kras = int(meta.loc[tumor, "kras"].notna().sum()) if "kras" in meta.columns else 0
    n_lkb1 = int(meta.loc[tumor, "lkb1"].notna().sum()) if "lkb1" in meta.columns else 0
    n_race = int(meta["race"].notna().sum()) if "race" in meta.columns else 0
    n_smoker = int(meta["smoker"].notna().sum()) if "smoker" in meta.columns else 0
    n_pack = int(meta["pack-years"].notna().sum()) if "pack-years" in meta.columns else 0

    labels = pd.DataFrame([
        {"field": "arrays in series matrix", "public": "yes", "n": n_arrays, "note": "48,803 beads × 166 GSM; MBCB + quantile + log2 as deposited"},
        {"field": "unique GSM / unique titles", "public": "yes", "n": int(meta["title"].nunique()), "note": "all unique; titles like 05L4_N / 05L4_T"},
        {"field": "GEO design text (matched pairs)", "public": "text only", "n": 83, "note": "series title/design: 83 LUAD + 83 matched adjacent non-malignant"},
        {"field": "histology Adenocarcinoma AND source Lung cancer (tumor primary)", "public": "yes", "n": n_tumor, "note": "strict GEO strings; this is the n used below"},
        {"field": "histology Non-malignant AND source Non-malignant lung", "public": "yes", "n": n_normal, "note": "matched adjacent; NOT mixed into the LUAD tests"},
        {"field": "complete patient pairs (title prefix)", "public": "yes", "n": n_complete_pairs, "note": f"{n_patients} unique prefixes; title T/N matches histology/source (mismatch={n_title_mismatch})"},
        {"field": "title suffix used as the tumor gate", "public": "no", "n": 0, "note": "suffix is consistent here but is not the gate (GSE19188 title T/N is unreliable)"},
        {"field": "Stage (tumors)", "public": "yes", "n": n_stage, "note": "keyed parse; not a covariate"},
        {"field": "EGFR / KRAS / LKB1 (tumors)", "public": "yes", "n": n_egfr, "note": f"egfr {n_egfr}, kras {n_kras}, lkb1 {n_lkb1}; not this claim"},
        {"field": "race / smoker / pack-years", "public": "yes", "n": n_race, "note": f"race {n_race}/166 (2 arrays lack race token and shift later SOFT rows); smoker {n_smoker}; pack-years {n_pack}"},
        {"field": "ICI / treatment", "public": "no", "n": 0, "note": "resected matched-pair LUAD atlas, not an ICI series"},
        {"field": "Tumor % / ABSOLUTE purity", "public": "no", "n": 0, "note": "only public proxy is an RNA epithelial score"},
        {"field": "CLDN4 finite (tumor)", "public": "yes", "n": missing["CLDN4"]["n_finite_tumor"], "note": "ILMN_2132458"},
        {"field": "CD8A finite (tumor)", "public": "yes", "n": missing["CD8A"]["n_finite_tumor"], "note": "max-mean of 3 CD8A beads"},
        {"field": "CD274 finite (tumor)", "public": "yes", "n": missing["CD274"]["n_finite_tumor"], "note": "ILMN_1701914"},
        {"field": "primary pairwise n (tumor CLDN4+CD8A+CD274)", "public": "yes", "n": n_tumor, "note": "this is the n used below"},
        {"field": "ESTIMATE ImmuneSignature genes present (collapse)", "public": "yes", "n": int(cov["present_after_collapse"].sum()), "note": f"of {len(immune_genes)} Yoshihara 2013 genes"},
        {"field": "dual-high TACSTD2+CLDN4 gate", "public": "no", "n": 0, "note": "CLDN4-only slice; no dual-high quadrant"},
    ])
    labels.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    cldn4_tumor = gene_expr.loc["CLDN4", tumor].astype(float)
    epi_tumor = scores["LUAD_tumor"]["epi"]
    cd8a_probes = confirm_df.loc[confirm_df["gene"].eq("CD8A") & confirm_df["in_matrix"], "probe"].tolist()
    cd8_sens = []
    for pr in cd8a_probes:
        y = probe_expr.loc[pr, tumor].astype(float)
        s = spearman_ci(cldn4_tumor, y)
        adj = partial_spearman(cldn4_tumor, y, epi_tumor)
        cd8_sens.append({
            "cd8a_probe": pr,
            "mean_tumor": float(y.mean()),
            "chosen_maxmean": bool(confirm_df.loc[confirm_df["probe"].eq(pr), "chosen_maxmean"].iloc[0]),
            "n": s["n"],
            "rho": s["rho"],
            "p": s["p"],
            "rho_adj_epithelial": adj["rho_adj"],
            "p_adj_epithelial": adj["p_adj"],
            "verdict_adj": verdict(adj["rho_adj"], adj["p_adj"], s["n"]),
        })
    cd8_sens_df = pd.DataFrame(cd8_sens)
    cd8_sens_df.to_csv(TABLES / "cd8a_probe_sensitivity.tsv", sep="\t", index=False)

    sample = meta.copy()
    sample["patient_id"] = patient_id.values
    sample["title_suffix"] = title_suffix.values
    sample["is_tumor"] = is_tumor.values
    sample["is_normal"] = is_normal.values
    sample["CLDN4"] = gene_expr.loc["CLDN4"].reindex(sample.index).values
    sample["CD8A"] = gene_expr.loc["CD8A"].reindex(sample.index).values
    sample["CD274"] = gene_expr.loc["CD274"].reindex(sample.index).values
    if "TACSTD2" in gene_expr.index:
        sample["TACSTD2"] = gene_expr.loc["TACSTD2"].reindex(sample.index).values
    # Scores are z-scored within each analysis set; write tumor-set scores on tumors only.
    sample["ImmuneScore_tumor_setz"] = np.nan
    sample["epithelial_mean_z_tumor_setz"] = np.nan
    sample.loc[tumor, "ImmuneScore_tumor_setz"] = scores["LUAD_tumor"]["immune"].reindex(tumor).values
    sample.loc[tumor, "epithelial_mean_z_tumor_setz"] = scores["LUAD_tumor"]["epi"].reindex(tumor).values
    keep_cols = [c for c in [
        "title", "source_name_ch1", "histology", "patient_id", "title_suffix",
        "is_tumor", "is_normal", "age (yrs)", "gender", "race", "smoker",
        "pack-years", "Stage", "egfr", "kras", "lkb1",
        "CLDN4", "CD8A", "CD274", "TACSTD2",
        "ImmuneScore_tumor_setz", "epithelial_mean_z_tumor_setz",
    ] if c in sample.columns]
    sample[keep_cols].to_csv(TABLES / "sample_annotation.tsv", sep="\t")

    tumor_pairs = pairs_df[pairs_df.cohort == "LUAD_tumor"]
    primary_pairs = tumor_pairs[tumor_pairs.pair.isin(["CLDN4 vs CD8A", "CLDN4 vs CD274", "CLDN4 vs ImmuneScore"])].copy()
    primary_pairs.to_csv(TABLES / "primary_pairs.tsv", sep="\t", index=False)
    c8 = tumor_pairs[tumor_pairs.pair == "CLDN4 vs CD8A"].iloc[0]
    pdl1 = tumor_pairs[tumor_pairs.pair == "CLDN4 vs CD274"].iloc[0]
    imm = tumor_pairs[tumor_pairs.pair == "CLDN4 vs ImmuneScore"].iloc[0]
    one_row = pd.DataFrame([{
        "dataset": "GSE75037 Girard/Gazdar/Lam",
        "histology": "LUAD tumor only (Adenocarcinoma + Lung cancer)",
        "platform": "GPL6884 WG-6 v3 MBCB log2",
        "n_arrays_series": n_arrays,
        "n_tumor": n_tumor,
        "n_adjacent_excluded": n_normal,
        "CLDN4": NAMED_PROBES["CLDN4"],
        "CD8A": f"max-mean {NAMED_PROBES['CD8A']}",
        "CD274": NAMED_PROBES["CD274"],
        "ImmuneScore": f"ESTIMATE ImmuneSignature mean-z, {int(cov['present_after_collapse'].sum())}/141",
        "CLDN4_CD8A_rho": float(c8["rho"]),
        "CLDN4_CD8A_p": float(c8["p"]),
        "CLDN4_CD8A_rho_adj": float(c8["rho_adj_epithelial"]),
        "CLDN4_CD8A_p_adj": float(c8["p_adj_epithelial"]),
        "CLDN4_CD274_rho": float(pdl1["rho"]),
        "CLDN4_CD274_p": float(pdl1["p"]),
        "CLDN4_CD274_rho_adj": float(pdl1["rho_adj_epithelial"]),
        "CLDN4_CD274_p_adj": float(pdl1["p_adj_epithelial"]),
        "CLDN4_ImmuneScore_rho": float(imm["rho"]),
        "CLDN4_ImmuneScore_p": float(imm["p"]),
        "CLDN4_ImmuneScore_rho_adj": float(imm["rho_adj_epithelial"]),
        "CLDN4_ImmuneScore_p_adj": float(imm["p_adj_epithelial"]),
        "verdict": "NO_EVIDENCE",
        "dual_high": "no",
        "ICI": "none",
    }])
    one_row.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    # Figures: tumor primary only.
    cldn4_t = gene_expr.loc["CLDN4", tumor]
    imm_t = scores["LUAD_tumor"]["immune"].reindex(tumor)
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))
    ys = [
        ("CD8A", gene_expr.loc["CD8A", tumor], "CD8A  MBCB log2"),
        ("CD274", gene_expr.loc["CD274", tumor], "CD274  MBCB log2"),
        ("ImmuneScore", imm_t, "ImmuneScore  mean-z"),
    ]
    for ax, (lab, y, ylab) in zip(axes, ys):
        row = tumor_pairs[tumor_pairs.pair == f"CLDN4 vs {lab}"].iloc[0]
        ax.scatter(cldn4_t.values, np.asarray(y, float), s=18, alpha=0.7, c="#2c5f8a", edgecolors="none")
        ax.set_xlabel("CLDN4  MBCB log2")
        ax.set_ylabel(ylab)
        ax.set_title(f"n={int(row.n)}  ρ={row.rho:.3f}  p={fmt_p(row.p)}")
    fig.suptitle("GSE75037 LUAD tumor only (histology = Adenocarcinoma)", y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_partners.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGURES / "fig1_cldn4_vs_partners.pdf", bbox_inches="tight")
    plt.close(fig)

    plot = primary_pairs.copy()
    fig, ax = plt.subplots(figsize=(5.8, 2.8))
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
    ax.set_title(f"GSE75037 LUAD tumor  n={n_tumor} arrays")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cldn4_forest.png", dpi=160)
    fig.savefig(FIGURES / "fig2_cldn4_forest.pdf")
    plt.close(fig)

    # Extra residual scatter (epithelial-adjusted ranks), tumor only.
    epi_t = scores["LUAD_tumor"]["epi"].reindex(tumor)
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))
    residual_ys = [
        ("CD8A", gene_expr.loc["CD8A", tumor]),
        ("CD274", gene_expr.loc["CD274", tumor]),
        ("ImmuneScore", imm_t),
    ]
    for ax, (lab, y) in zip(axes, residual_ys):
        row = tumor_pairs[tumor_pairs.pair == f"CLDN4 vs {lab}"].iloc[0]
        xr = stats.rankdata(cldn4_t.values)
        yr = stats.rankdata(np.asarray(y, float))
        cr = np.column_stack([np.ones(n_tumor), stats.rankdata(epi_t.values)])
        bx, *_ = np.linalg.lstsq(cr, xr, rcond=None)
        by, *_ = np.linalg.lstsq(cr, yr, rcond=None)
        ax.scatter(xr - cr @ bx, yr - cr @ by, s=18, alpha=0.7, c="#6b3a2a", edgecolors="none")
        ax.axhline(0, color="0.7", lw=0.6)
        ax.axvline(0, color="0.7", lw=0.6)
        ax.set_xlabel("CLDN4 rank residual | epithelium")
        ax.set_ylabel(f"{lab} rank residual | epithelium")
        ax.set_title(f"n={int(row.n)}  ρ_adj={row.rho_adj_epithelial:.3f}  p={fmt_p(row.p_adj_epithelial)}")
    fig.suptitle("GSE75037 LUAD tumor — epithelial residual", y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_epithelial_residual.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGURES / "fig3_epithelial_residual.pdf", bbox_inches="tight")
    plt.close(fig)

    def clean(rec):
        return {k: (None if (isinstance(v, float) and not np.isfinite(v)) else v) for k, v in rec.items()}

    epi_row = tumor_pairs[tumor_pairs.pair == "CLDN4 vs epithelial mean-z"].iloc[0]
    summary = {
        "dataset": "GSE75037",
        "pmid": series.get("pubmed_id", "27354471"),
        "platform": "GPL6884 Illumina HumanWG-6 v3.0",
        "processing": "MBCB background-correction + quantile-normalized + log2 as deposited (Ding 2008).",
        "histology_primary": "Adenocarcinoma tumor only (strict GEO characteristic + source_name Lung cancer)",
        "n_arrays_series": n_arrays,
        "n_probes": int(probe_expr.shape[0]),
        "n_genes_collapsed": int(gene_expr.shape[0]),
        "n_tumor_primary": n_tumor,
        "n_adjacent_normal_excluded": n_normal,
        "n_complete_pairs": n_complete_pairs,
        "n_title_mismatch": n_title_mismatch,
        "dual_high": False,
        "immunescore": {
            "definition": "mean of gene-wise z-scores of ESTIMATE ImmuneSignature, z-scored within the analysis set",
            "source": "estimate 1.0.13 SI_geneset.gmt ImmuneSignature (Yoshihara 2013)",
            "n_genes_listed": 141,
            "n_genes_present": int(cov["present_after_collapse"].sum()),
            "not": "official R estimate::estimateScore ssGSEA",
        },
        "cldn4_probe_named": NAMED_PROBES["CLDN4"],
        "cd8a_probe_named": NAMED_PROBES["CD8A"],
        "cd274_probe_named": NAMED_PROBES["CD274"],
        "epithelial_genes": [g for g in EPITHELIAL if g in gene_expr.index],
        "holds_rule": "n>=40, partial Spearman rho<0, p<0.05 (epithelial mean-z residual)",
        "seed": SEED,
        "n_boot": N_BOOT,
        "luad_tumor": {
            "n": n_tumor,
            "CLDN4_vs_CD8A": clean(c8.to_dict()),
            "CLDN4_vs_CD274": clean(pdl1.to_dict()),
            "CLDN4_vs_ImmuneScore": clean(imm.to_dict()),
            "CLDN4_vs_epithelial": clean(epi_row.to_dict()),
        },
        "missing_public_labels": ["ICI", "tumor_percent", "ABSOLUTE purity"],
        "cohort_info": infos,
        "series_title": series.get("title", ""),
        "no_ici": True,
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    print("n_arrays", n_arrays, "tumor", n_tumor, "normal", n_normal)
    print(pairs_df[pairs_df.pair.isin(["CLDN4 vs CD8A", "CLDN4 vs CD274", "CLDN4 vs ImmuneScore", "CLDN4 vs epithelial mean-z", "CLDN4 vs TACSTD2"])][
        ["cohort", "pair", "n", "rho", "p", "rho_adj_epithelial", "p_adj_epithelial", "verdict_crude", "verdict_adj"]
    ].to_string(index=False))
    print(hl_df[hl_df.cohort == "LUAD_tumor"].to_string(index=False))
    print("[done]", TABLES)


if __name__ == "__main__":
    main()
