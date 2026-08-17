#!/usr/bin/env python3
"""GSE41271 LUAD Illumina WG-6: CLDN4 vs CD8A / CD274.

Public MD Anderson / Girard NSCLC series (GEO GSE41271, GPL6884).
Primary slice is LUAD arrays labelled histology: Adenocarcinoma.
Downloads stay under $GSE41271_CLDN4_DATA (default /tmp/gse41271_cldn4)
and are not committed.
"""

from __future__ import annotations

import gzip
import io
import json
import os
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
DATA = Path(os.environ.get("GSE41271_CLDN4_DATA", "/tmp/gse41271_cldn4"))
SEED = 20260817
N_BOOT = 2000
HOLDS_N = 40

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41271/"
    "matrix/GSE41271_series_matrix.txt.gz"
)
GPL_ANNOT = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6884/"
    "annot/GPL6884.annot.gz"
)

# Primary pairs requested. TACSTD2 is a same-run companion only.
PRIMARY_GENES = ["CLDN4", "CD8A", "CD274"]
COMPANION_GENES = ["TACSTD2"]
EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]

# Official GEO annot (2016-08-09). CD8A has three beads; max-mean is
# chosen at run time (this matrix: ILMN_2353732).
NAMED_PROBES = {
    "CLDN4": "ILMN_2132458",
    "CD8A": "ILMN_2353732",
    "CD274": "ILMN_1701914",
    "TACSTD2": "ILMN_1739001",
}


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "gse41271-cldn4/1.0"})
    with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as out:
        out.write(r.read())
    return dest


def parse_soft_table(text: str) -> pd.DataFrame:
    lines = text.splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_begin")) + 1
    end = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_end"))
    return pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", dtype=str)


def parse_series_meta(path: Path) -> tuple[pd.DataFrame, dict]:
    """Parse sample fields. Characteristics are keyed by prefix, not row index.

    Three arrays have a missing tobacco-history token, which shifts later
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


def scatter(x, y, xlabel, ylabel, title, out: Path, n: int, rho: float, p: float) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    ax.scatter(x, y, s=18, alpha=0.65, c="#1f4e79", edgecolors="none")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.text(
        0.03,
        0.97,
        f"n={n}\nρ={rho:+.3f}\np={fmt_p(p)}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="#ccc"),
    )
    fig.tight_layout()
    fig.savefig(out.with_suffix(".png"), dpi=160)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    matrix_path = dl(DATA / "GSE41271_series_matrix.txt.gz", GEO_MATRIX)
    annot_path = dl(DATA / "GPL6884.annot.gz", GPL_ANNOT)

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
    is_luad = hist.eq("Adenocarcinoma")
    n_arrays = int(meta.shape[0])
    n_luad = int(is_luad.sum())
    hist_counts = Counter(hist.tolist())

    # Finite-value n for the three requested genes on LUAD.
    luad_cols = meta.index[is_luad]
    missing = {}
    for g in PRIMARY_GENES + COMPANION_GENES + EPITHELIAL:
        if g not in gene_expr.index:
            missing[g] = {"present": False, "n_finite_luad": 0}
            continue
        vals = gene_expr.loc[g, luad_cols].astype(float)
        missing[g] = {
            "present": True,
            "n_finite_luad": int(np.isfinite(vals.to_numpy()).sum()),
            "n_nan_luad": int((~np.isfinite(vals.to_numpy())).sum()),
        }

    epi_present = [g for g in EPITHELIAL if g in gene_expr.index]
    if len(epi_present) < 3:
        raise RuntimeError(f"epithelial genes missing: {epi_present}")
    z = gene_expr.loc[epi_present, luad_cols].astype(float)
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    epi_score = z.mean(axis=0)

    cldn4 = gene_expr.loc["CLDN4", luad_cols].astype(float)
    cd8a = gene_expr.loc["CD8A", luad_cols].astype(float)
    cd274 = gene_expr.loc["CD274", luad_cols].astype(float)
    tacstd2 = gene_expr.loc["TACSTD2", luad_cols].astype(float) if "TACSTD2" in gene_expr.index else None

    pairs = [
        ("CLDN4", "CD8A", cldn4, cd8a),
        ("CLDN4", "CD274", cldn4, cd274),
    ]
    if tacstd2 is not None:
        pairs.append(("CLDN4", "TACSTD2", cldn4, tacstd2))
        pairs.append(("TACSTD2", "CD8A", tacstd2, cd8a))
        pairs.append(("TACSTD2", "CD274", tacstd2, cd274))
    pairs.append(("CLDN4", "epithelial_mean_z", cldn4, epi_score))
    pairs.append(("CD8A", "epithelial_mean_z", cd8a, epi_score))
    pairs.append(("CD274", "epithelial_mean_z", cd274, epi_score))

    rows = []
    for a, b, x, y in pairs:
        s = spearman_ci(x, y)
        adj = partial_spearman(x, y, epi_score) if b != "epithelial_mean_z" else dict(n=s["n"], rho_adj=np.nan, p_adj=np.nan)
        v = verdict(adj["rho_adj"], adj["p_adj"], s["n"]) if b != "epithelial_mean_z" else "—"
        rows.append({
            "pair": f"{a} vs {b}",
            "predictor": a,
            "endpoint": b,
            "subset": "LUAD Adenocarcinoma",
            "n": s["n"],
            "rho": s["rho"],
            "p": s["p"],
            "ci_low": s["ci_low"],
            "ci_high": s["ci_high"],
            "rho_adj_epithelial": adj["rho_adj"],
            "p_adj_epithelial": adj["p_adj"],
            "verdict": v,
        })
    stats_df = pd.DataFrame(rows)

    # Q4 vs Q1 on CD8A and CD274.
    q = cldn4.quantile([0.25, 0.75])
    lo = cldn4 <= q.iloc[0]
    hi = cldn4 >= q.iloc[1]
    hl_rows = []
    for name, y in [("CD8A", cd8a), ("CD274", cd274)]:
        a = y[hi].to_numpy(float)
        b = y[lo].to_numpy(float)
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        # rank-biserial: 1 - 2U/(n1 n2), sign so Q4>Q1 is positive
        rrb = (2 * u) / (len(a) * len(b)) - 1
        hl_rows.append({
            "endpoint": name,
            "n_q4": int(a.size),
            "n_q1": int(b.size),
            "median_q4": float(np.median(a)),
            "median_q1": float(np.median(b)),
            "mwu_U": float(u),
            "mwu_p": float(p),
            "rank_biserial_q4_minus_q1": float(rrb),
        })
    hl_df = pd.DataFrame(hl_rows)

    # Probe confirmation.
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

    # CD8A has three beads. Report each vs CLDN4 so max-mean is not hidden.
    cd8a_probes = confirm_df.loc[confirm_df["gene"].eq("CD8A") & confirm_df["in_matrix"], "probe"].tolist()
    cd8_sens = []
    for pr in cd8a_probes:
        y = probe_expr.loc[pr, luad_cols].astype(float)
        s = spearman_ci(cldn4, y)
        adj = partial_spearman(cldn4, y, epi_score)
        cd8_sens.append({
            "cd8a_probe": pr,
            "mean_luad": float(y.mean()),
            "chosen_maxmean": bool(confirm_df.loc[confirm_df["probe"].eq(pr), "chosen_maxmean"].iloc[0]),
            "n": s["n"],
            "rho": s["rho"],
            "p": s["p"],
            "rho_adj_epithelial": adj["rho_adj"],
            "p_adj_epithelial": adj["p_adj"],
            "verdict": verdict(adj["rho_adj"], adj["p_adj"], s["n"]),
        })
    cd8_sens_df = pd.DataFrame(cd8_sens)

    # Honest-n inventory.
    n_tobacco = int(meta["tobacco history"].notna().sum()) if "tobacco history" in meta.columns else 0
    n_stage = int(meta["final patient stage"].notna().sum()) if "final patient stage" in meta.columns else 0
    n_vital = int(meta["vital statistics"].isin(["A", "D"]).sum()) if "vital statistics" in meta.columns else 0
    n_rec = int(meta["recurrence"].isin(["Y", "N"]).sum()) if "recurrence" in meta.columns else 0
    inventory = pd.DataFrame([
        {"item": "arrays in series matrix", "public": "yes", "n": n_arrays, "note": "275 GSM × 48,803 Illumina beads"},
        {"item": "unique GSM / unique titles", "public": "yes", "n": int(meta["title"].nunique()), "note": "all unique; titles like 15-T"},
        {"item": "histology Adenocarcinoma (LUAD)", "public": "yes", "n": n_luad, "note": "strict label; not adenosquamous / sarcomatoid-adeno"},
        {"item": "histology Squamous", "public": "yes", "n": int(hist.eq("Squamous").sum()), "note": "dropped from primary"},
        {"item": "other histology", "public": "yes", "n": int((~hist.isin(["Adenocarcinoma", "Squamous"])).sum()), "note": "LCC-NE / adenosquamous / sarcomatoid / NSCLC / mixed"},
        {"item": "GEO design text ADC", "public": "text only", "n": 183, "note": "series summary: mainly adenocarcinomas (n = 183)"},
        {"item": "CLDN4 finite on LUAD", "public": "yes", "n": missing["CLDN4"]["n_finite_luad"], "note": "ILMN_2132458"},
        {"item": "CD8A finite on LUAD", "public": "yes", "n": missing["CD8A"]["n_finite_luad"], "note": "max-mean of 3 CD8A beads"},
        {"item": "CD274 finite on LUAD", "public": "yes", "n": missing["CD274"]["n_finite_luad"], "note": "ILMN_1701914"},
        {"item": "primary pairwise n (LUAD, CLDN4+CD8A+CD274)", "public": "yes", "n": n_luad, "note": "this is the n used below"},
        {"item": "tobacco history parsed", "public": "yes", "n": n_tobacco, "note": "3 arrays lack this token; later SOFT rows shift"},
        {"item": "final patient stage parsed", "public": "yes", "n": n_stage, "note": "not used in the correlations"},
        {"item": "vital A/D", "public": "yes", "n": n_vital, "note": "dates present; no OS time computed here"},
        {"item": "recurrence Y/N", "public": "yes", "n": n_rec, "note": "not an ICI series"},
        {"item": "ICI / PD-1 treatment", "public": "no", "n": 0, "note": "surgical 1997–2005 MD Anderson tumors"},
    ])

    per_sample = meta.loc[luad_cols, [c for c in ["title", "histology", "gender", "race", "tobacco history", "final patient stage", "vital statistics", "recurrence"] if c in meta.columns]].copy()
    per_sample["CLDN4"] = cldn4
    per_sample["CD8A"] = cd8a
    per_sample["CD274"] = cd274
    if tacstd2 is not None:
        per_sample["TACSTD2"] = tacstd2
    per_sample["epithelial_mean_z"] = epi_score

    one_row = stats_df[stats_df["endpoint"].isin(["CD8A", "CD274"]) & stats_df["predictor"].eq("CLDN4")].copy()

    inventory.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)
    stats_df.to_csv(TABLES / "spearman_cldn4_vs_cd8a_cd274.tsv", sep="\t", index=False)
    one_row.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)
    hl_df.to_csv(TABLES / "highlow_cldn4.tsv", sep="\t", index=False)
    confirm_df.to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)
    cd8_sens_df.to_csv(TABLES / "cd8a_probe_sensitivity.tsv", sep="\t", index=False)
    per_sample.to_csv(TABLES / "sample_annotation.tsv", sep="\t")
    pd.DataFrame(
        [{"histology": k, "n": v} for k, v in sorted(hist_counts.items(), key=lambda kv: (-kv[1], kv[0]))]
    ).to_csv(TABLES / "histology_counts.tsv", sep="\t", index=False)
    coverage = pd.DataFrame([
        {"gene": g, **missing[g]} for g in PRIMARY_GENES + COMPANION_GENES + EPITHELIAL
    ])
    coverage.to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)

    primary = {r["endpoint"]: r for r in one_row.to_dict(orient="records")}
    summary = {
        "dataset": "GSE41271",
        "platform": "GPL6884 Illumina HumanWG-6 v3.0",
        "n_arrays": n_arrays,
        "n_luad_adenocarcinoma": n_luad,
        "histology_counts": dict(hist_counts),
        "primary_n": n_luad,
        "seed": SEED,
        "n_boot": N_BOOT,
        "epithelial_genes": epi_present,
        "holds_rule": "n>=40, partial Spearman rho<0, p<0.05 (epithelial mean-z residual)",
        "cldn4_vs_cd8a": primary.get("CD8A"),
        "cldn4_vs_cd274": primary.get("CD274"),
        "q4q1": hl_df.to_dict(orient="records"),
        "cd8a_probe_sensitivity": cd8_sens_df.to_dict(orient="records"),
        "series_title": series.get("title", ""),
        "pubmed": series.get("pubmed_id", ""),
        "no_ici": True,
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    scatter(
        cldn4, cd8a,
        "CLDN4 (max-mean, deposited)",
        "CD8A (max-mean, deposited)",
        "GSE41271 LUAD — CLDN4 vs CD8A",
        FIGURES / "fig1_cldn4_vs_cd8a",
        n_luad, float(primary["CD8A"]["rho"]), float(primary["CD8A"]["p"]),
    )
    scatter(
        cldn4, cd274,
        "CLDN4 (max-mean, deposited)",
        "CD274 (ILMN_1701914, deposited)",
        "GSE41271 LUAD — CLDN4 vs CD274",
        FIGURES / "fig2_cldn4_vs_cd274",
        n_luad, float(primary["CD274"]["rho"]), float(primary["CD274"]["p"]),
    )

    print("n_arrays", n_arrays, "n_luad", n_luad)
    print(stats_df.to_string(index=False))
    print(hl_df.to_string(index=False))
    print("histology", dict(hist_counts))


if __name__ == "__main__":
    main()
