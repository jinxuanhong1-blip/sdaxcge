#!/usr/bin/env python3
"""Additive GSE182328 extra cuts: CLDN4 vs CD274/HLA and Q4 vs Q1 CD8.

Continuous CLDN4 vs CD8A is already in methods/lung_ici_bulk_cldn4_leftover
(ρ=−0.40, p=0.0067, n=44). This script does not re-claim that row.

Public GEO processed counts only. Patient = sample. Honest pairwise-complete n.
"""
from __future__ import annotations

import gzip
import json
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
CACHE = Path("/tmp/gse182328_cldn4_extra")
FIG = HERE / "figures"
TAB = HERE / "tables"
for d in (CACHE, FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

COUNTS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182328/"
    "suppl/GSE182328_Gene_counts_matrix.txt.gz"
)
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182328/"
    "matrix/GSE182328_series_matrix.txt.gz"
)

# Extra-cut genes. CD8A is scored only for the Q4 vs Q1 cut and the known-row reprint.
HLA_I = ["HLA-A", "HLA-B", "HLA-C"]
HLA_EXTRA = ["HLA-E", "HLA-F", "HLA-G", "HLA-DRA", "B2M"]
GENES = ["CLDN4", "CD8A", "CD8B", "CD274", "TACSTD2"] + HLA_I + HLA_EXTRA


def fetch(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    return dest


def parse_series(path: Path) -> pd.DataFrame:
    stored: dict[str, list[str]] = {}
    char_rows: list[list[str]] = []

    def fields(line: str) -> list[str]:
        return [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]

    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_geo_accession"):
                stored["gsm"] = fields(line)
            elif line.startswith("!Sample_title"):
                stored["title"] = fields(line)
            elif line.startswith("!Sample_source_name"):
                stored["source"] = fields(line)
            elif line.startswith("!Sample_characteristics"):
                char_rows.append(fields(line))
    geo = stored.get("gsm") or []
    n = len(geo)
    meta = pd.DataFrame({"gsm": geo})
    for key in ("title", "source"):
        vals = stored.get(key) or [None] * n
        meta[key] = (vals + [None] * n)[:n]
    chars: dict[str, dict[str, str]] = {g: {} for g in geo}
    for vals in char_rows:
        for g, v in zip(geo, vals):
            if not v:
                continue
            if ":" in v:
                k, val = v.split(":", 1)
                chars[g][k.strip().lower()] = val.strip()
            else:
                chars[g]["characteristic"] = v
    keys = sorted({k for d in chars.values() for k in d})
    for k in keys:
        meta[k] = [chars[g].get(k) for g in geo]
    return meta


def log2cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts.div(lib, axis=1) * 1e6 + 1)


def mean_z(log: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    found = [g for g in genes if g in log.index]
    arr = log.loc[found].to_numpy(dtype=float)
    mu = np.nanmean(arr, axis=1, keepdims=True)
    sd = np.nanstd(arr, axis=1, keepdims=True)
    sd = np.where(sd == 0, np.nan, sd)
    z = (arr - mu) / sd
    return pd.Series(np.nanmean(z, axis=0), index=log.columns), found


def spearman(x, y) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    n = int(m.sum())
    if n < 4:
        return {"n": n, "rho": np.nan, "p": np.nan}
    r, p = stats.spearmanr(x[m], y[m])
    return {"n": n, "rho": float(r), "p": float(p)}


def quartiles(s: pd.Series) -> pd.Series:
    out = pd.Series(index=s.index, dtype=object)
    ok = s.dropna()
    if ok.empty:
        return out
    labels = pd.qcut(ok.rank(method="first"), 4, labels=["Q1", "Q2", "Q3", "Q4"])
    out.loc[ok.index] = labels.astype(str)
    return out


def mwu(a, b) -> dict:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    n4, n1 = len(a), len(b)
    rec = {
        "n_q4": n4,
        "n_q1": n1,
        "median_q4": float(np.median(a)) if n4 else np.nan,
        "median_q1": float(np.median(b)) if n1 else np.nan,
        "delta_median": np.nan,
        "U": np.nan,
        "rank_biserial_q4_gt_q1": np.nan,
        "cliffs_delta_q4_minus_q1": np.nan,
        "p": np.nan,
    }
    if n4 < 2 or n1 < 2:
        return rec
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    r = 2.0 * float(U) / (n4 * n1) - 1.0
    rec.update(
        {
            "median_q4": float(np.median(a)),
            "median_q1": float(np.median(b)),
            "delta_median": float(np.median(a) - np.median(b)),
            "U": float(U),
            "rank_biserial_q4_gt_q1": float(r),
            "cliffs_delta_q4_minus_q1": float(r),
            "p": float(p),
        }
    )
    return rec


def fmt_p(p: float) -> str:
    if p != p:
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_r(x: float, digits: int = 3) -> str:
    if x != x:
        return "NA"
    return f"{x:.{digits}f}"


def save_scatter(path: Path, x, y, title: str, xlab: str, ylab: str) -> None:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    fig, ax = plt.subplots(figsize=(4.0, 4.0))
    ax.scatter(x[m], y[m], s=28, c="#3d5a80", alpha=0.85)
    if m.sum() >= 3:
        lr = stats.linregress(x[m], y[m])
        xs = np.linspace(float(np.nanmin(x[m])), float(np.nanmax(x[m])), 50)
        ax.plot(xs, lr.intercept + lr.slope * xs, color="#e07a5f", lw=1.4)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def save_box(path: Path, q4, q1, title: str, ylab: str) -> None:
    fig, ax = plt.subplots(figsize=(3.6, 4.0))
    data = [np.asarray(q1, float), np.asarray(q4, float)]
    data = [d[~np.isnan(d)] for d in data]
    bp = ax.boxplot(
        data,
        tick_labels=[f"Q1\nn={len(data[0])}", f"Q4\nn={len(data[1])}"],
        widths=0.55,
        showfliers=False,
        patch_artist=True,
    )
    for patch, c in zip(bp["boxes"], ["#c9d6df", "#f6b26b"]):
        patch.set_facecolor(c)
        patch.set_alpha(0.85)
    rng = np.random.default_rng(0)
    for i, d in enumerate(data, start=1):
        ax.scatter(i + rng.uniform(-0.08, 0.08, size=len(d)), d, s=22, c="k", alpha=0.7, zorder=3)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    counts_path = fetch(COUNTS_URL, CACHE / "GSE182328_Gene_counts_matrix.txt.gz")
    matrix_path = fetch(MATRIX_URL, CACHE / "GSE182328_series_matrix.txt.gz")

    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    counts = counts.apply(pd.to_numeric, errors="coerce")
    counts.index = counts.index.astype(str)
    # one row per HUGO symbol in this author matrix
    dup = counts.index.duplicated()
    if dup.any():
        counts = counts.groupby(counts.index).mean()

    log = log2cpm(counts)
    meta = parse_series(matrix_path)
    meta["sample"] = meta["title"]

    present = {g: g in log.index for g in GENES}
    missing = [g for g, ok in present.items() if not ok]
    if missing:
        raise SystemExit(f"missing genes in author matrix: {missing}")

    scores = pd.DataFrame({g: log.loc[g].astype(float) for g in GENES})
    scores.index.name = "sample"
    scores = scores.reset_index()
    df = meta.merge(scores, on="sample", how="inner")
    df["patient"] = df["title"]
    assert df["patient"].is_unique
    assert len(df) == 44

    hla_i, hla_i_used = mean_z(log, HLA_I)
    cd8sig, cd8_used = mean_z(log, ["CD8A", "CD8B"])
    df = df.merge(
        pd.DataFrame({"sample": hla_i.index, "HLA_I": hla_i.values, "CD8sig": cd8sig.values}),
        on="sample",
        how="left",
    )
    df["CLDN4_q"] = quartiles(df["CLDN4"])
    df["CD8A_q"] = quartiles(df["CD8A"])

    # --- continuous extra cuts (CD274 / HLA). CD8A reprint is labeled known. ---
    continuous_rows = []
    known_rows = []
    for endpoint, col, kind in [
        ("CD8A", "CD8A", "known"),
        ("CD8B", "CD8B", "extra"),
        ("CD8sig (CD8A+CD8B mean-z)", "CD8sig", "extra"),
        ("CD274", "CD274", "extra"),
        ("HLA-A", "HLA-A", "extra"),
        ("HLA-B", "HLA-B", "extra"),
        ("HLA-C", "HLA-C", "extra"),
        ("HLA_I (HLA-A/B/C mean-z)", "HLA_I", "extra"),
        ("HLA-E", "HLA-E", "extra"),
        ("HLA-F", "HLA-F", "extra"),
        ("HLA-G", "HLA-G", "extra"),
        ("HLA-DRA", "HLA-DRA", "extra"),
        ("B2M", "B2M", "extra"),
        ("TACSTD2 (companion)", "TACSTD2", "companion"),
    ]:
        sp = spearman(df["CLDN4"], df[col])
        rec = {
            "cut": "continuous_Spearman",
            "predictor": "CLDN4",
            "endpoint": endpoint,
            "kind": kind,
            "n": sp["n"],
            "n_q4": "",
            "n_q1": "",
            "metric": "Spearman_rho",
            "effect": sp["rho"],
            "p": sp["p"],
            "median_q4": "",
            "median_q1": "",
            "delta_median": "",
            "U": "",
            "rank_biserial_q4_gt_q1": "",
            "note": "already known leftover row" if kind == "known" else "",
        }
        (known_rows if kind == "known" else continuous_rows).append(rec)

    # --- Q4 vs Q1 extra cuts ---
    q_rows = []
    # Primary requested: Q4 vs Q1 CD8. Two honest readings, both reported.
    # 1) CLDN4 quartiles, CD8A as endpoint (standard extra cut when continuous CD8A is known)
    # 2) CD8A quartiles, CLDN4 as endpoint ("Q4 vs Q1 CD8")
    for predictor_q, predictor, endpoint, kind, note in [
        ("CLDN4_q", "CLDN4", "CD8A", "extra_primary", "CLDN4 Q4 vs Q1; endpoint CD8A"),
        ("CLDN4_q", "CLDN4", "CD8sig", "extra", "CLDN4 Q4 vs Q1; endpoint CD8A+CD8B mean-z"),
        ("CLDN4_q", "CLDN4", "CD274", "extra", "CLDN4 Q4 vs Q1; endpoint CD274"),
        ("CLDN4_q", "CLDN4", "HLA_I", "extra", "CLDN4 Q4 vs Q1; endpoint HLA-A/B/C mean-z"),
        ("CLDN4_q", "CLDN4", "HLA-A", "extra", "CLDN4 Q4 vs Q1; endpoint HLA-A"),
        ("CLDN4_q", "CLDN4", "HLA-B", "extra", "CLDN4 Q4 vs Q1; endpoint HLA-B"),
        ("CLDN4_q", "CLDN4", "HLA-C", "extra", "CLDN4 Q4 vs Q1; endpoint HLA-C"),
        ("CD8A_q", "CD8A", "CLDN4", "extra_primary", "CD8A Q4 vs Q1; endpoint CLDN4"),
    ]:
        q4 = df.loc[df[predictor_q].eq("Q4"), endpoint]
        q1 = df.loc[df[predictor_q].eq("Q1"), endpoint]
        blk = mwu(q4, q1)
        q_rows.append(
            {
                "cut": f"{predictor}_Q4_vs_Q1",
                "predictor": predictor,
                "endpoint": endpoint,
                "kind": kind,
                "n": int(blk["n_q4"] + blk["n_q1"]),
                "n_q4": blk["n_q4"],
                "n_q1": blk["n_q1"],
                "metric": "MWU_rank_biserial_Q4_gt_Q1",
                "effect": blk["rank_biserial_q4_gt_q1"],
                "p": blk["p"],
                "median_q4": blk["median_q4"],
                "median_q1": blk["median_q1"],
                "delta_median": blk["delta_median"],
                "U": blk["U"],
                "rank_biserial_q4_gt_q1": blk["rank_biserial_q4_gt_q1"],
                "note": note,
            }
        )

    all_rows = known_rows + continuous_rows + q_rows
    tests = pd.DataFrame(all_rows)
    tests.to_csv(TAB / "extra_cuts.tsv", sep="\t", index=False)

    n_table = pd.DataFrame(
        [
            {"item": "GEO series samples (lung tumor RNA-seq)", "n": 44},
            {"item": "Author count-matrix columns", "n": int(counts.shape[1])},
            {"item": "Joined patient table (title = sample ID)", "n": int(len(df))},
            {"item": "Unique patients", "n": int(df["patient"].nunique())},
            {"item": "CLDN4 / CD8A / CD274 / HLA-A/B/C non-NA", "n": 44},
            {"item": "CLDN4 Q1 / Q2 / Q3 / Q4", "n": "11 / 11 / 11 / 11"},
            {"item": "CD8A Q1 / Q2 / Q3 / Q4", "n": "11 / 11 / 11 / 11"},
            {"item": "Q4 vs Q1 used", "n": "11 vs 11"},
            {"item": "RECIST / DCB / PFS / MPR on GEO", "n": 0},
            {"item": "HLA-I genes used", "n": ",".join(hla_i_used)},
            {"item": "CD8sig genes used", "n": ",".join(cd8_used)},
        ]
    )
    n_table.to_csv(TAB / "n_table.tsv", sep="\t", index=False)
    df.to_csv(TAB / "samples.tsv", sep="\t", index=False)

    # figures
    save_scatter(
        FIG / "cldn4_vs_cd274.png",
        df["CLDN4"],
        df["CD274"],
        "GSE182328 CLDN4 vs CD274",
        "CLDN4 log2(CPM+1)",
        "CD274 log2(CPM+1)",
    )
    save_scatter(
        FIG / "cldn4_vs_hla_i.png",
        df["CLDN4"],
        df["HLA_I"],
        "GSE182328 CLDN4 vs HLA-A/B/C mean-z",
        "CLDN4 log2(CPM+1)",
        "HLA-I mean-z",
    )
    save_box(
        FIG / "cldn4_q4q1_cd8a.png",
        df.loc[df["CLDN4_q"].eq("Q4"), "CD8A"],
        df.loc[df["CLDN4_q"].eq("Q1"), "CD8A"],
        "GSE182328 CD8A in CLDN4 Q4 vs Q1",
        "CD8A log2(CPM+1)",
    )
    save_box(
        FIG / "cd8a_q4q1_cldn4.png",
        df.loc[df["CD8A_q"].eq("Q4"), "CLDN4"],
        df.loc[df["CD8A_q"].eq("Q1"), "CLDN4"],
        "GSE182328 CLDN4 in CD8A Q4 vs Q1",
        "CLDN4 log2(CPM+1)",
    )

    def pick(endpoint: str, cut: str = "continuous_Spearman") -> dict:
        hit = tests[(tests["endpoint"] == endpoint) & (tests["cut"] == cut)]
        return hit.iloc[0].to_dict()

    def pick_q(cut: str, endpoint: str) -> dict:
        hit = tests[(tests["cut"] == cut) & (tests["endpoint"] == endpoint)]
        return hit.iloc[0].to_dict()

    summary = {
        "cohort": "GSE182328",
        "n_patients": 44,
        "already_known": {
            "CLDN4_vs_CD8A_Spearman": pick("CD8A"),
        },
        "extra_continuous": {
            "CLDN4_vs_CD274": pick("CD274"),
            "CLDN4_vs_HLA_A": pick("HLA-A"),
            "CLDN4_vs_HLA_B": pick("HLA-B"),
            "CLDN4_vs_HLA_C": pick("HLA-C"),
            "CLDN4_vs_HLA_I": pick("HLA_I (HLA-A/B/C mean-z)"),
        },
        "extra_quartile": {
            "CLDN4_Q4_vs_Q1_CD8A": pick_q("CLDN4_Q4_vs_Q1", "CD8A"),
            "CD8A_Q4_vs_Q1_CLDN4": pick_q("CD8A_Q4_vs_Q1", "CLDN4"),
        },
        "counts_url": COUNTS_URL,
        "hla_i_genes": hla_i_used,
        "cd8sig_genes": cd8_used,
        "quartile_rule": "pd.qcut(rank(method='first'), 4) among the 44 tumors",
        "normalization": "log2(CPM+1) from author gene-count matrix",
        "unit": "patient = GEO title / count-matrix column (1:1, n=44)",
        "no_ici_endpoint_on_geo": True,
    }
    (TAB / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    # FINDING.md is written by the caller from these tables so the PR table
    # cannot drift from the computed TSV. Still emit a machine-readable block.
    print(json.dumps(summary, indent=2, default=str))
    print("wrote", TAB / "extra_cuts.tsv")


if __name__ == "__main__":
    main()
