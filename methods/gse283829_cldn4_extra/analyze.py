#!/usr/bin/env python3
"""GSE283829 leftover extra cuts: CLDN4 vs CD274/HLA and Q4 vs ImmuneScore.

Continuous CLDN4 vs ESTIMATE ImmuneScore is already reported on the leftover
OPEN lung ICI bulk page (ρ=−0.57, p=0.0020, n=27). This script does not
re-claim that test. It adds:

1. CLDN4 vs CD274 and HLA class I/II genes (Spearman + partial | ImmuneScore)
2. CLDN4 Q4 vs Q1 vs ImmuneScore (and the same extra genes)

Public GEO processed counts only. Patient is the unit. Honest n.
"""
from __future__ import annotations

import gzip
import json
import math
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
CACHE = Path("/tmp/gse283829_cldn4_extra")
FIG = HERE / "figures"
TAB = HERE / "tables"
PROC = HERE / "processed"
for d in (CACHE, FIG, TAB, PROC):
    d.mkdir(parents=True, exist_ok=True)

COUNTS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/GSE283829/suppl/"
    "GSE283829_raw_express_matrix_all_samples.txt.gz"
)
SERIES_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/GSE283829/matrix/"
    "GSE283829_series_matrix.txt.gz"
)

EST = pd.read_csv(HERE / "estimate_gene_sets.csv")
STROMAL = [g for g in EST["stromal_signature"].dropna().astype(str) if g]
IMMUNE = [g for g in EST["immune_signature"].dropna().astype(str) if g]
SYM2ENS = (
    pd.read_csv(HERE / "hgnc_symbol_ensembl.tsv", sep="\t")
    .dropna()
    .drop_duplicates("symbol")
    .set_index("symbol")["ensembl_gene_id"]
    .to_dict()
)
ENS2SYM = {v: k for k, v in SYM2ENS.items()}

# Extra genes (not scored as individual Spearman rows on the leftover page)
EXTRA_GENES = ["CD274", "HLA-A", "HLA-B", "HLA-C", "HLA-DRA"]
MHC_I_GENES = ["HLA-A", "HLA-B", "HLA-C"]
SEED = 20260817
N_BOOT = 2000
RNG = np.random.default_rng(SEED)

# Already-known leftover continuous ImmuneScore (do not re-claim as new)
KNOWN_IMMUNE = {"n": 27, "rho": -0.5671550671550671, "p": 0.0020357432118934676}


def fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    return dest


def _split_fields(line: str) -> list[str]:
    return [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]


def parse_series(path: Path) -> pd.DataFrame:
    stored: dict[str, list[str]] = {}
    char_rows = []
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_geo_accession"):
                stored["geo"] = _split_fields(line)
            elif line.startswith("!Sample_title"):
                stored["title"] = _split_fields(line)
            elif line.startswith("!Sample_source_name"):
                stored["source"] = _split_fields(line)
            elif line.startswith("!Sample_characteristics"):
                char_rows.append(_split_fields(line))
    geo = stored.get("geo") or []
    meta = pd.DataFrame({"gsm": geo})
    n = len(geo)
    for key in ("title", "source"):
        vals = stored.get(key) or [None] * n
        meta[key] = (vals + [None] * n)[:n]
    chars: dict[str, dict[str, str]] = defaultdict(dict)
    for vals in char_rows:
        for g, v in zip(geo, vals):
            if not v:
                continue
            if ":" in v:
                k, val = v.split(":", 1)
                chars[g][k.strip().lower()] = val.strip()
    keys = sorted({k for d in chars.values() for k in d})
    for k in keys:
        meta[k] = [chars[g].get(k) for g in geo]
    return meta


def map_index_to_symbol(idx: pd.Index) -> pd.Index:
    out = []
    for x in idx.astype(str):
        x0 = x.split(".")[0]
        if x0.startswith("ENSG") and x0 in ENS2SYM:
            out.append(ENS2SYM[x0])
        else:
            out.append(x0)
    return pd.Index(out)


def counts_to_log2cpm(expr: pd.DataFrame) -> pd.DataFrame:
    lib = expr.sum(axis=0).replace(0, np.nan)
    cpm = expr.div(lib, axis=1) * 1e6
    return np.log2(cpm + 1)


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    arr = df.to_numpy(dtype=float)
    mu = np.nanmean(arr, axis=1, keepdims=True)
    sd = np.nanstd(arr, axis=1, keepdims=True)
    sd = np.where(sd == 0, np.nan, sd)
    return pd.DataFrame((arr - mu) / sd, index=df.index, columns=df.columns)


def mean_z(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series | None, list[str]]:
    found = [g for g in genes if g in expr.index]
    if len(found) < 2:
        return None, found
    return zscore_rows(expr.loc[found]).mean(axis=0), found


def spearman(x, y) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    if m.sum() < 4:
        return {"n": int(m.sum()), "rho": np.nan, "p": np.nan}
    r, p = stats.spearmanr(x[m], y[m])
    return {"n": int(m.sum()), "rho": float(r), "p": float(p)}


def partial_spearman(x, y, z) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    z = np.asarray(z, float)
    m = ~(np.isnan(x) | np.isnan(y) | np.isnan(z))
    n = int(m.sum())
    if n < 5:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rx, ry, rz = stats.rankdata(x[m]), stats.rankdata(y[m]), stats.rankdata(z[m])
    Z = np.column_stack([np.ones(n), rz])
    bx, *_ = np.linalg.lstsq(Z, rx, rcond=None)
    by, *_ = np.linalg.lstsq(Z, ry, rcond=None)
    ex, ey = rx - Z @ bx, ry - Z @ by
    if np.std(ex) == 0 or np.std(ey) == 0:
        return {"n": n, "rho": np.nan, "p": np.nan}
    r = float(np.corrcoef(ex, ey)[0, 1])
    df = n - 3
    t = r * math.sqrt(df / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), df))
    return {"n": n, "rho": r, "p": p}


def extreme_quartile(a: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Equal-count Q1 / Q4 via rank qcut (ties broken by first appearance)."""
    valid = a.dropna()
    ranks = valid.rank(method="first")
    q = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"])
    q1 = pd.Series(False, index=a.index)
    q4 = pd.Series(False, index=a.index)
    q1.loc[q.index[q == "Q1"]] = True
    q4.loc[q.index[q == "Q4"]] = True
    return q1, q4


def rank_biserial(x_q4: np.ndarray, x_q1: np.ndarray) -> float:
    n4, n1 = len(x_q4), len(x_q1)
    if n4 < 2 or n1 < 2:
        return float("nan")
    u = stats.mannwhitneyu(x_q4, x_q1, alternative="two-sided").statistic
    return float(2 * u / (n4 * n1) - 1)


def bootstrap_rb(x_q4: np.ndarray, x_q1: np.ndarray, n_boot: int = N_BOOT) -> tuple[float, float]:
    if len(x_q4) < 3 or len(x_q1) < 3:
        return float("nan"), float("nan")
    vals = []
    for _ in range(n_boot):
        a = RNG.choice(x_q4, size=len(x_q4), replace=True)
        b = RNG.choice(x_q1, size=len(x_q1), replace=True)
        vals.append(rank_biserial(a, b))
    lo, hi = np.nanpercentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


def mwu_q4_q1(y: pd.Series, q1: pd.Series, q4: pd.Series) -> dict:
    y1 = pd.to_numeric(y[q1], errors="coerce").to_numpy(float)
    y4 = pd.to_numeric(y[q4], errors="coerce").to_numpy(float)
    y1, y4 = y1[~np.isnan(y1)], y4[~np.isnan(y4)]
    n1, n4 = len(y1), len(y4)
    if n1 < 2 or n4 < 2:
        return {
            "n_q1": n1,
            "n_q4": n4,
            "median_q1": np.nan,
            "median_q4": np.nan,
            "delta_median": np.nan,
            "U": np.nan,
            "rank_biserial": np.nan,
            "rb_lo": np.nan,
            "rb_hi": np.nan,
            "p": np.nan,
        }
    U, p = stats.mannwhitneyu(y4, y1, alternative="two-sided")
    rb = rank_biserial(y4, y1)
    lo, hi = bootstrap_rb(y4, y1)
    return {
        "n_q1": n1,
        "n_q4": n4,
        "median_q1": float(np.median(y1)),
        "median_q4": float(np.median(y4)),
        "delta_median": float(np.median(y4) - np.median(y1)),
        "U": float(U),
        "rank_biserial": rb,
        "rb_lo": lo,
        "rb_hi": hi,
        "p": float(p),
    }


def fmt_p(p: float) -> str:
    if p != p:
        return "NA"
    if p < 1e-3:
        return f"{p:.1e}"
    return f"{p:.3f}"


def fmt_num(x: float, digits: int = 3) -> str:
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
        xs = np.linspace(np.nanmin(x[m]), np.nanmax(x[m]), 50)
        ax.plot(xs, lr.intercept + lr.slope * xs, color="#e07a5f", lw=1.4)
        rho, p = stats.spearmanr(x[m], y[m])
        ax.text(0.04, 0.96, f"ρ={rho:.2f}\np={p:.3g}\nn={int(m.sum())}", transform=ax.transAxes, va="top", fontsize=8)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def save_q4_box(path: Path, y_q1, y_q4, title: str, ylab: str) -> None:
    fig, ax = plt.subplots(figsize=(3.6, 4.0))
    data = [np.asarray(y_q1, float), np.asarray(y_q4, float)]
    data = [d[~np.isnan(d)] for d in data]
    bp = ax.boxplot(data, tick_labels=[f"Q1\nn={len(data[0])}", f"Q4\nn={len(data[1])}"], widths=0.55, showfliers=False, patch_artist=True)
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
    counts_path = fetch(COUNTS_URL, CACHE / "GSE283829_raw_express_matrix_all_samples.txt.gz")
    series_path = fetch(SERIES_URL, CACHE / "GSE283829_series_matrix.txt.gz")

    meta = parse_series(series_path)
    expr = pd.read_csv(counts_path, sep="\t", index_col=0)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr.index = map_index_to_symbol(expr.index)
    expr = expr.groupby(expr.index).mean()
    log = counts_to_log2cpm(expr)

    immune, immune_used = mean_z(log, IMMUNE)
    stromal, stromal_used = mean_z(log, STROMAL)
    mhc_i, mhc_used = mean_z(log, MHC_I_GENES)

    meta["sample"] = "S" + meta["title"].astype(str)
    rows = []
    missing_genes = []
    for g in ["CLDN4", "TACSTD2", "CD8A"] + EXTRA_GENES:
        if g not in log.index:
            missing_genes.append(g)
    if missing_genes:
        raise SystemExit(f"Missing genes in mapped matrix: {missing_genes}")

    for _, rec in meta.iterrows():
        s = rec["sample"]
        if s not in log.columns:
            raise SystemExit(f"Sample {s} not in count matrix columns")
        row = rec.to_dict()
        row["CLDN4"] = float(log.loc["CLDN4", s])
        row["TACSTD2"] = float(log.loc["TACSTD2", s])
        row["CD8A"] = float(log.loc["CD8A", s])
        for g in EXTRA_GENES:
            row[g] = float(log.loc[g, s])
        row["MHC_I"] = float(mhc_i[s]) if mhc_i is not None else np.nan
        row["ImmuneScore"] = float(immune[s]) if immune is not None else np.nan
        row["StromalScore"] = float(stromal[s]) if stromal is not None else np.nan
        rows.append(row)
    df = pd.DataFrame(rows)
    if len(df) != 27:
        raise SystemExit(f"Expected 27 samples, got {len(df)}")
    if not df["gsm"].is_unique:
        raise SystemExit("Duplicate GSM")

    q1, q4 = extreme_quartile(df["CLDN4"])
    df["CLDN4_quartile"] = "Q2Q3"
    df.loc[q1, "CLDN4_quartile"] = "Q1"
    df.loc[q4, "CLDN4_quartile"] = "Q4"

    # Reproduce leftover continuous ImmuneScore (sanity; not a new claim)
    known_check = spearman(df["CLDN4"], df["ImmuneScore"])
    if abs(known_check["rho"] - KNOWN_IMMUNE["rho"]) > 0.02:
        raise SystemExit(
            f"ImmuneScore Spearman drifted from leftover known value: {known_check} vs {KNOWN_IMMUNE}"
        )

    tests = []
    # Extra continuous: CD274 / HLA (new). ImmuneScore continuous is known — record as given.
    tests.append(
        {
            "family": "already_known",
            "test": "CLDN4_vs_ImmuneScore_continuous",
            "endpoint": "ImmuneScore",
            "metric": "Spearman_rho",
            "n": KNOWN_IMMUNE["n"],
            "n_q1": "",
            "n_q4": "",
            "effect": KNOWN_IMMUNE["rho"],
            "p": KNOWN_IMMUNE["p"],
            "note": "already on leftover OPEN lung ICI bulk page; not a new claim",
        }
    )
    extra_endpoints = EXTRA_GENES + ["MHC_I"]
    for ep in extra_endpoints:
        sp = spearman(df["CLDN4"], df[ep])
        tests.append(
            {
                "family": "extra_continuous",
                "test": f"CLDN4_vs_{ep}",
                "endpoint": ep,
                "metric": "Spearman_rho",
                "n": sp["n"],
                "n_q1": "",
                "n_q4": "",
                "effect": sp["rho"],
                "p": sp["p"],
                "note": "log2(CPM+1); extra cut (not on leftover continuous ImmuneScore row)",
            }
        )
        psp = partial_spearman(df["CLDN4"], df[ep], df["ImmuneScore"])
        tests.append(
            {
                "family": "extra_partial",
                "test": f"CLDN4_vs_{ep}_partial_ImmuneScore",
                "endpoint": ep,
                "metric": "partial_Spearman_rho",
                "n": psp["n"],
                "n_q1": "",
                "n_q4": "",
                "effect": psp["rho"],
                "p": psp["p"],
                "note": "partial Spearman after ESTIMATE ImmuneScore mean-z",
            }
        )

    q_endpoints = ["ImmuneScore"] + extra_endpoints
    for ep in q_endpoints:
        blk = mwu_q4_q1(df[ep], q1, q4)
        tests.append(
            {
                "family": "extra_q4q1",
                "test": f"CLDN4_Q4_vs_Q1_{ep}",
                "endpoint": ep,
                "metric": "MWU_Q4_minus_Q1",
                "n": blk["n_q1"] + blk["n_q4"],
                "n_q1": blk["n_q1"],
                "n_q4": blk["n_q4"],
                "median_q1": blk["median_q1"],
                "median_q4": blk["median_q4"],
                "delta_median": blk["delta_median"],
                "U": blk["U"],
                "rank_biserial": blk["rank_biserial"],
                "rb_lo": blk["rb_lo"],
                "rb_hi": blk["rb_hi"],
                "effect": blk["rank_biserial"],
                "p": blk["p"],
                "note": "CLDN4 rank-qcut Q4 vs Q1; middle quartiles excluded",
            }
        )

    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(TAB / "extra_cuts.tsv", sep="\t", index=False)
    df.to_csv(PROC / "sample_scores.tsv", sep="\t", index=False)

    ranks = df["CLDN4"].rank(method="first")
    q_all = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"])
    q_counts = q_all.value_counts().sort_index()
    n_table = pd.DataFrame(
        [
            {"item": "GEO RNA-seq samples (matrix n)", "n": 27, "note": "do not write a larger n"},
            {"item": "CR / SD / PD", "n": "7 / 10 / 10", "note": "GEO disease stage; no PR"},
            {"item": "continuous extra tests (CLDN4 vs CD274/HLA)", "n": 27, "note": "pairwise-complete; all genes present"},
            {"item": "ESTIMATE ImmuneScore genes used", "n": len(immune_used), "note": f"{len(immune_used)}/{len(IMMUNE)} Yoshihara immune genes mapped"},
            {"item": "MHC-I mean-z genes", "n": f"{len(mhc_used)}/3", "note": ",".join(mhc_used)},
            {
                "item": "CLDN4 Q1 / Q2 / Q3 / Q4",
                "n": " / ".join(str(int(q_counts[lab])) for lab in ["Q1", "Q2", "Q3", "Q4"]),
                "note": "rank(method=first) then qcut; ties do not inflate a bin",
            },
            {
                "item": "Q4 vs Q1 used",
                "n": f"{int(q_counts['Q4'])} vs {int(q_counts['Q1'])}",
                "note": "middle 13 samples excluded; do not write n=27 for Q4 vs Q1",
            },
        ]
    )
    n_table.to_csv(TAB / "honest_n.tsv", sep="\t", index=False)

    save_scatter(
        FIG / "CLDN4_vs_CD274.png",
        df["CLDN4"],
        df["CD274"],
        "GSE283829 CLDN4 vs CD274 (extra)",
        "CLDN4 log2(CPM+1)",
        "CD274 log2(CPM+1)",
    )
    save_scatter(
        FIG / "CLDN4_vs_HLA_A.png",
        df["CLDN4"],
        df["HLA-A"],
        "GSE283829 CLDN4 vs HLA-A (extra)",
        "CLDN4 log2(CPM+1)",
        "HLA-A log2(CPM+1)",
    )
    save_scatter(
        FIG / "CLDN4_vs_MHC_I.png",
        df["CLDN4"],
        df["MHC_I"],
        "GSE283829 CLDN4 vs MHC-I mean-z (HLA-A/B/C)",
        "CLDN4 log2(CPM+1)",
        "MHC-I mean-z",
    )
    save_q4_box(
        FIG / "CLDN4_Q4_vs_Q1_ImmuneScore.png",
        df.loc[q1, "ImmuneScore"],
        df.loc[q4, "ImmuneScore"],
        "GSE283829 CLDN4 Q4 vs Q1 ImmuneScore",
        "ESTIMATE ImmuneScore (mean-z)",
    )

    summary = {
        "cohort": "GSE283829",
        "n_geo": 27,
        "n_q1": int(q_counts["Q1"]),
        "n_q4": int(q_counts["Q4"]),
        "immune_genes_used": len(immune_used),
        "stromal_genes_used": len(stromal_used),
        "mhc_i_genes": mhc_used,
        "known_immune_spearman": KNOWN_IMMUNE,
        "reproduced_immune_spearman": known_check,
        "seed": SEED,
        "bootstrap": N_BOOT,
    }
    (TAB / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(tests_df.to_string(index=False))
    print(n_table.to_string(index=False))


if __name__ == "__main__":
    main()
