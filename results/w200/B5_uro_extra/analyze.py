#!/usr/bin/env python3
"""TACSTD2 / CLDN4 versus atezolizumab response in two public urothelial cohorts.

Primary question (prespecified):
    Do pretreatment tumor TACSTD2 or CLDN4 transcript levels separate
    RECIST responders (CR/PR) from nonresponders (SD/PD)?

Cohorts:
    1. IMvigor210 — Mariathasan et al., Nature 2018 (atezolizumab; n=348 RNA)
    2. Snyder 2017 — Snyder et al., PLoS Med 2017 (atezolizumab; n=25 RNA)

Expression scale: log2(CPM + 1) from raw/estimated counts and library size.
This is a within-gene, between-sample rank comparison. It is not a
genome-wide differential-expression analysis.

The script reads the small processed tables in processed/. Rebuild those
tables from the public sources with:  python3 analyze.py --rebuild-processed
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats
from sklearn.metrics import roc_auc_score
from statsmodels.stats.contingency_tables import Table2x2

HERE = Path(__file__).resolve().parent
PROCESSED = HERE / "processed"
FIGDIR = HERE / "figures"
TABLEDIR = HERE / "tables"
DOWNLOADS = HERE / "downloads"

IMVIGOR_URL = (
    "http://research-pub.gene.com/IMvigor210CoreBiologies/"
    "packageVersions/IMvigor210CoreBiologies_1.0.1.tar.gz"
)
SNYDER_KALLISTO_URL = (
    "https://raw.githubusercontent.com/hammerlab/"
    "multi-omic-urothelial-anti-pdl1/master/data_kallisto.csv"
)
SNYDER_CLINICAL_URL = (
    "https://raw.githubusercontent.com/hammerlab/"
    "multi-omic-urothelial-anti-pdl1/master/data_clinical.csv"
)

TARGET_GENES = [
    "TACSTD2",
    "CLDN4",
    "CD274",
    "CD8A",
    "GZMA",
    "PRF1",
    "CXCL9",
    "CXCL13",
    "EPCAM",
    "KRT5",
    "KRT20",
    "GATA3",
    "FOXA1",
    "UPK2",
    "CDH1",
    "TGFB1",
]
PRIMARY_GENES = ["TACSTD2", "CLDN4"]
CONTROL_GENES = ["CD8A", "CXCL9"]
RNG = np.random.default_rng(1)
N_BOOT = 2000

RESPONDER_COLOR = "#2c7bb6"
NONRESP_COLOR = "#b2182b"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def log2cpm(count, libsize) -> np.ndarray:
    count = np.asarray(count, dtype=float)
    libsize = np.asarray(libsize, dtype=float)
    return np.log2(count / libsize * 1e6 + 1.0)


def bootstrap_auc(y: np.ndarray, x: np.ndarray, n: int = N_BOOT) -> tuple[float, float, float]:
    y = np.asarray(y)
    x = np.asarray(x)
    auc = float(roc_auc_score(y, x))
    boots = []
    for _ in range(n):
        idx = RNG.integers(0, len(y), len(y))
        if y[idx].min() == y[idx].max():
            continue
        boots.append(float(roc_auc_score(y[idx], x[idx])))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return auc, float(lo), float(hi)


def mannwhitney(x_pos: np.ndarray, x_neg: np.ndarray) -> float:
    return float(stats.mannwhitneyu(x_pos, x_neg, alternative="two-sided").pvalue)


def median_or(expr: np.ndarray, resp: np.ndarray) -> dict:
    """Response odds in expression-high vs expression-low (cohort median)."""
    med = float(np.median(expr))
    high = expr >= med
    a = int((high & (resp == 1)).sum())
    b = int((high & (resp == 0)).sum())
    c = int((~high & (resp == 1)).sum())
    d = int((~high & (resp == 0)).sum())
    table = [[a, b], [c, d]]
    fisher_or, fisher_p = stats.fisher_exact(table)
    try:
        t = Table2x2(np.array(table, dtype=float))
        or_ci = t.oddsratio_confint()
        or_pt = float(t.oddsratio)
    except Exception:
        # Haldane-Anscombe if a cell is zero
        or_pt = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        se = np.sqrt(1 / (a + 0.5) + 1 / (b + 0.5) + 1 / (c + 0.5) + 1 / (d + 0.5))
        or_ci = (float(np.exp(np.log(or_pt) - 1.96 * se)), float(np.exp(np.log(or_pt) + 1.96 * se)))
    return {
        "median_cut": med,
        "high_R": a,
        "high_NR": b,
        "low_R": c,
        "low_NR": d,
        "OR": float(or_pt),
        "OR_fisher": float(fisher_or),
        "OR_lo": float(or_ci[0]),
        "OR_hi": float(or_ci[1]),
        "fisher_p": float(fisher_p),
    }


def gene_vs_response(expr: np.ndarray, resp: np.ndarray) -> dict:
    x_pos = expr[resp == 1]
    x_neg = expr[resp == 0]
    auc, auc_lo, auc_hi = bootstrap_auc(resp, expr)
    out = {
        "n_R": int(x_pos.size),
        "n_NR": int(x_neg.size),
        "median_R": float(np.median(x_pos)),
        "median_NR": float(np.median(x_neg)),
        "mean_R": float(np.mean(x_pos)),
        "mean_NR": float(np.mean(x_neg)),
        "delta_median": float(np.median(x_pos) - np.median(x_neg)),
        "auc_R_higher": auc,
        "auc_lo": auc_lo,
        "auc_hi": auc_hi,
        "mw_p": mannwhitney(x_pos, x_neg),
    }
    out.update(median_or(expr, resp))
    return out


def cox_os(expr: np.ndarray, time: np.ndarray, event: np.ndarray) -> dict:
    df = pd.DataFrame({"g": expr, "t": time, "e": event}).dropna()
    if df["e"].sum() < 5 or df["g"].nunique() < 3:
        return {"n": int(len(df)), "events": int(df["e"].sum()), "HR": np.nan, "HR_lo": np.nan, "HR_hi": np.nan, "p": np.nan}
    cph = CoxPHFitter()
    cph.fit(df, duration_col="t", event_col="e")
    s = cph.summary.loc["g"]
    rho, rho_p = stats.spearmanr(df["g"], df["t"])
    return {
        "n": int(len(df)),
        "events": int(df["e"].sum()),
        "HR": float(s["exp(coef)"]),
        "HR_lo": float(s["exp(coef) lower 95%"]),
        "HR_hi": float(s["exp(coef) upper 95%"]),
        "p": float(s["p"]),
        "spearman_os": float(rho),
        "spearman_os_p": float(rho_p),
    }


def download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    import urllib.request

    print(f"Downloading {url} -> {dest}", flush=True)
    urllib.request.urlretrieve(url, dest)
    return dest


def rebuild_processed() -> None:
    """Download public sources and write the small processed tables."""
    PROCESSED.mkdir(parents=True, exist_ok=True)
    DOWNLOADS.mkdir(parents=True, exist_ok=True)

    tar_path = download(IMVIGOR_URL, DOWNLOADS / "IMvigor210CoreBiologies_1.0.1.tar.gz")
    extract_dir = DOWNLOADS / "IMvigor210CoreBiologies"
    cds = extract_dir / "data" / "cds.RData"
    if not cds.exists():
        subprocess.check_call(["tar", "xzf", str(tar_path), "-C", str(DOWNLOADS)])
        # tarball root is IMvigor210CoreBiologies/
        cds = DOWNLOADS / "IMvigor210CoreBiologies" / "data" / "cds.RData"
    subprocess.check_call(
        ["Rscript", str(HERE / "extract_imvigor.R"), str(cds), str(PROCESSED)]
    )

    kal = download(SNYDER_KALLISTO_URL, DOWNLOADS / "snyder_data_kallisto.csv")
    clin = download(SNYDER_CLINICAL_URL, DOWNLOADS / "snyder_data_clinical.csv")

    genes = set(TARGET_GENES)
    lib: dict[str, float] = {}
    rows: list[tuple[str, str, float]] = []
    with kal.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            pid = str(row["patient_id"]).zfill(4)
            count = float(row["est_counts"])
            lib[pid] = lib.get(pid, 0.0) + count
            if row["gene_name"] in genes:
                rows.append((pid, row["gene_name"], count))
    with (PROCESSED / "snyder_libsize.csv").open("w", newline="") as handle:
        w = csv.writer(handle)
        w.writerow(["patient_id", "libsize"])
        for pid, value in sorted(lib.items()):
            w.writerow([pid, value])
    with (PROCESSED / "snyder_gene_counts.csv").open("w", newline="") as handle:
        w = csv.writer(handle)
        w.writerow(["patient_id", "gene", "est_counts"])
        w.writerows(rows)

    raw_clin = pd.read_csv(clin)
    keep_cols = [
        c
        for c in raw_clin.columns
        if not str(c).startswith("Unnamed")
    ]
    raw_clin[keep_cols].to_csv(PROCESSED / "snyder_clinical.csv", index=False)

    prov = [
        ("imvigor210_package", IMVIGOR_URL, sha256(tar_path), tar_path.stat().st_size),
        ("snyder_kallisto", SNYDER_KALLISTO_URL, sha256(kal), kal.stat().st_size),
        ("snyder_clinical", SNYDER_CLINICAL_URL, sha256(clin), clin.stat().st_size),
    ]
    with (PROCESSED / "provenance.tsv").open("w") as handle:
        handle.write("name\turl\tsha256\tbytes\n")
        for row in prov:
            handle.write("\t".join(map(str, row)) + "\n")
    print("Rebuilt processed/ tables.", flush=True)


def load_imvigor() -> pd.DataFrame:
    pdata = pd.read_csv(PROCESSED / "imvigor210_pdata.csv")
    genes = pd.read_csv(PROCESSED / "imvigor210_gene_counts.csv")
    df = pdata.merge(genes, on="sample_id", how="inner", validate="one_to_one")
    if len(df) != 348:
        raise SystemExit(f"IMvigor210 expected 348 RNA samples, got {len(df)}")
    recist = df["Best Confirmed Overall Response"].astype(str)
    df["recist"] = recist
    df["evaluable"] = recist.isin(["CR", "PR", "SD", "PD"])
    df["responder"] = recist.isin(["CR", "PR"])
    # package binaryResponse is CR/PR vs SD/PD; NE is NA
    pkg = df["binaryResponse"]
    derived = pd.Series(pd.NA, index=df.index, dtype="object")
    derived.loc[df["evaluable"] & df["responder"]] = "CR/PR"
    derived.loc[df["evaluable"] & ~df["responder"]] = "SD/PD"
    mismatch = df["evaluable"] & pkg.notna() & (pkg != derived)
    if mismatch.any():
        raise SystemExit("IMvigor210 binaryResponse does not match RECIST grouping")
    for gene in PRIMARY_GENES + CONTROL_GENES + ["GZMA", "PRF1", "CD274"]:
        df[f"{gene}_log2cpm"] = log2cpm(df[gene], df["libsize"])
    df["CYT_log2cpm"] = (df["GZMA_log2cpm"] + df["PRF1_log2cpm"]) / 2.0
    df["cohort"] = "IMvigor210"
    return df


def load_snyder() -> pd.DataFrame:
    clin = pd.read_csv(PROCESSED / "snyder_clinical.csv")
    clin["patient_id"] = clin["patient_id"].astype(str).str.zfill(4)
    genes = pd.read_csv(PROCESSED / "snyder_gene_counts.csv")
    genes["patient_id"] = genes["patient_id"].astype(str).str.zfill(4)
    lib = pd.read_csv(PROCESSED / "snyder_libsize.csv")
    lib["patient_id"] = lib["patient_id"].astype(str).str.zfill(4)
    wide = genes.pivot(index="patient_id", columns="gene", values="est_counts")
    df = clin.merge(wide, on="patient_id", how="inner").merge(lib, on="patient_id", how="inner")
    if len(df) != 25:
        raise SystemExit(f"Snyder expected 25 RNA patients, got {len(df)}")
    recist = df["Best Response RECIST 1.1"].astype(str).str.strip()
    recist = recist.replace(
        {
            "Only scanned baseline": "NE",
            "Only scanned at baseline": "NE",
        }
    )
    df["recist"] = recist
    df["evaluable"] = recist.isin(["CR", "PR", "SD", "PD"])
    df["responder"] = recist.isin(["CR", "PR"])
    df["benefit"] = df["is_benefit"].astype(bool)
    for gene in PRIMARY_GENES + CONTROL_GENES + ["GZMA", "PRF1", "CD274"]:
        df[f"{gene}_log2cpm"] = log2cpm(df[gene], df["libsize"])
    df["CYT_log2cpm"] = (df["GZMA_log2cpm"] + df["PRF1_log2cpm"]) / 2.0
    df["cohort"] = "Snyder2017"
    return df


def add_row(rows: list[dict], cohort: str, endpoint: str, gene: str, stats_d: dict) -> None:
    rec = {"cohort": cohort, "endpoint": endpoint, "gene": gene}
    rec.update(stats_d)
    rows.append(rec)


def run_tests(imv: pd.DataFrame, sny: pd.DataFrame) -> tuple[list[dict], dict]:
    rows: list[dict] = []
    genes_imv = PRIMARY_GENES + CONTROL_GENES + ["CYT"]

    ev = imv[imv["evaluable"]].copy()
    y = ev["responder"].astype(int).to_numpy()
    for gene in genes_imv:
        col = "CYT_log2cpm" if gene == "CYT" else f"{gene}_log2cpm"
        add_row(rows, "IMvigor210", "RECIST_CRPR_vs_SDPD", gene, gene_vs_response(ev[col].to_numpy(), y))

    pd_only = imv[imv["recist"].isin(["CR", "PR", "PD"])].copy()
    y2 = pd_only["responder"].astype(int).to_numpy()
    for gene in PRIMARY_GENES + ["CD8A"]:
        add_row(
            rows,
            "IMvigor210",
            "RECIST_CRPR_vs_PD",
            gene,
            gene_vs_response(pd_only[f"{gene}_log2cpm"].to_numpy(), y2),
        )

    # sizeFactor sensitivity (package DESeq size factors)
    for gene in PRIMARY_GENES:
        x = np.log2(ev[gene].to_numpy() / ev["sizeFactor"].to_numpy() + 1.0)
        add_row(rows, "IMvigor210", "RECIST_sizeFactor_log2", gene, gene_vs_response(x, y))

    for gene in PRIMARY_GENES + ["CD8A"]:
        os = cox_os(imv[f"{gene}_log2cpm"].to_numpy(), imv["os"].to_numpy(), imv["censOS"].to_numpy())
        add_row(rows, "IMvigor210", "OS_cox_per_log2cpm", gene, os)

    # Snyder primary: RECIST-evaluable
    sev = sny[sny["evaluable"]].copy()
    ys = sev["responder"].astype(int).to_numpy()
    for gene in genes_imv:
        col = "CYT_log2cpm" if gene == "CYT" else f"{gene}_log2cpm"
        add_row(rows, "Snyder2017", "RECIST_CRPR_vs_SDPD", gene, gene_vs_response(sev[col].to_numpy(), ys))

    yb = sny["benefit"].astype(int).to_numpy()
    for gene in PRIMARY_GENES + ["CD8A"]:
        add_row(rows, "Snyder2017", "paper_clinical_benefit", gene, gene_vs_response(sny[f"{gene}_log2cpm"].to_numpy(), yb))

    for gene in PRIMARY_GENES + ["CD8A"]:
        os = cox_os(
            sny[f"{gene}_log2cpm"].to_numpy(),
            sny["os"].to_numpy(),
            sny["is_deceased"].astype(int).to_numpy(),
        )
        add_row(rows, "Snyder2017", "OS_cox_per_log2cpm", gene, os)

    # Context: IMvigor210 subtype / immune phenotype (exploratory)
    for stratum_col in ("TCGA Subtype", "Immune phenotype"):
        for level, sub in ev.groupby(stratum_col):
            if pd.isna(level) or len(sub) < 20:
                continue
            yy = sub["responder"].astype(int).to_numpy()
            if yy.min() == yy.max():
                continue
            for gene in PRIMARY_GENES:
                add_row(
                    rows,
                    "IMvigor210",
                    f"exploratory_within_{stratum_col.replace(' ', '_')}_{level}",
                    gene,
                    gene_vs_response(sub[f"{gene}_log2cpm"].to_numpy(), yy),
                )

    headline = {
        "imvigor_n_rna": int(len(imv)),
        "imvigor_n_evaluable": int(ev["evaluable"].sum()) if "evaluable" in ev else int(len(ev)),
        "imvigor_n_R": int(y.sum()),
        "imvigor_n_NR": int((1 - y).sum()),
        "snyder_n_rna": int(len(sny)),
        "snyder_n_evaluable": int(len(sev)),
        "snyder_n_R": int(ys.sum()),
        "snyder_n_NR": int((1 - ys).sum()),
        "snyder_n_benefit": int(yb.sum()),
        "tacstd2_cldn4_spearman_imvigor": float(stats.spearmanr(imv["TACSTD2_log2cpm"], imv["CLDN4_log2cpm"]).statistic),
        "tacstd2_cldn4_spearman_snyder": float(stats.spearmanr(sny["TACSTD2_log2cpm"], sny["CLDN4_log2cpm"]).statistic),
    }
    return rows, headline


def write_sample_tables(imv: pd.DataFrame, sny: pd.DataFrame) -> None:
    TABLEDIR.mkdir(parents=True, exist_ok=True)
    imv_out = imv[
        [
            "sample_id",
            "ANONPT_ID",
            "recist",
            "binaryResponse",
            "evaluable",
            "responder",
            "os",
            "censOS",
            "IC Level",
            "TCGA Subtype",
            "Immune phenotype",
            "FMOne mutation burden per MB",
            "libsize",
            "TACSTD2_log2cpm",
            "CLDN4_log2cpm",
            "CD8A_log2cpm",
            "CXCL9_log2cpm",
            "CYT_log2cpm",
        ]
    ].copy()
    imv_out.to_csv(TABLEDIR / "sample_imvigor210.csv", index=False)

    sny_out = sny[
        [
            "patient_id",
            "recist",
            "evaluable",
            "responder",
            "benefit",
            "os",
            "is_deceased",
            "pfs",
            "PD-L1",
            "libsize",
            "TACSTD2_log2cpm",
            "CLDN4_log2cpm",
            "CD8A_log2cpm",
            "CXCL9_log2cpm",
            "CYT_log2cpm",
        ]
    ].copy()
    sny_out.to_csv(TABLEDIR / "sample_snyder2017.csv", index=False)


def plot_boxplots(imv: pd.DataFrame, sny: pd.DataFrame) -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 7.2))
    panels = [
        (axes[0, 0], imv[imv["evaluable"]], "TACSTD2_log2cpm", "IMvigor210 TACSTD2"),
        (axes[0, 1], imv[imv["evaluable"]], "CLDN4_log2cpm", "IMvigor210 CLDN4"),
        (axes[1, 0], sny[sny["evaluable"]], "TACSTD2_log2cpm", "Snyder 2017 TACSTD2"),
        (axes[1, 1], sny[sny["evaluable"]], "CLDN4_log2cpm", "Snyder 2017 CLDN4"),
    ]
    for ax, df, col, title in panels:
        r = df.loc[df["responder"], col].to_numpy()
        n = df.loc[~df["responder"], col].to_numpy()
        data = [n, r]
        bp = ax.boxplot(
            data,
            tick_labels=["SD/PD", "CR/PR"],
            widths=0.55,
            patch_artist=True,
            medianprops={"color": "black", "linewidth": 1.4},
            flierprops={"marker": "o", "markersize": 3, "alpha": 0.5},
        )
        bp["boxes"][0].set_facecolor("#f4cccc")
        bp["boxes"][1].set_facecolor("#cfe2f3")
        # jittered points
        rng = np.random.default_rng(2)
        for i, vals, color in ((1, n, NONRESP_COLOR), (2, r, RESPONDER_COLOR)):
            x = i + rng.uniform(-0.12, 0.12, size=len(vals))
            ax.scatter(x, vals, s=14 if len(vals) < 40 else 8, c=color, alpha=0.65, zorder=3, linewidths=0)
        p = mannwhitney(r, n)
        ax.set_title(f"{title}\nMW p={p:.3g}", fontsize=10)
        ax.set_ylabel("log2(CPM + 1)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlabel(f"n={len(n)} / {len(r)}")
    fig.suptitle("Pretreatment expression vs RECIST response (atezolizumab, urothelial)", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGDIR / "response_boxplots.png", dpi=160)
    fig.savefig(FIGDIR / "response_boxplots.svg")
    plt.close(fig)


def plot_km(imv: pd.DataFrame) -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    for ax, gene in zip(axes, PRIMARY_GENES):
        col = f"{gene}_log2cpm"
        med = float(imv[col].median())
        high = imv[col] >= med
        km = KaplanMeierFitter()
        km.fit(imv.loc[~high, "os"], imv.loc[~high, "censOS"], label=f"{gene}-low (n={(~high).sum()})")
        km.plot(ax=ax, ci_show=False, color="#636363", linewidth=1.8)
        km.fit(imv.loc[high, "os"], imv.loc[high, "censOS"], label=f"{gene}-high (n={high.sum()})")
        km.plot(ax=ax, ci_show=False, color="#e6550d", linewidth=1.8)
        res = logrank_test(
            imv.loc[high, "os"],
            imv.loc[~high, "os"],
            imv.loc[high, "censOS"],
            imv.loc[~high, "censOS"],
        )
        ax.set_title(f"IMvigor210 OS by median {gene}\nlog-rank p={res.p_value:.3g}")
        ax.set_xlabel("Overall survival (months)")
        ax.set_ylabel("Survival probability")
        ax.set_ylim(0, 1.02)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGDIR / "imvigor210_os_km.png", dpi=160)
    fig.savefig(FIGDIR / "imvigor210_os_km.svg")
    plt.close(fig)


def plot_or_forest(rows: list[dict]) -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    want = [
        r
        for r in rows
        if r["endpoint"] == "RECIST_CRPR_vs_SDPD" and r["gene"] in PRIMARY_GENES + CONTROL_GENES
    ]
    gene_order = PRIMARY_GENES + CONTROL_GENES
    want = sorted(want, key=lambda r: (r["cohort"], gene_order.index(r["gene"])))
    labels = [f"{r['cohort']} {r['gene']}" for r in want]
    ors = [r["OR"] for r in want]
    lo = [r["OR_lo"] for r in want]
    hi = [r["OR_hi"] for r in want]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    y = np.arange(len(want))[::-1]
    ax.axvline(1.0, color="0.4", linewidth=1, linestyle="--")
    ax.errorbar(
        ors,
        y,
        xerr=[np.array(ors) - np.array(lo), np.array(hi) - np.array(ors)],
        fmt="o",
        color="#222222",
        ecolor="#222222",
        capsize=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xscale("log")
    ax.set_xlabel("Odds ratio for CR/PR (expression-high vs low, median split)")
    ax.set_xlim(0.15, 8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("Median-split response OR (OR<1 = high expression, worse response)")
    fig.tight_layout()
    fig.savefig(FIGDIR / "median_split_or.png", dpi=160)
    fig.savefig(FIGDIR / "median_split_or.svg")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rebuild-processed",
        action="store_true",
        help="Re-download public sources and rebuild processed/ tables",
    )
    args = parser.parse_args()
    if args.rebuild_processed or not (PROCESSED / "imvigor210_pdata.csv").exists():
        rebuild_processed()

    imv = load_imvigor()
    sny = load_snyder()
    rows, headline = run_tests(imv, sny)

    TABLEDIR.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(TABLEDIR / "tests.csv", index=False)
    write_sample_tables(imv, sny)

    # subtype context table
    ev = imv[imv["evaluable"]]
    subtype_rows = []
    for level, sub in ev.groupby("TCGA Subtype"):
        subtype_rows.append(
            {
                "TCGA_subtype": level,
                "n_evaluable": int(len(sub)),
                "n_CRPR": int(sub["responder"].sum()),
                "ORR": float(sub["responder"].mean()),
                "TACSTD2_median_log2cpm": float(sub["TACSTD2_log2cpm"].median()),
                "CLDN4_median_log2cpm": float(sub["CLDN4_log2cpm"].median()),
            }
        )
    pd.DataFrame(subtype_rows).to_csv(TABLEDIR / "imvigor210_subtype_context.csv", index=False)

    plot_boxplots(imv, sny)
    plot_km(imv)
    plot_or_forest(rows)

    primary = [
        r
        for r in rows
        if r["endpoint"] == "RECIST_CRPR_vs_SDPD" and r["gene"] in PRIMARY_GENES
    ]
    summary = {"headline": headline, "primary": primary}
    (TABLEDIR / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print("=== PRIMARY (RECIST CR/PR vs SD/PD, log2CPM) ===")
    for r in primary:
        print(
            f"{r['cohort']:12} {r['gene']:8}  "
            f"n={r['n_R']}/{r['n_NR']}  "
            f"Δmed={r['delta_median']:+.3f}  "
            f"AUC={r['auc_R_higher']:.3f} [{r['auc_lo']:.3f},{r['auc_hi']:.3f}]  "
            f"MW p={r['mw_p']:.3g}  "
            f"OR_high={r['OR']:.2f} [{r['OR_lo']:.2f},{r['OR_hi']:.2f}]  "
            f"Fisher p={r['fisher_p']:.3g}"
        )
    print("Wrote tables/ and figures/ under", HERE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
