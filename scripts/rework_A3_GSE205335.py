#!/usr/bin/env python3
"""REWORK A3 — GSE205335 malignant TACSTD2 after MPR check, 0-cell drop, %positive cutoffs.

Self-contained. Public GEO processed files only.
Writes results/rework/A3_GSE205335/.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "GSE205335"
OUT = ROOT / "results" / "rework" / "A3_GSE205335"
FIG = OUT / "figures"

GEO_CELL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
    "GSE205335_Lung_IO_CellIdentity.txt.gz"
)
GEO_RDS = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
    "GSE205335_Lung_IO_UMI_matrix.rds.gz"
)
GEO_SERIES = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335"

# Locked before inspecting cutoff p-values.
CELL_CUTOFFS = [
    {"id": "umi_ge_1", "label": "UMI ≥ 1 (detection)", "kind": "umi", "thr": 1},
    {"id": "umi_ge_2", "label": "UMI ≥ 2", "kind": "umi", "thr": 2},
    {"id": "umi_ge_3", "label": "UMI ≥ 3", "kind": "umi", "thr": 3},
    {"id": "log1p_ge_0.5", "label": "log1p(CP10K) ≥ 0.5", "kind": "log1p", "thr": 0.5},
    {"id": "log1p_ge_1.0", "label": "log1p(CP10K) ≥ 1.0", "kind": "log1p", "thr": 1.0},
]
PATIENT_PCT_CUTOFFS = [10.0, 25.0, 50.0, 75.0]
PRIMARY_CELL_CUTOFF = "umi_ge_1"
RESP_MAP = {"PR": "R", "CR": "R", "SD": "NR", "PD": "NR", "NE": "NE"}
NSCLC = {"ADC", "SQ"}


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    print(f"Downloading {url} -> {dest}")
    urllib.request.urlretrieve(url, dest)


def fetch_gsm_metadata() -> pd.DataFrame:
    cache = DATA / "gsm_sample_metadata.csv"
    if cache.exists():
        return pd.read_csv(cache)
    import re
    import time

    rows = []
    for n in range(6210624, 6210657):
        gsm = f"GSM{n}"
        url = (
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
            f"?acc={gsm}&targ=self&form=text&view=quick"
        )
        with urllib.request.urlopen(url, timeout=60) as r:
            txt = r.read().decode("utf-8", "replace")
        title = re.search(r"!Sample_title = (.+)", txt)
        chars = dict(
            p.split(": ", 1)
            for p in re.findall(r"!Sample_characteristics_ch1 = (.+)", txt)
            if ": " in p
        )
        rows.append(
            {
                "gsm": gsm,
                "title": title.group(1).strip() if title else "",
                "patient": chars.get("patient", ""),
                "tissue": chars.get("tissue", ""),
                "tumor_stage": chars.get("tumor stage", ""),
                "cancer_subtype": chars.get("cancer subtype", ""),
                "recist": chars.get("recist", ""),
                "platform": chars.get("platform", ""),
            }
        )
        time.sleep(0.1)
    gsm = pd.DataFrame(rows)
    gsm.to_csv(cache, index=False)
    return gsm


def ensure_extract() -> Path:
    percell = DATA / "percell_tacstd2.csv.gz"
    if percell.exists() and percell.stat().st_size > 0:
        return percell
    rds_gz = DATA / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    rds = DATA / "GSE205335_Lung_IO_UMI_matrix.rds"
    download(GEO_RDS, rds_gz)
    if not rds.exists():
        # GEO file is double-gzipped; gunzip once, leave inner gzip for readRDS.
        subprocess.check_call(["gzip", "-dc", str(rds_gz)], stdout=open(rds, "wb"))
    script = ROOT / "scripts" / "extract_gse205335_tacstd2.R"
    subprocess.check_call(["Rscript", str(script), str(rds), str(percell)])
    return percell


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """Cliff's delta: P(x>y) - P(x<y). Positive => x stochastically larger."""
    if len(x) == 0 or len(y) == 0:
        return float("nan")
    u, _ = stats.mannwhitneyu(x, y, alternative="two-sided")
    return float((2.0 * u) / (len(x) * len(y)) - 1.0)


def mwu(x: np.ndarray, y: np.ndarray) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    out = {
        "n_R": int(len(x)),
        "n_NR": int(len(y)),
        "median_R": float(np.median(x)) if len(x) else float("nan"),
        "median_NR": float(np.median(y)) if len(y) else float("nan"),
        "U": float("nan"),
        "p": float("nan"),
        "cliffs_delta_R_minus_NR": float("nan"),
    }
    if len(x) >= 2 and len(y) >= 2:
        u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
        out["U"] = float(u)
        out["p"] = float(p)
        out["cliffs_delta_R_minus_NR"] = cliffs_delta(x, y)
    return out


def fisher_high(r_high, r_low, nr_high, nr_low) -> dict:
    table = np.array([[r_high, r_low], [nr_high, nr_low]], dtype=int)
    if table.min() < 0 or table.sum() == 0:
        return {"odds_ratio": float("nan"), "p": float("nan")}
    # If a row/col is all zero, fisher is undefined.
    if (table.sum(axis=0) == 0).any() or (table.sum(axis=1) == 0).any():
        return {"odds_ratio": float("nan"), "p": float("nan"), "note": "degenerate_table"}
    or_, p = stats.fisher_exact(table, alternative="two-sided")
    return {"odds_ratio": float(or_), "p": float(p)}


def fmt_p(p: float) -> str:
    if p != p:
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_num(x: float, nd=3) -> str:
    if x != x:
        return "NA"
    return f"{x:.{nd}f}"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    download(GEO_CELL, DATA / "GSE205335_Lung_IO_CellIdentity.txt.gz")
    percell_path = ensure_extract()
    gsm = fetch_gsm_metadata()
    # normalize column names if cache came from the hunt file
    rename = {
        "tumor stage": "tumor_stage",
        "cancer subtype": "cancer_subtype",
    }
    gsm = gsm.rename(columns=rename)
    if "cancer_subtype" not in gsm.columns and "cancer subtype" in gsm.columns:
        gsm["cancer_subtype"] = gsm["cancer subtype"]

    ident = pd.read_csv(DATA / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    percell = pd.read_csv(percell_path)
    df = ident.merge(percell, on="barcode", how="inner", validate="1:1")
    if len(df) != len(ident) or len(df) != len(percell):
        raise SystemExit(f"barcode merge mismatch: ident={len(ident)} percell={len(percell)} merged={len(df)}")

    gsm["token"] = gsm["title"].str.split().str[1].str.replace("_", "-", regex=False)
    gsm["suffix"] = np.where(gsm["platform"].str.contains("3", regex=False), "3P", "5P")
    # hunt used 3' -> 3P; platform strings are "Single Cell 3'" / "Single Cell 5'"
    gsm["orig.ident"] = gsm["token"] + "-" + gsm["suffix"]
    unmapped = set(df["orig.ident"]) - set(gsm["orig.ident"])
    if unmapped:
        raise SystemExit(f"unmapped orig.ident: {sorted(unmapped)}")

    df = df.merge(
        gsm[
            [
                "orig.ident",
                "gsm",
                "patient",
                "tissue",
                "recist",
                "platform",
                "cancer_subtype",
                "tumor_stage",
            ]
        ],
        on="orig.ident",
        how="left",
    )
    df["response"] = df["recist"].map(RESP_MAP)
    df["tacstd2_cp10k"] = df["tacstd2_umi"] / df["total_umi"] * 1e4
    df["tacstd2_log1p_cp10k"] = np.log1p(df["tacstd2_cp10k"])
    df["is_malignant"] = df["lineage.sub"] == "Malignant cells"
    df["is_tnk"] = df["lineage.total"] == "T/NK cells"
    df["is_core_cell"] = df["core.patient"] == "Core"
    df["is_nsclc"] = df["cancer_subtype"].isin(NSCLC)

    # Label scan: MPR / NMPR / pathologic response
    mpr_hits = []
    for col in ["recist", "tissue", "cancer_subtype", "tumor_stage", "title"]:
        if col in gsm.columns:
            s = gsm[col].astype(str).str.upper()
            if s.str.contains("MPR|NMPR|PATHOLOG", regex=True).any():
                mpr_hits.append(col)
    mpr_labeled = len(mpr_hits) > 0

    # Patient-level coverage (all 26 patients)
    all_pat = (
        gsm.groupby("patient", as_index=False)
        .agg(
            recist=("recist", "first"),
            cancer_subtype=("cancer_subtype", "first"),
            tumor_stage=("tumor_stage", "first"),
            n_samples=("gsm", "nunique"),
            tissues=("tissue", lambda s: ";".join(sorted(set(s)))),
        )
    )
    all_pat["response"] = all_pat["recist"].map(RESP_MAP)
    mal_n = df.loc[df["is_malignant"]].groupby("patient").size().rename("n_malignant")
    tnk_n = df.loc[df["is_tnk"]].groupby("patient").size().rename("n_tnk")
    core_mal_n = (
        df.loc[df["is_malignant"] & df["is_core_cell"]]
        .groupby("patient")
        .size()
        .rename("n_malignant_core")
    )
    all_pat = all_pat.merge(mal_n, on="patient", how="left")
    all_pat = all_pat.merge(tnk_n, on="patient", how="left")
    all_pat = all_pat.merge(core_mal_n, on="patient", how="left")
    for c in ["n_malignant", "n_tnk", "n_malignant_core"]:
        all_pat[c] = all_pat[c].fillna(0).astype(int)
    all_pat["zero_malignant"] = all_pat["n_malignant"] == 0
    all_pat["in_primary"] = (~all_pat["zero_malignant"]) & all_pat["response"].isin(["R", "NR"])
    all_pat["is_nsclc"] = all_pat["cancer_subtype"].isin(NSCLC)
    all_pat["has_core_malignant"] = all_pat["n_malignant_core"] > 0

    # Per-patient malignant metrics at each cell cutoff
    mal = df.loc[df["is_malignant"]].copy()
    tnk = df.loc[df["is_tnk"]].copy()

    def patient_metrics(cells: pd.DataFrame, prefix: str) -> pd.DataFrame:
        rows = []
        for pid, g in cells.groupby("patient"):
            rec = {
                "patient": pid,
                f"{prefix}n_cells": int(len(g)),
                f"{prefix}mean_log1p_cp10k": float(g["tacstd2_log1p_cp10k"].mean()),
                f"{prefix}pseudobulk_cpm": float(
                    g["tacstd2_umi"].sum() / g["total_umi"].sum() * 1e6
                ),
            }
            for spec in CELL_CUTOFFS:
                if spec["kind"] == "umi":
                    pos = g["tacstd2_umi"] >= spec["thr"]
                else:
                    pos = g["tacstd2_log1p_cp10k"] >= spec["thr"]
                rec[f"{prefix}pct_{spec['id']}"] = float(100.0 * pos.mean())
            rows.append(rec)
        return pd.DataFrame(rows)

    mal_pat = patient_metrics(mal, "")
    tnk_pat = patient_metrics(tnk, "tnk_")
    core_mal = mal.loc[mal["is_core_cell"]]
    core_pat = patient_metrics(core_mal, "core_") if len(core_mal) else pd.DataFrame()

    sample = all_pat.merge(mal_pat, on="patient", how="left")
    sample = sample.merge(tnk_pat, on="patient", how="left")
    if len(core_pat):
        sample = sample.merge(core_pat, on="patient", how="left")
    for spec in CELL_CUTOFFS:
        col = f"pct_{spec['id']}"
        if col in sample.columns:
            sample[col] = sample[col].where(sample["n_malignant"] > 0)
    sample["mean_log1p_cp10k"] = sample["mean_log1p_cp10k"].where(sample["n_malignant"] > 0)
    sample["pseudobulk_cpm"] = sample["pseudobulk_cpm"].where(sample["n_malignant"] > 0)

    # Analysis sets
    sets = {
        "primary_drop0_RvsNR": sample["in_primary"],
        "min20_RvsNR": sample["in_primary"] & (sample["n_malignant"] >= 20),
        "nsclc_drop0_RvsNR": sample["in_primary"] & sample["is_nsclc"],
        "core_malignant_RvsNR": sample["in_primary"] & sample["has_core_malignant"],
        "PRvsPD_drop0": (sample["n_malignant"] > 0) & sample["recist"].isin(["PR", "PD"]),
        "nsclc_PRvsPD_drop0": (
            (sample["n_malignant"] > 0)
            & sample["is_nsclc"]
            & sample["recist"].isin(["PR", "PD"])
        ),
    }

    def split_rn(sub: pd.DataFrame, recist_only: bool = False):
        if recist_only:
            r = sub[sub["recist"] == "PR"]
            nr = sub[sub["recist"] == "PD"]
        else:
            r = sub[sub["response"] == "R"]
            nr = sub[sub["response"] == "NR"]
        return r, nr

    continuous_rows = []
    for set_name, mask in sets.items():
        sub = sample.loc[mask].copy()
        recist_only = "PRvsPD" in set_name
        r, nr = split_rn(sub, recist_only=recist_only)
        for metric, label in [
            ("mean_log1p_cp10k", "mean_log1p_cp10k"),
            ("pseudobulk_cpm", "pseudobulk_cpm"),
        ]:
            st = mwu(r[metric].dropna().values, nr[metric].dropna().values)
            continuous_rows.append(
                {
                    "analysis_set": set_name,
                    "metric": label,
                    "cutoff_family": "continuous",
                    **st,
                }
            )
        for spec in CELL_CUTOFFS:
            col = f"pct_{spec['id']}"
            st = mwu(r[col].dropna().values, nr[col].dropna().values)
            continuous_rows.append(
                {
                    "analysis_set": set_name,
                    "metric": f"pct_pos_{spec['id']}",
                    "cutoff_family": "cell_pct_positive",
                    "cell_cutoff": spec["id"],
                    "cell_cutoff_label": spec["label"],
                    **st,
                }
            )
    continuous = pd.DataFrame(continuous_rows)

    fisher_rows = []
    for set_name, mask in sets.items():
        sub = sample.loc[mask].copy()
        recist_only = "PRvsPD" in set_name
        r, nr = split_rn(sub, recist_only=recist_only)
        pct = r["pct_umi_ge_1"].dropna()
        pct_nr = nr["pct_umi_ge_1"].dropna()
        for thr in PATIENT_PCT_CUTOFFS:
            r_high = int((pct >= thr).sum())
            r_low = int((pct < thr).sum())
            nr_high = int((pct_nr >= thr).sum())
            nr_low = int((pct_nr < thr).sum())
            ft = fisher_high(r_high, r_low, nr_high, nr_low)
            fisher_rows.append(
                {
                    "analysis_set": set_name,
                    "patient_pct_cutoff": thr,
                    "definition": f"TROP2-high if %UMI≥1 ≥ {thr:.0f}%",
                    "R_high": r_high,
                    "R_low": r_low,
                    "NR_high": nr_high,
                    "NR_low": nr_low,
                    **ft,
                }
            )
    fisher = pd.DataFrame(fisher_rows)

    # Paired malignant vs T/NK among patients with both >0
    paired = sample[(sample["n_malignant"] > 0) & (sample["n_tnk"] > 0)].copy()
    paired_stats = {}
    for metric_m, metric_t, name in [
        ("mean_log1p_cp10k", "tnk_mean_log1p_cp10k", "mean_log1p_cp10k"),
        ("pct_umi_ge_1", "tnk_pct_umi_ge_1", "pct_umi_ge_1"),
    ]:
        a = paired[metric_m].values
        b = paired[metric_t].values
        w, p = stats.wilcoxon(a, b, alternative="two-sided")
        paired_stats[name] = {
            "n_pairs": int(len(paired)),
            "median_malignant": float(np.median(a)),
            "median_tnk": float(np.median(b)),
            "W": float(w),
            "p": float(p),
        }

    # Spearman malignant TACSTD2 vs T/NK fraction (A3-analog, exploratory)
    sample["tnk_fraction"] = sample["n_tnk"] / (sample["n_tnk"] + sample["n_malignant"])
    spe_rows = []
    for set_name, mask in {
        "drop0_any_response": sample["n_malignant"] > 0,
        "primary_drop0_RvsNR": sample["in_primary"],
        "nsclc_drop0": (sample["n_malignant"] > 0) & sample["is_nsclc"],
    }.items():
        sub = sample.loc[mask & (sample["n_tnk"] > 0)].copy()
        if len(sub) >= 5:
            rho, p = stats.spearmanr(sub["pct_umi_ge_1"], sub["tnk_fraction"])
            spe_rows.append(
                {
                    "analysis_set": set_name,
                    "x": "malignant_pct_umi_ge_1",
                    "y": "TNK_fraction_among_mal_plus_TNK",
                    "n": int(len(sub)),
                    "spearman_rho": float(rho),
                    "p": float(p),
                }
            )
    spearman = pd.DataFrame(spe_rows)

    # Sanity
    sanity = {
        "n_cells": int(len(df)),
        "n_malignant": int(df["is_malignant"].sum()),
        "n_tnk": int(df["is_tnk"].sum()),
        "epcam_pos_malignant": float((mal["epcam_umi"] > 0).mean() * 100),
        "epcam_pos_tnk": float((tnk["epcam_umi"] > 0).mean() * 100),
        "ptprc_pos_malignant": float((mal["ptprc_umi"] > 0).mean() * 100),
        "ptprc_pos_tnk": float((tnk["ptprc_umi"] > 0).mean() * 100),
        "tacstd2_pos_all_cells": float((df["tacstd2_umi"] > 0).mean() * 100),
    }

    # BH-FDR within primary cell-cutoff family only (5 tests)
    prim = continuous[
        (continuous["analysis_set"] == "primary_drop0_RvsNR")
        & (continuous["cutoff_family"] == "cell_pct_positive")
    ].copy()
    pvals = prim["p"].to_numpy(dtype=float)
    if hasattr(stats, "false_discovery_control"):
        fdr = np.asarray(stats.false_discovery_control(pvals, method="bh"), dtype=float)
    else:
        order = np.argsort(pvals)
        ranked = np.empty_like(pvals)
        n = len(pvals)
        prev = 1.0
        for i, idx in enumerate(order[::-1]):
            rank = n - i
            val = min(prev, pvals[idx] * n / rank)
            ranked[idx] = val
            prev = val
        fdr = ranked
    prim = prim.assign(bh_fdr_5cellcutoffs=fdr)
    continuous = continuous.merge(
        prim[["analysis_set", "metric", "bh_fdr_5cellcutoffs"]],
        on=["analysis_set", "metric"],
        how="left",
    )

    # Write tables
    sample_out = sample.copy()
    sample_out.to_csv(OUT / "sample_table.tsv", sep="\t", index=False)
    all_pat.to_csv(OUT / "patient_coverage.tsv", sep="\t", index=False)
    continuous.to_csv(OUT / "cutoff_tests.tsv", sep="\t", index=False)
    fisher.to_csv(OUT / "patient_highlow_fisher.tsv", sep="\t", index=False)
    spearman.to_csv(OUT / "tnk_fraction_spearman.tsv", sep="\t", index=False)
    gsm.to_csv(OUT / "gsm_sample_metadata.csv", index=False)

    # Key pulls
    def row(set_name, metric):
        hit = continuous[
            (continuous["analysis_set"] == set_name) & (continuous["metric"] == metric)
        ]
        return hit.iloc[0].to_dict() if len(hit) else {}

    p_mean = row("primary_drop0_RvsNR", "mean_log1p_cp10k")
    p_pct = row("primary_drop0_RvsNR", "pct_pos_umi_ge_1")
    p_cpm = row("primary_drop0_RvsNR", "pseudobulk_cpm")
    p_nsclc_pct = row("nsclc_drop0_RvsNR", "pct_pos_umi_ge_1")
    p_nsclc_mean = row("nsclc_drop0_RvsNR", "mean_log1p_cp10k")
    p_prpd = row("PRvsPD_drop0", "pct_pos_umi_ge_1")
    p_core = row("core_malignant_RvsNR", "pct_pos_umi_ge_1")

    dropped = all_pat[all_pat["zero_malignant"]]
    n_r_drop = int(dropped["response"].eq("R").sum())
    n_nr_drop = int(dropped["response"].eq("NR").sum())
    n_ne_drop = int(dropped["response"].eq("NE").sum())
    n_between = int(
        ((all_pat["n_malignant"] > 0) & (all_pat["n_malignant"] < 20)).sum()
    )

    # Verdict
    # Locked rule: primary is drop-0 R vs NR on %pos UMI>=1 and mean log1p.
    # Do not promote a cutoff because it is smaller.
    primary_null = (p_pct.get("p", 1) >= 0.05) and (p_mean.get("p", 1) >= 0.05)
    any_cell_sig = bool((prim["p"] < 0.05).any())
    any_fisher_prim = bool(
        (
            fisher.loc[fisher["analysis_set"] == "primary_drop0_RvsNR", "p"] < 0.05
        ).any()
    )
    if mpr_labeled:
        verdict_short = "MPR labeled — see tests"
    elif primary_null and not any_cell_sig and not any_fisher_prim:
        verdict_short = (
            "NO — MPR not labeled; after dropping 0-malignant patients, "
            "malignant TACSTD2 is still null at every pre-specified %positive cutoff"
        )
    elif primary_null:
        verdict_short = (
            "NO on the primary metrics; a secondary cutoff/set is nominally p<0.05 "
            "and is not the verdict"
        )
    else:
        verdict_short = "YES — primary drop-0 R vs NR is significant"

    # Figures
    prim_sample = sample.loc[sample["n_malignant"] > 0].copy()
    colors = {"R": "#2166ac", "NR": "#b2182b", "NE": "#777777"}
    rng = np.random.default_rng(1)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6))
    # A: %pos UMI>=1
    ax = axes[0]
    order = ["R", "NR", "NE"]
    for i, g in enumerate(order):
        v = prim_sample.loc[prim_sample["response"] == g, "pct_umi_ge_1"].dropna().values
        if len(v) == 0:
            continue
        ax.scatter(rng.normal(i, 0.07, len(v)), v, c=colors[g], s=55, zorder=3, alpha=0.9)
        ax.hlines(np.median(v), i - 0.22, i + 0.22, color="black", lw=2, zorder=4)
    ax.set_xticks(range(3))
    ax.set_xticklabels(
        [f"{g}\n(n={(prim_sample['response']==g).sum()})" for g in order]
    )
    ax.set_ylabel("% TACSTD2+ malignant cells (UMI ≥ 1)")
    ax.set_title("Primary: drop 0-malignant patients\nRECIST R vs NR (not MPR)")
    ax.set_ylim(-5, 105)

    # B: mean log1p
    ax = axes[1]
    for i, g in enumerate(order):
        v = prim_sample.loc[prim_sample["response"] == g, "mean_log1p_cp10k"].dropna().values
        if len(v) == 0:
            continue
        ax.scatter(rng.normal(i, 0.07, len(v)), v, c=colors[g], s=55, zorder=3, alpha=0.9)
        ax.hlines(np.median(v), i - 0.22, i + 0.22, color="black", lw=2, zorder=4)
    ax.set_xticks(range(3))
    ax.set_xticklabels(
        [f"{g}\n(n={(prim_sample['response']==g).sum()})" for g in order]
    )
    ax.set_ylabel("TACSTD2 mean log1p(CP10K)")
    ax.set_title(f"Same patients, mean expression\n(original hunt p={p_mean.get('p', float('nan')):.2f})")

    # C: cutoff p-values
    ax = axes[2]
    prim_cut = continuous[
        (continuous["analysis_set"] == "primary_drop0_RvsNR")
        & (continuous["cutoff_family"] == "cell_pct_positive")
    ]
    ys = prim_cut["p"].values
    labels = [s.replace("pct_pos_", "") for s in prim_cut["metric"]]
    ax.barh(range(len(ys)), ys, color="#4d4d4d")
    ax.axvline(0.05, color="#b2182b", ls="--", lw=1)
    ax.set_yticks(range(len(ys)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Mann–Whitney p (R vs NR)")
    ax.set_title("Pre-specified cell %positive cutoffs\n(primary set)")
    ax.set_xlim(0, 1)
    fig.tight_layout()
    fig.savefig(FIG / "primary_r_vs_nr_and_cutoffs.png", dpi=200)
    plt.close(fig)

    # Per-patient bars
    fig, ax = plt.subplots(figsize=(12, 4.8))
    mo = prim_sample.sort_values(["response", "pct_umi_ge_1"], ascending=[True, False])
    ax.bar(
        range(len(mo)),
        mo["pct_umi_ge_1"],
        color=[colors[r] for r in mo["response"]],
    )
    ax.set_xticks(range(len(mo)))
    ax.set_xticklabels(
        [f"{p}\n{rec}/{sub}" for p, rec, sub in zip(mo["patient"], mo["recist"], mo["cancer_subtype"])],
        fontsize=7,
    )
    ax.set_ylabel("% TACSTD2+ malignant cells (UMI ≥ 1)")
    ax.set_title("Per-patient malignant %positive (blue=R, red=NR, grey=NE). SCLC/NUT labeled.")
    ax.set_ylim(0, 100)
    fig.tight_layout()
    fig.savefig(FIG / "per_patient_pctpos_bars.png", dpi=200)
    plt.close(fig)

    # High/low mosaic at 50%
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    f50 = fisher[
        (fisher["analysis_set"] == "primary_drop0_RvsNR")
        & (fisher["patient_pct_cutoff"] == 50.0)
    ].iloc[0]
    mat = np.array(
        [[f50["R_high"], f50["R_low"]], [f50["NR_high"], f50["NR_low"]]], dtype=float
    )
    im = ax.imshow(mat, cmap="Blues")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["TROP2-high\n(≥50% +)", "TROP2-low\n(<50% +)"])
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["R (PR)", "NR (SD/PD)"])
    for i in range(2):
        for j in range(2):
            ax.text(j, i, int(mat[i, j]), ha="center", va="center", fontsize=14)
    ax.set_title(f"Patient-level 50% cutoff\nFisher p={fmt_p(f50['p'])}")
    fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(FIG / "patient_highlow_50pct.png", dpi=200)
    plt.close(fig)

    # Provenance
    files_hashed = {}
    for p in [
        DATA / "GSE205335_Lung_IO_CellIdentity.txt.gz",
        DATA / "percell_tacstd2.csv.gz",
    ]:
        if p.exists():
            files_hashed[p.name] = {"bytes": p.stat().st_size, "md5": md5sum(p)}
    rds = DATA / "GSE205335_Lung_IO_UMI_matrix.rds"
    rds_gz = DATA / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    for p in [rds, rds_gz]:
        if p.exists():
            files_hashed[p.name] = {"bytes": p.stat().st_size, "md5": md5sum(p)}

    provenance = {
        "accession": "GSE205335",
        "series_url": GEO_SERIES,
        "paper": (
            "Ahn/Lee et al., eLife reviewed preprint 98366 "
            "(Unveiling the influence of tumor and immune signatures "
            "on immune checkpoint therapy in advanced lung cancer)"
        ),
        "processed_files": {
            "cell_identity": GEO_CELL,
            "umi_matrix": GEO_RDS,
        },
        "raw_not_used": "EGAD00001008703 (controlled; not downloaded)",
        "hashes": files_hashed,
        "cell_cutoffs_locked": CELL_CUTOFFS,
        "patient_pct_cutoffs_locked": PATIENT_PCT_CUTOFFS,
        "primary_cell_cutoff": PRIMARY_CELL_CUTOFF,
        "response_map": RESP_MAP,
        "mpr_labeled": mpr_labeled,
        "mpr_hit_columns": mpr_hits,
        "malignant_definition": "authors lineage.sub == 'Malignant cells' (no CNV re-inference)",
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")

    summary = {
        "verdict": verdict_short,
        "mpr_labeled": mpr_labeled,
        "n_patients_geo": int(all_pat["patient"].nunique()),
        "n_dropped_zero_malignant": int(all_pat["zero_malignant"].sum()),
        "dropped_zero_malignant": dropped["patient"].tolist(),
        "dropped_asymmetric": {
            "R": n_r_drop,
            "NR": n_nr_drop,
            "NE": n_ne_drop,
        },
        "n_patients_1_to_19_malignant": n_between,
        "primary": {
            "set": "patients with >0 author-annotated malignant cells; RECIST R vs NR; NE excluded",
            "n_R": p_pct.get("n_R"),
            "n_NR": p_pct.get("n_NR"),
            "pct_pos_umi_ge_1": p_pct,
            "mean_log1p_cp10k": p_mean,
            "pseudobulk_cpm": p_cpm,
        },
        "nsclc_only": {
            "pct_pos_umi_ge_1": p_nsclc_pct,
            "mean_log1p_cp10k": p_nsclc_mean,
        },
        "core_malignant": p_core,
        "PRvsPD": p_prpd,
        "paired_malignant_vs_tnk": paired_stats,
        "spearman_tnk_fraction": spe_rows,
        "sanity": sanity,
        "any_primary_cell_cutoff_p_lt_0.05": any_cell_sig,
        "any_primary_fisher_p_lt_0.05": any_fisher_prim,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # REPORT
    drop_list = ", ".join(
        f"{r.patient} ({r.recist}, {r.cancer_subtype}, {r.tissues})"
        for r in dropped.itertuples()
    )
    cell_lines = []
    for spec in CELL_CUTOFFS:
        st = row("primary_drop0_RvsNR", f"pct_pos_{spec['id']}")
        fdrv = prim.loc[prim["metric"] == f"pct_pos_{spec['id']}", "bh_fdr_5cellcutoffs"]
        fdr_s = fmt_p(float(fdrv.iloc[0])) if len(fdrv) else "NA"
        cell_lines.append(
            f"| {spec['label']} | {fmt_num(st.get('median_R', float('nan')), 1)} | "
            f"{fmt_num(st.get('median_NR', float('nan')), 1)} | "
            f"{fmt_p(st.get('p', float('nan')))} | {fdr_s} | "
            f"{fmt_num(st.get('cliffs_delta_R_minus_NR', float('nan')))} |"
        )
    fish_lines = []
    for thr in PATIENT_PCT_CUTOFFS:
        ft = fisher[
            (fisher["analysis_set"] == "primary_drop0_RvsNR")
            & (fisher["patient_pct_cutoff"] == thr)
        ].iloc[0]
        fish_lines.append(
            f"| ≥ {thr:.0f}% UMI≥1 | {int(ft.R_high)}/{int(ft.R_high+ft.R_low)} | "
            f"{int(ft.NR_high)}/{int(ft.NR_high+ft.NR_low)} | "
            f"{fmt_num(ft.odds_ratio)} | {fmt_p(ft.p)} |"
        )
    set_lines = []
    for set_name, label in [
        ("primary_drop0_RvsNR", "All histologies, drop 0 malignant, R vs NR"),
        ("min20_RvsNR", "≥20 malignant cells, R vs NR"),
        ("nsclc_drop0_RvsNR", "ADC+SQ only, drop 0, R vs NR"),
        ("core_malignant_RvsNR", "Author Core malignant cells, R vs NR"),
        ("PRvsPD_drop0", "All histologies, drop 0, PR vs PD"),
        ("nsclc_PRvsPD_drop0", "ADC+SQ only, drop 0, PR vs PD"),
    ]:
        st = row(set_name, "pct_pos_umi_ge_1")
        set_lines.append(
            f"| {label} | {st.get('n_R','')} vs {st.get('n_NR','')} | "
            f"{fmt_num(st.get('median_R', float('nan')), 1)} | "
            f"{fmt_num(st.get('median_NR', float('nan')), 1)} | "
            f"{fmt_p(st.get('p', float('nan')))} |"
        )
    spe_lines = []
    for r in spe_rows:
        spe_lines.append(
            f"| {r['analysis_set']} | {r['n']} | {fmt_num(r['spearman_rho'])} | {fmt_p(r['p'])} |"
        )

    report = f"""# REWORK A3 — GSE205335 malignant TACSTD2: MPR vs NMPR, 0-cell drop, %positive cutoffs

**Self-contained. Public processed GEO files only. Written to be read without the rest of the repo.**

**Verdict: {verdict_short}**

GSE205335 is advanced / palliative ICI, labeled with **RECIST 1.1 (PR / SD / PD / NE) only**.
There is **no MPR or NMPR field** in GEO sample characteristics or in the eLife paper.
The requested MPR vs NMPR contrast **cannot be run**. After dropping the 4 patients
with 0 author-annotated malignant cells, malignant TACSTD2 is still null for
responders vs non-responders: mean log1p(CP10K) p = {fmt_p(p_mean.get('p', float('nan')))}
(this is the published hunt p ≈ 0.96), %positive (UMI ≥ 1) p = {fmt_p(p_pct.get('p', float('nan')))}.
Every pre-specified cell-level %positive cutoff and every patient-level high/low
cutoff is also null on the primary set. An ADC+SQ-only sensitivity at UMI ≥ 1
is nominally p = {fmt_p(p_nsclc_pct.get('p', float('nan')))} (n = 4 vs 8) and
is **not** used as the verdict. This is an **underpowered primary null
(n = {p_pct.get('n_R')} R vs {p_pct.get('n_NR')} NR), not proof of no association**.

## Why this rework exists

A prior hunt (`results/hunt_gse205335/`) reported malignant-cell TACSTD2
R vs NR p = 0.96 (mean log1p CP10K; 6 vs 10 patients). Claim A3 in this project
is actually about **GSE207422** (neoadjuvant NSCLC, MPR vs NMPR). This folder
asks whether the GSE205335 null changes if we:

1. use **MPR vs NMPR if labeled**,
2. **drop patients with 0 malignant cells**,
3. score TACSTD2 as **%positive at pre-specified cutoffs** instead of only the mean.

## MPR vs NMPR — not labeled

| Source | Pathologic response (MPR / NMPR)? | What is labeled |
|---|---|---|
| GEO `!Sample_characteristics_ch1` (33 GSM, fetched for this run) | **No** | `recist`: PR / SD / PD / NE |
| eLife preprint 98366 (Ahn / Lee) | **No** | RECIST 1.1; R = PR, NR = SD+PD; no CR |
| Author cell table `GSE205335_Lung_IO_CellIdentity.txt` | **No** | cell lineage only |

This is a **stage III–IV / extensive-disease biopsy cohort on palliative ICI**,
not a neoadjuvant resection cohort. MPR requires a resected primary and a
percent residual viable tumor. That endpoint does not exist here. Substituting
RECIST for MPR would be a silent endpoint swap; it is not done.

Paper coding, reused: **R = PR** (no CR), **NR = SD + PD**, **NE excluded**
from response tests. Sensitivity: PR vs PD (drop SD).

## Analysis set and the 0-malignant drop

Author malignant call: `lineage.sub == "Malignant cells"` (n = {sanity['n_malignant']:,} cells).
CNV was **not** re-inferred. Patient is the unit; samples from the same patient
are pooled.

| | n |
|---|---:|
| Patients in GEO | {int(all_pat['patient'].nunique())} |
| Dropped (0 malignant cells) | {int(all_pat['zero_malignant'].sum())} |
| Patients with 1–19 malignant cells | {n_between} |
| Remaining with malignant cells | {int((all_pat['n_malignant']>0).sum())} |
| Of those, R / NR / NE | {int(((all_pat['n_malignant']>0)&(all_pat['response']=='R')).sum())} / {int(((all_pat['n_malignant']>0)&(all_pat['response']=='NR')).sum())} / {int(((all_pat['n_malignant']>0)&(all_pat['response']=='NE')).sum())} |
| Primary test (drop 0, R vs NR) | {p_pct.get('n_R')} vs {p_pct.get('n_NR')} |

Dropped zero-malignant patients: {drop_list}.

Dropout is **asymmetric**: {n_r_drop} responder vs {n_nr_drop} non-responder
({n_ne_drop} NE). Three of the four are **normal tissue** (normal LN / normal brain)
from early-stage or metastatic patients — they have no tumor cells by sampling,
not by a failed annotation. Absence of captured tumor in some PR patients can
also be treatment or sampling; it shrinks the R arm and is not imputed.

**Drop-0 and the original hunt's ≥20-cell floor are the same set** in this
cohort: nobody has 1–19 malignant cells (P4001 is the smallest remaining, 27 cells).
Rework item 2 therefore **does not change n or the p = 0.96 mean-expression test**.
It only makes the exclusion rule explicit.

## Pre-specified %positive cutoffs

Locked before looking at cutoff p-values. A malignant cell is TACSTD2+ if it
meets the threshold. Patient score = 100 × (positive malignant cells / malignant cells).
Primary cutoff = **UMI ≥ 1** (detection). The other four are robustness, not a
search for a significant threshold.

Cell-level thresholds → patient %positive → two-sided Mann–Whitney U, R vs NR:

| Cell cutoff | median %pos R | median %pos NR | p | BH-FDR (5 cutoffs) | Cliff's δ (R−NR) |
|---|---:|---:|---:|---:|---:|
{chr(10).join(cell_lines)}

Patient-level TROP2-high (using % UMI ≥ 1) → two-sided Fisher exact:

| High if | R high / R | NR high / NR | OR (high in R vs NR) | p |
|---|---:|---:|---:|---:|
{chr(10).join(fish_lines)}

Continuous metrics on the same primary patients:

| Metric | median R | median NR | p | Cliff's δ |
|---|---:|---:|---:|---:|
| mean log1p(CP10K) | {fmt_num(p_mean.get('median_R', float('nan')))} | {fmt_num(p_mean.get('median_NR', float('nan')))} | {fmt_p(p_mean.get('p', float('nan')))} | {fmt_num(p_mean.get('cliffs_delta_R_minus_NR', float('nan')))} |
| %pos UMI ≥ 1 | {fmt_num(p_pct.get('median_R', float('nan')), 1)} | {fmt_num(p_pct.get('median_NR', float('nan')), 1)} | {fmt_p(p_pct.get('p', float('nan')))} | {fmt_num(p_pct.get('cliffs_delta_R_minus_NR', float('nan')))} |
| pseudobulk CPM | {fmt_num(p_cpm.get('median_R', float('nan')), 1)} | {fmt_num(p_cpm.get('median_NR', float('nan')), 1)} | {fmt_p(p_cpm.get('p', float('nan')))} | {fmt_num(p_cpm.get('cliffs_delta_R_minus_NR', float('nan')))} |

No primary cell-cutoff p < 0.05. No primary Fisher p < 0.05. Cliff's δ values
are small. Do not quote a "best" cutoff.

## Sensitivity (not used for the verdict)

| Set | n R vs NR | median %pos R | median %pos NR | p |
|---|---|---:|---:|---:|
{chr(10).join(set_lines)}

NSCLC-only (ADC+SQ) removes two SCLC PR patients with very low malignant
TACSTD2 (P1016, P1115) plus SCLC PD (P1025) and NUT SD (P1056). That
**raises** the R-arm median %positive. At UMI ≥ 1 this sensitivity is
**nominally p = {fmt_p(p_nsclc_pct.get('p', float('nan')))}**
(n = {p_nsclc_pct.get('n_R')} vs {p_nsclc_pct.get('n_NR')};
Cliff's δ = {fmt_num(p_nsclc_pct.get('cliffs_delta_R_minus_NR', float('nan')))}).
It is **not the verdict**: n = 4 vs 8, uncorrected across six analysis sets
and five cutoffs; UMI ≥ 2 in the same NSCLC set is p = 0.214; NSCLC PR vs PD
is p = {fmt_p(row('nsclc_PRvsPD_drop0', 'pct_pos_umi_ge_1').get('p', float('nan')))};
mean log1p in NSCLC is p = {fmt_p(p_nsclc_mean.get('p', float('nan')))}.
A single n=4 Mann–Whitney that kisses 0.05 after dropping neuroendocrine
tumors is a hypothesis, not a finding.

Author `core.patient == Core` is the paper's 14-sample / 11-patient clinical
core. It is a sensitivity, not the primary, because the rework asked for all
patients after a 0-malignant drop.

## Exploratory: T/NK (not the rework question)

Paired malignant vs T/NK TACSTD2 remains large in every patient with both
compartments (n = {paired_stats['pct_umi_ge_1']['n_pairs']}): median %pos
{fmt_num(paired_stats['pct_umi_ge_1']['median_malignant'], 1)}% vs
{fmt_num(paired_stats['pct_umi_ge_1']['median_tnk'], 1)}%, Wilcoxon
p = {fmt_p(paired_stats['pct_umi_ge_1']['p'])}. TROP2 is tumor-restricted
here regardless of ICI response.

Claim A3's GSE207422 number was ρ ≈ −0.40 to −0.50 between malignant TACSTD2
and T/NK presence. In GSE205335, Spearman of malignant %pos (UMI ≥ 1) vs
T/NK fraction among (malignant + T/NK):

| Set | n | ρ | p |
|---|---:|---:|---:|
{chr(10).join(spe_lines)}

That is **not** the −0.40 to −0.50 claim, and this cohort is not GSE207422.

## Marker sanity

| Marker | % positive malignant | % positive T/NK |
|---|---:|---:|
| EPCAM (UMI > 0) | {sanity['epcam_pos_malignant']:.1f} | {sanity['epcam_pos_tnk']:.1f} |
| PTPRC / CD45 (UMI > 0) | {sanity['ptprc_pos_malignant']:.1f} | {sanity['ptprc_pos_tnk']:.1f} |

Consistent with the authors' malignant vs lymphocyte labels.

## Honest interpretation

1. **MPR vs NMPR is impossible here.** Do not write "MPR" on GSE205335 figures.
2. **Dropping 0-malignant patients does not move p = 0.96.** Those four patients
   already contributed no tumor-cell score; no one sits between 1 and 19 cells.
3. **On the primary set, %positive cutoffs do not rescue a response
   association.** Detection, UMI ≥ 2/3, and two log1p thresholds are all
   null; so are 10/25/50/75% patient high/low splits. This was pre-specified,
   not p-hacked.
4. **n = 6 vs 10 (4 vs 8 if NSCLC-only) is underpowered.** Only a large effect
   could have been seen. A null here is inconclusive, not "TROP2 is unrelated
   to ICI".
5. **SCLC in the R arm is a real confounder.** Two PR SCLC tumors are
   TACSTD2-low, as expected for neuroendocrine histology. Restricting to
   ADC+SQ makes UMI ≥ 1 nominally p = 0.048 (4 vs 8). That is a sensitivity,
   not a confirmed NSCLC effect: other cutoffs and PR-vs-PD in the same
   subset are not significant, and n = 4 cannot carry a claim.
6. **Asymmetric tumor-cell dropout.** 3 PR vs 1 PD had zero malignant cells.
   Two PR dropouts are normal LN. Do not treat the remaining R arm as a
   random sample of responders.
7. **mRNA ≠ protein ≠ ADC target occupancy.** UMI ≥ 1 is a transcript detection
   call, not an IHC H-score.
8. **Mixed sites and 3′/5′ chemistry** are confounded with patient; n is too
   small to adjust.
9. **Malignant labels are the authors'.** No independent inferCNV in this rework.

## Reproduce

```
pip install -r requirements.txt
# R with Matrix is required for the RDS extract
Rscript scripts/extract_gse205335_tacstd2.R
python scripts/rework_A3_GSE205335.py
```

The Python script will download the two GEO processed supplements into
`data/GSE205335/` (~500 MB UMI matrix, gitignored) if they are missing, gunzip
the RDS once (GEO file is double-gzipped), and call the R extract.

## Files

- `sample_table.tsv` — one row per patient: RECIST, histology, n malignant, all %pos cutoffs
- `patient_coverage.tsv` — the 0-cell drop and who is in the primary set
- `cutoff_tests.tsv` — Mann–Whitney for every set × metric
- `patient_highlow_fisher.tsv` — Fisher exact at 10/25/50/75%
- `tnk_fraction_spearman.tsv` — exploratory ρ
- `gsm_sample_metadata.csv` — GEO characteristics as fetched
- `summary.json` / `provenance.json`
- `figures/primary_r_vs_nr_and_cutoffs.png`
- `figures/per_patient_pctpos_bars.png`
- `figures/patient_highlow_50pct.png`

## Data

- GEO **GSE205335** processed: `GSE205335_Lung_IO_UMI_matrix.rds.gz` (dgCMatrix,
  33,714 genes × 96,505 cells, raw UMI) and `GSE205335_Lung_IO_CellIdentity.txt.gz`.
- Sample characteristics: NCBI GEO GSM6210624–GSM6210656 (`recist`, subtype, tissue).
- Raw FASTQ is EGA **EGAD00001008703** (controlled). Not used.
"""
    (OUT / "REPORT.md").write_text(report)
    print(json.dumps({"verdict": verdict_short, "primary_pct": p_pct, "primary_mean": p_mean}, indent=2))
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
