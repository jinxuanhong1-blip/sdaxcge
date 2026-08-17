#!/usr/bin/env python3
"""IMvigor210 urothelial analog: CLDN4 vs CD274 and vs atezolizumab response.

Additive. Public processed RNA from IMvigor210CoreBiologies 1.0.0
(Mariathasan et al., Nature 2018). Metastatic urothelial carcinoma
(anti-PD-L1), not lung. Label as analog.

Usage:
  python3 methods/imvigor210_cldn4/analyze.py
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import tarfile
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
CACHE = Path("/tmp/imvigor210_cldn4")

PKG_URL = (
    "http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/"
    "IMvigor210CoreBiologies_1.0.0.tar.gz"
)
PKG_NAME = "IMvigor210CoreBiologies_1.0.0.tar.gz"
EXPECTED_SIZE = 122_127_298
EXPECTED_SHA256 = "cfdd3176d7b34de5b04fb9416bfd2b20fa4b6e238aaad5f20b048a34329ea178"
N_BOOT = 2000
SEED = 20260817

GENES = ["CLDN4", "CD274", "CD8A", "CXCL9"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def find_tarball() -> Path | None:
    candidates = [
        CACHE / PKG_NAME,
        Path("/tmp/imvigor210") / PKG_NAME,
        Path("/tmp") / PKG_NAME,
    ]
    for p in candidates:
        if p.is_file() and p.stat().st_size == EXPECTED_SIZE:
            return p
    return None


def download_tarball() -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / PKG_NAME
    existing = find_tarball()
    if existing is not None:
        if existing.resolve() != dest.resolve():
            dest.write_bytes(existing.read_bytes())
        digest = sha256(dest)
        if digest != EXPECTED_SHA256:
            raise RuntimeError(f"Unexpected sha256 {digest}")
        return dest
    print(f"Downloading {PKG_URL}")
    urllib.request.urlretrieve(PKG_URL, dest)
    if dest.stat().st_size != EXPECTED_SIZE:
        raise RuntimeError(f"Unexpected size {dest.stat().st_size}")
    digest = sha256(dest)
    if digest != EXPECTED_SHA256:
        raise RuntimeError(f"Unexpected sha256 {digest}")
    return dest


def extract_sample_table(tarball: Path) -> pd.DataFrame:
    TABLES.mkdir(parents=True, exist_ok=True)
    out_tsv = TABLES / "sample_level.tsv"
    rdata = CACHE / "cds.RData"
    CACHE.mkdir(parents=True, exist_ok=True)
    if not rdata.is_file() or rdata.stat().st_size == 0:
        member = "IMvigor210CoreBiologies/data/cds.RData"
        with tarfile.open(tarball, "r:gz") as tf:
            src = tf.extractfile(member)
            if src is None:
                raise RuntimeError(f"{member} missing")
            rdata.write_bytes(src.read())
    rscript = HERE / "extract_cds.R"
    cmd = ["Rscript", str(rscript), str(rdata), str(out_tsv)]
    print("Running", " ".join(cmd))
    subprocess.check_call(cmd)
    df = pd.read_csv(out_tsv, sep="\t")
    if df.shape[0] != 348:
        raise RuntimeError(f"Expected 348 RNA samples, got {df.shape[0]}")
    return df


def fisher_or(high: np.ndarray, resp: np.ndarray) -> dict:
    # 2x2: rows = low/high, cols = SD/PD / CR/PR
    a = int(((~high) & (~resp)).sum())  # low, NR
    b = int(((~high) & resp).sum())  # low, R
    c = int((high & (~resp)).sum())  # high, NR
    d = int((high & resp).sum())  # high, R
    table = np.array([[a, b], [c, d]], dtype=int)
    or_, p = stats.fisher_exact(table, alternative="two-sided")
    # Woolf logit CI
    if min(a, b, c, d) == 0:
        lo = hi = float("nan")
    else:
        se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        lo = math.exp(math.log(or_) - 1.96 * se)
        hi = math.exp(math.log(or_) + 1.96 * se)
    return {
        "n": int(len(high)),
        "n_low": int((~high).sum()),
        "n_high": int(high.sum()),
        "r_low": b,
        "nr_low": a,
        "r_high": d,
        "nr_high": c,
        "orr_low": b / (a + b) if (a + b) else float("nan"),
        "orr_high": d / (c + d) if (c + d) else float("nan"),
        "or": float(or_),
        "or_lo": float(lo),
        "or_hi": float(hi),
        "p": float(p),
    }


def logistic_or(x: np.ndarray, y: np.ndarray) -> dict:
    z = (x - np.mean(x)) / np.std(x, ddof=1)
    X = sm.add_constant(z)
    fit = sm.Logit(y, X).fit(disp=False)
    est = float(fit.params[1])
    se = float(fit.bse[1])
    return {
        "n": int(len(y)),
        "or_per_sd": math.exp(est),
        "or_lo": math.exp(est - 1.96 * se),
        "or_hi": math.exp(est + 1.96 * se),
        "p": float(fit.pvalues[1]),
    }


def auc_mwu(values: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    pos = values[labels == 1]
    neg = values[labels == 0]
    res = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    auc = float(res.statistic / (len(pos) * len(neg)))
    return auc, float(res.pvalue)


def bootstrap_auc_ci(
    values: np.ndarray, labels: np.ndarray, n_boot: int = N_BOOT, seed: int = SEED
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    pos_idx = np.where(labels == 1)[0]
    neg_idx = np.where(labels == 0)[0]
    aucs = np.empty(n_boot)
    for i in range(n_boot):
        idx = np.concatenate(
            [
                rng.choice(pos_idx, len(pos_idx), replace=True),
                rng.choice(neg_idx, len(neg_idx), replace=True),
            ]
        )
        aucs[i], _ = auc_mwu(values[idx], labels[idx])
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return float(lo), float(hi)


def spearman_row(x: np.ndarray, y: np.ndarray, pair: str) -> dict:
    mask = np.isfinite(x) & np.isfinite(y)
    xx, yy = x[mask], y[mask]
    rho, p = stats.spearmanr(xx, yy)
    n = int(len(xx))
    # Fisher-z CI
    if n > 3 and abs(rho) < 1:
        z = np.arctanh(rho)
        se = 1.0 / math.sqrt(n - 3)
        lo = math.tanh(z - 1.96 * se)
        hi = math.tanh(z + 1.96 * se)
    else:
        lo = hi = float("nan")
    return {
        "pair": pair,
        "n": n,
        "rho": float(rho),
        "p": float(p),
        "rho_lo": float(lo),
        "rho_hi": float(hi),
    }


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> dict:
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    xx, yy, zz = x[mask], y[mask], z[mask]
    rx = stats.rankdata(xx)
    ry = stats.rankdata(yy)
    rz = stats.rankdata(zz)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    rz = rz - rz.mean()
    # residualise on z ranks
    bxz = np.dot(rx, rz) / np.dot(rz, rz)
    byz = np.dot(ry, rz) / np.dot(rz, rz)
    ex = rx - bxz * rz
    ey = ry - byz * rz
    rho, p = stats.pearsonr(ex, ey)
    n = int(len(xx))
    if n > 4 and abs(rho) < 1:
        zf = np.arctanh(rho)
        se = 1.0 / math.sqrt(n - 4)
        lo = math.tanh(zf - 1.96 * se)
        hi = math.tanh(zf + 1.96 * se)
    else:
        lo = hi = float("nan")
    return {
        "n": n,
        "rho": float(rho),
        "p": float(p),
        "rho_lo": float(lo),
        "rho_hi": float(hi),
    }


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def fmt_num(x: float, d: int = 3) -> str:
    if not np.isfinite(x):
        return "NA"
    return f"{x:.{d}f}"


def dedup_patients(df: pd.DataFrame) -> pd.DataFrame:
    """One row per ANONPT_ID. Keep larger sizeFactor if duplicated."""
    d = df.sort_values("sizeFactor", ascending=False).drop_duplicates(
        "ANONPT_ID", keep="first"
    )
    return d.reset_index(drop=True)


def add_flags(d: pd.DataFrame) -> pd.DataFrame:
    out = d.copy()
    out["orr_evaluable"] = out["binaryResponse"].isin(["CR/PR", "SD/PD"])
    out["responder"] = np.where(
        out["binaryResponse"] == "CR/PR",
        1,
        np.where(out["binaryResponse"] == "SD/PD", 0, np.nan),
    )
    out["tissue_bladder"] = out["Tissue"].eq("bladder")
    for g in GENES:
        med = out[f"{g}_log2TPM1"].median()
        out[f"{g}_high"] = out[f"{g}_log2TPM1"] >= med
        out[f"{g}_median"] = med
    return out


def orr_block(sub: pd.DataFrame, gene: str, subset: str) -> dict:
    x = sub[f"{gene}_log2TPM1"].to_numpy(dtype=float)
    y = sub["responder"].to_numpy(dtype=int)
    auc, p_mwu = auc_mwu(x, y)
    lo, hi = bootstrap_auc_ci(x, y)
    high = sub[f"{gene}_high"].to_numpy(dtype=bool)
    fish = fisher_or(high, y.astype(bool))
    logi = logistic_or(x, y)
    q1 = sub[f"{gene}_log2TPM1"] <= sub[f"{gene}_log2TPM1"].quantile(0.25)
    q4 = sub[f"{gene}_log2TPM1"] >= sub[f"{gene}_log2TPM1"].quantile(0.75)
    q = sub.loc[q1 | q4]
    q_high = q[f"{gene}_log2TPM1"] >= q[f"{gene}_log2TPM1"].median()
    # Q4 vs Q1 using the subset quantiles
    q_high = q[f"{gene}_log2TPM1"] >= sub[f"{gene}_log2TPM1"].quantile(0.75)
    q_y = q["responder"].to_numpy(dtype=int)
    q_fish = fisher_or(q_high.to_numpy(dtype=bool), q_y.astype(bool))
    return {
        "subset": subset,
        "gene": gene,
        "n": int(len(sub)),
        "n_CRPR": int(y.sum()),
        "n_SDPD": int((y == 0).sum()),
        "median_CRPR": float(np.median(x[y == 1])),
        "median_SDPD": float(np.median(x[y == 0])),
        "auc": auc,
        "auc_lo": lo,
        "auc_hi": hi,
        "mwu_p": p_mwu,
        "fisher_or": fish["or"],
        "fisher_or_lo": fish["or_lo"],
        "fisher_or_hi": fish["or_hi"],
        "fisher_p": fish["p"],
        "orr_high": fish["orr_high"],
        "orr_low": fish["orr_low"],
        "r_high": fish["r_high"],
        "n_high": fish["n_high"],
        "r_low": fish["r_low"],
        "n_low": fish["n_low"],
        "logit_or_per_sd": logi["or_per_sd"],
        "logit_or_lo": logi["or_lo"],
        "logit_or_hi": logi["or_hi"],
        "logit_p": logi["p"],
        "q4q1_or": q_fish["or"],
        "q4q1_or_lo": q_fish["or_lo"],
        "q4q1_or_hi": q_fish["or_hi"],
        "q4q1_p": q_fish["p"],
        "q4q1_n": q_fish["n"],
        "q4q1_r_high": q_fish["r_high"],
        "q4q1_n_high": q_fish["n_high"],
        "q4q1_r_low": q_fish["r_low"],
        "q4q1_n_low": q_fish["n_low"],
    }


def write_figures(d: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    ev = d[d["orr_evaluable"]].copy()
    ev["resp_lab"] = np.where(ev["responder"] == 1, "CR/PR", "SD/PD")

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    ax.scatter(
        ev["CLDN4_log2TPM1"],
        ev["CD274_log2TPM1"],
        c=np.where(ev["responder"] == 1, "#c0392b", "#2980b9"),
        s=18,
        alpha=0.7,
        linewidths=0,
    )
    rho, p = stats.spearmanr(ev["CLDN4_log2TPM1"], ev["CD274_log2TPM1"])
    ax.set_xlabel("CLDN4 log2(TPM+1)")
    ax.set_ylabel("CD274 log2(TPM+1)")
    ax.set_title(
        f"IMvigor210 analog (urothelial, not lung)\n"
        f"ORR-evaluable n={len(ev)}; Spearman ρ={rho:.3f} p={fmt_p(p)}"
    )
    ax.scatter([], [], c="#c0392b", s=18, label="CR/PR")
    ax.scatter([], [], c="#2980b9", s=18, label="SD/PD")
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd274.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.2), sharey=False)
    for ax, gene in zip(axes, ["CLDN4", "CD274"]):
        data = [
            ev.loc[ev["responder"] == 0, f"{gene}_log2TPM1"],
            ev.loc[ev["responder"] == 1, f"{gene}_log2TPM1"],
        ]
        bp = ax.boxplot(data, tick_labels=["SD/PD", "CR/PR"], patch_artist=True, widths=0.55)
        for patch, color in zip(bp["boxes"], ["#85c1e9", "#f5b7b1"]):
            patch.set_facecolor(color)
        auc, p = auc_mwu(
            ev[f"{gene}_log2TPM1"].to_numpy(), ev["responder"].to_numpy(dtype=int)
        )
        ax.set_title(f"{gene} vs response\nAUC={auc:.3f} p={fmt_p(p)}")
        ax.set_ylabel(f"{gene} log2(TPM+1)")
    fig.suptitle("IMvigor210 analog — urothelial atezolizumab (not lung)", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cldn4_cd274_by_response.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    rows = []
    for gene in ["CLDN4", "CD274", "CD8A", "CXCL9"]:
        high = ev[f"{gene}_high"].to_numpy(dtype=bool)
        y = ev["responder"].to_numpy(dtype=bool)
        fish = fisher_or(high, y)
        rows.append((gene, fish))
    ys = np.arange(len(rows))
    ors = [r[1]["or"] for r in rows]
    los = [r[1]["or_lo"] for r in rows]
    his = [r[1]["or_hi"] for r in rows]
    ax.errorbar(
        ors,
        ys,
        xerr=[np.array(ors) - np.array(los), np.array(his) - np.array(ors)],
        fmt="o",
        color="#1a5276",
        capsize=3,
    )
    ax.axvline(1.0, color="0.4", lw=1)
    ax.set_yticks(ys)
    ax.set_yticklabels(
        [
            f"{g} high vs low\nORR {r['orr_high']*100:.1f}% vs {r['orr_low']*100:.1f}%"
            for g, r in rows
        ]
    )
    ax.set_xlabel("Fisher OR (median-high vs low) for CR/PR")
    ax.set_title(f"IMvigor210 analog ORR n={len(ev)} (not lung)")
    ax.set_xlim(0.2, max(his) * 1.15)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_orr_forest.png", dpi=160)
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    tarball = download_tarball()
    raw = extract_sample_table(tarball)
    n_libraries = int(len(raw))
    n_patients_raw = int(raw["ANONPT_ID"].nunique())
    d = add_flags(dedup_patients(raw))
    n_patients = int(len(d))
    n_dropped = n_libraries - n_patients

    ev = d[d["orr_evaluable"]].copy()
    bladder = d[d["tissue_bladder"]].copy()
    bladder_ev = ev[ev["tissue_bladder"]].copy()

    coverage = [
        {"item": "RNA libraries in cds", "n": n_libraries, "rule": "CountDataSet columns"},
        {"item": "unique ANONPT_ID", "n": n_patients_raw, "rule": "package patient id"},
        {
            "item": "patients after sizeFactor dedup",
            "n": n_patients,
            "rule": "drop 1 extra library for ANONPT_ID 10285 (both NE)",
        },
        {"item": "libraries dropped as duplicate", "n": n_dropped, "rule": "keep larger sizeFactor"},
        {
            "item": "ORR-evaluable (CR/PR vs SD/PD)",
            "n": int(len(ev)),
            "rule": "package binaryResponse; NE excluded",
        },
        {"item": "CR/PR", "n": int(ev["responder"].sum()), "rule": "binaryResponse"},
        {"item": "SD/PD", "n": int((ev["responder"] == 0).sum()), "rule": "binaryResponse"},
        {"item": "NE / missing binaryResponse", "n": int((~d["orr_evaluable"]).sum()), "rule": "excluded from ORR"},
        {
            "item": "Tissue == bladder",
            "n": int(d["tissue_bladder"].sum()),
            "rule": "package Tissue; other sites are urothelial mets",
        },
        {
            "item": "Tissue == bladder and ORR-evaluable",
            "n": int(len(bladder_ev)),
            "rule": "sensitivity, not primary",
        },
        {
            "item": "Tissue == lung (urothelial met, not NSCLC)",
            "n": int((d["Tissue"] == "lung").sum()),
            "rule": "metastatic deposit; still urothelial analog",
        },
        {
            "item": "CLDN4 / CD274 / CD8A / CXCL9 finite",
            "n": int(d[ [f"{g}_log2TPM1" for g in GENES] ].notna().all(axis=1).sum()),
            "rule": "one row each in cds featureData",
        },
        {"item": "ICI / PD-L1 labels", "n": n_patients, "rule": "single-arm atezolizumab"},
        {"item": "lung-primary NSCLC RNA", "n": 0, "rule": "not this cohort; analog only"},
    ]
    pd.DataFrame(coverage).to_csv(TABLES / "coverage.tsv", sep="\t", index=False)

    tissue_counts = (
        d["Tissue"].fillna("NA").value_counts(dropna=False).rename_axis("Tissue").reset_index(name="n")
    )
    tissue_counts.to_csv(TABLES / "tissue_counts.tsv", sep="\t", index=False)

    spearman_rows = []
    for subset_name, sub in [
        ("all_patients", d),
        ("ORR_evaluable", ev),
        ("bladder_tissue", bladder),
        ("bladder_ORR_evaluable", bladder_ev),
    ]:
        for a, b in [
            ("CLDN4", "CD274"),
            ("CLDN4", "CD8A"),
            ("CLDN4", "CXCL9"),
            ("CD274", "CD8A"),
            ("CD274", "CXCL9"),
            ("CD8A", "CXCL9"),
        ]:
            row = spearman_row(
                sub[f"{a}_log2TPM1"].to_numpy(),
                sub[f"{b}_log2TPM1"].to_numpy(),
                f"{a} vs {b}",
            )
            row["subset"] = subset_name
            row["scale"] = "log2(TPM+1)"
            part = partial_spearman(
                sub[f"{a}_log2TPM1"].to_numpy(),
                sub[f"{b}_log2TPM1"].to_numpy(),
                sub["CD8A_log2TPM1"].to_numpy(),
            )
            if a != "CD8A" and b != "CD8A":
                row["partial_rho_CD8A"] = part["rho"]
                row["partial_p_CD8A"] = part["p"]
            else:
                row["partial_rho_CD8A"] = float("nan")
                row["partial_p_CD8A"] = float("nan")
            # DESeq scale check for primary pair
            if (a, b) == ("CLDN4", "CD274"):
                sf_row = spearman_row(
                    sub["CLDN4_log2SF"].to_numpy(),
                    sub["CD274_log2SF"].to_numpy(),
                    "CLDN4 vs CD274",
                )
                row["rho_log2SF"] = sf_row["rho"]
                row["p_log2SF"] = sf_row["p"]
            spearman_rows.append(row)
    spearman_df = pd.DataFrame(spearman_rows)
    spearman_df.to_csv(TABLES / "spearman.tsv", sep="\t", index=False)

    orr_rows = []
    for subset_name, sub in [
        ("ORR_evaluable", ev),
        ("bladder_ORR_evaluable", bladder_ev),
    ]:
        for gene in GENES:
            orr_rows.append(orr_block(sub, gene, subset_name))
    orr_df = pd.DataFrame(orr_rows)
    orr_df.to_csv(TABLES / "orr.tsv", sep="\t", index=False)

    # one-row headline
    p_all = spearman_df[
        (spearman_df["subset"] == "all_patients") & (spearman_df["pair"] == "CLDN4 vs CD274")
    ].iloc[0]
    p_orr = spearman_df[
        (spearman_df["subset"] == "ORR_evaluable") & (spearman_df["pair"] == "CLDN4 vs CD274")
    ].iloc[0]
    cldn4_orr = orr_df[(orr_df["subset"] == "ORR_evaluable") & (orr_df["gene"] == "CLDN4")].iloc[0]
    cd274_orr = orr_df[(orr_df["subset"] == "ORR_evaluable") & (orr_df["gene"] == "CD274")].iloc[0]
    cxcl9_orr = orr_df[(orr_df["subset"] == "ORR_evaluable") & (orr_df["gene"] == "CXCL9")].iloc[0]
    cd8a_orr = orr_df[(orr_df["subset"] == "ORR_evaluable") & (orr_df["gene"] == "CD8A")].iloc[0]
    bladder_sp = spearman_df[
        (spearman_df["subset"] == "bladder_tissue") & (spearman_df["pair"] == "CLDN4 vs CD274")
    ].iloc[0]
    bladder_cldn4 = orr_df[
        (orr_df["subset"] == "bladder_ORR_evaluable") & (orr_df["gene"] == "CLDN4")
    ].iloc[0]

    one = pd.DataFrame(
        [
            {
                "dataset": "IMvigor210 urothelial analog (not lung)",
                "n_patients": n_patients,
                "n_ORR": int(cldn4_orr["n"]),
                "n_CRPR": int(cldn4_orr["n_CRPR"]),
                "n_SDPD": int(cldn4_orr["n_SDPD"]),
                "CLDN4_CD274_rho_all": p_all["rho"],
                "CLDN4_CD274_p_all": p_all["p"],
                "CLDN4_CD274_rho_ORR": p_orr["rho"],
                "CLDN4_CD274_p_ORR": p_orr["p"],
                "CLDN4_CD274_partial_CD8A_ORR": p_orr["partial_rho_CD8A"],
                "CLDN4_CD274_partial_p_ORR": p_orr["partial_p_CD8A"],
                "CLDN4_ORR_AUC": cldn4_orr["auc"],
                "CLDN4_ORR_mwu_p": cldn4_orr["mwu_p"],
                "CLDN4_ORR_fisher_OR": cldn4_orr["fisher_or"],
                "CLDN4_ORR_fisher_p": cldn4_orr["fisher_p"],
                "CD274_ORR_AUC": cd274_orr["auc"],
                "CD274_ORR_mwu_p": cd274_orr["mwu_p"],
                "CD274_ORR_fisher_OR": cd274_orr["fisher_or"],
                "CD274_ORR_fisher_p": cd274_orr["fisher_p"],
                "CXCL9_ORR_mwu_p": cxcl9_orr["mwu_p"],
                "CD8A_ORR_mwu_p": cd8a_orr["mwu_p"],
                "bladder_n": int(len(bladder)),
                "bladder_ORR_n": int(len(bladder_ev)),
                "bladder_CLDN4_CD274_rho": bladder_sp["rho"],
                "bladder_CLDN4_CD274_p": bladder_sp["p"],
                "bladder_CLDN4_ORR_fisher_OR": bladder_cldn4["fisher_or"],
                "bladder_CLDN4_ORR_fisher_p": bladder_cldn4["fisher_p"],
            }
        ]
    )
    one.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    write_figures(d)

    summary = {
        "slice": "methods/imvigor210_cldn4",
        "cohort": "IMvigor210 metastatic urothelial carcinoma, atezolizumab",
        "analog_not_lung": True,
        "reachable": True,
        "package": "IMvigor210CoreBiologies 1.0.0",
        "sha256": EXPECTED_SHA256,
        "url": PKG_URL,
        "n_libraries": n_libraries,
        "n_patients": n_patients,
        "n_ORR": int(cldn4_orr["n"]),
        "n_CRPR": int(cldn4_orr["n_CRPR"]),
        "n_SDPD": int(cldn4_orr["n_SDPD"]),
        "tissue_counts": tissue_counts.set_index("Tissue")["n"].to_dict(),
        "CLDN4_vs_CD274_all": {
            "n": int(p_all["n"]),
            "rho": float(p_all["rho"]),
            "p": float(p_all["p"]),
            "rho_lo": float(p_all["rho_lo"]),
            "rho_hi": float(p_all["rho_hi"]),
        },
        "CLDN4_vs_CD274_ORR": {
            "n": int(p_orr["n"]),
            "rho": float(p_orr["rho"]),
            "p": float(p_orr["p"]),
            "partial_rho_CD8A": float(p_orr["partial_rho_CD8A"]),
            "partial_p_CD8A": float(p_orr["partial_p_CD8A"]),
        },
        "CLDN4_vs_ORR": {k: (None if (isinstance(v, float) and not np.isfinite(v)) else (int(v) if isinstance(v, (np.integer,)) else (float(v) if isinstance(v, (np.floating, float)) else v))) for k, v in cldn4_orr.to_dict().items()},
        "CD274_vs_ORR": {k: (None if (isinstance(v, float) and not np.isfinite(v)) else (int(v) if isinstance(v, (np.integer,)) else (float(v) if isinstance(v, (np.floating, float)) else v))) for k, v in cd274_orr.to_dict().items()},
        "CXCL9_vs_ORR_mwu_p": float(cxcl9_orr["mwu_p"]),
        "CD8A_vs_ORR_mwu_p": float(cd8a_orr["mwu_p"]),
        "bladder_tissue_n": int(len(bladder)),
        "bladder_ORR_n": int(len(bladder_ev)),
        "bladder_CLDN4_CD274_rho": float(bladder_sp["rho"]),
        "bladder_CLDN4_ORR_fisher_OR": float(bladder_cldn4["fisher_or"]),
        "bladder_CLDN4_ORR_fisher_p": float(bladder_cldn4["fisher_p"]),
        "note": "Urothelial analog. Do not quote as lung / NSCLC.",
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print("Wrote", TABLES)
    print("Wrote", FIGURES)


if __name__ == "__main__":
    main()
