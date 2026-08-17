#!/usr/bin/env python3
"""ADDITIVE leftover: TCGA-LUAD CLDN4 vs CXCL9 / CXCL10 / CXCL13 / CXCR3.

Single-gene extra. Predictor is CLDN4 RNA only. No TACSTD2 gate.
Does NOT re-audit Q4 vs CD8 (PR #352) or TJ7 vs CD8 (PR #69 / #90).

Axis: continuous Spearman + CLDN4 Q4 vs Q1 + ESTIMATE residual
(partial Spearman | ESTIMATE_score; OLS residual of each chemokine on
ESTIMATE_score as a matching residual column).

Primary matrix: UCSC Xena HiSeqV2 log2(RSEM norm_count+1), primary -01.
ESTIMATE: official MDACC RNAseqV2 (Yoshihara 2013). We do not recompute it.

These four genes track official ImmuneScore tightly (ρ ≈ 0.60–0.72 here).
Residualising on ImmuneScore is circular for infiltrate-tracking chemokines;
residualising on ESTIMATE_score is still a within-infiltrate leftover
(ImmuneScore is a term in ESTIMATE_score).
"""
from __future__ import annotations

import argparse
import gzip
import json
from math import atanh, tanh
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

GENES = ("CLDN4", "CXCL9", "CXCL10", "CXCL13", "CXCR3")
ENDPOINTS = ("CXCL9", "CXCL10", "CXCL13", "CXCR3")


def sample01(barcode: str) -> str | None:
    b = str(barcode).replace(".", "-")
    parts = b.split("-")
    if len(parts) < 4 or not parts[3].startswith("01"):
        return None
    return "-".join(parts[:3]) + "-01"


def extract_hiseqv2(gz_path: Path, genes: tuple[str, ...]) -> pd.DataFrame:
    with gzip.open(gz_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[1:]
        want = set(genes)
        collected: dict[str, np.ndarray] = {}
        for line in f:
            g = line.split("\t", 1)[0]
            if g not in want:
                continue
            vals = np.array(line.rstrip("\n").split("\t")[1:], dtype=float)
            collected[g] = vals
            if len(collected) == len(want):
                break
    missing = [g for g in genes if g not in collected]
    if missing:
        raise SystemExit(f"HiSeqV2 missing genes: {missing}")
    mat = pd.DataFrame(collected, index=samples)
    sid = [sample01(s) for s in mat.index]
    mat = mat.loc[[s is not None for s in sid]].copy()
    mat.index = [sample01(s) for s in mat.index]
    mat = mat.groupby(level=0).mean()
    mat.index.name = "sample"
    return mat


def load_estimate(path: Path) -> pd.DataFrame:
    est = pd.read_csv(path, sep="\t")
    est["sample"] = est["ID"].map(sample01)
    est = est.dropna(subset=["sample"]).groupby("sample").mean(numeric_only=True)
    return est.rename(
        columns={
            "Immune_score": "ImmuneScore",
            "Stromal_score": "StromalScore",
            "ESTIMATE_score": "ESTIMATE_score",
        }
    )


def quartiles(s: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=s.index, dtype=object)
    ok = s.dropna()
    if len(ok) < 8:
        return out
    labels = pd.qcut(ok.rank(method="first"), 4, labels=["Q1", "Q2", "Q3", "Q4"])
    out.loc[ok.index] = labels.astype(str)
    return out


def spearman(x: pd.Series, y: pd.Series) -> dict:
    d = pd.concat([x, y], axis=1).dropna()
    n = len(d)
    rec = {"n": int(n), "rho": np.nan, "p": np.nan, "ci_lo": np.nan, "ci_hi": np.nan}
    if n < 4 or d.iloc[:, 0].nunique() < 2 or d.iloc[:, 1].nunique() < 2:
        return rec
    rho, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    rec["rho"] = float(rho)
    rec["p"] = float(p)
    z = atanh(np.clip(rho, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(max(n - 3, 1))
    rec["ci_lo"] = float(tanh(z - 1.96 * se))
    rec["ci_hi"] = float(tanh(z + 1.96 * se))
    return rec


def partial_spearman(x: pd.Series, y: pd.Series, z: pd.Series) -> dict:
    d = pd.concat([x, y, z], axis=1).dropna()
    n = len(d)
    rec = {"n": int(n), "rho": np.nan, "p": np.nan}
    if n < 6:
        return rec
    rxy = float(stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1]).statistic)
    rxz = float(stats.spearmanr(d.iloc[:, 0], d.iloc[:, 2]).statistic)
    ryz = float(stats.spearmanr(d.iloc[:, 1], d.iloc[:, 2]).statistic)
    denom = np.sqrt(max((1 - rxz**2) * (1 - ryz**2), 1e-12))
    r = float(np.clip((rxy - rxz * ryz) / denom, -0.999999, 0.999999))
    df = n - 3
    t = r * np.sqrt(df / max(1 - r**2, 1e-12))
    p = float(2 * stats.t.sf(abs(t), df))
    rec.update({"rho": r, "p": p})
    return rec


def ols_residual(y: pd.Series, x: pd.Series) -> pd.Series:
    d = pd.concat([y, x], axis=1).dropna()
    out = pd.Series(np.nan, index=y.index, dtype=float)
    if len(d) < 4 or d.iloc[:, 1].nunique() < 2:
        return out
    xv = d.iloc[:, 1].to_numpy(float)
    yv = d.iloc[:, 0].to_numpy(float)
    slope, intercept = np.polyfit(xv, yv, 1)
    out.loc[d.index] = yv - (intercept + slope * xv)
    return out


def mannwhitney(q4: pd.Series, q1: pd.Series) -> dict:
    a = pd.to_numeric(q4, errors="coerce").dropna()
    b = pd.to_numeric(q1, errors="coerce").dropna()
    rec = {
        "n_q4": int(len(a)),
        "n_q1": int(len(b)),
        "median_q4": float(a.median()) if len(a) else np.nan,
        "median_q1": float(b.median()) if len(b) else np.nan,
        "delta_median_q4_minus_q1": np.nan,
        "U": np.nan,
        "p": np.nan,
        "rank_biserial_q4_gt_q1": np.nan,
    }
    if len(a) and len(b):
        rec["delta_median_q4_minus_q1"] = float(a.median() - b.median())
    if len(a) < 3 or len(b) < 3:
        return rec
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    r = (2.0 * U) / (len(a) * len(b)) - 1.0
    rec.update({"U": float(U), "p": float(p), "rank_biserial_q4_gt_q1": float(r)})
    return rec


def bh(pvals: list[float]) -> list[float]:
    n = len(pvals)
    order = np.argsort(pvals)
    q = np.empty(n, dtype=float)
    prev = 1.0
    for rank, i in enumerate(order[::-1], start=1):
        adj = pvals[i] * n / (n - rank + 1)
        prev = min(prev, adj)
        q[i] = prev
    return [float(min(v, 1.0)) for v in q]


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_rho_p(rho: float, p: float) -> str:
    if not np.isfinite(rho):
        return "NA"
    sign = "+" if rho >= 0 else ""
    return f"{sign}{rho:.3f} ({fmt_p(p)})"


def label_rho(rho: float, q: float) -> str:
    if not np.isfinite(rho) or not np.isfinite(q):
        return "NA"
    if q >= 0.05:
        return "NULL"
    mag = abs(rho)
    direction = "POSITIVE" if rho > 0 else "NEGATIVE"
    if mag < 0.15:
        return f"WEAK_{direction}"
    if mag < 0.30:
        return f"SMALL_{direction}"
    return f"ASSOCIATED_{direction}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--hiseqv2",
        default="data/tcga_luad_cldn4_cxcl/TCGA.LUAD.HiSeqV2.gz",
        type=Path,
    )
    ap.add_argument(
        "--estimate",
        default="data/tcga_luad_cldn4_cxcl/MDACC_estimate_LUAD_RNAseqV2.txt",
        type=Path,
    )
    ap.add_argument(
        "--slim",
        default="methods/tcga_luad_cldn4_cxcl/harvested/hiseqv2_cldn4_cxcl.tsv",
        type=Path,
    )
    ap.add_argument("--outdir", default="methods/tcga_luad_cldn4_cxcl", type=Path)
    args = ap.parse_args()
    out = args.outdir
    for sub in ("tables", "figures", "harvested", "results"):
        (out / sub).mkdir(parents=True, exist_ok=True)

    if args.slim.exists():
        expr = pd.read_csv(args.slim, sep="\t", index_col=0)
    elif args.hiseqv2.exists():
        expr = extract_hiseqv2(args.hiseqv2, GENES)
        args.slim.parent.mkdir(parents=True, exist_ok=True)
        expr.to_csv(args.slim, sep="\t")
    else:
        raise SystemExit("Need HiSeqV2.gz or the slim harvested gene matrix")

    n_hiseq_primary = int(len(expr))
    est = load_estimate(args.estimate)
    n_est_primary = int(len(est))

    core = expr.join(est[["ImmuneScore", "StromalScore", "ESTIMATE_score"]], how="inner")
    for g in GENES:
        if core[g].isna().any():
            raise SystemExit(f"NA values in {g}")
    core["quartile"] = quartiles(core["CLDN4"])
    for col in ENDPOINTS:
        core[f"{col}_resid_ESTIMATE"] = ols_residual(core[col], core["ESTIMATE_score"])
        core[f"{col}_resid_Immune"] = ols_residual(core[col], core["ImmuneScore"])
    core["CLDN4_resid_ESTIMATE"] = ols_residual(core["CLDN4"], core["ESTIMATE_score"])

    q_counts = core["quartile"].value_counts(dropna=False).to_dict()
    q1_max = float(core.loc[core["quartile"] == "Q1", "CLDN4"].max())
    q4_min = float(core.loc[core["quartile"] == "Q4", "CLDN4"].min())

    presence_rows = []
    for g in GENES:
        presence_rows.append(
            {
                "gene": g,
                "on_HiSeqV2": "yes",
                "n_finite": int(core[g].notna().sum()),
                "role": "predictor" if g == "CLDN4" else "endpoint",
            }
        )
    presence = pd.DataFrame(presence_rows)

    n_rows = [
        {
            "item": "HiSeqV2 primary tumors (-01)",
            "n": n_hiseq_primary,
            "note": "Xena TCGA.LUAD.sampleMap/HiSeqV2; one column per -01",
        },
        {
            "item": "MDACC ESTIMATE primaries (-01)",
            "n": n_est_primary,
            "note": "lung_adenocarcinoma_RNAseqV2.txt",
        },
        {
            "item": "HiSeqV2 ∩ ESTIMATE (analysis set)",
            "n": int(len(core)),
            "note": "complete-case CLDN4 / CXCL9 / CXCL10 / CXCL13 / CXCR3 / ESTIMATE",
        },
        {
            "item": "ABSOLUTE not required",
            "n": int(len(core)),
            "note": "do not write n=502; that drop is purity-only (PR #188 / #311)",
        },
        {
            "item": "Q1 / Q2 / Q3 / Q4",
            "n": "{}/{}/{}/{}".format(
                int(q_counts.get("Q1", 0)),
                int(q_counts.get("Q2", 0)),
                int(q_counts.get("Q3", 0)),
                int(q_counts.get("Q4", 0)),
            ),
            "note": "pd.qcut(rank(method='first'), 4) on CLDN4",
        },
        {
            "item": "Q4 vs Q1 used (all four CXCL endpoints)",
            "n": f"{int(q_counts.get('Q4', 0))} vs {int(q_counts.get('Q1', 0))}",
            "note": "same n for CXCL9 / CXCL10 / CXCL13 / CXCR3",
        },
        {
            "item": "endpoints tracking ImmuneScore (ρ>0.5)",
            "n": 4,
            "note": "CXCL9/10/13/CXCR3 vs ImmuneScore ρ=0.60–0.72; ESTIMATE residual is within-infiltrate leftover",
        },
        {
            "item": "TJ7 genes in this test",
            "n": 0,
            "note": "not recomputed; TJ7 vs CD8 already known (PR #69 / #90)",
        },
        {
            "item": "CD8 / CD274 re-audit",
            "n": 0,
            "note": "out of scope; see PR #352",
        },
        {
            "item": "TACSTD2 gate",
            "n": 0,
            "note": "none",
        },
    ]
    n_table = pd.DataFrame(n_rows)

    sp_rows = []
    q4_rows = []
    resid_rows = []
    for col in ENDPOINTS:
        sp = spearman(core["CLDN4"], core[col])
        sp.update({"endpoint": col, "test": "unadjusted_spearman"})
        sp_rows.append(sp)

        rec = mannwhitney(
            core.loc[core["quartile"] == "Q4", col],
            core.loc[core["quartile"] == "Q1", col],
        )
        rec.update({"endpoint": col, "scale": "raw_HiSeqV2"})
        q4_rows.append(rec)

        part_est = partial_spearman(core["CLDN4"], core[col], core["ESTIMATE_score"])
        part_imm = partial_spearman(core["CLDN4"], core[col], core["ImmuneScore"])
        resid_sp = spearman(core["CLDN4"], core[f"{col}_resid_ESTIMATE"])
        both_resid = spearman(core["CLDN4_resid_ESTIMATE"], core[f"{col}_resid_ESTIMATE"])
        resid_rows.append(
            {
                "endpoint": col,
                "n": part_est["n"],
                "partial_rho_ESTIMATE": part_est["rho"],
                "partial_p_ESTIMATE": part_est["p"],
                "partial_rho_ImmuneScore": part_imm["rho"],
                "partial_p_ImmuneScore": part_imm["p"],
                "ols_resid_rho_CLDN4_vs_Yresid": resid_sp["rho"],
                "ols_resid_p_CLDN4_vs_Yresid": resid_sp["p"],
                "ols_both_resid_rho": both_resid["rho"],
                "ols_both_resid_p": both_resid["p"],
                "circularity": "endpoint tracks ImmuneScore; ESTIMATE residual is within-infiltrate leftover",
            }
        )

        rec_r = mannwhitney(
            core.loc[core["quartile"] == "Q4", f"{col}_resid_ESTIMATE"],
            core.loc[core["quartile"] == "Q1", f"{col}_resid_ESTIMATE"],
        )
        rec_r.update({"endpoint": col, "scale": "OLS_residual_on_ESTIMATE_score"})
        q4_rows.append(rec_r)

    q_unadj = bh([r["p"] for r in sp_rows])
    q_part = bh([r["partial_p_ESTIMATE"] for r in resid_rows])
    for r, q in zip(sp_rows, q_unadj):
        r["bh_q"] = q
        r["label"] = label_rho(r["rho"], q)
    for r, q in zip(resid_rows, q_part):
        r["partial_bh_q_ESTIMATE"] = q
        r["partial_label"] = label_rho(r["partial_rho_ESTIMATE"], q)

    # CLDN4 vs ESTIMATE scores (context, not a chemokine claim)
    context = [
        {"pair": "CLDN4 vs ImmuneScore", **spearman(core["CLDN4"], core["ImmuneScore"])},
        {"pair": "CLDN4 vs ESTIMATE_score", **spearman(core["CLDN4"], core["ESTIMATE_score"])},
        {"pair": "CLDN4 vs StromalScore", **spearman(core["CLDN4"], core["StromalScore"])},
    ]
    for col in ENDPOINTS:
        context.append(
            {
                "pair": f"{col} vs ESTIMATE_score",
                **spearman(core[col], core["ESTIMATE_score"]),
            }
        )
        context.append(
            {
                "pair": f"{col} vs ImmuneScore",
                **spearman(core[col], core["ImmuneScore"]),
            }
        )

    sp_df = pd.DataFrame(sp_rows)
    q4_df = pd.DataFrame(q4_rows)
    resid_df = pd.DataFrame(resid_rows)
    context_df = pd.DataFrame(context)

    one = {
        "dataset": "TCGA-LUAD HiSeqV2 ∩ MDACC ESTIMATE",
        "n": int(len(core)),
        "n_Q4_vs_Q1": f"{int(q_counts.get('Q4', 0))} vs {int(q_counts.get('Q1', 0))}",
        "n_rule": "HiSeqV2 ∩ official ESTIMATE primaries (-01); complete pairs; not ABSOLUTE n=502",
        "CLDN4_CXCL9_rho": sp_rows[0]["rho"],
        "CLDN4_CXCL9_p": sp_rows[0]["p"],
        "CLDN4_CXCL9_q": sp_rows[0]["bh_q"],
        "CLDN4_CXCL10_rho": sp_rows[1]["rho"],
        "CLDN4_CXCL10_p": sp_rows[1]["p"],
        "CLDN4_CXCL10_q": sp_rows[1]["bh_q"],
        "CLDN4_CXCL13_rho": sp_rows[2]["rho"],
        "CLDN4_CXCL13_p": sp_rows[2]["p"],
        "CLDN4_CXCL13_q": sp_rows[2]["bh_q"],
        "CLDN4_CXCR3_rho": sp_rows[3]["rho"],
        "CLDN4_CXCR3_p": sp_rows[3]["p"],
        "CLDN4_CXCR3_q": sp_rows[3]["bh_q"],
        "CXCL9_ESTIMATE_partial_rho": resid_rows[0]["partial_rho_ESTIMATE"],
        "CXCL9_ESTIMATE_partial_p": resid_rows[0]["partial_p_ESTIMATE"],
        "CXCL9_ESTIMATE_partial_q": resid_rows[0]["partial_bh_q_ESTIMATE"],
        "CXCL10_ESTIMATE_partial_rho": resid_rows[1]["partial_rho_ESTIMATE"],
        "CXCL10_ESTIMATE_partial_p": resid_rows[1]["partial_p_ESTIMATE"],
        "CXCL10_ESTIMATE_partial_q": resid_rows[1]["partial_bh_q_ESTIMATE"],
        "CXCL13_ESTIMATE_partial_rho": resid_rows[2]["partial_rho_ESTIMATE"],
        "CXCL13_ESTIMATE_partial_p": resid_rows[2]["partial_p_ESTIMATE"],
        "CXCL13_ESTIMATE_partial_q": resid_rows[2]["partial_bh_q_ESTIMATE"],
        "CXCR3_ESTIMATE_partial_rho": resid_rows[3]["partial_rho_ESTIMATE"],
        "CXCR3_ESTIMATE_partial_p": resid_rows[3]["partial_p_ESTIMATE"],
        "CXCR3_ESTIMATE_partial_q": resid_rows[3]["partial_bh_q_ESTIMATE"],
        "CD8_reaudit": "no (PR #352 out of scope)",
        "TJ7_reaudit": "no (PR #69 / #90 out of scope)",
        "TACSTD2_gate": "none",
        "ICI_labels": "none (treatment-naive surgical RNA)",
    }
    one_df = pd.DataFrame([one])

    core.to_csv(out / "tables" / "sample_scores.tsv", sep="\t")
    core.to_csv(out / "results" / "sample_scores.tsv", sep="\t")
    n_table.to_csv(out / "tables" / "n_table.tsv", sep="\t", index=False)
    n_table.to_csv(out / "results" / "n_table.tsv", sep="\t", index=False)
    presence.to_csv(out / "tables" / "gene_presence.tsv", sep="\t", index=False)
    sp_df.to_csv(out / "tables" / "spearman.tsv", sep="\t", index=False)
    q4_df.to_csv(out / "tables" / "q4_vs_q1.tsv", sep="\t", index=False)
    resid_df.to_csv(out / "tables" / "estimate_residual.tsv", sep="\t", index=False)
    context_df.to_csv(out / "tables" / "estimate_context.tsv", sep="\t", index=False)
    one_df.to_csv(out / "tables" / "one_row.tsv", sep="\t", index=False)

    summary = {
        "cohort": "TCGA-LUAD",
        "predictor": "CLDN4 RNA (single gene)",
        "axis": "CXCL9 / CXCL10 / CXCL13 / CXCR3",
        "q4_vs_cd8_reaudit": False,
        "tj7_vs_cd8_reaudit": False,
        "tacstd2_gate": None,
        "n_hiseqv2_primary": n_hiseq_primary,
        "n_estimate_primary": n_est_primary,
        "n_analysis": int(len(core)),
        "n_q1": int(q_counts.get("Q1", 0)),
        "n_q4": int(q_counts.get("Q4", 0)),
        "cldn4_q1_max": q1_max,
        "cldn4_q4_min": q4_min,
        "infiltrate_tracking_endpoints": list(ENDPOINTS),
        "spearman": {r["endpoint"]: r for r in sp_rows},
        "q4_vs_q1": {r["endpoint"] + "_" + r["scale"]: r for r in q4_rows},
        "estimate_residual": {r["endpoint"]: r for r in resid_rows},
        "one_row": one,
        "matrix": "UCSC Xena HiSeqV2 log2(RSEM norm_count+1)",
        "estimate": "MDACC ESTIMATE RNAseqV2 (ImmuneScore / StromalScore / ESTIMATE_score)",
    }
    (out / "tables" / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out / "results" / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # Extra figures
    fig, ax = plt.subplots(figsize=(7.2, 3.8), constrained_layout=True)
    ys = np.arange(len(ENDPOINTS))
    rho_u = [r["rho"] for r in sp_rows]
    lo_u = [r["ci_lo"] for r in sp_rows]
    hi_u = [r["ci_hi"] for r in sp_rows]
    rho_r = [r["partial_rho_ESTIMATE"] for r in resid_rows]
    ax.axvline(0, color="0.6", lw=0.8)
    ax.errorbar(
        rho_u,
        ys + 0.12,
        xerr=[np.array(rho_u) - np.array(lo_u), np.array(hi_u) - np.array(rho_u)],
        fmt="o",
        color="#2c5f8a",
        label="unadjusted Spearman",
        capsize=3,
    )
    ax.plot(rho_r, ys - 0.12, "s", color="#b85c38", label="partial | ESTIMATE_score")
    ax.set_yticks(ys, ENDPOINTS)
    ax.set_xlabel("Spearman ρ vs CLDN4")
    ax.set_title(
        f"TCGA-LUAD CLDN4 vs CXCL axis  (n={len(core)}; Q4 vs Q1 "
        f"{q_counts.get('Q4', 0)} vs {q_counts.get('Q1', 0)})"
    )
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    ax.set_xlim(-0.45, 0.35)
    fig.savefig(out / "figures" / "fig1_spearman_forest.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.2), constrained_layout=True)
    for ax, col, rec in zip(axes.ravel(), ENDPOINTS, sp_rows):
        ax.scatter(core["CLDN4"], core[col], s=8, alpha=0.35, c="#3d6b8a", linewidths=0)
        ax.set_xlabel("CLDN4")
        ax.set_ylabel(col)
        ax.set_title(f"{col}  ρ={rec['rho']:+.3f}  p={fmt_p(rec['p'])}  n={rec['n']}")
    fig.suptitle("TCGA-LUAD CLDN4 vs CXCL9 / CXCL10 / CXCL13 / CXCR3 (HiSeqV2 primaries)", fontsize=10)
    fig.savefig(out / "figures" / "fig2_scatter.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 4, figsize=(11.2, 3.4), constrained_layout=True)
    for ax, col in zip(axes, ENDPOINTS):
        d = core[[col, "quartile"]].dropna()
        data = [d.loc[d["quartile"] == q, col].to_numpy() for q in ("Q1", "Q4")]
        parts = ax.violinplot(data, positions=[1, 2], showextrema=False, widths=0.7)
        for pc in parts["bodies"]:
            pc.set_facecolor("#6b8cae")
            pc.set_alpha(0.55)
        ax.boxplot(
            data,
            positions=[1, 2],
            widths=0.18,
            showfliers=False,
            medianprops={"color": "black"},
        )
        rec = next(r for r in q4_rows if r["endpoint"] == col and r["scale"] == "raw_HiSeqV2")
        ax.set_xticks([1, 2], [f"Q1\nn={rec['n_q1']}", f"Q4\nn={rec['n_q4']}"])
        ax.set_title(f"{col}\nΔ={rec['delta_median_q4_minus_q1']:.3g}  p={fmt_p(rec['p'])}")
        ax.set_ylabel(col)
    fig.suptitle(
        f"CLDN4 Q4 vs Q1 CXCL axis  (analysis n={len(core)}; MWU "
        f"{q_counts.get('Q4', 0)} vs {q_counts.get('Q1', 0)})",
        fontsize=10,
    )
    fig.savefig(out / "figures" / "fig3_q4_vs_q1.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.2), constrained_layout=True)
    for ax, col, rec in zip(axes.ravel(), ENDPOINTS, resid_rows):
        ax.scatter(
            core["CLDN4"],
            core[f"{col}_resid_ESTIMATE"],
            s=8,
            alpha=0.35,
            c="#8a4b2c",
            linewidths=0,
        )
        ax.axhline(0, color="0.7", lw=0.6, ls=":")
        ax.set_xlabel("CLDN4")
        ax.set_ylabel(f"{col} residual | ESTIMATE_score")
        ax.set_title(
            f"{col} residual  partial ρ={rec['partial_rho_ESTIMATE']:+.3f}  "
            f"p={fmt_p(rec['partial_p_ESTIMATE'])}"
        )
    fig.suptitle(
        "ESTIMATE residual leftover (OLS chemokine ~ ESTIMATE_score; infiltrate-tracking genes)",
        fontsize=10,
    )
    fig.savefig(out / "figures" / "fig4_estimate_residual.png", dpi=160)
    plt.close(fig)

    print(n_table.to_string(index=False))
    print()
    print(sp_df.to_string(index=False))
    print()
    print(q4_df.to_string(index=False))
    print()
    print(resid_df.to_string(index=False))
    print()
    print(one_df.to_string(index=False))


if __name__ == "__main__":
    main()
