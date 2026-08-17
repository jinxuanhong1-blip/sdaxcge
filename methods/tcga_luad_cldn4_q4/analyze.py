#!/usr/bin/env python3
"""TCGA-LUAD: single-gene CLDN4 Q4 vs Q1 vs ImmuneScore / CD8 / CD274.

Additive extra. TJ7 vs CD8 is already known (PR #69 / #90) and is not
recomputed. Predictor is CLDN4 RNA only. No TACSTD2 gate.

Primary matrix: UCSC Xena HiSeqV2 log2(RSEM norm_count+1), primary tumors
(-01). ImmuneScore = MD Anderson ESTIMATE Immune_score on the matching
RNAseqV2 freeze. Honest n is the HiSeqV2 ∩ ESTIMATE primary intersection
(not the ABSOLUTE-restricted n=502 used elsewhere).
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

GENES = ("CLDN4", "CD8A", "CD8B", "CD274")
ENDPOINTS = (
    ("ImmuneScore", "ESTIMATE Immune_score (MDACC RNAseqV2)"),
    ("CD8A", "CD8A HiSeqV2 log2(norm_count+1)"),
    ("CD8_score", "mean(CD8A, CD8B) — B3-style CD8"),
    ("CD274", "CD274 / PD-L1 HiSeqV2 log2(norm_count+1)"),
)


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
    """Equal-count Q1–Q4. rank-first so ties still fill four bins."""
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
    if n < 4 or d.iloc[:, 0].nunique() < 2 or d.iloc[:, 1].nunique() < 2:
        return {"n": int(n), "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    return {"n": int(n), "rho": float(rho), "p": float(p)}


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


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--hiseqv2",
        default="data/tcga_luad_cldn4_q4/TCGA.LUAD.HiSeqV2.gz",
        type=Path,
    )
    ap.add_argument(
        "--estimate",
        default="data/tcga_luad_cldn4_q4/MDACC_estimate_LUAD_RNAseqV2.txt",
        type=Path,
    )
    ap.add_argument(
        "--slim",
        default="methods/tcga_luad_cldn4_q4/harvested/hiseqv2_cldn4_cd8_cd274.tsv",
        type=Path,
    )
    ap.add_argument("--outdir", default="methods/tcga_luad_cldn4_q4", type=Path)
    args = ap.parse_args()
    out = args.outdir
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "harvested").mkdir(parents=True, exist_ok=True)

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

    core = expr.join(est[["ImmuneScore"]], how="inner")
    core["CD8_score"] = core[["CD8A", "CD8B"]].mean(axis=1)
    core["quartile"] = quartiles(core["CLDN4"])
    core["tacstd2_gate"] = "none"
    core["tj7_used"] = "no"

    q_counts = core["quartile"].value_counts(dropna=False).to_dict()
    cuts = core["CLDN4"]
    q1_max = float(cuts[core["quartile"] == "Q1"].max())
    q4_min = float(cuts[core["quartile"] == "Q4"].min())

    n_rows = [
        {
            "item": "HiSeqV2 primary tumors (-01)",
            "n": n_hiseq_primary,
            "note": "Xena TCGA.LUAD.sampleMap/HiSeqV2; one column per -01",
        },
        {
            "item": "MDACC ESTIMATE primaries (-01)",
            "n": n_est_primary,
            "note": "lung_adenocarcinoma_RNAseqV2.txt Immune_score",
        },
        {
            "item": "HiSeqV2 ∩ ESTIMATE (analysis set)",
            "n": int(len(core)),
            "note": "complete-case for CLDN4 / CD8A / CD8B / CD274 / ImmuneScore",
        },
        {
            "item": "ABSOLUTE not required",
            "n": int(len(core)),
            "note": "do not write n=502; that drop is purity-only",
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
            "item": "Q4 vs Q1 used (all four endpoints)",
            "n": f"{int(q_counts.get('Q4', 0))} vs {int(q_counts.get('Q1', 0))}",
            "note": "same n for ImmuneScore, CD8A, CD8_score, CD274",
        },
        {
            "item": "TJ7 genes in this test",
            "n": 0,
            "note": "single-gene CLDN4 extra; TJ7 vs CD8 already known",
        },
        {
            "item": "TACSTD2 gate",
            "n": 0,
            "note": "none",
        },
    ]
    n_table = pd.DataFrame(n_rows)

    q4_rows = []
    sp_rows = []
    for col, desc in ENDPOINTS:
        d = core[[col, "quartile"]].dropna()
        rec = mannwhitney(d.loc[d["quartile"] == "Q4", col], d.loc[d["quartile"] == "Q1", col])
        rec.update({"endpoint": col, "definition": desc, "n_complete": int(len(d))})
        q4_rows.append(rec)
        sp = spearman(core["CLDN4"], core[col])
        sp.update({"endpoint": col, "definition": desc})
        sp_rows.append(sp)

    q4_df = pd.DataFrame(q4_rows)
    sp_df = pd.DataFrame(sp_rows)

    core.to_csv(out / "results" / "sample_scores.tsv", sep="\t")
    n_table.to_csv(out / "results" / "n_table.tsv", sep="\t", index=False)
    q4_df.to_csv(out / "results" / "q4_vs_q1.tsv", sep="\t", index=False)
    sp_df.to_csv(out / "results" / "spearman.tsv", sep="\t", index=False)

    summary = {
        "cohort": "TCGA-LUAD",
        "predictor": "CLDN4 RNA (single gene)",
        "tj7_recomputed": False,
        "tacstd2_gate": None,
        "n_hiseqv2_primary": n_hiseq_primary,
        "n_estimate_primary": n_est_primary,
        "n_analysis": int(len(core)),
        "n_q1": int(q_counts.get("Q1", 0)),
        "n_q4": int(q_counts.get("Q4", 0)),
        "cldn4_q1_max": q1_max,
        "cldn4_q4_min": q4_min,
        "q4_vs_q1": {r["endpoint"]: r for r in q4_rows},
        "spearman": {r["endpoint"]: r for r in sp_rows},
        "matrix": "UCSC Xena HiSeqV2 log2(RSEM norm_count+1)",
        "immunescore": "MDACC ESTIMATE Immune_score RNAseqV2",
    }
    (out / "results" / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # Three-panel primary figure: ImmuneScore, CD8A, CD274
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.4), constrained_layout=True)
    panel = [("ImmuneScore", "ImmuneScore"), ("CD8A", "CD8A"), ("CD274", "CD274")]
    for ax, (col, title) in zip(axes, panel):
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
        rec = next(r for r in q4_rows if r["endpoint"] == col)
        ax.set_xticks([1, 2], [f"Q1\nn={rec['n_q1']}", f"Q4\nn={rec['n_q4']}"])
        ax.set_title(f"{title}\nΔ={rec['delta_median_q4_minus_q1']:.3g}  p={fmt_p(rec['p'])}")
        ax.set_ylabel(title)
        ax.axhline(d[col].median(), color="0.7", lw=0.6, ls=":")
    fig.suptitle(
        f"TCGA-LUAD CLDN4 Q4 vs Q1  (n={int(len(core))} primaries; Q4 vs Q1 "
        f"{int(q_counts.get('Q4', 0))} vs {int(q_counts.get('Q1', 0))})",
        fontsize=10,
    )
    fig.savefig(out / "figures" / "fig_q4_vs_q1.png", dpi=160)
    plt.close(fig)

    print(n_table.to_string(index=False))
    print()
    print(q4_df.to_string(index=False))
    print()
    print(sp_df.to_string(index=False))


if __name__ == "__main__":
    main()
