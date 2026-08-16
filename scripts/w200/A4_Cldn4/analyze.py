#!/usr/bin/env python3
"""Score TISMO Cldn4 after ICB on the official Tacstd2 64-group universe."""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("results/w200/A4_Cldn4")
DATA = ROOT / "data"
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"
STEM_RE = re.compile(r"\(n=\d+\)$")


def stem(label: str) -> str:
    return STEM_RE.sub("", str(label))


def line_of(group: str) -> str:
    return group.split("_", 1)[0]


def load_export(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    need = {"Samples", "geneID", "value", "cell_line", "Responder", "GSE_ID", "Mouse_treatment"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")
    df = df.copy()
    df["stem"] = df["cell_line"].map(stem)
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    return df


def cancer_map(meta_path: Path) -> dict[str, str]:
    payload = json.loads(meta_path.read_text())
    out: dict[str, str] = {}
    for row in payload.get("data") or []:
        line = row.get("cellLine")
        ctype = row.get("cancerType")
        if line and ctype:
            out[str(line)] = str(ctype)
    return out


def per_group(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group, sub in df.groupby("stem", sort=False):
        base = sub.loc[sub["Responder"] == "Baseline", "value"].dropna()
        treated = sub.loc[sub["Responder"] != "Baseline", "value"].dropna()
        r = sub.loc[sub["Responder"] == "Responders", "value"].dropna()
        nr = sub.loc[sub["Responder"] == "Non-responders", "value"].dropna()
        if len(base) == 0 or len(treated) == 0:
            dmean = np.nan
            dmed = np.nan
        else:
            dmean = float(treated.mean() - base.mean())
            dmed = float(treated.median() - base.median())
        if pd.isna(dmean):
            direction = "NA"
        elif dmean > 0:
            direction = "up"
        elif dmean < 0:
            direction = "down"
        else:
            direction = "tie"
        rows.append(
            {
                "group": group,
                "tismo_label": sub["cell_line"].iloc[0],
                "cell_line": line_of(group),
                "gse_id": sub["GSE_ID"].iloc[0],
                "mouse_treatment": sub["Mouse_treatment"].iloc[0],
                "n_baseline": int(len(base)),
                "n_treated": int(len(treated)),
                "n_responders": int(len(r)),
                "n_nonresponders": int(len(nr)),
                "treated_classes": ",".join(
                    sorted(sub.loc[sub["Responder"] != "Baseline", "Responder"].dropna().unique())
                ),
                "mean_baseline": float(base.mean()) if len(base) else np.nan,
                "mean_treated": float(treated.mean()) if len(treated) else np.nan,
                "median_baseline": float(base.median()) if len(base) else np.nan,
                "median_treated": float(treated.median()) if len(treated) else np.nan,
                "delta_mean": dmean,
                "delta_median": dmed,
                "direction": direction,
            }
        )
    return pd.DataFrame(rows)


def sign_counts(series: pd.Series) -> dict[str, int]:
    return {
        "up": int((series == "up").sum()),
        "down": int((series == "down").sum()),
        "tie": int((series == "tie").sum()),
        "n": int(series.notna().sum()),
    }


def wilcoxon_and_binom(deltas: pd.Series) -> dict:
    d = pd.to_numeric(deltas, errors="coerce").dropna()
    n_up = int((d > 0).sum())
    n_down = int((d < 0).sum())
    n_tie = int((d == 0).sum())
    out = {
        "n": int(len(d)),
        "n_up": n_up,
        "n_down": n_down,
        "n_tie": n_tie,
        "mean_delta": float(d.mean()) if len(d) else None,
        "median_delta": float(d.median()) if len(d) else None,
        "wilcoxon_stat": None,
        "wilcoxon_p": None,
        "binom_k": n_up,
        "binom_n": n_up + n_down,
        "binom_p": None,
    }
    if len(d) >= 2 and (d != 0).any():
        w = stats.wilcoxon(d.to_numpy(), alternative="two-sided", zero_method="wilcox")
        out["wilcoxon_stat"] = float(w.statistic)
        out["wilcoxon_p"] = float(w.pvalue)
    if n_up + n_down >= 1:
        b = stats.binomtest(n_up, n=n_up + n_down, p=0.5, alternative="two-sided")
        out["binom_p"] = float(b.pvalue)
    return out


def mammary_bucket(cancer: str) -> str:
    if cancer.startswith("Mammary"):
        return "Mammary (pooled)"
    return cancer


def write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def waterfall(df: pd.DataFrame, title: str, dest: Path) -> None:
    plot = df.sort_values("delta_mean", ascending=True).reset_index(drop=True)
    colors = [
        "#2ca02c" if v > 0 else ("#7f7f7f" if v == 0 else "#d62728") for v in plot["delta_mean"]
    ]
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.bar(np.arange(len(plot)), plot["delta_mean"], color=colors, width=0.85)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel("Cldn4 mean log(TPM) treated − baseline")
    ax.set_xlabel("TISMO ICB comparison group (same 64 as Tacstd2)")
    ax.set_title(title)
    ax.set_xticks([])
    n_up = int((plot["delta_mean"] > 0).sum())
    n_down = int((plot["delta_mean"] < 0).sum())
    n_tie = int((plot["delta_mean"] == 0).sum())
    ax.text(
        0.02,
        0.95,
        f"{n_up} up / {n_down} down / {n_tie} tie",
        transform=ax.transAxes,
        va="top",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(dest, dpi=150)
    plt.close(fig)


def concordance_scatter(merged: pd.DataFrame, dest: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 5.4))
    ax.axhline(0, color="0.6", lw=0.7)
    ax.axvline(0, color="0.6", lw=0.7)
    ax.scatter(
        merged["delta_mean_tacstd2"],
        merged["delta_mean_cldn4"],
        s=28,
        c="#1f77b4",
        alpha=0.85,
        edgecolors="none",
    )
    ax.set_xlabel("Tacstd2 Δ mean log(TPM)")
    ax.set_ylabel("Cldn4 Δ mean log(TPM)")
    ax.set_title("Same 64 TISMO ICB groups")
    rho, p = stats.spearmanr(merged["delta_mean_tacstd2"], merged["delta_mean_cldn4"])
    ax.text(
        0.04,
        0.96,
        f"Spearman ρ = {rho:.2f}\np = {p:.3g}\nn = {len(merged)}",
        transform=ax.transAxes,
        va="top",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(dest, dpi=150)
    plt.close(fig)


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    cldn4 = load_export(DATA / "cldn4_tismo_vivo.csv")
    tacstd2 = load_export(DATA / "tacstd2_tismo_vivo.csv")
    cmap = cancer_map(DATA / "tismo_vivo_meta.json")

    tac_groups = per_group(tacstd2)
    cld_groups = per_group(cldn4)
    tac_stems = set(tac_groups["group"])
    extra = cld_groups.loc[~cld_groups["group"].isin(tac_stems)].copy()
    cld64 = cld_groups.loc[cld_groups["group"].isin(tac_stems)].copy()
    if len(cld64) != 64 or len(tac_groups) != 64:
        raise RuntimeError(
            f"expected 64/64 groups, got Tacstd2={len(tac_groups)} Cldn4-matched={len(cld64)}"
        )

    cld64["cancer_type"] = cld64["cell_line"].map(cmap)
    cld64["cancer_bucket"] = cld64["cancer_type"].map(mammary_bucket)
    tac_groups["cancer_type"] = tac_groups["cell_line"].map(cmap)

    merged = tac_groups[["group", "direction", "delta_mean"]].merge(
        cld64[["group", "direction", "delta_mean", "cancer_type", "cell_line"]],
        on="group",
        suffixes=("_tacstd2", "_cldn4"),
    )
    same_dir = int((merged["direction_tacstd2"] == merged["direction_cldn4"]).sum())
    spearman = stats.spearmanr(merged["delta_mean_tacstd2"], merged["delta_mean_cldn4"])

    primary = wilcoxon_and_binom(cld64["delta_mean"])
    tac_stats = wilcoxon_and_binom(tac_groups["delta_mean"])
    median_stats = wilcoxon_and_binom(cld64["delta_median"])
    abs01 = wilcoxon_and_binom(cld64.loc[cld64["delta_mean"].abs() >= 0.1, "delta_mean"])
    abs025 = wilcoxon_and_binom(cld64.loc[cld64["delta_mean"].abs() >= 0.25, "delta_mean"])

    line = (
        cld64.groupby("cell_line", as_index=False)
        .agg(
            n_groups=("delta_mean", "size"),
            mean_delta=("delta_mean", "mean"),
            cancer_type=("cancer_type", "first"),
        )
        .sort_values("mean_delta", ascending=False)
    )
    line["direction"] = np.where(
        line["mean_delta"] > 0, "up", np.where(line["mean_delta"] < 0, "down", "tie")
    )
    line_stats = wilcoxon_and_binom(line["mean_delta"])

    cancer_rows = []
    for cancer, sub in cld64.groupby("cancer_bucket"):
        s = wilcoxon_and_binom(sub["delta_mean"])
        s["cancer_bucket"] = cancer
        s["n_groups"] = int(len(sub))
        cancer_rows.append(s)
    cancer_df = pd.DataFrame(cancer_rows).sort_values("n_groups", ascending=False)

    lung = cld64.loc[cld64["cancer_type"] == "Lung carcinoma"].copy()
    extra_stats = wilcoxon_and_binom(cld_groups["delta_mean"])

    both = cld64[(cld64["n_responders"] > 0) & (cld64["n_nonresponders"] > 0)]
    rnr = []
    for group, sub in cldn4.loc[cldn4["stem"].isin(tac_stems)].groupby("stem"):
        r = sub.loc[sub["Responder"] == "Responders", "value"].dropna()
        nr = sub.loc[sub["Responder"] == "Non-responders", "value"].dropna()
        if len(r) and len(nr):
            rnr.append({"group": group, "delta_R_minus_NR": float(r.mean() - nr.mean())})
    rnr_df = pd.DataFrame(rnr)
    rnr_stats = wilcoxon_and_binom(rnr_df["delta_R_minus_NR"]) if len(rnr_df) else {}

    cld64.sort_values("delta_mean", ascending=False).to_csv(
        TABLES / "cldn4_on_tacstd2_64.tsv", sep="\t", index=False
    )
    tac_groups.sort_values("delta_mean", ascending=False).to_csv(
        TABLES / "tacstd2_64.tsv", sep="\t", index=False
    )
    extra.to_csv(TABLES / "cldn4_extra_groups_not_in_64.tsv", sep="\t", index=False)
    line.to_csv(TABLES / "cldn4_by_cell_line.tsv", sep="\t", index=False)
    cancer_df.to_csv(TABLES / "cldn4_by_cancer.tsv", sep="\t", index=False)
    merged.to_csv(TABLES / "cldn4_vs_tacstd2_64.tsv", sep="\t", index=False)
    if len(rnr_df):
        rnr_df.to_csv(TABLES / "cldn4_R_vs_NR.tsv", sep="\t", index=False)

    summary = {
        "task": "A4 analog: TISMO Cldn4 after ICB on the same 64 groups as the Tacstd2 claim",
        "verdict": "NOT_SUPPORTED",
        "verdict_sentence": (
            "Cldn4 is not up after ICB in 49/64 TISMO groups. On the official Tacstd2 "
            "64-group universe, Cldn4 mean treated−baseline is up in 34, down in 28, "
            "and tied in 2 (Wilcoxon p=0.078; sign-test p=0.53)."
        ),
        "universe": {
            "definition": (
                "TISMO gene-module ICB export (All treatments × All tumors, type=3). "
                "One group = cell_line stem with Baseline vs treated (Responders and/or "
                "Non-responders pooled). Same 64 stems as Tacstd2."
            ),
            "n_groups_tacstd2": 64,
            "n_groups_cldn4_export": int(cld_groups["group"].nunique()),
            "n_groups_cldn4_matched": 64,
            "n_cell_lines": int(cld64["cell_line"].nunique()),
            "n_studies": int(cld64["gse_id"].nunique()),
            "extra_cldn4_only_groups": extra["group"].tolist(),
        },
        "tacstd2_lock": {
            "up_down_tie": sign_counts(tac_groups["direction"]),
            "stats": tac_stats,
            "claim_match": (
                tac_stats["n_up"] == 49
                and tac_stats["n_down"] == 15
                and tac_stats["n"] == 64
                and tac_stats["wilcoxon_p"] is not None
                and abs(tac_stats["wilcoxon_p"] - 5.8e-5) < 5e-6
            ),
        },
        "cldn4_primary_same_64": {
            "up_down_tie": sign_counts(cld64["direction"]),
            "stats": primary,
            "matches_tacstd2_49_of_64": False,
        },
        "cldn4_sensitivity": {
            "median_delta": median_stats,
            "abs_delta_ge_0.1": abs01,
            "abs_delta_ge_0.25": abs025,
            "collapse_17_cell_lines": {
                "up_down_tie": sign_counts(line["direction"]),
                "stats": line_stats,
            },
            "all_cldn4_export_groups_including_MOC22": extra_stats,
        },
        "concordance_with_tacstd2": {
            "same_direction": same_dir,
            "n": int(len(merged)),
            "spearman_rho": float(spearman.statistic),
            "spearman_p": float(spearman.pvalue),
        },
        "lung": {
            "n_groups": int(len(lung)),
            "groups": lung["group"].tolist(),
            "up_down_tie": sign_counts(lung["direction"]),
            "note": "Only LLC GSE155972. n=2 cannot carry a 49/64 claim.",
        },
        "response_R_vs_NR": {
            "n_groups_with_both": int(len(both)),
            "stats": rnr_stats,
            "note": "Almost untestable. Do not interpret as response association.",
        },
        "honest_limits": [
            "Not 64 independent models: 17 cell lines, 22 studies; 29/64 groups are mammary.",
            "Not paired longitudinal mice; TISMO pairs Baseline vs treated groups.",
            "Direction is the sign of the mean difference, not DESeq2 significance.",
            "Cldn4 export has one extra group (MOC22_RU31562203_antiPD1, down) excluded to keep the Tacstd2 64.",
            "R vs NR is not the A4 contrast and has only 2 groups with both classes.",
        ],
    }
    write_json(ROOT / "summary.json", summary)
    write_json(
        ROOT / "audit.json",
        {
            "source": {
                "database": "TISMO",
                "url": "https://tismo.pku-genomics.org/",
                "endpoint": "POST /rtismo/gene/downVivoExprn",
                "params": {
                    "filename": "genetreatment_vivo.csv",
                    "type": "3",
                    "icbList": '["All"]',
                    "tumorList": '["All"]',
                    "genes": ["Cldn4", "Tacstd2"],
                },
                "citation": "Zeng et al. Nucleic Acids Research 2022; PMID 34534350",
            },
            "n_rows_cldn4": int(len(cldn4)),
            "n_rows_tacstd2": int(len(tacstd2)),
            "n_samples_cldn4": int(cldn4["Samples"].nunique()),
            "n_samples_tacstd2": int(tacstd2["Samples"].nunique()),
            "sample_id_overlap": int(len(set(cldn4["Samples"]) & set(tacstd2["Samples"]))),
            "files": {
                "cldn4_csv": "data/cldn4_tismo_vivo.csv",
                "tacstd2_csv": "data/tacstd2_tismo_vivo.csv",
                "vivo_meta": "data/tismo_vivo_meta.json",
            },
        },
    )

    pd.DataFrame(
        [
            {
                "analysis": "Tacstd2 lock (official 64)",
                "n_up": tac_stats["n_up"],
                "n_down": tac_stats["n_down"],
                "n_tie": tac_stats["n_tie"],
                "wilcoxon_p": tac_stats["wilcoxon_p"],
                "binom_p": tac_stats["binom_p"],
            },
            {
                "analysis": "Cldn4 on same 64 (primary)",
                "n_up": primary["n_up"],
                "n_down": primary["n_down"],
                "n_tie": primary["n_tie"],
                "wilcoxon_p": primary["wilcoxon_p"],
                "binom_p": primary["binom_p"],
            },
            {
                "analysis": "Cldn4 median-delta on 64",
                "n_up": median_stats["n_up"],
                "n_down": median_stats["n_down"],
                "n_tie": median_stats["n_tie"],
                "wilcoxon_p": median_stats["wilcoxon_p"],
                "binom_p": median_stats["binom_p"],
            },
            {
                "analysis": "Cldn4 |Δ|≥0.1",
                "n_up": abs01["n_up"],
                "n_down": abs01["n_down"],
                "n_tie": abs01["n_tie"],
                "wilcoxon_p": abs01["wilcoxon_p"],
                "binom_p": abs01["binom_p"],
            },
            {
                "analysis": "Cldn4 collapse to 17 cell lines",
                "n_up": line_stats["n_up"],
                "n_down": line_stats["n_down"],
                "n_tie": line_stats["n_tie"],
                "wilcoxon_p": line_stats["wilcoxon_p"],
                "binom_p": line_stats["binom_p"],
            },
        ]
    ).to_csv(TABLES / "summary_stats.tsv", sep="\t", index=False)

    waterfall(
        cld64,
        "Cldn4 after ICB · TISMO official 64 groups (Tacstd2 universe)",
        FIGS / "cldn4_waterfall_64.png",
    )
    concordance_scatter(merged, FIGS / "cldn4_vs_tacstd2_scatter.png")

    print(json.dumps(summary["cldn4_primary_same_64"], indent=2))
    print("verdict:", summary["verdict"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
