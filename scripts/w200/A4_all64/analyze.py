#!/usr/bin/env python3
"""Honest recompute of A4: Tacstd2 after ICB across TISMO's 64 slices."""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon

ROOT = Path("results/w200/A4_all64")
SRC = ROOT / "tacstd2_vivo_icb.csv"
FIG = ROOT / "figures"

# Official TISMO cellLineMeta cancerType, plus a collapsed group for splits.
CANCER = {
    "4T1": ("Mammary carcinoma", "Mammary"),
    "E0771": ("Mammary adenocarcinoma", "Mammary"),
    "EMT6": ("Mammary carcinoma", "Mammary"),
    "KPB25L": ("Mammary cancer, NOS", "Mammary"),
    "T11": ("Mammary cancer, NOS", "Mammary"),
    "p53-2225L": ("Mammary cancer, NOS", "Mammary"),
    "p53-2336R": ("Mammary cancer, NOS", "Mammary"),
    "B16": ("Melanoma", "Melanoma"),
    "YUMM1.7": ("Melanoma", "Melanoma"),
    "D3UV2": ("Melanoma", "Melanoma"),
    "D4M.3A.3": ("Melanoma", "Melanoma"),
    "CT26": ("Colorectal carcinoma", "Colorectal"),
    "MC38": ("Colorectal carcinoma", "Colorectal"),
    "LLC": ("Lung carcinoma", "Lung"),
    "402230": ("Sarcoma", "Sarcoma"),
    "BNL-MEA": ("Hepatocellular carcinoma", "Liver"),
    "YTN16": ("Gastric adenocarcinoma", "Gastric"),
}

# Tokens that are ICB-agent labels, not a separate biological model.
ICB_TOKEN_RE = re.compile(
    r"antiCTLA4&antiPD[L]?1|antiCTLA4|antiPD1_FcS|antiPD1|antiPDL1|antiPDL2",
    re.I,
)
TIME_RE = re.compile(r"(?:^|_)(day\d+|end)(?:_|$)")
N_RE = re.compile(r"\(n=\d+\)$")


def sign_counts(values: pd.Series | np.ndarray) -> dict:
    v = pd.Series(values, dtype=float).dropna()
    n_up = int((v > 0).sum())
    n_down = int((v < 0).sum())
    n_tie = int((v == 0).sum())
    n = n_up + n_down
    if n == 0:
        p_two = p_greater = np.nan
    else:
        p_two = float(binomtest(n_up, n, 0.5, alternative="two-sided").pvalue)
        p_greater = float(binomtest(n_up, n, 0.5, alternative="greater").pvalue)
    return {
        "n_models": int(len(v)),
        "n_up": n_up,
        "n_down": n_down,
        "n_tie": n_tie,
        "n_tested": n,
        "sign_p_two_sided": p_two,
        "sign_p_greater": p_greater,
    }


def wilcoxon_delta(values: pd.Series, zero_method: str = "wilcox") -> dict:
    v = pd.Series(values, dtype=float).dropna()
    if zero_method == "wilcox":
        v = v[v != 0]
    if len(v) < 1:
        return {"n": int(len(v)), "statistic": np.nan, "p_two_sided": np.nan}
    res = wilcoxon(v, alternative="two-sided", zero_method=zero_method)
    return {
        "n": int(len(v)),
        "statistic": float(res.statistic),
        "p_two_sided": float(res.pvalue),
    }


def parse_line(model: str) -> str:
    return model.split("_")[0]


def parse_gse(model: str) -> str:
    m = re.search(r"(GSE\d+|ERP\d+)", model)
    return m.group(1) if m else "NA"


def strip_time(model: str) -> str:
    return TIME_RE.sub("_", N_RE.sub("", model)).strip("_")


def time_rank(model: str) -> float | None:
    m = re.search(r"_(day(\d+)|end)_", model)
    if not m:
        return None
    if m.group(1) == "end":
        return 999.0
    return float(m.group(2))


def independent_key(model: str) -> str:
    """Collapse ICB agent + timepoint; keep genotype / co-treatment / site."""
    rest = N_RE.sub("", model)
    rest = TIME_RE.sub("_", rest)
    rest = ICB_TOKEN_RE.sub("", rest)
    rest = re.sub(r"_+", "_", rest).strip("_")
    return rest or parse_line(model)


def per_model_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for model, g in df.groupby("cell_line", sort=True):
        line = parse_line(model)
        cancer, group = CANCER[line]
        base = g.loc[g["Baseline"] == 1, "value"].astype(float)
        trt = g.loc[g["Baseline"] == 0, "value"].astype(float)
        resp = g.loc[g["Responder"] == "Responders", "value"].astype(float)
        nr = g.loc[g["Responder"] == "Non-responders", "value"].astype(float)
        d_med = float(trt.median() - base.median())
        d_mean = float(trt.mean() - base.mean())
        rows.append(
            {
                "model": model,
                "cell_line": line,
                "gse_id": g["GSE_ID"].iloc[0],
                "cancer_type": cancer,
                "cancer_group": group,
                "independent_key": independent_key(model),
                "n_baseline": int(len(base)),
                "n_icb": int(len(trt)),
                "n_responder": int(len(resp)),
                "n_nonresponder": int(len(nr)),
                "median_baseline": float(base.median()),
                "median_icb": float(trt.median()),
                "mean_baseline": float(base.mean()),
                "mean_icb": float(trt.mean()),
                "delta_median": d_med,
                "delta_mean": d_mean,
                "dir_median": int(np.sign(d_med)),
                "dir_mean": int(np.sign(d_mean)) if d_mean != 0 else 0,
                "treatments": ";".join(sorted(g["Mouse_treatment"].dropna().unique())),
                "tismo_pvalue": (
                    float(g["pvalue"].iloc[0]) if pd.notna(g["pvalue"].iloc[0]) else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def last_timepoint(models: pd.DataFrame) -> pd.DataFrame:
    keep = []
    for _, sub in models.groupby(models["model"].map(strip_time)):
        ranks = sub["model"].map(time_rank)
        if ranks.notna().any():
            keep.extend(sub.loc[ranks == ranks.max(), "model"])
        else:
            keep.extend(sub["model"])
    return models[models["model"].isin(keep)].copy()


def collapse_unique_samples(df: pd.DataFrame, models: pd.DataFrame, key_col: str) -> pd.DataFrame:
    """One pair per key using unique SRX IDs (shared baselines counted once)."""
    rows = []
    lookup = models.set_index("model")
    for key, sub in models.groupby(key_col, sort=True):
        names = set(sub["model"])
        g = df[df["cell_line"].isin(names)].drop_duplicates("Samples")
        base = g.loc[g["Baseline"] == 1, "value"].astype(float)
        trt = g.loc[g["Baseline"] == 0, "value"].astype(float)
        if len(base) == 0 or len(trt) == 0:
            continue
        line = sub["cell_line"].iloc[0]
        rows.append(
            {
                "unit": key,
                "cell_line": line,
                "cancer_group": lookup.loc[sub["model"].iloc[0], "cancer_group"],
                "n_slices": int(len(sub)),
                "n_baseline": int(len(base)),
                "n_icb": int(len(trt)),
                "delta_median": float(trt.median() - base.median()),
                "delta_mean": float(trt.mean() - base.mean()),
            }
        )
    return pd.DataFrame(rows)


def row_from_counts(label: str, counts: dict, extra: dict | None = None) -> dict:
    out = {"analysis": label, **counts}
    if extra:
        out.update(extra)
    return out


def _json_clean(obj):
    if isinstance(obj, dict):
        return {k: _json_clean(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_json_clean(v) for v in obj]
    if isinstance(obj, float) and (np.isnan(obj) or np.isinf(obj)):
        return None
    return obj


def write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(_json_clean(obj), indent=2) + "\n")


def plot_waterfall(models: pd.DataFrame, path: Path) -> None:
    d = models.sort_values(["delta_median", "model"]).reset_index(drop=True)
    colors = {
        "Mammary": "#4C78A8",
        "Melanoma": "#F58518",
        "Colorectal": "#54A24B",
        "Gastric": "#E45756",
        "Lung": "#72B7B2",
        "Sarcoma": "#B279A2",
        "Liver": "#FF9DA6",
    }
    fig, ax = plt.subplots(figsize=(11, 4.2))
    ax.axhline(0, color="0.4", lw=0.8)
    ax.bar(
        np.arange(len(d)),
        d["delta_median"],
        color=[colors[g] for g in d["cancer_group"]],
        width=0.9,
    )
    ax.set_ylabel("median Tacstd2 (ICB − baseline)")
    ax.set_xlabel("TISMO ICB slice (n = 64), sorted")
    ax.set_xticks([])
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=c, label=k) for k, c in colors.items() if k in set(d["cancer_group"])
    ]
    ax.legend(handles=handles, frameon=False, ncol=4, fontsize=8, loc="upper left")
    ax.set_title("Tacstd2 after ICB: per-slice median difference")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_cancer_split(split: pd.DataFrame, path: Path) -> None:
    sub = split[split["location"] == "median"].copy()
    order = ["Mammary", "Melanoma", "Colorectal", "Gastric", "Lung", "Sarcoma", "Liver"]
    sub = sub.set_index("cancer_group").loc[[g for g in order if g in set(sub["cancer_group"])]]
    x = np.arange(len(sub))
    fig, ax = plt.subplots(figsize=(8.5, 4.0))
    ax.bar(x - 0.2, sub["n_up"], 0.4, label="up", color="#4C78A8")
    ax.bar(x + 0.2, sub["n_down"], 0.4, label="down", color="#E45756")
    ax.set_xticks(x)
    ax.set_xticklabels(sub.index, rotation=20, ha="right")
    ax.set_ylabel("TISMO slices (ties omitted)")
    ax.legend(frameon=False)
    ax.set_title("Cancer-type split of median ICB vs baseline sign")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    if not SRC.exists():
        raise SystemExit(f"missing {SRC}; run download.py first")

    raw = pd.read_csv(SRC)
    expected = {
        "Samples",
        "geneID",
        "value",
        "cell_line",
        "Responder",
        "Baseline",
        "GSE_ID",
        "Mouse_treatment",
    }
    missing = expected - set(raw.columns)
    if missing:
        raise SystemExit(f"unexpected TISMO table columns, missing {missing}")

    models = per_model_table(raw)
    if len(models) != 64:
        raise SystemExit(f"expected 64 TISMO slices, got {len(models)}")

    FIG.mkdir(parents=True, exist_ok=True)

    median_sign = sign_counts(models["delta_median"])
    mean_sign = sign_counts(models["delta_mean"])
    w_med = wilcoxon_delta(models["delta_median"])
    w_mean = wilcoxon_delta(models["delta_mean"])

    last = last_timepoint(models)
    last_med = sign_counts(last["delta_median"])
    last_mean = sign_counts(last["delta_mean"])

    no_gse124821 = models[models["gse_id"] != "GSE124821"]
    drop_zero = models[~((models["median_baseline"] == 0) & (models["median_icb"] == 0))]

    by_line = collapse_unique_samples(raw, models, "cell_line")
    by_indep = collapse_unique_samples(raw, models, "independent_key")
    by_gse = models.groupby("gse_id", as_index=False)["delta_median"].median()

    split_rows = []
    for loc, col in (("median", "delta_median"), ("mean", "delta_mean")):
        for grp, sub in models.groupby("cancer_group"):
            c = sign_counts(sub[col])
            split_rows.append(
                {
                    "cancer_group": grp,
                    "location": loc,
                    "n_slices": int(len(sub)),
                    "n_cell_lines": int(sub["cell_line"].nunique()),
                    **c,
                }
            )
    split = pd.DataFrame(split_rows).sort_values(["location", "n_slices"], ascending=[True, False])

    sensitivity = pd.DataFrame(
        [
            row_from_counts("sign_test_median_64_slices", median_sign, {"note": "primary honest sign test"}),
            row_from_counts("sign_test_mean_64_slices", mean_sign, {"note": "reproduces user 49/64 count"}),
            row_from_counts(
                "wilcoxon_signed_rank_median_deltas",
                {"n_models": w_med["n"], "n_up": median_sign["n_up"], "n_down": median_sign["n_down"], "n_tie": median_sign["n_tie"], "n_tested": w_med["n"], "sign_p_two_sided": w_med["p_two_sided"], "sign_p_greater": np.nan},
                {"note": "Wilcoxon on median deltas; zeros dropped"},
            ),
            row_from_counts(
                "wilcoxon_signed_rank_mean_deltas",
                {"n_models": w_mean["n"], "n_up": mean_sign["n_up"], "n_down": mean_sign["n_down"], "n_tie": mean_sign["n_tie"], "n_tested": w_mean["n"], "sign_p_two_sided": w_mean["p_two_sided"], "sign_p_greater": np.nan},
                {"note": "reproduces user p=5.8e-5 (this is Wilcoxon, not a sign test)"},
            ),
            row_from_counts("sign_test_median_last_timepoint", last_med, {"note": f"n_slices={len(last)}"}),
            row_from_counts("sign_test_mean_last_timepoint", last_mean, {"note": f"n_slices={len(last)}"}),
            row_from_counts(
                "sign_test_median_drop_GSE124821",
                sign_counts(no_gse124821["delta_median"]),
                {"note": "GSE124821 is 19/64 time-course slices"},
            ),
            row_from_counts(
                "sign_test_median_drop_both_zero",
                sign_counts(drop_zero["delta_median"]),
                {"note": "drop slices with median Tacstd2 = 0 in both arms"},
            ),
            row_from_counts(
                "sign_test_median_unique_cell_line",
                sign_counts(by_line["delta_median"]),
                {"note": "unique SRX pooled per cell line"},
            ),
            row_from_counts(
                "sign_test_median_independent_key",
                sign_counts(by_indep["delta_median"]),
                {"note": "collapse ICB agent + timepoint; unique SRX"},
            ),
            row_from_counts(
                "sign_test_median_unique_GSE",
                sign_counts(by_gse["delta_median"]),
                {"note": "median of slice-level median deltas per GSE"},
            ),
        ]
    )

    shared_base = (
        raw[raw["Baseline"] == 1]
        .groupby("Samples")["cell_line"]
        .nunique()
    )
    audit = {
        "source": "TISMO /gene/downVivoExprn gene=Tacstd2 type=3 (table)",
        "url": "https://tismo.pku-genomics.org/",
        "n_rows": int(len(raw)),
        "n_unique_samples": int(raw["Samples"].nunique()),
        "n_tismo_slices": int(models["model"].nunique()),
        "n_cell_lines": int(models["cell_line"].nunique()),
        "n_gse": int(models["gse_id"].nunique()),
        "n_gse124821_slices": int((models["gse_id"] == "GSE124821").sum()),
        "n_baseline_samples_reused_across_slices": int((shared_base > 1).sum()),
        "gene": "Tacstd2",
        "pairing": "within-slice Baseline==1 vs Baseline==0 (ICB-treated, R+NR pooled)",
        "primary_test": "two-sided exact binomial sign test on median(ICB)-median(baseline); ties dropped",
        "user_claim": {"n_up": 49, "n_models": 64, "p": 5.8e-5, "label": "sign test"},
        "user_claim_reproduced_as": {
            "n_up_mean": mean_sign["n_up"],
            "n_down_mean": mean_sign["n_down"],
            "sign_test_p_two_sided_on_means": mean_sign["sign_p_two_sided"],
            "wilcoxon_p_two_sided_on_mean_deltas": w_mean["p_two_sided"],
        },
    }

    summary = {
        "claim": "A4 TISMO Tacstd2 after ICB, all 64 slices",
        "n_slices": 64,
        "primary_median_sign_test": median_sign,
        "mean_sign_test_user_count": mean_sign,
        "wilcoxon_mean_deltas_user_p": w_mean,
        "wilcoxon_median_deltas": w_med,
        "cancer_groups": {
            g: int((models["cancer_group"] == g).sum()) for g in sorted(models["cancer_group"].unique())
        },
        "verdict": (
            "User 49/64 is the mean-based sign count. User p=5.8e-5 is Wilcoxon "
            "signed-rank on the 64 mean differences, not a sign test. Honest "
            "median sign test is 40 up / 18 down / 6 ties (p=5.35e-3). No cancer "
            "type is significant on its own; mammary (29 slices) is 16/12/1, p=0.57."
        ),
    }

    models.sort_values(["cancer_group", "cell_line", "gse_id", "model"]).to_csv(
        ROOT / "per_model.tsv", sep="\t", index=False
    )
    split.to_csv(ROOT / "cancer_type_split.tsv", sep="\t", index=False)
    sensitivity.to_csv(ROOT / "sensitivity.tsv", sep="\t", index=False)
    by_indep.sort_values("unit").to_csv(ROOT / "independent_units.tsv", sep="\t", index=False)
    write_json(ROOT / "summary.json", summary)
    write_json(ROOT / "audit.json", audit)

    plot_waterfall(models, FIG / "fig_waterfall_median_delta.png")
    plot_cancer_split(split, FIG / "fig_cancer_type_sign.png")

    print("=== primary (median sign test, 64 slices) ===")
    print(median_sign)
    print("=== user-style (mean sign count / Wilcoxon on means) ===")
    print(mean_sign)
    print("wilcoxon means", w_mean)
    print("=== cancer split (median) ===")
    print(split[split["location"] == "median"].to_string(index=False))
    print("wrote", ROOT)


if __name__ == "__main__":
    main()
