#!/usr/bin/env python3
"""Median-split odds ratios for high CLDN4 versus objective response.

The slide claim is OR = 0.42 for objective response in CLDN4-high versus
CLDN4-low (OR < 1 is that direction). This script puts one primary OR on
every open cohort that has CLDN4 and a responder label.

Extra cuts run only on mixed or underpowered cohorts, judged from the
primary rank tests already in open_gene_tests.tsv:

- Cho (GSE126044): n_responder = 4 and the GEO response string disagrees
  on GSM3589680. Label sensitivity, continuous CLDN4, TACSTD2, and a
  KRT18/KRT19 residual are reported. The Data S9 label stays the forest point.
- Gide pretreatment: CLDN4 AUC p = 0.070. Continuous CLDN4, TACSTD2, and
  the keratin residual are reported. Pretreatment is the forest point.
  Pre+on is not a second cohort.

Hugo, Jung, VanAllen, Liu, and Mariathasan keep the primary median split
only. OAK/POPLAR are not in the per-sample table and are not imputed.
"""

from __future__ import annotations

import gzip
import json
import os
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/ici_11_cohorts"
FIG = OUT / "figures"
CACHE = Path(os.environ.get("ICI11_CACHE", "/tmp/ici11"))
PER = OUT / "per_sample_open.tsv"
SLIDE_OR = 0.42

# Forest order. Gide uses pretreatment sample ids (suffix _Pre).
FOREST = [
    ("Hugo", "Hugo GSE78220", None),
    ("Jung", "Jung GSE135222", None),
    ("Cho", "Cho GSE126044", None),
    ("Liu", "Liu DFCI 2019", None),
    ("VanAllen", "VanAllen DFCI 2015", None),
    ("Gide_pre", "Gide pretreatment", None),
    ("Mariathasan", "Mariathasan IMvigor210", None),
]


def fetch_json(url: str, payload: dict | None = None):
    data = None if payload is None else json.dumps(payload).encode()
    headers = {"Accept": "application/json", "User-Agent": "Mozilla/5.0"}
    if payload is not None:
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=headers)
    with urllib.request.urlopen(req, timeout=180) as resp:
        return json.load(resp)


def load_frames() -> dict[str, pd.DataFrame]:
    raw = pd.read_csv(PER, sep="\t")
    frames = {}
    for cohort, sub in raw.groupby("cohort"):
        frames[cohort] = sub.copy()
    gide = frames["Gide_all"]
    frames["Gide_pre"] = gide[gide.sample_id.str.endswith("_Pre")].copy()
    frames["Gide_on"] = gide[gide.sample_id.str.endswith("_On")].copy()
    return frames


def split_counts(y: np.ndarray, x: np.ndarray) -> dict | None:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask].astype(int)
    if len(x) < 4 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return None
    med = float(np.median(x))
    high = x > med
    # Values tied at the median stay low, so the high arm is strictly above.
    a = int(((high) & (y == 1)).sum())
    b = int(((high) & (y == 0)).sum())
    c = int(((~high) & (y == 1)).sum())
    d = int(((~high) & (y == 0)).sum())
    return {"n": int(len(y)), "n_high": int(high.sum()), "n_low": int((~high).sum()),
            "high_R": a, "high_NR": b, "low_R": c, "low_NR": d, "median": med}


def or_from_counts(counts: dict) -> dict:
    a, b, c, d = counts["high_R"], counts["high_NR"], counts["low_R"], counts["low_NR"]
    zero = min(a, b, c, d) == 0
    corr = 0.5 if zero else 0.0
    aa, bb, cc, dd = a + corr, b + corr, c + corr, d + corr
    lor = float(np.log((aa * dd) / (bb * cc)))
    se = float(np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd))
    _, p = fisher_exact([[a, b], [c, d]])
    z = (lor - np.log(SLIDE_OR)) / se
    # two-sided normal test of H0: OR = 0.42
    from math import erfc

    p_slide = float(erfc(abs(z) / np.sqrt(2)))
    out = dict(counts)
    out.update(
        {
            "or_response_high_vs_low": float(np.exp(lor)),
            "ci_low": float(np.exp(lor - 1.96 * se)),
            "ci_high": float(np.exp(lor + 1.96 * se)),
            "log_or": lor,
            "se_log_or": se,
            "haldane": zero,
            "fisher_p": float(p),
            "z_vs_0.42": float(z),
            "p_vs_0.42": p_slide,
        }
    )
    return out


def median_or(frame: pd.DataFrame, gene: str) -> dict | None:
    if gene not in frame.columns:
        return None
    counts = split_counts(frame["y"].to_numpy(float), frame[gene].to_numpy(float))
    if counts is None:
        return None
    return or_from_counts(counts)


def logit_per_sd(y: np.ndarray, x: np.ndarray) -> dict | None:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask].astype(float)
    y = y[mask].astype(float)
    if len(x) < 8 or np.unique(y).size < 2 or np.std(x) == 0:
        return None
    z = (x - x.mean()) / x.std(ddof=0)
    try:
        fit = sm.Logit(y, sm.add_constant(z)).fit(disp=False, maxiter=100)
    except Exception as exc:  # separation or singular
        return {"ok": False, "reason": str(exc), "n": int(len(y))}
    b = float(fit.params[1])
    se = float(fit.bse[1])
    return {
        "ok": True,
        "n": int(len(y)),
        "or_per_sd": float(np.exp(b)),
        "ci_low": float(np.exp(b - 1.96 * se)),
        "ci_high": float(np.exp(b + 1.96 * se)),
        "p": float(fit.pvalues[1]),
        "log_or": b,
        "se": se,
    }


def log_expr(series: pd.Series, already_log: bool) -> np.ndarray:
    x = series.to_numpy(float)
    if already_log:
        return x
    return np.log2(np.clip(x, 0, None) + 1)


def keratin_residual(y_gene: np.ndarray, k18: np.ndarray, k19: np.ndarray) -> np.ndarray:
    mask = np.isfinite(y_gene) & np.isfinite(k18) & np.isfinite(k19)
    resid = np.full(len(y_gene), np.nan)
    X = np.column_stack([np.ones(mask.sum()), k18[mask], k19[mask]])
    beta, *_ = np.linalg.lstsq(X, y_gene[mask], rcond=None)
    resid[mask] = y_gene[mask] - X @ beta
    return resid


def cho_keratins() -> pd.DataFrame:
    """GSM-level KRT18/KRT19 raw counts for GSE126044."""
    series = CACHE / "GSE126044_series_matrix.txt.gz"
    counts = CACHE / "GSE126044_counts.txt.gz"
    gsm, title = [], []
    with gzip.open(series, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_geo_accession"):
                gsm = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_title"):
                title = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
    col_of = {g: t.replace("RNA-seq_", "") for g, t in zip(gsm, title)}
    wanted = {"KRT18": {}, "KRT19": {}}
    with gzip.open(counts, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        samples = header[1:]
        for line in handle:
            gene, *vals = line.rstrip("\n").split("\t")
            if gene in wanted:
                wanted[gene] = {s: float(v) for s, v in zip(samples, vals)}
    rows = []
    for sample_id, col in col_of.items():
        rows.append(
            {
                "sample_id": sample_id,
                "KRT18": wanted["KRT18"].get(col, np.nan),
                "KRT19": wanted["KRT19"].get(col, np.nan),
            }
        )
    return pd.DataFrame(rows)


def gide_keratins() -> pd.DataFrame:
    mol = fetch_json(
        "https://www.cbioportal.org/api/molecular-profiles/mel_iatlas_gide_2019_rna_seq_mrna/molecular-data/fetch",
        {"entrezGeneIds": [3875, 3880], "sampleListId": "mel_iatlas_gide_2019_all"},
    )
    gene_of = {3875: "KRT18", 3880: "KRT19"}
    rows = [
        {"sample_id": rec["sampleId"], "gene": gene_of[rec["entrezGeneId"]], "value": rec["value"]}
        for rec in mol
    ]
    wide = pd.DataFrame(rows).pivot_table(index="sample_id", columns="gene", values="value", aggfunc="first")
    return wide.reset_index()


def pool(rows: list[dict]) -> dict:
    lor = np.array([r["log_or"] for r in rows], float)
    se = np.array([r["se_log_or"] for r in rows], float)
    w = 1 / se**2
    lor_fe = float(np.sum(w * lor) / np.sum(w))
    se_fe = float(np.sqrt(1 / np.sum(w)))
    q = float(np.sum(w * (lor - lor_fe) ** 2))
    k = len(rows)
    df = k - 1
    c = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    w_re = 1 / (se**2 + tau2)
    lor_re = float(np.sum(w_re * lor) / np.sum(w_re))
    se_re = float(np.sqrt(1 / np.sum(w_re)))
    i2 = max(0.0, (q - df) / q) if q > 0 else 0.0

    def pack(name, lor_hat, se_hat):
        z = (lor_hat - np.log(SLIDE_OR)) / se_hat
        from math import erfc

        return {
            "model": name,
            "k": k,
            "or": float(np.exp(lor_hat)),
            "ci_low": float(np.exp(lor_hat - 1.96 * se_hat)),
            "ci_high": float(np.exp(lor_hat + 1.96 * se_hat)),
            "log_or": lor_hat,
            "se": se_hat,
            "z_vs_0.42": float(z),
            "p_vs_0.42": float(erfc(abs(z) / np.sqrt(2))),
            "Q": q,
            "I2": float(i2),
            "tau2": float(tau2),
        }

    return {"fixed": pack("fixed", lor_fe, se_fe), "random": pack("random", lor_re, se_re)}


def draw_forest(primary: pd.DataFrame, sensitivity: dict | None, pool_re: dict) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    rows = list(primary.itertuples(index=False))
    labels = [r.label for r in rows]
    ors = [r.or_response_high_vs_low for r in rows]
    lo = [r.ci_low for r in rows]
    hi = [r.ci_high for r in rows]
    if sensitivity is not None:
        labels.append("Cho GEO label (not pooled)")
        ors.append(sensitivity["or_response_high_vs_low"])
        lo.append(sensitivity["ci_low"])
        hi.append(sensitivity["ci_high"])
    labels.append(f"Random-effects pool (k={pool_re['k']})")
    ors.append(pool_re["or"])
    lo.append(pool_re["ci_low"])
    hi.append(pool_re["ci_high"])

    ypos = np.arange(len(labels))[::-1]
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    for i, (o, a, b) in enumerate(zip(ors, lo, hi)):
        y = ypos[i]
        color = "#333333"
        marker = "o"
        if i == len(labels) - 1:
            color = "#1f4e79"
            marker = "D"
        elif sensitivity is not None and i == len(labels) - 2:
            color = "#888888"
            marker = "o"
        ax.plot([a, b], [y, y], color=color, lw=1.4)
        ax.plot(o, y, marker, color=color, ms=7 if marker == "D" else 6)
    ax.axvline(1.0, color="#444444", lw=0.8)
    ax.axvline(SLIDE_OR, color="#a33b32", lw=1.0, ls="--")
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels)
    ax.set_xscale("log")
    ax.set_xlabel("OR of objective response, CLDN4 high vs low")
    # Cho's Haldane interval starts near 0.002; keep that whisker on the axis.
    ax.set_xlim(0.0015, 20)
    ax.text(SLIDE_OR, ypos[0] + 0.7, "slide 0.42", color="#a33b32", ha="center", fontsize=8)
    ax.set_title("Open ICI cohorts: high CLDN4 and response")
    fig.tight_layout()
    fig.savefig(FIG / "cldn4_or_forest.png", dpi=160)
    fig.savefig(FIG / "cldn4_or_forest.pdf")
    plt.close(fig)


def main() -> None:
    frames = load_frames()
    primary_rows = []
    for key, label, _ in FOREST:
        est = median_or(frames[key], "CLDN4")
        if est is None:
            raise SystemExit(f"no primary OR for {key}")
        est["cohort"] = key
        est["label"] = label
        est["spec"] = "primary_median_split"
        est["in_pool"] = True
        primary_rows.append(est)

    # Footnote only: Gide pre+on mixes on-treatment RNA. Not pooled.
    gide_all = median_or(frames["Gide_all"], "CLDN4")
    gide_all["cohort"] = "Gide_all"
    gide_all["label"] = "Gide pre+on (not pooled)"
    gide_all["spec"] = "footnote_not_pooled"
    gide_all["in_pool"] = False

    pooled = pool(primary_rows)

    # Cho label sensitivity: GEO calls GSM3589680 a responder.
    cho = frames["Cho"].copy()
    cho_geo = cho.copy()
    flip = cho_geo.sample_id == "GSM3589680"
    if int(flip.sum()) != 1:
        raise SystemExit("GSM3589680 missing from Cho")
    cho_geo.loc[flip, "y"] = 1
    cho_geo.loc[flip, "response"] = "responder"
    cho_sens = median_or(cho_geo, "CLDN4")
    cho_sens["cohort"] = "Cho"
    cho_sens["label"] = "Cho GEO label (not pooled)"
    cho_sens["spec"] = "label_sensitivity_GSM3589680"
    cho_sens["in_pool"] = False

    # Continuous and keratin cuts: Cho and Gide pretreatment only.
    extras = []

    def add_extra(cohort, spec, gene, result):
        if result is None:
            return
        result = dict(result)
        result.update({"cohort": cohort, "spec": spec, "gene": gene})
        extras.append(result)

    # Cho counts are linear. Gide cBio values are already log2 UQ.
    cho_log = log_expr(cho["CLDN4"], already_log=False)
    cho_log_t = log_expr(cho["TACSTD2"], already_log=False)
    add_extra("Cho", "continuous_log2_per_sd", "CLDN4", logit_per_sd(cho["y"].to_numpy(float), cho_log))
    add_extra("Cho", "continuous_log2_per_sd", "TACSTD2", logit_per_sd(cho["y"].to_numpy(float), cho_log_t))
    add_extra("Cho", "median_split", "TACSTD2", median_or(cho, "TACSTD2"))
    add_extra(
        "Cho",
        "median_split_GEO_label",
        "CLDN4",
        cho_sens,
    )
    add_extra(
        "Cho",
        "continuous_log2_per_sd_GEO_label",
        "CLDN4",
        logit_per_sd(cho_geo["y"].to_numpy(float), cho_log),
    )

    krt = cho_keratins()
    cho_k = cho.merge(krt, on="sample_id", how="left")
    resid = keratin_residual(
        log_expr(cho_k["CLDN4"], False),
        log_expr(cho_k["KRT18"], False),
        log_expr(cho_k["KRT19"], False),
    )
    cho_k = cho_k.copy()
    cho_k["CLDN4_krt"] = resid
    add_extra("Cho", "keratin_residual_KRT18_KRT19_median_split", "CLDN4", median_or(cho_k, "CLDN4_krt"))
    add_extra(
        "Cho",
        "keratin_residual_KRT18_KRT19_per_sd",
        "CLDN4",
        logit_per_sd(cho_k["y"].to_numpy(float), cho_k["CLDN4_krt"].to_numpy(float)),
    )

    gpre = frames["Gide_pre"]
    add_extra("Gide_pre", "continuous_deposited_per_sd", "CLDN4", logit_per_sd(gpre["y"].to_numpy(float), gpre["CLDN4"].to_numpy(float)))
    add_extra("Gide_pre", "continuous_deposited_per_sd", "TACSTD2", logit_per_sd(gpre["y"].to_numpy(float), gpre["TACSTD2"].to_numpy(float)))
    add_extra("Gide_pre", "median_split", "TACSTD2", median_or(gpre, "TACSTD2"))
    gk = gide_keratins()
    gpre_k = gpre.merge(gk, on="sample_id", how="left")
    gresid = keratin_residual(
        gpre_k["CLDN4"].to_numpy(float),
        gpre_k["KRT18"].to_numpy(float),
        gpre_k["KRT19"].to_numpy(float),
    )
    gpre_k = gpre_k.copy()
    gpre_k["CLDN4_krt"] = gresid
    add_extra("Gide_pre", "keratin_residual_KRT18_KRT19_median_split", "CLDN4", median_or(gpre_k, "CLDN4_krt"))
    add_extra(
        "Gide_pre",
        "keratin_residual_KRT18_KRT19_per_sd",
        "CLDN4",
        logit_per_sd(gpre_k["y"].to_numpy(float), gpre_k["CLDN4_krt"].to_numpy(float)),
    )

    primary_df = pd.DataFrame(primary_rows)
    footnote = pd.DataFrame([gide_all, {**cho_sens}])
    extra_df = pd.DataFrame(extras)
    pool_df = pd.DataFrame([pooled["fixed"], pooled["random"]])

    primary_df.to_csv(OUT / "forest_primary_or.tsv", sep="\t", index=False)
    footnote.to_csv(OUT / "forest_not_pooled.tsv", sep="\t", index=False)
    extra_df.to_csv(OUT / "mixed_underpowered_extras.tsv", sep="\t", index=False)
    pool_df.to_csv(OUT / "forest_pool.tsv", sep="\t", index=False)
    draw_forest(primary_df, cho_sens, pooled["random"])

    show = primary_df[
        ["label", "n", "high_R", "high_NR", "low_R", "low_NR", "or_response_high_vs_low", "ci_low", "ci_high", "fisher_p", "p_vs_0.42"]
    ]
    print(show.to_string(index=False))
    print("--- pool ---")
    print(pool_df.to_string(index=False))
    print("--- extras ---")
    cols = [c for c in ["cohort", "spec", "gene", "or_response_high_vs_low", "or_per_sd", "ci_low", "ci_high", "fisher_p", "p", "ok"] if c in extra_df.columns]
    print(extra_df[cols].to_string(index=False))
    print("--- gide all footnote ---")
    print(gide_all["or_response_high_vs_low"], gide_all["ci_low"], gide_all["ci_high"], gide_all["fisher_p"])


if __name__ == "__main__":
    main()
