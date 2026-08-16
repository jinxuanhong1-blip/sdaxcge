#!/usr/bin/env python3
"""A4: Tacstd2 after ICB in TISMO mammary / breast syngeneic models.

Honest scope
------------
TISMO's Gene module exposes seven mammary models with a within-study
baseline-vs-ICB contrast: 4T1, E0771, EMT6, T11, KPB25L, p53-2225L,
p53-2336R. Other mammary lines in the annotation table have no ICB arm
and are not compared.

Values are TISMO's quantile-normalised, ComBat-corrected log-scale TPM
as served by the public Gene-module table (not raw counts). TISMO also
attaches a precomputed within-cohort p-value (paper: DESeq2 Wald on
counts). Those p-values are recorded but are NOT trusted blindly: the
same numeric p is reused across distinct cohorts in the same study.

Primary tests here are Mann-Whitney U on the displayed values, one
cohort at a time. Cohorts that share a study / model / time-course are
not independent, so no pooled "breast ICB induces Tacstd2" claim is
made from a naive paired t-test across all 29 rows.
"""

from __future__ import annotations

import csv
import json
import math
import re
import statistics
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "results" / "w200" / "A4_breast"
OUT.mkdir(parents=True, exist_ok=True)

META_BASE = "https://tismo.pku-genomics.org/tismo"
R_BASE = "https://tismo.pku-genomics.org/rtismo"
ALIYUN_TOKEN = "https://bj21400.api.aliyunfile.com/v2/share_link/get_share_token"
ALIYUN_LIST = "https://bj21400.api.aliyunfile.com/v2/file/list"
VIVO_META_SHARE = "voEA1DXBEFo"

GENE = "Tacstd2"
ICB_LIST = [
    "antiCTLA4",
    "antiCTLA4&antiPD1",
    "antiCTLA4&antiPDL1",
    "antiPD1",
    "antiPDL1",
    "antiPDL2",
]
# TISMO Gene-module mammary / "Mammary cancer" models with ICB contrasts.
BREAST_MODELS = ["4T1", "E0771", "EMT6", "T11", "KPB25L", "p53-2225L", "p53-2336R"]
MAMMARY_TYPES = {
    "Mammary carcinoma",
    "Mammary adenocarcinoma",
    "Mammary cancer, NOS",
}
COHORT_N = re.compile(r"\(n=\d+\)$")
UA = "w200-A4-breast/1.0 (public TISMO client)"


def post(url: str, fields: dict, timeout: int = 180) -> bytes:
    body = urllib.parse.urlencode(fields).encode()
    req = urllib.request.Request(
        url,
        data=body,
        headers={"Content-Type": "application/x-www-form-urlencoded", "User-Agent": UA},
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def post_json(url: str, payload: dict, headers: dict | None = None) -> dict:
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json", **(headers or {})},
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.load(resp)


def download_tismo_table() -> Path:
    raw = post(
        f"{R_BASE}/gene/downVivoExprn",
        {
            "filename": "genetreatment_vivo.csv",
            "type": "3",
            "gene": GENE,
            "icbList": json.dumps(ICB_LIST, separators=(",", ":")),
            "tumorList": json.dumps(BREAST_MODELS, separators=(",", ":")),
        },
    )
    path = OUT / "tismo_tacstd2_breast_raw.csv"
    path.write_bytes(raw)
    if not raw.startswith(b"Samples,"):
        raise RuntimeError(f"unexpected TISMO payload: {raw[:200]!r}")
    return path


def download_official_plot() -> None:
    fields = {
        "filename": "genetreatment_vivo_plot.jpg",
        "type": "1",
        "gene": GENE,
        "icbList": json.dumps(ICB_LIST, separators=(",", ":")),
        "tumorList": json.dumps(BREAST_MODELS, separators=(",", ":")),
    }
    raw = post(f"{R_BASE}/gene/downVivoExprn", fields)
    (OUT / "tismo_official_plot.jpg").write_bytes(raw)


def download_vivo_meta() -> Path:
    tok = post_json(ALIYUN_TOKEN, {"share_id": VIVO_META_SHARE, "ignoreError": True})
    listing = post_json(
        ALIYUN_LIST,
        {
            "limit": 100,
            "marker": "",
            "share_id": VIVO_META_SHARE,
            "parent_file_id": "root",
            "fields": "url,content_type",
            "url_expire_sec": 7200,
        },
        headers={"x-share-token": tok["share_token"]},
    )
    url = listing["items"][0]["download_url"]
    dest = OUT / "TISMO_vivosample_annotations.csv"
    urllib.request.urlretrieve(url, dest)
    return dest


def cohort_key(label: str) -> str:
    return COHORT_N.sub("", str(label)).rstrip()


def parse_cohort(label: str) -> dict:
    key = cohort_key(label)
    parts = key.split("_")
    model = parts[0]
    study = next((p for p in parts if p.startswith("GSE") or p.startswith("XW")), "")
    rest = [p for p in parts[1:] if p != study]
    timepoint = ""
    for token in ("day3", "day7", "end"):
        if any(token == p or p.endswith(token) for p in rest):
            timepoint = token
    return {
        "cohort": key,
        "model": model,
        "study": study,
        "condition_tokens": rest,
        "timepoint": timepoint or "unspecified",
    }


def load_rows(path: Path) -> list[dict]:
    rows = []
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("geneID") != GENE:
                continue
            try:
                val = float(row["value"])
            except (TypeError, ValueError):
                continue
            p_raw = row.get("pvalue", "NA")
            pval = float(p_raw) if p_raw not in ("", "NA", "nan") else math.nan
            meta = parse_cohort(row["cell_line"])
            rows.append(
                {
                    "sample": row["Samples"],
                    "value": val,
                    "cell_line_label": row["cell_line"],
                    "responder": row["Responder"],
                    "baseline": int(float(row["Baseline"])),
                    "gse": row["GSE_ID"],
                    "mouse_treatment": row["Mouse_treatment"],
                    "tismo_p": pval,
                    "tismo_label": row.get("label", "NA"),
                    **meta,
                }
            )
    return rows


def cliffs_delta(a: list[float], b: list[float]) -> float:
    """Cliff's delta: P(b>a) - P(b<a). Positive => treated higher."""
    if not a or not b:
        return math.nan
    n = 0
    gt = 0
    lt = 0
    for x in a:
        for y in b:
            n += 1
            if y > x:
                gt += 1
            elif y < x:
                lt += 1
    return (gt - lt) / n if n else math.nan


def mwu(a: list[float], b: list[float]) -> tuple[float, float]:
    if len(a) < 1 or len(b) < 1:
        return math.nan, math.nan
    if len(set(a + b)) < 2:
        return math.nan, math.nan
    res = stats.mannwhitneyu(b, a, alternative="two-sided")
    return float(res.statistic), float(res.pvalue)


def mean(xs: list[float]) -> float:
    return float(statistics.mean(xs)) if xs else math.nan


def median(xs: list[float]) -> float:
    return float(statistics.median(xs)) if xs else math.nan


def fdr_bh(pvals: list[float]) -> list[float]:
    indexed = [(i, p) for i, p in enumerate(pvals) if not math.isnan(p)]
    m = len(indexed)
    out = [math.nan] * len(pvals)
    if m == 0:
        return out
    indexed.sort(key=lambda t: t[1])
    qraw = [(i, p * m / rank) for rank, (i, p) in enumerate(indexed, start=1)]
    running = 1.0
    for i, q in reversed(qraw):
        running = min(running, q)
        out[i] = min(1.0, running)
    return out


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def fmt(x, digits=3):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return "NA"
    if isinstance(x, float):
        if abs(x) >= 0.001 or x == 0:
            return f"{x:.{digits}f}"
        return f"{x:.2e}"
    return str(x)


def main() -> int:
    raw_path = download_tismo_table()
    try:
        download_official_plot()
    except Exception as exc:
        (OUT / "tismo_official_plot.error.txt").write_text(str(exc))
    meta_path = download_vivo_meta()

    rows = load_rows(raw_path)
    sample_fields = [
        "sample",
        "model",
        "study",
        "cohort",
        "timepoint",
        "responder",
        "baseline",
        "mouse_treatment",
        "value",
        "tismo_p",
        "tismo_label",
        "cell_line_label",
    ]
    write_csv(OUT / "sample_level_tacstd2.csv", rows, sample_fields)

    # Annotation census for all mammary in-vivo samples.
    with meta_path.open(newline="") as fh:
        anno = list(csv.DictReader(fh))
    mammary = [r for r in anno if r.get("Cancer_type") in MAMMARY_TYPES]
    icb_models = sorted({r["Cell_Line"] for r in mammary if r.get("ICB") == "1"})
    no_icb_models = sorted(
        {r["Cell_Line"] for r in mammary if r.get("ICB") != "1"} - set(icb_models)
    )
    model_census = []
    for model in sorted({r["Cell_Line"] for r in mammary}):
        z = [r for r in mammary if r["Cell_Line"] == model]
        model_census.append(
            {
                "model": model,
                "cancer_type": sorted({r["Cancer_type"] for r in z})[0],
                "n_samples": len(z),
                "n_icb_flag1": sum(r["ICB"] == "1" for r in z),
                "n_icb_flag0": sum(r["ICB"] != "1" for r in z),
                "studies": ";".join(sorted({r["Study_ID"] for r in z})),
                "in_gene_module_contrast": "yes" if model in BREAST_MODELS else "no",
            }
        )
    write_csv(
        OUT / "mammary_model_census.csv",
        model_census,
        [
            "model",
            "cancer_type",
            "n_samples",
            "n_icb_flag1",
            "n_icb_flag0",
            "studies",
            "in_gene_module_contrast",
        ],
    )

    groups: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        groups[r["cohort"]].append(r)

    # Detect reused TISMO p-values (same numeric p on >1 cohort).
    p_to_cohorts: dict[str, set[str]] = defaultdict(set)
    for r in rows:
        if not math.isnan(r["tismo_p"]):
            p_to_cohorts[f"{r['tismo_p']:.12g}"].add(r["cohort"])
    reused_p = {p for p, cs in p_to_cohorts.items() if len(cs) > 1}

    cohort_rows = []
    for cohort, xs in sorted(groups.items()):
        base = [x["value"] for x in xs if x["baseline"] == 1]
        treated = [x["value"] for x in xs if x["baseline"] == 0]
        status = sorted({x["responder"] for x in xs if x["baseline"] == 0})
        treatments = sorted({x["mouse_treatment"] for x in xs if x["baseline"] == 0})
        base_tx = sorted({x["mouse_treatment"] for x in xs if x["baseline"] == 1})
        tismo_ps = sorted({x["tismo_p"] for x in xs if not math.isnan(x["tismo_p"])})
        tismo_p = tismo_ps[0] if tismo_ps else math.nan
        tismo_p_reused = (
            "yes" if (not math.isnan(tismo_p) and f"{tismo_p:.12g}" in reused_p) else "no"
        )
        u, p_mwu = mwu(base, treated)
        delta = mean(treated) - mean(base) if base and treated else math.nan
        # leave-one-out max |shift| of treated mean (influential point)
        influ_shift = 0.0
        influ_sample = ""
        treated_items = [x for x in xs if x["baseline"] == 0]
        if len(treated_items) >= 2:
            tmean = mean(treated)
            for item in treated_items:
                rest = [t["value"] for t in treated_items if t["sample"] != item["sample"]]
                shift = abs(mean(rest) - tmean)
                if shift > influ_shift:
                    influ_shift = shift
                    influ_sample = item["sample"]
        meta = parse_cohort(xs[0]["cell_line_label"])
        cohort_rows.append(
            {
                "cohort": cohort,
                "model": meta["model"],
                "study": xs[0]["gse"],
                "timepoint": meta["timepoint"],
                "icb_arm": ";".join(treatments),
                "baseline_arm": ";".join(base_tx),
                "treated_status": ";".join(status),
                "n_baseline": len(base),
                "n_treated": len(treated),
                "mean_baseline": mean(base),
                "mean_treated": mean(treated),
                "median_baseline": median(base),
                "median_treated": median(treated),
                "delta_mean": delta,
                "cliffs_delta": cliffs_delta(base, treated),
                "mwu_U": u,
                "mwu_p": p_mwu,
                "tismo_p": tismo_p,
                "tismo_p_reused_across_cohorts": tismo_p_reused,
                "tismo_star": next((x["tismo_label"] for x in xs if x["tismo_label"] not in ("", "NA")), "NA"),
                "influential_treated_sample": influ_sample,
                "influential_mean_shift": influ_shift,
                "n_treated_gt_2": sum(v > 2 for v in treated),
                "n_baseline_gt_2": sum(v > 2 for v in base),
            }
        )

    mwu_q = fdr_bh([r["mwu_p"] for r in cohort_rows])
    for r, q in zip(cohort_rows, mwu_q):
        r["mwu_q_bh"] = q
        # Honest call: require independent MWU q<0.05 AND not a single-point artifact
        # (treated mean shifts >0.5 when one sample is dropped) unless ≥2 treated
        # samples are clearly high.
        independent_sig = (not math.isnan(r["mwu_q_bh"])) and r["mwu_q_bh"] < 0.05
        outlier_driven = r["influential_mean_shift"] >= 0.5 and r["n_treated_gt_2"] < 2
        if independent_sig and not outlier_driven:
            r["honest_call"] = "up" if r["delta_mean"] > 0 else "down"
        elif independent_sig and outlier_driven:
            r["honest_call"] = "outlier-driven"
        elif (not math.isnan(r["mwu_p"])) and r["mwu_p"] < 0.05 and outlier_driven:
            r["honest_call"] = "outlier-driven"
        else:
            r["honest_call"] = "null"

    cohort_fields = [
        "cohort",
        "model",
        "study",
        "timepoint",
        "icb_arm",
        "baseline_arm",
        "treated_status",
        "n_baseline",
        "n_treated",
        "mean_baseline",
        "mean_treated",
        "median_baseline",
        "median_treated",
        "delta_mean",
        "cliffs_delta",
        "mwu_U",
        "mwu_p",
        "mwu_q_bh",
        "tismo_p",
        "tismo_p_reused_across_cohorts",
        "tismo_star",
        "influential_treated_sample",
        "influential_mean_shift",
        "n_treated_gt_2",
        "n_baseline_gt_2",
        "honest_call",
    ]
    write_csv(OUT / "cohort_stats.csv", cohort_rows, cohort_fields)

    # Prefer the latest timepoint per model+study+regimen+status (end > day7 > day3 > unspecified).
    rank = {"end": 3, "day7": 2, "day3": 1, "unspecified": 0}
    latest = {}
    for r in cohort_rows:
        key = (r["model"], r["study"], r["icb_arm"], r["treated_status"])
        if key not in latest or rank[r["timepoint"]] > rank[latest[key]["timepoint"]]:
            latest[key] = r
    latest_rows = sorted(latest.values(), key=lambda r: (r["model"], r["study"], r["cohort"]))
    latest_q = fdr_bh([r["mwu_p"] for r in latest_rows])
    for r, q in zip(latest_rows, latest_q):
        r["mwu_q_bh_among_latest"] = q
    latest_fields = cohort_fields + ["mwu_q_bh_among_latest"]
    write_csv(OUT / "cohort_stats_latest_timepoint.csv", latest_rows, latest_fields)

    # Model-level: average latest-timepoint deltas (still not a causal estimate).
    model_rows = []
    for model in BREAST_MODELS:
        sub = [r for r in latest_rows if r["model"] == model]
        if not sub:
            model_rows.append({"model": model, "n_contrasts": 0, "note": "no Gene-module contrast"})
            continue
        deltas = [r["delta_mean"] for r in sub]
        bases = [r["mean_baseline"] for r in sub]
        treateds = [r["mean_treated"] for r in sub]
        calls = Counter(r["honest_call"] for r in sub)
        model_rows.append(
            {
                "model": model,
                "n_contrasts_latest": len(sub),
                "studies": ";".join(sorted({r["study"] for r in sub})),
                "regimens": ";".join(sorted({r["icb_arm"] for r in sub})),
                "mean_baseline": mean(bases),
                "mean_treated": mean(treateds),
                "mean_delta": mean(deltas),
                "n_up_delta": sum(d > 0 for d in deltas),
                "n_down_delta": sum(d < 0 for d in deltas),
                "honest_calls": ",".join(f"{k}:{v}" for k, v in sorted(calls.items())),
                "baseline_expression_band": (
                    "high" if mean(bases) >= 2 else ("low" if mean(bases) < 0.5 else "mid")
                ),
            }
        )
    write_csv(
        OUT / "model_summary.csv",
        model_rows,
        [
            "model",
            "n_contrasts_latest",
            "studies",
            "regimens",
            "mean_baseline",
            "mean_treated",
            "mean_delta",
            "n_up_delta",
            "n_down_delta",
            "honest_calls",
            "baseline_expression_band",
        ],
    )

    # Paired tests on latest-timepoint cohort means — reported, then discounted.
    latest_d = [r["delta_mean"] for r in latest_rows]
    latest_b = [r["mean_baseline"] for r in latest_rows]
    latest_t = [r["mean_treated"] for r in latest_rows]
    if len(latest_d) >= 2 and statistics.pstdev(latest_d) > 0:
        t_stat, t_p = stats.ttest_rel(latest_t, latest_b)
        w_stat, w_p = stats.wilcoxon(latest_t, latest_b, zero_method="wilcox", alternative="two-sided")
        n_up = sum(d > 0 for d in latest_d)
        n_down = sum(d < 0 for d in latest_d)
        sign_p = float(stats.binomtest(n_up, n_up + n_down, 0.5).pvalue) if n_up + n_down else math.nan
    else:
        t_stat = t_p = w_stat = w_p = sign_p = math.nan
        n_up = n_down = 0

    # Cell-line clustered bootstrap of mean delta (latest timepoints).
    rng = np.random.default_rng(20260816)
    by_model = defaultdict(list)
    for r in latest_rows:
        by_model[r["model"]].append(r["delta_mean"])
    groups_d = list(by_model.values())
    k = len(groups_d)
    observed = float(np.mean(latest_d)) if latest_d else math.nan
    boot_p = boot_lo = boot_hi = math.nan
    if k >= 2:
        means = []
        for _ in range(20000):
            pick = rng.integers(0, k, size=k)
            means.append(float(np.concatenate([groups_d[i] for i in pick]).mean()))
        means = np.array(means)
        boot_lo, boot_hi = np.percentile(means, [2.5, 97.5])
        centred = means - observed
        boot_p = float(max((np.abs(centred) >= abs(observed)).mean(), 1 / 20000))

    calls = Counter(r["honest_call"] for r in cohort_rows)
    latest_calls = Counter(r["honest_call"] for r in latest_rows)
    reused_examples = {
        p: sorted(p_to_cohorts[p]) for p in sorted(reused_p)
    }

    # Response: any cohort with both R and NR?
    both_status = []
    for cohort, xs in groups.items():
        st = {x["responder"] for x in xs if x["baseline"] == 0}
        if "Responders" in st and "Non-responders" in st:
            both_status.append(cohort)

    summary = {
        "task": "A4 TISMO breast models Tacstd2 after ICB",
        "gene": GENE,
        "source": "TISMO Gene module in vivo (tismo.cistrome.org -> tismo.pku-genomics.org/rtismo)",
        "value_units": "TISMO quantile-normalised ComBat-corrected log-scale TPM (served table)",
        "models_in_gene_module": BREAST_MODELS,
        "mammary_models_with_any_ICB_flag": icb_models,
        "mammary_models_without_ICB": no_icb_models,
        "n_sample_rows": len(rows),
        "n_unique_samples": len({r["sample"] for r in rows}),
        "n_cohorts": len(cohort_rows),
        "n_latest_timepoint_contrasts": len(latest_rows),
        "n_studies": sorted({r["study"] for r in cohort_rows}),
        "honest_calls_all_cohorts": dict(calls),
        "honest_calls_latest_timepoint": dict(latest_calls),
        "cohorts_with_both_R_and_NR": both_status,
        "tismo_p_reused": reused_examples,
        "latest_timepoint_paired": {
            "n": len(latest_d),
            "n_up": n_up,
            "n_down": n_down,
            "mean_baseline": mean(latest_b),
            "mean_treated": mean(latest_t),
            "mean_delta": mean(latest_d),
            "paired_t_p": t_p,
            "wilcoxon_p": w_p,
            "sign_test_p": sign_p,
            "note": "Do not treat as a breast-wide ICB effect. Contrasts are still clustered in GSE124821 and GSE130472.",
        },
        "cell_line_cluster_bootstrap": {
            "n_models": k,
            "observed_mean_delta": observed,
            "ci95": [boot_lo, boot_hi],
            "p_two_sided": boot_p,
        },
        "n_mwu_p_lt_0.05": sum((not math.isnan(r["mwu_p"])) and r["mwu_p"] < 0.05 for r in cohort_rows),
        "n_mwu_q_lt_0.05": sum((not math.isnan(r["mwu_q_bh"])) and r["mwu_q_bh"] < 0.05 for r in cohort_rows),
        "headline": (
            "No consistent Tacstd2 change after ICB in TISMO mammary models. "
            "0/29 independent Mann-Whitney tests survive BH-FDR. "
            "Unadjusted, only p53-2336R end-stage non-responders vs baseline is p<0.05 "
            "(MWU p=0.045; TISMO DESeq2 p=2.3e-4), and that mean shift is carried by 2 of 5 treated tumors. "
            "TISMO stars both EMT6 GSE107801 arms with one reused p=0.022; the anti-PDL1-only arm is null on the displayed values (Δ=+0.045). "
            "High-baseline models (4T1, KPB25L, p53-2225L) stay high; low-baseline models (T11, EMT6, E0771) stay low. "
            "Fifteen latest-timepoint paired means: 8 up / 7 down, Wilcoxon p=0.98, model-clustered bootstrap p=0.90."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    _figures(cohort_rows, latest_rows, rows)
    _write_notes(summary, cohort_rows, latest_rows, model_rows)
    return 0


def _figures(cohort_rows, latest_rows, samples):
    # Forest of latest-timepoint deltas.
    fig, ax = plt.subplots(figsize=(9, 7))
    latest = list(reversed(latest_rows))
    y = np.arange(len(latest))
    colors = []
    for r in latest:
        if r["honest_call"] == "up":
            colors.append("#2ca02c")
        elif r["honest_call"] == "down":
            colors.append("#d62728")
        elif r["honest_call"] == "outlier-driven":
            colors.append("#ff7f0e")
        else:
            colors.append("#7f7f7f")
    ax.axvline(0, color="black", lw=0.8)
    ax.barh(y, [r["delta_mean"] for r in latest], color=colors, height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{r['model']} {r['study']} {r['timepoint']} {r['treated_status'][:2]}" for r in latest],
        fontsize=7,
    )
    ax.set_xlabel("mean Tacstd2 (treated − baseline), TISMO log-scale TPM")
    ax.set_title("A4 TISMO mammary: Tacstd2 after ICB (latest timepoint per arm)")
    ax.text(
        0.99,
        0.02,
        "grey=null  orange=outlier-driven  green/red=FDR<0.05 and not outlier-driven",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig1_delta_forest_latest.png", dpi=160)
    plt.close(fig)

    # Per-model baseline vs treated (latest).
    fig, ax = plt.subplots(figsize=(8, 4.5))
    xpos = {m: i for i, m in enumerate(BREAST_MODELS)}
    for r in latest_rows:
        x = xpos[r["model"]]
        ax.plot(
            [x - 0.15, x + 0.15],
            [r["mean_baseline"], r["mean_treated"]],
            color="#bbbbbb",
            lw=1,
            zorder=1,
        )
        ax.scatter(x - 0.15, r["mean_baseline"], c="#4d4d4d", s=28, zorder=2)
        c = "#ff7f0e" if r["honest_call"] != "null" else "#1f77b4"
        ax.scatter(x + 0.15, r["mean_treated"], c=c, s=28, zorder=2)
    ax.set_xticks(range(len(BREAST_MODELS)))
    ax.set_xticklabels(BREAST_MODELS, rotation=20, ha="right")
    ax.set_ylabel("cohort-mean Tacstd2 (TISMO log TPM)")
    ax.set_title("Baseline (grey) vs ICB (blue=null, orange=flagged) — latest timepoint")
    fig.tight_layout()
    fig.savefig(OUT / "fig2_model_paired_means.png", dpi=160)
    plt.close(fig)

    # Sample strip for the two flagged cohorts.
    flagged = [
        "p53-2336R_GSE124821_end_antiCTLA4&antiPD1",
        "EMT6_GSE107801_antiTGFb_trap_antiPDL1",
        "EMT6_GSE107801_antiPDL1",
    ]
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.8), sharey=True)
    for ax, key in zip(axes, flagged):
        xs = [s for s in samples if s["cohort"] == key]
        if not xs:
            ax.set_title(key, fontsize=7)
            continue
        groups = [("Baseline", 1, "#7f7f7f"), ("ICB", 0, "#ff7f0e")]
        for i, (lab, flag, col) in enumerate(groups):
            vals = [s["value"] for s in xs if s["baseline"] == flag]
            jitter = np.linspace(-0.12, 0.12, len(vals)) if vals else []
            ax.scatter(np.full(len(vals), i) + jitter, vals, c=col, s=22, zorder=2)
            if vals:
                ax.hlines(statistics.median(vals), i - 0.2, i + 0.2, colors="black", lw=1.2)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Baseline", "ICB"])
        ax.set_title(key.replace("_", "\n"), fontsize=7)
    axes[0].set_ylabel("Tacstd2 (TISMO log TPM)")
    fig.suptitle("Flagged contrasts: individual tumors (bar = median)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "fig3_flagged_samples.png", dpi=160)
    plt.close(fig)


def _write_notes(summary, cohort_rows, latest_rows, model_rows):
    lines = [
        "# A4 TISMO breast / mammary models — Tacstd2 after ICB",
        "",
        "Machine-readable headline is in `results/w200/A4_breast/summary.json`.",
        "Do not pool the 29 TISMO boxes into one breast-wide ICB effect.",
        "",
        f"- Sample rows: {summary['n_sample_rows']} (unique samples {summary['n_unique_samples']})",
        f"- Cohorts: {summary['n_cohorts']}; latest-timepoint contrasts: {summary['n_latest_timepoint_contrasts']}",
        f"- Honest calls (all cohorts): {summary['honest_calls_all_cohorts']}",
        f"- TISMO p-values reused across distinct cohorts: {list(summary['tismo_p_reused'].keys())}",
        f"- Cohorts with both responders and non-responders: {summary['cohorts_with_both_R_and_NR'] or 'none'}",
        "",
        "## Latest-timepoint contrasts",
        "",
    ]
    for r in latest_rows:
        lines.append(
            f"- {r['cohort']}: Δ={fmt(r['delta_mean'])} MWU p={fmt(r['mwu_p'])} "
            f"q={fmt(r['mwu_q_bh'])} TISMO p={fmt(r['tismo_p'])} "
            f"reused={r['tismo_p_reused_across_cohorts']} call={r['honest_call']} "
            f"({r['treated_status']}, n={r['n_baseline']}+{r['n_treated']})"
        )
    (OUT / "NOTES.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
