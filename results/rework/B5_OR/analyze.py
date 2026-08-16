#!/usr/bin/env python3
"""B5 rework: CLDN4-high vs ICI response at ALL locked cutoffs.

Locked before looking at pooled numbers
--------------------------------------
1. Every cutoff is computed and written: median, tertile T3-vs-T1,
   quartile Q4-vs-Q1, and continuous (logistic OR per 1 SD of log expression).
2. Every binary endpoint that a cohort supports is computed:
   curated R vs NR; RECIST CR/PR vs PD (SD dropped); RECIST CR/PR vs SD+PD.
3. GSE135222_DCB180 (PFS>=180d) is a sensitivity endpoint on the SAME
   patients as GSE135222_Jung. It is never counted as an extra study in
   the primary meta-analysis.
4. Primary pooled comparison to the user number uses median + curated
   R-vs-NR (or cohort-native binary) + DerSimonian–Laird random effects.
   That choice is locked here. Other cutoff×endpoint combinations are
   still reported in full. We do not select the combination closest to
   OR=0.42.
5. Cohorts with a missing gene or non-estimable 2×2 are kept in the
   inventory and are not dropped because they are null.

User-reported claim (not a result of this script)
    CLDN4-high ICI OR = 0.42 [0.18–0.95], k=11.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact, norm

HERE = Path(__file__).resolve().parent

# Independent studies in the primary meta (one row per study).
PRIMARY_COHORTS = [
    "GSE126044",
    "GSE135222_Jung",
    "GSE166449",
    "GSE207422",
    "IMvigor210_Mariathasan",
    "Gide",
    "Liu",
    "Riaz",
    "Braun",
]
LUNG_COHORTS = ["GSE126044", "GSE135222_Jung", "GSE166449", "GSE207422"]
# Same patients as GSE135222_Jung; excluded from primary k.
SENSITIVITY_ONLY = {"GSE135222_DCB180"}

USER_OR = 0.42
USER_LO = 0.18
USER_HI = 0.95
USER_K = 11


def log_expr(raw: pd.Series, assay: str) -> tuple[pd.Series, str]:
    """Return (expression used for cutoffs, scale label).

    Zenodo ICB TPM/FPKM objects are already log-scale (negative values).
    Applying log2(x+1) would clip those to 0 and collapse the median split.
    """
    x = pd.to_numeric(raw, errors="coerce")
    assay = (assay or "").lower()
    if assay in {"log2tpm", "log2(tpm)", "log2_tpm"}:
        return x, "already_log2tpm"
    if x.notna().any() and float(x.min()) < 0:
        return x, "already_log_detected"
    if assay == "counts":
        return np.log2(x.clip(lower=0) + 1.0), "log2_count_plus1"
    return np.log2(x.clip(lower=0) + 1.0), "log2_plus1"


def recode_recist(s: pd.Series) -> pd.Series:
    u = s.astype(str).str.upper().str.strip()
    out = pd.Series(pd.NA, index=s.index, dtype=object)
    out[u.isin(["CR", "PR", "PRCR"])] = "R"
    out[u.isin(["PD"])] = "NR"
    out[u.isin(["SD"])] = "SD"
    return out


def binary_from_cohort(df: pd.DataFrame, endpoint: str) -> pd.Series:
    if endpoint == "curated_R_vs_NR":
        r = df["response"].astype(str).str.upper()
        out = pd.Series(pd.NA, index=df.index, dtype=object)
        out[r.eq("R")] = "R"
        out[r.eq("NR")] = "NR"
        return out
    rec = recode_recist(df["recist"])
    if endpoint == "CRPR_vs_PD":
        out = rec.copy()
        out[out.eq("SD")] = pd.NA
        return out
    if endpoint == "CRPR_vs_SDPD":
        out = rec.copy()
        out[out.eq("SD")] = "NR"
        return out
    raise ValueError(endpoint)


def assign_cutoff(x: pd.Series, cutoff: str) -> pd.Series:
    """Return High/Low/NA. Cutoffs are cohort-internal on the analysis subset."""
    x = pd.to_numeric(x, errors="coerce")
    lab = pd.Series(pd.NA, index=x.index, dtype=object)
    ok = x.notna()
    if ok.sum() < 4:
        return lab
    v = x[ok]
    if cutoff == "median":
        thr = float(v.median())
        lab[ok & (x >= thr)] = "High"
        lab[ok & (x < thr)] = "Low"
        return lab
    if cutoff == "tertile":
        q1, q2 = v.quantile([1 / 3, 2 / 3])
        lab[ok & (x >= q2)] = "High"
        lab[ok & (x <= q1)] = "Low"
        return lab
    if cutoff == "quartile":
        q1, q3 = v.quantile([0.25, 0.75])
        lab[ok & (x >= q3)] = "High"
        lab[ok & (x <= q1)] = "Low"
        return lab
    raise ValueError(cutoff)


def or_from_2x2(a: int, b: int, c: int, d: int) -> dict:
    """OR of response in CLDN4-high vs CLDN4-low.

    a=High-R, b=High-NR, c=Low-R, d=Low-NR.
    Haldane–Anscombe 0.5 correction if any cell is 0.
    """
    cells = [a, b, c, d]
    used_haldane = any(x == 0 for x in cells)
    aa, bb, cc, dd = ([x + 0.5 for x in cells] if used_haldane else cells)
    logor = math.log((aa * dd) / (bb * cc))
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    lo = math.exp(logor - 1.96 * se)
    hi = math.exp(logor + 1.96 * se)
    table = np.array([[a, b], [c, d]], dtype=int)
    if table.sum() == 0 or table.min() < 0:
        fisher_p = np.nan
    else:
        _, fisher_p = fisher_exact(table, alternative="two-sided")
    return {
        "n_high_R": a,
        "n_high_NR": b,
        "n_low_R": c,
        "n_low_NR": d,
        "n_used": int(a + b + c + d),
        "or": math.exp(logor),
        "logor": logor,
        "se_logor": se,
        "or_lo": lo,
        "or_hi": hi,
        "fisher_p": float(fisher_p) if fisher_p == fisher_p else np.nan,
        "haldane": used_haldane,
        "estimable": True,
    }


def logistic_or_per_sd(y: pd.Series, x: pd.Series) -> dict:
    mask = y.notna() & x.notna()
    yy = (y[mask] == "R").astype(int)
    xx = pd.to_numeric(x[mask], errors="coerce")
    if yy.nunique() < 2 or xx.notna().sum() < 8:
        return {"estimable": False, "reason": "too_few_or_one_class", "n_used": int(mask.sum())}
    z = (xx - xx.mean()) / (xx.std(ddof=1) if xx.std(ddof=1) else 1.0)
    X = sm.add_constant(z.to_numpy())
    try:
        fit = sm.Logit(yy.to_numpy(), X).fit(disp=False, maxiter=100)
    except Exception as e:
        return {"estimable": False, "reason": str(e), "n_used": int(mask.sum())}
    logor = float(fit.params[1])
    se = float(fit.bse[1])
    return {
        "n_used": int(mask.sum()),
        "n_R": int(yy.sum()),
        "n_NR": int((1 - yy).sum()),
        "or": math.exp(logor),
        "logor": logor,
        "se_logor": se,
        "or_lo": math.exp(logor - 1.96 * se),
        "or_hi": math.exp(logor + 1.96 * se),
        "wald_p": float(fit.pvalues[1]),
        "estimable": True,
    }


def dl_meta(rows: pd.DataFrame) -> dict:
    """DerSimonian–Laird random-effects + inverse-variance fixed-effect on logOR."""
    d = rows.dropna(subset=["logor", "se_logor"]).copy()
    d = d[d["se_logor"] > 0]
    k = int(len(d))
    if k == 0:
        return {"k": 0, "estimable": False}
    w = 1.0 / (d["se_logor"] ** 2)
    fe_logor = float(np.sum(w * d["logor"]) / np.sum(w))
    q = float(np.sum(w * (d["logor"] - fe_logor) ** 2))
    df = k - 1
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w)) if k > 1 else 0.0
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    w_re = 1.0 / (d["se_logor"] ** 2 + tau2)
    re_logor = float(np.sum(w_re * d["logor"]) / np.sum(w_re))
    re_se = float(math.sqrt(1.0 / np.sum(w_re)))
    fe_se = float(math.sqrt(1.0 / np.sum(w)))
    i2 = max(0.0, (q - df) / q) * 100 if q > 0 and df > 0 else 0.0
    z = re_logor / re_se if re_se else np.nan
    p = float(2 * (1 - norm.cdf(abs(z)))) if z == z else np.nan
    return {
        "k": k,
        "n_total": int(d["n_used"].sum()) if "n_used" in d else None,
        "fe_or": math.exp(fe_logor),
        "fe_or_lo": math.exp(fe_logor - 1.96 * fe_se),
        "fe_or_hi": math.exp(fe_logor + 1.96 * fe_se),
        "re_or": math.exp(re_logor),
        "re_or_lo": math.exp(re_logor - 1.96 * re_se),
        "re_or_hi": math.exp(re_logor + 1.96 * re_se),
        "re_logor": re_logor,
        "re_se": re_se,
        "re_p": p,
        "tau2": tau2,
        "Q": q,
        "I2": i2,
        "estimable": True,
        "cohorts": ",".join(d["cohort"].astype(str)),
    }


def fmt_or(or_, lo, hi) -> str:
    if or_ != or_:
        return "NA"
    return f"{or_:.2f} [{lo:.2f}–{hi:.2f}]"


def load_processed(proc: Path) -> dict[str, pd.DataFrame]:
    out = {}
    for p in sorted(proc.glob("*.csv")):
        if p.name in {"unusable_no_CLDN4.csv"}:
            continue
        df = pd.read_csv(p)
        if "cohort" not in df.columns:
            continue
        out[str(df["cohort"].iloc[0])] = df
    return out


def analyze_one(df: pd.DataFrame, cutoff: str, endpoint: str) -> dict:
    assay = str(df["assay"].iloc[0]) if "assay" in df.columns else ""
    work = df.copy()
    work["cldn4_log"], scale = log_expr(work["cldn4_raw"], assay)
    work["y"] = binary_from_cohort(work, endpoint)
    work = work.loc[work["y"].isin(["R", "NR"]) & work["cldn4_log"].notna()].copy()
    base = {
        "cohort": str(df["cohort"].iloc[0]),
        "cancer": str(df["cancer"].iloc[0]) if "cancer" in df.columns else "",
        "cutoff": cutoff,
        "endpoint": endpoint,
        "n_with_expr": int(df["cldn4_raw"].notna().sum()),
        "n_with_binary": int(work.shape[0]),
        "n_R": int((work["y"] == "R").sum()),
        "n_NR": int((work["y"] == "NR").sum()),
        "expr_scale": scale,
    }
    if cutoff == "continuous":
        stats = logistic_or_per_sd(work["y"], work["cldn4_log"])
        base.update(stats)
        return base
    if work.shape[0] < 6 or work["y"].nunique() < 2:
        base.update({"estimable": False, "reason": "too_few"})
        return base
    grp = assign_cutoff(work["cldn4_log"], cutoff)
    a = int(((grp == "High") & (work["y"] == "R")).sum())
    b = int(((grp == "High") & (work["y"] == "NR")).sum())
    c = int(((grp == "Low") & (work["y"] == "R")).sum())
    d = int(((grp == "Low") & (work["y"] == "NR")).sum())
    if (a + b) == 0 or (c + d) == 0:
        base.update({"estimable": False, "reason": "empty_cutoff_arm", "n_high_R": a, "n_high_NR": b, "n_low_R": c, "n_low_NR": d})
        return base
    base.update(or_from_2x2(a, b, c, d))
    return base


def forest(df: pd.DataFrame, meta: dict, title: str, out: Path) -> None:
    d = df.dropna(subset=["or", "or_lo", "or_hi"]).copy()
    if d.empty:
        return
    d = d.sort_values("cohort")
    labels = [
        f"{row['cohort']}  {int(row['n_R'])}/{int(row['n_NR'])}  "
        f"{row['or']:.2f} [{row['or_lo']:.2f}–{row['or_hi']:.2f}]"
        for _, row in d.iterrows()
    ]
    ys = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(9, 0.45 * (len(d) + 3) + 1.5))
    ax.errorbar(
        d["or"],
        ys,
        xerr=[d["or"] - d["or_lo"], d["or_hi"] - d["or"]],
        fmt="o",
        color="black",
        capsize=3,
    )
    if meta.get("estimable"):
        ax.axvline(meta["re_or"], color="tab:blue", ls="--", label=f"RE {fmt_or(meta['re_or'], meta['re_or_lo'], meta['re_or_hi'])}")
        ax.axvspan(meta["re_or_lo"], meta["re_or_hi"], color="tab:blue", alpha=0.08)
    ax.axvline(1.0, color="grey", lw=1)
    ax.axvline(USER_OR, color="tab:red", ls=":", lw=1, label=f"user 0.42 [{USER_LO}–{USER_HI}]")
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("Odds ratio (response, CLDN4-high vs low; continuous = per 1 SD)")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def write_writeup(out_dir: Path, per: pd.DataFrame, pooled: pd.DataFrame, unusable: pd.DataFrame | None) -> None:
    def row(cut, end, subset):
        m = pooled[(pooled.cutoff == cut) & (pooled.endpoint == end) & (pooled.subset == subset)]
        return m.iloc[0] if len(m) else None

    prim = row("median", "curated_R_vs_NR", "all_independent")
    lung = row("median", "curated_R_vs_NR", "lung_only")
    lines = []
    lines.append("# B5 rework: open ICI RNA, CLDN4 vs response, all cutoffs")
    lines.append("")
    lines.append("## User-reported number (not produced here)")
    lines.append("")
    lines.append(f"CLDN4-high ICI **OR = 0.42 [0.18–0.95], k=11**.")
    lines.append("Open lung ICI was previously NS. This slice recomputes the claim")
    lines.append("on public matrices and **reports every locked cutoff**. It does not")
    lines.append("choose the cutoff whose pooled OR is closest to 0.42.")
    lines.append("")
    lines.append("## Locked methods")
    lines.append("")
    lines.append("- Gene: CLDN4 (`ENSG00000189143` / symbol).")
    lines.append("- Expression: `log2(x+1)` for raw counts/TPM/FPKM. If the public matrix")
    lines.append("  already contains negative values (Zenodo ICB log-scale), it is used as-is.")
    lines.append("  GSE207422 is already log2TPM. Scale is in `tables/per_cohort_or.csv`.")
    lines.append("- Cutoffs (cohort-internal, on the analysis subset):")
    lines.append("  - median: High = ≥ median")
    lines.append("  - tertile: High = T3 vs Low = T1 (middle tertile dropped)")
    lines.append("  - quartile: High = Q4 vs Low = Q1 (Q2–Q3 dropped)")
    lines.append("  - continuous: logistic OR per 1 SD of log expression")
    lines.append("- Endpoints (all reported when the cohort has the labels):")
    lines.append("  - `curated_R_vs_NR`: ICB/GEO native R vs NR (for IMvigor210/Gide/Liu/Riaz/Braun this is CR/PR vs PD; SD is already excluded in the curated field)")
    lines.append("  - `CRPR_vs_PD`: RECIST CR/PR vs PD, SD dropped")
    lines.append("  - `CRPR_vs_SDPD`: RECIST CR/PR vs SD+PD")
    lines.append("- OR: Woolf interval; Haldane–Anscombe 0.5 if any 2×2 cell is 0; Fisher exact p.")
    lines.append("- Meta: DerSimonian–Laird random effects and inverse-variance fixed effect on logOR.")
    lines.append("- GSE135222_DCB180 is the same n as GSE135222_Jung with DCB = PFS ≥ 180 days. It is **not** an extra independent study.")
    lines.append("- GSE207422 is neoadjuvant ICI + chemotherapy, pathologic MPR vs NMPR. RECIST has no PD.")
    lines.append("- Liu raw data are dbGaP; Braun raw WES are EGA. Processed expression+response used here are the public Zenodo ICB TSVs (CC-BY-4.0).")
    lines.append("")
    lines.append("## Inventory")
    lines.append("")
    lines.append("| Cohort | Cancer | Matrix | CLDN4 | Binary n (R/NR) | Notes |")
    lines.append("|---|---|---|---|---|---|")
    notes = {
        "GSE126044": "GEO NSCLC anti-PD-1; GEO responder/non-responder",
        "GSE135222_Jung": "GEO GSE135222 / ICB_Jung; curated R/NR",
        "GSE135222_DCB180": "SAME patients as Jung; DCB=PFS≥180d; not independent",
        "GSE166449": "GEO NSCLC ICI; title responder/non-responder",
        "GSE207422": "GEO neoadjuvant ICI+chemo; MPR vs NMPR",
        "IMvigor210_Mariathasan": "Zenodo ICB_Mariathasan; atezolizumab urothelial",
        "Gide": "Zenodo ICB_Gide; melanoma anti-PD-1; n=41 public TPM",
        "Liu": "Zenodo ICB_Liu processed FPKM (raw dbGaP)",
        "Riaz": "Zenodo ICB_Riaz; melanoma nivolumab; patients with RNA",
        "Braun": "Zenodo ICB_Braun; ccRCC nivolumab; patients with RNA",
    }
    inv = per[(per.cutoff == "median") & (per.endpoint == "curated_R_vs_NR")]
    for _, r in inv.iterrows():
        lines.append(
            f"| {r.cohort} | {r.cancer} | public | yes | {r.n_with_binary} ({r.n_R}/{r.n_NR}) | {notes.get(r.cohort, '')} |"
        )
    if unusable is not None and len(unusable):
        for _, r in unusable.iterrows():
            lines.append(f"| {r.cohort} | {r.cancer} | public panel | **no** | — | {r.reason} |")
    lines.append("")
    lines.append(f"Independent studies with CLDN4 + binary response: **k = {len(PRIMARY_COHORTS)}**, not 11.")
    lines.append("")
    lines.append("## All pooled ORs (do not cherry-pick)")
    lines.append("")
    lines.append("Every cutoff × endpoint × subset is listed. The user number is")
    lines.append("shown only as a reference line. `matches_user_0.42_at_2dp` is")
    lines.append("true only when the random-effects point estimate rounds to 0.42.")
    lines.append("")
    lines.append("| Subset | Cutoff | Endpoint | k | RE OR [95% CI] | p | I² | matches 0.42 at 2dp | 0.42 in CI |")
    lines.append("|---|---|---|---:|---|---:|---:|---|---|")
    for _, r in pooled.sort_values(["subset", "endpoint", "cutoff"]).iterrows():
        if not r.estimable:
            continue
        match = abs(r.re_or - USER_OR) < 0.005
        inside = (r.re_or_lo <= USER_OR <= r.re_or_hi)
        lines.append(
            f"| {r.subset} | {r.cutoff} | {r.endpoint} | {r.k} | "
            f"{fmt_or(r.re_or, r.re_or_lo, r.re_or_hi)} | {r.re_p:.3g} | {r.I2:.0f}% | "
            f"{'YES' if match else 'no'} | {'yes' if inside else 'no'} |"
        )
    lines.append("")
    lines.append("## Primary locked comparison (median, curated R vs NR, all independent)")
    lines.append("")
    if prim is not None and prim.estimable:
        lines.append(
            f"Random-effects OR = **{fmt_or(prim.re_or, prim.re_or_lo, prim.re_or_hi)}**, "
            f"k={int(prim.k)}, p={prim.re_p:.3g}, I²={prim.I2:.0f}%."
        )
        match = abs(prim.re_or - USER_OR) < 0.005
        inside = prim.re_or_lo <= USER_OR <= prim.re_or_hi
        lines.append("")
        if match:
            lines.append("The primary point estimate rounds to 0.42. That is reported because median + R-vs-NR was locked, not because it was selected after seeing the table.")
        else:
            lines.append(
                f"The primary point estimate does **not** match the user 0.42 "
                f"(observed {prim.re_or:.2f})."
            )
        if inside:
            lines.append("The user 0.42 does fall inside the primary 95% CI (wide CI, not confirmation).")
        else:
            lines.append("The user 0.42 is outside the primary 95% CI.")
        if int(prim.k) != USER_K:
            lines.append(f"k={int(prim.k)} public independent matrices, not the claimed k=11.")
    else:
        lines.append("Primary meta was not estimable.")
    lines.append("")
    lines.append("## Open lung ICI only (median, curated R vs NR)")
    lines.append("")
    if lung is not None and lung.estimable:
        lines.append(
            f"Random-effects OR = **{fmt_or(lung.re_or, lung.re_or_lo, lung.re_or_hi)}**, "
            f"k={int(lung.k)}, p={lung.re_p:.3g}, I²={lung.I2:.0f}%."
        )
        if lung.re_p >= 0.05:
            lines.append("This lung-only pool is **not significant** (NS), consistent with the prior open-lung note.")
        else:
            lines.append("This lung-only pool is statistically significant at α=0.05 — that is reported, not hidden.")
    else:
        lines.append("Lung-only meta not estimable (small n / empty cells). Per-cohort rows are in `tables/per_cohort_or.csv`.")
    lines.append("")
    lines.append("## Per-cohort median / curated R vs NR")
    lines.append("")
    lines.append("| Cohort | n (R/NR) | High R/NR | Low R/NR | OR [95% CI] | Fisher p | Haldane |")
    lines.append("|---|---|---|---|---|---:|---|")
    sub = per[(per.cutoff == "median") & (per.endpoint == "curated_R_vs_NR") & (~per.cohort.isin(SENSITIVITY_ONLY))]
    for _, r in sub.sort_values("cohort").iterrows():
        if not r.get("estimable", False):
            lines.append(f"| {r.cohort} | {r.n_with_binary} ({r.n_R}/{r.n_NR}) | — | — | not estimable | — | — |")
            continue
        lines.append(
            f"| {r.cohort} | {r.n_used} ({r.n_R}/{r.n_NR}) | {int(r.n_high_R)}/{int(r.n_high_NR)} | "
            f"{int(r.n_low_R)}/{int(r.n_low_NR)} | {fmt_or(r['or'], r.or_lo, r.or_hi)} | {r.fisher_p:.3g} | "
            f"{'yes' if r.haldane else 'no'} |"
        )
    lines.append("")
    lines.append("## Honest verdict")
    lines.append("")
    lines.append("- The user number is a **user-reported** statistic. This slice does not treat it as ground truth.")
    lines.append("- All cutoffs are in `tables/pooled_or.csv` and `tables/per_cohort_or.csv`.")
    lines.append("- We did not drop a cohort because CLDN4 was null or opposite.")
    lines.append("- We did not add extra melanoma/HNSCC series after seeing results in order to reach k=11.")
    lines.append("- Small lung n (16–27) makes lung-only ORs unstable; tertile/quartile drop the middle and are even smaller.")
    lines.append("- Cross-cancer pooling (bladder + melanoma + RCC + lung) is not a lung-specific test.")
    lines.append("- Neoadjuvant GSE207422 is ICI+chemo pathologic response, not metastatic ORR.")
    lines.append("- Median splits that land everyone on one side (empty High or Low arm) are reported as not estimable, not recoded after the fact to force an OR.")
    lines.append("")
    lines.append("## Reproduction")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 -m pip install -r results/rework/B5_OR/requirements.txt")
    lines.append("python3 results/rework/B5_OR/download.py")
    lines.append("python3 results/rework/B5_OR/analyze.py")
    lines.append("```")
    lines.append("")
    (out_dir / "WRITEUP.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--processed", default=str(HERE / "processed"))
    ap.add_argument("--out-dir", default=str(HERE))
    args = ap.parse_args()
    proc = Path(args.processed)
    out_dir = Path(args.out_dir)
    tables = out_dir / "tables"
    figs = out_dir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    cohorts = load_processed(proc)
    unusable_path = proc / "unusable_no_CLDN4.csv"
    unusable = pd.read_csv(unusable_path) if unusable_path.exists() else None

    cutoffs = ["median", "tertile", "quartile", "continuous"]
    endpoints = ["curated_R_vs_NR", "CRPR_vs_PD", "CRPR_vs_SDPD"]
    rows = []
    for name, df in cohorts.items():
        for cutoff in cutoffs:
            for endpoint in endpoints:
                rows.append(analyze_one(df, cutoff, endpoint))
    per = pd.DataFrame(rows)
    per.to_csv(tables / "per_cohort_or.csv", index=False)

    pooled_rows = []
    subsets = {
        "all_independent": [c for c in PRIMARY_COHORTS if c in set(per.cohort)],
        "lung_only": [c for c in LUNG_COHORTS if c in set(per.cohort)],
        "non_lung": [c for c in PRIMARY_COHORTS if c not in LUNG_COHORTS and c in set(per.cohort)],
    }
    for subset_name, keep in subsets.items():
        for cutoff in cutoffs:
            for endpoint in endpoints:
                sub = per[
                    (per.cohort.isin(keep))
                    & (per.cutoff == cutoff)
                    & (per.endpoint == endpoint)
                    & (per.get("estimable") == True)  # noqa: E712
                ]
                meta = dl_meta(sub)
                meta.update({"subset": subset_name, "cutoff": cutoff, "endpoint": endpoint})
                pooled_rows.append(meta)
                if meta.get("estimable") and subset_name in {"all_independent", "lung_only"}:
                    title = f"{subset_name} | {cutoff} | {endpoint} | k={meta['k']}"
                    safe = f"forest_{subset_name}_{cutoff}_{endpoint}.png"
                    forest(sub, meta, title, figs / safe)
    pooled = pd.DataFrame(pooled_rows)
    pooled.to_csv(tables / "pooled_or.csv", index=False)

    # compact primary-facing table
    view_cols = [
        "subset",
        "cutoff",
        "endpoint",
        "k",
        "n_total",
        "re_or",
        "re_or_lo",
        "re_or_hi",
        "re_p",
        "I2",
        "fe_or",
        "fe_or_lo",
        "fe_or_hi",
    ]
    compact = pooled.loc[pooled.get("estimable") == True, [c for c in view_cols if c in pooled.columns]]  # noqa: E712
    compact = compact.copy()
    compact["matches_user_0.42_at_2dp"] = (compact["re_or"] - USER_OR).abs() < 0.005
    compact["user_0.42_in_re_ci"] = (compact["re_or_lo"] <= USER_OR) & (compact["re_or_hi"] >= USER_OR)
    compact.to_csv(tables / "pooled_or_compact.csv", index=False)

    prim = pooled[
        (pooled.subset == "all_independent")
        & (pooled.cutoff == "median")
        & (pooled.endpoint == "curated_R_vs_NR")
    ]
    prim = prim.iloc[0].to_dict() if len(prim) else {}
    lung = pooled[
        (pooled.subset == "lung_only")
        & (pooled.cutoff == "median")
        & (pooled.endpoint == "curated_R_vs_NR")
    ]
    lung = lung.iloc[0].to_dict() if len(lung) else {}

    # nearest-to-0.42 among all pooled (for honesty: show it, do not adopt it)
    est = pooled[pooled.get("estimable") == True].copy()  # noqa: E712
    if len(est):
        est["abs_diff_user"] = (est["re_or"] - USER_OR).abs()
        nearest = est.sort_values("abs_diff_user").iloc[0]
        nearest_rec = {
            "subset": nearest.subset,
            "cutoff": nearest.cutoff,
            "endpoint": nearest.endpoint,
            "re_or": float(nearest.re_or),
            "k": int(nearest.k),
            "abs_diff": float(nearest.abs_diff_user),
            "adopted": False,
            "note": "Listed only to document that we saw the closest cell and did not select it.",
        }
    else:
        nearest_rec = {}

    summary = {
        "user_claim": {"or": USER_OR, "lo": USER_LO, "hi": USER_HI, "k": USER_K},
        "primary_locked": {
            "cutoff": "median",
            "endpoint": "curated_R_vs_NR",
            "model": "DerSimonian-Laird random effects",
            "cohorts": PRIMARY_COHORTS,
            "result": {k: prim.get(k) for k in ["k", "n_total", "re_or", "re_or_lo", "re_or_hi", "re_p", "I2", "fe_or"]},
        },
        "lung_only_locked": {
            "cutoff": "median",
            "endpoint": "curated_R_vs_NR",
            "result": {k: lung.get(k) for k in ["k", "n_total", "re_or", "re_or_lo", "re_or_hi", "re_p", "I2"]},
        },
        "closest_to_user_0.42_not_adopted": nearest_rec,
        "independent_k": len(PRIMARY_COHORTS),
        "did_we_pick_cutoff_to_hit_0.42": False,
        "honest_verdict": {
            "stat_compared": "RE OR, median, curated R vs NR, all independent public matrices",
            "user_claimed_or": USER_OR,
            "observed_re_or": prim.get("re_or"),
            "observed_k": prim.get("k"),
            "matches_user_0.42_at_2dp": bool(prim.get("re_or") is not None and abs(prim.get("re_or", 99) - USER_OR) < 0.005),
            "matches_user_k11": bool(prim.get("k") == USER_K),
            "label": "COMPARE_ALL_CUTOFFS",
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=float), encoding="utf-8")
    write_writeup(out_dir, per, pooled, unusable)
    print("wrote", out_dir / "WRITEUP.md")
    print("primary", json.dumps(summary["primary_locked"]["result"], default=float))
    print("lung", json.dumps(summary["lung_only_locked"]["result"], default=float))
    print("closest_not_adopted", json.dumps(nearest_rec, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
