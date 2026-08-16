"""GSE271689 GeoMx WTA: tumour-segment TACSTD2/CLDN4 vs ICI OS and CD8-ish signatures.

Primary cohort = Yale first-line ICI (GSE271689 DCC reprocessed here; outcomes from
Aung et al. Nat Genet 2025 Source Data Fig. 6b). Replication uses the paper's
public tumour-compartment tables for UQ (CTA; no CLDN4) and Greece (WTA).

Honest constraints are written into the report, not papered over:
  - Yale n=37 / 17 OS events. Underpowered for multivariable claims.
  - Cytotoxic genes other than CD8A/CXCL9 sit near the GeoMx LOQ in tumour AOIs.
  - UQ is CTA (~1.8k genes), not WTA; CLDN4 is absent.
  - Greek source-data values are already centred/scaled (can be negative).
"""

from __future__ import annotations

import json
import os
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats

warnings.filterwarnings("ignore", category=UserWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IN = os.path.join(ROOT, "results", "hunt_gse271689", "data")
OUT = os.path.join(ROOT, "results", "hunt_gse271689")
os.makedirs(OUT, exist_ok=True)

# Pre-specified signatures. Genes missing from a matrix are dropped, not imputed.
SIGS = {
    "cd8_cytotoxic": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7"],
    "ayers_ifng": ["CXCL9", "CXCL10", "STAT1", "HLA-DRA", "HLA.DRA", "IDO1", "IFNG"],
    "tcell": ["CD3D", "CD3E", "CD2", "CD8A"],
}

FOCUS = [
    "TACSTD2",
    "CLDN4",
    "CD8A",
    "CD8B",
    "GZMA",
    "GZMB",
    "PRF1",
    "NKG7",
    "IFNG",
    "CXCL9",
    "CXCL10",
    "STAT1",
    "HLA-DRA",
    "HLA.DRA",
    "IDO1",
    "CD3D",
    "CD3E",
    "CD2",
    "EPCAM",
    "KRT19",
    "PTPRC",
]


def zscore(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    sd = s.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return s * np.nan
    return (s - s.mean()) / sd


def pick_col(df: pd.DataFrame, *names) -> str | None:
    for n in names:
        if n in df.columns:
            return n
    return None


def signature_score(df: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in df.columns]
    # collapse HLA-DRA aliases
    if "HLA-DRA" in present and "HLA.DRA" in present:
        present.remove("HLA.DRA")
    if not present:
        return pd.Series(np.nan, index=df.index), []
    z = pd.concat([zscore(df[g]) for g in present], axis=1)
    return z.mean(axis=1), present


def cox_continuous(df: pd.DataFrame, time: str, event: str, feature: str) -> dict:
    sub = df[[time, event, feature]].dropna()
    n, ev = len(sub), int(sub[event].sum())
    out = {
        "n": n,
        "events": ev,
        "hr": np.nan,
        "hr_lo": np.nan,
        "hr_hi": np.nan,
        "p": np.nan,
        "cindex": np.nan,
    }
    if n < 8 or ev < 4 or sub[feature].std(ddof=0) == 0:
        return out
    work = sub.copy()
    work[feature] = zscore(work[feature])
    cph = CoxPHFitter()
    try:
        cph.fit(work, duration_col=time, event_col=event)
        s = cph.summary.loc[feature]
        out.update(
            {
                "hr": float(s["exp(coef)"]),
                "hr_lo": float(s["exp(coef) lower 95%"]),
                "hr_hi": float(s["exp(coef) upper 95%"]),
                "p": float(s["p"]),
                "cindex": float(cph.concordance_index_),
            }
        )
    except Exception as exc:  # singular / non-convergence
        out["error"] = str(exc)
    return out


def km_median(df: pd.DataFrame, time: str, event: str, feature: str) -> dict:
    sub = df[[time, event, feature]].dropna()
    med = float(sub[feature].median())
    hi = sub[feature] > med
    out = {
        "n_hi": int(hi.sum()),
        "n_lo": int((~hi).sum()),
        "events_hi": int(sub.loc[hi, event].sum()),
        "events_lo": int(sub.loc[~hi, event].sum()),
        "median_cut": med,
        "logrank_p": np.nan,
        "hr_hi_vs_lo": np.nan,
        "hr_lo": np.nan,
        "hr_hi": np.nan,
    }
    if out["n_hi"] < 4 or out["n_lo"] < 4:
        return out
    res = logrank_test(
        sub.loc[hi, time],
        sub.loc[~hi, time],
        event_observed_A=sub.loc[hi, event],
        event_observed_B=sub.loc[~hi, event],
    )
    out["logrank_p"] = float(res.p_value)
    work = sub.copy()
    work["hi"] = hi.astype(int)
    cph = CoxPHFitter()
    try:
        cph.fit(work[[time, event, "hi"]], duration_col=time, event_col=event)
        s = cph.summary.loc["hi"]
        out["hr_hi_vs_lo"] = float(s["exp(coef)"])
        out["hr_lo"] = float(s["exp(coef) lower 95%"])
        out["hr_hi"] = float(s["exp(coef) upper 95%"])
    except Exception:
        pass
    return out


def km_plot(df, time, event, feature, title, path, xlabel="Time (days)"):
    sub = df[[time, event, feature]].dropna()
    med = sub[feature].median()
    hi = sub[feature] > med
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    km = KaplanMeierFitter()
    km.fit(sub.loc[~hi, time], sub.loc[~hi, event], label=f"low (n={int((~hi).sum())})")
    km.plot(ax=ax, ci_show=False, color="#4C78A8")
    km.fit(sub.loc[hi, time], sub.loc[hi, event], label=f"high (n={int(hi.sum())})")
    km.plot(ax=ax, ci_show=False, color="#E45756")
    lr = logrank_test(
        sub.loc[hi, time],
        sub.loc[~hi, time],
        event_observed_A=sub.loc[hi, event],
        event_observed_B=sub.loc[~hi, event],
    )
    ax.set_title(f"{title}\nmedian split, log-rank p={lr.p_value:.3g}")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Survival")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def forest_plot(rows: pd.DataFrame, path: str, title: str):
    r = rows.dropna(subset=["hr"]).copy()
    if r.empty:
        return
    r = r.iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.2, 0.42 * len(r) + 1.4))
    y = np.arange(len(r))
    ax.errorbar(
        r["hr"],
        y,
        xerr=[r["hr"] - r["hr_lo"], r["hr_hi"] - r["hr"]],
        fmt="o",
        color="#222",
        ecolor="#666",
        capsize=2,
    )
    ax.axvline(1, color="#888", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(r["label"].tolist(), fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("HR per +1 SD (Cox)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def load_yale() -> pd.DataFrame:
    df = pd.read_csv(os.path.join(IN, "yale_tumour_patient_level.csv.gz"))
    df["cohort"] = "Yale"
    df["time_os"] = df["OS_Days"]
    df["event_os"] = df["OS_Index"].astype(int)
    df["time_pfs"] = df["PFS_Days"]
    df["event_pfs"] = df["PFS_Index"].astype(int)
    df["histology"] = df["Histotype"]
    df["line"] = "1st"
    return df


def load_uq() -> pd.DataFrame:
    df = pd.read_csv(os.path.join(IN, "uq_tumour_patient_level.csv.gz"))
    df["cohort"] = "UQ"
    # Paper OS for UQ is 2-year censored (OS_Days_2Yrs_byIT).
    df["time_os"] = df["OS_Days_2Yrs_byIT"]
    df["event_os"] = df["OS_Index_2Yrs_byIT"].astype(int)
    df["time_pfs"] = df["PFS_Days_byIT"]
    df["event_pfs"] = df["PFS_Index_byIT"].astype(int)
    df["line"] = df["Line.of.therapy.Itx"].astype(str)
    df["histology"] = np.nan
    return df


def load_greek() -> pd.DataFrame:
    df = pd.read_csv(os.path.join(IN, "greek_tumour_roi_level.csv.gz"), low_memory=False)
    df["cohort"] = "Greek"
    df["time_os"] = pd.to_numeric(df["OS_1"], errors="coerce")
    df["event_os"] = df["Death"].map({"Yes": 1, "No": 0}).astype(int)
    df["time_pfs"] = pd.to_numeric(df["PFS_1"], errors="coerce")
    # PFS_2Years_index is the paper's 2-year PFS event flag
    df["event_pfs"] = pd.to_numeric(df["PFS_2Years_index"], errors="coerce").fillna(0).astype(int)
    df["histology"] = df["Histotype"]
    df["line"] = np.nan
    return df


def add_scores(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    used = {}
    for name, genes in SIGS.items():
        score, present = signature_score(out, genes)
        out[name] = score
        used[name] = present
    out["_sig_genes"] = json.dumps(used)
    return out


def analyse_block(df: pd.DataFrame, cohort: str, features: list[str]) -> list[dict]:
    rows = []
    for feat in features:
        if feat not in df.columns or df[feat].notna().sum() < 8:
            continue
        for endpoint, t, e in [("OS", "time_os", "event_os"), ("PFS", "time_pfs", "event_pfs")]:
            cont = cox_continuous(df, t, e, feat)
            med = km_median(df, t, e, feat)
            rows.append(
                {
                    "cohort": cohort,
                    "n": cont["n"],
                    "events": cont["events"],
                    "endpoint": endpoint,
                    "feature": feat,
                    "model": "continuous_z",
                    **{k: cont[k] for k in ["hr", "hr_lo", "hr_hi", "p", "cindex"]},
                    "logrank_p_median": med["logrank_p"],
                    "hr_median_hi_vs_lo": med["hr_hi_vs_lo"],
                    "n_hi": med["n_hi"],
                    "n_lo": med["n_lo"],
                }
            )
    return rows


def scatter_plot(df, x, y, title, path):
    a = pd.to_numeric(df[x], errors="coerce")
    b = pd.to_numeric(df[y], errors="coerce")
    m = a.notna() & b.notna()
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    ax.scatter(a[m], b[m], s=28, c="#4C78A8", alpha=0.85, edgecolors="none")
    if m.sum() >= 5:
        rho, p = stats.spearmanr(a[m], b[m])
        ax.set_title(f"{title}\nSpearman ρ={rho:.2f}, p={p:.3g}")
    else:
        ax.set_title(title)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    yale = add_scores(load_yale())
    uq = add_scores(load_uq())
    greek = add_scores(load_greek())

    features = ["TACSTD2", "CLDN4", "CD8A", "cd8_cytotoxic", "ayers_ifng", "tcell", "CXCL9"]
    rows = []
    rows += analyse_block(yale, "Yale", features)
    rows += analyse_block(yale[yale["histology"] == "Adenocarcinoma"], "Yale_LUAD", features)
    rows += analyse_block(uq, "UQ", features)
    rows += analyse_block(uq[uq["line"].str.contains("1st", case=False, na=False)], "UQ_1st", features)
    rows += analyse_block(greek, "Greek", features)
    rows += analyse_block(greek[greek["histology"] == "Adenocarcinoma"], "Greek_LUAD", features)
    cox_tbl = pd.DataFrame(rows)
    cox_tbl.to_csv(os.path.join(OUT, "cox_os_pfs.csv"), index=False)

    # correlations
    corr_rows = []
    for name, df in [("Yale", yale), ("UQ", uq), ("Greek", greek)]:
        pairs = [
            ("TACSTD2", "CLDN4"),
            ("TACSTD2", "CD8A"),
            ("TACSTD2", "cd8_cytotoxic"),
            ("TACSTD2", "ayers_ifng"),
            ("TACSTD2", "CXCL9"),
            ("CLDN4", "CD8A"),
            ("CLDN4", "cd8_cytotoxic"),
            ("CLDN4", "ayers_ifng"),
        ]
        for a, b in pairs:
            if a not in df.columns or b not in df.columns:
                continue
            x = pd.to_numeric(df[a], errors="coerce")
            y = pd.to_numeric(df[b], errors="coerce")
            m = x.notna() & y.notna()
            if m.sum() < 8:
                continue
            rho, p = stats.spearmanr(x[m], y[m])
            corr_rows.append(
                {"cohort": name, "x": a, "y": b, "n": int(m.sum()), "spearman_rho": float(rho), "p": float(p)}
            )
    corr_tbl = pd.DataFrame(corr_rows)
    corr_tbl.to_csv(os.path.join(OUT, "correlations.csv"), index=False)

    # multivariable: TACSTD2 + CD8A + histology (Yale only; n is small)
    mv_rows = []
    y = yale[["time_os", "event_os", "TACSTD2", "CD8A", "histology"]].dropna()
    y["TACSTD2_z"] = zscore(y["TACSTD2"])
    y["CD8A_z"] = zscore(y["CD8A"])
    y["adeno"] = (y["histology"] == "Adenocarcinoma").astype(int)
    if len(y) >= 20 and y["event_os"].sum() >= 8:
        cph = CoxPHFitter()
        try:
            cph.fit(y[["time_os", "event_os", "TACSTD2_z", "CD8A_z", "adeno"]], "time_os", "event_os")
            for feat in ["TACSTD2_z", "CD8A_z", "adeno"]:
                s = cph.summary.loc[feat]
                mv_rows.append(
                    {
                        "cohort": "Yale",
                        "model": "OS ~ TACSTD2_z + CD8A_z + adeno",
                        "feature": feat,
                        "n": len(y),
                        "events": int(y["event_os"].sum()),
                        "hr": float(s["exp(coef)"]),
                        "hr_lo": float(s["exp(coef) lower 95%"]),
                        "hr_hi": float(s["exp(coef) upper 95%"]),
                        "p": float(s["p"]),
                    }
                )
        except Exception as exc:
            mv_rows.append({"cohort": "Yale", "error": str(exc)})
    pd.DataFrame(mv_rows).to_csv(os.path.join(OUT, "cox_multivariable_yale.csv"), index=False)

    # compact patient tables (no 18k-gene dump)
    def compact(df, extra):
        keep = ["cohort", "time_os", "event_os", "time_pfs", "event_pfs", "histology", "line"]
        keep += [c for c in extra if c in df.columns]
        keep += [c for c in ["TACSTD2", "CLDN4", "CD8A", "CXCL9", "cd8_cytotoxic", "ayers_ifng", "tcell"] if c in df.columns]
        return df[keep].copy()

    compact(yale, ["spot_id", "n_tumour_aoi", "Agent_ITx", "Stage_at_ITx"]).to_csv(
        os.path.join(OUT, "yale_patient_compact.csv"), index=False
    )
    compact(uq, ["new_ID"]).to_csv(os.path.join(OUT, "uq_patient_compact.csv"), index=False)
    compact(greek, ["ROILabel", "Histotype"]).to_csv(os.path.join(OUT, "greek_roi_compact.csv"), index=False)

    # plots
    km_plot(yale, "time_os", "event_os", "TACSTD2", "Yale tumour TACSTD2 — OS", os.path.join(OUT, "km_yale_os_tacstd2.png"))
    km_plot(yale, "time_os", "event_os", "CLDN4", "Yale tumour CLDN4 — OS", os.path.join(OUT, "km_yale_os_cldn4.png"))
    km_plot(yale, "time_os", "event_os", "CD8A", "Yale tumour CD8A — OS", os.path.join(OUT, "km_yale_os_cd8a.png"))
    km_plot(yale, "time_pfs", "event_pfs", "TACSTD2", "Yale tumour TACSTD2 — PFS", os.path.join(OUT, "km_yale_pfs_tacstd2.png"))
    km_plot(uq, "time_os", "event_os", "TACSTD2", "UQ tumour TACSTD2 — OS (2y-censored)", os.path.join(OUT, "km_uq_os_tacstd2.png"))
    km_plot(greek, "time_os", "event_os", "TACSTD2", "Greek tumour TACSTD2 — OS", os.path.join(OUT, "km_greek_os_tacstd2.png"))
    km_plot(greek, "time_os", "event_os", "CLDN4", "Greek tumour CLDN4 — OS", os.path.join(OUT, "km_greek_os_cldn4.png"))

    scatter_plot(yale, "TACSTD2", "CLDN4", "Yale tumour", os.path.join(OUT, "scatter_yale_tacstd2_cldn4.png"))
    scatter_plot(yale, "TACSTD2", "CD8A", "Yale tumour", os.path.join(OUT, "scatter_yale_tacstd2_cd8a.png"))
    scatter_plot(yale, "TACSTD2", "cd8_cytotoxic", "Yale tumour", os.path.join(OUT, "scatter_yale_tacstd2_cd8sig.png"))
    scatter_plot(yale, "TACSTD2", "ayers_ifng", "Yale tumour", os.path.join(OUT, "scatter_yale_tacstd2_ifng.png"))

    os_cont = cox_tbl[(cox_tbl["endpoint"] == "OS") & (cox_tbl["model"] == "continuous_z")].copy()
    os_cont["label"] = os_cont.apply(lambda r: f"{r['cohort']}  {r['feature']}  (n={int(r['n'])}, ev={int(r['events'])})", axis=1)
    forest_plot(
        os_cont[os_cont["feature"].isin(["TACSTD2", "CLDN4", "CD8A", "cd8_cytotoxic"])],
        os.path.join(OUT, "forest_os_continuous.png"),
        "OS Cox HR per +1 SD (pre-specified features)",
    )

    # stats.txt
    def fmt(row):
        if row is None or (isinstance(row, float) and np.isnan(row)):
            return "NA"
        return (
            f"HR={row['hr']:.3f} [{row['hr_lo']:.3f}-{row['hr_hi']:.3f}] p={row['p']:.3g} "
            f"n={int(row['n'])} ev={int(row['events'])} c={row['cindex']:.3f}"
        )

    def grab(tbl, cohort, feat, endpoint="OS"):
        hit = tbl[(tbl.cohort == cohort) & (tbl.feature == feat) & (tbl.endpoint == endpoint)]
        return None if hit.empty else hit.iloc[0]

    yale_sig = json.loads(yale["_sig_genes"].iloc[0])
    uq_sig = json.loads(uq["_sig_genes"].iloc[0])
    greek_sig = json.loads(greek["_sig_genes"].iloc[0])

    lines = []
    lines.append("GSE271689 GeoMx WTA tumour-segment TACSTD2/CLDN4 vs ICI OS")
    lines.append("")
    lines.append("PRIMARY Yale (GSE271689 DCC reprocessed; Q3-log2; mean of tumour AOIs)")
    lines.append(f"  n={len(yale)} OS events={int(yale.event_os.sum())} PFS events={int(yale.event_pfs.sum())}")
    lines.append(f"  all first-line ICI, pre-treatment biopsy; adeno={int((yale.histology=='Adenocarcinoma').sum())} sq={int((yale.histology=='Squamous').sum())}")
    lines.append(f"  TACSTD2 OS  {fmt(grab(cox_tbl,'Yale','TACSTD2'))}")
    lines.append(f"  CLDN4   OS  {fmt(grab(cox_tbl,'Yale','CLDN4'))}")
    lines.append(f"  CD8A    OS  {fmt(grab(cox_tbl,'Yale','CD8A'))}")
    lines.append(f"  cd8sig  OS  {fmt(grab(cox_tbl,'Yale','cd8_cytotoxic'))}")
    lines.append(f"  ayers   OS  {fmt(grab(cox_tbl,'Yale','ayers_ifng'))}")
    lines.append(f"  TACSTD2 PFS {fmt(grab(cox_tbl,'Yale','TACSTD2','PFS'))}")
    lines.append(f"  CLDN4   PFS {fmt(grab(cox_tbl,'Yale','CLDN4','PFS'))}")
    ymed = km_median(yale, "time_os", "event_os", "TACSTD2")
    lines.append(
        f"  TACSTD2 OS median-split HR_hi/lo={ymed['hr_hi_vs_lo']:.3f} logrank_p={ymed['logrank_p']:.3g} "
        f"n_hi={ymed['n_hi']} n_lo={ymed['n_lo']}"
    )
    lines.append(f"  Yale LUAD-only TACSTD2 OS {fmt(grab(cox_tbl,'Yale_LUAD','TACSTD2'))}")
    lines.append("")
    lines.append("CORRELATIONS Yale")
    for _, r in corr_tbl[corr_tbl.cohort == "Yale"].iterrows():
        lines.append(f"  {r['x']} vs {r['y']}: rho={r['spearman_rho']:.3f} p={r['p']:.3g} n={int(r['n'])}")
    lines.append(f"  signature genes used: {yale_sig}")
    lines.append("")
    lines.append("REPLICATION UQ (paper Source Data Fig.6c; CTA panel; CLDN4 absent)")
    lines.append(f"  n={len(uq)} OS events={int(uq.event_os.sum())}  1st-line n={int(uq.line.str.contains('1st', case=False, na=False).sum())}")
    lines.append(f"  TACSTD2 OS  {fmt(grab(cox_tbl,'UQ','TACSTD2'))}")
    lines.append(f"  CD8A    OS  {fmt(grab(cox_tbl,'UQ','CD8A'))}")
    lines.append(f"  cd8sig  OS  {fmt(grab(cox_tbl,'UQ','cd8_cytotoxic'))}")
    lines.append(f"  UQ 1st-line TACSTD2 OS {fmt(grab(cox_tbl,'UQ_1st','TACSTD2'))}")
    for _, r in corr_tbl[corr_tbl.cohort == "UQ"].iterrows():
        lines.append(f"  {r['x']} vs {r['y']}: rho={r['spearman_rho']:.3f} p={r['p']:.3g}")
    lines.append(f"  signature genes used: {uq_sig}")
    lines.append("")
    lines.append("REPLICATION Greek (paper Source Data Fig.6d; WTA tumour ROI ≈ patient)")
    lines.append(f"  n={len(greek)} OS events={int(greek.event_os.sum())}")
    lines.append(f"  TACSTD2 OS  {fmt(grab(cox_tbl,'Greek','TACSTD2'))}")
    lines.append(f"  CLDN4   OS  {fmt(grab(cox_tbl,'Greek','CLDN4'))}")
    lines.append(f"  CD8A    OS  {fmt(grab(cox_tbl,'Greek','CD8A'))}")
    lines.append(f"  Greek LUAD TACSTD2 OS {fmt(grab(cox_tbl,'Greek_LUAD','TACSTD2'))}")
    for _, r in corr_tbl[corr_tbl.cohort == "Greek"].iterrows():
        lines.append(f"  {r['x']} vs {r['y']}: rho={r['spearman_rho']:.3f} p={r['p']:.3g}")
    lines.append(f"  signature genes used: {greek_sig}")
    lines.append("")
    if mv_rows and "hr" in mv_rows[0]:
        lines.append("YALE multivariable OS ~ TACSTD2_z + CD8A_z + adeno")
        for r in mv_rows:
            lines.append(f"  {r['feature']}: HR={r['hr']:.3f} [{r['hr_lo']:.3f}-{r['hr_hi']:.3f}] p={r['p']:.3g}")
    lines.append("")
    lines.append("DIRECTION vs user prior (high TROP2 worse):")
    for cohort, feat in [("Yale", "TACSTD2"), ("Yale", "CLDN4"), ("UQ", "TACSTD2"), ("Greek", "TACSTD2"), ("Greek", "CLDN4")]:
        r = grab(cox_tbl, cohort, feat)
        if r is None or not np.isfinite(r["hr"]):
            lines.append(f"  {cohort} {feat}: not estimable")
        else:
            direction = "WORSE (HR>1)" if r["hr"] > 1 else "BETTER (HR<1)"
            sig = "nominal p<0.05" if r["p"] < 0.05 else "not significant"
            lines.append(f"  {cohort} {feat}: {direction}, {sig}, {fmt(r)}")
    open(os.path.join(OUT, "stats.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
