#!/usr/bin/env python3
"""TACSTD2 / CLDN4 vs ICI response in the selected lung datasets.

Datasets (all downloaded by 05_download.py):
  - GSE261345 (ES-SCLC, GeoMx DSP CTA)  : TACSTD2 vs RECIST + PFS
  - GSE261348 (ES-SCLC, GeoMx DSP CTA)  : TACSTD2 vs RECIST + PFS
  - GSE233203 (NSCLC, 10x scRNA pseudobulk): TACSTD2 & CLDN4 vs response

CLDN4 is NOT on the GeoMx CTA panel, so it is only assessable in GSE233203.

Outputs -> results/w200/GEO_2026/:
  analysis_summary.json, analysis_results.tsv,
  per-dataset value tables, and boxplot PNGs.
"""
import gzip
import io
import json
import math
import tarfile
import urllib.request
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RES = Path(__file__).resolve().parents[2] / "results" / "w200" / "GEO_2026"
DATA = RES / "data"
GENES = ["TACSTD2", "CLDN4"]
RESULTS = []  # collected rows for analysis_results.tsv


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def full_matrix(gse: str) -> str:
    stub = gse[:-3] + "nnn"
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/matrix/"
           f"{gse}_series_matrix.txt.gz")
    raw = urllib.request.urlopen(
        urllib.request.Request(url, headers={"User-Agent": "geo"}),
        timeout=90).read()
    return gzip.GzipFile(fileobj=io.BytesIO(raw)).read().decode(
        "utf-8", "replace")


def parse_series(gse: str) -> pd.DataFrame:
    """Return per-GSM dataframe: geo_accession, title, description, and one
    column per characteristic key."""
    txt = full_matrix(gse)
    fields = {}

    def row(prefix):
        for ln in txt.splitlines():
            if ln.startswith(prefix):
                return [x.strip().strip('"') for x in ln.split("\t")[1:]]
        return None

    fields["geo_accession"] = row("!Sample_geo_accession")
    fields["title"] = row("!Sample_title")
    fields["description"] = row("!Sample_description")
    n = len(fields["geo_accession"])
    # Parse per-CELL: each characteristic cell is its own "key: value" pair,
    # so we assign by the cell's own key. This is robust to GEO's common
    # column-misalignment where different samples list characteristics in a
    # different order.
    per_sample = [dict() for _ in range(n)]
    for ln in txt.splitlines():
        if ln.startswith("!Sample_characteristics"):
            vals = [x.strip().strip('"') for x in ln.split("\t")[1:]]
            for i, v in enumerate(vals):
                if i < n and ":" in v:
                    k, val = v.split(":", 1)
                    k = k.strip().lower()
                    val = val.strip()
                    if val and (k not in per_sample[i] or not per_sample[i][k]):
                        per_sample[i][k] = val
    all_keys = sorted({k for d in per_sample for k in d})
    df = pd.DataFrame({"geo_accession": fields["geo_accession"],
                       "title": fields["title"],
                       "description": fields["description"]})
    for k in all_keys:
        df[k] = [d.get(k, "") for d in per_sample]
    return df


def mannwhitney(a, b, method="auto"):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan")
    u, p = stats.mannwhitneyu(
        a, b, alternative="two-sided", method=method
    )
    return u, p


def cliffs_delta(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    gt = sum(1 for x in a for y in b if x > y)
    lt = sum(1 for x in a for y in b if x < y)
    return (gt - lt) / (len(a) * len(b))


def boxplot(groups, labels, title, ylab, path):
    fig, ax = plt.subplots(figsize=(4.2, 4))
    data = [g for g in groups]
    bp = ax.boxplot(data, tick_labels=labels, showfliers=False,
                    patch_artist=True, widths=0.6)
    for patch, c in zip(bp["boxes"], ["#4C9F70", "#C0504D", "#4472C4"]):
        patch.set_facecolor(c)
        patch.set_alpha(0.45)
    rng = np.random.default_rng(0)
    for i, g in enumerate(data, 1):
        x = rng.normal(i, 0.05, size=len(g))
        ax.scatter(x, g, s=18, color="black", zorder=3, alpha=0.7)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylab)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def marker_scatter(df, path):
    fig, ax = plt.subplots(figsize=(4.2, 4))
    for label, color in (
        ("Response", "#4C9F70"), ("Non-response", "#C0504D")
    ):
        group = df[df["response"] == label]
        ax.scatter(
            group["CLDN4"], group["TACSTD2"],
            color=color, label=label, s=38, alpha=0.85,
        )
    for _, row in df.iterrows():
        ax.annotate(
            row["sample_file"].split("_", 1)[1],
            (row["CLDN4"], row["TACSTD2"]),
            xytext=(3, 3), textcoords="offset points", fontsize=6,
        )
    ax.set_xlabel("CLDN4 log2(pseudobulk CPM+1)")
    ax.set_ylabel("TACSTD2 log2(pseudobulk CPM+1)")
    ax.set_title("GSE233203: whole-sample marker concordance", fontsize=10)
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------
# GSE261345 / GSE261348  (GeoMx DSP, ES-SCLC)
# ---------------------------------------------------------------------------
RESP_MAP = {
    "complete response": "R", "partial response": "R",
    "stable disease": "NR", "progressive disease": "NR",
}


def parse_date(s):
    s = (s or "").strip()
    for fmt in ("%m/%d/%Y", "%Y-%m-%d", "%d/%m/%Y"):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            continue
    return None


def analyze_dsp(gse: str, xlsx: str):
    print("\n" + "=" * 70, "\n", gse)
    meta = parse_series(gse)
    # normalized target counts: genes x segments (SegmentDisplayName cols)
    tcm = pd.read_excel(DATA / gse / xlsx, sheet_name="TargetCountMatrix")
    tcm = tcm.set_index(tcm.columns[0])
    tcm.index = tcm.index.astype(str).str.upper()
    seg_cols = list(tcm.columns)

    recist_col = "best recist response to treatment"
    prog_col = "disease progression or death"
    pid_col = "patient id"

    # per-segment table keyed by SegmentDisplayName == series 'description'
    rows = []
    for _, r in meta.iterrows():
        seg = r["description"]
        if seg not in tcm.columns:
            continue
        val = float(tcm.loc["TACSTD2", seg]) if "TACSTD2" in tcm.index else np.nan
        rows.append({
            "segment": seg, "patient": r.get(pid_col, ""),
            "recist": r.get(recist_col, "").lower(),
            "prog_event": r.get(prog_col, ""),
            "first_dose": r.get("date of first dose of treatment", ""),
            "prog_date": r.get("date of disease progression or death", ""),
            "last_fu": r.get("date of last follow-up", ""),
            "tissue": r.get("tissue", ""),
            "TACSTD2_norm": val,
        })
    seg_df = pd.DataFrame(rows)
    seg_df["TACSTD2_log2"] = np.log2(seg_df["TACSTD2_norm"] + 1)
    matched = seg_df["segment"].notna().sum()
    print(f"  segments matched to metadata: {matched}/{len(seg_cols)}")

    # aggregate to patient level (mean across that patient's ROIs)
    pat = seg_df.groupby("patient").agg(
        TACSTD2_log2=("TACSTD2_log2", "mean"),
        n_roi=("segment", "count"),
        recist=("recist", "first"),
        prog_event=("prog_event", "first"),
        first_dose=("first_dose", "first"),
        prog_date=("prog_date", "first"),
        last_fu=("last_fu", "first"),
    ).reset_index()
    pat["resp_group"] = pat["recist"].map(RESP_MAP)
    print("  RECIST value counts:", seg_df["recist"].value_counts().to_dict())
    print(f"  patients: {len(pat)}  (R={sum(pat.resp_group=='R')}, "
          f"NR={sum(pat.resp_group=='NR')})")

    pat.to_csv(RES / f"{gse}_patient_TACSTD2.tsv", sep="\t", index=False)

    # RECIST responder vs non-responder (patient level)
    r = pat.loc[pat.resp_group == "R", "TACSTD2_log2"].dropna()
    nr = pat.loc[pat.resp_group == "NR", "TACSTD2_log2"].dropna()
    u, p = mannwhitney(r, nr)
    d = cliffs_delta(r, nr)
    print(f"  TACSTD2 R vs NR: median R={r.median():.3f} (n={len(r)}), "
          f"NR={nr.median():.3f} (n={len(nr)}), MWU p={p:.3f}, cliff={d:.3f}")
    RESULTS.append({
        "dataset": gse, "cohort": "ES-SCLC (DSP)", "gene": "TACSTD2",
        "comparison": "RECIST responder(CR/PR) vs non(SD/PD), patient-level",
        "test": "Mann-Whitney U (two-sided)",
        "n_group1": len(r), "n_group2": len(nr),
        "median_group1": round(float(r.median()), 4) if len(r) else None,
        "median_group2": round(float(nr.median()), 4) if len(nr) else None,
        "statistic": round(float(u), 4) if not math.isnan(u) else None,
        "p_value": round(float(p), 4) if not math.isnan(p) else None,
        "effect_type": "Cliff's delta (R minus NR)",
        "effect_value": round(float(d), 4) if not math.isnan(d) else None,
    })
    if len(r) >= 2 and len(nr) >= 2:
        boxplot([r.values, nr.values], [f"R\n(n={len(r)})", f"NR\n(n={len(nr)})"],
                f"{gse}: TACSTD2 vs RECIST", "log2(Q3-norm counts+1)",
                RES / f"{gse}_TACSTD2_RECIST.png")

    # PFS association (event = progression or death Yes)
    def pfs_days(row):
        d0 = parse_date(row["first_dose"])
        if d0 is None:
            return None, None
        event = str(row["prog_event"]).strip().lower() in ("yes", "1", "true")
        if event:
            d1 = parse_date(row["prog_date"])
        else:
            d1 = parse_date(row["last_fu"])
        if d1 is None:
            return None, event
        return (d1 - d0).days, event

    pfs = pat.apply(lambda x: pfs_days(x), axis=1)
    pat["pfs_days"] = [x[0] for x in pfs]
    pat["pfs_event"] = [x[1] for x in pfs]
    surv = pat.dropna(subset=["pfs_days", "TACSTD2_log2"])
    if len(surv) >= 6:
        rho, prho = stats.spearmanr(surv["TACSTD2_log2"], surv["pfs_days"])
        print(f"  Spearman(TACSTD2, PFS days) = {rho:.3f} (p={prho:.3f}, "
              f"n={len(surv)})")
        RESULTS.append({
            "dataset": gse, "cohort": "ES-SCLC (DSP)", "gene": "TACSTD2",
            "comparison": "Spearman TACSTD2 vs PFS days (patient-level)",
            "test": "Spearman rank correlation (censoring ignored)",
            "n_group1": len(surv), "n_group2": None,
            "median_group1": None, "median_group2": None,
            "statistic": round(float(rho), 4),
            "p_value": round(float(prho), 4),
            "effect_type": "Spearman rho",
            "effect_value": round(float(rho), 4),
        })
        # median-split logrank if lifelines present
        try:
            from lifelines.statistics import logrank_test
            med = surv["TACSTD2_log2"].median()
            hi = surv[surv.TACSTD2_log2 > med]
            lo = surv[surv.TACSTD2_log2 <= med]
            lr = logrank_test(hi.pfs_days, lo.pfs_days,
                              event_observed_A=hi.pfs_event.astype(int),
                              event_observed_B=lo.pfs_event.astype(int))
            print(f"  PFS logrank TACSTD2 hi vs lo: p={lr.p_value:.3f}")
            RESULTS.append({
                "dataset": gse, "cohort": "ES-SCLC (DSP)", "gene": "TACSTD2",
                "comparison": "PFS logrank TACSTD2 hi vs lo (median split)",
                "test": "Log-rank (median split)",
                "n_group1": len(hi), "n_group2": len(lo),
                "median_group1": None, "median_group2": None,
                "statistic": round(float(lr.test_statistic), 4),
                "p_value": round(float(lr.p_value), 4),
                "effect_type": None, "effect_value": None,
            })
        except Exception as e:  # noqa: BLE001
            print("  (lifelines unavailable, skipping logrank:", e, ")")
    pat.to_csv(RES / f"{gse}_patient_TACSTD2.tsv", sep="\t", index=False)


# ---------------------------------------------------------------------------
# GSE233203  (NSCLC 10x scRNA pseudobulk)
# ---------------------------------------------------------------------------
def analyze_scrna(gse: str):
    print("\n" + "=" * 70, "\n", gse)
    from scipy.io import mmread
    meta = parse_series(gse)
    resp_col = [c for c in meta.columns if "therapeutic response" in c][0]
    meta["sample"] = meta["title"]
    resp = dict(zip(meta["title"], meta[resp_col]))

    tar = DATA / gse / "GSE233203_RAW.tar"
    rows = []
    with tarfile.open(tar) as t:
        names = t.getnames()
        samples = sorted({n.split("_matrix.mtx")[0] for n in names
                          if n.endswith("matrix.mtx.gz")})
        for s in samples:
            feats = gzip.decompress(
                t.extractfile(s + "_features.tsv.gz").read()).decode().splitlines()
            syms = [ln.split("\t")[1] if "\t" in ln else ln for ln in feats]
            mtx = mmread(io.BytesIO(gzip.decompress(
                t.extractfile(s + "_matrix.mtx.gz").read())), spmatrix=True)
            mtx = mtx.tocsr()  # genes x cells
            total = mtx.sum()
            gene_sum = np.asarray(mtx.sum(axis=1)).ravel()
            title = s.split("_", 1)[1]  # e.g. NCCLu_162
            # GSM prefix in filename -> match to meta by title suffix
            row = {"sample_file": s, "n_cells": mtx.shape[1],
                   "total_counts": float(total)}
            sym_idx = {sy.upper(): i for i, sy in enumerate(syms)}
            for g in GENES:
                if g in sym_idx:
                    cpm = gene_sum[sym_idx[g]] / total * 1e6
                    row[g] = math.log2(cpm + 1)
                else:
                    row[g] = float("nan")
            # response from meta (match by the NCCLu_xxx token in title)
            token = title
            match = [v for k, v in resp.items() if token in k or k in token]
            row["response"] = match[0] if match else ""
            rows.append(row)
            print(f"  {s}: cells={mtx.shape[1]}, resp={row['response']}, "
                  f"TACSTD2={row['TACSTD2']:.2f}, CLDN4={row['CLDN4']:.2f}")
    df = pd.DataFrame(rows)
    df.to_csv(RES / f"{gse}_pseudobulk.tsv", sep="\t", index=False)

    for g in GENES:
        r = df.loc[df.response == "Response", g].dropna()
        nr = df.loc[df.response == "Non-response", g].dropna()
        u, p = mannwhitney(r, nr, method="exact")
        d = cliffs_delta(r, nr)
        print(f"  {g}: R median={r.median():.3f} (n={len(r)}), "
              f"NR median={nr.median():.3f} (n={len(nr)}), MWU p={p:.3f}, "
              f"cliff={d:.3f}")
        RESULTS.append({
            "dataset": gse,
            "cohort": "NSCLC (whole-sample scRNA pseudobulk)", "gene": g,
            "comparison": "Response vs Non-response, one whole-sample pseudobulk per patient",
            "test": "Mann-Whitney U (two-sided, exact)",
            "n_group1": len(r), "n_group2": len(nr),
            "median_group1": round(float(r.median()), 4) if len(r) else None,
            "median_group2": round(float(nr.median()), 4) if len(nr) else None,
            "statistic": round(float(u), 4) if not math.isnan(u) else None,
            "p_value": round(float(p), 4) if not math.isnan(p) else None,
            "effect_type": "Cliff's delta (R minus NR)",
            "effect_value": round(float(d), 4) if not math.isnan(d) else None,
        })
        if len(r) >= 2 and len(nr) >= 2:
            boxplot([r.values, nr.values],
                    [f"R\n(n={len(r)})", f"NR\n(n={len(nr)})"],
                    f"{gse}: {g} vs response", "log2(pseudobulk CPM+1)",
                    RES / f"{gse}_{g}_response.png")

    rho, p = stats.spearmanr(df["TACSTD2"], df["CLDN4"])
    RESULTS.append({
        "dataset": gse,
        "cohort": "NSCLC (whole-sample scRNA pseudobulk)",
        "gene": "TACSTD2~CLDN4",
        "comparison": "Marker concordance across patients",
        "test": "Spearman rank correlation",
        "n_group1": len(df), "n_group2": None,
        "median_group1": None, "median_group2": None,
        "statistic": round(float(rho), 4),
        "p_value": round(float(p), 4),
        "effect_type": "Spearman rho",
        "effect_value": round(float(rho), 4),
    })
    marker_scatter(df, RES / f"{gse}_TACSTD2_CLDN4_scatter.png")


def main():
    analyze_dsp("GSE261345", "GSE261345_CANTABRICO_DSP_normalizedcounts.xlsx")
    analyze_dsp("GSE261348", "GSE261348_IMfirst_DSP_normalizedcounts.xlsx")
    analyze_scrna("GSE233203")

    out = pd.DataFrame(RESULTS)
    out.to_csv(RES / "analysis_results.tsv", sep="\t", index=False)
    (RES / "analysis_summary.json").write_text(json.dumps(RESULTS, indent=2))
    print("\nWrote", RES / "analysis_results.tsv")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
