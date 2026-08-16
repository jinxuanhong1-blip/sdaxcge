#!/usr/bin/env python3
"""
TACSTD2 (TROP2) and CLDN4 vs Durable Clinical Benefit (DCB)
in anti-PD1 treated advanced NSCLC (GSE190265 / GSE190266).

Honest constraints baked in:
  - DCB cutoff is the conventional 6-month PFS rule, not tuned to these genes
  - two-sided tests only; no one-sided fishing
  - platforms differ (HiSeq/NextSeq vs NovaSeq); no naive TPM merge
  - GSE190266 GEO TPM matrix is truncated at 16,384 columns (ends at MTMR14),
    so TACSTD2 is unavailable in that cohort
  - GSE190266 PFS is deposited already capped at 6 months
"""
from __future__ import annotations

import gzip
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, "data")
OUTDIR = os.path.join(REPO, "results", "w200", "GSE190265")
os.makedirs(OUTDIR, exist_ok=True)

GENES = ["TACSTD2", "CLDN4"]
DCB_MONTHS = 6.0


def classify_dcb(time_pfs, evt, cutoff=DCB_MONTHS):
    if pd.isna(time_pfs):
        return None
    if float(time_pfs) >= cutoff:
        return "DCB"
    if evt == 1:
        return "NDB"
    return None


def parse_decimal(x):
    if pd.isna(x):
        return np.nan
    if isinstance(x, (int, float, np.floating, np.integer)):
        return float(x)
    s = str(x).strip()
    if s == "" or s.lower() in {"na", "nan", "none"}:
        return np.nan
    if "," in s and "." not in s:
        s = s.replace(",", ".")
    try:
        return float(s)
    except ValueError:
        return np.nan


def load_clinical_265():
    p = os.path.join(DATA, "GSE190265", "GSE190265_samples_info_France3.csv.gz")
    with gzip.open(p, "rt") as fh:
        df = pd.read_csv(fh, sep=";")
    df["sample"] = df["sample"].astype(str)
    df["pfs_time"] = pd.to_numeric(df["time_PFS"], errors="coerce")
    df["pfs_evt"] = pd.to_numeric(df["evtPFS"], errors="coerce")
    df["DCB"] = [classify_dcb(t, e) for t, e in zip(df["pfs_time"], df["pfs_evt"])]
    df["pfs_capped"] = False
    hist = load_histology("GSE190265")
    df = df.merge(hist, on="sample", how="left")
    return df[["sample", "pfs_time", "pfs_evt", "DCB", "histology", "pfs_capped"]]


def load_series_matrix_rows(path):
    with gzip.open(path, "rt") as fh:
        lines = fh.read().splitlines()
    titles = None
    chars = []
    for ln in lines:
        if ln.startswith("!Sample_title"):
            titles = [x.strip().strip('"') for x in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_characteristics_ch1"):
            chars.append([x.strip().strip('"') for x in ln.split("\t")[1:]])
    return titles, chars


def load_histology(acc):
    p = os.path.join(DATA, acc, f"{acc}_series_matrix.txt.gz")
    titles, chars = load_series_matrix_rows(p)
    hist = []
    for vals in chars:
        if vals and vals[0].lower().startswith("disease state"):
            hist = [re.sub(r"^disease state:\s*", "", v, flags=re.I) for v in vals]
            break
    out = pd.DataFrame({"sample": titles, "histology": hist})
    out["histology"] = out["histology"].replace(
        {"NA": np.nan, "na": np.nan, "5": np.nan, "": np.nan}
    )
    return out


def load_clinical_266():
    p = os.path.join(DATA, "GSE190266", "GSE190266_series_matrix.txt.gz")
    titles, chars = load_series_matrix_rows(p)
    pfs_time = pfs_evt = hist = None
    for vals in chars:
        joined = " ".join(vals)
        if "pfs_time" in joined:
            pfs_time = [parse_decimal(re.sub(r".*:\s*", "", v)) for v in vals]
        elif "pfs_evt" in joined:
            pfs_evt = [parse_decimal(re.sub(r".*:\s*", "", v)) for v in vals]
        elif vals and vals[0].lower().startswith("disease state"):
            hist = [re.sub(r"^disease state:\s*", "", v, flags=re.I) for v in vals]
    df = pd.DataFrame(
        {
            "sample": titles,
            "pfs_time": pfs_time,
            "pfs_evt": pfs_evt,
            "histology": hist,
        }
    )
    df["sample"] = df["sample"].astype(str)
    df["histology"] = df["histology"].replace(
        {"NA": np.nan, "na": np.nan, "5": np.nan, "": np.nan}
    )
    df["DCB"] = [classify_dcb(t, e) for t, e in zip(df["pfs_time"], df["pfs_evt"])]
    df["pfs_capped"] = True
    return df[["sample", "pfs_time", "pfs_evt", "DCB", "histology", "pfs_capped"]]


def load_expression(path):
    """Return samples x genes TPM for TACSTD2/CLDN4 (whichever exist).

    GSE190265 header has no leading blank: data has one extra leading sample field.
    GSE190266 header has a leading blank and uses comma decimals.
    """
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n")
        body = fh.readlines()
    cols = header.split(";")
    rows = [ln.rstrip("\n").split(";") for ln in body]
    n_data = len(rows[0])
    if cols[0] == "":
        genes = cols[1:]
        samples = [r[0] for r in rows]
        mat = [r[1:] for r in rows]
        if len(genes) != len(mat[0]):
            raise ValueError(f"{path}: gene/value length mismatch")
    elif n_data == len(cols) + 1:
        genes = cols
        samples = [r[0] for r in rows]
        mat = [r[1:] for r in rows]
    else:
        raise ValueError(
            f"{path}: cannot align header ({len(cols)}) to data ({n_data})"
        )

    keep_idx = {g: genes.index(g) for g in GENES if g in genes}
    data = {"sample": [str(s) for s in samples]}
    for g, j in keep_idx.items():
        data[g] = [parse_decimal(r[j]) for r in mat]
    df = pd.DataFrame(data).set_index("sample")
    df.attrs["n_gene_columns"] = len(genes)
    df.attrs["last_gene"] = genes[-1]
    df.attrs["missing_genes"] = [g for g in GENES if g not in df.columns]
    return df


def rank_biserial_auc(x_hi, x_lo):
    n1, n2 = len(x_hi), len(x_lo)
    if n1 == 0 or n2 == 0:
        return np.nan, np.nan
    U, _ = stats.mannwhitneyu(x_hi, x_lo, alternative="two-sided")
    auc = U / (n1 * n2)
    return auc, 2 * auc - 1


def logrank_two_group(time, event, group):
    """Mantel-Haenszel log-rank for a binary group (1 vs 0)."""
    t = np.asarray(time, dtype=float)
    e = np.asarray(event, dtype=int)
    g = np.asarray(group, dtype=int)
    ok = np.isfinite(t) & np.isfinite(e) & np.isfinite(g)
    t, e, g = t[ok], e[ok], g[ok]
    if g.min() == g.max() or e.sum() < 1:
        return np.nan, np.nan
    times = np.unique(t[e == 1])
    o1 = e1 = v = 0.0
    for ti in times:
        at_risk = t >= ti
        n = at_risk.sum()
        n1 = ((g == 1) & at_risk).sum()
        n0 = n - n1
        if n1 == 0 or n0 == 0:
            continue
        d = ((t == ti) & (e == 1)).sum()
        d1 = ((t == ti) & (e == 1) & (g == 1)).sum()
        o1 += d1
        e1 += d * n1 / n
        v += (d * (n - d) * n1 * n0) / (n * n * (n - 1)) if n > 1 else 0.0
    if v <= 0:
        return np.nan, np.nan
    z = (o1 - e1) / np.sqrt(v)
    p = 2 * stats.norm.sf(abs(z))
    return float(z), float(p)


def per_cohort_stats(expr, clin, cohort):
    m = expr.join(clin.set_index("sample"), how="inner")
    labeled = m[m["DCB"].isin(["DCB", "NDB"])].copy()
    rows = []
    for g in GENES:
        if g not in labeled.columns:
            rows.append(
                dict(
                    cohort=cohort,
                    gene=g,
                    n_DCB=0,
                    n_NDB=0,
                    wilcoxon_p=np.nan,
                    note="gene absent from deposited TPM matrix",
                )
            )
            continue
        log = np.log2(labeled[g].astype(float) + 1.0)
        dcb = log[labeled["DCB"] == "DCB"].to_numpy()
        ndb = log[labeled["DCB"] == "NDB"].to_numpy()
        p = np.nan
        if len(dcb) >= 3 and len(ndb) >= 3:
            _, p = stats.mannwhitneyu(dcb, ndb, alternative="two-sided")
        auc, rbc = rank_biserial_auc(dcb, ndb)
        # continuous PFS association (flag if time is 6-month-capped)
        rho, prho = stats.spearmanr(log, labeled["pfs_time"])
        med = np.median(log)
        zlr, plr = logrank_two_group(
            labeled["pfs_time"], labeled["pfs_evt"], (log >= med).astype(int)
        )
        rows.append(
            dict(
                cohort=cohort,
                gene=g,
                n_DCB=int(len(dcb)),
                n_NDB=int(len(ndb)),
                median_TPM_DCB=float(np.median(labeled.loc[labeled["DCB"] == "DCB", g])),
                median_TPM_NDB=float(np.median(labeled.loc[labeled["DCB"] == "NDB", g])),
                median_log2_DCB=float(np.median(dcb)),
                median_log2_NDB=float(np.median(ndb)),
                delta_median_log2=float(np.median(dcb) - np.median(ndb)),
                AUC_DCB_gt_NDB=auc,
                rank_biserial=rbc,
                wilcoxon_p=p,
                spearman_rho_vs_PFS=float(rho) if np.isfinite(rho) else np.nan,
                spearman_p=float(prho) if np.isfinite(prho) else np.nan,
                logrank_high_vs_low_Z=zlr,
                logrank_high_vs_low_p=plr,
                pfs_time_capped=bool(labeled["pfs_capped"].iloc[0]),
                note="",
            )
        )
    return pd.DataFrame(rows), labeled


def van_elteren(pooled, gene):
    S = E = V = 0.0
    used = 0
    for _, sub in pooled[pooled["gene"] == gene].groupby("cohort"):
        vals = sub["log2"].to_numpy()
        grp = (sub["DCB"] == "DCB").to_numpy()
        n, n1 = len(vals), int(grp.sum())
        n2 = n - n1
        if n1 < 1 or n2 < 1:
            continue
        ranks = stats.rankdata(vals)
        W = ranks[grp].sum()
        EW = n1 * (n + 1) / 2.0
        _, counts = np.unique(vals, return_counts=True)
        tie = (counts**3 - counts).sum()
        VW = n1 * n2 / 12.0 * ((n + 1) - (tie / (n * (n - 1)) if n > 1 else 0))
        S += W
        E += EW
        V += VW
        used += 1
    if used == 0 or V <= 0:
        return np.nan, np.nan, used
    Z = (S - E) / np.sqrt(V)
    return float(Z), float(2 * stats.norm.sf(abs(Z))), used


def pooled_z_test(pooled, gene):
    sub = pooled[pooled["gene"] == gene].copy()
    if sub.empty:
        return dict(gene=gene, n_DCB=0, n_NDB=0, wilcoxon_p=np.nan, note="no data")
    sub["z"] = sub.groupby("cohort")["log2"].transform(
        lambda v: (v - v.mean()) / v.std(ddof=0) if v.std(ddof=0) > 0 else 0.0
    )
    dcb = sub.loc[sub["DCB"] == "DCB", "z"].to_numpy()
    ndb = sub.loc[sub["DCB"] == "NDB", "z"].to_numpy()
    p = np.nan
    if len(dcb) >= 3 and len(ndb) >= 3:
        _, p = stats.mannwhitneyu(dcb, ndb, alternative="two-sided")
    auc, rbc = rank_biserial_auc(dcb, ndb)
    return dict(
        gene=gene,
        n_DCB=int(len(dcb)),
        n_NDB=int(len(ndb)),
        n_cohorts=int(sub["cohort"].nunique()),
        median_z_DCB=float(np.median(dcb)) if len(dcb) else np.nan,
        median_z_NDB=float(np.median(ndb)) if len(ndb) else np.nan,
        AUC_DCB_gt_NDB=auc,
        rank_biserial=rbc,
        wilcoxon_p=p,
    )


def bh_adjust(pvals):
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    mask = np.isfinite(p)
    pv = p[mask]
    n = len(pv)
    if n == 0:
        return out
    order = np.argsort(pv)
    adj = pv[order] * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    res = np.empty(n)
    res[order] = adj
    out[mask] = res
    return out


def md_table(df, cols=None, nd=4):
    if cols is None:
        cols = list(df.columns)
    def fmt(v):
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            return "NA"
        if isinstance(v, (float, np.floating)):
            return f"{float(v):.{nd}f}"
        return str(v)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(fmt(r.get(c)) for c in cols) + " |")
    return "\n".join(lines)


def plot_box(pooled, missing, outpath):
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.4))
    panels = ["GSE190265", "GSE190266", "POOLED (z-scored)"]
    for r, g in enumerate(GENES):
        for c, panel in enumerate(panels):
            ax = axes[r][c]
            if panel.startswith("POOLED"):
                sub = pooled[pooled["gene"] == g].copy()
                if sub.empty:
                    ax.text(0.5, 0.5, f"{g} unavailable\n{missing.get(g, '')}",
                            ha="center", va="center", transform=ax.transAxes)
                    ax.set_axis_off()
                    continue
                sub["val"] = sub.groupby("cohort")["log2"].transform(
                    lambda v: (v - v.mean()) / v.std(ddof=0) if v.std(ddof=0) > 0 else 0.0
                )
                ylab = "within-cohort z(log2 TPM+1)"
            else:
                sub = pooled[(pooled["gene"] == g) & (pooled["cohort"] == panel)].copy()
                if sub.empty:
                    ax.text(0.5, 0.5, f"{g} not in {panel}\nGEO TPM matrix",
                            ha="center", va="center", transform=ax.transAxes, fontsize=9)
                    ax.set_title(f"{g} — {panel}", fontsize=10)
                    ax.set_axis_off()
                    continue
                sub["val"] = sub["log2"]
                ylab = "log2(TPM + 1)"
            data = [sub.loc[sub["DCB"] == lab, "val"].to_numpy() for lab in ("NDB", "DCB")]
            kw = dict(widths=0.55, showfliers=False, patch_artist=True)
            try:
                bp = ax.boxplot(data, tick_labels=["NDB", "DCB"], **kw)
            except TypeError:
                bp = ax.boxplot(data, labels=["NDB", "DCB"], **kw)
            for patch, col in zip(bp["boxes"], ["#d0d0d0", "#7fb3d5"]):
                patch.set_facecolor(col)
            rng = np.random.default_rng(1)
            for i, arr in enumerate(data, start=1):
                if len(arr):
                    ax.scatter(
                        rng.normal(i, 0.055, len(arr)),
                        arr,
                        s=18,
                        color="#2c3e50",
                        alpha=0.75,
                        zorder=3,
                    )
            ax.set_title(f"{g} — {panel}\nNDB n={len(data[0])}, DCB n={len(data[1])}", fontsize=10)
            ax.set_ylabel(ylab, fontsize=9)
            ax.grid(axis="y", ls=":", alpha=0.45)
    fig.suptitle(
        "TACSTD2 / CLDN4 vs DCB (PFS ≥ 6 mo), anti-PD1 NSCLC (GSE190265/266)",
        fontsize=12,
        y=1.01,
    )
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    clin265 = load_clinical_265()
    clin266 = load_clinical_266()
    expr265 = load_expression(os.path.join(DATA, "GSE190265", "GSE190265_TPM_France3.csv.gz"))
    expr266 = load_expression(os.path.join(DATA, "GSE190266", "GSE190266_TPM_France4.csv.gz"))

    stats265, m265 = per_cohort_stats(expr265, clin265, "GSE190265")
    stats266, m266 = per_cohort_stats(expr266, clin266, "GSE190266")
    per_cohort = pd.concat([stats265, stats266], ignore_index=True)
    per_cohort["wilcoxon_p_BH"] = bh_adjust(per_cohort["wilcoxon_p"].to_numpy())

    def to_long(m, cohort):
        recs = []
        for s, r in m.iterrows():
            for g in GENES:
                if g in m.columns:
                    recs.append(
                        dict(
                            cohort=cohort,
                            sample=s,
                            gene=g,
                            DCB=r["DCB"],
                            histology=r.get("histology"),
                            pfs_time=r["pfs_time"],
                            pfs_evt=r["pfs_evt"],
                            tpm=float(r[g]),
                            log2=float(np.log2(float(r[g]) + 1.0)),
                        )
                    )
        return pd.DataFrame(recs)

    pooled = pd.concat([to_long(m265, "GSE190265"), to_long(m266, "GSE190266")], ignore_index=True)

    strat_rows = []
    pz_rows = []
    for g in GENES:
        z, p, ncoh = van_elteren(pooled, g)
        strat_rows.append(dict(gene=g, n_cohorts=ncoh, van_elteren_Z=z, van_elteren_p=p))
        pz_rows.append(pooled_z_test(pooled, g))
    strat = pd.DataFrame(strat_rows)
    pooledz = pd.DataFrame(pz_rows)

    # histology x DCB (confounder check)
    hist_rows = []
    for cohort, clin in (("GSE190265", clin265), ("GSE190266", clin266)):
        sub = clin[clin["DCB"].isin(["DCB", "NDB"])]
        ct = pd.crosstab(sub["histology"].fillna("unknown"), sub["DCB"])
        for h, r in ct.iterrows():
            hist_rows.append(
                dict(
                    cohort=cohort,
                    histology=h,
                    n_DCB=int(r.get("DCB", 0)),
                    n_NDB=int(r.get("NDB", 0)),
                )
            )
    hist_tab = pd.DataFrame(hist_rows)

    # gene vs histology (GSE190265 both genes; 266 CLDN4 only)
    hist_expr_rows = []
    for cohort, m in (("GSE190265", m265), ("GSE190266", m266)):
        mm = m[m["histology"].isin(["squamous", "non-squamous"])].copy()
        for g in GENES:
            if g not in mm.columns:
                continue
            a = np.log2(mm.loc[mm["histology"] == "non-squamous", g].astype(float) + 1)
            b = np.log2(mm.loc[mm["histology"] == "squamous", g].astype(float) + 1)
            p = np.nan
            if len(a) >= 3 and len(b) >= 3:
                _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            hist_expr_rows.append(
                dict(
                    cohort=cohort,
                    gene=g,
                    n_nonsquamous=int(len(a)),
                    n_squamous=int(len(b)),
                    median_log2_nonsquamous=float(np.median(a)) if len(a) else np.nan,
                    median_log2_squamous=float(np.median(b)) if len(b) else np.nan,
                    wilcoxon_p=p,
                )
            )
    hist_expr = pd.DataFrame(hist_expr_rows)

    clin_all = pd.concat(
        [clin265.assign(cohort="GSE190265"), clin266.assign(cohort="GSE190266")],
        ignore_index=True,
    )
    clin_all.to_csv(os.path.join(OUTDIR, "clinical_dcb_labels.csv"), index=False)
    per_cohort.to_csv(os.path.join(OUTDIR, "per_cohort_stats.csv"), index=False)
    strat.to_csv(os.path.join(OUTDIR, "stratified_van_elteren.csv"), index=False)
    pooledz.to_csv(os.path.join(OUTDIR, "pooled_zscored_stats.csv"), index=False)
    pooled.to_csv(os.path.join(OUTDIR, "expression_long.csv"), index=False)
    hist_tab.to_csv(os.path.join(OUTDIR, "histology_by_dcb.csv"), index=False)
    hist_expr.to_csv(os.path.join(OUTDIR, "expression_by_histology.csv"), index=False)

    # exploratory: DCB vs NDB inside non-squamous only (histology confounder check)
    ns_rows = []
    for cohort, m in (("GSE190265", m265), ("GSE190266", m266)):
        mm = m[(m["histology"] == "non-squamous") & (m["DCB"].isin(["DCB", "NDB"]))]
        for g in GENES:
            if g not in mm.columns:
                continue
            dcb = np.log2(mm.loc[mm["DCB"] == "DCB", g].astype(float) + 1).to_numpy()
            ndb = np.log2(mm.loc[mm["DCB"] == "NDB", g].astype(float) + 1).to_numpy()
            p = np.nan
            if len(dcb) >= 3 and len(ndb) >= 3:
                _, p = stats.mannwhitneyu(dcb, ndb, alternative="two-sided")
            auc, rbc = rank_biserial_auc(dcb, ndb)
            ns_rows.append(
                dict(
                    cohort=cohort,
                    gene=g,
                    n_DCB=int(len(dcb)),
                    n_NDB=int(len(ndb)),
                    median_log2_DCB=float(np.median(dcb)) if len(dcb) else np.nan,
                    median_log2_NDB=float(np.median(ndb)) if len(ndb) else np.nan,
                    AUC_DCB_gt_NDB=auc,
                    wilcoxon_p=p,
                )
            )
    ns_tab = pd.DataFrame(ns_rows)
    ns_tab.to_csv(os.path.join(OUTDIR, "nonsquamous_only_stats.csv"), index=False)

    missing = {
        "TACSTD2": (
            f"GSE190266 deposited matrix has {expr266.attrs['n_gene_columns']} gene columns "
            f"and ends at {expr266.attrs['last_gene']}"
        )
    }
    plot_box(pooled, missing, os.path.join(OUTDIR, "TACSTD2_CLDN4_DCB_boxplots.png"))

    def counts(clin, expr):
        overlap = set(clin["sample"]) & set(expr.index)
        dcb = clin.loc[clin["sample"].isin(overlap), "DCB"]
        return dict(
            n_clin=len(clin),
            n_expr=len(expr),
            n_overlap=len(overlap),
            n_DCB=int((dcb == "DCB").sum()),
            n_NDB=int((dcb == "NDB").sum()),
            n_excluded=int(dcb.isna().sum()),
        )

    c265, c266 = counts(clin265, expr265), counts(clin266, expr266)

    # write a numbers-only summary; REPORT.md is authored after inspection
    lines = [
        "# Machine summary (numbers only)",
        "",
        f"GSE190265 TPM genes in file: {expr265.attrs['n_gene_columns']} (last={expr265.attrs['last_gene']}); missing={expr265.attrs['missing_genes']}",
        f"GSE190266 TPM genes in file: {expr266.attrs['n_gene_columns']} (last={expr266.attrs['last_gene']}); missing={expr266.attrs['missing_genes']}",
        "",
        f"GSE190265 match: {c265}",
        f"GSE190266 match: {c266}",
        "",
        "## per-cohort",
        md_table(per_cohort),
        "",
        "## van Elteren",
        md_table(strat),
        "",
        "## pooled z",
        md_table(pooledz),
        "",
        "## histology x DCB",
        md_table(hist_tab, nd=0),
        "",
        "## expression vs histology",
        md_table(hist_expr),
    ]
    with open(os.path.join(OUTDIR, "machine_summary.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
