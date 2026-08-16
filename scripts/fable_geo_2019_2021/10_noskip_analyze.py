#!/usr/bin/env python3
"""Compute TACSTD2/CLDN4 on leftover 2019-2021 GEO series (no skip for size/tissue).

Only statistics that can be computed from deposited files are written.
No IDs or p-values are invented. Series without the genes, without a processed
matrix, or without a deposited outcome are recorded as such.
"""
import csv
import gzip
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "noskip" / "GEO_2019_2021"
DL = OUT / "downloads"
TAB = OUT / "tables"
FIG = OUT / "figures"
TAB.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

rows = []


def add(**kw):
    rows.append(kw)
    print("{dataset} | {gene} | {outcome} | n={n} | {effect} | p={p_value}".format(
        dataset=kw.get("dataset"), gene=kw.get("gene"), outcome=kw.get("outcome"),
        n=kw.get("n"), effect=kw.get("effect"), p_value=kw.get("p_value")))


def auc_mw(pos, neg):
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan"), float("nan")
    u, p = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    return u / (len(pos) * len(neg)), p


def parse_series_chars(path):
    table = None
    with gzip.open(path, "rt", errors="replace") as f:
        titles = gsms = None
        chars = []
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                gsms = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
            elif line.startswith("!Sample_characteristics"):
                vals = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
                field = None
                for v in vals:
                    if ":" in v:
                        field = v.split(":", 1)[0].strip()
                        break
                if field:
                    parsed = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
                    chars.append((field, parsed))
    n = len(gsms or titles or [])
    recs = []
    for i in range(n):
        r = {"gsm": gsms[i] if gsms else "", "title": titles[i] if titles else ""}
        for field, parsed in chars:
            r[field] = parsed[i]
        recs.append(r)
    return pd.DataFrame(recs)


def boxplot(df, group, genes, gse, ylab, order=None):
    order = order or list(df[group].dropna().unique())
    fig, axes = plt.subplots(1, len(genes), figsize=(4 * len(genes), 4))
    if len(genes) == 1:
        axes = [axes]
    for ax, g in zip(axes, genes):
        data = [df.loc[df[group] == lvl, g].dropna().values for lvl in order]
        ax.boxplot(data, tick_labels=order, showfliers=False)
        for i, lvl in enumerate(order, start=1):
            y = df.loc[df[group] == lvl, g].dropna().values
            x = np.random.normal(i, 0.05, size=len(y))
            ax.scatter(x, y, s=16, alpha=0.7, color="#377eb8")
        ax.set_title(f"{gse}: {g}")
        ax.set_ylabel(ylab)
        ax.tick_params(axis="x", rotation=15)
    fig.tight_layout()
    fig.savefig(FIG / f"{gse}_{group}_boxplot.png", dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------
# GSE190266 — leftover ICI outcome: 6-month PFS. CLDN4 present; TACSTD2 absent.
# ---------------------------------------------------------------------------
def gse190266():
    gse = "GSE190266"
    clin = parse_series_chars(DL / gse / f"{gse}_series_matrix.txt.gz")
    tpm = pd.read_csv(DL / gse / f"{gse}_TPM_France4.csv.gz", sep=";", decimal=",")
    tpm = tpm.rename(columns={tpm.columns[0]: "sample"})
    # European CSV: first col is sample id
    if "CLDN4" not in tpm.columns:
        add(dataset=gse, compartment="tumor biopsy", gene="CLDN4",
            outcome="PFS 6-month (deposited)", n=len(clin),
            n_group1="", n_group2="", effect="CLDN4 column missing after parse",
            p_value="", direction="n/a", note="parse error")
        return
    df = clin.copy()
    df["sample"] = df["title"]
    df = df.merge(tpm[["sample", "CLDN4"]], on="sample", how="inner")
    df["CLDN4_log"] = np.log2(df["CLDN4"].astype(float) + 1)
    df["pfs_time"] = pd.to_numeric(df["pfs_time (6 months)"], errors="coerce")
    df["pfs_evt"] = pd.to_numeric(df["pfs_evt (6 months)"], errors="coerce")
    df = df.dropna(subset=["CLDN4_log", "pfs_time", "pfs_evt"])
    df.to_csv(TAB / f"{gse}_CLDN4_PFS.csv", index=False)
    cph = CoxPHFitter()
    cph.fit(df[["CLDN4_log", "pfs_time", "pfs_evt"]],
            duration_col="pfs_time", event_col="pfs_evt")
    hr = math.exp(cph.params_["CLDN4_log"])
    p_cox = float(cph.summary.loc["CLDN4_log", "p"])
    med = df["CLDN4_log"].median()
    high = df[df["CLDN4_log"] > med]
    low = df[df["CLDN4_log"] <= med]
    lr = logrank_test(high["pfs_time"], low["pfs_time"], high["pfs_evt"], low["pfs_evt"])
    add(dataset=gse, compartment="tumor biopsy", gene="CLDN4",
        outcome="6-month PFS (Cox per log2(TPM+1) / KM median-split)",
        n=len(df), n_group1=len(high), n_group2=len(low),
        effect=f"HR={hr:.3f}; logrank p={lr.p_value:.4f}",
        p_value=round(p_cox, 4),
        direction="higher expr -> shorter PFS" if hr > 1 else "higher expr -> longer PFS",
        note="TACSTD2 not in deposited TPM matrix (16383 genes; CLDN4 present). PFS is 6-month truncated as deposited.")
    add(dataset=gse, compartment="tumor biopsy", gene="TACSTD2",
        outcome="6-month PFS", n=len(clin), n_group1="", n_group2="",
        effect="gene absent from deposited TPM matrix", p_value="",
        direction="n/a", note="not computed")
    kmf = KaplanMeierFitter()
    fig, ax = plt.subplots(figsize=(5, 4))
    kmf.fit(high["pfs_time"], high["pfs_evt"], label=f"CLDN4 high n={len(high)}")
    kmf.plot_survival_function(ax=ax, ci_show=False)
    kmf.fit(low["pfs_time"], low["pfs_evt"], label=f"CLDN4 low n={len(low)}")
    kmf.plot_survival_function(ax=ax, ci_show=False)
    ax.set_title("GSE190266: 6-month PFS by CLDN4")
    ax.set_xlabel("PFS time (months, truncated at 6)")
    fig.tight_layout()
    fig.savefig(FIG / f"{gse}_CLDN4_KM.png", dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------
# GSE190265 — same study, no PFS deposited. Both genes present. Histology only.
# ---------------------------------------------------------------------------
def gse190265():
    gse = "GSE190265"
    clin = parse_series_chars(DL / gse / f"{gse}_series_matrix.txt.gz")
    # Header is genes-only; each data row is sample_id + values (one extra field).
    with gzip.open(DL / gse / f"{gse}_TPM_France3.csv.gz", "rt") as f:
        genes = f.readline().rstrip("\n").split(";")
    tpm = pd.read_csv(DL / gse / f"{gse}_TPM_France3.csv.gz", sep=";", header=None,
                      skiprows=1, names=["sample"] + genes)
    tpm["sample"] = tpm["sample"].astype(str)
    tpm = tpm.set_index("sample")
    clin["sample"] = clin["title"].astype(str)
    keep = [s for s in clin["sample"] if s in tpm.index]
    df = clin[clin["sample"].isin(keep)].copy()
    for g in ["TACSTD2", "CLDN4"]:
        df[g] = np.log2(tpm.loc[df["sample"], g].astype(float).values + 1)
    df.to_csv(TAB / f"{gse}_expr_histology.csv", index=False)
    hist = df["disease state"].replace({"NA": np.nan})
    df = df.assign(histology=hist)
    sub = df[df["histology"].isin(["squamous", "non-squamous"])]
    for g in ["TACSTD2", "CLDN4"]:
        a = sub.loc[sub.histology == "squamous", g]
        b = sub.loc[sub.histology == "non-squamous", g]
        auc, p = auc_mw(a, b)
        add(dataset=gse, compartment="tumor biopsy", gene=g,
            outcome="histology (squamous vs non-squamous) — NOT an ICI endpoint",
            n=len(sub), n_group1=len(a), n_group2=len(b),
            effect=f"AUC_sq_vs_nsq={auc:.3f}; med_sq={a.median():.3f}; med_nsq={b.median():.3f}",
            p_value=round(p, 4),
            direction="higher in squamous" if a.median() > b.median() else "higher in non-squamous",
            note="No response/PFS deposited in GEO for this series.")
    boxplot(sub, "histology", ["TACSTD2", "CLDN4"], gse, "log2(TPM+1)",
            order=["non-squamous", "squamous"])


# ---------------------------------------------------------------------------
# GSE181820 — NSCLC RNA-seq, GEO groups A/B/C (meaning not deposited).
# ---------------------------------------------------------------------------
def gse181820():
    gse = "GSE181820"
    clin = parse_series_chars(DL / gse / f"{gse}_series_matrix.txt.gz")
    expr = pd.read_csv(DL / gse / f"{gse}_RNAseq_NSCLC.txt.gz", sep="\t", index_col=0)
    # columns are patient numbers matching "NSCLC {n}"
    clin["pt"] = clin["title"].str.replace("NSCLC ", "", regex=False)
    # expression columns may be int-like
    expr.columns = [str(c) for c in expr.columns]
    recs = []
    for _, r in clin.iterrows():
        pt = r["pt"]
        if pt not in expr.columns:
            continue
        recs.append({
            "sample": r["title"], "group": r["group"],
            "TACSTD2": float(expr.loc["TACSTD2", pt]),
            "CLDN4": float(expr.loc["CLDN4", pt]),
        })
    df = pd.DataFrame(recs)
    df["TACSTD2_log"] = np.log2(df["TACSTD2"] + 1)
    df["CLDN4_log"] = np.log2(df["CLDN4"] + 1)
    df.to_csv(TAB / f"{gse}_expr_group.csv", index=False)
    for g, col in [("TACSTD2", "TACSTD2_log"), ("CLDN4", "CLDN4_log")]:
        groups = [df.loc[df.group == lvl, col].values for lvl in ["A", "B", "C"]]
        try:
            h, p = stats.kruskal(*groups)
        except Exception:  # noqa: BLE001
            h, p = float("nan"), float("nan")
        meds = {lvl: float(df.loc[df.group == lvl, col].median()) for lvl in ["A", "B", "C"]}
        add(dataset=gse, compartment="tumor", gene=g,
            outcome="GEO group A/B/C (Kruskal–Wallis) — NOT a deposited ICI response",
            n=len(df), n_group1=int((df.group == "A").sum()),
            n_group2=int((df.group == "B").sum()) + int((df.group == "C").sum()),
            effect=f"H={h:.3f}; medA={meds['A']:.3f} medB={meds['B']:.3f} medC={meds['C']:.3f}",
            p_value=round(p, 4) if p == p else "",
            direction="see medians",
            note="Group labels A/B/C are deposited; their biological meaning is not in series_matrix. Do not treat as ICI outcome.")
    boxplot(df, "group", ["TACSTD2_log", "CLDN4_log"], gse, "log2(expr+1)",
            order=["A", "B", "C"])


# ---------------------------------------------------------------------------
# GSE146100 — 3 nodules / 1 patient, pembrolizumab. scRNA means only.
# ---------------------------------------------------------------------------
def gse146100():
    gse = "GSE146100"
    path = DL / gse / f"{gse}_NormData.txt.gz"
    found = {}
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip().split("\t")
        cols = header[1:]
        prefixes = [c.split("_", 1)[0] for c in cols]
        for line in f:
            gene = line.split("\t", 1)[0]
            if gene in ("TACSTD2", "CLDN4"):
                vals = np.array(line.rstrip().split("\t")[1:], dtype=float)
                found[gene] = vals
            if len(found) == 2:
                break
    recs = []
    for nod in ["W1", "W2", "W3"]:
        idx = [i for i, p in enumerate(prefixes) if p == nod]
        rec = {"nodule": nod, "n_cells": len(idx),
               "response_label": "responder" if nod == "W2" else "non-responder"}
        for g, vals in found.items():
            rec[f"{g}_mean"] = float(vals[idx].mean()) if idx else float("nan")
            rec[f"{g}_pct_pos"] = float((vals[idx] > 0).mean()) if idx else float("nan")
        recs.append(rec)
    df = pd.DataFrame(recs)
    df.to_csv(TAB / f"{gse}_per_nodule_means.csv", index=False)
    for g in found:
        add(dataset=gse, compartment="scRNA nodules (1 patient)", gene=g,
            outcome="per-nodule mean (W2 responded vs W1/W3 non-responded)",
            n=3, n_group1=1, n_group2=2,
            effect="; ".join(f"{r['nodule']}({r['response_label']}) mean={r[f'{g}_mean']:.4f} pct>0={r[f'{g}_pct_pos']:.3f}"
                             for _, r in df.iterrows()),
            p_value="not tested (n=3 nodules, 1 patient)",
            direction="descriptive only",
            note="No inferential test. Single-patient scRNA.")


# ---------------------------------------------------------------------------
# GSE145896 — PD-1+ TIL CD28+ vs CD28- (not a patient ICI endpoint)
# ---------------------------------------------------------------------------
def gse145896():
    gse = "GSE145896"
    expr = pd.read_csv(DL / gse / f"{gse}_countData.csv.gz", index_col=0)
    plus = [c for c in expr.columns if "CD28_plus" in c]
    minus = [c for c in expr.columns if "CD28_minus" in c]
    recs = []
    for c in expr.columns:
        recs.append({
            "sample": c,
            "subset": "CD28+" if "CD28_plus" in c else "CD28-",
            "TACSTD2": float(expr.loc["TACSTD2", c]),
            "CLDN4": float(expr.loc["CLDN4", c]),
        })
    df = pd.DataFrame(recs)
    df.to_csv(TAB / f"{gse}_TIL_subsets.csv", index=False)
    for g in ["TACSTD2", "CLDN4"]:
        a = df.loc[df.subset == "CD28+", g]
        b = df.loc[df.subset == "CD28-", g]
        auc, p = auc_mw(a, b)
        add(dataset=gse, compartment="PD-1+ TILs", gene=g,
            outcome="CD28+ vs CD28- subset (NOT patient ICI response)",
            n=len(df), n_group1=len(a), n_group2=len(b),
            effect=f"AUC(+ vs -)={auc:.3f}; med+={a.median():.1f}; med-={b.median():.1f}",
            p_value=round(p, 4),
            direction="higher in CD28+" if a.median() > b.median() else "higher in CD28-",
            note="Paired subset RNA from same tumors; not a patient-level ICI endpoint.")
    boxplot(df, "subset", ["TACSTD2", "CLDN4"], gse, "raw counts",
            order=["CD28+", "CD28-"])


# ---------------------------------------------------------------------------
# GSE131933 — explant scRNA, durvalumab. Tumor-cell gene means by treatment.
# ---------------------------------------------------------------------------
def gse131933():
    gse = "GSE131933"
    files = [
        "GSE131933_T1T3D0_gene_count.txt.gz",
        "GSE131933_T3D11_gene_count.txt.gz",
    ]
    # T1D7 (958 MB) has no TACSTD2/CLDN4 gene symbols (verified by scan).
    add(dataset=gse, compartment="ex-vivo NSCLC explant scRNA", gene="TACSTD2/CLDN4",
        outcome="T1D7 file (958 MB)", n="", n_group1="", n_group2="",
        effect="gene symbols TACSTD2 and CLDN4 absent from this file",
        p_value="", direction="n/a",
        note="Downloaded (no size skip). File uses gene symbols; targets not present.")
    collected = []
    for fname in files:
        path = DL / gse / fname
        with gzip.open(path, "rt") as f:
            header = f.readline().rstrip()
            # quoted tokens
            cols = re.findall(r'"([^"]+)"', header)
            if not cols:
                cols = header.split()
            gene_rows = {}
            for line in f:
                tok = line.split(None, 1)[0].strip('"')
                if tok in ("TACSTD2", "CLDN4"):
                    rest = line.split(None, 1)[1] if " " in line or "\t" in line else ""
                    vals = np.array(rest.split(), dtype=float)
                    gene_rows[tok] = vals
                if len(gene_rows) == 2:
                    break
        for i, col in enumerate(cols):
            # names like T1D0_BARCODE_Tumor / T1D7_..._1_NK
            parts = col.split("_")
            treat = parts[0] if parts else col
            ctype = parts[-1]
            rec = {"file": fname, "cell": col, "sample_treat": treat, "cell_type": ctype}
            for g, vals in gene_rows.items():
                rec[g] = float(vals[i]) if i < len(vals) else float("nan")
            collected.append(rec)
    if not collected:
        return
    df = pd.DataFrame(collected)
    tumor = df[df["cell_type"].str.contains("Tumor", case=False, na=False)]
    summary = (tumor.groupby("sample_treat")[["TACSTD2", "CLDN4"]]
               .agg(["mean", "median", "count"]))
    # flatten
    flat = tumor.groupby("sample_treat").agg(
        n=("cell", "size"),
        TACSTD2_mean=("TACSTD2", "mean"),
        TACSTD2_pct=("TACSTD2", lambda s: float((s > 0).mean())),
        CLDN4_mean=("CLDN4", "mean"),
        CLDN4_pct=("CLDN4", lambda s: float((s > 0).mean())),
    ).reset_index()
    flat.to_csv(TAB / f"{gse}_tumor_cells_by_treatment.csv", index=False)
    add(dataset=gse, compartment="ex-vivo explant Tumor cells", gene="TACSTD2",
        outcome="descriptive mean by sample_treat (D0 / Durva / combo) — not patient ICI survival",
        n=int(flat["n"].sum()), n_group1="", n_group2="",
        effect="; ".join(f"{r.sample_treat} n={int(r.n)} mean={r.TACSTD2_mean:.4f} pct>0={r.TACSTD2_pct:.3f}"
                         for r in flat.itertuples()),
        p_value="not tested (explant timepoints, not independent patients)",
        direction="descriptive only",
        note="Durvalumab/tremelimumab applied ex vivo. Not a clinical ICI endpoint.")
    add(dataset=gse, compartment="ex-vivo explant Tumor cells", gene="CLDN4",
        outcome="descriptive mean by sample_treat — not patient ICI survival",
        n=int(flat["n"].sum()), n_group1="", n_group2="",
        effect="; ".join(f"{r.sample_treat} n={int(r.n)} mean={r.CLDN4_mean:.4f} pct>0={r.CLDN4_pct:.3f}"
                         for r in flat.itertuples()),
        p_value="not tested (explant timepoints, not independent patients)",
        direction="descriptive only",
        note="Durvalumab/tremelimumab applied ex vivo. Not a clinical ICI endpoint.")


# ---------------------------------------------------------------------------
# GSE184053 — Temra tumor vs blood. Epithelial genes ~0.
# ---------------------------------------------------------------------------
def gse184053():
    gse = "GSE184053"
    frames = []
    for fname in ["GSE184053_raw_gene_counts_1.txt.gz",
                  "GSE184053_raw_gene_counts_2.txt.gz",
                  "GSE184053_raw_gene_counts_3.txt.gz"]:
        df = pd.read_csv(DL / gse / fname, sep="\t", index_col=0)
        frames.append(df)
    expr = pd.concat(frames, axis=1)
    # drop duplicate columns if any
    expr = expr.loc[:, ~expr.columns.duplicated()]
    genes = {}
    for ens, name in [("ENSG00000184292", "TACSTD2"), ("ENSG00000189143", "CLDN4")]:
        if ens in expr.index:
            genes[name] = expr.loc[ens]
    recs = []
    for col in expr.columns:
        rec = {"sample": col,
               "compartment": "tumor_TIL" if col.startswith("til") else "blood"}
        for name, s in genes.items():
            rec[name] = float(s[col]) if col in s.index else float("nan")
        recs.append(rec)
    df = pd.DataFrame(recs)
    df.to_csv(TAB / f"{gse}_Temra.csv", index=False)
    for g in genes:
        add(dataset=gse, compartment="Temra (tumor vs blood)", gene=g,
            outcome="expression magnitude (epithelial-null check)",
            n=len(df), n_group1=int((df.compartment == "tumor_TIL").sum()),
            n_group2=int((df.compartment == "blood").sum()),
            effect=f"mean={df[g].mean():.4f}; max={df[g].max():.4f}",
            p_value="not an ICI endpoint",
            direction="near-zero as expected in T cells",
            note="No patient ICI response deposited.")


# ---------------------------------------------------------------------------
# Honest non-computes
# ---------------------------------------------------------------------------
def record_absences():
    absences = [
        dict(dataset="GSE152590", compartment="PBMC CD8+ T cells (blood)",
             gene="TACSTD2/CLDN4", outcome="n/a", n=8, n_group1="", n_group2="",
             effect="both genes absent from deposited TPM matrix (24530 genes; CLDN4 not among CLDN* listed)",
             p_value="", direction="n/a",
             note="Blood leftover revisited. Processed xlsx downloaded. No fabricate."),
        dict(dataset="GSE141479", compartment="PBMC CD8+ (blood, PD-1 treated)",
             gene="TACSTD2/CLDN4", outcome="treatment timepoints deposited; no response/PFS",
             n=74, n_group1="", n_group2="",
             effect="series_matrix uses probe IDs (TC*.hg.1); TACSTD2/CLDN4 symbols not present",
             p_value="", direction="n/a",
             note="Blood leftover revisited (57 MB series_matrix downloaded). No symbol mapping invented."),
        dict(dataset="GSE179994", compartment="T-cell scRNA (anti-PD-1 lung)",
             gene="TACSTD2", outcome="pre vs on-treatment only (no response/PFS in GEO)",
             n=47, n_group1="", n_group2="",
             effect="TACSTD2 string present in 441 MB RDS; patient-level matrix not extracted (R object). No stats invented.",
             p_value="", direction="n/a",
             note="Large leftover downloaded (no size skip). T-cell compartment; epithelial genes not a valid ICI-outcome test without cell-type aggregation we did not fabricate."),
        dict(dataset="GSE180347", compartment="LUAD tumor",
             gene="TACSTD2/CLDN4", outcome="PD-L1 group deposited; no processed expression matrix",
             n=144, n_group1="", n_group2="",
             effect="suppl is RAW.tar only; series_matrix has metadata, no gene rows",
             p_value="", direction="n/a",
             note="No open processed matrix. Not computed."),
        dict(dataset="GSE173351", compartment="T cells / multi-site",
             gene="TACSTD2/CLDN4", outcome="no processed expression matrix",
             n=212, n_group1="", n_group2="",
             effect="RAW.tar only", p_value="", direction="n/a",
             note="No open processed matrix. Not computed."),
        dict(dataset="GSE176021", compartment="scRNA annotations",
             gene="TACSTD2/CLDN4", outcome="response status IS deposited; no processed gene matrix",
             n=110, n_group1="", n_group2="",
             effect="suppl has annotation RDS + RAW.tar, no gene-expression matrix",
             p_value="", direction="n/a",
             note="Would have been an ICI-response leftover if a processed matrix existed. Honest skip: no processed expression."),
        dict(dataset="GSE136961", compartment="tumor (Oncomine 395-gene panel)",
             gene="TACSTD2/CLDN4", outcome="DCB/NDB + survival deposited",
             n=21, n_group1="", n_group2="",
             effect="target genes absent from panel (first-pass verified)",
             p_value="", direction="n/a",
             note="Reconfirmed: cannot compute."),
    ]
    for a in absences:
        add(**a)


def main():
    gse190266()
    gse190265()
    gse181820()
    gse146100()
    gse145896()
    gse131933()
    gse184053()
    record_absences()
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "leftover_TACSTD2_CLDN4_computes.csv", index=False)
    print("\nWrote", TAB / "leftover_TACSTD2_CLDN4_computes.csv")


if __name__ == "__main__":
    main()
