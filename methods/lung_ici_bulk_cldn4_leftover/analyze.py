#!/usr/bin/env python3
"""CLDN4 (TACSTD2 companion) vs ICI endpoints and immune axes in leftover OPEN lung ICI bulk.

Public GEO processed matrices only. Patient is the unit. Honest n.
Does not re-audit GSE135222 / GSE190266 / GSE285029 / GSE253564 / GSE248378.
"""
from __future__ import annotations

import gzip
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = Path("/tmp/lung_ici_leftover")
FIG = HERE / "figures"
TAB = HERE / "tables"
PROC = HERE / "processed"
for d in (FIG, TAB, PROC):
    d.mkdir(parents=True, exist_ok=True)

EST = pd.read_csv(HERE / "estimate_gene_sets.csv")
STROMAL = [g for g in EST["stromal_signature"].dropna().astype(str) if g]
IMMUNE = [g for g in EST["immune_signature"].dropna().astype(str) if g]
SYM2ENS = (
    pd.read_csv(HERE / "hgnc_symbol_ensembl.tsv", sep="\t")
    .dropna()
    .drop_duplicates("symbol")
    .set_index("symbol")["ensembl_gene_id"]
    .to_dict()
)
ENS2SYM = {v: k for k, v in SYM2ENS.items()}

IFN_GENES = ["IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA"]
MHC_GENES = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2"]
TARGETS = ["CLDN4", "TACSTD2"]


def _split_fields(line: str) -> list[str]:
    return [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]


def parse_series(path: Path):
    stored = {}
    char_rows = []
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_geo_accession"):
                stored["geo"] = _split_fields(line)
            elif line.startswith("!Sample_title"):
                stored["title"] = _split_fields(line)
            elif line.startswith("!Sample_source_name"):
                stored["source"] = _split_fields(line)
            elif line.startswith("!Sample_description"):
                stored["description"] = _split_fields(line)
            elif line.startswith("!Sample_characteristics"):
                char_rows.append(_split_fields(line))
    geo = stored.get("geo") or []
    meta = pd.DataFrame({"gsm": geo})
    n = len(geo)
    for key in ("title", "source", "description"):
        vals = stored.get(key) or [None] * n
        meta[key] = (vals + [None] * n)[:n]
    chars = defaultdict(dict)
    for vals in char_rows:
        for g, v in zip(geo, vals):
            if not v:
                continue
            if ":" in v:
                k, val = v.split(":", 1)
                chars[g][k.strip().lower()] = val.strip()
            else:
                chars[g]["characteristic"] = v
    keys = sorted({k for d in chars.values() for k in d})
    for k in keys:
        meta[k] = [chars[g].get(k) for g in geo]
    return meta


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    arr = df.to_numpy(dtype=float)
    mu = np.nanmean(arr, axis=1, keepdims=True)
    sd = np.nanstd(arr, axis=1, keepdims=True)
    sd = np.where(sd == 0, np.nan, sd)
    return pd.DataFrame((arr - mu) / sd, index=df.index, columns=df.columns)


def mean_z(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series | None, list[str]]:
    found = [g for g in genes if g in expr.index]
    if len(found) < 2:
        return None, found
    return zscore_rows(expr.loc[found]).mean(axis=0), found


def log2p1(s: pd.Series) -> pd.Series:
    v = pd.to_numeric(s, errors="coerce")
    if v.dropna().empty:
        return v
    mx = float(v.max())
    # already log-like
    if mx <= 20:
        return v
    return np.log2(v.clip(lower=0) + 1)


def counts_to_log2cpm(expr: pd.DataFrame) -> pd.DataFrame:
    lib = expr.sum(axis=0).replace(0, np.nan)
    cpm = expr.div(lib, axis=1) * 1e6
    return np.log2(cpm + 1)


def spearman(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    if m.sum() < 4:
        return {"n": int(m.sum()), "rho": np.nan, "p": np.nan}
    r, p = stats.spearmanr(x[m], y[m])
    return {"n": int(m.sum()), "rho": float(r), "p": float(p)}


def partial_spearman(x, y, z):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    z = np.asarray(z, float)
    m = ~(np.isnan(x) | np.isnan(y) | np.isnan(z))
    n = int(m.sum())
    if n < 5:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rx, ry, rz = stats.rankdata(x[m]), stats.rankdata(y[m]), stats.rankdata(z[m])
    Z = np.column_stack([np.ones(n), rz])
    bx, *_ = np.linalg.lstsq(Z, rx, rcond=None)
    by, *_ = np.linalg.lstsq(Z, ry, rcond=None)
    ex, ey = rx - Z @ bx, ry - Z @ by
    if np.std(ex) == 0 or np.std(ey) == 0:
        return {"n": n, "rho": np.nan, "p": np.nan}
    r = float(np.corrcoef(ex, ey)[0, 1])
    df = n - 3
    t = r * math.sqrt(df / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), df))
    return {"n": n, "rho": r, "p": p}


def residualize(y, covar):
    y = np.asarray(y, float)
    c = np.asarray(covar, float)
    m = ~(np.isnan(y) | np.isnan(c))
    out = np.full_like(y, np.nan, dtype=float)
    if m.sum() < 4:
        return out
    Z = np.column_stack([np.ones(m.sum()), c[m]])
    b, *_ = np.linalg.lstsq(Z, y[m], rcond=None)
    out[m] = y[m] - Z @ b
    return out


def mwu_block(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return {
            "n_R": n1,
            "n_NR": n2,
            "median_R": float(np.median(a)) if n1 else np.nan,
            "median_NR": float(np.median(b)) if n2 else np.nan,
            "U": np.nan,
            "auc_R_gt_NR": np.nan,
            "cliffs_delta": np.nan,
            "p": np.nan,
        }
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    auc = float(U / (n1 * n2))
    return {
        "n_R": n1,
        "n_NR": n2,
        "median_R": float(np.median(a)),
        "median_NR": float(np.median(b)),
        "U": float(U),
        "auc_R_gt_NR": auc,
        "cliffs_delta": float(2 * auc - 1),
        "p": float(p),
    }


def map_index_to_symbol(idx: pd.Index) -> pd.Index:
    out = []
    for x in idx.astype(str):
        x0 = x.split(".")[0]
        if x0.startswith("ENSG") and x0 in ENS2SYM:
            out.append(ENS2SYM[x0])
        elif "_" in x and not x.startswith("ENSG"):
            out.append(x.split("_")[0])
        else:
            out.append(x)
    return pd.Index(out)


def collapse_symbols(expr: pd.DataFrame) -> pd.DataFrame:
    expr = expr.copy()
    expr.index = map_index_to_symbol(expr.index)
    expr = expr.groupby(expr.index).mean()
    return expr


def score_matrix(expr: pd.DataFrame, already_log: bool, is_counts: bool) -> pd.DataFrame:
    expr = collapse_symbols(expr)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    if is_counts:
        log = counts_to_log2cpm(expr)
    elif already_log:
        log = expr
    else:
        log = np.log2(expr.clip(lower=0) + 1)
    return log


def attach_scores(log: pd.DataFrame) -> pd.DataFrame:
    rows = []
    immune, immune_used = mean_z(log, IMMUNE)
    stromal, stromal_used = mean_z(log, STROMAL)
    ifn, ifn_used = mean_z(log, IFN_GENES)
    mhc, mhc_used = mean_z(log, MHC_GENES)
    for s in log.columns:
        rec = {"sample": s}
        for g in TARGETS + ["CD8A"]:
            rec[g] = float(log.loc[g, s]) if g in log.index else np.nan
        rec["IFN"] = float(ifn[s]) if ifn is not None else np.nan
        rec["MHC"] = float(mhc[s]) if mhc is not None else np.nan
        rec["ImmuneScore"] = float(immune[s]) if immune is not None else np.nan
        rec["StromalScore"] = float(stromal[s]) if stromal is not None else np.nan
        if immune is not None and stromal is not None:
            rec["ESTIMATEScore"] = rec["ImmuneScore"] + rec["StromalScore"]
        else:
            rec["ESTIMATEScore"] = np.nan
        rows.append(rec)
    out = pd.DataFrame(rows)
    if out["ImmuneScore"].notna().sum() >= 4:
        out["CLDN4_resid_Immune"] = residualize(out["CLDN4"], out["ImmuneScore"])
        out["TACSTD2_resid_Immune"] = residualize(out["TACSTD2"], out["ImmuneScore"])
    else:
        out["CLDN4_resid_Immune"] = np.nan
        out["TACSTD2_resid_Immune"] = np.nan
    out.attrs["immune_n"] = len(immune_used)
    out.attrs["stromal_n"] = len(stromal_used)
    out.attrs["ifn_n"] = len(ifn_used)
    out.attrs["mhc_n"] = len(mhc_used)
    return out


def box_strip(ax, groups, labels, colors):
    data = [np.asarray(g, float) for g in groups]
    data = [d[~np.isnan(d)] for d in data]
    bp = ax.boxplot(data, tick_labels=labels, widths=0.55, showfliers=False, patch_artist=True)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.85)
    rng = np.random.default_rng(0)
    for i, d in enumerate(data, start=1):
        ax.scatter(i + rng.uniform(-0.08, 0.08, size=len(d)), d, s=22, c="k", alpha=0.7, zorder=3)


def save_response_box(path, gene, r, nr, title, ylab):
    fig, ax = plt.subplots(figsize=(3.6, 4.0))
    box_strip(ax, [nr, r], [f"NR\nn={np.isfinite(nr).sum()}", f"R\nn={np.isfinite(r).sum()}"], ["#c9d6df", "#f6b26b"])
    ax.set_ylabel(ylab)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def save_scatter(path, x, y, title, xlab, ylab):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    fig, ax = plt.subplots(figsize=(4.0, 4.0))
    ax.scatter(x[m], y[m], s=28, c="#3d5a80", alpha=0.85)
    if m.sum() >= 3:
        lr = stats.linregress(x[m], y[m])
        xs = np.linspace(np.nanmin(x[m]), np.nanmax(x[m]), 50)
        ax.plot(xs, lr.intercept + lr.slope * xs, color="#e07a5f", lw=1.4)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


ALL_TESTS = []
PRIMARY = []
INVENTORY = []


def add_test(**kw):
    ALL_TESTS.append(kw)


def add_primary(**kw):
    PRIMARY.append(kw)


def add_inv(**kw):
    INVENTORY.append(kw)


def immune_tests(cohort, df, note=""):
    for gene in ["CLDN4", "TACSTD2"]:
        for axis in ["CD8A", "IFN", "MHC", "ImmuneScore"]:
            if gene not in df or axis not in df:
                continue
            sp = spearman(df[gene], df[axis])
            add_test(
                cohort=cohort,
                gene=gene,
                endpoint=f"vs_{axis}",
                n=sp["n"],
                n_R="",
                n_NR="",
                metric="Spearman_rho",
                effect=sp["rho"],
                p=sp["p"],
                note=note,
            )
        if "CLDN4_resid_Immune" in df and gene == "CLDN4":
            for axis in ["CD8A", "IFN", "MHC"]:
                if axis not in df or df["ImmuneScore"].notna().sum() < 5:
                    continue
                psp = partial_spearman(df[gene], df[axis], df["ImmuneScore"])
                add_test(
                    cohort=cohort,
                    gene=gene,
                    endpoint=f"vs_{axis}_partial_ImmuneScore",
                    n=psp["n"],
                    n_R="",
                    n_NR="",
                    metric="partial_Spearman_rho",
                    effect=psp["rho"],
                    p=psp["p"],
                    note="residual/partial after ESTIMATE ImmuneScore",
                )


def binary_gene_tests(cohort, df, mask_r, mask_nr, endpoint, note=""):
    for gene in ["CLDN4", "TACSTD2"]:
        if gene not in df:
            continue
        blk = mwu_block(df.loc[mask_r, gene], df.loc[mask_nr, gene])
        add_test(
            cohort=cohort,
            gene=gene,
            endpoint=endpoint,
            n=blk["n_R"] + blk["n_NR"],
            n_R=blk["n_R"],
            n_NR=blk["n_NR"],
            metric="cliffs_delta_R_minus_NR",
            effect=blk["cliffs_delta"],
            p=blk["p"],
            note=f"{note}; AUC={blk['auc_R_gt_NR']}; medR={blk['median_R']}; medNR={blk['median_NR']}",
        )
        return_blk = blk if gene == "CLDN4" else None
    return mwu_block(df.loc[mask_r, "CLDN4"], df.loc[mask_nr, "CLDN4"]) if "CLDN4" in df else None


# ---------------------------------------------------------------------------
# GSE126044 — NSCLC tumor, anti-PD-1, R vs NR
# ---------------------------------------------------------------------------
meta = parse_series(DATA / "series/GSE126044_series_matrix.txt.gz")
expr = pd.read_csv(DATA / "expr/GSE126044_counts.txt.gz", sep="\t", index_col=0)
log = score_matrix(expr, already_log=False, is_counts=True)
sc = attach_scores(log)
meta["sample"] = meta["title"].str.replace("RNA-seq_", "", regex=False)
df = meta.merge(sc, on="sample", how="inner")
df["patient"] = df["gsm"]
df["resp"] = df["patient response"].str.lower()
assert df["patient"].is_unique
df.to_csv(PROC / "GSE126044_patient.tsv", sep="\t", index=False)
r = df["resp"].eq("responder")
nr = df["resp"].eq("non-responder")
blk = binary_gene_tests("GSE126044", df, r, nr, "ORR_R_vs_NR", "tumor bulk counts; GEO patient response")
immune_tests("GSE126044", df, "ESTIMATE lists used")
add_inv(
    cohort="GSE126044",
    tissue="NSCLC tumor biopsy",
    assay="RNA-seq counts",
    n_geo=16,
    n_patients=int(df["patient"].nunique()),
    CLDN4="yes",
    TACSTD2="yes",
    ICI_label="responder / non-responder (GEO)",
    estimate_immune_genes=sc.attrs.get("immune_n"),
)
add_primary(
    cohort="GSE126044",
    tissue="NSCLC tumor",
    n=int(r.sum() + nr.sum()),
    n_R=int(r.sum()),
    n_NR=int(nr.sum()),
    endpoint="ORR (GEO responder vs non-responder)",
    gene="CLDN4",
    metric="Cliff δ (R−NR)",
    effect=blk["cliffs_delta"],
    p=blk["p"],
    note=f"log2CPM; med R {blk['median_R']:.2f} vs NR {blk['median_NR']:.2f}; TACSTD2 companion in all_tests",
)
save_response_box(
    FIG / "GSE126044_CLDN4_ORR.png",
    "CLDN4",
    df.loc[r, "CLDN4"],
    df.loc[nr, "CLDN4"],
    "GSE126044 CLDN4 vs anti-PD-1 response",
    "CLDN4 log2(CPM+1)",
)
save_scatter(
    FIG / "GSE126044_CLDN4_vs_IFN.png",
    df["CLDN4"],
    df["IFN"],
    "GSE126044 CLDN4 vs IFN (Ayers-6 z)",
    "CLDN4 log2(CPM+1)",
    "IFN mean-z",
)
save_scatter(
    FIG / "GSE126044_CLDN4_vs_MHC.png",
    df["CLDN4"],
    df["MHC"],
    "GSE126044 CLDN4 vs MHC-I mean-z",
    "CLDN4 log2(CPM+1)",
    "MHC-I mean-z",
)
save_scatter(
    FIG / "GSE126044_CLDN4_vs_CD8A.png",
    df["CLDN4"],
    df["CD8A"],
    "GSE126044 CLDN4 vs CD8A",
    "CLDN4 log2(CPM+1)",
    "CD8A log2(CPM+1)",
)

# ---------------------------------------------------------------------------
# GSE136961 — Oncomine Immune Response; DCB in titles (D vs N)
# ---------------------------------------------------------------------------
meta = parse_series(DATA / "series/GSE136961_series_matrix.txt.gz")
expr = pd.read_csv(DATA / "expr/GSE136961_TPM.tsv.gz", sep="\t")
expr = expr.set_index(expr["Symbol_ID"].astype(str).str.split("_").str[0])
expr = expr.drop(columns=["Symbol_ID"])
log = score_matrix(expr, already_log=False, is_counts=False)
sc = attach_scores(log)
meta["sample"] = meta["title"]
df = meta.merge(sc, on="sample", how="inner")
df["patient"] = df["title"]
df["DCB"] = np.where(df["title"].str.startswith("D"), "DCB", np.where(df["title"].str.startswith("N"), "NDB", None))


def parse_pfs(s):
    if not isinstance(s, str):
        return np.nan, np.nan
    pfs = np.nan
    rec = np.nan
    for part in s.split(","):
        part = part.strip()
        if part.startswith("progression-free-survival_days"):
            pfs = float(part.split(":")[1])
        if part.startswith("recurrence"):
            rec = float(part.split(":")[1])
    return pfs, rec


pfs_os = df["survival"].map(lambda s: parse_pfs(s) if isinstance(s, str) else (np.nan, np.nan))
df["pfs_days"] = [a for a, _ in pfs_os]
df["pfs_event"] = [b for _, b in pfs_os]
df.to_csv(PROC / "GSE136961_patient.tsv", sep="\t", index=False)
add_inv(
    cohort="GSE136961",
    tissue="NSCLC tumor (Oncomine Immune Response 395)",
    assay="targeted panel TPM",
    n_geo=21,
    n_patients=21,
    CLDN4="no",
    TACSTD2="no",
    ICI_label="DCB vs NDB in titles (D/N) + PFS days",
    estimate_immune_genes=sc.attrs.get("immune_n"),
)
add_primary(
    cohort="GSE136961",
    tissue="NSCLC tumor (Oncomine 395)",
    n=21,
    n_R=int((df["DCB"] == "DCB").sum()),
    n_NR=int((df["DCB"] == "NDB").sum()),
    endpoint="DCB (title D vs N); PFS in GEO",
    gene="CLDN4",
    metric="not on panel",
    effect=np.nan,
    p=np.nan,
    note="Oncomine Immune Response panel has CD8A/IFNG/HLA-A/B/C; CLDN4 and TACSTD2 absent. No CLDN4–DCB test.",
)
# immune-only DCB check is not the question; do not invent CLDN4

# ---------------------------------------------------------------------------
# GSE166449 — LUAD pembro, titles Responder / nonResponder
# ---------------------------------------------------------------------------
meta = parse_series(DATA / "series/GSE166449_series_matrix.txt.gz")
expr = pd.read_csv(DATA / "expr/GSE166449_Raw_gene_TPM_matrix.txt.gz", sep="\t", index_col=0)
log = score_matrix(expr, already_log=False, is_counts=False)
sc = attach_scores(log)
meta["sample"] = meta["description"]
meta["resp"] = np.where(
    meta["title"].str.contains("nonResponder", case=False),
    "NR",
    np.where(meta["title"].str.contains("Responder", case=False), "R", None),
)
df = meta.merge(sc, on="sample", how="inner")
df["patient"] = df["gsm"]
assert df["patient"].is_unique
df.to_csv(PROC / "GSE166449_patient.tsv", sep="\t", index=False)
r = df["resp"].eq("R")
nr = df["resp"].eq("NR")
blk = binary_gene_tests("GSE166449", df, r, nr, "ORR_R_vs_NR", "tumor TPM; GEO titles Responder/nonResponder")
immune_tests("GSE166449", df)
add_inv(
    cohort="GSE166449",
    tissue="advanced LUAD tumor",
    assay="RNA-seq TPM",
    n_geo=22,
    n_patients=22,
    CLDN4="yes",
    TACSTD2="yes",
    ICI_label="Responder 7 / nonResponder 15 (titles)",
    estimate_immune_genes=sc.attrs.get("immune_n"),
)
add_primary(
    cohort="GSE166449",
    tissue="LUAD tumor (pembro)",
    n=22,
    n_R=int(r.sum()),
    n_NR=int(nr.sum()),
    endpoint="ORR (GEO title Responder vs nonResponder)",
    gene="CLDN4",
    metric="Cliff δ (R−NR)",
    effect=blk["cliffs_delta"],
    p=blk["p"],
    note=f"log2(TPM+1); med R {blk['median_R']:.2f} vs NR {blk['median_NR']:.2f}",
)
save_response_box(
    FIG / "GSE166449_CLDN4_ORR.png",
    "CLDN4",
    df.loc[r, "CLDN4"],
    df.loc[nr, "CLDN4"],
    "GSE166449 CLDN4 vs pembrolizumab response",
    "CLDN4 log2(TPM+1)",
)
save_scatter(FIG / "GSE166449_CLDN4_vs_IFN.png", df["CLDN4"], df["IFN"], "GSE166449 CLDN4 vs IFN", "CLDN4", "IFN mean-z")
save_scatter(FIG / "GSE166449_CLDN4_vs_MHC.png", df["CLDN4"], df["MHC"], "GSE166449 CLDN4 vs MHC-I", "CLDN4", "MHC-I mean-z")

# ---------------------------------------------------------------------------
# GSE182328 — tumor RNA-seq; GEO has Akkermansia only (not RECIST/PFS)
# ---------------------------------------------------------------------------
meta = parse_series(DATA / "series/GSE182328_series_matrix.txt.gz")
expr = pd.read_csv(DATA / "expr/GSE182328_Gene_counts_matrix.txt.gz", sep="\t", index_col=0)
log = score_matrix(expr, already_log=False, is_counts=True)
sc = attach_scores(log)
meta["sample"] = meta["title"]
df = meta.merge(sc, on="sample", how="inner")
df["patient"] = df["title"]
df.to_csv(PROC / "GSE182328_patient.tsv", sep="\t", index=False)
# deposited grouping: Akkermansia detectable = better-prognosis surrogate in the paper, not a GEO RECIST field
akk_pos = df["akk_metaominer"].eq("detectable")
akk_neg = df["akk_metaominer"].eq("not_detectable")
# treat detectable as "better" group for direction (paper: Akk+ better ICI)
blk = binary_gene_tests(
    "GSE182328",
    df,
    akk_pos,
    akk_neg,
    "Akkermansia_detectable_vs_not",
    "NOT DCB/ORR/PFS/MPR in GEO; Akkermansia grouping only",
)
immune_tests("GSE182328", df)
add_inv(
    cohort="GSE182328",
    tissue="lung tumor",
    assay="RNA-seq counts",
    n_geo=44,
    n_patients=44,
    CLDN4="yes",
    TACSTD2="yes",
    ICI_label="none in GEO (Akkermansia detectable/not only)",
    estimate_immune_genes=sc.attrs.get("immune_n"),
)
add_primary(
    cohort="GSE182328",
    tissue="NSCLC tumor (ICI-treated; GEO has no RECIST/PFS)",
    n=44,
    n_R=int(akk_pos.sum()),
    n_NR=int(akk_neg.sum()),
    endpoint="Akkermansia detectable vs not (GEO surrogate; not DCB/ORR/PFS/MPR)",
    gene="CLDN4",
    metric="Cliff δ (Akk+ − Akk−)",
    effect=blk["cliffs_delta"],
    p=blk["p"],
    note=f"log2CPM; med Akk+ {blk['median_R']:.2f} vs Akk− {blk['median_NR']:.2f}. No DCB/MPR/PFS/ORR field on GEO.",
)
save_response_box(
    FIG / "GSE182328_CLDN4_Akkermansia.png",
    "CLDN4",
    df.loc[akk_pos, "CLDN4"],
    df.loc[akk_neg, "CLDN4"],
    "GSE182328 CLDN4 vs Akkermansia (GEO grouping)",
    "CLDN4 log2(CPM+1)",
)
save_scatter(FIG / "GSE182328_CLDN4_vs_IFN.png", df["CLDN4"], df["IFN"], "GSE182328 CLDN4 vs IFN", "CLDN4", "IFN mean-z")
save_scatter(FIG / "GSE182328_CLDN4_vs_MHC.png", df["CLDN4"], df["MHC"], "GSE182328 CLDN4 vs MHC-I", "CLDN4", "MHC-I mean-z")
save_scatter(FIG / "GSE182328_CLDN4_vs_CD8A.png", df["CLDN4"], df["CD8A"], "GSE182328 CLDN4 vs CD8A", "CLDN4", "CD8A")

# ---------------------------------------------------------------------------
# GSE161537 — HTG OBP; ICI RECIST+PFS; CLDN4/TACSTD2 absent
# ---------------------------------------------------------------------------
meta = parse_series(DATA / "series/GSE161537_series_matrix.txt.gz")
raw = pd.read_csv(DATA / "expr/GSE161537_nivobio_log2cpm.csv.gz", sep=";")
raw = raw.set_index(raw.columns[0])
raw = raw.apply(lambda col: pd.to_numeric(col.astype(str).str.replace(",", ".", regex=False), errors="coerce"))
log = score_matrix(raw, already_log=True, is_counts=False)
sc = attach_scores(log)
meta["sample"] = meta["patient id"].astype(str)
sc["sample"] = sc["sample"].astype(str)
df = meta.merge(sc, on="sample", how="inner")
df["patient"] = df["patient id"]
df["recist"] = df["best response on immunotherapy (recist)"]
df.to_csv(PROC / "GSE161537_patient.tsv", sep="\t", index=False)
n_eval = int(df["recist"].isin(["CR", "PR", "SD", "PD"]).sum())
add_inv(
    cohort="GSE161537",
    tissue="advanced NSCLC tumor (HTG OBP ~2560)",
    assay="HTG EdgeSeq log2CPM",
    n_geo=82,
    n_patients=int(df["patient"].nunique()),
    CLDN4="no",
    TACSTD2="no",
    ICI_label="RECIST + PFS/OS (2L PD-1/PD-L1)",
    estimate_immune_genes=sc.attrs.get("immune_n"),
)
add_primary(
    cohort="GSE161537",
    tissue="advanced NSCLC (HTG OBP)",
    n=n_eval,
    n_R=int(df["recist"].isin(["CR", "PR"]).sum()),
    n_NR=int(df["recist"].isin(["SD", "PD"]).sum()),
    endpoint="ORR (GEO RECIST CR/PR vs SD/PD); PFS present",
    gene="CLDN4",
    metric="not on panel",
    effect=np.nan,
    p=np.nan,
    note="HTG Oncology Biomarker Panel has CLDN3/EPCAM/CD8A/IFNG, not CLDN4 or TACSTD2. RECIST n=76 evaluable (6 NA).",
)

# ---------------------------------------------------------------------------
# GSE162520 — HTG OBP; surgically treated early NSCLC; NOT ICI
# ---------------------------------------------------------------------------
meta = parse_series(DATA / "series/GSE162520_series_matrix.txt.gz")
raw = pd.read_csv(DATA / "expr/GSE162520_GEO_data_TUMADOR_log2cpm.csv.gz", sep=";")
raw = raw.set_index(raw.columns[0])
raw = raw.apply(lambda col: pd.to_numeric(col.astype(str).str.replace(",", ".", regex=False), errors="coerce"))
log = score_matrix(raw, already_log=True, is_counts=False)
sc = attach_scores(log)
meta["sample"] = meta["description"]
df = meta.merge(sc, on="sample", how="inner")
df["patient"] = df["description"]
df.to_csv(PROC / "GSE162520_patient.tsv", sep="\t", index=False)
add_inv(
    cohort="GSE162520",
    tissue="early NSCLC tumor (surgical, HTG OBP)",
    assay="HTG EdgeSeq log2CPM",
    n_geo=92,
    n_patients=92,
    CLDN4="no",
    TACSTD2="no",
    ICI_label="none — surgery cohort; PFS/OS are not ICI endpoints",
    estimate_immune_genes=sc.attrs.get("immune_n"),
)
add_primary(
    cohort="GSE162520",
    tissue="early NSCLC (surgery, not ICI)",
    n=92,
    n_R="",
    n_NR="",
    endpoint="no ICI DCB/MPR/PFS/ORR (surgical PFS/OS only)",
    gene="CLDN4",
    metric="not on panel; no ICI label",
    effect=np.nan,
    p=np.nan,
    note="GEO overall design: 92 surgically treated NSCLC. HTG OBP lacks CLDN4/TACSTD2. Not an ICI-response cohort.",
)

# ---------------------------------------------------------------------------
# GSE93157 — NanoString 730; NSCLC subset; CLDN4/TACSTD2 absent
# ---------------------------------------------------------------------------
meta = parse_series(DATA / "series/GSE93157_series_matrix.txt.gz")
raw = pd.read_csv(DATA / "expr/GSE93157_raw_data_values.txt.gz", sep="\t", header=None)
hdr_i = int(raw.index[raw.iloc[:, 0].astype(str).eq("ID_REF")][0])
n_samp = len(meta)
mat = raw.iloc[hdr_i + 1 :, : n_samp + 1].copy()
mat.columns = ["gene"] + [f"SAMPLE_{i}" for i in range(1, n_samp + 1)]
mat = mat.set_index("gene")
mat = mat.apply(pd.to_numeric, errors="coerce")
log = score_matrix(mat, already_log=False, is_counts=False)
sc = attach_scores(log)
meta = meta.reset_index(drop=True)
meta["sample"] = [f"SAMPLE_{i}" for i in range(1, n_samp + 1)]
df = meta.merge(sc, on="sample", how="inner")
df["nsclc"] = df["source"].str.contains("LUNG", case=False, na=False)
lung = df[df["nsclc"]].copy()
lung["patient"] = lung["title"]
lung["orr"] = np.where(lung["best.resp"].isin(["CR", "PR"]), "R", np.where(lung["best.resp"].isin(["SD", "PD"]), "NR", None))
lung["dcb"] = np.where(lung["best.resp"].isin(["CR", "PR", "SD"]), "DCB", np.where(lung["best.resp"].eq("PD"), "NDB", None))
lung["pfs"] = pd.to_numeric(lung["pfs"], errors="coerce")
lung["pfse"] = pd.to_numeric(lung["pfse"], errors="coerce")
lung.to_csv(PROC / "GSE93157_NSCLC_patient.tsv", sep="\t", index=False)
add_inv(
    cohort="GSE93157_NSCLC",
    tissue="NSCLC tumor (nCounter PanCancer Immune 730)",
    assay="NanoString",
    n_geo=65,
    n_patients=int(lung["patient"].nunique()),
    CLDN4="no",
    TACSTD2="no",
    ICI_label="best.resp + PFS (nivo/pembro)",
    estimate_immune_genes=sc.attrs.get("immune_n"),
)
add_primary(
    cohort="GSE93157_NSCLC",
    tissue="NSCLC tumor (NanoString 730)",
    n=int(lung.shape[0]),
    n_R=int((lung["orr"] == "R").sum()),
    n_NR=int((lung["orr"] == "NR").sum()),
    endpoint="ORR CR/PR vs SD/PD; DCB and PFS also in GEO",
    gene="CLDN4",
    metric="not on panel",
    effect=np.nan,
    p=np.nan,
    note=f"NSCLC subset {int(lung.shape[0])} (22 non-squamous + 13 squamous). Panel has CD8A/IFNG/HLA; no CLDN4/TACSTD2.",
)

# ---------------------------------------------------------------------------
# GSE280232 — 10x scRNA RAW.tar only; no processed bulk; no MPR/RECIST
# ---------------------------------------------------------------------------
add_inv(
    cohort="GSE280232",
    tissue="lung tumor / adjacent (10x GEX+TCR)",
    assay="10x MTX in RAW.tar only",
    n_geo=42,
    n_patients="",
    CLDN4="not scored (no processed bulk)",
    TACSTD2="not scored (no processed bulk)",
    ICI_label="none — genotype (KRAS/STK11) + regimen only",
    estimate_immune_genes="",
)
add_primary(
    cohort="GSE280232",
    tissue="resectable KRAS-mutant lung (neoadjuvant ICI)",
    n="",
    n_R="",
    n_NR="",
    endpoint="no MPR/ORR/PFS/DCB in GEO",
    gene="CLDN4",
    metric="no processed bulk matrix",
    effect=np.nan,
    p=np.nan,
    note="GEO FTP has RAW.tar 10x MTX/TCR only (no counts/TPM table). Characteristics: treatment + KRAS/STK11, no response field. Public processed bulk not available.",
)

# ---------------------------------------------------------------------------
# GSE283829 — NSCLC tumor counts; CR / SD / PD
# ---------------------------------------------------------------------------
meta = parse_series(DATA / "series/GSE283829_series_matrix.txt.gz")
expr = pd.read_csv(DATA / "expr/GSE283829_raw_express_matrix_all_samples.txt.gz", sep="\t", index_col=0)
log = score_matrix(expr, already_log=False, is_counts=True)
sc = attach_scores(log)
meta["sample"] = "S" + meta["title"].astype(str)
df = meta.merge(sc, on="sample", how="inner")
df["patient"] = df["title"]
df["recist"] = df["disease stage"]
assert df["patient"].is_unique
df.to_csv(PROC / "GSE283829_patient.tsv", sep="\t", index=False)
cr = df["recist"].eq("CR")
pd_ = df["recist"].eq("PD")
dcr = df["recist"].isin(["CR", "SD"])
blk_crpd = binary_gene_tests("GSE283829", df, cr, pd_, "CR_vs_PD", "tumor counts; GEO disease stage")
blk_dcr = binary_gene_tests("GSE283829", df, dcr, pd_, "DCR_CRorSD_vs_PD", "tumor counts")
immune_tests("GSE283829", df)
add_inv(
    cohort="GSE283829",
    tissue="NSCLC tumor",
    assay="RNA-seq counts (Ensembl)",
    n_geo=27,
    n_patients=27,
    CLDN4="yes",
    TACSTD2="yes",
    ICI_label="disease stage CR/SD/PD",
    estimate_immune_genes=sc.attrs.get("immune_n"),
)
add_primary(
    cohort="GSE283829",
    tissue="NSCLC tumor",
    n=int(cr.sum() + pd_.sum()),
    n_R=int(cr.sum()),
    n_NR=int(pd_.sum()),
    endpoint="ORR-extreme CR vs PD (GEO disease stage); DCR also tested",
    gene="CLDN4",
    metric="Cliff δ (CR−PD)",
    effect=blk_crpd["cliffs_delta"],
    p=blk_crpd["p"],
    note=f"log2CPM; CR n={int(cr.sum())} PD n={int(pd_.sum())} SD n={int(df['recist'].eq('SD').sum())}; DCR δ={blk_dcr['cliffs_delta']:.3f} p={blk_dcr['p']:.3g}",
)
save_response_box(
    FIG / "GSE283829_CLDN4_CR_vs_PD.png",
    "CLDN4",
    df.loc[cr, "CLDN4"],
    df.loc[pd_, "CLDN4"],
    "GSE283829 CLDN4 CR vs PD",
    "CLDN4 log2(CPM+1)",
)
save_scatter(FIG / "GSE283829_CLDN4_vs_IFN.png", df["CLDN4"], df["IFN"], "GSE283829 CLDN4 vs IFN", "CLDN4", "IFN mean-z")
save_scatter(FIG / "GSE283829_CLDN4_vs_MHC.png", df["CLDN4"], df["MHC"], "GSE283829 CLDN4 vs MHC-I", "CLDN4", "MHC-I mean-z")
save_scatter(FIG / "GSE283829_CLDN4_vs_CD8A.png", df["CLDN4"], df["CD8A"], "GSE283829 CLDN4 vs CD8A", "CLDN4", "CD8A")

# ---------------------------------------------------------------------------
# GSE260770 — peripheral blood FPKM; sintilimab response
# ---------------------------------------------------------------------------
meta = parse_series(DATA / "series/GSE260770_series_matrix.txt.gz")
expr = pd.read_csv(DATA / "expr/GSE260770_mRNA_FPKM.txt.gz", sep="\t")
expr = expr.set_index("Symbol")
expr = expr.drop(columns=[c for c in expr.columns if c.startswith("#") or c == "Symbol"], errors="ignore")
log = score_matrix(expr, already_log=False, is_counts=False)
sc = attach_scores(log)
meta["sample"] = meta["title"]
df = meta.merge(sc, on="sample", how="inner")
df["patient"] = df["title"]
df["resp"] = np.where(
    df["group"].str.contains("Non-responsed", case=False, na=False),
    "NR",
    np.where(df["group"].str.contains("Responsed", case=False, na=False), "R", None),
)
df.to_csv(PROC / "GSE260770_patient.tsv", sep="\t", index=False)
r = df["resp"].eq("R")
nr = df["resp"].eq("NR")
blk = binary_gene_tests(
    "GSE260770",
    df,
    r,
    nr,
    "sintilimab_R_vs_NR",
    "peripheral blood FPKM; epithelial genes at floor (not a tumor test)",
)
immune_tests("GSE260770", df, "blood; CLDN4 median 0")
add_inv(
    cohort="GSE260770",
    tissue="peripheral blood (not tumor)",
    assay="RNA-seq FPKM",
    n_geo=50,
    n_patients=50,
    CLDN4="yes (floor)",
    TACSTD2="yes (floor)",
    ICI_label="Responsed / Non-responsed to sintilimab",
    estimate_immune_genes=sc.attrs.get("immune_n"),
)
add_primary(
    cohort="GSE260770",
    tissue="peripheral blood (sintilimab, GGO/MPLC)",
    n=50,
    n_R=int(r.sum()),
    n_NR=int(nr.sum()),
    endpoint="GEO group Responsed vs Non-responsed",
    gene="CLDN4",
    metric="Cliff δ (R−NR)",
    effect=blk["cliffs_delta"],
    p=blk["p"],
    note=f"blood FPKM; CLDN4 nonzero {int((df['CLDN4']>0).sum())}/50 after log2(FPKM+1) still near floor. Not tumor epithelium.",
)
save_response_box(
    FIG / "GSE260770_CLDN4_blood_response.png",
    "CLDN4",
    df.loc[r, "CLDN4"],
    df.loc[nr, "CLDN4"],
    "GSE260770 blood CLDN4 vs sintilimab (floor)",
    "CLDN4 log2(FPKM+1)",
)

# ---------------------------------------------------------------------------
# GSE293591 — pan-cancer TPM vs IHC; no ICI labels. Lung subset immune axis only.
# ---------------------------------------------------------------------------
meta = parse_series(DATA / "series/GSE293591_series_matrix.txt.gz")
wanted = set(TARGETS + ["CD8A"] + IFN_GENES + MHC_GENES + IMMUNE + STROMAL)
chunks = []
for chunk in pd.read_csv(DATA / "expr/GSE293591_TPM_all_samples.tsv.gz", sep="\t", index_col=0, chunksize=4000):
    keep = [i for i in chunk.index if str(i) in wanted]
    if keep:
        chunks.append(chunk.loc[keep])
expr = pd.concat(chunks) if chunks else pd.DataFrame()
log = score_matrix(expr, already_log=False, is_counts=False)
sc = attach_scores(log)
meta["sample"] = meta["title"].str.split(",").str[0].str.replace("Sample_", "Sample_", regex=False)
# titles are "Sample_001, Breast cancer, biopsy"
meta["sample"] = meta["title"].str.extract(r"(Sample_\d+)", expand=False)
df = meta.merge(sc, on="sample", how="inner")
df["patient"] = df["sample"]
lung = df[df["diagnosis"].isin(["Lung Adenocarcinoma", "Squamous Cell Carcinoma of Lung"])].copy()
lung.to_csv(PROC / "GSE293591_lung_patient.tsv", sep="\t", index=False)
df.to_csv(PROC / "GSE293591_all_patient.tsv", sep="\t", index=False)
immune_tests("GSE293591_lung", lung, "no ICI labels; lung tumors only")
add_inv(
    cohort="GSE293591",
    tissue="FFPE/FF solid tumor (pan-cancer)",
    assay="RNA-seq TPM",
    n_geo=365,
    n_patients=int(lung["patient"].nunique()) if len(lung) else 0,
    CLDN4="yes",
    TACSTD2="yes",
    ICI_label="none",
    estimate_immune_genes=sc.attrs.get("immune_n"),
)
sp_ifn = spearman(lung["CLDN4"], lung["IFN"]) if len(lung) else {"n": 0, "rho": np.nan, "p": np.nan}
sp_mhc = spearman(lung["CLDN4"], lung["MHC"]) if len(lung) else {"n": 0, "rho": np.nan, "p": np.nan}
add_primary(
    cohort="GSE293591_lung",
    tissue="LUAD+LUSC tumor (no ICI)",
    n=int(len(lung)),
    n_R="",
    n_NR="",
    endpoint="no ICI DCB/MPR/PFS/ORR; CLDN4 vs IFN (immune axis)",
    gene="CLDN4",
    metric="Spearman ρ vs IFN",
    effect=sp_ifn["rho"],
    p=sp_ifn["p"],
    note=f"Lung ADC n={(lung['diagnosis']=='Lung Adenocarcinoma').sum() if len(lung) else 0}, LUSC n={(lung['diagnosis']=='Squamous Cell Carcinoma of Lung').sum() if len(lung) else 0}. vs MHC ρ={sp_mhc['rho']} p={sp_mhc['p']}. No ICI field on GEO.",
)
if len(lung):
    save_scatter(
        FIG / "GSE293591_lung_CLDN4_vs_IFN.png",
        lung["CLDN4"],
        lung["IFN"],
        "GSE293591 lung CLDN4 vs IFN (no ICI labels)",
        "CLDN4 log2(TPM+1)",
        "IFN mean-z",
    )
    save_scatter(
        FIG / "GSE293591_lung_CLDN4_vs_MHC.png",
        lung["CLDN4"],
        lung["MHC"],
        "GSE293591 lung CLDN4 vs MHC-I (no ICI labels)",
        "CLDN4 log2(TPM+1)",
        "MHC-I mean-z",
    )

# ---------------------------------------------------------------------------
# write tables
# ---------------------------------------------------------------------------
prim = pd.DataFrame(PRIMARY)
inv = pd.DataFrame(INVENTORY)
tests = pd.DataFrame(ALL_TESTS)
prim.to_csv(TAB / "cohort_primary.tsv", sep="\t", index=False)
inv.to_csv(TAB / "inventory.tsv", sep="\t", index=False)
tests.to_csv(TAB / "all_tests.tsv", sep="\t", index=False)

# compact JSON for FINDING
summary = {
    "primary": prim.to_dict(orient="records"),
    "inventory": inv.to_dict(orient="records"),
    "n_tests": int(len(tests)),
}
with open(TAB / "summary.json", "w") as f:
    json.dump(summary, f, indent=2, default=str)

print("PRIMARY")
print(prim.to_string(index=False))
print("\nWrote", TAB / "cohort_primary.tsv", "n_tests", len(tests))

# Extra 4-panel: barrier / ICI-poor vs CLDN4-high MHC-I (rank partial)
a = pd.read_csv(PROC / "GSE126044_patient.tsv", sep="\t")
b = pd.read_csv(PROC / "GSE182328_patient.tsv", sep="\t")
c = pd.read_csv(PROC / "GSE283829_patient.tsv", sep="\t")
d = pd.read_csv(PROC / "GSE293591_lung_patient.tsv", sep="\t")
fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.6))
r = a["resp"].eq("responder")
nr = a["resp"].eq("non-responder")
box_strip(
    axes[0, 0],
    [a.loc[nr, "CLDN4"], a.loc[r, "CLDN4"]],
    [f"NR n={int(nr.sum())}", f"R n={int(r.sum())}"],
    ["#c9d6df", "#f6b26b"],
)
axes[0, 0].set_title("A  GSE126044 CLDN4 vs anti-PD-1 ORR\n(higher in NR; ICI-poor axis)", fontsize=10)
axes[0, 0].set_ylabel("CLDN4 log2(CPM+1)")
for ax, x, y, title, xlab, ylab in (
    (
        axes[0, 1],
        b["CLDN4"],
        b["CD8A"],
        "B  GSE182328 CLDN4 vs CD8A\n(immune-low / barrier axis)",
        "CLDN4 log2(CPM+1)",
        "CD8A log2(CPM+1)",
    ),
    (
        axes[1, 0],
        c["CLDN4"],
        c["ImmuneScore"],
        "C  GSE283829 CLDN4 vs ESTIMATE ImmuneScore\n(immune-low axis)",
        "CLDN4 log2(CPM+1)",
        "ImmuneScore (mean-z)",
    ),
):
    xx, yy = np.asarray(x, float), np.asarray(y, float)
    m = ~(np.isnan(xx) | np.isnan(yy))
    ax.scatter(xx[m], yy[m], s=22, c="#3d5a80", alpha=0.85)
    if m.sum() >= 3:
        lr = stats.linregress(xx[m], yy[m])
        xs = np.linspace(np.nanmin(xx[m]), np.nanmax(xx[m]), 40)
        ax.plot(xs, lr.intercept + lr.slope * xs, color="#e07a5f", lw=1.3)
        rho, p = stats.spearmanr(xx[m], yy[m])
        ax.text(0.04, 0.96, f"ρ={rho:.2f}\np={p:.3g}", transform=ax.transAxes, va="top", fontsize=8)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlab, fontsize=9)
    ax.set_ylabel(ylab, fontsize=9)
psp = partial_spearman(d["CLDN4"], d["MHC"], d["ImmuneScore"])
# plot rank residuals so the panel matches the partial Spearman
x = np.asarray(d["CLDN4"], float)
y = np.asarray(d["MHC"], float)
z = np.asarray(d["ImmuneScore"], float)
m = ~(np.isnan(x) | np.isnan(y) | np.isnan(z))
rx, ry, rz = stats.rankdata(x[m]), stats.rankdata(y[m]), stats.rankdata(z[m])
Z = np.column_stack([np.ones(m.sum()), rz])
bx, *_ = np.linalg.lstsq(Z, rx, rcond=None)
by, *_ = np.linalg.lstsq(Z, ry, rcond=None)
axes[1, 1].scatter(rx - Z @ bx, ry - Z @ by, s=22, c="#2a9d8f", alpha=0.85)
axes[1, 1].text(
    0.04,
    0.96,
    f"partial ρ={psp['rho']:.2f}\np={psp['p']:.3g}",
    transform=axes[1, 1].transAxes,
    va="top",
    fontsize=8,
)
axes[1, 1].set_title("D  GSE293591 lung CLDN4 vs MHC-I\nafter ImmuneScore (ADC+ICI axis)", fontsize=10)
axes[1, 1].set_xlabel("CLDN4 rank residual | ImmuneScore", fontsize=9)
axes[1, 1].set_ylabel("MHC-I rank residual | ImmuneScore", fontsize=9)
fig.suptitle("Leftover OPEN lung ICI bulk — two CLDN4 axes (not pooled)", fontsize=12, y=1.01)
fig.tight_layout()
fig.savefig(FIG / "extra_CLDN4_two_axes.png", dpi=170, bbox_inches="tight")
plt.close(fig)
print("Wrote", FIG / "extra_CLDN4_two_axes.png")
