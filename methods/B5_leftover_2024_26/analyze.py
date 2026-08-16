#!/usr/bin/env python3
"""Additive 2023–2026 open lung ICI leftover for claim B5.

User B5 (11-cohort CLDN4-high ICI OR=0.42) is taken as given.
This script does not re-cut those cohorts and does not pool new rows into that OR.

Must-try:
  GSE328294  H23 KRAS-mutant NSCLC ± IBI351 (n=6 FPKM; cell line, no ICI OR)
  GSE244944  SCLC ChIP-seq (no expression-vs-ICB test)
  GSE244945  SCLC cell-line bulk RNA-seq (no patient ICB labels)
  GSE244946  SCLC cell-line scRNA (no patient ICB labels)

Downloadable tumor bulk with response labels (not in the original 11):
  GSE274975  Poddubskaya 2024 NSCLC ICI, n=61 counts + Table S1 RECIST/PFS
  GSE283829  Lindberg 2025 NSCLC ICI, n=27 counts + GEO RECIST-like labels

Public GEO / EuropePMC only. No EGA/dbGaP.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
import shutil
import tarfile
import urllib.request
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "B5_leftover_2024_26"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
PROC = OUT / "processed"
CACHE = Path("/tmp/geo_b5_leftover")
CACHE.mkdir(parents=True, exist_ok=True)
for d in (TABLES, FIGS, PROC):
    d.mkdir(parents=True, exist_ok=True)

GENES = {
    "TACSTD2": "ENSG00000184292",
    "CLDN4": "ENSG00000189143",
}
# locked cutoffs; do not pick the one closest to OR=0.42
CUTOFFS = ("median", "tertile")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, min_bytes: int = 200) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size >= min_bytes:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "B5-leftover/1.0"})
    with urllib.request.urlopen(req, timeout=180) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f)
    if tmp.stat().st_size < min_bytes:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"download too small: {url}")
    tmp.replace(dest)
    return dest


def parse_soft_samples(soft_gz: Path) -> list[dict]:
    samples = []
    cur = None
    with gzip.open(soft_gz, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("^SAMPLE"):
                if cur:
                    samples.append(cur)
                cur = {"gsm": line.split("=", 1)[1].strip(), "chars": {}}
            elif cur is None:
                continue
            elif line.startswith("!Sample_title"):
                cur["title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_source_name"):
                cur["source"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_description"):
                cur["description"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics"):
                val = line.split("=", 1)[1].strip()
                if ":" in val:
                    k, v = val.split(":", 1)
                    cur["chars"][k.strip().lower()] = v.strip()
    if cur:
        samples.append(cur)
    return samples


def log2_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts.div(lib, axis=1) * 1e6 + 1.0)


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """δ > 0 if x tends to be larger than y."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if len(x) == 0 or len(y) == 0:
        return float("nan")
    u = stats.mannwhitneyu(x, y, alternative="two-sided").statistic
    return float(2.0 * u / (len(x) * len(y)) - 1.0)


def rank_biserial(x: np.ndarray, y: np.ndarray) -> float:
    """Positive if x > y (same sign as Cliff's δ)."""
    d = cliffs_delta(x, y)
    return d  # Glass rank-biserial equals Cliff's δ for two samples


def haldane_or(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    """OR for [[a,b],[c,d]] = high-R, high-NR / low-R, low-NR. Woolf CI."""
    aa, bb, cc, dd = a, b, c, d
    if min(a, b, c, d) == 0:
        aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    or_ = (aa * dd) / (bb * cc)
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    lo, hi = math.exp(math.log(or_) - 1.96 * se), math.exp(math.log(or_) + 1.96 * se)
    return float(or_), float(lo), float(hi)


def fisher_p(a: int, b: int, c: int, d: int) -> float:
    return float(stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")[1])


def assign_cutoff(values: pd.Series, how: str) -> pd.Series:
    v = values.astype(float)
    if how == "median":
        med = v.median()
        # ties at the median go to low so both arms are defined
        return pd.Series(np.where(v > med, "high", "low"), index=v.index)
    if how == "tertile":
        q1, q2 = v.quantile([1 / 3, 2 / 3])
        out = pd.Series("mid", index=v.index)
        out[v <= q1] = "low"
        out[v >= q2] = "high"
        return out
    raise ValueError(how)


def binary_tests(df: pd.DataFrame, gene: str, ycol: str, cohort: str, endpoint: str) -> list[dict]:
    rows = []
    sub = df.dropna(subset=[gene, ycol]).copy()
    y = sub[ycol].astype(int)
    x = sub[gene].astype(float)
    n_r, n_nr = int(y.sum()), int((1 - y).sum())
    if n_r < 1 or n_nr < 1:
        return [
            {
                "cohort": cohort,
                "gene": gene,
                "endpoint": endpoint,
                "cutoff": "unestimable",
                "n": int(len(sub)),
                "n_high_or_R": n_r,
                "n_low_or_NR": n_nr,
                "effect": "NA",
                "effect_value": np.nan,
                "ci_low": np.nan,
                "ci_high": np.nan,
                "p": np.nan,
                "note": "one class empty",
            }
        ]
    xr, xnr = x[y == 1].to_numpy(), x[y == 0].to_numpy()
    u = stats.mannwhitneyu(xr, xnr, alternative="two-sided")
    rows.append(
        {
            "cohort": cohort,
            "gene": gene,
            "endpoint": endpoint,
            "cutoff": "continuous",
            "n": int(len(sub)),
            "n_high_or_R": n_r,
            "n_low_or_NR": n_nr,
            "effect": "rank_biserial_R_minus_NR",
            "effect_value": rank_biserial(xr, xnr),
            "ci_low": np.nan,
            "ci_high": np.nan,
            "p": float(u.pvalue),
            "note": f"Cliff_delta={cliffs_delta(xr, xnr):.4f}; MWU_U={u.statistic:.1f}",
        }
    )
    for how in CUTOFFS:
        lab = assign_cutoff(x, how)
        use = sub.copy()
        use["_bin"] = lab.values
        if how == "tertile":
            use = use[use["_bin"].isin(["high", "low"])]
        a = int(((use["_bin"] == "high") & (use[ycol] == 1)).sum())
        b = int(((use["_bin"] == "high") & (use[ycol] == 0)).sum())
        c = int(((use["_bin"] == "low") & (use[ycol] == 1)).sum())
        d = int(((use["_bin"] == "low") & (use[ycol] == 0)).sum())
        or_, lo, hi = haldane_or(a, b, c, d)
        rows.append(
            {
                "cohort": cohort,
                "gene": gene,
                "endpoint": endpoint,
                "cutoff": how if how == "median" else "tertile_T3_vs_T1",
                "n": int(a + b + c + d),
                "n_high_or_R": a + b,
                "n_low_or_NR": c + d,
                "effect": "OR_high_vs_low",
                "effect_value": or_,
                "ci_low": lo,
                "ci_high": hi,
                "p": fisher_p(a, b, c, d),
                "note": f"table high-R/high-NR/low-R/low-NR={a}/{b}/{c}/{d}"
                + ("; Haldane-Anscombe +0.5" if min(a, b, c, d) == 0 else ""),
            }
        )
    return rows


def logrank_p(time: np.ndarray, event: np.ndarray, group: np.ndarray) -> float:
    """Two-sample log-rank p (high=1 vs low=0)."""
    t = np.asarray(time, float)
    e = np.asarray(event, int)
    g = np.asarray(group, int)
    times = np.unique(t[e == 1])
    o1 = v = 0.0
    for ti in times:
        at_risk = t >= ti
        d = (t == ti) & (e == 1)
        n = at_risk.sum()
        n1 = ((g == 1) & at_risk).sum()
        d_tot = d.sum()
        d1 = ((g == 1) & d).sum()
        if n <= 1 or n1 == 0 or n1 == n:
            continue
        e1 = d_tot * n1 / n
        o1 += d1 - e1
        v += (n1 / n) * (1 - n1 / n) * ((n - d_tot) / (n - 1)) * d_tot
    if v <= 0:
        return float("nan")
    z = o1 / math.sqrt(v)
    return float(2 * stats.norm.sf(abs(z)))


def cox_hr_simple(time: np.ndarray, event: np.ndarray, x: np.ndarray) -> tuple[float, float]:
    """One-covariate Cox HR by Newton on partial likelihood. x is the predictor."""
    t = np.asarray(time, float)
    e = np.asarray(event, int)
    z = np.asarray(x, float)
    z = (z - np.nanmean(z)) / (np.nanstd(z) + 1e-12)
    beta = 0.0
    for _ in range(40):
        ll1 = ll2 = 0.0
        for i in np.where(e == 1)[0]:
            risk = t >= t[i]
            ez = np.exp(beta * z[risk])
            s0 = ez.sum()
            s1 = (z[risk] * ez).sum()
            s2 = ((z[risk] ** 2) * ez).sum()
            ll1 += z[i] - s1 / s0
            ll2 += -(s2 / s0 - (s1 / s0) ** 2)
        if abs(ll2) < 1e-12:
            break
        step = ll1 / (-ll2) if ll2 < 0 else ll1
        beta += float(np.clip(step, -2, 2))
        if abs(step) < 1e-8:
            break
    se = math.sqrt(1 / max(-ll2, 1e-12)) if ll2 < 0 else float("nan")
    return float(math.exp(beta)), float(2 * stats.norm.sf(abs(beta / se))) if se == se else float("nan")


def boxplot_by_group(df: pd.DataFrame, gene: str, group: str, title: str, path: Path) -> None:
    cats = [c for c in ["CR", "PR", "SD", "PD", "R", "NR", "DMSO", "IBI351"] if c in set(df[group].astype(str))]
    if not cats:
        cats = sorted(df[group].astype(str).unique())
    data = [df.loc[df[group].astype(str) == c, gene].dropna().to_numpy() for c in cats]
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.boxplot(data, tick_labels=cats, widths=0.55)
    for i, ys in enumerate(data, start=1):
        ax.scatter(np.random.default_rng(0).normal(i, 0.06, size=len(ys)), ys, s=16, alpha=0.7, c="#333")
    ax.set_ylabel(gene)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------------------
# GSE274975
# ---------------------------------------------------------------------------

def analyze_gse274975() -> tuple[pd.DataFrame, list[dict]]:
    soft = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274975/soft/GSE274975_family.soft.gz",
        CACHE / "GSE274975_family.soft.gz",
    )
    counts_p = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274975/suppl/GSE274975_raw_counts.tsv.gz",
        CACHE / "GSE274975_raw_counts.tsv.gz",
    )
    supp_zip = download(
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11669362/supplementaryFiles",
        CACHE / "PMC11669362_SupplementaryFiles.zip",
        min_bytes=1000,
    )
    with zipfile.ZipFile(supp_zip) as z:
        z.extract("Table1.xlsx", CACHE / "supp")
    clin = pd.read_excel(CACHE / "supp" / "Table1.xlsx")
    clin.columns = [c.strip() for c in clin.columns]
    samples = parse_soft_samples(soft)
    meta = pd.DataFrame(
        [
            {
                "gsm": s["gsm"],
                "title": s.get("title"),
                "source": s.get("source", ""),
                "histology_geo": s.get("source", ""),
            }
            for s in samples
        ]
    )
    counts = pd.read_csv(counts_p, sep="\t", index_col=0)
    # map count columns -> Sample_ID (OB_pat_LuC_N)
    def col_to_sid(col: str) -> str | None:
        m = re.search(r"Lu[Cc][-_]?(\d+)", col)
        if not m:
            return None
        return f"OB_pat_LuC_{int(m.group(1))}"

    colmap = {c: col_to_sid(c) for c in counts.columns}
    if any(v is None for v in colmap.values()):
        raise RuntimeError(f"unmapped count columns: {[c for c,v in colmap.items() if v is None]}")
    counts = counts.rename(columns=colmap)
    if counts.columns.duplicated().any():
        raise RuntimeError("duplicate sample map")
    expr = log2_cpm(counts).loc[list(GENES.values())].T
    expr.columns = list(GENES)
    expr["Sample_ID"] = expr.index
    df = clin.merge(expr, on="Sample_ID", how="inner")
    df = df.merge(meta, left_on="Sample_ID", right_on="title", how="left")
    if len(df) != 61:
        raise RuntimeError(f"GSE274975 join n={len(df)} expected 61")
    recist_raw = df["RECIST Response"]
    recist = recist_raw.astype(str).str.upper()
    recist = recist.where(recist_raw.notna() & ~recist.isin(["NAN", "NONE", "NA", ""]), other=pd.NA)
    df["recist"] = recist
    df["orr"] = np.where(recist.isin(["CR", "PR"]), 1, np.where(recist.isin(["SD", "PD"]), 0, np.nan))
    df["orr_vs_pd"] = np.where(recist.isin(["CR", "PR"]), 1, np.where(recist.eq("PD"), 0, np.nan))
    df["dcb6"] = (pd.to_numeric(df["PFS time, months"], errors="coerce") >= 6).astype(int)
    df["pfs"] = pd.to_numeric(df["PFS time, months"], errors="coerce")
    # Table S1 "Response status available" tracks progression events (1 in PD / some SD)
    df["pfs_event"] = pd.to_numeric(df["Response status available"], errors="coerce").fillna(0).astype(int)
    df["is_sclc"] = df["Histotype"].astype(str).str.contains("small cell", case=False) | df[
        "histology_geo"
    ].astype(str).str.contains("Small cell", case=False)
    nsclc = df.loc[~df["is_sclc"]].copy()
    df.to_csv(PROC / "GSE274975_patient_level.tsv", sep="\t", index=False)
    rows = []
    for gene in GENES:
        rows += binary_tests(nsclc.dropna(subset=["orr", gene]), gene, "orr", "GSE274975_NSCLC", "ORR_CRPR_vs_SDPD")
        rows += binary_tests(nsclc.dropna(subset=["orr_vs_pd", gene]), gene, "orr_vs_pd", "GSE274975_NSCLC", "ORR_CRPR_vs_PD")
        rows += binary_tests(nsclc.dropna(subset=["dcb6", gene]), gene, "dcb6", "GSE274975_NSCLC", "DCB_PFS_ge_6mo")
        # PFS: median split log-rank + continuous Cox (per 1 SD)
        surv = nsclc.dropna(subset=[gene, "pfs", "pfs_event"])
        x = surv[gene].astype(float)
        high = (x > x.median()).astype(int)
        p_lr = logrank_p(surv["pfs"].to_numpy(), surv["pfs_event"].to_numpy(), high.to_numpy())
        hr, p_cox = cox_hr_simple(surv["pfs"].to_numpy(), surv["pfs_event"].to_numpy(), x.to_numpy())
        rows.append(
            {
                "cohort": "GSE274975_NSCLC",
                "gene": gene,
                "endpoint": "PFS_logrank_median",
                "cutoff": "median",
                "n": int(len(surv)),
                "n_high_or_R": int(high.sum()),
                "n_low_or_NR": int((1 - high).sum()),
                "effect": "logrank_p_only",
                "effect_value": np.nan,
                "ci_low": np.nan,
                "ci_high": np.nan,
                "p": p_lr,
                "note": "event=TableS1 Response_status_available; not pooled into B5 OR",
            }
        )
        rows.append(
            {
                "cohort": "GSE274975_NSCLC",
                "gene": gene,
                "endpoint": "PFS_cox_per_SD",
                "cutoff": "continuous",
                "n": int(len(surv)),
                "n_high_or_R": int(surv["pfs_event"].sum()),
                "n_low_or_NR": int((1 - surv["pfs_event"]).sum()),
                "effect": "HR_per_SD",
                "effect_value": hr,
                "ci_low": np.nan,
                "ci_high": np.nan,
                "p": p_cox,
                "note": "univariable Cox; event column as above",
            }
        )
    boxplot_by_group(
        nsclc, "CLDN4", "recist", "GSE274975 NSCLC CLDN4 log2(CPM+1) vs RECIST", FIGS / "GSE274975_CLDN4_RECIST.png"
    )
    boxplot_by_group(
        nsclc, "TACSTD2", "recist", "GSE274975 NSCLC TACSTD2 log2(CPM+1) vs RECIST", FIGS / "GSE274975_TACSTD2_RECIST.png"
    )
    return df, rows


# ---------------------------------------------------------------------------
# GSE283829
# ---------------------------------------------------------------------------

def analyze_gse283829() -> tuple[pd.DataFrame, list[dict]]:
    soft = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/GSE283829/soft/GSE283829_family.soft.gz",
        CACHE / "GSE283829_family.soft.gz",
    )
    counts_p = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/GSE283829/suppl/GSE283829_raw_express_matrix_all_samples.txt.gz",
        CACHE / "GSE283829_raw_express_matrix_all_samples.txt.gz",
    )
    samples = parse_soft_samples(soft)
    meta = pd.DataFrame(
        [
            {
                "title": s["title"],
                "recist": s["chars"].get("disease stage"),
                "tumor_type": s["chars"].get("tumor type"),
                "pla": s["chars"].get("pla"),
                "sex": s["chars"].get("sex"),
            }
            for s in samples
        ]
    )
    counts = pd.read_csv(counts_p, sep="\t", index_col=0)
    counts.columns = [re.sub(r"^S", "", c.strip().strip('"')) for c in counts.columns]
    expr = log2_cpm(counts).loc[list(GENES.values())].T
    expr.columns = list(GENES)
    expr["title"] = expr.index
    df = meta.merge(expr, on="title", how="inner")
    if len(df) != 27:
        raise RuntimeError(f"GSE283829 join n={len(df)} expected 27")
    recist = df["recist"].astype(str).str.upper()
    df["recist"] = recist
    df["orr"] = recist.eq("CR").astype(int)  # no PR labels in GEO
    df["cr_vs_pd"] = np.where(recist.eq("CR"), 1, np.where(recist.eq("PD"), 0, np.nan))
    df["dcr"] = recist.isin(["CR", "SD"]).astype(int)
    df.to_csv(PROC / "GSE283829_patient_level.tsv", sep="\t", index=False)
    rows = []
    for gene in GENES:
        rows += binary_tests(df.dropna(subset=["cr_vs_pd"]), gene, "cr_vs_pd", "GSE283829_NSCLC", "CR_vs_PD")
        rows += binary_tests(df, gene, "dcr", "GSE283829_NSCLC", "DCR_CRSD_vs_PD")
    boxplot_by_group(df, "CLDN4", "recist", "GSE283829 CLDN4 log2(CPM+1) vs GEO disease stage", FIGS / "GSE283829_CLDN4_RECIST.png")
    boxplot_by_group(df, "TACSTD2", "recist", "GSE283829 TACSTD2 log2(CPM+1) vs GEO disease stage", FIGS / "GSE283829_TACSTD2_RECIST.png")
    return df, rows


# ---------------------------------------------------------------------------
# GSE328294 cell line
# ---------------------------------------------------------------------------

def analyze_gse328294() -> tuple[pd.DataFrame, list[dict]]:
    xlsx = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE328nnn/GSE328294/suppl/GSE328294_genes_fpkm_expression.xlsx",
        CACHE / "GSE328294_genes_fpkm_expression.xlsx",
        min_bytes=1000,
    )
    raw = pd.read_excel(xlsx)
    sub = raw[raw["gene_name"].isin(list(GENES) + ["CD274", "STC1"])].copy()
    fpkm_cols = [c for c in sub.columns if c.startswith("FPKM.")]
    long = sub.melt(id_vars=["gene_name"], value_vars=fpkm_cols, var_name="sample", value_name="fpkm")
    long["sample"] = long["sample"].str.replace("FPKM.", "", regex=False)
    long["arm"] = np.where(long["sample"].str.startswith("IBI351"), "IBI351", "DMSO")
    long["log2_fpkm1"] = np.log2(long["fpkm"].astype(float) + 1)
    wide = long.pivot_table(index="sample", columns="gene_name", values="log2_fpkm1")
    wide["arm"] = np.where(wide.index.str.startswith("IBI351"), "IBI351", "DMSO")
    wide.to_csv(PROC / "GSE328294_cellline_log2fpkm.tsv", sep="\t")
    long.to_csv(PROC / "GSE328294_gene_fpkm_long.tsv", sep="\t", index=False)
    rows = []
    for gene in GENES:
        a = wide.loc[wide["arm"] == "IBI351", gene].to_numpy()
        b = wide.loc[wide["arm"] == "DMSO", gene].to_numpy()
        u = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append(
            {
                "cohort": "GSE328294_H23_cellline",
                "gene": gene,
                "endpoint": "IBI351_vs_DMSO",
                "cutoff": "continuous",
                "n": 6,
                "n_high_or_R": 3,
                "n_low_or_NR": 3,
                "effect": "rank_biserial_IBI351_minus_DMSO",
                "effect_value": rank_biserial(a, b),
                "ci_low": np.nan,
                "ci_high": np.nan,
                "p": float(u.pvalue),
                "note": "NOT an ICI response OR. H23 ± IBI351 48h. No patient labels in GEO.",
            }
        )
    boxplot_by_group(wide.reset_index(), "CLDN4", "arm", "GSE328294 H23 CLDN4 log2(FPKM+1)", FIGS / "GSE328294_CLDN4_IBI351.png")
    boxplot_by_group(wide.reset_index(), "TACSTD2", "arm", "GSE328294 H23 TACSTD2 log2(FPKM+1)", FIGS / "GSE328294_TACSTD2_IBI351.png")
    return wide.reset_index(), rows


# ---------------------------------------------------------------------------
# GSE244945 SCLC cell lines
# ---------------------------------------------------------------------------

def analyze_gse244945() -> tuple[pd.DataFrame, list[dict]]:
    csv_p = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244945/suppl/GSE244945_Human_SCLC_RSEM_RPKM_log2.csv.gz",
        CACHE / "GSE244945_Human_SCLC_RSEM_RPKM_log2.csv.gz",
    )
    mat = pd.read_csv(csv_p, index_col=0)
    keep = [g for g in list(GENES) + ["NOTCH1", "CD274", "TMEM173"] if g in mat.index]
    expr = mat.loc[keep].T
    expr["sample"] = expr.index
    # H82 parental vs H82 NOTCH1 +dox (exclude STING KO)
    def arm(name: str) -> str | None:
        n = name
        if "STING_KO" in n:
            return None
        if re.search(r"H82_parental", n):
            return "H82_parental"
        if re.search(r"H82_NOTCH1_dox", n) and "no_dox" not in n:
            return "H82_NOTCH1_dox"
        if re.search(r"H69_parental", n):
            return "H69_parental"
        if re.search(r"H69_NOTCH1_dox", n) and "no_dox" not in n:
            return "H69_NOTCH1_dox"
        return None

    expr["arm"] = expr["sample"].map(arm)
    use = expr.dropna(subset=["arm"]).copy()
    use.to_csv(PROC / "GSE244945_SCLC_cellline_selected.tsv", sep="\t", index=False)
    rows = []
    for line, a0, a1 in (("H82", "H82_parental", "H82_NOTCH1_dox"), ("H69", "H69_parental", "H69_NOTCH1_dox")):
        for gene in GENES:
            if gene not in use.columns:
                continue
            x = use.loc[use["arm"] == a1, gene].to_numpy()
            y = use.loc[use["arm"] == a0, gene].to_numpy()
            if len(x) < 2 or len(y) < 2:
                continue
            u = stats.mannwhitneyu(x, y, alternative="two-sided")
            rows.append(
                {
                    "cohort": f"GSE244945_SCLC_{line}_cellline",
                    "gene": gene,
                    "endpoint": "NOTCH1_dox_vs_parental",
                    "cutoff": "continuous",
                    "n": int(len(x) + len(y)),
                    "n_high_or_R": int(len(x)),
                    "n_low_or_NR": int(len(y)),
                    "effect": "rank_biserial_NOTCH1_minus_parental",
                    "effect_value": rank_biserial(x, y),
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "p": float(u.pvalue),
                    "note": "SCLC cell line, not patient ICB. Do not mix into NSCLC OR.",
                }
            )
    if "CLDN4" in use.columns:
        boxplot_by_group(use, "CLDN4", "arm", "GSE244945 SCLC cell-line CLDN4 log2(RPKM+1)", FIGS / "GSE244945_CLDN4_NOTCH1.png")
        boxplot_by_group(use, "TACSTD2", "arm", "GSE244945 SCLC cell-line TACSTD2 log2(RPKM+1)", FIGS / "GSE244945_TACSTD2_NOTCH1.png")
    return use, rows


def inventory() -> pd.DataFrame:
    rows = [
        {
            "accession": "GSE328294",
            "year_public": 2026,
            "histology": "NSCLC cell line (H23 KRAS G12C)",
            "assay": "bulk RNA-seq FPKM n=6",
            "response_labels": "none (DMSO vs IBI351)",
            "disposition": "must-try scored as treatment contrast only; no ICI OR",
            "in_original_B5_11": "no",
        },
        {
            "accession": "GSE244944",
            "year_public": 2025,
            "histology": "SCLC cell lines",
            "assay": "ChIP-seq",
            "response_labels": "none",
            "disposition": "must-try; no expression-vs-ICB test. IMpower133 RNA not in GEO",
            "in_original_B5_11": "no",
        },
        {
            "accession": "GSE244945",
            "year_public": 2025,
            "histology": "SCLC cell lines",
            "assay": "bulk RNA-seq log2 RPKM",
            "response_labels": "none (NOTCH1 activation)",
            "disposition": "must-try; SCLC table only; no patient ICB labels",
            "in_original_B5_11": "no",
        },
        {
            "accession": "GSE244946",
            "year_public": 2025,
            "histology": "SCLC cell line COR-L88",
            "assay": "scRNA-seq n=3",
            "response_labels": "none (DMSO / LSD1i / GSI)",
            "disposition": "must-try; not patient ICB; not mixed into NSCLC OR",
            "in_original_B5_11": "no",
        },
        {
            "accession": "GSE274975",
            "year_public": 2024,
            "histology": "NSCLC tumor FFPE (1 SCLC excluded from NSCLC OR)",
            "assay": "bulk RNA-seq STAR counts n=61",
            "response_labels": "Table S1 RECIST + PFS (PMID 39723204)",
            "disposition": "INCLUDED extra NSCLC rows",
            "in_original_B5_11": "no",
        },
        {
            "accession": "GSE283829",
            "year_public": 2025,
            "histology": "NSCLC tumor",
            "assay": "bulk RNA-seq counts n=27",
            "response_labels": "GEO disease stage CR/SD/PD (no PR)",
            "disposition": "INCLUDED extra NSCLC rows",
            "in_original_B5_11": "no",
        },
        {
            "accession": "GSE202417",
            "year_public": 2025,
            "histology": "NSCLC",
            "assay": "Clariom D array of PBMC CD8 T cells n=28 (14 pts)",
            "response_labels": "R/NR in titles",
            "disposition": "SKIPPED for tumor OR: blood CD8, not tumor bulk",
            "in_original_B5_11": "no",
        },
        {
            "accession": "GSE317309",
            "year_public": 2026,
            "histology": "NSCLC",
            "assay": "scRNA/TCR tumor+PBMC",
            "response_labels": "not in GEO sample fields",
            "disposition": "SKIPPED: not bulk tumor matrix",
            "in_original_B5_11": "no",
        },
        {
            "accession": "GSE291670",
            "year_public": 2025,
            "histology": "NSCLC",
            "assay": "scRNA-seq n=6 (3 MPR / 3 non-MPR)",
            "response_labels": "MPR in titles",
            "disposition": "SKIPPED: scRNA, not bulk",
            "in_original_B5_11": "no",
        },
        {
            "accession": "GSE207422",
            "year_public": 2023,
            "histology": "NSCLC",
            "assay": "bulk + scRNA",
            "response_labels": "MPR",
            "disposition": "SKIPPED: original B5 / already used; not recut",
            "in_original_B5_11": "yes",
        },
        {
            "accession": "GSE253564 / GSE248378",
            "year_public": "2023–2024",
            "histology": "NSCLC neoadjuvant durvalumab ± SBRT",
            "assay": "bulk FPKM",
            "response_labels": "MPR/PFS in papers, not GEO matrix",
            "disposition": "not recut here (sibling leftover already exists); not original 11",
            "in_original_B5_11": "no",
        },
        {
            "accession": "IMpower133 RNA",
            "year_public": "closed",
            "histology": "SCLC",
            "assay": "trial RNA (Roche/Genentech)",
            "response_labels": "OS with atezo+chemo",
            "disposition": "SKIPPED: not in GEO; EGA/dbGaP-class / sponsor-controlled",
            "in_original_B5_11": "no",
        },
    ]
    return pd.DataFrame(rows)


def extra_row_view(stats_df: pd.DataFrame) -> pd.DataFrame:
    """Compact supplement rows: n, OR or rank-biserial, p."""
    keep_end = {
        "ORR_CRPR_vs_SDPD",
        "CR_vs_PD",
        "DCB_PFS_ge_6mo",
        "IBI351_vs_DMSO",
        "NOTCH1_dox_vs_parental",
        "PFS_cox_per_SD",
    }
    sub = stats_df[stats_df["endpoint"].isin(keep_end)].copy()
    # primary locked cutoff for OR rows is median; keep continuous rank-biserial too
    sub = sub[sub["cutoff"].isin(["median", "continuous", "tertile_T3_vs_T1"])]
    return sub


def main() -> None:
    inv = inventory()
    inv.to_csv(TABLES / "series_inventory.tsv", sep="\t", index=False)

    all_rows: list[dict] = []
    g274, r274 = analyze_gse274975()
    g283, r283 = analyze_gse283829()
    g328, r328 = analyze_gse328294()
    g244, r244 = analyze_gse244945()
    all_rows.extend(r274 + r283 + r328 + r244)
    stats_df = pd.DataFrame(all_rows)
    stats_df.to_csv(TABLES / "all_tests.tsv", sep="\t", index=False)

    nsclc = stats_df[stats_df["cohort"].str.contains("GSE274975_NSCLC|GSE283829_NSCLC")].copy()
    sclc = stats_df[stats_df["cohort"].str.contains("SCLC|GSE24494")].copy()
    extra_n = extra_row_view(nsclc)
    extra_s = extra_row_view(sclc)
    extra_n.to_csv(TABLES / "extra_rows_NSCLC.tsv", sep="\t", index=False)
    extra_s.to_csv(TABLES / "extra_rows_SCLC.tsv", sep="\t", index=False)

    # primary locked extra rows for the supplement (median OR + continuous rank-biserial)
    primary = []
    for cohort, endpoint in (
        ("GSE274975_NSCLC", "ORR_CRPR_vs_SDPD"),
        ("GSE274975_NSCLC", "DCB_PFS_ge_6mo"),
        ("GSE283829_NSCLC", "CR_vs_PD"),
    ):
        for gene in GENES:
            for cutoff in ("median", "continuous"):
                hit = stats_df[
                    (stats_df.cohort == cohort)
                    & (stats_df.gene == gene)
                    & (stats_df.endpoint == endpoint)
                    & (stats_df.cutoff == cutoff)
                ]
                primary.append(hit)
    prim = pd.concat(primary, ignore_index=True)
    prim.to_csv(TABLES / "extra_rows_NSCLC_primary.tsv", sep="\t", index=False)

    manifest = []
    for p in sorted(CACHE.glob("*")):
        if p.is_file() and p.stat().st_size < 80_000_000:
            manifest.append({"file": p.name, "bytes": p.stat().st_size, "sha256": sha256(p)})
    pd.DataFrame(manifest).to_csv(TABLES / "input_manifest.tsv", sep="\t", index=False)

    summary = {
        "claim_B5_taken_as_given": {
            "statement": "11-cohort CLDN4-high ICI OR=0.42",
            "recut": False,
            "pooled_with_new_rows": False,
        },
        "n_GSE274975_joined": int(len(g274)),
        "n_GSE274975_NSCLC": int((~g274["is_sclc"]).sum()),
        "n_GSE274975_SCLC_excluded": int(g274["is_sclc"].sum()),
        "GSE274975_RECIST": g274.loc[~g274["is_sclc"], "recist"].value_counts().to_dict(),
        "n_GSE283829": int(len(g283)),
        "GSE283829_RECIST": g283["recist"].value_counts().to_dict(),
        "GSE328294_note": "H23 cell line ± IBI351; no patient ICI labels",
        "GSE24494x_note": "SCLC cell line / ChIP / scRNA; IMpower133 RNA not public on GEO",
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(prim.to_string(index=False))


if __name__ == "__main__":
    main()
