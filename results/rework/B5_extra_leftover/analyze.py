#!/usr/bin/env python3
"""Extra 2023–2026 open lung ICI rows for a future expanded B5 meta.

The user's 11-cohort CLDN4-high ICI meta (OR=0.42) is taken as given.
This script does not re-cut those cohorts and does not pool new rows into that OR.

Named leftover series (score if DCB/ORR + tumor CLDN4 exist):
  GSE328294  H23 ± IBI351 (cell line; no patient ICI OR)
  GSE308745  PBMC scRNA/TCR before neoadjuvant PD-1 (no tumor-bulk OR)
  GSE302284  TROP2 CAR-T / EGFR DTP scRNA (no ICI/ADC response labels)
  GSE244944  SCLC ChIP-seq (separate SCLC table)
  GSE244945  SCLC cell-line RNA (separate SCLC table; IMpower133 not on GEO)
  GSE244946  SCLC COR-L88 scRNA (separate SCLC table)

Downloadable tumor-bulk leftovers with public DCB/ORR (not in the original 11):
  GSE274975  Poddubskaya 2024 NSCLC ICI, STAR counts + Table S1 RECIST/PFS
  GSE283829  Lindberg 2025 NSCLC ICI, counts + GEO CR/SD/PD

Also inventoried: GSE309652 (NanoString R/NR, CLDN4/TACSTD2 absent),
GSE309446 (TCR-seq only).

Public GEO / EuropePMC only. No EGA / dbGaP / FASTQ.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
import shutil
import urllib.request
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "rework" / "B5_extra_leftover"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
PROC = OUT / "processed"
CACHE = Path("/tmp/b5_extra_raw")
CACHE.mkdir(parents=True, exist_ok=True)
for d in (TABLES, FIGS, PROC):
    d.mkdir(parents=True, exist_ok=True)

GENES = {"CLDN4": "ENSG00000189143", "TACSTD2": "ENSG00000184292"}
UA = {"User-Agent": "B5-extra-leftover/1.0 (public GEO/EuropePMC)"}


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
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=180) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f)
    if tmp.stat().st_size < min_bytes:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"download too small: {url} -> {tmp.stat().st_size if tmp.exists() else 0}")
    tmp.replace(dest)
    return dest


def parse_soft_samples(soft_gz: Path) -> list[dict]:
    samples: list[dict] = []
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
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if len(x) == 0 or len(y) == 0:
        return float("nan")
    u = stats.mannwhitneyu(x, y, alternative="two-sided").statistic
    return float(2.0 * u / (len(x) * len(y)) - 1.0)


def haldane_or(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    aa, bb, cc, dd = a, b, c, d
    if min(a, b, c, d) == 0:
        aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    or_ = (aa * dd) / (bb * cc)
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    lo = math.exp(math.log(or_) - 1.96 * se)
    hi = math.exp(math.log(or_) + 1.96 * se)
    return float(or_), float(lo), float(hi)


def fisher_p(a: int, b: int, c: int, d: int) -> float:
    return float(stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")[1])


def assign_cutoff(values: pd.Series, how: str) -> pd.Series:
    v = values.astype(float)
    if how == "median":
        med = v.median()
        return pd.Series(np.where(v > med, "high", "low"), index=v.index)
    if how == "tertile":
        q1, q2 = v.quantile([1 / 3, 2 / 3])
        out = pd.Series(np.where(v >= q2, "high", np.where(v <= q1, "low", "mid")), index=v.index)
        return out
    raise ValueError(how)


def binary_tests(df: pd.DataFrame, gene: str, ycol: str, cohort: str, endpoint: str, note: str = "") -> list[dict]:
    rows = []
    sub = df.dropna(subset=[gene, ycol]).copy()
    sub[ycol] = pd.to_numeric(sub[ycol], errors="coerce")
    sub = sub.dropna(subset=[ycol])
    if sub.empty:
        return rows
    y = sub[ycol].astype(int)
    # continuous
    xr = sub.loc[y == 1, gene].to_numpy(float)
    xn = sub.loc[y == 0, gene].to_numpy(float)
    if len(xr) and len(xn):
        rb = cliffs_delta(xr, xn)
        p_mw = float(stats.mannwhitneyu(xr, xn, alternative="two-sided").pvalue)
        rows.append(
            {
                "cohort": cohort,
                "gene": gene,
                "endpoint": endpoint,
                "cutoff": "continuous",
                "n": int(len(sub)),
                "n_R": int(y.sum()),
                "n_NR": int((1 - y).sum()),
                "n_high_R": np.nan,
                "n_high_NR": np.nan,
                "n_low_R": np.nan,
                "n_low_NR": np.nan,
                "effect": "rank_biserial_R_minus_NR",
                "effect_value": rb,
                "ci_low": np.nan,
                "ci_high": np.nan,
                "p": p_mw,
                "note": note,
            }
        )
    for how in ("median", "tertile"):
        lab = assign_cutoff(sub[gene], how)
        use = sub.loc[lab.isin(["high", "low"])]
        hi = lab.loc[use.index].eq("high")
        yy = y.loc[use.index]
        a = int(((hi) & (yy == 1)).sum())
        b = int(((hi) & (yy == 0)).sum())
        c = int((~hi & (yy == 1)).sum())
        d = int((~hi & (yy == 0)).sum())
        if a + b + c + d < 4 or min(a + c, b + d) == 0:
            continue
        or_, lo, hi_ci = haldane_or(a, b, c, d)
        rows.append(
            {
                "cohort": cohort,
                "gene": gene,
                "endpoint": endpoint,
                "cutoff": "median" if how == "median" else "tertile_T3_vs_T1",
                "n": int(len(use)),
                "n_R": int(yy.sum()),
                "n_NR": int((1 - yy).sum()),
                "n_high_R": a,
                "n_high_NR": b,
                "n_low_R": c,
                "n_low_NR": d,
                "effect": "OR_high_vs_low",
                "effect_value": or_,
                "ci_low": lo,
                "ci_high": hi_ci,
                "p": fisher_p(a, b, c, d),
                "note": note,
            }
        )
    return rows


def boxplot_by_group(df: pd.DataFrame, gene: str, group: str, title: str, dest: Path) -> None:
    sub = df.dropna(subset=[gene, group]).copy()
    if sub.empty:
        return
    order = [g for g in ["CR", "PR", "SD", "PD", "DMSO", "IBI351"] if g in set(sub[group].astype(str))]
    if not order:
        order = sorted(sub[group].astype(str).unique())
    data = [sub.loc[sub[group].astype(str).eq(g), gene].to_numpy(float) for g in order]
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.boxplot(data, tick_labels=order, patch_artist=True)
    for i, arr in enumerate(data, start=1):
        ax.scatter(np.random.default_rng(1).normal(i, 0.06, size=len(arr)), arr, s=12, color="k", alpha=0.55)
    ax.set_ylabel(f"{gene}")
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def cox_hr_per_sd(time: np.ndarray, event: np.ndarray, x: np.ndarray) -> tuple[float, float]:
    """Univariable Cox via statsmodels if available; else nan."""
    try:
        import statsmodels.duration.hazard_regression as hr
    except Exception:
        return float("nan"), float("nan")
    t = np.asarray(time, float)
    e = np.asarray(event, float)
    z = np.asarray(x, float)
    ok = np.isfinite(t) & np.isfinite(e) & np.isfinite(z) & (t > 0)
    t, e, z = t[ok], e[ok], z[ok]
    if e.sum() < 2 or len(t) < 8:
        return float("nan"), float("nan")
    z = (z - z.mean()) / (z.std(ddof=0) or 1.0)
    try:
        model = hr.PHReg(t, z[:, None], status=e)
        res = model.fit(disp=0)
        beta = float(np.asarray(res.params).ravel()[0])
        p = float(np.asarray(res.pvalues).ravel()[0])
        return math.exp(beta), p
    except Exception:
        return float("nan"), float("nan")


# ---------------------------------------------------------------------------
# GSE274975 — extra NSCLC tumor-bulk ICI row
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
            }
            for s in samples
        ]
    )
    counts = pd.read_csv(counts_p, sep="\t", index_col=0)

    def col_to_sid(col: str) -> str | None:
        m = re.search(r"Lu[Cc][-_]?(\d+)", str(col))
        if not m:
            return None
        return f"OB_pat_LuC_{int(m.group(1))}"

    colmap = {c: col_to_sid(c) for c in counts.columns}
    missing = [c for c, v in colmap.items() if v is None]
    if missing:
        raise RuntimeError(f"unmapped count columns: {missing}")
    counts = counts.rename(columns=colmap)
    if counts.columns.duplicated().any():
        raise RuntimeError("duplicate sample map")
    expr = log2_cpm(counts).loc[list(GENES.values())].T
    expr.columns = list(GENES)
    expr["Sample_ID"] = expr.index
    df = clin.merge(expr, on="Sample_ID", how="inner")
    if not meta.empty:
        df = df.merge(meta, left_on="Sample_ID", right_on="title", how="left")
    if len(df) != 61:
        raise RuntimeError(f"GSE274975 join n={len(df)} expected 61")
    recist = df["RECIST Response"].astype(str).str.upper()
    recist = recist.where(df["RECIST Response"].notna() & ~recist.isin(["NAN", "NONE", "NA", ""]), other=pd.NA)
    df["recist"] = recist
    df["orr"] = np.where(recist.isin(["CR", "PR"]), 1, np.where(recist.isin(["SD", "PD"]), 0, np.nan))
    df["pfs"] = pd.to_numeric(df["PFS time, months"], errors="coerce")
    df["pfs_event"] = pd.to_numeric(df["Response status available"], errors="coerce").fillna(0).astype(int)
    # Locked DCB: PFS≥6 mo; early-censored before 6 mo = NA
    df["dcb6"] = np.where(
        df["pfs"] >= 6,
        1,
        np.where((df["pfs"] < 6) & (df["pfs_event"] == 1), 0, np.nan),
    )
    src = df.get("source", pd.Series("", index=df.index)).astype(str)
    histo = df["Histotype"].astype(str)
    df["is_sclc"] = src.str.contains("Small cell", case=False) | histo.str.contains("small cell", case=False)
    df["histo_bin"] = np.where(
        histo.str.contains("Squamous", case=False),
        "LUSC",
        np.where(histo.str.contains("Adeno", case=False), "LUAD", "other"),
    )
    nsclc = df.loc[~df["is_sclc"]].copy()
    df.to_csv(PROC / "GSE274975_patient_level.tsv", sep="\t", index=False)
    note = (
        "PMID 39723204 Table S1; DCB=PFS>=6mo, early-censored <6mo excluded; "
        "not pooled into the 11-cohort meta"
    )
    rows: list[dict] = []
    for gene in GENES:
        rows += binary_tests(nsclc, gene, "orr", "GSE274975_NSCLC", "ORR_CRPR_vs_SDPD", note)
        rows += binary_tests(nsclc, gene, "dcb6", "GSE274975_NSCLC", "DCB_PFS_ge_6mo", note)
        for histo_name, sub in (("LUAD", nsclc[nsclc["histo_bin"].eq("LUAD")]), ("LUSC", nsclc[nsclc["histo_bin"].eq("LUSC")])):
            rows += binary_tests(sub, gene, "orr", f"GSE274975_{histo_name}", "ORR_CRPR_vs_SDPD", note + f"; {histo_name} only")
            rows += binary_tests(sub, gene, "dcb6", f"GSE274975_{histo_name}", "DCB_PFS_ge_6mo", note + f"; {histo_name} only")
        surv = nsclc.dropna(subset=[gene, "pfs", "pfs_event"])
        hr, p_cox = cox_hr_per_sd(surv["pfs"].to_numpy(), surv["pfs_event"].to_numpy(), surv[gene].to_numpy())
        rows.append(
            {
                "cohort": "GSE274975_NSCLC",
                "gene": gene,
                "endpoint": "PFS_cox_per_SD",
                "cutoff": "continuous",
                "n": int(len(surv)),
                "n_R": int(surv["pfs_event"].sum()),
                "n_NR": int((1 - surv["pfs_event"]).sum()),
                "n_high_R": np.nan,
                "n_high_NR": np.nan,
                "n_low_R": np.nan,
                "n_low_NR": np.nan,
                "effect": "HR_per_SD",
                "effect_value": hr,
                "ci_low": np.nan,
                "ci_high": np.nan,
                "p": p_cox,
                "note": "univariable Cox; event=Table S1 Response status available",
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
# GSE283829 — leftover 2025 NSCLC tumor bulk (slide did not use)
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
                "gsm": s["gsm"],
                "title": s.get("title"),
                "recist": s["chars"].get("disease stage"),
                "tumor_type": s["chars"].get("tumor type"),
                "sex": s["chars"].get("sex"),
            }
            for s in samples
        ]
    )
    counts = pd.read_csv(counts_p, sep="\t", index_col=0)
    # columns S105691_047 <-> title 105691_047
    def col_to_title(col: str) -> str:
        return str(col)[1:] if str(col).startswith("S") else str(col)

    counts = counts.rename(columns=col_to_title)
    expr = log2_cpm(counts)
    missing = [g for g in GENES.values() if g not in expr.index]
    if missing:
        raise RuntimeError(f"GSE283829 missing genes {missing}")
    expr = expr.loc[list(GENES.values())].T
    expr.columns = list(GENES)
    expr["title"] = expr.index
    df = meta.merge(expr, on="title", how="inner")
    recist = df["recist"].astype(str).str.upper()
    df["recist"] = recist
    # GEO has CR/SD/PD only (no PR). ORR = CR vs SD+PD; also CR vs PD.
    df["orr"] = np.where(recist.eq("CR"), 1, np.where(recist.isin(["SD", "PD"]), 0, np.nan))
    df["cr_vs_pd"] = np.where(recist.eq("CR"), 1, np.where(recist.eq("PD"), 0, np.nan))
    df.to_csv(PROC / "GSE283829_patient_level.tsv", sep="\t", index=False)
    note = "GEO disease stage is best response CR/SD/PD (no PR, no PFS); not pooled into the 11-cohort meta"
    rows: list[dict] = []
    for gene in GENES:
        rows += binary_tests(df, gene, "orr", "GSE283829_NSCLC", "ORR_CR_vs_SDPD", note)
        rows += binary_tests(df, gene, "cr_vs_pd", "GSE283829_NSCLC", "CR_vs_PD", note)
    boxplot_by_group(df, "CLDN4", "recist", "GSE283829 NSCLC CLDN4 log2(CPM+1) vs GEO response", FIGS / "GSE283829_CLDN4_RECIST.png")
    boxplot_by_group(df, "TACSTD2", "recist", "GSE283829 NSCLC TACSTD2 log2(CPM+1) vs GEO response", FIGS / "GSE283829_TACSTD2_RECIST.png")
    return df, rows


# ---------------------------------------------------------------------------
# GSE328294 — cell line, not a patient OR
# ---------------------------------------------------------------------------

def analyze_gse328294() -> tuple[pd.DataFrame, list[dict]]:
    xlsx = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE328nnn/GSE328294/suppl/GSE328294_genes_fpkm_expression.xlsx",
        CACHE / "GSE328294_genes_fpkm_expression.xlsx",
        min_bytes=1000,
    )
    raw = pd.read_excel(xlsx)
    name_col = "gene_name" if "gene_name" in raw.columns else None
    if name_col is None:
        raise RuntimeError(f"GSE328294 columns: {list(raw.columns)[:12]}")
    fpkm_cols = [c for c in raw.columns if str(c).lower().startswith("fpkm.") or str(c).endswith("_fpkm") or "fpkm" in str(c).lower()]
    # file has count.* and likely fpkm.*
    if not fpkm_cols:
        fpkm_cols = [c for c in raw.columns if str(c).lower().startswith("fpkm")]
    if not fpkm_cols:
        # fallback: columns named like FPKM.CTL1
        fpkm_cols = [c for c in raw.columns if re.search(r"(CTL|IBI351)", str(c), re.I) and "count" not in str(c).lower()]
    sub = raw.loc[raw[name_col].isin(list(GENES)), [name_col] + fpkm_cols].copy()
    if sub.empty:
        # try gene_id
        sub = raw.loc[raw.get("gene_id", pd.Series()).isin(list(GENES.values())), ["gene_id"] + fpkm_cols].copy()
        sub = sub.rename(columns={"gene_id": name_col})
        sub[name_col] = sub[name_col].map({v: k for k, v in GENES.items()})
    long_rows = []
    for _, r in sub.iterrows():
        gene = r[name_col]
        for c in fpkm_cols:
            val = r[c]
            cs = str(c)
            treat = "IBI351" if "IBI351" in cs.upper() else ("DMSO" if "CTL" in cs.upper() or "DMSO" in cs.upper() else "other")
            long_rows.append({"gene": gene, "sample": c, "fpkm": float(val), "treatment": treat, "log2fpkm1": math.log2(float(val) + 1.0)})
    long = pd.DataFrame(long_rows)
    long.to_csv(PROC / "GSE328294_cellline_fpkm.tsv", sep="\t", index=False)
    rows = []
    for gene in GENES:
        g = long[long["gene"].eq(gene)]
        x = g.loc[g["treatment"].eq("IBI351"), "log2fpkm1"].to_numpy(float)
        y = g.loc[g["treatment"].eq("DMSO"), "log2fpkm1"].to_numpy(float)
        if len(x) and len(y):
            rb = cliffs_delta(x, y)
            p = float(stats.mannwhitneyu(x, y, alternative="two-sided").pvalue)
            rows.append(
                {
                    "cohort": "GSE328294_H23_cellline",
                    "gene": gene,
                    "endpoint": "IBI351_vs_DMSO",
                    "cutoff": "continuous",
                    "n": int(len(x) + len(y)),
                    "n_R": int(len(x)),
                    "n_NR": int(len(y)),
                    "n_high_R": np.nan,
                    "n_high_NR": np.nan,
                    "n_low_R": np.nan,
                    "n_low_NR": np.nan,
                    "effect": "rank_biserial_IBI351_minus_DMSO",
                    "effect_value": rb,
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "p": p,
                    "note": "H23 cell line; not a patient ICI DCB/ORR row",
                }
            )
    wide = long.pivot_table(index="sample", columns="gene", values="log2fpkm1", aggfunc="first").reset_index()
    wide["treatment"] = wide["sample"].map(lambda s: "IBI351" if "IBI351" in str(s).upper() else "DMSO")
    if "CLDN4" in wide.columns:
        boxplot_by_group(wide, "CLDN4", "treatment", "GSE328294 H23 CLDN4 log2(FPKM+1)", FIGS / "GSE328294_CLDN4_IBI351.png")
    if "TACSTD2" in wide.columns:
        boxplot_by_group(wide, "TACSTD2", "treatment", "GSE328294 H23 TACSTD2 log2(FPKM+1)", FIGS / "GSE328294_TACSTD2_IBI351.png")
    return long, rows


# ---------------------------------------------------------------------------
# GSE244945 — SCLC cell-line RNA (separate table)
# ---------------------------------------------------------------------------

def analyze_gse244945() -> tuple[pd.DataFrame, list[dict]]:
    csv_p = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244945/suppl/GSE244945_Human_SCLC_RSEM_RPKM_log2.csv.gz",
        CACHE / "GSE244945_Human_SCLC_RSEM_RPKM_log2.csv.gz",
        min_bytes=1000,
    )
    mat = pd.read_csv(csv_p, index_col=0)
    have = [g for g in GENES if g in mat.index]
    if not have:
        return pd.DataFrame(), []
    expr = mat.loc[have].T
    expr["sample"] = expr.index.astype(str)

    def parse_sample(s: str) -> dict:
        # e.g. 7_H82_NOTCH1_dox , 10_CORL88_DMSO
        body = re.sub(r"^\d+_", "", s)
        body = re.sub(r"\.fastq\.gz$", "", body, flags=re.I)
        line = "unknown"
        for key in ("H82", "H524", "H69", "CORL88", "COR-L88", "H82"):
            if key.replace("-", "") in body.replace("-", "").upper() or key in body.upper():
                line = key.replace("-", "")
                break
        if "REST" in body.upper() or "_KO" in body.upper() or "KO_" in body.upper():
            treat = "KO_related"
        elif "NOTCH1" in body.upper() and "dox" in body and "no_dox" not in body:
            treat = "NOTCH1_dox"
        elif "no_dox" in body:
            treat = "NOTCH1_no_dox"
        elif "parental" in body:
            treat = "parental"
        elif "TAS1440_GSI" in body.upper() or "TAS1440_GSI" in body:
            treat = "TAS1440_GSI"
        elif "TAS1440" in body.upper():
            treat = "TAS1440"
        elif "DMSO" in body.upper():
            treat = "DMSO"
        else:
            treat = body
        return {"line": line, "treatment": treat}

    parsed = expr["sample"].map(parse_sample).apply(pd.Series)
    df = pd.concat([expr.reset_index(drop=True), parsed.reset_index(drop=True)], axis=1)
    df.to_csv(PROC / "GSE244945_SCLC_cellline.tsv", sep="\t", index=False)
    rows = []
    for line in sorted(df["line"].dropna().unique()):
        sub = df[df["line"].eq(line)]
        if {"parental", "NOTCH1_dox"}.issubset(set(sub["treatment"])):
            for gene in have:
                x = sub.loc[sub["treatment"].eq("NOTCH1_dox"), gene].to_numpy(float)
                y = sub.loc[sub["treatment"].eq("parental"), gene].to_numpy(float)
                if len(x) >= 2 and len(y) >= 2:
                    rows.append(
                        {
                            "cohort": f"GSE244945_{line}",
                            "gene": gene,
                            "endpoint": "NOTCH1_dox_vs_parental",
                            "cutoff": "continuous",
                            "n": int(len(x) + len(y)),
                            "n_R": int(len(x)),
                            "n_NR": int(len(y)),
                            "n_high_R": np.nan,
                            "n_high_NR": np.nan,
                            "n_low_R": np.nan,
                            "n_low_NR": np.nan,
                            "effect": "rank_biserial_dox_minus_parental",
                            "effect_value": cliffs_delta(x, y),
                            "ci_low": np.nan,
                            "ci_high": np.nan,
                            "p": float(stats.mannwhitneyu(x, y, alternative="two-sided").pvalue),
                            "note": "SCLC cell line; no patient ICB labels; IMpower133 RNA not on GEO",
                        }
                    )
    return df, rows


def inventory_rows() -> list[dict]:
    return [
        {
            "accession": "GSE328294",
            "year_public": 2026,
            "what_it_is": "H23 KRAS-mutant NSCLC cell line, DMSO vs IBI351, n=6 FPKM",
            "tumor_bulk": "no (cell line)",
            "CLDN4": "yes",
            "DCB_ORR_labels": "no patient ICI/DCB/ORR",
            "extra_OR_row": "no",
            "table": "cell-line note only",
        },
        {
            "accession": "GSE308745",
            "year_public": 2026,
            "what_it_is": "PBMC scRNA/scTCR/ADT, CD3/CD4 T cells, before neoadjuvant PD-1 (PB cohort2); PMID 41904165",
            "tumor_bulk": "no (PBMC scRNA)",
            "CLDN4": "not scored (blood T cells, not tumor epithelium)",
            "DCB_ORR_labels": "none in GEO SOFT (titles RE/I/G are not a published codebook)",
            "extra_OR_row": "no",
            "table": "inventory",
        },
        {
            "accession": "GSE302284",
            "year_public": 2025,
            "what_it_is": "scRNA of residual EGFR-mutant tumor/LN + DFCI282/PC9 osi vs vehicle; TROP2 CAR-T paper",
            "tumor_bulk": "no (scRNA, n=6 libraries)",
            "CLDN4": "not scored",
            "DCB_ORR_labels": "no ICI or ADC response labels (treatment = osimertinib vs vehicle)",
            "extra_OR_row": "no",
            "table": "inventory",
        },
        {
            "accession": "GSE244944",
            "year_public": 2025,
            "what_it_is": "SCLC H3K27ac ChIP-seq, COR-L88 / NCI-H82, DMSO vs TAS1440 ± GSI",
            "tumor_bulk": "no (ChIP-seq cell line)",
            "CLDN4": "no expression matrix",
            "DCB_ORR_labels": "no patient ICB; IMpower133 RNA not on GEO",
            "extra_OR_row": "no",
            "table": "SCLC separate",
        },
        {
            "accession": "GSE244945",
            "year_public": 2025,
            "what_it_is": "SCLC cell-line bulk RNA-seq (H82/H524/COR-L88/H69), NOTCH1/LSD1/REST",
            "tumor_bulk": "no (cell line)",
            "CLDN4": "yes",
            "DCB_ORR_labels": "no patient ICB; IMpower133 RNA not on GEO",
            "extra_OR_row": "no",
            "table": "SCLC separate (mechanistic)",
        },
        {
            "accession": "GSE244946",
            "year_public": 2025,
            "what_it_is": "COR-L88 scRNA, DMSO / TAS1440 / TAS1440+GSI, n=3",
            "tumor_bulk": "no (scRNA cell line)",
            "CLDN4": "not scored as patient OR",
            "DCB_ORR_labels": "no patient ICB",
            "extra_OR_row": "no",
            "table": "SCLC separate",
        },
        {
            "accession": "GSE274975",
            "year_public": 2024,
            "what_it_is": "61 FFPE lung tumors, STAR counts; ICI RECIST+PFS in PMID 39723204 Table S1",
            "tumor_bulk": "yes",
            "CLDN4": "yes (ENSG00000189143)",
            "DCB_ORR_labels": "yes (RECIST + PFS months)",
            "extra_OR_row": "yes",
            "table": "NSCLC extra",
        },
        {
            "accession": "GSE283829",
            "year_public": 2025,
            "what_it_is": "27 NSCLC tumors, ICI; GEO disease stage = CR/SD/PD",
            "tumor_bulk": "yes",
            "CLDN4": "yes",
            "DCB_ORR_labels": "ORR-like CR vs SD+PD; no PFS/DCB",
            "extra_OR_row": "yes",
            "table": "NSCLC extra",
        },
        {
            "accession": "GSE309652",
            "year_public": 2025,
            "what_it_is": "NanoString nCounter, n=72 stage-IV NSCLC anti-PD-(L)1; author R/NR 24/48",
            "tumor_bulk": "yes (FFPE NanoString)",
            "CLDN4": "no (not on 768-gene panel; TACSTD2 also absent)",
            "DCB_ORR_labels": "author R/NR present; cannot score CLDN4",
            "extra_OR_row": "no",
            "table": "inventory",
        },
        {
            "accession": "GSE309446",
            "year_public": 2025,
            "what_it_is": "Phase II sitravatinib + tislelizumab + docetaxel; TCR-seq blood/FFPE",
            "tumor_bulk": "no RNA expression matrix (TCR clones only)",
            "CLDN4": "no",
            "DCB_ORR_labels": "no per-patient DCB/ORR in GEO",
            "extra_OR_row": "no",
            "table": "inventory",
        },
    ]


def fmt_or(row: dict) -> str:
    if row["effect"] != "OR_high_vs_low" or not np.isfinite(row.get("effect_value", np.nan)):
        return ""
    return f"{row['effect_value']:.2f} ({row['ci_low']:.2f}–{row['ci_high']:.2f})"


def main() -> None:
    all_rows: list[dict] = []
    g274, r274 = analyze_gse274975()
    all_rows += r274
    g283, r283 = analyze_gse283829()
    all_rows += r283
    g328, r328 = analyze_gse328294()
    all_rows += r328
    g244, r244 = analyze_gse244945()
    all_rows += r244

    stats_df = pd.DataFrame(all_rows)
    stats_df.to_csv(TABLES / "all_tests.tsv", sep="\t", index=False)

    inv = pd.DataFrame(inventory_rows())
    inv.to_csv(TABLES / "series_inventory.tsv", sep="\t", index=False)

    # Extra NSCLC rows: locked median OR + continuous complement
    nsclc = stats_df[stats_df["cohort"].str.contains("GSE274975_NSCLC|GSE283829_NSCLC")].copy()
    primary = nsclc[
        nsclc["cutoff"].isin(["median", "continuous"])
        & nsclc["endpoint"].isin(["ORR_CRPR_vs_SDPD", "DCB_PFS_ge_6mo", "ORR_CR_vs_SDPD", "CR_vs_PD", "PFS_cox_per_SD"])
        & nsclc["cohort"].isin(["GSE274975_NSCLC", "GSE283829_NSCLC"])
    ].copy()
    primary.to_csv(TABLES / "extra_rows_NSCLC.tsv", sep="\t", index=False)

    histo = stats_df[stats_df["cohort"].str.contains("GSE274975_LUAD|GSE274975_LUSC")].copy()
    histo.to_csv(TABLES / "extra_rows_NSCLC_histology.tsv", sep="\t", index=False)

    sclc = stats_df[stats_df["cohort"].str.contains("GSE244945|GSE328294")].copy()
    # keep 328294 out of SCLC table
    sclc_only = stats_df[stats_df["cohort"].str.contains("GSE244945")].copy()
    # add empty-label rows for 244944/6
    extra_sclc = [
        {
            "cohort": "GSE244944_ChIPseq",
            "gene": "—",
            "endpoint": "none",
            "cutoff": "—",
            "n": 18,
            "n_R": np.nan,
            "n_NR": np.nan,
            "n_high_R": np.nan,
            "n_high_NR": np.nan,
            "n_low_R": np.nan,
            "n_low_NR": np.nan,
            "effect": "—",
            "effect_value": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "p": np.nan,
            "note": "SCLC ChIP-seq cell line; no patient ICB; not mixed into NSCLC",
        },
        {
            "cohort": "GSE244946_scRNA",
            "gene": "—",
            "endpoint": "none",
            "cutoff": "—",
            "n": 3,
            "n_R": np.nan,
            "n_NR": np.nan,
            "n_high_R": np.nan,
            "n_high_NR": np.nan,
            "n_low_R": np.nan,
            "n_low_NR": np.nan,
            "effect": "—",
            "effect_value": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "p": np.nan,
            "note": "COR-L88 scRNA; no patient ICB; not mixed into NSCLC",
        },
    ]
    sclc_tab = pd.concat([sclc_only, pd.DataFrame(extra_sclc)], ignore_index=True)
    sclc_tab.to_csv(TABLES / "extra_rows_SCLC.tsv", sep="\t", index=False)

    cell = stats_df[stats_df["cohort"].str.contains("GSE328294")].copy()
    cell.to_csv(TABLES / "GSE328294_cellline_note.tsv", sep="\t", index=False)

    summary = {
        "scope": "extra leftover 2023-2026 open lung ICI rows; 11-cohort OR=0.42 taken as given; not pooled",
        "n_GSE274975_joined": int(len(g274)),
        "n_GSE274975_NSCLC": int((~g274["is_sclc"]).sum()),
        "n_GSE274975_SCLC_excluded": int(g274["is_sclc"].sum()),
        "GSE274975_RECIST_NSCLC": g274.loc[~g274["is_sclc"], "recist"].value_counts(dropna=False).to_dict(),
        "GSE274975_DCB_NSCLC": g274.loc[~g274["is_sclc"], "dcb6"].value_counts(dropna=False).to_dict(),
        "n_GSE283829": int(len(g283)),
        "GSE283829_RECIST": g283["recist"].value_counts().to_dict(),
        "n_GSE328294": int(g328["sample"].nunique()) if len(g328) else 0,
        "n_GSE244945_libraries": int(len(g244)) if len(g244) else 0,
        "inputs": {
            "GSE274975_counts_sha256": sha256(CACHE / "GSE274975_raw_counts.tsv.gz"),
            "TableS1_sha256": sha256(CACHE / "supp" / "Table1.xlsx"),
            "GSE283829_counts_sha256": sha256(CACHE / "GSE283829_raw_express_matrix_all_samples.txt.gz"),
        },
    }
    # cast numpy types
    def _cast(o):
        if isinstance(o, dict):
            return {str(k): _cast(v) for k, v in o.items()}
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        return o

    (TABLES / "summary.json").write_text(json.dumps(_cast(summary), indent=2) + "\n")
    print(json.dumps(_cast(summary), indent=2))
    print(stats_df.to_string(index=False))


if __name__ == "__main__":
    main()
