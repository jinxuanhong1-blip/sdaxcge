#!/usr/bin/env python3
"""B5 wave 2: lung-first then all-open ICI, CLDN4 and TACSTD2.

Locked before looking at pooled numbers
--------------------------------------
Genes: CLDN4 and TACSTD2. Both raw and ESTIMATEScore-residualized.
Cutoffs: median; tertile T3 vs T1; continuous logistic OR per 1 SD.
Endpoints kept separate (never mixed in one meta):
  DCB  = PFS >= 6 months (180 days). Early-censored before 6 mo = NA.
  ORR  = RECIST CR/PR vs SD+PD.
  author_R_vs_NR = GEO/PredictIO native label only (not DCB, not RECIST).
Histology split only from pathology fields (LUAD vs LUSC). No inferred histology
in the primary tables.
Meta: DerSimonian–Laird random effects on logOR / logHR.
Lung-first, then all independent open ICI. Same patients are never counted twice.
User claim CLDN4-high OR=0.42, k=11 is a reference line, not a target.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import math
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact, norm
from statsmodels.duration.hazard_regression import PHReg

HERE = Path(__file__).resolve().parent
GMT_PATH = HERE / "SI_geneset.gmt"

USER_OR = 0.42
USER_LO = 0.18
USER_HI = 0.95
USER_K = 11
DCB_MONTHS = 6.0
DCB_DAYS = 180.0

ENSG = {
    "TACSTD2": "ENSG00000184292",
    "CLDN4": "ENSG00000189143",
    "CD8A": "ENSG00000153563",
}
ALIAS = {
    "WISP1": "CCN4",
    "CCN4": "WISP1",
    "GPR124": "ADGRA2",
    "ADGRA2": "GPR124",
    "LPPR4": "PLPPR4",
    "PLPPR4": "LPPR4",
    "ODZ4": "TENM4",
    "TENM4": "ODZ4",
    "TXNDC3": "NME8",
    "NME8": "TXNDC3",
}

# Same patients as GSE135222; excluded from every meta.
DUPLICATE_OF = {"ICB_Jung": "GSE135222"}

# Locked independent sets (filled after load with whoever is actually usable).
LUNG_DCB_CANDIDATES = ["GSE135222", "GSE190265", "GSE190266"]
LUNG_ORR_CANDIDATES = ["GSE207422", "GSE283829"]
LUNG_AUTHOR_CANDIDATES = ["GSE126044", "GSE166449"]
# All-open adds Zenodo ICB studies that are not Jung and not lung-GEO duplicates.
NONLUNG_ICB = [
    "IMvigor210",
    "Gide",
    "Hugo",
    "Liu",
    "Riaz",
    "Puch",
    "Braun",
    "Shiuan",
    "Miao1",
    "Snyder",
    "Kim",
]
CTLA4_ONLY = ["Van_Allen", "Nathanson"]


def load_gmt(path: Path) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    for line in path.read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def parse_series_matrix(path: Path) -> pd.DataFrame:
    titles = geo = descriptions = None
    char_rows: list[list[str]] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            parts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
            key, vals = parts[0], parts[1:]
            if key == "!Sample_title":
                titles = vals
            elif key == "!Sample_geo_accession":
                geo = vals
            elif key == "!Sample_description":
                descriptions = vals
            elif key == "!Sample_characteristics_ch1":
                char_rows.append(vals)
    rows = []
    for i, title in enumerate(titles or []):
        d = {"title": title}
        if geo and i < len(geo):
            d["gsm"] = geo[i]
        if descriptions and i < len(descriptions):
            d["description"] = descriptions[i]
        for row in char_rows:
            if i < len(row) and row[i] and ":" in row[i]:
                k, v = row[i].split(":", 1)
                d[k.strip().lower()] = v.strip()
        rows.append(d)
    return pd.DataFrame(rows)


def norm_histology(value) -> str:
    if not isinstance(value, str):
        return "NA"
    v = value.strip().lower()
    if v in {
        "squamous",
        "sqcc",
        "lscc",
        "lusc",
        "squamous cell carcinoma",
        "lung squamous cell carcinoma",
        "squamous cell lung carcinoma",
    }:
        return "LUSC"
    if v in {
        "adeno",
        "adenocarcinoma",
        "adc",
        "ac",
        "luad",
        "lung adenocarcinoma",
        "adenocarcinoma of lung",
    }:
        return "LUAD"
    if v in {"non-squamous", "nonsquamous", "non squamous"}:
        return "nonsquamous"
    return "NA"


def collapse_index(expr: pd.DataFrame) -> pd.DataFrame:
    idx = pd.Index([str(i).split(".")[0] for i in expr.index])
    expr = expr.copy()
    expr.index = idx
    expr = expr[~expr.index.duplicated(keep="first")]
    expr = expr.apply(pd.to_numeric, errors="coerce")
    return expr


def map_gene(expr: pd.DataFrame, symbol: str) -> str | None:
    idx_u = {str(g).upper(): g for g in expr.index}
    if symbol.upper() in idx_u:
        return idx_u[symbol.upper()]
    ensg = ENSG.get(symbol)
    if ensg and ensg in idx_u:
        return idx_u[ensg]
    # Ensembl already stripped of version in collapse_index
    if ensg and ensg.upper() in idx_u:
        return idx_u[ensg.upper()]
    alt = ALIAS.get(symbol)
    if alt and alt.upper() in idx_u:
        return idx_u[alt.upper()]
    return None


def resolve_set(expr: pd.DataFrame, genes: list[str]) -> dict[str, str]:
    resolved = {}
    idx_u = {str(g).upper(): g for g in expr.index}
    for g in genes:
        if g.upper() in idx_u:
            resolved[g] = idx_u[g.upper()]
            continue
        alt = ALIAS.get(g)
        if alt and alt.upper() in idx_u:
            resolved[g] = idx_u[alt.upper()]
    return resolved


def estimate_scores(expr: pd.DataFrame, gene_sets: dict[str, list[str]]):
    genes = expr.index.to_numpy()
    n_genes = expr.shape[0]
    m = expr.rank(axis=0, method="average").to_numpy() * (10000.0 / n_genes)
    rows = {}
    overlaps = {}
    for name, out_col in (("StromalSignature", "StromalScore"), ("ImmuneSignature", "ImmuneScore")):
        mapped = [resolve_set(expr, gene_sets[name])[g] for g in gene_sets[name] if g in resolve_set(expr, gene_sets[name])]
        # resolve once
        resolved = resolve_set(expr, gene_sets[name])
        mapped = list(resolved.values())
        in_set = np.isin(genes, mapped)
        overlaps[name] = int(in_set.sum())
        es = np.empty(expr.shape[1])
        for j in range(expr.shape[1]):
            order = np.argsort(-m[:, j], kind="stable")
            tag = in_set[order].astype(float)
            w = np.abs(m[order, j]) ** 0.25
            nh = tag.sum()
            nm = n_genes - nh
            if nh == 0 or (tag * w).sum() == 0 or nm == 0:
                es[j] = np.nan
                continue
            pn = np.cumsum(tag * w) / (tag * w).sum()
            p0 = np.cumsum((1.0 - tag) / nm)
            es[j] = float(np.sum(pn - p0))
        rows[out_col] = es
    df = pd.DataFrame(rows, index=expr.columns)
    df["ESTIMATEScore"] = df["StromalScore"] + df["ImmuneScore"]
    return df, overlaps


def residualize(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    y = np.asarray(y, float)
    x = np.asarray(x, float)
    out = np.full_like(y, np.nan, dtype=float)
    m = ~(np.isnan(y) | np.isnan(x))
    if m.sum() < 3:
        return out
    X = np.column_stack([np.ones(m.sum()), x[m]])
    beta, *_ = np.linalg.lstsq(X, y[m], rcond=None)
    out[m] = y[m] - X @ beta
    return out


def recode_recist(s: pd.Series) -> pd.Series:
    u = s.astype(str).str.upper().str.strip()
    out = pd.Series(pd.NA, index=s.index, dtype=object)
    out[u.isin(["CR", "PR", "PRCR", "CRPR"])] = "R"
    out[u.isin(["PD"])] = "NR"
    out[u.isin(["SD"])] = "SD"
    return out


def dcb_from_pfs(time, event, unit: str) -> pd.Series:
    t = pd.to_numeric(time, errors="coerce")
    e = pd.to_numeric(event, errors="coerce")
    unit = (unit or "").lower()
    if "day" in unit:
        thr = DCB_DAYS
    elif "month" in unit:
        thr = DCB_MONTHS
    else:
        thr = DCB_DAYS if float(np.nanmax(t)) > 36 else DCB_MONTHS
    out = pd.Series(pd.NA, index=t.index if hasattr(t, "index") else range(len(t)), dtype=object)
    t = pd.Series(t).reset_index(drop=True)
    e = pd.Series(e).reset_index(drop=True)
    out = pd.Series(pd.NA, index=t.index, dtype=object)
    out[t >= thr] = "R"
    out[(t < thr) & (e == 1)] = "NR"
    # event 0 and t < thr remains NA (early censor)
    return out


def assign_cutoff(x: pd.Series, cutoff: str) -> pd.Series:
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
    raise ValueError(cutoff)


def or_from_2x2(a: int, b: int, c: int, d: int) -> dict:
    cells = [a, b, c, d]
    used_haldane = any(x == 0 for x in cells)
    aa, bb, cc, dd = ([x + 0.5 for x in cells] if used_haldane else cells)
    logor = math.log((aa * dd) / (bb * cc))
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    _, fisher_p = fisher_exact(np.array([[a, b], [c, d]], dtype=int), alternative="two-sided")
    return {
        "n_high_R": a,
        "n_high_NR": b,
        "n_low_R": c,
        "n_low_NR": d,
        "n_used": int(a + b + c + d),
        "or": math.exp(logor),
        "logor": logor,
        "se_logor": se,
        "or_lo": math.exp(logor - 1.96 * se),
        "or_hi": math.exp(logor + 1.96 * se),
        "p": float(fisher_p),
        "haldane": used_haldane,
        "estimable": True,
        "effect": "OR",
    }


def logistic_or_per_sd(y: pd.Series, x: pd.Series) -> dict:
    mask = y.notna() & x.notna()
    yy = (y[mask] == "R").astype(int)
    xx = pd.to_numeric(x[mask], errors="coerce")
    if yy.nunique() < 2 or xx.notna().sum() < 8:
        return {"estimable": False, "reason": "too_few_or_one_class", "n_used": int(mask.sum()), "effect": "OR"}
    z = (xx - xx.mean()) / (xx.std(ddof=1) if xx.std(ddof=1) else 1.0)
    X = sm.add_constant(z.to_numpy())
    try:
        fit = sm.Logit(yy.to_numpy(), X).fit(disp=False, maxiter=100)
    except Exception as e:
        return {"estimable": False, "reason": str(e), "n_used": int(mask.sum()), "effect": "OR"}
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
        "p": float(fit.pvalues[1]),
        "estimable": True,
        "effect": "OR",
    }


def cox_hr(time, event, x, binary: bool = False) -> dict:
    t = pd.to_numeric(time, errors="coerce")
    e = pd.to_numeric(event, errors="coerce")
    xx = pd.to_numeric(x, errors="coerce")
    m = t.notna() & e.notna() & xx.notna() & (t > 0)
    t, e, xx = t[m], e[m], xx[m]
    if len(t) < 8 or int(e.sum()) < 3 or xx.nunique() < 2:
        return {"estimable": False, "reason": "too_few_events", "n_used": int(m.sum()), "n_events": int(e.sum()) if len(e) else 0, "effect": "HR"}
    if not binary:
        z = (xx - xx.mean()) / (xx.std(ddof=1) if xx.std(ddof=1) else 1.0)
        exog = z.to_numpy()
    else:
        exog = xx.to_numpy()
    try:
        res = PHReg(t.to_numpy(), exog, status=e.to_numpy()).fit(disp=0)
        loghr = float(np.asarray(res.params).reshape(-1)[0])
        se = float(np.asarray(res.bse).reshape(-1)[0])
        p = float(np.asarray(res.pvalues).reshape(-1)[0])
    except Exception as exc:
        return {"estimable": False, "reason": str(exc), "n_used": int(len(t)), "n_events": int(e.sum()), "effect": "HR"}
    return {
        "n_used": int(len(t)),
        "n_events": int(e.sum()),
        "or": math.exp(loghr),  # HR stored in or fields for meta reuse
        "logor": loghr,
        "se_logor": se,
        "or_lo": math.exp(loghr - 1.96 * se),
        "or_hi": math.exp(loghr + 1.96 * se),
        "p": p,
        "estimable": True,
        "effect": "HR",
    }


def dl_meta(rows: pd.DataFrame) -> dict:
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


# ---------------------------------------------------------------------------
# Cohort loaders
# ---------------------------------------------------------------------------

def load_gse126044(raw: Path):
    meta = parse_series_matrix(raw / "GSE126044_series_matrix.txt.gz")
    counts = pd.read_csv(raw / "GSE126044_counts.txt.gz", sep="\t", index_col=0)
    counts = collapse_index(counts)
    cpm = counts.divide(counts.sum(axis=0), axis=1) * 1e6
    expr = np.log2(cpm.clip(lower=0) + 1)
    recs = []
    for _, s in meta.iterrows():
        sid = str(s["title"]).replace("RNA-seq_", "")
        resp = str(s.get("patient response", "")).lower()
        y = "R" if resp == "responder" else ("NR" if resp in {"non-responder", "nonresponder"} else pd.NA)
        recs.append({"sample": sid, "author": y, "histology": "NA", "histology_source": "none_in_GEO"})
    clin = pd.DataFrame(recs).drop_duplicates("sample").set_index("sample")
    keep = [c for c in expr.columns if c in clin.index]
    return expr[keep], clin.loc[keep], "NSCLC", "lung", "author_R_vs_NR", "log2CPM"


def load_gse135222(raw: Path):
    meta = parse_series_matrix(raw / "GSE135222_series_matrix.txt.gz")
    tpm = pd.read_csv(raw / "GSE135222_exp.tsv.gz", sep="\t", index_col=0)
    tpm = collapse_index(tpm)
    expr = np.log2(tpm.clip(lower=0) + 1)
    recs = []
    for _, s in meta.iterrows():
        sid = str(s["title"]).replace(" ", "")
        recs.append(
            {
                "sample": sid,
                "pfs_time": pd.to_numeric(s.get("pfs.time"), errors="coerce"),
                "pfs_event": pd.to_numeric(s.get("progression-free survival (pfs)"), errors="coerce"),
                "pfs_unit": "day",
                "histology": "NA",
                "histology_source": "none_in_GEO",
            }
        )
    clin = pd.DataFrame(recs).drop_duplicates("sample").set_index("sample")
    clin["dcb"] = dcb_from_pfs(clin["pfs_time"], clin["pfs_event"], "day").values
    keep = [c for c in expr.columns if c in clin.index]
    return expr[keep], clin.loc[keep], "NSCLC", "lung", "DCB", "log2(TPM+1)"


def load_gse166449(raw: Path):
    meta = parse_series_matrix(raw / "GSE166449_series_matrix.txt.gz")
    tpm = pd.read_csv(raw / "GSE166449_TPM.txt.gz", sep="\t", index_col=0)
    tpm = collapse_index(tpm)
    expr = np.log2(tpm.clip(lower=0) + 1)
    recs = []
    for _, s in meta.iterrows():
        sid = s.get("description")
        title = str(s.get("title", "")).lower()
        if "nonresponder" in title or "non-responder" in title:
            y = "NR"
        elif "responder" in title:
            y = "R"
        else:
            y = pd.NA
        recs.append({"sample": sid, "author": y, "histology": "NA", "histology_source": "none_in_GEO"})
    clin = pd.DataFrame(recs).drop_duplicates("sample").set_index("sample")
    keep = [c for c in expr.columns if c in clin.index]
    return expr[keep], clin.loc[keep], "NSCLC", "lung", "author_R_vs_NR", "log2(TPM+1)"


def load_gse190265(raw: Path):
    path = raw / "GSE190265_TPM_France3.csv.gz"
    with gzip.open(path, "rt") as fh:
        genes = fh.readline().strip().split(";")
    tpm = pd.read_csv(path, sep=";", header=None, skiprows=1, index_col=0)
    tpm.columns = genes
    tpm = collapse_index(tpm.T)
    expr = np.log2(tpm.clip(lower=0).astype(float) + 1)
    info = pd.read_csv(raw / "GSE190265_samples_info_France3.csv.gz", sep=";").set_index("sample")
    meta = parse_series_matrix(raw / "GSE190265_series_matrix.txt.gz").set_index("title")
    common = [c for c in expr.columns if c in info.index]
    clin = pd.DataFrame(index=common)
    clin["pfs_time"] = pd.to_numeric(info.loc[common, "time_PFS"], errors="coerce").values
    clin["pfs_event"] = pd.to_numeric(info.loc[common, "evtPFS"], errors="coerce").values
    clin["pfs_unit"] = "month"
    clin["dcb"] = dcb_from_pfs(clin["pfs_time"], clin["pfs_event"], "month").values
    clin["histology"] = [norm_histology(meta["disease state"].get(s) if s in meta.index else None) for s in common]
    clin["histology_source"] = "pathology_GEO"
    return expr[common], clin, "NSCLC", "lung", "DCB", "log2(TPM+1)"


def load_gse190266(raw: Path):
    path = raw / "GSE190266_TPM_France4.csv.gz"
    with gzip.open(path, "rt", encoding="utf8", errors="replace") as fh:
        text = fh.read()
    frame = pd.read_csv(io.StringIO(text), sep=";", decimal=",", index_col=0, low_memory=False)
    expr = np.log2(collapse_index(frame.T.apply(pd.to_numeric, errors="coerce")).clip(lower=0) + 1)
    meta = parse_series_matrix(raw / "GSE190266_series_matrix.txt.gz")
    meta["sample_id"] = meta["title"].astype(str)
    meta = meta.set_index("sample_id")
    common = [c for c in expr.columns if c in meta.index]
    clin = pd.DataFrame(index=common)
    pfs_t = None
    pfs_e = None
    for col in meta.columns:
        if "pfs_time" in col.lower():
            pfs_t = pd.to_numeric(meta.loc[common, col], errors="coerce")
        if "pfs_evt" in col.lower() or "pfs_event" in col.lower():
            pfs_e = pd.to_numeric(meta.loc[common, col], errors="coerce")
    # series-matrix keys from opus_lusc: pfs_time (6 months) / pfs_evt (6 months)
    if pfs_t is None:
        for key in meta.columns:
            if "pfs" in key.lower() and "time" in key.lower():
                pfs_t = pd.to_numeric(meta.loc[common, key], errors="coerce")
            if "pfs" in key.lower() and ("evt" in key.lower() or "event" in key.lower()):
                pfs_e = pd.to_numeric(meta.loc[common, key], errors="coerce")
    clin["pfs_time"] = pfs_t.values if pfs_t is not None else np.nan
    clin["pfs_event"] = pfs_e.values if pfs_e is not None else np.nan
    clin["pfs_unit"] = "month"
    clin["dcb"] = dcb_from_pfs(clin["pfs_time"], clin["pfs_event"], "month").values
    hist_col = "disease state" if "disease state" in meta.columns else None
    clin["histology"] = [norm_histology(meta.loc[s, hist_col] if hist_col else None) for s in common]
    clin["histology_source"] = "pathology_GEO"
    return expr[common], clin, "NSCLC", "lung", "DCB", "log2(TPM+1)"


def load_gse207422(raw: Path):
    md = pd.read_excel(raw / "GSE207422_metadata.xlsx")
    md = md.dropna(subset=["Sample"]).copy()
    expr = collapse_index(pd.read_csv(raw / "GSE207422_log2TPM.txt.gz", sep="\t", index_col=0))
    # Bulk xlsx, not the BD Rhapsody series-matrix titles. Prefer pre-treatment biopsy.
    if "Resource" in md.columns:
        pre = md["Resource"].astype(str).str.contains("Pre", case=False, na=False)
        if pre.any():
            md = md.loc[pre].copy()
    md = md[md["Sample"].astype(str).isin(expr.columns)]
    md = md.drop_duplicates("Sample")
    recist = md["RECIST"].astype(str).str.upper()
    rec = recode_recist(recist)
    orr = rec.copy()
    orr[orr.eq("SD")] = "NR"
    mpr = md["Pathologic Response"].astype(str).map(
        lambda x: "R" if str(x).upper().startswith("MPR") else ("NR" if str(x).upper().startswith("NMPR") else pd.NA)
    )
    clin = pd.DataFrame(index=md["Sample"].astype(str).values)
    clin["orr"] = orr.values
    clin["recist_raw"] = recist.values
    clin["author"] = mpr.values  # MPR stored separately; not mixed into ORR meta
    path_col = "Pathology" if "Pathology" in md.columns else None
    clin["histology"] = [norm_histology(md.set_index(md["Sample"].astype(str)).loc[s, path_col] if path_col else None) for s in clin.index]
    clin["histology_source"] = "pathology_xlsx"
    keep = [c for c in clin.index if c in expr.columns]
    return expr[keep], clin.loc[keep], "NSCLC", "lung", "ORR", "already_log2TPM"


def load_gse283829(raw: Path):
    counts = pd.read_csv(raw / "GSE283829_raw_express_matrix_all_samples.txt.gz", sep="\t", index_col=0)
    counts.columns = [str(c).lstrip("S") for c in counts.columns]
    counts = collapse_index(counts)
    cpm = counts.divide(counts.sum(axis=0), axis=1) * 1e6
    expr = np.log2(cpm.clip(lower=0) + 1)
    meta = parse_series_matrix(raw / "GSE283829_series_matrix.txt.gz")
    meta["key"] = meta["title"].astype(str)
    # titles may already match stripped columns
    recs = []
    for _, s in meta.iterrows():
        sid = str(s["title"])
        sid2 = sid.lstrip("S")
        raw_r = str(s.get("disease stage", s.get("response", "")))
        rec = recode_recist(pd.Series([raw_r])).iloc[0]
        orr = "R" if rec == "R" else ("NR" if rec in {"NR", "SD"} else pd.NA)
        recs.append(
            {
                "sample": sid if sid in expr.columns else sid2,
                "orr": orr,
                "recist_raw": raw_r,
                "histology": norm_histology(s.get("tumor type")),
                "histology_source": "pathology_GEO",
            }
        )
    clin = pd.DataFrame(recs).drop_duplicates("sample").set_index("sample")
    keep = [c for c in expr.columns if c in clin.index]
    return expr[keep], clin.loc[keep], "NSCLC", "lung", "ORR", "log2CPM"


def load_gse253564(raw: Path):
    frame = pd.read_csv(raw / "GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz", sep="\t")
    gene_col = "gene" if "gene" in frame.columns else frame.columns[0]
    frame = frame.drop(columns=[c for c in ("Entrez.ID",) if c in frame.columns])
    frame = frame.set_index(gene_col)
    expr = np.log2(collapse_index(frame).clip(lower=0) + 1)
    meta = parse_series_matrix(raw / "GSE253564_series_matrix.txt.gz")
    meta["key"] = meta["title"].astype(str)
    recs = []
    for _, s in meta.iterrows():
        recs.append(
            {
                "sample": str(s["title"]),
                "histology": norm_histology(s.get("cell type") or s.get("histology") or s.get("tumor type")),
                "histology_source": "pathology_GEO",
            }
        )
    clin = pd.DataFrame(recs).drop_duplicates("sample").set_index("sample")
    keep = [c for c in expr.columns if c in clin.index]
    if not keep:
        # titles may not match matrix columns; keep expression samples with NA clinical
        clin = pd.DataFrame(index=expr.columns)
        clin["histology"] = "NA"
        clin["histology_source"] = "unmatched_titles"
        keep = list(expr.columns)
    return expr[keep], clin.loc[keep], "NSCLC", "lung", "none_public", "log2(FPKM+1)"


def _icb_find_members(z: zipfile.ZipFile) -> tuple[str, str]:
    names = z.namelist()
    meta = next(n for n in names if n.endswith("metadata.tsv") or n.endswith("_metadata.csv"))
    expr_cands = [n for n in names if "expr" in n.lower() and n.endswith(".tsv")]
    # prefer gene-level tpm/fpkm over isoform
    expr_cands.sort(key=lambda n: (("gene" not in n.lower()), ("tpm" not in n.lower()), len(n)))
    if not expr_cands:
        raise FileNotFoundError(f"no expr tsv in {z.filename}")
    return meta, expr_cands[0]


def load_icb(raw: Path, zip_name: str, cohort: str, cancer: str, tissue: str):
    zpath = raw / zip_name
    with zipfile.ZipFile(zpath) as z:
        meta_n, expr_n = _icb_find_members(z)
        with z.open(meta_n) as fh:
            meta = pd.read_csv(io.TextIOWrapper(fh, encoding="utf-8"), sep="\t")
        with z.open(expr_n) as fh:
            expr = pd.read_csv(io.TextIOWrapper(fh, encoding="utf-8"), sep="\t", index_col=0)
    expr = collapse_index(expr)
    # already-log Zenodo objects keep sign; do not log2 again
    if float(np.nanmin(expr.to_numpy())) >= 0:
        expr = np.log2(expr.clip(lower=0) + 1)
        scale = "log2_plus1_from_nonneg"
    else:
        scale = "already_log_zenodo"
    pid_col = "patientid" if "patientid" in meta.columns else meta.columns[0]
    meta[pid_col] = meta[pid_col].astype(str)
    expr.columns = expr.columns.astype(str)
    common = [p for p in meta[pid_col] if p in expr.columns]
    # some ICB objects are genes x samples with patient ids
    if not common:
        # try samples as rows
        if set(meta[pid_col]).intersection(set(expr.index.astype(str))):
            expr = expr.T
            common = [p for p in meta[pid_col] if p in expr.columns]
    meta = meta[meta[pid_col].isin(common)].copy()
    meta = meta.drop_duplicates(pid_col)
    meta = meta.set_index(pid_col)
    common = [c for c in common if c in meta.index]
    clin = pd.DataFrame(index=common)
    if "response" in meta.columns:
        r = meta.loc[common, "response"].astype(str).str.upper()
        clin["author"] = r.map(lambda x: "R" if x == "R" else ("NR" if x == "NR" else pd.NA))
    if "recist" in meta.columns:
        rec = recode_recist(meta.loc[common, "recist"])
        orr = rec.copy()
        orr[orr.eq("SD")] = "NR"
        clin["orr"] = orr.values
        clin["recist_raw"] = meta.loc[common, "recist"].astype(str).values
        clin["orr_crpr_vs_pd"] = rec.replace({"SD": pd.NA}).values
    pfs_t = meta.loc[common, "survival_time_pfs"] if "survival_time_pfs" in meta.columns else None
    pfs_e = meta.loc[common, "event_occurred_pfs"] if "event_occurred_pfs" in meta.columns else None
    unit = ""
    if "survival_unit" in meta.columns:
        unit = str(meta.loc[common, "survival_unit"].dropna().astype(str).iloc[0]) if meta.loc[common, "survival_unit"].notna().any() else ""
    if pfs_t is not None:
        clin["pfs_time"] = pd.to_numeric(pfs_t, errors="coerce").values
        clin["pfs_event"] = pd.to_numeric(pfs_e, errors="coerce").values if pfs_e is not None else np.nan
        clin["pfs_unit"] = unit or "unknown"
        clin["dcb"] = dcb_from_pfs(clin["pfs_time"], clin["pfs_event"], unit).values
    clin["histology"] = "NA"
    clin["histology_source"] = "not_lung_or_absent"
    if "histo" in meta.columns:
        clin["histology"] = [norm_histology(v) for v in meta.loc[common, "histo"]]
        clin["histology_source"] = "icb_metadata"
    return expr[common], clin, cancer, tissue, "icb", scale


LOADERS = {
    "GSE126044": lambda raw: load_gse126044(raw),
    "GSE135222": lambda raw: load_gse135222(raw),
    "GSE166449": lambda raw: load_gse166449(raw),
    "GSE190265": lambda raw: load_gse190265(raw),
    "GSE190266": lambda raw: load_gse190266(raw),
    "GSE207422": lambda raw: load_gse207422(raw),
    "GSE283829": lambda raw: load_gse283829(raw),
    "GSE253564": lambda raw: load_gse253564(raw),
}

ICB_SPECS = [
    ("ICB_Mariathasan.zip", "IMvigor210", "urothelial", "nonlung"),
    ("ICB_Gide.zip", "Gide", "melanoma", "nonlung"),
    ("ICB_Liu.zip", "Liu", "melanoma", "nonlung"),
    ("ICB_Riaz.zip", "Riaz", "melanoma", "nonlung"),
    ("ICB_Braun.zip", "Braun", "ccRCC", "nonlung"),
    ("ICB_Jung.zip", "ICB_Jung", "NSCLC", "lung"),
    ("ICB_Hugo.zip", "Hugo", "melanoma", "nonlung"),
    ("ICB_Snyder.zip", "Snyder", "urothelial", "nonlung"),
    ("ICB_Shiuan.zip", "Shiuan", "RCC", "nonlung"),
    ("ICB_Puch.zip", "Puch", "melanoma", "nonlung"),
    ("ICB_Miao1.zip", "Miao1", "RCC", "nonlung"),
    ("ICB_Kim.zip", "Kim", "gastric", "nonlung"),
    ("ICB_Van_Allen.zip", "Van_Allen", "melanoma", "nonlung"),
    ("ICB_Nathanson.zip", "Nathanson", "melanoma", "nonlung"),
]


def analyze_binary(cohort, cancer, tissue, gene, measure, cutoff, endpoint, y, x) -> dict:
    base = {
        "cohort": cohort,
        "cancer": cancer,
        "tissue": tissue,
        "gene": gene,
        "measure": measure,
        "cutoff": cutoff,
        "endpoint": endpoint,
        "n_R": int((y == "R").sum()),
        "n_NR": int((y == "NR").sum()),
        "n_with_binary": int(y.isin(["R", "NR"]).sum()),
    }
    mask = y.isin(["R", "NR"]) & pd.to_numeric(x, errors="coerce").notna()
    yy, xx = y[mask], pd.to_numeric(x[mask], errors="coerce")
    if cutoff == "continuous":
        stats = logistic_or_per_sd(yy, xx)
        base.update(stats)
        if "n_R" not in stats:
            base["n_R"] = int((yy == "R").sum())
            base["n_NR"] = int((yy == "NR").sum())
        return base
    if yy.shape[0] < 6 or yy.nunique() < 2:
        base.update({"estimable": False, "reason": "too_few", "effect": "OR"})
        return base
    grp = assign_cutoff(xx, cutoff)
    a = int(((grp == "High") & (yy == "R")).sum())
    b = int(((grp == "High") & (yy == "NR")).sum())
    c = int(((grp == "Low") & (yy == "R")).sum())
    d = int(((grp == "Low") & (yy == "NR")).sum())
    if (a + b) == 0 or (c + d) == 0:
        base.update({"estimable": False, "reason": "empty_cutoff_arm", "n_high_R": a, "n_high_NR": b, "n_low_R": c, "n_low_NR": d, "effect": "OR"})
        return base
    base.update(or_from_2x2(a, b, c, d))
    return base


def forest(df: pd.DataFrame, meta: dict, title: str, out: Path, xlab: str) -> None:
    d = df.dropna(subset=["or", "or_lo", "or_hi"]).copy()
    if d.empty:
        return
    d = d.sort_values("cohort")
    labels = [
        f"{row['cohort']}  {int(row.get('n_R', 0))}/{int(row.get('n_NR', 0))}  "
        f"{row['or']:.2f} [{row['or_lo']:.2f}–{row['or_hi']:.2f}]"
        for _, row in d.iterrows()
    ]
    ys = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(9.5, 0.45 * (len(d) + 3) + 1.6))
    ax.errorbar(d["or"], ys, xerr=[d["or"] - d["or_lo"], d["or_hi"] - d["or"]], fmt="o", color="black", capsize=3)
    if meta.get("estimable"):
        ax.axvline(meta["re_or"], color="tab:blue", ls="--", label=f"RE {fmt_or(meta['re_or'], meta['re_or_lo'], meta['re_or_hi'])}")
        ax.axvspan(meta["re_or_lo"], meta["re_or_hi"], color="tab:blue", alpha=0.08)
    ax.axvline(1.0, color="grey", lw=1)
    ax.axvline(USER_OR, color="tab:red", ls=":", lw=1, label=f"user 0.42 [{USER_LO}–{USER_HI}]")
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel(xlab)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def write_writeup(out_dir: Path, inv: pd.DataFrame, per: pd.DataFrame, pooled: pd.DataFrame, cox: pd.DataFrame) -> None:
    def get(subset, gene, measure, cutoff, endpoint):
        m = pooled[
            (pooled.subset == subset)
            & (pooled.gene == gene)
            & (pooled.measure == measure)
            & (pooled.cutoff == cutoff)
            & (pooled.endpoint == endpoint)
        ]
        return m.iloc[0] if len(m) else None

    lines = []
    lines.append("# B5 wave 2: lung-first then all-open ICI, CLDN4 and TACSTD2")
    lines.append("")
    lines.append("## User-reported number (not produced here)")
    lines.append("")
    lines.append("CLDN4-high ICI **OR = 0.42 [0.18–0.95], k=11**.")
    lines.append("Prior open tests: IMvigor210 OR=1.47 p=0.214 (opposite); open lung ICI bulk")
    lines.append("all NS after purity residual (PR #114); LUSC-only ICI g=−0.29 p=0.46 (PR #205).")
    lines.append("This wave recomputes **lung first**, then all independent open ICI, for")
    lines.append("**CLDN4 and TACSTD2**, at median / tertile / continuous logistic, on raw")
    lines.append("expression **and** ESTIMATEScore residuals, with DCB and RECIST ORR kept")
    lines.append("as **separate** endpoints. It does not pick the cell closest to 0.42.")
    lines.append("")
    lines.append("## Locked methods")
    lines.append("")
    lines.append("- Public processed matrices only. No EGA. Liu/Braun raw are controlled;")
    lines.append("  processed Zenodo ICB TSVs (CC-BY-4.0, 10.5281/zenodo.7058399) are used.")
    lines.append("- ICB_Jung is the same patients as GSE135222 and is **not** an extra study.")
    lines.append("- DCB = PFS ≥ 6 months (180 days). Early-censored before 6 months = NA.")
    lines.append("- ORR = RECIST CR/PR vs SD+PD. Author/PredictIO R vs NR is a third label.")
    lines.append("- ESTIMATE: Yoshihara 2013 rank-ssGSEA (R estimate v1.0.13 port). Residual")
    lines.append("  = OLS of the gene on ESTIMATEScore. Residual skipped if a signature has")
    lines.append("  <30 genes overlapping the matrix.")
    lines.append("- Histology split only from pathology fields. GSE126044 / GSE135222 /")
    lines.append("  GSE166449 have no public histology and are not inferred here.")
    lines.append("- OR: Woolf + Haldane–Anscombe 0.5 if any 2×2 cell is 0; Fisher p.")
    lines.append("- Continuous: logistic OR per 1 SD. Survival: Cox HR per 1 SD and median.")
    lines.append("- Meta: DerSimonian–Laird RE + IVW FE on logOR / logHR.")
    lines.append("- GSE207422 is neoadjuvant ICI+chemo (RECIST ORR; MPR is not ORR).")
    lines.append("- GSE253564 has no public DCB/ORR in GEO; inventory only.")
    lines.append("- Van_Allen / Nathanson are CTLA-4 and are not in the PD-1/PD-L1 metas.")
    lines.append("")
    lines.append("## Inventory")
    lines.append("")
    lines.append("| Cohort | Cancer | CLDN4 | TACSTD2 | DCB n (R/NR) | ORR n (R/NR) | Author n | Histology | ESTIMATE overlap |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for _, r in inv.sort_values(["tissue", "cohort"]).iterrows():
        lines.append(
            f"| {r.cohort} | {r.cancer} | {r.cldn4} | {r.tacstd2} | {r.dcb_n} | {r.orr_n} | {r.author_n} | {r.histology_note} | {r.estimate_overlap} |"
        )
    lines.append("")
    lines.append("## Locked pooled results (residual, median)")
    lines.append("")
    lines.append("| Pool | Gene | Endpoint | k | n | RE OR [95% CI] | p | I² | matches 0.42 | 0.42 in CI |")
    lines.append("|---|---|---|---:|---:|---|---:|---:|---|---|")
    locked = [
        ("lung_DCB", "CLDN4", "DCB"),
        ("lung_ORR", "CLDN4", "ORR"),
        ("all_open_DCB", "CLDN4", "DCB"),
        ("all_open_ORR", "CLDN4", "ORR"),
        ("lung_DCB", "TACSTD2", "DCB"),
        ("lung_ORR", "TACSTD2", "ORR"),
        ("all_open_DCB", "TACSTD2", "DCB"),
        ("all_open_ORR", "TACSTD2", "ORR"),
    ]
    for subset, gene, end in locked:
        r = get(subset, gene, "residual", "median", end)
        if r is None or not r.estimable:
            lines.append(f"| {subset} | {gene} | {end} | — | — | not estimable | — | — | no | no |")
            continue
        match = abs(r.re_or - USER_OR) < 0.005
        inside = r.re_or_lo <= USER_OR <= r.re_or_hi
        lines.append(
            f"| {subset} | {gene} | {end} | {int(r.k)} | {int(r.n_total)} | "
            f"{fmt_or(r.re_or, r.re_or_lo, r.re_or_hi)} | {r.re_p:.3g} | {r.I2:.0f}% | "
            f"{'YES' if match else 'no'} | {'yes' if inside else 'no'} |"
        )
    lines.append("")
    lines.append("## All pooled cells (do not cherry-pick)")
    lines.append("")
    lines.append("Every gene × raw/residual × cutoff × endpoint × subset that was estimable.")
    lines.append("`matches_user_0.42_at_2dp` is true only when the RE point estimate rounds to 0.42.")
    lines.append("")
    lines.append("| Subset | Gene | Measure | Cutoff | Endpoint | k | RE [95% CI] | p | I² | 0.42@2dp | 0.42 in CI |")
    lines.append("|---|---|---|---|---|---:|---|---:|---:|---|---|")
    view = pooled[pooled.estimable == True].sort_values(["subset", "gene", "endpoint", "measure", "cutoff"])  # noqa: E712
    for _, r in view.iterrows():
        match = abs(r.re_or - USER_OR) < 0.005
        inside = r.re_or_lo <= USER_OR <= r.re_or_hi
        lines.append(
            f"| {r.subset} | {r.gene} | {r.measure} | {r.cutoff} | {r.endpoint} | {int(r.k)} | "
            f"{fmt_or(r.re_or, r.re_or_lo, r.re_or_hi)} | {r.re_p:.3g} | {r.I2:.0f}% | "
            f"{'YES' if match else 'no'} | {'yes' if inside else 'no'} |"
        )
    lines.append("")
    lines.append("## Per-cohort residual median (CLDN4 and TACSTD2)")
    lines.append("")
    lines.append("| Cohort | Gene | Endpoint | n (R/NR) | High R/NR | Low R/NR | OR [95% CI] | p |")
    lines.append("|---|---|---|---|---|---|---|---:|")
    sub = per[
        (per.measure == "residual")
        & (per.cutoff == "median")
        & (per.endpoint.isin(["DCB", "ORR", "author_R_vs_NR"]))
        & (~per.cohort.isin(DUPLICATE_OF))
    ]
    for _, r in sub.sort_values(["endpoint", "gene", "cohort"]).iterrows():
        if not r.get("estimable", False):
            lines.append(f"| {r.cohort} | {r.gene} | {r.endpoint} | {r.n_with_binary} ({r.n_R}/{r.n_NR}) | — | — | {r.get('reason', 'not estimable')} | — |")
            continue
        lines.append(
            f"| {r.cohort} | {r.gene} | {r.endpoint} | {int(r.n_used)} ({int(r.n_R)}/{int(r.n_NR)}) | "
            f"{int(r.n_high_R)}/{int(r.n_high_NR)} | {int(r.n_low_R)}/{int(r.n_low_NR)} | "
            f"{fmt_or(r['or'], r.or_lo, r.or_hi)} | {r.p:.3g} |"
        )
    lines.append("")
    lines.append("## Cox PFS (residual, continuous per 1 SD)")
    lines.append("")
    if cox is not None and len(cox):
        lines.append("| Cohort | Gene | n | events | HR [95% CI] | p |")
        lines.append("|---|---|---:|---:|---|---:|")
        csub = cox[(cox.measure == "residual") & (cox.cutoff == "continuous") & (~cox.cohort.isin(DUPLICATE_OF))]
        for _, r in csub.sort_values(["gene", "cohort"]).iterrows():
            if not r.get("estimable", False):
                lines.append(f"| {r.cohort} | {r.gene} | {r.n_used} | {r.get('n_events', '')} | not estimable | — |")
                continue
            lines.append(f"| {r.cohort} | {r.gene} | {int(r.n_used)} | {int(r.n_events)} | {fmt_or(r['or'], r.or_lo, r.or_hi)} | {r.p:.3g} |")
    else:
        lines.append("No Cox rows.")
    lines.append("")
    lines.append("## Honest verdict")
    lines.append("")
    prim = get("all_open_ORR", "CLDN4", "residual", "median", "ORR")
    lung = get("lung_DCB", "CLDN4", "residual", "median", "DCB")
    if prim is not None and prim.estimable:
        lines.append(
            f"- All-open ORR, CLDN4 residual median: RE OR = **{fmt_or(prim.re_or, prim.re_or_lo, prim.re_or_hi)}**, "
            f"k={int(prim.k)}, p={prim.re_p:.3g}. "
            + ("Rounds to 0.42." if abs(prim.re_or - USER_OR) < 0.005 else f"Does **not** match user 0.42 (observed {prim.re_or:.2f}).")
        )
        if int(prim.k) != USER_K:
            lines.append(f"- Public independent ORR k={int(prim.k)}, not the claimed k=11.")
    if lung is not None and lung.estimable:
        lines.append(
            f"- Lung-first DCB, CLDN4 residual median: RE OR = **{fmt_or(lung.re_or, lung.re_or_lo, lung.re_or_hi)}**, "
            f"k={int(lung.k)}, p={lung.re_p:.3g}."
        )
    lines.append("- We did not drop a cohort because it was null or opposite.")
    lines.append("- We did not add extra series after seeing results in order to reach k=11.")
    lines.append("- DCB and ORR were never pooled as if they were the same endpoint.")
    lines.append("- Cross-cancer all-open pooling is not a lung-specific test.")
    lines.append("- Small lung n makes tertile splits unstable; empty arms are not recoded.")
    lines.append("")
    lines.append("## Reproduction")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 -m pip install -r results/rework/B5_wave2/requirements.txt")
    lines.append("python3 results/rework/B5_wave2/download.py")
    lines.append("python3 results/rework/B5_wave2/analyze.py")
    lines.append("```")
    lines.append("")
    (out_dir / "WRITEUP.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def process_cohort(name, expr, clin, cancer, tissue, scale, gene_sets, processed: Path):
    est, overlaps = estimate_scores(expr, gene_sets)
    min_overlap = min(overlaps.values()) if overlaps else 0
    residual_ok = min_overlap >= 30
    genes_found = {}
    for g in ["CLDN4", "TACSTD2", "CD8A"]:
        rid = map_gene(expr, g)
        genes_found[g] = rid
    rows = []
    for g, rid in genes_found.items():
        if rid is None:
            continue
        raw = expr.loc[rid].reindex(clin.index)
        resid = pd.Series(residualize(raw.to_numpy(), est["ESTIMATEScore"].reindex(clin.index).to_numpy()), index=clin.index) if residual_ok else pd.Series(np.nan, index=clin.index)
        rows.append((g, raw, resid))
        clin[f"{g}_raw"] = raw.values
        clin[f"{g}_resid"] = resid.values
    clin["ESTIMATEScore"] = est["ESTIMATEScore"].reindex(clin.index).values
    clin["StromalScore"] = est["StromalScore"].reindex(clin.index).values
    clin["ImmuneScore"] = est["ImmuneScore"].reindex(clin.index).values
    clin["cohort"] = name
    clin["cancer"] = cancer
    clin["tissue"] = tissue
    clin["expr_scale"] = scale
    clin = clin.reset_index().rename(columns={"index": "sample"})
    clin.to_csv(processed / f"{name}.csv", index=False)

    def n_pair(col):
        if col not in clin.columns:
            return "—"
        s = clin[col]
        return f"{int(s.isin(['R','NR']).sum())} ({int((s=='R').sum())}/{int((s=='NR').sum())})"

    hist = clin["histology"] if "histology" in clin.columns else pd.Series(["NA"] * len(clin))
    hist_note = ",".join(f"{k}:{int(v)}" for k, v in hist.value_counts(dropna=False).items())
    inv = {
        "cohort": name,
        "cancer": cancer,
        "tissue": tissue,
        "n_samples": int(len(clin)),
        "cldn4": "yes" if genes_found["CLDN4"] else "no",
        "tacstd2": "yes" if genes_found["TACSTD2"] else "no",
        "dcb_n": n_pair("dcb"),
        "orr_n": n_pair("orr"),
        "author_n": n_pair("author"),
        "histology_note": hist_note,
        "histology_source": str(clin["histology_source"].iloc[0]) if "histology_source" in clin.columns else "",
        "estimate_overlap": f"stromal={overlaps.get('StromalSignature', 0)},immune={overlaps.get('ImmuneSignature', 0)}",
        "residual_ok": residual_ok,
        "expr_scale": scale,
        "duplicate_of": DUPLICATE_OF.get(name, ""),
    }

    or_rows = []
    cox_rows = []
    endpoints = []
    if "dcb" in clin.columns and clin["dcb"].isin(["R", "NR"]).sum() >= 6:
        endpoints.append(("DCB", clin["dcb"]))
    if "orr" in clin.columns and clin["orr"].isin(["R", "NR"]).sum() >= 6:
        endpoints.append(("ORR", clin["orr"]))
    if "author" in clin.columns and clin["author"].isin(["R", "NR"]).sum() >= 6:
        endpoints.append(("author_R_vs_NR", clin["author"]))

    histologies = [("all", pd.Series(True, index=clin.index))]
    if "histology" in clin.columns:
        for h in ["LUAD", "LUSC"]:
            m = clin["histology"] == h
            if m.sum() >= 8:
                histologies.append((h, m))
        non_lusc = clin["histology"].isin(["LUAD", "nonsquamous"])
        if non_lusc.sum() >= 8 and (clin["histology"] == "LUSC").sum() >= 1:
            histologies.append(("nonLUSC", non_lusc))

    for g, raw, resid in rows:
        if g == "CD8A":
            continue
        measures = [("raw", raw)]
        if residual_ok:
            measures.append(("residual", resid))
        for measure, vec in measures:
            x = pd.Series(vec.values, index=clin.index)
            for hname, hmask in histologies:
                for end_name, y in endpoints:
                    yy = y.where(hmask, other=pd.NA)
                    for cutoff in ["median", "tertile", "continuous"]:
                        rec = analyze_binary(name, cancer, tissue, g, measure, cutoff, end_name, yy, x)
                        rec["histology_split"] = hname
                        or_rows.append(rec)
            if "pfs_time" in clin.columns:
                for cutoff, xx, binary in [
                    ("continuous", x, False),
                    ("median", (assign_cutoff(x, "median") == "High").astype(float).where(assign_cutoff(x, "median").notna()), True),
                ]:
                    rec = cox_hr(clin["pfs_time"], clin.get("pfs_event"), xx, binary=binary)
                    rec.update({"cohort": name, "cancer": cancer, "gene": g, "measure": measure, "cutoff": cutoff, "endpoint": "PFS_Cox"})
                    cox_rows.append(rec)
    return inv, or_rows, cox_rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", default="/tmp/b5_wave2_raw")
    ap.add_argument("--out-dir", default=str(HERE))
    args = ap.parse_args()
    raw = Path(args.raw_dir)
    out_dir = Path(args.out_dir)
    tables = out_dir / "tables"
    figs = out_dir / "figures"
    processed = out_dir / "processed"
    for d in (tables, figs, processed):
        d.mkdir(parents=True, exist_ok=True)

    gene_sets = load_gmt(GMT_PATH)
    inventory = []
    or_rows = []
    cox_rows = []
    unusable = []

    for name, loader in LOADERS.items():
        print(f"===== {name} =====")
        try:
            expr, clin, cancer, tissue, _kind, scale = loader(raw)
        except Exception as e:
            print(f"FAIL {name}: {e}")
            unusable.append({"cohort": name, "reason": str(e)})
            continue
        print(f"  expr {expr.shape} scale={scale} clin={clin.shape}")
        inv, ors, coxs = process_cohort(name, expr, clin, cancer, tissue, scale, gene_sets, processed)
        inventory.append(inv)
        or_rows.extend(ors)
        cox_rows.extend(coxs)
        print(f"  CLDN4={inv['cldn4']} TACSTD2={inv['tacstd2']} residual_ok={inv['residual_ok']}")

    for zip_name, cohort, cancer, tissue in ICB_SPECS:
        print(f"===== {cohort} =====")
        try:
            expr, clin, cancer, tissue, _kind, scale = load_icb(raw, zip_name, cohort, cancer, tissue)
        except Exception as e:
            print(f"FAIL {cohort}: {e}")
            unusable.append({"cohort": cohort, "reason": str(e)})
            continue
        print(f"  expr {expr.shape} scale={scale} clin={clin.shape}")
        inv, ors, coxs = process_cohort(cohort, expr, clin, cancer, tissue, scale, gene_sets, processed)
        inventory.append(inv)
        or_rows.extend(ors)
        cox_rows.extend(coxs)
        print(f"  CLDN4={inv['cldn4']} TACSTD2={inv['tacstd2']} residual_ok={inv['residual_ok']}")

    # panel objects: document gene absence without loading as studies
    for zip_name, cohort in [("ICB_Hwang.zip", "Hwang"), ("ICB_Jerby_Arnon.zip", "Jerby_Arnon"), ("ICB_Roh.zip", "Roh")]:
        zpath = raw / zip_name
        if not zpath.exists():
            unusable.append({"cohort": cohort, "reason": "zip missing"})
            continue
        try:
            with zipfile.ZipFile(zpath) as z:
                _, expr_n = _icb_find_members(z)
                with z.open(expr_n) as fh:
                    expr = pd.read_csv(io.TextIOWrapper(fh, encoding="utf-8"), sep="\t", index_col=0, nrows=0)
                # read only index
                with z.open(expr_n) as fh:
                    idx = pd.read_csv(io.TextIOWrapper(fh, encoding="utf-8"), sep="\t", usecols=[0])
            idx = [str(x).split(".")[0].upper() for x in idx.iloc[:, 0]]
            has_c = ("CLDN4" in idx) or ("ENSG00000189143" in idx)
            has_t = ("TACSTD2" in idx) or ("ENSG00000184292" in idx)
            unusable.append({"cohort": cohort, "reason": f"panel; CLDN4={has_c} TACSTD2={has_t}"})
            print(f"{cohort} panel CLDN4={has_c} TACSTD2={has_t}")
        except Exception as e:
            unusable.append({"cohort": cohort, "reason": str(e)})

    inv_df = pd.DataFrame(inventory)
    per = pd.DataFrame(or_rows)
    cox = pd.DataFrame(cox_rows) if cox_rows else pd.DataFrame()
    inv_df.to_csv(tables / "cohort_inventory.csv", index=False)
    per.to_csv(tables / "per_cohort_or.csv", index=False)
    if len(cox):
        cox.to_csv(tables / "per_cohort_hr.csv", index=False)
    pd.DataFrame(unusable).to_csv(tables / "unusable.csv", index=False)

    # Primary rows are histology_split == all
    per_all = per[per.get("histology_split", "all") == "all"].copy() if "histology_split" in per.columns else per

    present = set(inv_df["cohort"])
    lung_dcb = [c for c in LUNG_DCB_CANDIDATES if c in present]
    lung_orr = [c for c in LUNG_ORR_CANDIDATES if c in present]
    lung_author = [c for c in LUNG_AUTHOR_CANDIDATES if c in present]
    nonlung = [c for c in NONLUNG_ICB if c in present]
    all_open_dcb = lung_dcb + [c for c in nonlung if c not in DUPLICATE_OF]
    all_open_orr = lung_orr + [c for c in nonlung if c not in DUPLICATE_OF]
    # IMvigor etc. contribute DCB only if they have DCB rows

    subsets = {
        "lung_DCB": (lung_dcb, "DCB"),
        "lung_ORR": (lung_orr, "ORR"),
        "lung_author": (lung_author, "author_R_vs_NR"),
        "all_open_DCB": (all_open_dcb, "DCB"),
        "all_open_ORR": (all_open_orr, "ORR"),
        "nonlung_ORR": (nonlung, "ORR"),
        "nonlung_DCB": (nonlung, "DCB"),
        "lung_LUAD_DCB": (lung_dcb, "DCB"),
        "lung_LUSC_DCB": (lung_dcb, "DCB"),
        "lung_nonLUSC_DCB": (lung_dcb, "DCB"),
        "lung_LUAD_ORR": (lung_orr, "ORR"),
        "lung_LUSC_ORR": (lung_orr, "ORR"),
        "lung_nonLUSC_ORR": (lung_orr, "ORR"),
    }

    pooled_rows = []
    for subset_name, (keep, endpoint) in subsets.items():
        hist_filter = "all"
        if "nonLUSC" in subset_name:
            hist_filter = "nonLUSC"
        elif "LUAD" in subset_name:
            hist_filter = "LUAD"
        elif "LUSC" in subset_name:
            hist_filter = "LUSC"
        src = per if hist_filter != "all" else per_all
        if hist_filter != "all" and "histology_split" in src.columns:
            src = src[src.histology_split == hist_filter]
        for gene in ["CLDN4", "TACSTD2"]:
            for measure in ["raw", "residual"]:
                for cutoff in ["median", "tertile", "continuous"]:
                    sub = src[
                        (src.cohort.isin(keep))
                        & (src.gene == gene)
                        & (src.measure == measure)
                        & (src.cutoff == cutoff)
                        & (src.endpoint == endpoint)
                        & (src.get("estimable") == True)  # noqa: E712
                    ]
                    meta = dl_meta(sub)
                    meta.update(
                        {
                            "subset": subset_name,
                            "gene": gene,
                            "measure": measure,
                            "cutoff": cutoff,
                            "endpoint": endpoint,
                            "histology_split": hist_filter,
                        }
                    )
                    pooled_rows.append(meta)
                    if (
                        meta.get("estimable")
                        and subset_name in {"lung_DCB", "lung_ORR", "all_open_DCB", "all_open_ORR"}
                        and measure == "residual"
                        and cutoff == "median"
                    ):
                        title = f"{subset_name} | {gene} residual median | {endpoint} | k={meta['k']}"
                        forest(
                            sub,
                            meta,
                            title,
                            figs / f"forest_{subset_name}_{gene}_residual_median_{endpoint}.png",
                            f"OR (response, {gene}-high vs low; residual)",
                        )

    # Cox metas
    if len(cox):
        for subset_name, keep in {
            "lung_PFS": lung_dcb,
            "all_open_PFS": lung_dcb + nonlung,
        }.items():
            for gene in ["CLDN4", "TACSTD2"]:
                for measure in ["raw", "residual"]:
                    for cutoff in ["continuous", "median"]:
                        sub = cox[
                            (cox.cohort.isin(keep))
                            & (cox.gene == gene)
                            & (cox.measure == measure)
                            & (cox.cutoff == cutoff)
                            & (cox.get("estimable") == True)  # noqa: E712
                        ]
                        meta = dl_meta(sub)
                        meta.update(
                            {
                                "subset": subset_name,
                                "gene": gene,
                                "measure": measure,
                                "cutoff": cutoff,
                                "endpoint": "PFS_Cox",
                                "histology_split": "all",
                            }
                        )
                        pooled_rows.append(meta)

    pooled = pd.DataFrame(pooled_rows)
    pooled.to_csv(tables / "pooled.csv", index=False)
    est = pooled[pooled.get("estimable") == True].copy()  # noqa: E712
    if len(est):
        est["matches_user_0.42_at_2dp"] = (est["re_or"] - USER_OR).abs() < 0.005
        est["user_0.42_in_re_ci"] = (est["re_or_lo"] <= USER_OR) & (est["re_or_hi"] >= USER_OR)
        compact_cols = [
            "subset",
            "gene",
            "measure",
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
            "matches_user_0.42_at_2dp",
            "user_0.42_in_re_ci",
            "cohorts",
        ]
        est[[c for c in compact_cols if c in est.columns]].to_csv(tables / "pooled_compact.csv", index=False)
        nearest = est.assign(abs_diff=(est["re_or"] - USER_OR).abs()).sort_values("abs_diff").iloc[0]
        nearest_rec = {
            "subset": nearest.subset,
            "gene": nearest.gene,
            "measure": nearest.measure,
            "cutoff": nearest.cutoff,
            "endpoint": nearest.endpoint,
            "re_or": float(nearest.re_or),
            "k": int(nearest.k),
            "abs_diff": float(nearest.abs_diff),
            "adopted": False,
            "note": "Listed only to document that we saw the closest cell and did not select it.",
        }
    else:
        nearest_rec = {}

    def pack(subset, gene, measure, cutoff, endpoint):
        m = pooled[
            (pooled.subset == subset)
            & (pooled.gene == gene)
            & (pooled.measure == measure)
            & (pooled.cutoff == cutoff)
            & (pooled.endpoint == endpoint)
        ]
        if not len(m):
            return {}
        r = m.iloc[0]
        return {k: r.get(k) for k in ["k", "n_total", "re_or", "re_or_lo", "re_or_hi", "re_p", "I2", "cohorts"]}

    summary = {
        "user_claim": {"or": USER_OR, "lo": USER_LO, "hi": USER_HI, "k": USER_K},
        "locked_primary": {
            "lung_DCB_CLDN4_residual_median": pack("lung_DCB", "CLDN4", "residual", "median", "DCB"),
            "lung_ORR_CLDN4_residual_median": pack("lung_ORR", "CLDN4", "residual", "median", "ORR"),
            "all_open_DCB_CLDN4_residual_median": pack("all_open_DCB", "CLDN4", "residual", "median", "DCB"),
            "all_open_ORR_CLDN4_residual_median": pack("all_open_ORR", "CLDN4", "residual", "median", "ORR"),
            "lung_DCB_TACSTD2_residual_median": pack("lung_DCB", "TACSTD2", "residual", "median", "DCB"),
            "lung_ORR_TACSTD2_residual_median": pack("lung_ORR", "TACSTD2", "residual", "median", "ORR"),
            "all_open_DCB_TACSTD2_residual_median": pack("all_open_DCB", "TACSTD2", "residual", "median", "DCB"),
            "all_open_ORR_TACSTD2_residual_median": pack("all_open_ORR", "TACSTD2", "residual", "median", "ORR"),
        },
        "closest_to_user_0.42_not_adopted": nearest_rec,
        "did_we_pick_cutoff_to_hit_0.42": False,
        "independent_note": "ICB_Jung excluded as duplicate of GSE135222; CTLA-4 cohorts excluded from PD-1 metas",
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=float), encoding="utf-8")
    write_writeup(out_dir, inv_df, per_all, pooled, cox)
    print("wrote", out_dir / "WRITEUP.md")
    print(json.dumps(summary["locked_primary"], indent=2, default=float))
    print("closest_not_adopted", json.dumps(nearest_rec, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
