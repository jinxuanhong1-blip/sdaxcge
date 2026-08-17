#!/usr/bin/env python3
"""GSE7670 paired LUAD Affy U133A: tumor CLDN4 vs CD8A / CD274 / ImmuneScore.

Additive public bulk slice. CLDN4-only (no dual-high). Unit is the array.
Primary n is individual LUAD tumors. Matched adjacent-normal arrays are a
companion only and are not mixed into the tumor pairwise tests. Downloads
stay under $GSE7670_CLDN4_DATA (default /tmp/gse7670_cldn4) and are not
committed.
"""

from __future__ import annotations

import gzip
import io
import json
import os
import re
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
DATA = Path(os.environ.get("GSE7670_CLDN4_DATA", "/tmp/gse7670_cldn4"))
SEED = 20260817
N_BOOT = 2000

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7670/"
    "matrix/GSE7670_series_matrix.txt.gz"
)
GPL_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz"

# A1 8-gene T-cell effector list (PR 139 / OncoSG). Not ESTIMATE.
IMMUNE8 = ["CD8A", "GZMA", "GZMB", "IFNG", "EOMES", "CXCL9", "CXCL10", "TBX21"]
EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]
NAMED_PROBES = {
    "CLDN4": "201428_at",
    "CD8A": "205758_at",
    "TACSTD2": "202286_s_at",
    "PDCD1LG2": "220049_s_at",
}
# Plus-2 CD274 probes are not on GPL96; listed only to document ABSENT.
CD274_PLUS2 = ["223834_at", "227458_at"]

# Same HOLDS rule as PR 229 / GSE4573: n>=40, adj ρ<0, p<0.05.
HOLDS_N = 40


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "gse7670-cldn4/1.0"})
    with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as out:
        out.write(r.read())
    return dest


def parse_soft_table(text: str) -> pd.DataFrame:
    lines = text.splitlines()
    start = next(
        i
        for i, l in enumerate(lines)
        if l.startswith("!platform_table_begin") or l.startswith("!series_matrix_table_begin")
    ) + 1
    end = next(
        i
        for i, l in enumerate(lines)
        if l.startswith("!platform_table_end") or l.startswith("!series_matrix_table_end")
    )
    return pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", dtype=str)


def parse_series_meta(path: Path) -> tuple[pd.DataFrame, dict]:
    series: dict[str, list[str]] = {}
    sample_fields: dict[str, list[list[str]]] = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Series_"):
                k = line.split("\t", 1)[0][8:]
                v = line.split("\t", 1)[1].strip().strip('"') if "\t" in line else ""
                series.setdefault(k, []).append(v)
            elif line.startswith("!Sample_"):
                key = line.split("\t", 1)[0][8:]
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                sample_fields.setdefault(key, []).append(vals)

    first = {k: rows[0] for k, rows in sample_fields.items()}
    meta = pd.DataFrame(first)
    if "geo_accession" not in meta.columns:
        raise KeyError(f"no geo_accession in sample fields: {list(meta.columns)}")
    meta.index = meta["geo_accession"].astype(str)
    meta.index.name = "gsm"

    char_rows = sample_fields.get("characteristics_ch1", [])
    if char_rows:
        meta["tissue_phrase"] = char_rows[0]
    if len(char_rows) > 1:
        meta["sex_raw"] = char_rows[1]
    else:
        meta["sex_raw"] = ""
    return meta, {k: " | ".join(v) for k, v in series.items()}


def first_symbol(s: str) -> str:
    if not isinstance(s, str) or not s or s == "nan":
        return ""
    return s.split(" /// ")[0].strip()


def classify_sample(row: pd.Series) -> dict:
    src = str(row.get("source_name_ch1", "")).strip()
    title = str(row.get("title", "")).strip()
    phrase = str(row.get("tissue_phrase", "")).strip()
    sex_raw = str(row.get("sex_raw", "")).strip()
    sex = ""
    msex = re.search(r"sex:\s*(.+)$", sex_raw, flags=re.I)
    if msex:
        sex = msex.group(1).strip().lower()

    pair_id = ""
    role = "other"
    histology = "unknown"
    blob = f"{src} {title} {phrase}".lower()

    mt = re.fullmatch(r"(\d+)T", src)
    mn = re.fullmatch(r"(\d+)N", src)
    if mt:
        pair_id = mt.group(1)
        role = "patient_tumor"
    elif mn:
        pair_id = mn.group(1)
        role = "patient_normal"
    elif "mixture" in src.lower() or "mixture" in title.lower():
        role = "mixture_tumor" if src.startswith("T.") or title.lower().startswith("t. mixture") else "mixture_normal"
    elif "stratagene" in src.lower() or "clontech" in src.lower():
        role = "commercial_normal"
    elif src in {"A549", "H661", "H1299", "NL-20", "CL1-0", "CL1-1", "CL1-5", "CL1-5-F4"}:
        role = "cell_line"
    else:
        role = "other"

    if "lage cell" in blob or "large cell" in blob:
        histology = "LCC"
    elif "adenocarcinoma" in blob:
        histology = "LUAD"
    elif role == "cell_line":
        histology = {
            "A549": "LUAD_line",
            "H661": "LCC_line",
            "H1299": "NSCLC_line",
            "NL-20": "bronchial_line",
            "CL1-0": "LUAD_line",
            "CL1-1": "LUAD_line",
            "CL1-5": "LUAD_line",
            "CL1-5-F4": "LUAD_line",
        }.get(src, "cell_line")
    elif role.startswith("commercial"):
        histology = "normal_commercial"
    elif role.startswith("mixture"):
        histology = "LUAD_mixture"

    primary_tumor = bool(role == "patient_tumor" and histology == "LUAD")
    companion_normal = bool(role == "patient_normal" and histology == "LUAD")
    return {
        "source": src,
        "pair_id": pair_id,
        "role": role,
        "histology": histology,
        "sex": sex,
        "primary_luad_tumor": primary_tumor,
        "companion_luad_normal": companion_normal,
    }


def collapse_maxmean(probe_expr: pd.DataFrame, id2gene: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = probe_expr.index.map(lambda x: id2gene.get(str(x), ""))
    df = probe_expr.copy()
    df["_gene"] = genes
    df = df[df["_gene"].astype(str).str.len() > 0]
    means = df.drop(columns="_gene").astype(float).mean(axis=1)
    df["_mean"] = means
    df = df.sort_values("_mean", ascending=False)
    kept = df.loc[~df["_gene"].duplicated(keep="first")]
    gene_expr = kept.drop(columns=["_gene", "_mean"]).astype(float)
    gene_expr.index = kept["_gene"].values
    gene_expr = gene_expr.sort_index()
    tmp = df.reset_index()
    probe_col = tmp.columns[0]
    probe_audit = (
        tmp.rename(columns={probe_col: "probe", "_gene": "gene", "_mean": "mean"})
        [["gene", "probe", "mean"]]
        .sort_values(["gene", "mean"], ascending=[True, False])
    )
    return gene_expr, probe_audit


def mean_z(expr: pd.DataFrame, genes: list[str], samples: pd.Index) -> tuple[pd.Series, list[str], list[str]]:
    present = [g for g in genes if g in expr.index]
    missing = [g for g in genes if g not in expr.index]
    if not present:
        return pd.Series(np.nan, index=samples), present, missing
    block = expr.loc[present, samples].astype(float)
    z = block.apply(lambda r: (r - r.mean()) / (r.std(ddof=0) if r.std(ddof=0) else np.nan), axis=1)
    return z.mean(axis=0), present, missing


def spearman_ci(x, y, n_boot: int = N_BOOT, seed: int = SEED):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = int(x.size)
    if n < 8:
        return dict(n=n, rho=np.nan, p=np.nan, ci_low=np.nan, ci_high=np.nan)
    rho, p = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        if np.unique(x[idx]).size < 3 or np.unique(y[idx]).size < 3:
            boots[i] = np.nan
            continue
        boots[i] = stats.spearmanr(x[idx], y[idx])[0]
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return dict(n=n, rho=float(rho), p=float(p), ci_low=float(lo), ci_high=float(hi))


def partial_spearman(x, y, covar):
    x = pd.Series(np.asarray(x, float))
    y = pd.Series(np.asarray(y, float))
    c = pd.Series(np.asarray(covar, float), name="cov")
    df = pd.concat([x.rename("x"), y.rename("y"), c], axis=1).dropna()
    n = int(df.shape[0])
    if n < 12:
        return dict(n=n, rho_adj=np.nan, p_adj=np.nan)
    xr = stats.rankdata(df["x"].values)
    yr = stats.rankdata(df["y"].values)
    cr = np.column_stack([np.ones(n), stats.rankdata(df["cov"].values)])
    bx, *_ = np.linalg.lstsq(cr, xr, rcond=None)
    by, *_ = np.linalg.lstsq(cr, yr, rcond=None)
    rx = xr - cr @ bx
    ry = yr - cr @ by
    if np.std(rx) == 0 or np.std(ry) == 0:
        return dict(n=n, rho_adj=np.nan, p_adj=np.nan)
    r = float(np.corrcoef(rx, ry)[0, 1])
    dof = n - 3
    t = r * np.sqrt(dof / max(1e-12, 1 - r ** 2))
    p = float(2 * stats.t.sf(abs(t), dof))
    return dict(n=n, rho_adj=r, p_adj=p)


def verdict(rho_adj, p_adj, n, present: bool = True) -> str:
    if not present:
        return "ABSENT"
    if n < HOLDS_N:
        return "UNDERPOWERED"
    if np.isfinite(rho_adj) and rho_adj < 0 and p_adj < 0.05:
        return "HOLDS"
    if np.isfinite(rho_adj) and rho_adj > 0 and p_adj < 0.05:
        return "OPPOSITE"
    return "NO_EVIDENCE"


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_rho(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:+.3f}"


def pair_row(x, y, covar, gene_x, gene_y, subset, n_arrays, extra=None, present=True):
    if not present:
        rec = {
            "gene_x": gene_x,
            "gene_y": gene_y,
            "subset": subset,
            "n": 0,
            "n_arrays": n_arrays,
            "rho": np.nan,
            "p": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "rho_adj_epithelial": np.nan,
            "p_adj_epithelial": np.nan,
            "verdict": "ABSENT",
        }
        if extra:
            rec.update(extra)
        return rec
    crude = spearman_ci(x, y)
    adj = partial_spearman(x, y, covar)
    rec = {
        "gene_x": gene_x,
        "gene_y": gene_y,
        "subset": subset,
        "n": crude["n"],
        "n_arrays": n_arrays,
        "rho": crude["rho"],
        "p": crude["p"],
        "ci_low": crude["ci_low"],
        "ci_high": crude["ci_high"],
        "rho_adj_epithelial": adj["rho_adj"],
        "p_adj_epithelial": adj["p_adj"],
        "verdict": verdict(adj["rho_adj"], adj["p_adj"], adj["n"], present=True),
    }
    if extra:
        rec.update(extra)
    return rec


def highlow(x, y, endpoint: str) -> dict:
    q = pd.qcut(pd.Series(x), 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    a = pd.Series(y)[q == "Q4"]
    b = pd.Series(y)[q == "Q1"]
    if a.size < 3 or b.size < 3:
        return {
            "endpoint": endpoint,
            "n_Q4": int(a.size),
            "n_Q1": int(b.size),
            "median_Q4": float(np.median(a)) if a.size else np.nan,
            "median_Q1": float(np.median(b)) if b.size else np.nan,
            "U": np.nan,
            "p": np.nan,
            "rank_biserial_Q4_minus_Q1": np.nan,
            "note": "too few per quartile",
        }
    U, p_mw = stats.mannwhitneyu(a.values, b.values, alternative="two-sided")
    rbc = 2 * U / (a.size * b.size) - 1
    return {
        "endpoint": endpoint,
        "n_Q4": int(a.size),
        "n_Q1": int(b.size),
        "median_Q4": float(np.median(a)),
        "median_Q1": float(np.median(b)),
        "U": float(U),
        "p": float(p_mw),
        "rank_biserial_Q4_minus_Q1": float(rbc),
        "note": "",
    }


def paired_wilcoxon(tumor: pd.Series, normal: pd.Series, pair_map: pd.DataFrame, endpoint: str) -> dict:
    rows = []
    for _, r in pair_map.iterrows():
        t = tumor.get(r["tumor_gsm"], np.nan)
        n = normal.get(r["normal_gsm"], np.nan)
        if np.isfinite(t) and np.isfinite(n):
            rows.append((t, n))
    if not rows:
        return {"endpoint": endpoint, "n_pairs": 0, "n_tumor_gt_normal": 0, "median_T": np.nan,
                "median_N": np.nan, "median_T_minus_N": np.nan, "wilcoxon_stat": np.nan, "p": np.nan}
    t = np.array([a for a, _ in rows], float)
    n = np.array([b for _, b in rows], float)
    d = t - n
    try:
        w, p = stats.wilcoxon(t, n, alternative="two-sided", zero_method="wilcox")
    except ValueError:
        w, p = np.nan, np.nan
    return {
        "endpoint": endpoint,
        "n_pairs": int(len(rows)),
        "n_tumor_gt_normal": int((d > 0).sum()),
        "median_T": float(np.median(t)),
        "median_N": float(np.median(n)),
        "median_T_minus_N": float(np.median(d)),
        "wilcoxon_stat": float(w) if np.isfinite(w) else np.nan,
        "p": float(p) if np.isfinite(p) else np.nan,
    }


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    matrix_path = dl(DATA / "GSE7670_series_matrix.txt.gz", GEO_MATRIX)
    annot_path = dl(DATA / "GPL96.annot.gz", GPL_ANNOT)

    with gzip.open(matrix_path, "rt", errors="replace") as fh:
        raw = fh.read()
    meta, series = parse_series_meta(matrix_path)
    lines = raw.split("\n")
    start = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin")) + 1
    end = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end"))
    probe_expr = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", index_col=0)
    probe_expr.index = probe_expr.index.astype(str).str.strip('"')
    probe_expr = probe_expr.astype(float)
    probe_expr = probe_expr.loc[:, [c for c in probe_expr.columns if c in meta.index]]
    meta = meta.loc[probe_expr.columns]

    parsed = meta.apply(classify_sample, axis=1, result_type="expand")
    meta = pd.concat([meta, parsed], axis=1)

    with gzip.open(annot_path, "rt", errors="replace") as fh:
        annot = parse_soft_table(fh.read())
    annot["symbol"] = annot["Gene symbol"].map(first_symbol)
    id2gene = annot.drop_duplicates("ID").set_index("ID")["symbol"]
    entrez = annot.drop_duplicates("ID").set_index("ID")["Gene ID"].astype(str)

    gene_expr, probe_audit = collapse_maxmean(probe_expr, id2gene)
    n_missing = int(probe_expr.isna().sum().sum())
    n_array = int(probe_expr.shape[1])

    tumors = meta.index[meta["primary_luad_tumor"]]
    normals = meta.index[meta["companion_luad_normal"]]
    all_patient_tumors = meta.index[meta["role"] == "patient_tumor"]
    n_luad = int(len(tumors))
    n_lcc_tumor = int(((meta["role"] == "patient_tumor") & (meta["histology"] == "LCC")).sum())
    n_pairs_luad = n_luad
    n_cell = int((meta["role"] == "cell_line").sum())
    n_mix = int(meta["role"].str.startswith("mixture").sum())
    n_comm = int((meta["role"] == "commercial_normal").sum())

    if n_luad != 26:
        raise SystemExit(f"expected 26 LUAD tumors, got {n_luad}")
    if n_lcc_tumor != 1:
        raise SystemExit(f"expected 1 LCC tumor (pair 19), got {n_lcc_tumor}")

    # Pair map for LUAD companion T vs N
    pair_rows = []
    for pid, sub in meta[meta["pair_id"] != ""].groupby("pair_id"):
        t = sub.index[sub["role"] == "patient_tumor"]
        n = sub.index[sub["role"] == "patient_normal"]
        hist = sub.loc[t[0], "histology"] if len(t) else "unknown"
        if len(t) == 1 and len(n) == 1:
            pair_rows.append({
                "pair_id": pid,
                "tumor_gsm": t[0],
                "normal_gsm": n[0],
                "histology": hist,
                "sex": sub.loc[t[0], "sex"],
            })
    pair_map = pd.DataFrame(pair_rows).sort_values("pair_id")
    luad_pairs = pair_map[pair_map["histology"] == "LUAD"].copy()

    # Scores on LUAD tumors only (z within the primary set)
    epi, epi_present, epi_missing = mean_z(gene_expr, EPITHELIAL, tumors)
    imm, imm_present, imm_missing = mean_z(gene_expr, IMMUNE8, tumors)
    # Companion normals: score genes with the tumor-defined z is wrong; for T vs N
    # use raw gene values / a separately computed normal+tumor z is not needed.
    # For paired CLDN4 / CD8A use raw MAS5.

    cldn4_all = gene_expr.loc["CLDN4"]
    cd8a_all = gene_expr.loc["CD8A"]
    cldn4 = cldn4_all.reindex(tumors)
    cd8a = cd8a_all.reindex(tumors)
    cd274_present = "CD274" in gene_expr.index
    tacstd2_present = "TACSTD2" in gene_expr.index

    # Probe confirm
    def probe_for(gene: str) -> str:
        hit = probe_audit.loc[probe_audit["gene"] == gene, "probe"]
        return str(hit.iloc[0]) if len(hit) else ""

    probe_rows = []
    for gene, named in NAMED_PROBES.items():
        on_chip = named in probe_expr.index
        entrez_id = entrez.get(named, "") if on_chip else ""
        probe_rows.append({
            "gene": gene,
            "named_probe": named,
            "on_GPL96": on_chip,
            "collapse_probe": probe_for(gene),
            "entrez": entrez_id,
            "finite_on_LUAD_tumors": int(probe_expr.loc[named, tumors].notna().all()) if on_chip else 0,
        })
    for p in CD274_PLUS2:
        probe_rows.append({
            "gene": "CD274",
            "named_probe": p,
            "on_GPL96": p in probe_expr.index,
            "collapse_probe": "",
            "entrez": "29126",
            "finite_on_LUAD_tumors": 0,
        })
    probe_confirm = pd.DataFrame(probe_rows)
    probe_confirm.to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)

    # Label inventory — honest n
    labels = pd.DataFrame([
        {"field": "arrays in series matrix", "public": "yes", "n": n_array,
         "note": "22,283 GPL96 probes × 66 GSM; 0 missing MAS5; do not use 66 for the LUAD tests"},
        {"field": "unique GSM / unique titles", "public": "yes", "n": n_array, "note": "all unique"},
        {"field": "patients (GEO overall_design text)", "public": "text only", "n": 27,
         "note": "pairwise samples from 27 patients at Taipei VGH; one pair is LCC not LUAD"},
        {"field": "patient tumor arrays", "public": "yes", "n": int(len(all_patient_tumors)),
         "note": "26 LUAD + 1 LCC (pair 19, title typo 'lage cell')"},
        {"field": "LUAD patient tumors (primary n)", "public": "yes", "n": n_luad,
         "note": "source *T and adenocarcinoma; excludes pair 19, mixtures, cell lines"},
        {"field": "LUAD matched adjacent-normal (companion)", "public": "yes", "n": int(len(normals)),
         "note": "not mixed into the tumor pairwise tests"},
        {"field": "LCC pair 19 tumor+normal", "public": "yes", "n": 2,
         "note": "inventoried; not LUAD; not in primary n"},
        {"field": "tissue mixtures (Taichung VGH)", "public": "yes", "n": n_mix,
         "note": "Adj. N. mixture + T. mixture; not an individual tumor"},
        {"field": "commercial normal lung", "public": "yes", "n": n_comm,
         "note": "Stratagene 735020 + Clontech 636524"},
        {"field": "cell lines", "public": "yes", "n": n_cell,
         "note": "A549 H661 H1299 NL-20 CL1-0 CL1-1 CL1-5 CL1-5-F4; not primary"},
        {"field": "stage", "public": "no", "n": 0,
         "note": "overall_design says early and late stages; not a per-array characteristic"},
        {"field": "ICI / treatment", "public": "no", "n": 0, "note": "2007 surgical atlas, not an ICI series"},
        {"field": "tumor % / ABSOLUTE purity", "public": "no", "n": 0,
         "note": "only public proxy used here is an RNA epithelial mean-z"},
        {"field": "ESTIMATE ImmuneScore", "public": "no", "n": 0,
         "note": "skipped — A1 8-gene mean-z used instead; not Yoshihara ESTIMATE"},
        {"field": "CLDN4 finite (201428_at)", "public": "yes", "n": n_luad, "note": "named GPL96 / Entrez 1364"},
        {"field": "CD8A finite (205758_at)", "public": "yes", "n": n_luad, "note": "named GPL96 / Entrez 925"},
        {"field": "CD274 finite (Entrez 29126)", "public": "no", "n": 0,
         "note": "no GPL96 probe; Plus-2 223834_at / 227458_at are not on U133A"},
        {"field": "ImmuneScore A1 genes present", "public": "yes", "n": len(imm_present),
         "note": f"{len(imm_present)}/8; missing {','.join(imm_missing) if imm_missing else 'none'}"},
        {"field": "epithelial genes present", "public": "yes", "n": len(epi_present),
         "note": f"{len(epi_present)}/6; missing {','.join(epi_missing) if epi_missing else 'none'}"},
        {"field": "primary pairwise n (CLDN4 + CD8A)", "public": "yes", "n": n_luad,
         "note": "this is the n used below"},
        {"field": "primary pairwise n (CLDN4 + ImmuneScore)", "public": "yes", "n": n_luad,
         "note": "A1 mean-z on genes present; EOMES absent"},
        {"field": "primary pairwise n (CLDN4 + CD274)", "public": "no", "n": 0,
         "note": "ABSENT — do not invent a PD-L1 ρ"},
    ])
    labels.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    for g in ["CLDN4", "CD8A"]:
        if g not in gene_expr.index:
            raise SystemExit(f"required gene absent after collapse: {g}")
    if cd274_present:
        raise SystemExit("CD274 unexpectedly present on GPL96; re-check annotation")

    rows = []
    subsets = {
        "LUAD_tumors": tumors,
        "all_patient_tumors_26LUAD_1LCC": all_patient_tumors,
    }
    # epithelial / immune for the mixed 27: recompute z within that subset
    epi_27, _, _ = mean_z(gene_expr, EPITHELIAL, all_patient_tumors)
    imm_27, _, _ = mean_z(gene_expr, IMMUNE8, all_patient_tumors)

    for subset_name, idx in subsets.items():
        epi_s = epi if subset_name == "LUAD_tumors" else epi_27
        imm_s = imm if subset_name == "LUAD_tumors" else imm_27
        role = "primary" if subset_name == "LUAD_tumors" else "sensitivity_include_LCC"
        rows.append(pair_row(
            cldn4_all.reindex(idx).values, cd8a_all.reindex(idx).values, epi_s.reindex(idx).values,
            "CLDN4", "CD8A", subset_name, n_arrays=int(len(idx)), extra={"role": role},
        ))
        rows.append(pair_row(
            cldn4_all.reindex(idx).values, imm_s.reindex(idx).values, epi_s.reindex(idx).values,
            "CLDN4", "ImmuneScore", subset_name, n_arrays=int(len(idx)),
            extra={"role": role, "note": f"A1 {len(imm_present)}/8 mean-z; missing {','.join(imm_missing)}"},
        ))
        rows.append(pair_row(
            None, None, None, "CLDN4", "CD274", subset_name, n_arrays=int(len(idx)),
            extra={"role": role, "note": "Entrez 29126 not on GPL96"}, present=False,
        ))
        rows.append(pair_row(
            cldn4_all.reindex(idx).values, epi_s.reindex(idx).values, epi_s.reindex(idx).values,
            "CLDN4", "epithelial_mean_z", subset_name, n_arrays=int(len(idx)),
            extra={"role": "covariate"},
        ))
        rows[-1]["rho_adj_epithelial"] = np.nan
        rows[-1]["p_adj_epithelial"] = np.nan
        rows[-1]["verdict"] = "COMPANION"
        if tacstd2_present:
            rows.append(pair_row(
                cldn4_all.reindex(idx).values, gene_expr.loc["TACSTD2"].reindex(idx).values,
                epi_s.reindex(idx).values, "CLDN4", "TACSTD2", subset_name,
                n_arrays=int(len(idx)), extra={"role": "companion", "note": "not dual-high"},
            ))
            rows.append(pair_row(
                gene_expr.loc["TACSTD2"].reindex(idx).values, cd8a_all.reindex(idx).values,
                epi_s.reindex(idx).values, "TACSTD2", "CD8A", subset_name,
                n_arrays=int(len(idx)), extra={"role": "companion", "note": "not this claim"},
            ))
        rows.append(pair_row(
            cd8a_all.reindex(idx).values, imm_s.reindex(idx).values, epi_s.reindex(idx).values,
            "CD8A", "ImmuneScore", subset_name, n_arrays=int(len(idx)),
            extra={"role": "positive_control", "note": "CD8A is one of the ImmuneScore genes"},
        ))

    # PDCD1LG2 is not CD274 — document only
    if "PDCD1LG2" in gene_expr.index:
        rows.append(pair_row(
            cldn4.values, gene_expr.loc["PDCD1LG2"].reindex(tumors).values, epi.values,
            "CLDN4", "PDCD1LG2", "LUAD_tumors", n_arrays=n_luad,
            extra={"role": "not_CD274", "note": "PD-L2; not substituted for CD274"},
        ))

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(TABLES / "spearman_cldn4_vs_genes.tsv", sep="\t", index=False)

    hl = pd.DataFrame([
        highlow(cldn4.values, cd8a.values, "CD8A"),
        highlow(cldn4.values, imm.values, "ImmuneScore"),
        {"endpoint": "CD274", "n_Q4": 0, "n_Q1": 0, "median_Q4": np.nan, "median_Q1": np.nan,
         "U": np.nan, "p": np.nan, "rank_biserial_Q4_minus_Q1": np.nan, "note": "ABSENT on GPL96"},
    ])
    hl.to_csv(TABLES / "highlow_cldn4.tsv", sep="\t", index=False)

    paired = pd.DataFrame([
        paired_wilcoxon(cldn4_all, cldn4_all, luad_pairs, "CLDN4"),
        paired_wilcoxon(cd8a_all, cd8a_all, luad_pairs, "CD8A"),
        {"endpoint": "CD274", "n_pairs": 0, "n_tumor_gt_normal": 0, "median_T": np.nan,
         "median_N": np.nan, "median_T_minus_N": np.nan, "wilcoxon_stat": np.nan, "p": np.nan,
         "note": "ABSENT on GPL96"},
    ])
    # ImmuneScore is tumor-z; for paired T vs N use mean of the same present genes after
    # z-scoring within the 52 LUAD T+N arrays so the score is defined on both sides.
    tn_idx = pd.Index(list(luad_pairs["tumor_gsm"]) + list(luad_pairs["normal_gsm"]))
    imm_tn, _, _ = mean_z(gene_expr, IMMUNE8, tn_idx)
    paired_imm = paired_wilcoxon(imm_tn, imm_tn, luad_pairs, "ImmuneScore")
    paired = pd.concat([paired, pd.DataFrame([paired_imm])], ignore_index=True)
    paired.to_csv(TABLES / "paired_tumor_vs_normal.tsv", sep="\t", index=False)

    cov = pd.DataFrame([
        {"set": f"gene:{g}", "n_listed": 1, "n_present": int(g in gene_expr.index),
         "probe_used": probe_for(g),
         "missing": "" if g in gene_expr.index else g}
        for g in ["CLDN4", "CD8A", "CD274", "TACSTD2", "PDCD1LG2"] + EPITHELIAL + IMMUNE8
    ])
    cov.to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)

    # Sample annotation
    ann_out = meta[[
        "title", "source", "pair_id", "role", "histology", "sex",
        "primary_luad_tumor", "companion_luad_normal", "tissue_phrase",
    ]].copy()
    ann_out["CLDN4"] = cldn4_all.reindex(ann_out.index)
    ann_out["CD8A"] = cd8a_all.reindex(ann_out.index)
    ann_out["TACSTD2"] = gene_expr.loc["TACSTD2"].reindex(ann_out.index) if tacstd2_present else np.nan
    ann_out["epithelial_mean_z_LUAD_tumors"] = epi.reindex(ann_out.index)
    ann_out["ImmuneScore_A1_LUAD_tumors"] = imm.reindex(ann_out.index)
    ann_out.to_csv(TABLES / "sample_annotation.tsv", sep="\t")
    luad_pairs.to_csv(TABLES / "luad_pairs.tsv", sep="\t", index=False)

    # Figures
    def scatter(ax, x, y, title, xlab, ylab, color):
        ax.scatter(np.log2(np.asarray(x) + 1), np.log2(np.asarray(y) + 1),
                   s=28, alpha=0.8, c=color, edgecolors="none")
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(title)

    def pick(gene_x, gene_y, subset):
        hit = stats_df[(stats_df.gene_x == gene_x) & (stats_df.gene_y == gene_y) & (stats_df.subset == subset)]
        return hit.iloc[0].to_dict() if len(hit) else {}

    r8 = pick("CLDN4", "CD8A", "LUAD_tumors")
    ri = pick("CLDN4", "ImmuneScore", "LUAD_tumors")
    re = pick("CLDN4", "epithelial_mean_z", "LUAD_tumors")

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.4))
    scatter(
        axes[0], cldn4.values, cd8a.values,
        f"GSE7670 LUAD tumors  n={n_luad}\nCLDN4 vs CD8A  ρ={r8['rho']:.3f}  p={fmt_p(r8['p'])}",
        "CLDN4  log2(MAS5+1)", "CD8A  log2(MAS5+1)", "#2c5f8a",
    )
    axes[1].scatter(np.log2(cldn4.values + 1), imm.values, s=28, alpha=0.8, c="#3d6b4f", edgecolors="none")
    axes[1].set_xlabel("CLDN4  log2(MAS5+1)")
    axes[1].set_ylabel(f"ImmuneScore  A1 {len(imm_present)}/8 mean-z")
    axes[1].set_title(
        f"GSE7670 LUAD tumors  n={n_luad}\nCLDN4 vs ImmuneScore  ρ={ri['rho']:.3f}  p={fmt_p(ri['p'])}"
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_immunescore.png", dpi=160)
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_immunescore.pdf")
    plt.close(fig)

    forest = stats_df[
        (stats_df.gene_x == "CLDN4")
        & (stats_df.gene_y.isin(["CD8A", "ImmuneScore", "epithelial_mean_z"]))
        & (stats_df.subset == "LUAD_tumors")
    ].copy()
    forest["label"] = forest["gene_y"].map({
        "CD8A": "CLDN4 vs CD8A",
        "ImmuneScore": "CLDN4 vs ImmuneScore",
        "epithelial_mean_z": "CLDN4 vs epithelial mean-z",
    })
    order = ["CLDN4 vs CD8A", "CLDN4 vs ImmuneScore", "CLDN4 vs epithelial mean-z"]
    forest["label"] = pd.Categorical(forest["label"], categories=order, ordered=True)
    forest = forest.sort_values("label")
    fig, ax = plt.subplots(figsize=(6.8, 3.2))
    y_pos = np.arange(len(forest))
    ax.axvline(0, color="0.5", lw=0.8)
    ax.errorbar(
        forest["rho"], y_pos,
        xerr=[forest["rho"] - forest["ci_low"], forest["ci_high"] - forest["rho"]],
        fmt="o", color="#2c5f8a", ecolor="#2c5f8a", capsize=2, ms=6,
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(forest["label"])
    ax.set_xlabel("Spearman ρ (bootstrap 95% CI)")
    ax.set_title(f"GSE7670 LUAD tumors  n={n_luad}  (CD274 ABSENT on U133A)")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cldn4_forest.png", dpi=160)
    fig.savefig(FIGURES / "fig2_cldn4_forest.pdf")
    plt.close(fig)

    # Companion: paired T vs N for CLDN4 (not the claim)
    t_vals = cldn4_all.reindex(luad_pairs["tumor_gsm"]).values
    n_vals = cldn4_all.reindex(luad_pairs["normal_gsm"]).values
    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    ax.scatter(np.log2(n_vals + 1), np.log2(t_vals + 1), s=28, alpha=0.8, c="#6b3a2a", edgecolors="none")
    lims = [
        min(np.log2(n_vals + 1).min(), np.log2(t_vals + 1).min()) - 0.2,
        max(np.log2(n_vals + 1).max(), np.log2(t_vals + 1).max()) + 0.2,
    ]
    ax.plot(lims, lims, color="0.5", lw=0.8)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    pw = paired.loc[paired.endpoint == "CLDN4"].iloc[0]
    ax.set_xlabel("Adjacent normal CLDN4  log2(MAS5+1)")
    ax.set_ylabel("Matched tumor CLDN4  log2(MAS5+1)")
    ax.set_title(
        f"Companion only  {int(pw.n_pairs)} LUAD pairs\n"
        f"CLDN4 T vs N  Wilcoxon p={fmt_p(pw.p)}  T>N {int(pw.n_tumor_gt_normal)}/{int(pw.n_pairs)}"
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_companion_paired_cldn4.png", dpi=160)
    fig.savefig(FIGURES / "fig3_companion_paired_cldn4.pdf")
    plt.close(fig)

    # Residual scatter
    xr = stats.rankdata(cldn4.values)
    cr = np.column_stack([np.ones(n_luad), stats.rankdata(epi.values)])
    bx, *_ = np.linalg.lstsq(cr, xr, rcond=None)
    rx = xr - cr @ bx
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.4))
    for ax, yraw, name, rec, color in [
        (axes[0], cd8a.values, "CD8A", r8, "#2c5f8a"),
        (axes[1], imm.values, "ImmuneScore", ri, "#3d6b4f"),
    ]:
        yr = stats.rankdata(yraw)
        by, *_ = np.linalg.lstsq(cr, yr, rcond=None)
        ry = yr - cr @ by
        ax.scatter(rx, ry, s=28, alpha=0.8, c=color, edgecolors="none")
        ax.axhline(0, color="0.5", lw=0.6)
        ax.axvline(0, color="0.5", lw=0.6)
        ax.set_xlabel("CLDN4 rank residual | epithelium")
        ax.set_ylabel(f"{name} rank residual | epithelium")
        ax.set_title(
            f"n={n_luad}  ρ_adj={fmt_rho(rec['rho_adj_epithelial'])}  p={fmt_p(rec['p_adj_epithelial'])}"
        )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4_epithelial_residual.png", dpi=160)
    fig.savefig(FIGURES / "fig4_epithelial_residual.pdf")
    plt.close(fig)

    one_row = pd.DataFrame([{
        "dataset": "GSE7670 Su",
        "platform": "GPL96 U133A MAS5",
        "histology": "LUAD paired (pair 19 LCC excluded)",
        "n_arrays_in_series": n_array,
        "n_LUAD_tumors": n_luad,
        "n_LUAD_pairs": int(len(luad_pairs)),
        "n_do_not_use": "66 series / 27 GEO-text patients / 1 LCC pair / 8 lines / 2 mixtures",
        "CLDN4_probe": probe_for("CLDN4") or NAMED_PROBES["CLDN4"],
        "CD8A_probe": probe_for("CD8A") or NAMED_PROBES["CD8A"],
        "CD274": "ABSENT (Entrez 29126 not on U133A)",
        "ImmuneScore": f"A1 {len(imm_present)}/8 mean-z (missing {','.join(imm_missing)}; not ESTIMATE)",
        "epithelial": f"{len(epi_present)}/6 mean-z",
        "CLDN4_CD8A_n": r8["n"],
        "CLDN4_CD8A_rho": r8["rho"],
        "CLDN4_CD8A_p": r8["p"],
        "CLDN4_CD8A_rho_adj": r8["rho_adj_epithelial"],
        "CLDN4_CD8A_p_adj": r8["p_adj_epithelial"],
        "CLDN4_CD8A_verdict": r8["verdict"],
        "CLDN4_ImmuneScore_n": ri["n"],
        "CLDN4_ImmuneScore_rho": ri["rho"],
        "CLDN4_ImmuneScore_p": ri["p"],
        "CLDN4_ImmuneScore_rho_adj": ri["rho_adj_epithelial"],
        "CLDN4_ImmuneScore_p_adj": ri["p_adj_epithelial"],
        "CLDN4_ImmuneScore_verdict": ri["verdict"],
        "CLDN4_CD274_n": 0,
        "CLDN4_CD274_verdict": "ABSENT",
        "CLDN4_epithelial_rho": re["rho"],
        "CLDN4_epithelial_p": re["p"],
        "ICI": "not deposited",
        "dual_high": "not done",
    }])
    one_row.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    summary = {
        "dataset": "GSE7670",
        "pmid": "17540040",
        "title": series.get("title", ""),
        "platform": "GPL96 Affymetrix HG-U133A",
        "processing": "MAS5 as deposited (target 500). Spearman is rank-based.",
        "n_arrays": n_array,
        "n_probes": int(probe_expr.shape[0]),
        "n_missing_mas5": n_missing,
        "n_LUAD_tumors_primary": n_luad,
        "n_LUAD_pairs": int(len(luad_pairs)),
        "n_LCC_tumors_excluded": n_lcc_tumor,
        "n_cell_lines_excluded": n_cell,
        "n_mixtures_excluded": n_mix,
        "cldn4_probe": probe_for("CLDN4"),
        "cd8a_probe": probe_for("CD8A"),
        "cd274": "ABSENT",
        "immune_genes_used": imm_present,
        "immune_genes_missing": imm_missing,
        "epithelial_genes_used": epi_present,
        "primary": {
            "CLDN4_vs_CD8A_LUAD": r8,
            "CLDN4_vs_ImmuneScore_LUAD": ri,
            "CLDN4_vs_CD274_LUAD": pick("CLDN4", "CD274", "LUAD_tumors"),
            "CLDN4_vs_epithelial_LUAD": re,
        },
        "paired_companion": paired.to_dict(orient="records"),
        "highlow": hl.to_dict(orient="records"),
        "missing_public_labels": ["stage", "ICI", "tumor_percent", "CD274", "EOMES", "ESTIMATE"],
        "holds_rule": "n>=40, rho_adj<0, p_adj<0.05 (PR 229 / GSE4573). n=26 is UNDERPOWERED.",
        "dual_high": "not done",
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    print("n_arrays", n_array, "LUAD tumors", n_luad, "LUAD pairs", len(luad_pairs), "missing", n_missing)
    print("ImmuneScore genes", imm_present, "missing", imm_missing)
    print(stats_df[stats_df.subset == "LUAD_tumors"][
        ["gene_x", "gene_y", "n", "rho", "p", "rho_adj_epithelial", "p_adj_epithelial", "verdict"]
    ].to_string(index=False))
    print(paired.to_string(index=False))
    print("[done]", TABLES)


if __name__ == "__main__":
    main()
