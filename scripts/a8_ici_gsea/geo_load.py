"""Load leftover public ICI lung bulk matrices and align GEO labels."""
from __future__ import annotations

import gzip
import os
import re

import numpy as np
import pandas as pd

TARGET = "TACSTD2"


def parse_series_matrix(path: str) -> pd.DataFrame:
    titles = geo = descriptions = None
    char_rows = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                _, titles = _split(line)
            elif line.startswith("!Sample_geo_accession"):
                _, geo = _split(line)
            elif line.startswith("!Sample_description"):
                _, descriptions = _split(line)
            elif line.startswith("!Sample_characteristics_ch1"):
                _, vals = _split(line)
                char_rows.append(vals)
            elif line.startswith("!series_matrix_table_begin"):
                break
    n = len(titles) if titles else 0
    rows = []
    for i in range(n):
        d = {
            "title": titles[i] if titles else "",
            "geo_accession": geo[i] if geo else "",
            "description": descriptions[i] if descriptions and i < len(descriptions) else "",
        }
        for row in char_rows:
            if i < len(row) and row[i] and ":" in row[i]:
                k, v = row[i].split(":", 1)
                d[k.strip().lower()] = v.strip()
        rows.append(d)
    return pd.DataFrame(rows)


def _split(line: str):
    parts = line.rstrip("\n").split("\t")
    return parts[0].strip().strip('"'), [p.strip().strip('"') for p in parts[1:]]


def collapse_symbols(expr: pd.DataFrame) -> pd.DataFrame:
    expr = expr.copy()
    expr.index = expr.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    expr.index = expr.index.str.strip().str.strip('"')
    expr = expr.apply(pd.to_numeric, errors="coerce")
    if expr.index.duplicated().any():
        means = expr.mean(axis=1)
        keep_idx = means.groupby(level=0).idxmax()
        expr = expr.loc[keep_idx]
    return expr.loc[~expr.index.duplicated(keep="first")]


def log2cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts.div(lib, axis=1) * 1e6 + 1.0)


def log2p1(x: pd.DataFrame) -> pd.DataFrame:
    return np.log2(x.clip(lower=0) + 1.0)


def looks_logged(x: pd.DataFrame) -> bool:
    arr = x.to_numpy(dtype=float)
    return float(np.nanmax(arr)) < 25 and float(np.nanmedian(arr)) < 12


def load_ensembl_map(data_dir: str) -> dict[str, str]:
    path = os.path.join(data_dir, "hgnc_complete_set.tsv")
    if not os.path.exists(path):
        return {}
    h = pd.read_csv(path, sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str)
    h = h.dropna()
    return dict(zip(h["ensembl_gene_id"], h["symbol"]))


def map_ensembl(expr: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    idx = expr.index.astype(str)
    if idx.str.startswith("ENSG").mean() < 0.3:
        return expr
    new = [mapping.get(i, i) for i in idx]
    expr = expr.copy()
    expr.index = new
    return collapse_symbols(expr)


def _norm_key(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "", str(s)).upper()


def align_meta(expr: pd.DataFrame, meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rename expression columns to geo_accession when a unique match exists."""
    meta = meta.copy()
    meta["geo_accession"] = meta["geo_accession"].astype(str)
    candidates = {}
    for _, row in meta.iterrows():
        g = row["geo_accession"]
        for raw in (row.get("title", ""), row.get("description", ""), g):
            candidates.setdefault(_norm_key(raw), g)
            # also last token (Dis_01 from RNA-seq_Dis_01; NSCLC990 from NSCLC 990)
            tok = re.sub(r"[\s_]+", "", str(raw).split()[-1] if str(raw).split() else "")
            if tok:
                candidates.setdefault(_norm_key(tok), g)
            # strip common prefixes
            for pref in ("RNASEQ", "S"):
                k = _norm_key(raw)
                if k.startswith(pref) and len(k) > len(pref) + 2:
                    candidates.setdefault(k[len(pref) :], g)

    rename = {}
    for c in expr.columns:
        k = _norm_key(c)
        if c in set(meta["geo_accession"]):
            rename[c] = c
        elif k in candidates:
            rename[c] = candidates[k]
        else:
            # strip leading S
            if k.startswith("S") and k[1:] in candidates:
                rename[c] = candidates[k[1:]]
    if len(rename) >= max(8, expr.shape[1] // 2):
        expr = expr.rename(columns=rename)
    common = [c for c in expr.columns if c in set(meta["geo_accession"])]
    if not common and expr.shape[1] == len(meta):
        expr = expr.copy()
        expr.columns = list(meta["geo_accession"])
        common = list(expr.columns)
    expr = expr.loc[:, common]
    meta = meta.drop_duplicates("geo_accession").set_index("geo_accession", drop=False).loc[common]
    return expr, meta


def norm_response(v):
    s = str(v).strip().lower()
    if s in {"", "nan", "none", "na"}:
        return np.nan
    if s in {"r", "responder", "response", "pr", "cr", "mpr", "dcb", "yes", "y", "benefit"}:
        return "R"
    if s in {"nr", "non-responder", "nonresponder", "non-response", "pd", "nmpr", "ndb", "no", "n"}:
        return "NR"
    if "nonresponder" in s.replace("-", "").replace(" ", "") or ("non" in s and "respond" in s):
        return "NR"
    if "respond" in s:
        return "R"
    if s.startswith("complete") or s == "cr":
        return "R"
    if "partial" in s or s == "pr":
        return "R"
    if "progress" in s or s == "pd":
        return "NR"
    return np.nan


def load_gse126044(data: str, mapping: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(os.path.join(data, "GSE126044_counts.txt.gz"), sep="\t", index_col=0)
    counts = collapse_symbols(counts)
    expr = log2cpm(counts)
    meta = parse_series_matrix(os.path.join(data, "GSE126044_series_matrix.txt.gz"))
    expr, meta = align_meta(expr, meta)
    raw = meta.get("patient response", pd.Series("", index=meta.index))
    meta["response_raw"] = raw.astype(str)
    meta["response"] = meta["response_raw"].map(norm_response)
    return expr, meta


def load_gse135222(data: str, mapping: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(os.path.join(data, "GSE135222_exp.tsv.gz"), sep="\t", index_col=0)
    expr = map_ensembl(collapse_symbols(expr), mapping)
    if not looks_logged(expr):
        expr = log2p1(expr)
    meta = parse_series_matrix(os.path.join(data, "GSE135222_series_matrix.txt.gz"))
    expr, meta = align_meta(expr, meta)
    pfs = pd.to_numeric(meta.get("pfs.time", np.nan), errors="coerce")
    meta["pfs_days"] = pfs
    meta["response"] = ["R" if pd.notna(v) and v >= 180 else ("NR" if pd.notna(v) else np.nan) for v in pfs]
    meta["response_raw"] = [f"PFS_days={v}" if pd.notna(v) else "" for v in pfs]
    return expr, meta


def load_gse166449(data: str, mapping: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(os.path.join(data, "GSE166449_TPM.txt.gz"), sep="\t", index_col=0)
    expr = collapse_symbols(expr)
    if not looks_logged(expr):
        expr = log2p1(expr)
    meta = parse_series_matrix(os.path.join(data, "GSE166449_series_matrix.txt.gz"))
    expr, meta = align_meta(expr, meta)
    # titles are Immunotherapy_ResponderN / Immunotherapy_nonResponderN
    meta["response_raw"] = meta["title"].astype(str)
    meta["response"] = meta["title"].map(norm_response)
    return expr, meta


def load_gse253564(data: str, mapping: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(os.path.join(data, "GSE253564_pre_FPKM.txt.gz"), sep="\t")
    # gene, Entrez.ID, samples...
    gene_col = raw.columns[0]
    drop = [c for c in raw.columns if c.lower() in {"entrez.id", "entrez", "entrez_id"}]
    expr = raw.set_index(gene_col).drop(columns=drop, errors="ignore")
    expr = collapse_symbols(expr)
    expr = log2p1(expr)
    meta = parse_series_matrix(os.path.join(data, "GSE253564_series_matrix.txt.gz"))
    expr, meta = align_meta(expr, meta)
    meta["response_raw"] = ""
    meta["response"] = np.nan
    return expr, meta


def load_gse190265(data: str, mapping: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = os.path.join(data, "GSE190265_TPM_France3.csv.gz")
    with gzip.open(path, "rt") as fh:
        genes = fh.readline().rstrip("\n").split(";")
    expr = pd.read_csv(path, sep=";", names=["sample"] + genes, skiprows=1)
    expr = expr.set_index("sample").T
    expr = collapse_symbols(expr)
    if not looks_logged(expr):
        expr = log2p1(expr)
    info = pd.read_csv(os.path.join(data, "GSE190265_samples_info_France3.csv.gz"), sep=";")
    info["sample"] = info["sample"].astype(str)
    meta_sm = parse_series_matrix(os.path.join(data, "GSE190265_series_matrix.txt.gz"))
    # expression columns are sample IDs from the TPM file
    keep = [c for c in expr.columns if str(c) in set(info["sample"])]
    expr = expr.loc[:, keep]
    info = info.drop_duplicates("sample").set_index("sample").loc[keep]
    meta = info.copy()
    meta["geo_accession"] = keep
    meta["title"] = keep
    pfs = pd.to_numeric(meta["time_PFS"], errors="coerce")
    meta["pfs_months"] = pfs
    meta["response"] = ["R" if pd.notna(v) and v >= 6 else ("NR" if pd.notna(v) else np.nan) for v in pfs]
    meta["response_raw"] = [f"PFS_months={v}" if pd.notna(v) else "" for v in pfs]
    # attach histology from series matrix if titles match
    if not meta_sm.empty:
        tmap = {_norm_key(t): g for t, g in zip(meta_sm["title"], meta_sm["geo_accession"])}
        meta["series_gsm"] = [tmap.get(_norm_key(s), "") for s in keep]
    return expr, meta


def load_gse283829(data: str, mapping: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(os.path.join(data, "GSE283829_raw_counts.txt.gz"), sep="\t", index_col=0)
    raw.columns = [str(c).strip().strip('"') for c in raw.columns]
    raw = map_ensembl(collapse_symbols(raw), mapping)
    expr = log2cpm(raw)
    meta = parse_series_matrix(os.path.join(data, "GSE283829_series_matrix.txt.gz"))
    expr, meta = align_meta(expr, meta)
    stage = meta.get("disease stage", pd.Series("", index=meta.index)).astype(str)
    meta["response_raw"] = stage
    mapped = []
    for v in stage:
        s = str(v).strip().lower()
        if s in {"cr", "complete response"}:
            mapped.append("R")
        elif s in {"pd", "progressive disease"}:
            mapped.append("NR")
        else:
            mapped.append(np.nan)  # SD held out of R vs NR
    meta["response"] = mapped
    return expr, meta


LOADERS = {
    "GSE126044": load_gse126044,
    "GSE135222": load_gse135222,
    "GSE166449": load_gse166449,
    "GSE253564": load_gse253564,
    "GSE190265": load_gse190265,
    "GSE283829": load_gse283829,
}
