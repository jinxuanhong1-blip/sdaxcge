#!/usr/bin/env python3
"""Public mouse lung KP/KL matrices beyond GSE137244.

Score Cldn4, an NHEJ module, and a cGAS–STING module. Each accession
is tested on its own. Private 8 KL mice are not read and are not merged
with any public matrix. GSE137244 is the locked cell-line baseline and
is not re-opened.
"""

from __future__ import annotations

import gzip
import json
import os
import re
from io import StringIO
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from gene_sets import MODULES, all_aliases

ROOT = Path(__file__).resolve().parents[2]
DATA = Path(os.environ.get("KPKL_DATA", "/tmp/kpkl"))
OUT = ROOT / "results" / "mouse_kpkl_nhej_sting"
TAB = OUT / "tables"
FIG = OUT / "figures"

WANTED = all_aliases()


def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    import urllib.request

    req = urllib.request.Request(url, headers={"User-Agent": "cursor-agent"})
    with urllib.request.urlopen(req, timeout=180) as fh, dest.open("wb") as out:
        out.write(fh.read())
    return dest


def bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return []
    order = np.argsort(np.where(np.isnan(p), np.inf, p))
    ranked = p[order]
    q = np.full(n, np.nan)
    valid = ~np.isnan(ranked)
    if valid.any():
        ranks = np.arange(1, n + 1, dtype=float)
        raw = ranked * n / ranks
        # NaNs stay at the end of argsort; fix only the finite prefix.
        finite = np.where(valid)[0]
        raw_f = raw[finite].copy()
        raw_f = np.minimum.accumulate(raw_f[::-1])[::-1]
        raw_f = np.clip(raw_f, 0, 1)
        q[finite] = raw_f
    out = np.full(n, np.nan)
    out[order] = q
    return [float(x) if not np.isnan(x) else float("nan") for x in out]


def hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    """Hedges' g for b minus a (KL minus comparator)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return float("nan")
    va, vb = a.var(ddof=1), b.var(ddof=1)
    sp2 = ((na - 1) * va + (nb - 1) * vb) / (na + nb - 2)
    if sp2 <= 0:
        return float("nan")
    g = (b.mean() - a.mean()) / np.sqrt(sp2)
    j = 1 - 3 / (4 * (na + nb) - 9)
    return float(j * g)


def mwu(a: np.ndarray, b: np.ndarray) -> float:
    """Two-sided Mann–Whitney, comparator (a) vs KL (b)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 1 or len(b) < 1:
        return float("nan")
    if np.all(a == a[0]) and np.all(b == b[0]) and a[0] == b[0]:
        return 1.0
    res = stats.mannwhitneyu(b, a, alternative="two-sided", method="auto")
    return float(res.pvalue)


def spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = ~(np.isnan(x) | np.isnan(y))
    if m.sum() < 3 or np.unique(x[m]).size < 2 or np.unique(y[m]).size < 2:
        return float("nan"), float("nan")
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p)


def log2p1(df: pd.DataFrame) -> pd.DataFrame:
    return np.log2(df.clip(lower=0) + 1.0)


def cpm_log(df: pd.DataFrame) -> pd.DataFrame:
    lib = df.sum(axis=0).replace(0, np.nan)
    cpm = df.divide(lib, axis=1) * 1e6
    return log2p1(cpm)


def collapse_symbols(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.index = df.index.astype(str)
    df = df.groupby(level=0).mean()
    return df


def resolve_genes(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return a canonical-symbol matrix and a presence table."""
    index = {str(i): i for i in df.index}
    lower = {str(i).lower(): str(i) for i in df.index}
    rows = []
    found = {}
    for module, genes in MODULES.items():
        for canon, aliases in genes.items():
            hit = None
            for alias in aliases:
                if alias in index:
                    hit = alias
                    break
                if alias.lower() in lower:
                    hit = lower[alias.lower()]
                    break
            rows.append(
                {
                    "module": module,
                    "canonical": canon,
                    "matched_symbol": hit or "",
                    "present": "yes" if hit else "no",
                }
            )
            if hit:
                found[canon] = df.loc[hit]
    mat = pd.DataFrame(found).T if found else pd.DataFrame(index=[], columns=df.columns)
    if len(mat):
        # If a symbol collision produced a DataFrame, average.
        mat = mat.groupby(level=0).mean()
    return mat, pd.DataFrame(rows)


def zmean(mat: pd.DataFrame, genes: list[str], samples: list[str]) -> pd.Series:
    use = [g for g in genes if g in mat.index]
    if not use:
        return pd.Series(np.nan, index=samples)
    sub = mat.loc[use, samples].astype(float)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=0).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def module_gene_list(name: str) -> list[str]:
    return list(MODULES[name])


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def read_table(path: Path, sep: str = "\t") -> pd.DataFrame:
    return pd.read_csv(path, sep=sep, compression="gzip")


def load_symbol_matrix_simple(path: Path, gene_col: str, sample_cols: list[str] | None = None) -> pd.DataFrame:
    df = read_table(path)
    if sample_cols is None:
        sample_cols = [c for c in df.columns if c != gene_col]
    out = df.set_index(gene_col)[sample_cols].apply(pd.to_numeric, errors="coerce")
    return collapse_symbols(out)


def load_gse137396() -> tuple[pd.DataFrame, str]:
    path = DATA / "GSE137396_Normalized_logtransformed_medcentred_genetable_GEMMnodule.txt.gz"
    df = read_table(path)
    df = df.set_index("Genes")
    df = df.apply(pd.to_numeric, errors="coerce")
    return collapse_symbols(df), "depositor_log2_median_centered"


def load_gse274351(ensembl_map: dict[str, str]) -> tuple[pd.DataFrame, str]:
    path = DATA / "GSE274351_expredata_TPM_gene.txt.gz"
    df = read_table(path)
    df["gene_id"] = df["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    df["symbol"] = df["gene_id"].map(ensembl_map)
    df = df.dropna(subset=["symbol"])
    sample_cols = [c for c in df.columns if c not in ("gene_id", "symbol")]
    out = df.groupby("symbol")[sample_cols].mean()
    return log2p1(out), "log2_TPM_plus1"


def _gse274352_wide(path: Path) -> pd.DataFrame:
    df = read_table(path)
    gene = "external_gene_name"
    if gene not in df.columns:
        raise SystemExit(f"GSE274352 missing gene name column: {df.columns[:5].tolist()}")
    cols = [
        c
        for c in df.columns
        if c not in ("", gene)
        and not str(c).startswith("ENSMUS")
        and not str(c).startswith("Unnamed")
    ]
    out = df.set_index(gene)[cols].apply(pd.to_numeric, errors="coerce")
    return collapse_symbols(out)


def load_gse274352_ifnb() -> tuple[pd.DataFrame, str]:
    """Within-file normalized counts. Do not paste onto the STING file."""
    return log2p1(_gse274352_wide(DATA / "GSE274352_normalizedcounts_genes_IFNB_vs_empty.tsv.gz")), "log2_normalized_counts_plus1"


def load_gse274352_sting() -> tuple[pd.DataFrame, str]:
    return log2p1(_gse274352_wide(DATA / "GSE274352_normalizedcounts_genes_STING_vs_emtpy.tsv.gz")), "log2_normalized_counts_plus1"


def load_gse175479() -> tuple[pd.DataFrame, str]:
    path = DATA / "GSE175479_Raw_gene_count_matrix.txt.gz"
    # Header omits a leading row-index column. Data are: index, symbol, 12 counts.
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        header = [h.strip('"') for h in header]
        samples = header[1:]
        rows = []
        symbols = []
        for line in fh:
            parts = [p.strip('"') for p in line.rstrip("\n").split("\t")]
            if len(parts) < 3:
                continue
            symbols.append(parts[1])
            rows.append([float(x) if x not in ("", "NA") else np.nan for x in parts[2:]])
    df = pd.DataFrame(rows, index=symbols, columns=samples[: len(rows[0])])
    df = collapse_symbols(df)
    return cpm_log(df), "log2_CPM_plus1"


def load_gse164758() -> tuple[pd.DataFrame, str]:
    path = DATA / "GSE164758_primary_tumors_fpkm.txt.gz"
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # Column names sit after the annotation fields. The cmd= string in
        # field 0 also contains "Aligned.out.sam" and must not be parsed.
        sams = []
        for field in header:
            m = re.match(r"^(\S+)Aligned\.out\.sam FPKM$", field.strip())
            if m:
                sams.append(m.group(1))
        if len(sams) != 41:
            raise SystemExit(f"GSE164758: expected 41 FPKM columns, got {len(sams)}: {sams[:5]}")
        records = []
        symbols = []
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            ann_i = next((i for i, p in enumerate(parts) if "|" in p), None)
            if ann_i is None:
                continue
            vals = parts[ann_i + 1 :]
            if len(vals) < len(sams):
                continue
            vals = vals[: len(sams)]
            try:
                nums = [float(x) for x in vals]
            except ValueError:
                continue
            symbols.append(parts[ann_i].split("|")[0])
            records.append(nums)
    df = pd.DataFrame(records, index=symbols, columns=sams)
    df = collapse_symbols(df)
    if df.shape[1] != len(sams):
        raise SystemExit("GSE164758 column mismatch")
    return log2p1(df), "log2_FPKM_plus1"


def load_gse193895(ensembl_map: dict[str, str]) -> tuple[pd.DataFrame, str]:
    path = DATA / "GSE193895_Tumor_KrasModel_GEMMs_RNAseq_GenewiseCounts.txt.gz"
    df = read_table(path)
    df["gene_id"] = df["GeneID"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    df["symbol"] = df["gene_id"].map(ensembl_map)
    df = df.dropna(subset=["symbol"])
    sample_cols = [c for c in df.columns if c not in ("GeneID", "Length", "gene_id", "symbol")]
    out = df.groupby("symbol")[sample_cols].mean()
    return cpm_log(out), "log2_CPM_plus1"


def load_gse338923() -> tuple[pd.DataFrame, str]:
    path = DATA / "GSE338923_Lacun3_STK11_RNAseq_raw_counts.tsv.gz"
    df = read_table(path)
    sample_cols = [c for c in df.columns if re.match(r"^[CS]\d", c)]
    out = df.groupby("gene_name")[sample_cols].mean()
    return cpm_log(out.apply(pd.to_numeric, errors="coerce")), "log2_CPM_plus1"


def load_gse322570() -> tuple[pd.DataFrame, str]:
    path = DATA / "GSE322570_raw_count.txt.gz"
    df = read_table(path)
    sample_cols = [c for c in df.columns if c.startswith("sg")]
    out = df.groupby("gene_name")[sample_cols].mean()
    return cpm_log(out.apply(pd.to_numeric, errors="coerce")), "log2_CPM_plus1"


def load_gse133895() -> tuple[pd.DataFrame, str]:
    path = DATA / "GSE133895_Sik_exvivo_Count_Matrix.txt.gz"
    # Title line, then a header of sample names only. Data rows are symbol + counts.
    with gzip.open(path, "rt") as fh:
        _title = fh.readline()
        header = fh.readline().rstrip("\n").split("\t")
        rows = []
        symbols = []
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            symbols.append(parts[0])
            rows.append([float(x) if x not in ("", "NA") else np.nan for x in parts[1:]])
    df = pd.DataFrame(rows, index=symbols, columns=header[: len(rows[0])])
    return cpm_log(collapse_symbols(df)), "log2_CPM_plus1"


def load_gse244452() -> tuple[pd.DataFrame, str]:
    path = DATA / "GSE244452_KPvsKL_deg_all.txt.gz"
    df = read_table(path)
    cols = ["KP1", "KP2", "KP3", "KL1", "KL2", "KL3"]
    out = df.groupby("gene_name")[cols].mean().apply(pd.to_numeric, errors="coerce")
    return log2p1(out), "log2_depositor_sample_values_plus1"


def _series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    titles = geos = chars = None
    table: list[str] = []
    in_table = False
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                geos = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") and chars is None:
                chars = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            elif line.startswith("!series_matrix_table_end"):
                break
            elif in_table:
                table.append(line)
    expr = pd.read_csv(StringIO("".join(table)), sep="\t")
    expr = expr.set_index("ID_REF")
    expr.columns = [c.strip('"') for c in expr.columns]
    expr = expr.apply(pd.to_numeric, errors="coerce")
    meta = pd.DataFrame({"gsm": geos, "title": titles, "characteristics": chars})
    return expr, meta


def _annot_map(path: Path) -> dict[str, list[str]]:
    """probe -> list of gene symbols."""
    out: dict[str, list[str]] = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        go = False
        header = None
        sym_i = id_i = None
        for line in fh:
            if line.startswith("!platform_table_begin"):
                go = True
                header = next(fh).rstrip("\n").split("\t")
                sym_i = header.index("Gene symbol")
                id_i = header.index("ID")
                continue
            if line.startswith("!platform_table_end"):
                break
            if not go:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(sym_i, id_i):
                continue
            raw = parts[sym_i].strip()
            if not raw or raw == "---":
                continue
            syms = [s.strip() for s in raw.split("///") if s.strip()]
            if any(s in WANTED or s.lower() in {w.lower() for w in WANTED} for s in syms):
                out[parts[id_i]] = syms
    return out


def probes_to_genes(expr: pd.DataFrame, amap: dict[str, list[str]]) -> pd.DataFrame:
    buckets: dict[str, list[str]] = {}
    wanted_lower = {w.lower(): w for w in WANTED}
    for probe, syms in amap.items():
        if probe not in expr.index:
            continue
        for s in syms:
            key = s if s in WANTED else wanted_lower.get(s.lower())
            if key:
                buckets.setdefault(key, []).append(probe)
    rows = {}
    for sym, probes in buckets.items():
        rows[sym] = expr.loc[probes].mean(axis=0)
    return pd.DataFrame(rows).T


def maybe_log_array(df: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    mx = float(np.nanmax(df.values))
    mn = float(np.nanmin(df.values))
    if mx > 40 and mn >= 0:
        return log2p1(df), "log2_array_plus1"
    return df, "array_already_log_or_centered"


def build_ensembl_map() -> dict[str, str]:
    """Ensembl gene id (no version) -> symbol, from matrices that carry both."""
    mapping: dict[str, str] = {}
    sting = read_table(DATA / "GSE274352_normalizedcounts_genes_STING_vs_emtpy.tsv.gz")
    idcol = sting.columns[0]
    def _add(gids, syms, overwrite: bool) -> None:
        for gid, sym in zip(gids, syms):
            if not isinstance(sym, str) or not sym or sym == "nan":
                continue
            key = re.sub(r"\.\d+$", "", str(gid))
            if overwrite or key not in mapping:
                mapping[key] = sym

    _add(sting[idcol].tolist(), sting["external_gene_name"].tolist(), True)
    for path, idc, symc in [
        (DATA / "GSE322570_raw_count.txt.gz", "gene_id", "gene_name"),
        (DATA / "GSE338923_Lacun3_STK11_RNAseq_raw_counts.tsv.gz", "gene_id", "gene_name"),
    ]:
        df = read_table(path)
        _add(df[idc].tolist(), df[symc].tolist(), False)
    # Fill any canonical alias still missing via Ensembl REST.
    have = {s for s in mapping.values() if isinstance(s, str)}
    have_lower = {s.lower() for s in have}
    missing = [a for a in sorted(WANTED) if a not in have and a.lower() not in have_lower]
    if missing:
        import json as _json
        import urllib.request

        url = "https://rest.ensembl.org/lookup/symbol/mus_musculus"
        req = urllib.request.Request(
            url,
            data=_json.dumps({"symbols": missing}).encode(),
            headers={"Content-Type": "application/json", "Accept": "application/json"},
        )
        try:
            with urllib.request.urlopen(req, timeout=60) as fh:
                res = _json.loads(fh.read().decode())
            for sym, rec in res.items():
                if rec and rec.get("id"):
                    mapping[rec["id"]] = sym
        except Exception as exc:
            print("Ensembl REST fallback failed:", exc)
    return mapping


# ---------------------------------------------------------------------------
# Contrasts
# ---------------------------------------------------------------------------

def genotype_from_prefix(name: str) -> str:
    n = name
    if n.startswith("KPL") or n.startswith("KPL"):
        return "KPL"
    if n.startswith("KP"):
        return "KP"
    if n.startswith("KL"):
        return "KL"
    if n.startswith("Kras") or n.startswith("K-") or n == "K":
        return "K"
    return "other"


def score_contrast(
    accession: str,
    design: str,
    family: str,
    mat: pd.DataFrame,
    scale: str,
    groups: dict[str, list[str]],
    kl_name: str,
    ref_name: str,
    unit: str,
    note: str,
) -> tuple[list[dict], pd.DataFrame, list[dict]]:
    """groups maps a label to sample ids already at the analysis unit."""
    presence_rows = []
    gene_mat, presence = resolve_genes(mat)
    presence.insert(0, "accession", accession)
    presence_rows = presence.to_dict(orient="records")

    kl_samples = groups[kl_name]
    ref_samples = groups[ref_name]
    samples = ref_samples + kl_samples
    # Drop samples missing from the matrix.
    samples = [s for s in samples if s in gene_mat.columns]
    kl_samples = [s for s in kl_samples if s in gene_mat.columns]
    ref_samples = [s for s in ref_samples if s in gene_mat.columns]

    endpoints = {
        "Cldn4": None,
        "Sting1": None,
        "Cgas": None,
        "Stk11": None,
        "NHEJ": module_gene_list("NHEJ"),
        "STING_core": module_gene_list("STING_core"),
        "STING_ISG": module_gene_list("STING_ISG"),
    }
    values = {}
    for ep, genes in endpoints.items():
        if genes is None:
            if ep in gene_mat.index:
                values[ep] = gene_mat.loc[ep, samples].astype(float)
            else:
                values[ep] = pd.Series(np.nan, index=samples)
        else:
            values[ep] = zmean(gene_mat, genes, samples)

    present_counts = {}
    for mod in ("NHEJ", "STING_core", "STING_ISG"):
        genes = module_gene_list(mod)
        present_counts[mod] = sum(1 for g in genes if g in gene_mat.index)

    extra_genes = [
        g
        for g in gene_mat.index
        if g not in ("Cldn4", "Sting1", "Cgas", "Stk11", "NHEJ", "STING_core", "STING_ISG")
    ]
    for g in extra_genes:
        values[g] = gene_mat.loc[g, samples].astype(float)

    rows = []
    for ep in list(endpoints) + extra_genes:
        vec = values[ep]
        a = vec.reindex(ref_samples).to_numpy(dtype=float)
        b = vec.reindex(kl_samples).to_numpy(dtype=float)
        p = mwu(a, b)
        if ep in present_counts:
            n_in, n_req = present_counts[ep], len(module_gene_list(ep))
        elif ep in gene_mat.index:
            n_in, n_req = 1, 1
        else:
            n_in, n_req = 0, 1
        n_a = int(np.sum(~np.isnan(a)))
        n_b = int(np.sum(~np.isnan(b)))
        rows.append(
            {
                "accession": accession,
                "design": design,
                "family": family,
                "contrast": f"{kl_name}_minus_{ref_name}",
                "endpoint": ep,
                "scale": scale,
                "unit": unit,
                "n_kl": n_b,
                "n_ref": n_a,
                "mean_kl": float(np.nanmean(b)) if n_b else float("nan"),
                "mean_ref": float(np.nanmean(a)) if n_a else float("nan"),
                "delta_kl_minus_ref": float(np.nanmean(b) - np.nanmean(a)) if n_a and n_b else float("nan"),
                "hedges_g": hedges_g(a, b),
                "mannwhitney_p": p,
                "n_genes_in_module": n_in,
                "n_genes_requested": n_req,
                "underpowered_mw_floor": "yes" if min(n_a, n_b) < 4 else "no",
                "note": note,
            }
        )

    score = pd.DataFrame({ep: values[ep] for ep in values})
    score.insert(0, "sample", score.index)
    score.insert(0, "accession", accession)
    score.insert(2, "group", [kl_name if s in kl_samples else ref_name for s in score["sample"]])
    score["design"] = design

    spearman_rows = []
    for other in ("NHEJ", "STING_core", "Sting1", "STING_ISG"):
        r, p = spearman(values["Cldn4"].to_numpy(dtype=float), values[other].to_numpy(dtype=float))
        spearman_rows.append(
            {
                "accession": accession,
                "design": design,
                "family": family,
                "x": "Cldn4",
                "y": other,
                "n": int(np.sum(~(values["Cldn4"].isna() | values[other].isna()))),
                "spearman_r": r,
                "spearman_p": p,
                "note": "within contrast samples only; not a cross-dataset merge",
            }
        )
    return rows, score, spearman_rows, presence_rows


def mouse_mean(df_scores_later: pd.DataFrame) -> pd.DataFrame:
    return df_scores_later


def gse6135_groups(meta: pd.DataFrame) -> dict[str, dict[str, list[str]]]:
    """Return sample-level and mouse-level group maps. Unit id is GSM or mouse."""

    def klass(title: str) -> str:
        t = title.lower()
        if "a549" in t or "h2126" in t or "oxygen" in t:
            return "human_line"
        if re.search(r"\bmet\b", t):
            return "metastasis"
        if "primary" not in t:
            return "other"
        if "lkb" in t or "stk11" in t:
            if "l/+" in t or "l +" in t:
                return "KL_het"
            if "l/l" in t or "l/-" in t or "l/l" in t:
                return "KL"
            return "KL_other"
        if "p53" in t:
            return "KP"
        if "p16" in t or "ink4" in t:
            return "K_Ink4a"
        if "kras" in t or "k-ras" in t:
            return "K"
        return "other"

    def mouse_id(title: str) -> str:
        head = title.split("||")[0].strip().replace(" ", "")
        m = re.match(r"(\d+)", head)
        return m.group(1) if m else head

    meta = meta.copy()
    meta["genotype"] = [klass(t) for t in meta["title"]]
    meta["mouse"] = [mouse_id(t) for t in meta["title"]]
    # Mouse-level ids are genotype + mouse number so a shared numeric id
    # across genotypes is not collapsed.
    return meta


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    urls = {
        "GSE137396_Normalized_logtransformed_medcentred_genetable_GEMMnodule.txt.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137396/suppl/GSE137396_Normalized_logtransformed_medcentred_genetable_GEMMnodule.txt.gz",
        "GSE274351_expredata_TPM_gene.txt.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274351/suppl/GSE274351_expredata_TPM_gene.txt.gz",
        "GSE274352_normalizedcounts_genes_STING_vs_emtpy.tsv.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274352/suppl/GSE274352_normalizedcounts_genes_STING_vs_emtpy.tsv.gz",
        "GSE274352_normalizedcounts_genes_IFNB_vs_empty.tsv.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274352/suppl/GSE274352_normalizedcounts_genes_IFNB_vs_empty.tsv.gz",
        "GSE175479_Raw_gene_count_matrix.txt.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175479/suppl/GSE175479_Raw_gene_count_matrix.txt.gz",
        "GSE164758_primary_tumors_fpkm.txt.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164758/suppl/GSE164758_primary_tumors_fpkm.txt.gz",
        "GSE193895_Tumor_KrasModel_GEMMs_RNAseq_GenewiseCounts.txt.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193895/suppl/GSE193895_Tumor_KrasModel_GEMMs_RNAseq_GenewiseCounts.txt.gz",
        "GSE338923_Lacun3_STK11_RNAseq_raw_counts.tsv.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE338nnn/GSE338923/suppl/GSE338923_Lacun3_STK11_RNAseq_raw_counts.tsv.gz",
        "GSE322570_raw_count.txt.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE322nnn/GSE322570/suppl/GSE322570_raw_count.txt.gz",
        "GSE133895_Sik_exvivo_Count_Matrix.txt.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133895/suppl/GSE133895_Sik_exvivo_Count_Matrix.txt.gz",
        "GSE244452_KPvsKL_deg_all.txt.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244452/suppl/GSE244452_KPvsKL_deg_all.txt.gz",
        "GSE6135-GPL8321_series_matrix.txt.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE6nnn/GSE6135/matrix/GSE6135-GPL8321_series_matrix.txt.gz",
        "GSE69552_series_matrix.txt.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE69nnn/GSE69552/matrix/GSE69552_series_matrix.txt.gz",
        "GSE21581_series_matrix.txt.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE21nnn/GSE21581/matrix/GSE21581_series_matrix.txt.gz",
        "GPL8321.annot.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL8nnn/GPL8321/annot/GPL8321.annot.gz",
        "GPL6887.annot.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6887/annot/GPL6887.annot.gz",
    }
    for name, url in urls.items():
        print("fetch", name)
        _download(url, DATA / name)

    print("ensembl map")
    emap = build_ensembl_map()
    (TAB / "ensembl_symbol_map_used.tsv").write_text(
        "ensembl_id\tsymbol\n" + "\n".join(f"{k}\t{v}" for k, v in sorted(emap.items()) if v in WANTED or any(v.lower()==w.lower() for w in WANTED))
    )

    contrast_rows: list[dict] = []
    spearman_rows: list[dict] = []
    presence_rows: list[dict] = []
    score_frames: list[pd.DataFrame] = []

    def add(rows, score, sp, pres):
        contrast_rows.extend(rows)
        score_frames.append(score)
        spearman_rows.extend(sp)
        presence_rows.extend(pres)

    # GSE137396 nodules, KL vs KP
    mat, scale = load_gse137396()
    groups = {
        "KL": [c for c in mat.columns if c.startswith("KL")],
        "KP": [c for c in mat.columns if c.startswith("KP")],
    }
    add(*score_contrast(
        "GSE137396", "autochthonous_nodule", "kl_vs_kp", mat, scale, groups, "KL", "KP",
        "nodule",
        "Same study family as GSE137244 but in vivo GEMM nodules, not the locked cell lines. Matrix is depositor log2, median-centered.",
    ))

    # GSE274351 LCM adenomas
    mat, scale = load_gse274351(emap)
    groups = {
        "KL": [c for c in mat.columns if c.startswith("KL")],
        "KP": [c for c in mat.columns if c.startswith("KP")],
        "K": [c for c in mat.columns if re.match(r"K\d", c)],
    }
    add(*score_contrast(
        "GSE274351", "lcm_adenoma", "kl_vs_kp", mat, scale,
        {"KL": groups["KL"], "KP": groups["KP"]}, "KL", "KP",
        "adenoma",
        "LCM lung adenomas, K/KP/KL. Normal lung excluded from this contrast.",
    ))
    add(*score_contrast(
        "GSE274351", "lcm_adenoma", "kl_vs_k", mat, scale,
        {"KL": groups["KL"], "K": groups["K"]}, "KL", "K",
        "adenoma",
        "Same LCM series, KL vs Kras-only. Not a p53 comparator.",
    ))

    # GSE274352. The two supplementary matrices are normalized separately.
    # Empty-vector KL vs KP uses the IFN-β file because it has KP1–3 and KL1–3.
    # STING-V154M pairs stay inside the STING file.
    mat, scale = load_gse274352_ifnb()
    empty = {
        "KL": [c for c in mat.columns if str(c).startswith("KL") and "empty" in str(c).lower()],
        "KP": [c for c in mat.columns if str(c).startswith("KP") and "empty" in str(c).lower()],
    }
    add(*score_contrast(
        "GSE274352", "cell_line_empty_vector", "kl_vs_kp", mat, scale, empty, "KL", "KP",
        "cell_line",
        "KP and KL GEMM-derived lines, empty-vector columns from the IFN-β matrix (KP1–3, KL1–3). Not GSE137244. Not pasted onto the STING matrix.",
    ))

    def _paired(mat_p: pd.DataFrame, kind: str, token: str) -> list[dict]:
        rows_p = []
        gene_mat, _ = resolve_genes(mat_p)
        treated = [c for c in mat_p.columns if token in str(c).upper()]
        for col in treated:
            line = str(col).split("_")[0]
            empties = [c for c in mat_p.columns if str(c).startswith(line + "_") and "empty" in str(c).lower()]
            if len(empties) != 1:
                continue
            e = empties[0]
            for ep, genes in {
                "Cldn4": None,
                "Sting1": None,
                "NHEJ": module_gene_list("NHEJ"),
                "STING_core": module_gene_list("STING_core"),
            }.items():
                if genes is None:
                    if ep not in gene_mat.index:
                        continue
                    ev, tv = float(gene_mat.loc[ep, e]), float(gene_mat.loc[ep, col])
                else:
                    use = [g for g in genes if g in gene_mat.index]
                    if not use:
                        continue
                    # Mean log2, not a z-score. A two-sample z-score is only a sign.
                    ev = float(gene_mat.loc[use, e].mean())
                    tv = float(gene_mat.loc[use, col].mean())
                rows_p.append(
                    {
                        "accession": "GSE274352",
                        "perturbation": kind,
                        "line": line,
                        "genotype": "KL" if line.startswith("KL") else "KP",
                        "endpoint": ep,
                        "empty": ev,
                        "treated": tv,
                        "delta_treated_minus_empty": tv - ev,
                    }
                )
        return rows_p

    sting_mat, _sting_scale = load_gse274352_sting()
    rescue_rows = _paired(sting_mat, "STING_V154M", "STING") + _paired(mat, "IFNB", "IFN")
    rescue = pd.DataFrame(rescue_rows)
    rescue.to_csv(TAB / "gse274352_paired_perturbation.tsv", sep="\t", index=False)
    # Paired summary
    paired_sum = []
    if len(rescue):
        for (pert, ep), sub in rescue.groupby(["perturbation", "endpoint"]):
            d = sub["delta_treated_minus_empty"].to_numpy(dtype=float)
            d = d[~np.isnan(d)]
            try:
                if len(d) >= 4 and np.any(d != 0):
                    w = stats.wilcoxon(d, alternative="two-sided", zero_method="wilcox")
                    p = float(w.pvalue)
                else:
                    p = float("nan")
            except ValueError:
                p = float("nan")
            paired_sum.append(
                {
                    "perturbation": pert,
                    "endpoint": ep,
                    "n_pairs": int(len(d)),
                    "n_kl": int((sub["genotype"] == "KL").sum()),
                    "n_kp": int((sub["genotype"] == "KP").sum()),
                    "median_delta": float(np.median(d)),
                    "mean_delta": float(np.mean(d)),
                    "n_up": int(np.sum(d > 0)),
                    "wilcoxon_p": p,
                    "note": "Paired within line on mean log2. Two-sided Wilcoxon at n=4 cannot go below 0.125.",
                }
            )
    pd.DataFrame(paired_sum).to_csv(TAB / "gse274352_paired_summary.tsv", sep="\t", index=False)

    # GSE164758 primary GEMM tumors
    mat, scale = load_gse164758()
    geno = {c: ("KPL" if c.startswith("KPL") else "KP" if c.startswith("KP") else "KL" if c.startswith("KL") else "K" if c.startswith("Kras") else "other") for c in mat.columns}
    groups = {}
    for gname in ("KL", "KP", "K", "KPL"):
        groups[gname] = [c for c, g in geno.items() if g == gname]
    add(*score_contrast(
        "GSE164758", "autochthonous_primary_tumor", "kl_vs_kp", mat, scale,
        {"KL": groups["KL"], "KP": groups["KP"]}, "KL", "KP",
        "tumor",
        "Untreated primary GEMM tumors. KPL is not pooled into KL.",
    ))
    add(*score_contrast(
        "GSE164758", "autochthonous_primary_tumor", "kl_vs_k", mat, scale,
        {"KL": groups["KL"], "K": groups["K"]}, "KL", "K",
        "tumor",
        "Same series, KL vs Kras-only. KPL excluded.",
    ))
    # KPL vs KP is a different genotype (triple mutant). Report, do not call it KL.
    add(*score_contrast(
        "GSE164758", "autochthonous_primary_tumor", "kpl_vs_kp", mat, scale,
        {"KPL": groups["KPL"], "KP": groups["KP"]}, "KPL", "KP",
        "tumor",
        "Triple-mutant KPL vs KP. Not a KL contrast. Shown so it is not silently pooled.",
    ))

    # GSE244452 syngeneic subcutaneous
    mat, scale = load_gse244452()
    add(*score_contrast(
        "GSE244452", "syngeneic_subcutaneous", "kl_vs_kp", mat, scale,
        {"KL": ["KL1", "KL2", "KL3"], "KP": ["KP1", "KP2", "KP3"]}, "KL", "KP",
        "tumor",
        "Subcutaneous syngeneic KP vs KL tumors. Per-sample columns inside the depositor all-gene table. n=3 cannot reach MW p<0.05.",
    ))

    # GSE322570 isogenic sgLkb1 on KP tumors, sorted tumor cells
    mat, scale = load_gse322570()
    add(*score_contrast(
        "GSE322570", "isogenic_sglkb1_on_KP", "kl_vs_kp", mat, scale,
        {
            "sgLkb1": [c for c in mat.columns if c.startswith("sgLkb1")],
            "sgNeo": [c for c in mat.columns if c.startswith("sgNeo")],
        },
        "sgLkb1", "sgNeo",
        "sorted_tumor_cells",
        "GFP+ CD45− tumor cells from KP tumors, sgLkb1 vs sgNeo. Isogenic, n=3.",
    ))

    # GSE338923 Lacun3 parental vs STK11 KO, in vitro
    mat, scale = load_gse338923()
    # Confirm KO by Stk11 itself: S* should be lower if those are the KO.
    gene_mat, _ = resolve_genes(mat)
    c_cols = [c for c in mat.columns if c.startswith("C")]
    s_cols = [c for c in mat.columns if c.startswith("S")]
    stk = gene_mat.loc["Stk11"] if "Stk11" in gene_mat.index else None
    if stk is not None and float(stk[s_cols].mean()) > float(stk[c_cols].mean()):
        # labels swapped relative to the C=control assumption
        ko_cols, wt_cols = c_cols, s_cols
        label_note = "Column groups were assigned by Stk11 expression because the S columns were not the lower-Stk11 group."
    else:
        ko_cols, wt_cols = s_cols, c_cols
        label_note = "C columns = parental (higher Stk11), S columns = Stk11 KO, confirmed by Stk11 expression."
    add(*score_contrast(
        "GSE338923", "isogenic_cellline_stk11ko", "isogenic_stk11", mat, scale,
        {"STK11_KO": ko_cols, "parental": wt_cols}, "STK11_KO", "parental",
        "cell_line",
        "Lacun3 murine LUAD, in vitro, parental vs CRISPR Stk11 KO. Not a KP-vs-KL GEMM. " + label_note,
    ))

    # GSE175479 K vs KL bulk lung
    mat, scale = load_gse175479()
    add(*score_contrast(
        "GSE175479", "autochthonous_bulk_lung", "kl_vs_k", mat, scale,
        {
            "KL": [c for c in mat.columns if c.startswith("KL")],
            "K": [c for c in mat.columns if c.startswith("K-")],
        },
        "KL", "K",
        "bulk_lung",
        "Kras vs Kras;Lkb1 bulk lung. AMPK-double-knockout (KAA) arm excluded. No KP arm.",
    ))

    # GSE133895 ex vivo KT vs KT;Lkb1
    mat, scale = load_gse133895()
    add(*score_contrast(
        "GSE133895", "exvivo_sorted", "kl_vs_k", mat, scale,
        {
            "KL": [c for c in mat.columns if "Lkb1" in c],
            "KT": [c for c in mat.columns if c.startswith("KT#")],
        },
        "KL", "KT",
        "sorted_tumor_cells",
        "Ex vivo sorted KT vs KT;Lkb1. Sik-knockout arms excluded. No KP arm.",
    ))

    # GSE193895 KL vs Kras/Keap1, mouse-level. Not KP.
    mat, scale = load_gse193895(emap)
    def _mouse_cols(prefix: str) -> dict[str, list[str]]:
        cols = [c for c in mat.columns if c.startswith(prefix + "_")]
        mice: dict[str, list[str]] = {}
        for c in cols:
            mid = c.split("_")[1]
            mice.setdefault(mid, []).append(c)
        return mice

    def _avg_mice(mice: dict[str, list[str]], tag: str) -> pd.DataFrame:
        pieces = []
        for mid, cols in mice.items():
            pieces.append(mat[cols].mean(axis=1).rename(f"{tag}{mid}"))
        return pd.concat(pieces, axis=1) if pieces else pd.DataFrame()

    kl_m = _avg_mice(_mouse_cols("KL"), "KL")
    kk_m = _avg_mice(_mouse_cols("KK"), "KK")
    mouse_mat = pd.concat([kl_m, kk_m], axis=1)
    add(*score_contrast(
        "GSE193895", "autochthonous_tumor_fragments", "kl_vs_keap1", mouse_mat, scale,
        {"KL": list(kl_m.columns), "KK": list(kk_m.columns)}, "KL", "KK",
        "mouse",
        "Kras/Lkb1 vs Kras/Keap1. Two tumor fragments per mouse were averaged. This is not a KP (p53) contrast. Normals and KKL excluded.",
    ))

    # Arrays
    amap_8321 = _annot_map(DATA / "GPL8321.annot.gz")
    amap_6887 = _annot_map(DATA / "GPL6887.annot.gz")

    expr, meta = _series_matrix(DATA / "GSE6135-GPL8321_series_matrix.txt.gz")
    expr.index = expr.index.astype(str).str.strip('"')
    genes = probes_to_genes(expr, amap_8321)
    genes, scale = maybe_log_array(genes)
    meta = gse6135_groups(meta)
    meta = meta.set_index("gsm")
    # Align columns
    genes = genes.loc[:, [c for c in genes.columns if c in meta.index]]
    # Mouse-level mean within genotype class
    mouse_frames = []
    mouse_groups = {"KL": [], "KP": [], "K": []}
    for geno in ("KL", "KP", "K"):
        sub = meta[meta["genotype"] == geno]
        for mid, rows in sub.groupby("mouse"):
            cols = [g for g in rows.index if g in genes.columns]
            if not cols:
                continue
            name = f"{geno}_{mid}"
            mouse_frames.append(genes[cols].mean(axis=1).rename(name))
            mouse_groups[geno].append(name)
    mouse_mat = pd.concat(mouse_frames, axis=1)
    add(*score_contrast(
        "GSE6135", "autochthonous_primary_tumor_microarray", "kl_vs_kp", mouse_mat, scale,
        {"KL": mouse_groups["KL"], "KP": mouse_groups["KP"]}, "KL", "KP",
        "mouse",
        "Ji et al. primary lung tumors. KL = Lkb1 L/L or L/- only. L/+ heterozygotes, metastases, and human lines excluded. Tumors from one mouse averaged.",
    ))
    add(*score_contrast(
        "GSE6135", "autochthonous_primary_tumor_microarray", "kl_vs_k", mouse_mat, scale,
        {"KL": mouse_groups["KL"], "K": mouse_groups["K"]}, "KL", "K",
        "mouse",
        "Same Ji series, KL vs Kras-only primary tumors.",
    ))

    expr, meta = _series_matrix(DATA / "GSE69552_series_matrix.txt.gz")
    expr.index = expr.index.astype(str).str.strip('"')
    genes = probes_to_genes(expr, amap_6887)
    genes, scale = maybe_log_array(genes)
    meta = meta.set_index("gsm")
    genes = genes.loc[:, [c for c in genes.columns if c in meta.index]]
    title = meta["title"]
    kl_all = [g for g in genes.columns if "Lkb1" in title[g] or "LKB1" in title[g]]
    kp_all = [g for g in genes.columns if "p53" in title[g]]
    kl_pap = [g for g in kl_all if "papillary" in title[g].lower()]
    kp_pap = [g for g in kp_all if "papillary" in title[g].lower()]
    add(*score_contrast(
        "GSE69552", "autochthonous_histology_microarray", "kl_vs_kp", genes, scale,
        {"KL": kl_all, "KP": kp_all}, "KL", "KP",
        "tumor",
        "Pre-specified primary: all Kras;Lkb1 (papillary + adenosquamous) vs Kras;p53 papillary.",
    ))
    add(*score_contrast(
        "GSE69552", "autochthonous_histology_microarray", "kl_vs_kp_papillary", genes, scale,
        {"KL": kl_pap, "KP": kp_pap}, "KL", "KP",
        "tumor",
        "Pre-specified sensitivity: papillary KL vs papillary KP only.",
    ))

    expr, meta = _series_matrix(DATA / "GSE21581_series_matrix.txt.gz")
    expr.index = expr.index.astype(str).str.strip('"')
    genes = probes_to_genes(expr, amap_8321)
    genes, scale = maybe_log_array(genes)
    meta = meta.set_index("gsm")
    genes = genes.loc[:, [c for c in genes.columns if c in meta.index]]
    title = meta["title"]
    kl_p = [g for g in genes.columns if "Lkb1 primary" in title[g]]
    k_p = [g for g in genes.columns if title[g].startswith("Kras primary")]
    add(*score_contrast(
        "GSE21581", "autochthonous_primary_tumor_microarray", "kl_vs_k", genes, scale,
        {"KL": kl_p, "K": k_p}, "KL", "K",
        "tumor",
        "Kras/Lkb1 primary vs Kras primary. Metastases excluded. No KP arm.",
    ))

    # FDR within pre-specified KL-vs-KP primary endpoints.
    contrasts = pd.DataFrame(contrast_rows)
    primary_eps = ["Cldn4", "Sting1", "NHEJ", "STING_core"]
    fdr_mask = (contrasts["family"] == "kl_vs_kp") & (contrasts["endpoint"].isin(primary_eps))
    # Split FDR: in vivo tumors vs cell/isogenic, so a cell-line test is not pooled with GEMMs.
    in_vivo_designs = {
        "autochthonous_nodule",
        "lcm_adenoma",
        "autochthonous_primary_tumor",
        "syngeneic_subcutaneous",
        "autochthonous_primary_tumor_microarray",
        "autochthonous_histology_microarray",
    }
    contrasts["fdr_family"] = ""
    contrasts["fdr_bh"] = np.nan
    for label, mask_extra in [
        ("in_vivo_kl_vs_kp", contrasts["design"].isin(in_vivo_designs)),
        ("cell_or_isogenic_kl_vs_kp", ~contrasts["design"].isin(in_vivo_designs)),
    ]:
        m = fdr_mask & mask_extra
        if m.any():
            qs = bh(contrasts.loc[m, "mannwhitney_p"].tolist())
            contrasts.loc[m, "fdr_bh"] = qs
            contrasts.loc[m, "fdr_family"] = label

    # Spearman FDR within family kl_vs_kp only, all Cldn4 correlations.
    spear = pd.DataFrame(spearman_rows)
    sm = spear["family"] == "kl_vs_kp"
    if sm.any():
        spear.loc[sm, "fdr_bh"] = bh(spear.loc[sm, "spearman_p"].tolist())
        spear.loc[sm, "fdr_family"] = "kl_vs_kp_spearman"
    else:
        spear["fdr_bh"] = np.nan

    contrasts.to_csv(TAB / "contrasts.tsv", sep="\t", index=False)
    spear.to_csv(TAB / "cldn4_spearman.tsv", sep="\t", index=False)
    pd.concat(score_frames, ignore_index=True).to_csv(TAB / "sample_scores.tsv", sep="\t", index=False)
    pd.DataFrame(presence_rows).to_csv(TAB / "genes_used.tsv", sep="\t", index=False)

    # Inventory (includes series that were not scored).
    inventory = [
        ["GSE137244", "KP vs KL GEMM cell lines", "yes", "yes", "yes", "no", "Locked baseline. Cldn4 KL−KP +5.57, MW p=0.00794, n=5 vs 5. Not re-opened."],
        ["GSE137396", "KP vs KL autochthonous nodules", "yes", "yes", "yes", "yes", "In vivo sister of the GSE137244 study. Scored."],
        ["GSE274351", "K/KP/KL LCM adenomas", "yes", "yes", "yes", "yes", "Scored KL vs KP and KL vs K. Normals excluded."],
        ["GSE274352", "KP/KL lines ± IFN-β or STING-V154M", "yes", "yes", "yes", "yes", "Empty-vector genotype contrast plus paired STING and IFN-β."],
        ["GSE164758", "K/KP/KL/KPL primary GEMM tumors", "yes", "yes", "yes", "yes", "Largest open in vivo n. KPL not pooled into KL."],
        ["GSE6135", "Ji et al. K/KP/KL primary tumors, microarray", "yes", "yes", "yes", "yes", "Mouse-level. Lkb1 L/+ and metastases excluded."],
        ["GSE69552", "KL papillary and adenosquamous vs KP papillary", "yes", "yes", "yes", "yes", "Primary = all KL vs KP. Sensitivity = papillary only."],
        ["GSE244452", "Syngeneic subcutaneous KP vs KL tumors", "yes", "yes", "yes", "yes", "3 vs 3. Per-sample values were in the all-gene table."],
        ["GSE322570", "sgLkb1 vs sgNeo sorted KP tumor cells", "KP background", "sgLkb1", "yes", "yes", "Isogenic on KP. n=3. CD45− tumor cells."],
        ["GSE338923", "Lacun3 parental vs Stk11 KO, in vitro", "parental LUAD line", "Stk11 KO", "yes", "yes", "Isogenic cell line. Not called a KP GEMM."],
        ["GSE175479", "K vs KL vs AMPK bulk lung", "no", "yes", "yes", "yes", "KL vs K only. KAA arm excluded."],
        ["GSE21581", "K vs KL primary and KL metastases", "no", "yes", "yes", "yes", "Primary tumors only. Metastases excluded."],
        ["GSE133895", "KT vs KT;Lkb1 ex vivo", "no", "yes", "yes", "yes", "Sik arms excluded. No KP."],
        ["GSE193895", "Kras/Lkb1 vs Kras/Keap1 tumors", "no", "yes", "yes", "yes", "Mouse-level. Comparator is Keap1, not p53."],
        ["GSE133714", "Kras ± Lkb1 ± Keap1 microarray", "no", "yes", "yes", "no", "No KP arm. GPL18233 has no standalone annot file (634 MB family soft). Not scored."],
        ["GSE165640", "Whole-lung KL/KP plus Il1-family knockouts", "yes", "yes", "no full matrix", "no", "Depositor files are pairwise compare tables, not a sample-by-gene matrix."],
        ["GSE118246", "Sox2;Lkb1 and Pten;Lkb1 squamous plus KP LUAD", "yes", "no classic KL", "yes", "no", "Lkb1 models are not Kras;Lkb1. Not a KL vs KP contrast."],
        ["GSE277929", "KL GEMM tumors ± entinostat/trametinib", "no", "yes", "yes", "no", "KL only. No KP comparator."],
        ["GSE253613", "KL ± G6pd", "no", "yes", "yes", "no", "KL only."],
        ["GSE182228", "Lkb1-deficient LUAD ± anti-PD-1", "no", "yes", "yes", "no", "KL only, and Cldn4 is near the floor in the prior slice."],
        ["GSE180963", "scRNA K vs KL", "no", "yes", "yes", "no", "n=1 per arm. Not merged with private 8 KL."],
        ["GSE165641", "scRNA KL only", "no", "yes", "yes", "no", "n=2 KL. Not merged with private 8 KL."],
        ["GSE154977", "scRNA KP ± cisplatin", "yes", "no", "yes", "no", "No KL. Not merged with private 8 KL."],
        ["GSE194166", "KL vs KP immune-sorted scRNA", "yes", "yes", "yes", "no", "CD45-sorted. Cldn4 is the wrong compartment."],
        ["GSE179500", "Lkb1-restorable bulk", "mixed", "restorable", "DE tables", "no", "Not a clean KP vs KL matrix. Not re-scored."],
        ["user_private_8KL", "Private 8 KL mice", "no", "private", "not public", "no", "Not used. Not merged with any public series."],
    ]
    inv = pd.DataFrame(inventory, columns=["accession", "what", "kp_arm", "kl_arm", "open_matrix", "scored_here", "note"])
    inv.to_csv(TAB / "inventory.tsv", sep="\t", index=False)

    # Forest of in vivo + cell KL vs KP for the four endpoints.
    plot_df = contrasts[
        (contrasts["family"] == "kl_vs_kp") & (contrasts["endpoint"].isin(primary_eps))
    ].copy()
    short = {
        "autochthonous_nodule": "nodules",
        "lcm_adenoma": "LCM adenomas",
        "cell_line_empty_vector": "cell line, empty",
        "autochthonous_primary_tumor": "primary tumors",
        "syngeneic_subcutaneous": "syngeneic s.c.",
        "isogenic_sglkb1_on_KP": "sgLkb1 on KP",
        "autochthonous_primary_tumor_microarray": "primary, array",
        "autochthonous_histology_microarray": "histology, array",
    }
    plot_df["label"] = plot_df["accession"] + "\n" + plot_df["design"].map(short).fillna(plot_df["design"])
    endpoints_order = primary_eps
    fig, axes = plt.subplots(1, len(endpoints_order), figsize=(12.5, 6.2), sharey=True)
    labels = list(dict.fromkeys(plot_df["label"]))
    ymap = {lab: i for i, lab in enumerate(labels)}
    for ax, ep in zip(axes, endpoints_order):
        sub = plot_df[plot_df["endpoint"] == ep]
        for _, r in sub.iterrows():
            y = ymap[r["label"]]
            ax.scatter(r["delta_kl_minus_ref"], y, s=36, color="#1f4e79")
            ax.text(
                r["delta_kl_minus_ref"],
                y + 0.18,
                f"p={r['mannwhitney_p']:.3g}" if pd.notna(r["mannwhitney_p"]) else "",
                fontsize=6,
                ha="center",
                color="#444",
            )
        ax.axvline(0, color="#999", lw=0.8)
        ax.set_title(ep, fontsize=10)
        unit = "z-score" if ep in ("NHEJ", "STING_core") else "log2"
        ax.set_xlabel(f"KL − comparator ({unit})")
    axes[0].set_yticks(list(ymap.values()))
    axes[0].set_yticklabels(list(ymap.keys()), fontsize=7)
    fig.suptitle("Public mouse lung, beyond GSE137244\nKL (or Lkb1-loss) minus KP / Lkb1-intact", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG / "forest_kl_minus_kp.png", dpi=140)
    plt.close()

    # Compact JSON for the write-up.
    keep = contrasts[contrasts["endpoint"].isin(primary_eps + ["Stk11", "Cgas", "STING_ISG"])]
    summary = {
        "private_8kl": "not used",
        "gse137244": "not re-opened; locked Cldn4 KL-KP +5.57, MW p=0.00794, n=5 vs 5",
        "n_contrasts": int(len(contrasts)),
        "primary_kl_vs_kp": keep[keep["family"] == "kl_vs_kp"][
            ["accession", "design", "endpoint", "n_kl", "n_ref", "delta_kl_minus_ref", "mannwhitney_p", "fdr_bh", "hedges_g"]
        ].to_dict(orient="records"),
    }
    (TAB / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(contrasts[contrasts["endpoint"].isin(primary_eps)][
        ["accession", "family", "endpoint", "n_kl", "n_ref", "delta_kl_minus_ref", "mannwhitney_p", "fdr_bh"]
    ].to_string(index=False))
    print("wrote", TAB)


if __name__ == "__main__":
    main()
