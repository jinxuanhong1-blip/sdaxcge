#!/usr/bin/env python3
"""ESTIMATE-purity residualization of TACSTD2, CLDN4, and a TJ score vs ICI DCB/ORR.

Rework of the open lung ICI bulk slice (raw TACSTD2/CLDN4 vs response was NS).
For every requested GEO series, this script:

  1. Downloads official processed matrices + metadata (no FASTQ).
  2. Computes ESTIMATE Stromal/Immune/ESTIMATEScore/TumorPurity
     (Python port of estimate R v1.0.13; Yoshihara et al. 2013).
  3. Residualizes TACSTD2, CLDN4, and a locked tight-junction (TJ) score
     on ESTIMATEScore by ordinary least squares (primary). TumorPurity
     residualization is a sensitivity when the Affymetrix purity formula
     is in-bounds.
  4. Tests raw and residualized values against the cohort-native binary
     endpoint (DCB or ORR / author responder label / RECIST / MPR).
  5. Reports every requested cohort, including those that cannot be tested.

Requested series (all reported):
  GSE126044, GSE135222, GSE166449, GSE190265, GSE207422.

Outputs: results/rework/ICI_purity/
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2] if HERE.parts[-2:] == ("rework", "ICI_purity") else Path("/workspace")
OUT = Path(os.environ.get("ICI_PURITY_OUT", REPO / "results" / "rework" / "ICI_purity"))
DATA = Path(os.environ.get("ICI_PURITY_DATA", "/tmp/ici_purity_data"))
GMT_PATH = HERE / "SI_geneset.gmt"
HGNC_PATH = DATA / "hgnc_complete_set.txt"

TABLES = OUT / "tables"
FIGURES = OUT / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

SEED = 20260816
N_BOOT = 2000
DCB_DAYS = 183
DCB_MONTHS = 6.0

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

FILES = {
    "GSE126044_series_matrix.txt.gz": f"{FTP}/GSE126nnn/GSE126044/matrix/GSE126044_series_matrix.txt.gz",
    "GSE126044_counts.txt.gz": f"{FTP}/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz",
    "GSE135222_series_matrix.txt.gz": f"{FTP}/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz",
    "GSE135222_exp.tsv.gz": f"{FTP}/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz",
    "GSE166449_series_matrix.txt.gz": f"{FTP}/GSE166nnn/GSE166449/matrix/GSE166449_series_matrix.txt.gz",
    "GSE166449_TPM.txt.gz": f"{FTP}/GSE166nnn/GSE166449/suppl/GSE166449_Raw_gene_TPM_matrix.txt.gz",
    "GSE190265_series_matrix.txt.gz": f"{FTP}/GSE190nnn/GSE190265/matrix/GSE190265_series_matrix.txt.gz",
    "GSE190265_TPM_France3.csv.gz": f"{FTP}/GSE190nnn/GSE190265/suppl/GSE190265_TPM_France3.csv.gz",
    "GSE190265_samples_info_France3.csv.gz": f"{FTP}/GSE190nnn/GSE190265/suppl/GSE190265_samples_info_France3.csv.gz",
    "GSE207422_series_matrix.txt.gz": f"{FTP}/GSE207nnn/GSE207422/matrix/GSE207422_series_matrix.txt.gz",
    "GSE207422_log2TPM.txt.gz": f"{FTP}/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
    "GSE207422_metadata.xlsx": f"{FTP}/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx",
}

# Locked TJ programme used across this project (not a single gene).
TJ_5 = ["CLDN1", "CLDN4", "CLDN7", "F11R", "PARD3"]
# Structural TJ set from the B3 purity slice (sensitivity).
TJ_STRUCT = [
    "CLDN1", "CLDN3", "CLDN4", "CLDN7",
    "OCLN",
    "TJP1", "TJP2", "TJP3",
    "F11R", "JAM2", "JAM3",
    "MARVELD2", "MARVELD3",
    "CGN", "CGNL1",
]
TARGETS = ["TACSTD2", "CLDN4"]
CD8_GENES = ["CD8A", "CD8B"]
ENSG = {
    "TACSTD2": "ENSG00000184292",
    "CLDN4": "ENSG00000189143",
    "CLDN1": "ENSG00000163347",
    "CLDN7": "ENSG00000181885",
    "F11R": "ENSG00000158769",
    "PARD3": "ENSG00000148498",
    "CD8A": "ENSG00000153563",
    "CD8B": "ENSG00000172116",
}

SYNONYMS = {
    "F11R": ["JAM1", "JAMA", "JAM-A"],
    "PATJ": ["INADL"],
    "LPPR4": ["PLPPR4"],
    "GPR124": ["ADGRA2"],
    "WISP1": ["CCN4"],
}


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_err = None
    for attempt in range(1, 5):
        try:
            urllib.request.urlretrieve(url, tmp)
            if tmp.exists() and tmp.stat().st_size > 0:
                tmp.replace(dest)
                return dest
        except Exception as exc:  # noqa: BLE001
            last_err = exc
        if tmp.exists():
            tmp.unlink(missing_ok=True)
    raise RuntimeError(f"download failed: {url} ({last_err})")


def parse_series_matrix(path: Path) -> list[dict]:
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
    samples = []
    for i, title in enumerate(titles or []):
        d = {"title": title}
        if geo and i < len(geo):
            d["geo_accession"] = geo[i]
        if descriptions and i < len(descriptions):
            d["description"] = descriptions[i]
        for row in char_rows:
            if i < len(row) and row[i] and ":" in row[i]:
                k, v = row[i].split(":", 1)
                d[k.strip().lower()] = v.strip()
        samples.append(d)
    return samples


def load_gmt(path: Path) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    for line in path.read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def load_hgnc() -> pd.DataFrame:
    fetch(HGNC_URL, HGNC_PATH)
    return pd.read_csv(HGNC_PATH, sep="\t", dtype=str, low_memory=False)


def build_symbol_maps(hgnc: pd.DataFrame) -> tuple[dict[str, str], dict[str, str]]:
    """alias/prev/current symbol -> current symbol; current symbol -> ensembl."""
    alias_to_current: dict[str, str] = {}
    symbol_to_ensg: dict[str, str] = {}
    for _, row in hgnc.iterrows():
        sym = row.get("symbol")
        if not isinstance(sym, str) or not sym:
            continue
        alias_to_current[sym.upper()] = sym
        ensg = row.get("ensembl_gene_id")
        if isinstance(ensg, str) and ensg.startswith("ENSG"):
            symbol_to_ensg[sym] = ensg
        for col in ("prev_symbol", "alias_symbol"):
            raw = row.get(col)
            if isinstance(raw, str) and raw:
                for a in raw.split("|"):
                    a = a.strip()
                    if a:
                        alias_to_current.setdefault(a.upper(), sym)
    for old, news in SYNONYMS.items():
        for n in news:
            alias_to_current.setdefault(old.upper(), n)
            alias_to_current.setdefault(n.upper(), n)
    return alias_to_current, symbol_to_ensg


# ---------------------------------------------------------------------------
# ESTIMATE (exact port of estimateScore, R package v1.0.13)
# ---------------------------------------------------------------------------
def estimate_scores(expr: pd.DataFrame, gene_sets: dict[str, list[str]], resolved: dict[str, str]):
    """expr = genes x samples on any monotone-equivalent scale.

    `resolved` maps GMT symbols to row names actually present in `expr`.
    """
    genes = expr.index.to_numpy()
    n_genes = expr.shape[0]
    m = expr.rank(axis=0, method="average").to_numpy() * (10000.0 / n_genes)
    rows = {}
    overlaps = {}
    for name, out_col in (("StromalSignature", "StromalScore"), ("ImmuneSignature", "ImmuneScore")):
        mapped = [resolved[g] for g in gene_sets[name] if g in resolved]
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
    arg = 0.6049872018 + 0.0001467884 * df["ESTIMATEScore"].to_numpy()
    purity = np.cos(arg)
    oob = (arg < 0) | (arg > np.pi) | (purity < 0)
    purity = purity.astype(float)
    purity[oob] = np.nan
    df["TumorPurity"] = purity
    return df, overlaps


def resolve_gmt_genes(gmt_genes: list[str], expr_index: pd.Index, alias_to_current: dict[str, str],
                      symbol_to_ensg: dict[str, str], index_is_ensembl: bool) -> dict[str, str]:
    """Map GMT symbols onto expr row names. Returns gmt_symbol -> row name."""
    idx_upper = {str(g).upper(): g for g in expr_index}
    idx_ensg = {}
    if index_is_ensembl:
        for g in expr_index:
            idx_ensg[str(g).split(".")[0]] = g
    resolved = {}
    for g in gmt_genes:
        current = alias_to_current.get(g.upper(), g)
        if index_is_ensembl:
            ensg = symbol_to_ensg.get(current) or symbol_to_ensg.get(g) or ENSG.get(current) or ENSG.get(g)
            if ensg and ensg in idx_ensg:
                resolved[g] = idx_ensg[ensg]
                continue
        if current.upper() in idx_upper:
            resolved[g] = idx_upper[current.upper()]
        elif g.upper() in idx_upper:
            resolved[g] = idx_upper[g.upper()]
    return resolved


# ---------------------------------------------------------------------------
# Stats
# ---------------------------------------------------------------------------
def residualize(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """OLS residual of y on x with intercept. NaNs in x or y become NaN residuals."""
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


def zscore(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, float)
    sd = np.nanstd(v)
    if sd == 0 or np.isnan(sd):
        return np.zeros_like(v)
    return (v - np.nanmean(v)) / sd


def signature_score(expr: pd.DataFrame, genes: list[str], row_map: dict[str, str]) -> tuple[pd.Series, list[str]]:
    found = [row_map[g] for g in genes if g in row_map]
    if not found:
        return pd.Series(np.nan, index=expr.columns), []
    z = expr.loc[found].apply(lambda r: zscore(r.to_numpy()), axis=1, result_type="expand")
    z.columns = expr.columns
    return z.mean(axis=0), found


def map_wanted(genes: list[str], expr_index: pd.Index, alias_to_current: dict[str, str],
               symbol_to_ensg: dict[str, str], index_is_ensembl: bool) -> dict[str, str]:
    return resolve_gmt_genes(genes, expr_index, alias_to_current, symbol_to_ensg, index_is_ensembl)


def mwu_block(values: np.ndarray, labels: np.ndarray, pos_name: str, neg_name: str) -> dict:
    """labels: 1 = benefit/response (pos), 0 = no benefit."""
    v = np.asarray(values, float)
    lab = np.asarray(labels, int)
    m = ~np.isnan(v) & np.isin(lab, [0, 1])
    v, lab = v[m], lab[m]
    pos, neg = v[lab == 1], v[lab == 0]
    n1, n0 = len(pos), len(neg)
    if n1 < 2 or n0 < 2:
        return {
            "n_pos": n1, "n_neg": n0,
            "median_pos": float(np.median(pos)) if n1 else np.nan,
            "median_neg": float(np.median(neg)) if n0 else np.nan,
            "delta_median_pos_minus_neg": np.nan,
            "mannwhitney_U": np.nan, "auc_pos_gt_neg": np.nan,
            "rank_biserial": np.nan, "p_two_sided": np.nan,
            "pos_class": pos_name, "neg_class": neg_name,
        }
    U, p = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    auc = float(U / (n1 * n0))
    return {
        "n_pos": n1, "n_neg": n0,
        "median_pos": float(np.median(pos)),
        "median_neg": float(np.median(neg)),
        "delta_median_pos_minus_neg": float(np.median(pos) - np.median(neg)),
        "mannwhitney_U": float(U),
        "auc_pos_gt_neg": auc,
        "rank_biserial": float(2 * auc - 1),
        "p_two_sided": float(p),
        "pos_class": pos_name,
        "neg_class": neg_name,
    }


def spearman_safe(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    if m.sum() < 3:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def bootstrap_auc_ci(values, labels, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    v = np.asarray(values, float)
    lab = np.asarray(labels, int)
    m = ~np.isnan(v) & np.isin(lab, [0, 1])
    v, lab = v[m], lab[m]
    pos_idx, neg_idx = np.where(lab == 1)[0], np.where(lab == 0)[0]
    if len(pos_idx) < 2 or len(neg_idx) < 2:
        return np.nan, np.nan
    aucs = []
    for _ in range(n_boot):
        idx = np.concatenate(
            [rng.choice(pos_idx, len(pos_idx), replace=True),
             rng.choice(neg_idx, len(neg_idx), replace=True)]
        )
        pos, neg = v[idx][lab[idx] == 1], v[idx][lab[idx] == 0]
        U = stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic
        aucs.append(U / (len(pos) * len(neg)))
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return float(lo), float(hi)


# ---------------------------------------------------------------------------
# Cohort loaders: return expr (genes x samples, log-like), labels DataFrame
# ---------------------------------------------------------------------------
def load_gse126044():
    meta = parse_series_matrix(DATA / "GSE126044_series_matrix.txt.gz")
    resp = {}
    for s in meta:
        col = s["title"].replace("RNA-seq_", "")
        resp[col] = s.get("patient response")
    counts = pd.read_csv(DATA / "GSE126044_counts.txt.gz", sep="\t", index_col=0)
    counts = counts[~counts.index.duplicated(keep="first")]
    cpm = counts.divide(counts.sum(axis=0), axis=1) * 1e6
    expr = np.log2(cpm + 1)
    samples = [c for c in expr.columns if resp.get(c) in ("responder", "non-responder")]
    lab = pd.DataFrame({
        "sample": samples,
        "endpoint": "ORR_author_response",
        "endpoint_note": "GEO characteristic 'patient response' (Cho et al.); not a derived DCB",
        "label": [1 if resp[c] == "responder" else 0 for c in samples],
        "label_name": [resp[c] for c in samples],
    }).set_index("sample")
    return expr[samples], lab, "log2CPM", False


def load_gse135222():
    meta = parse_series_matrix(DATA / "GSE135222_series_matrix.txt.gz")
    recs = []
    for s in meta:
        col = s["title"].replace(" ", "")
        ev, t = s.get("progression-free survival (pfs)"), s.get("pfs.time")
        if ev is None or t is None:
            continue
        recs.append({"sample": col, "pfs_event": int(ev), "pfs_days": float(t)})
    clin = pd.DataFrame(recs).set_index("sample")
    early = clin[(clin.pfs_event == 0) & (clin.pfs_days < DCB_DAYS)]
    tpm = pd.read_csv(DATA / "GSE135222_exp.tsv.gz", sep="\t", index_col=0)
    tpm.index = [str(i).split(".")[0] for i in tpm.index]
    tpm = tpm[~pd.Index(tpm.index).duplicated(keep="first")]
    expr = np.log2(tpm.clip(lower=0) + 1)
    common = [c for c in expr.columns if c in clin.index]
    clin = clin.loc[common]
    clin["label"] = (clin.pfs_days >= DCB_DAYS).astype(int)
    clin["label_name"] = np.where(clin["label"] == 1, "DCB", "NDB")
    clin["endpoint"] = "DCB_PFS_ge_183d"
    clin["endpoint_note"] = (
        f"DCB = PFS >= {DCB_DAYS} days; "
        f"early-censored before cutoff: {len(early)} (none expected)"
    )
    return expr[common], clin, "log2(TPM+1)", True


def load_gse166449():
    meta = parse_series_matrix(DATA / "GSE166449_series_matrix.txt.gz")
    resp = {}
    for s in meta:
        col = s.get("description")
        title = (s.get("title") or "").lower()
        if "nonresponder" in title:
            resp[col] = "non-responder"
        elif "responder" in title:
            resp[col] = "responder"
    tpm = pd.read_csv(DATA / "GSE166449_TPM.txt.gz", sep="\t", index_col=0)
    tpm = tpm[~tpm.index.duplicated(keep="first")]
    expr = np.log2(tpm.clip(lower=0) + 1)
    samples = [c for c in expr.columns if resp.get(c) in ("responder", "non-responder")]
    lab = pd.DataFrame({
        "sample": samples,
        "endpoint": "ORR_author_response",
        "endpoint_note": "responder/nonresponder parsed from GEO sample title",
        "label": [1 if resp[c] == "responder" else 0 for c in samples],
        "label_name": [resp[c] for c in samples],
    }).set_index("sample")
    return expr[samples], lab, "log2(TPM+1)", False


def load_gse190265():
    path = DATA / "GSE190265_TPM_France3.csv.gz"
    with gzip.open(path, "rt") as fh:
        genes = fh.readline().strip().split(";")
    tpm = pd.read_csv(path, sep=";", header=None, skiprows=1, index_col=0)
    tpm.columns = genes
    tpm = tpm.T
    tpm = tpm[~tpm.index.duplicated(keep="first")]
    expr = np.log2(tpm.clip(lower=0).astype(float) + 1)
    info = pd.read_csv(DATA / "GSE190265_samples_info_France3.csv.gz", sep=";").set_index("sample")
    series = parse_series_matrix(DATA / "GSE190265_series_matrix.txt.gz")
    series_ids = {s["title"] for s in series}
    common = [c for c in expr.columns if c in info.index]
    clin = info.loc[common].copy()
    early = clin[(clin.evtPFS == 0) & (clin.time_PFS < DCB_MONTHS)]
    clin["label"] = (clin.time_PFS >= DCB_MONTHS).astype(int)
    clin["label_name"] = np.where(clin["label"] == 1, "DCB", "NDB")
    clin["endpoint"] = "DCB_PFS_ge_6mo"
    clin["endpoint_note"] = (
        f"DCB = time_PFS >= {DCB_MONTHS} months from GEO samples_info; "
        f"early-censored: {len(early)}; "
        f"n_in_series_matrix={sum(c in series_ids for c in common)}/{len(common)}"
    )
    clin["in_series_matrix"] = [c in series_ids for c in clin.index]
    return expr[common], clin, "log2(TPM+1)", False


def load_gse207422():
    md = pd.read_excel(DATA / "GSE207422_metadata.xlsx")
    md = md.dropna(subset=["Sample"]).copy()
    md = md[md["Sample"].astype(str).str.startswith("R")]
    expr_raw = pd.read_csv(DATA / "GSE207422_log2TPM.txt.gz", sep="\t", index_col=0)
    expr_raw = expr_raw[~expr_raw.index.duplicated(keep="first")]
    # Primary endpoint: RECIST ORR (CR/PR vs SD/PD). MPR kept as secondary.
    recist_map = {"CR": 1, "PR": 1, "SD": 0, "PD": 0}
    md["orr_label"] = md["RECIST"].map(recist_map)
    def mpr_lab(x):
        x = str(x)
        if x.startswith("NMPR"):
            return 0
        if "MPR" in x:
            return 1
        return np.nan
    md["mpr_label"] = md["Pathologic Response"].map(mpr_lab)
    md = md.set_index("Sample")
    samples = [c for c in expr_raw.columns if c in md.index and pd.notna(md.loc[c, "orr_label"])]
    lab = md.loc[samples, ["orr_label", "mpr_label", "RECIST", "Pathologic Response", "Residual Tumor"]].copy()
    lab["label"] = lab["orr_label"].astype(int)
    lab["label_name"] = np.where(lab["label"] == 1, "ORR_CR_PR", "nonORR_SD_PD")
    lab["endpoint"] = "ORR_RECIST"
    lab["endpoint_note"] = "RECIST CR/PR vs SD/PD from GEO bulk metadata xlsx; MPR stored as secondary"
    return expr_raw[samples], lab, "log2TPM", False


LOADERS = {
    "GSE126044": load_gse126044,
    "GSE135222": load_gse135222,
    "GSE166449": load_gse166449,
    "GSE190265": load_gse190265,
    "GSE207422": load_gse207422,
}


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def box_panel(ax, values, labels, title, ylabel, pos_name, neg_name):
    rng = np.random.default_rng(SEED)
    v = np.asarray(values, float)
    lab = np.asarray(labels, int)
    groups = [v[lab == 0], v[lab == 1]]
    names = [f"{neg_name}\n(n={np.isfinite(groups[0]).sum()})",
             f"{pos_name}\n(n={np.isfinite(groups[1]).sum()})"]
    ax.boxplot(groups, tick_labels=names, showfliers=False, widths=0.6)
    for i, g in enumerate(groups, start=1):
        g = g[np.isfinite(g)]
        ax.scatter(np.full(len(g), i) + rng.uniform(-0.08, 0.08, len(g)),
                   g, s=22, alpha=0.8, color="#1f77b4" if i == 1 else "#d62728", zorder=3)
    ax.set_title(title, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=8)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("Downloading missing files…")
    manifest = []
    for name, url in FILES.items():
        dest = DATA / name
        fetch(url, dest)
        manifest.append({"file": name, "bytes": dest.stat().st_size, "sha256": sha256(dest), "url": url})
    if not GMT_PATH.exists():
        raise FileNotFoundError(f"ESTIMATE GMT missing: {GMT_PATH}")
    gene_sets = load_gmt(GMT_PATH)
    hgnc = load_hgnc()
    alias_to_current, symbol_to_ensg = build_symbol_maps(hgnc)

    inventory_rows = []
    test_rows = []
    purity_rows = []
    sample_rows = []
    gene_cov_rows = []

    primary_features = ["TACSTD2", "CLDN4", "TJ_5"]

    for gse, loader in LOADERS.items():
        print(f"\n===== {gse} =====")
        try:
            expr, lab, unit, ensembl = loader()
            status = "analyzed"
            skip_reason = ""
        except Exception as exc:  # noqa: BLE001
            inventory_rows.append({
                "dataset": gse, "status": "failed_to_load", "skip_reason": str(exc),
                "n_expr_samples": 0, "n_labeled": 0, "n_pos": 0, "n_neg": 0,
                "endpoint": "", "unit": "", "bulk_present": False,
            })
            print("LOAD FAIL", exc)
            continue

        n_expr = expr.shape[1]
        n_lab = int(lab["label"].notna().sum()) if "label" in lab else 0
        n_pos = int((lab["label"] == 1).sum()) if n_lab else 0
        n_neg = int((lab["label"] == 0).sum()) if n_lab else 0
        endpoint = lab["endpoint"].iloc[0] if n_lab else ""
        note = lab["endpoint_note"].iloc[0] if n_lab else ""

        # Map GMT + analysis genes
        gmt_all = gene_sets["StromalSignature"] + gene_sets["ImmuneSignature"]
        gmt_map = resolve_gmt_genes(gmt_all, expr.index, alias_to_current, symbol_to_ensg, ensembl)
        tgt_map = map_wanted(TARGETS + TJ_5 + TJ_STRUCT + CD8_GENES, expr.index,
                             alias_to_current, symbol_to_ensg, ensembl)

        missing_targets = [g for g in TARGETS if g not in tgt_map]
        if missing_targets:
            status = "genes_missing"
            skip_reason = f"missing {missing_targets}"
        if n_pos < 2 or n_neg < 2:
            status = "insufficient_labels" if status == "analyzed" else status
            skip_reason = (skip_reason + "; " if skip_reason else "") + f"n_pos={n_pos} n_neg={n_neg}"

        inventory_rows.append({
            "dataset": gse, "status": status, "skip_reason": skip_reason,
            "n_expr_samples": n_expr, "n_genes": int(expr.shape[0]),
            "n_labeled": n_lab, "n_pos": n_pos, "n_neg": n_neg,
            "endpoint": endpoint, "endpoint_note": note, "unit": unit,
            "bulk_present": True, "index_is_ensembl": ensembl,
            "estimate_stromal_overlap": sum(g in gmt_map for g in gene_sets["StromalSignature"]),
            "estimate_immune_overlap": sum(g in gmt_map for g in gene_sets["ImmuneSignature"]),
            "TACSTD2_present": "TACSTD2" in tgt_map,
            "CLDN4_present": "CLDN4" in tgt_map,
            "TJ_5_found": sum(g in tgt_map for g in TJ_5),
            "TJ_struct_found": sum(g in tgt_map for g in TJ_STRUCT),
        })
        for g in TARGETS + TJ_5 + TJ_STRUCT + CD8_GENES:
            gene_cov_rows.append({
                "dataset": gse, "gene": g, "found": g in tgt_map,
                "row_name": tgt_map.get(g, ""),
            })

        if status != "analyzed":
            print("SKIP", status, skip_reason)
            continue

        est, overlaps = estimate_scores(expr, gene_sets, gmt_map)
        print("ESTIMATE overlap", overlaps, "purity OOB", int(est["TumorPurity"].isna().sum()))

        tj5, tj5_found = signature_score(expr, TJ_5, tgt_map)
        tjs, tjs_found = signature_score(expr, TJ_STRUCT, tgt_map)
        cd8, cd8_found = signature_score(expr, CD8_GENES, tgt_map)

        features = {}
        for g in TARGETS:
            features[g] = expr.loc[tgt_map[g]]
        features["TJ_5"] = tj5
        features["TJ_structural"] = tjs
        features["CD8"] = cd8

        residuals_score = {k: pd.Series(residualize(v.to_numpy(), est["ESTIMATEScore"].to_numpy()),
                                        index=v.index) for k, v in features.items()}
        residuals_pur = {k: pd.Series(residualize(v.to_numpy(), est["TumorPurity"].to_numpy()),
                                      index=v.index) for k, v in features.items()}

        labels = lab["label"].reindex(expr.columns)
        pos_name = "benefit"
        neg_name = "no_benefit"
        if endpoint.startswith("DCB"):
            pos_name, neg_name = "DCB", "NDB"
        elif endpoint.startswith("ORR"):
            pos_name, neg_name = "ORR+", "ORR-"
        elif "responder" in endpoint:
            pos_name, neg_name = "responder", "non-responder"

        for feat, raw in features.items():
            rho_s, p_s, n_s = spearman_safe(raw, est["ESTIMATEScore"])
            rho_p, p_p, n_p = spearman_safe(raw, est["TumorPurity"])
            rho_i, p_i, _ = spearman_safe(raw, est["ImmuneScore"])
            purity_rows.append({
                "dataset": gse, "feature": feat, "unit": unit,
                "spearman_vs_ESTIMATEScore": rho_s, "p_vs_ESTIMATEScore": p_s, "n_vs_score": n_s,
                "spearman_vs_TumorPurity": rho_p, "p_vs_TumorPurity": p_p, "n_vs_purity": n_p,
                "spearman_vs_ImmuneScore": rho_i, "p_vs_ImmuneScore": p_i,
                "n_purity_OOB": int(est["TumorPurity"].isna().sum()),
            })
            for measure, series in (
                ("raw", raw),
                ("residual_ESTIMATEScore", residuals_score[feat]),
                ("residual_TumorPurity", residuals_pur[feat]),
            ):
                blk = mwu_block(series.reindex(labels.index).to_numpy(),
                                labels.to_numpy(), pos_name, neg_name)
                lo, hi = bootstrap_auc_ci(series.reindex(labels.index).to_numpy(),
                                          labels.to_numpy())
                blk.update({
                    "dataset": gse, "feature": feat, "measure": measure,
                    "unit": unit if measure == "raw" else f"{unit}_resid",
                    "endpoint": endpoint, "auc_95ci_lo": lo, "auc_95ci_hi": hi,
                    "primary": feat in primary_features and measure == "residual_ESTIMATEScore",
                })
                test_rows.append(blk)

        # GSE207422 secondary: MPR
        if gse == "GSE207422" and "mpr_label" in lab.columns:
            mpr = lab["mpr_label"].reindex(expr.columns)
            for feat, raw in features.items():
                for measure, series in (
                    ("raw", raw),
                    ("residual_ESTIMATEScore", residuals_score[feat]),
                ):
                    blk = mwu_block(series.reindex(mpr.index).to_numpy(),
                                    mpr.to_numpy(), "MPR", "NMPR")
                    lo, hi = bootstrap_auc_ci(series.reindex(mpr.index).to_numpy(),
                                              mpr.to_numpy())
                    blk.update({
                        "dataset": gse, "feature": feat, "measure": measure,
                        "unit": unit if measure == "raw" else f"{unit}_resid",
                        "endpoint": "MPR_pathologic",
                        "auc_95ci_lo": lo, "auc_95ci_hi": hi, "primary": False,
                    })
                    test_rows.append(blk)

        # per-sample table
        for s in expr.columns:
            row = {
                "dataset": gse, "sample": s, "unit": unit,
                "endpoint": endpoint,
                "label": int(lab.loc[s, "label"]) if s in lab.index and pd.notna(lab.loc[s, "label"]) else np.nan,
                "label_name": lab.loc[s, "label_name"] if s in lab.index else "",
                "StromalScore": float(est.loc[s, "StromalScore"]),
                "ImmuneScore": float(est.loc[s, "ImmuneScore"]),
                "ESTIMATEScore": float(est.loc[s, "ESTIMATEScore"]),
                "TumorPurity": float(est.loc[s, "TumorPurity"]) if pd.notna(est.loc[s, "TumorPurity"]) else np.nan,
            }
            for feat in features:
                row[f"{feat}_raw"] = float(features[feat].loc[s])
                row[f"{feat}_resid_ESTIMATEScore"] = float(residuals_score[feat].loc[s])
                row[f"{feat}_resid_TumorPurity"] = float(residuals_pur[feat].loc[s])
            sample_rows.append(row)

        # figures: raw vs residual boxplots for primary features
        fig, axes = plt.subplots(2, 3, figsize=(11, 7))
        for j, feat in enumerate(primary_features):
            box_panel(axes[0, j], features[feat].reindex(labels.index), labels,
                      f"{gse} {feat} raw", unit, pos_name, neg_name)
            box_panel(axes[1, j], residuals_score[feat].reindex(labels.index), labels,
                      f"{gse} {feat} residual|ESTIMATEScore", "OLS residual", pos_name, neg_name)
        fig.suptitle(f"{gse}: raw vs ESTIMATEScore-residualized vs {endpoint}", fontsize=11)
        fig.tight_layout()
        fig.savefig(FIGURES / f"{gse}_raw_vs_residual_box.png", dpi=150)
        plt.close(fig)

        fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
        for ax, feat in zip(axes, primary_features):
            x = est["ESTIMATEScore"].reindex(labels.index)
            y = features[feat].reindex(labels.index)
            for val, col, name in ((0, "#1f77b4", neg_name), (1, "#d62728", pos_name)):
                m = labels == val
                ax.scatter(x[m], y[m], s=28, color=col, label=name, alpha=0.85,
                           edgecolor="k", linewidth=0.3)
            r, p, n = spearman_safe(x, y)
            ax.set_xlabel("ESTIMATEScore")
            ax.set_ylabel(feat)
            ax.set_title(f"{feat} vs purity proxy\nρ={r:.2f} p={p:.3f} n={n}", fontsize=9)
            ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(FIGURES / f"{gse}_feature_vs_ESTIMATEScore.png", dpi=150)
        plt.close(fig)

        print(f"n={n_lab} pos={n_pos} neg={n_neg} endpoint={endpoint}")

    # ---- assemble tables ----
    inv = pd.DataFrame(inventory_rows)
    tests = pd.DataFrame(test_rows)
    pur = pd.DataFrame(purity_rows)
    samples = pd.DataFrame(sample_rows)
    gcov = pd.DataFrame(gene_cov_rows)

    if not tests.empty:
        prim = tests["primary"].fillna(False).astype(bool)
        q = np.full(len(tests), np.nan)
        if prim.any():
            p = tests.loc[prim, "p_two_sided"].to_numpy()
            _, q_prim, _, _ = multipletests(p, method="fdr_bh")
            q[np.where(prim)[0]] = q_prim
        tests["q_BH_primary_residual"] = q
        # also within-cohort FDR on the 3 primary residual tests
        q2 = np.full(len(tests), np.nan)
        for ds, sub in tests[prim].groupby("dataset"):
            _, qq, _, _ = multipletests(sub["p_two_sided"], method="fdr_bh")
            q2[sub.index] = qq
        tests["q_BH_within_cohort"] = q2

    inv.to_csv(TABLES / "cohort_inventory.csv", index=False)
    tests.to_csv(TABLES / "response_tests.csv", index=False)
    pur.to_csv(TABLES / "feature_vs_purity.csv", index=False)
    samples.to_csv(TABLES / "sample_scores.csv", index=False)
    gcov.to_csv(TABLES / "gene_coverage.csv", index=False)
    pd.DataFrame(manifest).to_csv(TABLES / "download_manifest.csv", index=False)

    # summary forest of primary residual AUCs
    if not tests.empty:
        prim_df = tests[tests["primary"] == True]  # noqa: E712
        if len(prim_df):
            fig, ax = plt.subplots(figsize=(8, 4.5))
            y = np.arange(len(prim_df))
            ax.errorbar(prim_df["auc_pos_gt_neg"], y,
                        xerr=[prim_df["auc_pos_gt_neg"] - prim_df["auc_95ci_lo"],
                              prim_df["auc_95ci_hi"] - prim_df["auc_pos_gt_neg"]],
                        fmt="o", color="#333", ecolor="#888", capsize=3)
            ax.axvline(0.5, color="k", ls="--", lw=0.8)
            labels_y = [f"{r.dataset} {r.feature}" for r in prim_df.itertuples()]
            ax.set_yticks(y)
            ax.set_yticklabels(labels_y, fontsize=8)
            ax.set_xlabel("AUC (benefit > no-benefit) after ESTIMATEScore residualization")
            ax.set_title("Primary residualized tests (bootstrap 95% CI)")
            ax.invert_yaxis()
            fig.tight_layout()
            fig.savefig(FIGURES / "primary_residual_auc_forest.png", dpi=150)
            plt.close(fig)

    key = {
        "task": "residualize TACSTD2/CLDN4/TJ_score on ESTIMATE purity, then test DCB/ORR",
        "requested_cohorts": list(LOADERS),
        "tj_5_genes": TJ_5,
        "tj_structural_genes": TJ_STRUCT,
        "purity_primary_covariate": "ESTIMATEScore (Yoshihara 2013; rank-ssGSEA weight 0.25)",
        "purity_formula": "TumorPurity = cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)",
        "residualization": "OLS residual of feature on ESTIMATEScore (intercept + slope)",
        "primary_features": primary_features,
        "primary_measure": "residual_ESTIMATEScore",
        "n_boot_auc": N_BOOT,
        "seed": SEED,
        "cohorts": inv.to_dict(orient="records"),
        "primary_tests": tests[tests.get("primary", False) == True].to_dict(orient="records")
        if not tests.empty and "primary" in tests.columns else [],
    }
    (TABLES / "key_stats.json").write_text(json.dumps(key, indent=2, default=str))

    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 40)
    print("\n=== COHORT INVENTORY (all requested) ===")
    print(inv.to_string(index=False))
    if not tests.empty:
        print("\n=== PRIMARY residual_ESTIMATEScore tests ===")
        cols = ["dataset", "feature", "endpoint", "n_pos", "n_neg",
                "median_pos", "median_neg", "auc_pos_gt_neg", "auc_95ci_lo",
                "auc_95ci_hi", "p_two_sided", "q_BH_within_cohort", "q_BH_primary_residual"]
        print(tests[tests["primary"] == True][cols].to_string(index=False))
        print("\n=== RAW (unadjusted) TACSTD2/CLDN4/TJ_5 ===")
        raw = tests[(tests.measure == "raw") & (tests.feature.isin(primary_features))
                    & (~tests.endpoint.str.contains("MPR"))]
        print(raw[["dataset", "feature", "endpoint", "n_pos", "n_neg",
                   "auc_pos_gt_neg", "p_two_sided"]].to_string(index=False))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
