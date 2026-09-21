#!/usr/bin/env python3
"""Reanalyze public TROP2-ADC and TACSTD2 knockdown/knockout matrices.

Question: in accessible public omics, do CLDN4 and interferon / effector
programs change after TACSTD2 loss or TROP2-ADC exposure, and do those
immune-program shifts move with CLDN4?

Numbers are computed from the deposited matrices. Nothing in the result
tables is copied from papers. Sacituzumab govitecan single-cell studies
on ftp.ebi.ac.uk are catalogued but not reanalyzed when that host cannot
be reached.
"""

from __future__ import annotations

import gzip
import json
import math
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/tmp/trop2_data")
RES = ROOT / "resources"
OUT = ROOT / "results"
PSEUDO = 1.0
EXPR_FLOOR = 1.0  # gene must exceed this in at least one contrast sample

HUMAN_EFFECTOR = [
    "CD8A", "CD8B", "CD3D", "CD3E", "GZMA", "GZMB", "GZMK",
    "PRF1", "NKG7", "GNLY", "IFNG", "CTSW",
]
MOUSE_EFFECTOR = [
    "Cd8a", "Cd8b1", "Cd3d", "Cd3e", "Gzma", "Gzmb", "Gzmk",
    "Prf1", "Nkg7", "Ifng", "Ctsw",
]
HUMAN_LEUKOCYTE = [
    "PTPRC", "CD3D", "CD3E", "CD4", "CD8A", "CD68", "ITGAM",
    "MS4A1", "FCGR3A", "NCR1",
]
MOUSE_LEUKOCYTE = [
    "Ptprc", "Cd3d", "Cd3e", "Cd4", "Cd8a", "Cd68", "Itgam",
    "Ms4a1", "Ncr1", "Adgre1",
]
HUMAN_CONTEXT = [
    "TACSTD2", "CLDN4", "CLDN1", "CLDN3", "CLDN7",
    "EPCAM", "CDH1", "KRT8", "KRT18", "KRT19",
]
MOUSE_CONTEXT = [
    "Tacstd2", "Cldn4", "Cldn1", "Cldn3", "Cldn7",
    "Epcam", "Cdh1", "Krt8", "Krt18", "Krt19",
]

DOWNLOADS = {
    "GSE245459_fpkm.anno.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245459/suppl/GSE245459_fpkm.anno.txt.gz",
    "GSE334497_normalized_counts.csv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE334nnn/GSE334497/suppl/GSE334497_normalized_counts.csv.gz",
    "GSE304294_gene_fpkm.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE304nnn/GSE304294/suppl/GSE304294_gene_fpkm.txt.gz",
    "GSE312098_gene_fpkm.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE312nnn/GSE312098/suppl/GSE312098_gene_fpkm.txt.gz",
    "GSE311016_gene_fpkm.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311016/suppl/GSE311016_gene_fpkm.txt.gz",
    "GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289287/suppl/GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz",
    "GSE289287_DESeq2-DSG2KO_tumors_vs_WT.tsv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289287/suppl/GSE289287_DESeq2-DSG2KO_tumors_vs_WT.tsv.gz",
}


def ensure_downloads() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    for name, url in DOWNLOADS.items():
        dest = DATA / name
        if dest.exists() and dest.stat().st_size > 0:
            continue
        print(f"downloading {name}")
        urllib.request.urlretrieve(url, dest)


def read_table(path: Path, encoding: str | None = None) -> pd.DataFrame:
    if encoding is None:
        with gzip.open(path, "rb") as handle:
            magic = handle.read(4)
        encoding = "utf-16" if magic[:2] in (b"\xff\xfe", b"\xfe\xff") else "utf-8"
    sep = "," if path.name.endswith(".csv.gz") else "\t"
    return pd.read_csv(path, sep=sep, encoding=encoding, low_memory=False)


def load_geneset(path: Path) -> list[str]:
    genes = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(">") or line.startswith("HALLMARK_"):
            continue
        genes.append(line)
    if len(genes) < 20:
        raise RuntimeError(f"gene set {path} looks truncated ({len(genes)} genes)")
    return genes


def collapse_max_mean(frame: pd.DataFrame, symbol_col: str, sample_cols: list[str]) -> pd.DataFrame:
    work = frame[[symbol_col, *sample_cols]].copy()
    work[symbol_col] = work[symbol_col].astype(str)
    work = work[work[symbol_col].notna() & (work[symbol_col] != "") & (work[symbol_col] != "nan")]
    for col in sample_cols:
        work[col] = pd.to_numeric(work[col], errors="coerce")
    work["_mean"] = work[sample_cols].mean(axis=1)
    work = work.sort_values("_mean", ascending=False)
    work = work.drop_duplicates(symbol_col, keep="first")
    out = work.set_index(symbol_col)[sample_cols].astype(float)
    out.index.name = "symbol"
    return out


def ensembl_symbol_map(symbols: list[str]) -> dict[str, str]:
    """Mouse symbol -> Ensembl gene id. Missing symbols are omitted."""
    cache_path = OUT / "mouse_symbol_to_ensembl.tsv"
    cached: dict[str, str] = {}
    if cache_path.exists():
        prev = pd.read_csv(cache_path, sep="\t")
        cached = dict(zip(prev["symbol"], prev["ensembl_id"]))
    needed = [s for s in symbols if s not in cached]
    for start in range(0, len(needed), 200):
        chunk = needed[start:start + 200]
        payload = json.dumps({"symbols": chunk}).encode()
        req = urllib.request.Request(
            "https://rest.ensembl.org/lookup/symbol/mus_musculus",
            data=payload,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=120) as resp:
            data = json.load(resp)
        for symbol, rec in data.items():
            if isinstance(rec, dict) and rec.get("id"):
                cached[symbol] = rec["id"].split(".")[0]
    rows = [{"symbol": s, "ensembl_id": cached.get(s, "")} for s in sorted(set(symbols))]
    pd.DataFrame(rows).to_csv(cache_path, sep="\t", index=False)
    return {s: i for s, i in cached.items() if i}


def mouse_expr_by_symbol(ensembl_expr: pd.DataFrame, symbols: list[str]) -> pd.DataFrame:
    mapping = ensembl_symbol_map(symbols)
    rows = []
    index = []
    missing = []
    for symbol in symbols:
        ensg = mapping.get(symbol)
        if ensg is None or ensg not in ensembl_expr.index:
            missing.append(symbol)
            continue
        rows.append(ensembl_expr.loc[ensg])
        index.append(symbol)
    if missing:
        print(f"  mouse symbols absent from matrix or Ensembl: {missing}")
    if not rows:
        return pd.DataFrame(columns=ensembl_expr.columns)
    out = pd.DataFrame(rows, index=index)
    out.index.name = "symbol"
    # If two symbols collided on one id, keep the first.
    out = out[~out.index.duplicated(keep="first")]
    return out.astype(float)


def log_mean(values: np.ndarray) -> float:
    return float(np.mean(np.log2(np.asarray(values, dtype=float) + PSEUDO)))


def safe_mw(treat: np.ndarray, ctrl: np.ndarray) -> tuple[float, str]:
    treat = np.asarray(treat, dtype=float)
    ctrl = np.asarray(ctrl, dtype=float)
    if len(treat) == 0 or len(ctrl) == 0:
        return math.nan, "empty"
    try:
        res = stats.mannwhitneyu(treat, ctrl, alternative="two-sided", method="exact")
        return float(res.pvalue), "exact"
    except ValueError:
        res = stats.mannwhitneyu(treat, ctrl, alternative="two-sided", method="asymptotic")
        return float(res.pvalue), "asymptotic"


def safe_wilcoxon(deltas: np.ndarray) -> tuple[float, str, int]:
    deltas = np.asarray(deltas, dtype=float)
    nonzero = deltas[np.abs(deltas) > 0]
    if len(nonzero) == 0:
        return math.nan, "all_zero", 0
    try:
        res = stats.wilcoxon(deltas, alternative="two-sided", method="exact", zero_method="wilcox")
        return float(res.pvalue), "exact_drop_zero", int(len(nonzero))
    except ValueError:
        res = stats.wilcoxon(deltas, alternative="two-sided", method="asymptotic", zero_method="wilcox")
        return float(res.pvalue), "asymptotic_drop_zero", int(len(nonzero))


def _rank_corr(rx: np.ndarray, ry: np.ndarray) -> float:
    rc = rx - rx.mean()
    yc = ry - ry.mean()
    denom = math.sqrt(float(np.sum(rc ** 2) * np.sum(yc ** 2)))
    if denom == 0:
        return math.nan
    return float(np.dot(rc, yc) / denom)


def safe_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int, str]:
    """Spearman rho with a two-sided p-value.

    SciPy's t approximation goes to ~0 when |rho| is 1. For n <= 8 the p-value
    is the exact permutation test. Larger n keeps the t approximation, except
    |rho| == 1, which is reported as 2/n!.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    n = int(x.size)
    if n < 3 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return math.nan, math.nan, n, "undefined"
    rx = stats.rankdata(x).astype(float)
    ry = stats.rankdata(y).astype(float)
    rho = _rank_corr(rx, ry)
    if n <= 8:
        perms = np.array(list(__import__("itertools").permutations(range(n))))
        ry_perm = ry[perms]
        rc = rx - rx.mean()
        yc = ry_perm - ry_perm.mean(axis=1, keepdims=True)
        denom = np.sqrt(np.sum(rc ** 2) * np.sum(yc ** 2, axis=1))
        rhos = (yc @ rc) / denom
        p = float(np.mean(np.abs(rhos) + 1e-12 >= abs(rho)))
        return rho, p, n, "exact_permutation"
    if abs(rho) >= 1 - 1e-12:
        return rho, min(1.0, 2 / math.factorial(n)), n, "exact_perfect"
    _rho, p = stats.spearmanr(x, y)
    return float(_rho), float(p), n, "student_t"


def score_samples(expr: pd.DataFrame, genes: list[str], samples: list[str]) -> tuple[pd.Series, list[str], int]:
    """Mean log2(expr+1) over genes present and above the expression floor."""
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=samples), [], 0
    sub = expr.loc[present, samples]
    keep = [g for g in sub.index if float(sub.loc[g].max()) > EXPR_FLOOR]
    if len(keep) == 0:
        return pd.Series(np.nan, index=samples), [], len(present)
    logged = np.log2(sub.loc[keep] + PSEUDO)
    return logged.mean(axis=0), keep, len(present)


def gene_stats(expr: pd.DataFrame, symbol: str, treat: list[str], ctrl: list[str], paired: list[tuple[str, str]] | None) -> dict:
    if symbol not in expr.index:
        return {
            "symbol": symbol,
            "present": False,
            "mean_ctrl": math.nan,
            "mean_treat": math.nan,
            "delta_mean_log2p1": math.nan,
            "mw_p": math.nan,
            "mw_method": "absent",
            "n_ctrl": len(ctrl),
            "n_treat": len(treat),
        }
    cvals = expr.loc[symbol, ctrl].to_numpy(dtype=float)
    tvals = expr.loc[symbol, treat].to_numpy(dtype=float)
    p, method = safe_mw(np.log2(tvals + PSEUDO), np.log2(cvals + PSEUDO))
    row = {
        "symbol": symbol,
        "present": True,
        "mean_ctrl": float(np.mean(cvals)),
        "mean_treat": float(np.mean(tvals)),
        "median_ctrl": float(np.median(cvals)),
        "median_treat": float(np.median(tvals)),
        "delta_mean_log2p1": log_mean(tvals) - log_mean(cvals),
        "mw_p": p,
        "mw_method": method,
        "n_ctrl": len(ctrl),
        "n_treat": len(treat),
    }
    if paired:
        deltas = []
        for t_col, c_col in paired:
            deltas.append(
                math.log2(float(expr.loc[symbol, t_col]) + PSEUDO)
                - math.log2(float(expr.loc[symbol, c_col]) + PSEUDO)
            )
        wp, wmethod, nnz = safe_wilcoxon(np.asarray(deltas))
        row["paired_mean_delta_log2p1"] = float(np.mean(deltas))
        row["paired_median_delta_log2p1"] = float(np.median(deltas))
        row["wilcoxon_p"] = wp
        row["wilcoxon_method"] = wmethod
        row["wilcoxon_n_nonzero"] = nnz
        row["n_pairs_same_sign_as_mean"] = int(np.sum(np.sign(deltas) == np.sign(np.mean(deltas))))
    return row


def contrast_record(
    expr: pd.DataFrame,
    *,
    dataset: str,
    contrast: str,
    role: str,
    compartment: str,
    species: str,
    treat: list[str],
    ctrl: list[str],
    sets: dict[str, list[str]],
    context_symbols: list[str],
    target_symbol: str,
    cldn_symbol: str,
    paired: list[tuple[str, str]] | None = None,
    note: str = "",
) -> tuple[dict, list[dict], list[dict]]:
    samples = ctrl + treat
    gene_rows = []
    for symbol in context_symbols:
        rec = gene_stats(expr, symbol, treat, ctrl, paired)
        rec.update({"dataset": dataset, "contrast": contrast, "role": role})
        gene_rows.append(rec)

    set_rows = []
    scores = {}
    used_genes = {}
    for set_name, genes in sets.items():
        score, keep, n_present = score_samples(expr, genes, samples)
        scores[set_name] = score
        used_genes[set_name] = keep
        n_used = len(keep)
        cvals = score.loc[ctrl].to_numpy(dtype=float)
        tvals = score.loc[treat].to_numpy(dtype=float)
        if n_used == 0 or np.all(~np.isfinite(cvals)) or np.all(~np.isfinite(tvals)):
            p, method = math.nan, "no_genes_above_floor"
        else:
            p, method = safe_mw(tvals, cvals)
        row = {
            "dataset": dataset,
            "contrast": contrast,
            "role": role,
            "compartment": compartment,
            "species": species,
            "score": set_name,
            "n_genes_in_set": len(genes),
            "n_genes_present": n_present,
            "n_genes_used": n_used,
            "genes_used": ";".join(keep) if n_used <= 20 else "",
            "interpretable_program": bool(n_used >= 5),
            "expr_floor": EXPR_FLOOR,
            "mean_score_ctrl": float(np.nanmean(cvals)) if np.isfinite(cvals).any() else math.nan,
            "mean_score_treat": float(np.nanmean(tvals)) if np.isfinite(tvals).any() else math.nan,
            "delta_mean_score": (
                float(np.nanmean(tvals) - np.nanmean(cvals))
                if np.isfinite(cvals).any() and np.isfinite(tvals).any() else math.nan
            ),
            "mw_p": p,
            "mw_method": method,
            "n_ctrl": len(ctrl),
            "n_treat": len(treat),
            "note": note,
        }
        if paired and n_used > 0:
            deltas = [float(score.loc[t] - score.loc[c]) for t, c in paired]
            wp, wmethod, nnz = safe_wilcoxon(np.asarray(deltas))
            row["paired_mean_delta"] = float(np.mean(deltas))
            row["wilcoxon_p"] = wp
            row["wilcoxon_method"] = wmethod
            row["wilcoxon_n_nonzero"] = nnz
        set_rows.append(row)

    # Sample-level tracking inside this contrast only.
    cldn = None
    if cldn_symbol in expr.index:
        cldn = np.log2(expr.loc[cldn_symbol, samples].to_numpy(dtype=float) + PSEUDO)
    corr_rows = []
    for set_name, score in scores.items():
        if cldn is None or used_genes[set_name] == []:
            rho, p, n, pmethod = math.nan, math.nan, 0 if cldn is None else len(samples), "undefined"
        else:
            rho, p, n, pmethod = safe_spearman(cldn, score.loc[samples].to_numpy(dtype=float))
        corr_rows.append({
            "dataset": dataset,
            "contrast": contrast,
            "role": role,
            "compartment": compartment,
            "x": f"log2({cldn_symbol}+1)",
            "y": set_name,
            "scope": "contrast_samples",
            "spearman_rho": rho,
            "spearman_p": p,
            "spearman_p_method": pmethod,
            "n": n,
            "n_genes_used": len(used_genes[set_name]),
            "interpretable_program": len(used_genes[set_name]) >= 5,
        })
    if paired and cldn_symbol in expr.index:
        cldn_d = np.asarray([
            math.log2(float(expr.loc[cldn_symbol, t]) + PSEUDO)
            - math.log2(float(expr.loc[cldn_symbol, c]) + PSEUDO)
            for t, c in paired
        ])
        for set_name, score in scores.items():
            if score.isna().all():
                rho, p, n, pmethod = math.nan, math.nan, len(paired), "undefined"
                deltas = np.full(len(paired), np.nan)
            else:
                deltas = np.asarray([float(score.loc[t] - score.loc[c]) for t, c in paired])
                rho, p, n, pmethod = safe_spearman(cldn_d, deltas)
            corr_rows.append({
                "dataset": dataset,
                "contrast": contrast,
                "role": role,
                "compartment": compartment,
                "x": f"paired_delta_log2({cldn_symbol}+1)",
                "y": f"paired_delta_{set_name}",
                "scope": "paired_deltas",
                "spearman_rho": rho,
                "spearman_p": p,
                "spearman_p_method": pmethod,
                "n": n,
                "n_genes_used": len(used_genes[set_name]),
                "interpretable_program": len(used_genes[set_name]) >= 5,
                "paired_score_deltas": ";".join(f"{d:.6g}" for d in deltas),
                "paired_cldn_deltas": ";".join(f"{d:.6g}" for d in cldn_d),
                "paired_labels": ";".join(f"{t}_minus_{c}" for t, c in paired),
            })

    target = next(r for r in gene_rows if r["symbol"] == target_symbol)
    cldn_row = next(r for r in gene_rows if r["symbol"] == cldn_symbol)
    summary = {
        "dataset": dataset,
        "contrast": contrast,
        "role": role,
        "compartment": compartment,
        "species": species,
        "n_ctrl": len(ctrl),
        "n_treat": len(treat),
        "target_symbol": target_symbol,
        "target_mean_ctrl": target["mean_ctrl"],
        "target_mean_treat": target["mean_treat"],
        "target_delta_log2p1": target["delta_mean_log2p1"],
        "target_mw_p": target["mw_p"],
        "cldn_symbol": cldn_symbol,
        "cldn_mean_ctrl": cldn_row["mean_ctrl"],
        "cldn_mean_treat": cldn_row["mean_treat"],
        "cldn_delta_log2p1": cldn_row["delta_mean_log2p1"],
        "cldn_mw_p": cldn_row["mw_p"],
        "note": note,
    }
    for row in set_rows:
        summary[f"{row['score']}_delta"] = row["delta_mean_score"]
        summary[f"{row['score']}_mw_p"] = row["mw_p"]
        summary[f"{row['score']}_n_genes_used"] = row["n_genes_used"]
        if "wilcoxon_p" in row:
            summary[f"{row['score']}_wilcoxon_p"] = row["wilcoxon_p"]
            summary[f"{row['score']}_paired_mean_delta"] = row["paired_mean_delta"]
    if "wilcoxon_p" in cldn_row:
        summary["cldn_wilcoxon_p"] = cldn_row["wilcoxon_p"]
        summary["cldn_paired_mean_delta"] = cldn_row["paired_mean_delta_log2p1"]
        summary["target_wilcoxon_p"] = target.get("wilcoxon_p", math.nan)
        summary["target_paired_mean_delta"] = target.get("paired_mean_delta_log2p1", math.nan)
    return summary, gene_rows, set_rows, corr_rows


def load_symbol_matrix_from_fpkm(path: Path, sample_cols: list[str], encoding: str | None = None) -> pd.DataFrame:
    raw = read_table(path, encoding=encoding)
    symbol_col = "gene_name" if "gene_name" in raw.columns else "GeneName"
    missing = [c for c in sample_cols if c not in raw.columns]
    if missing:
        raise RuntimeError(f"{path.name} missing columns {missing}; have {list(raw.columns)[:20]}")
    return collapse_max_mean(raw, symbol_col, sample_cols)


def load_4t1() -> pd.DataFrame:
    raw = read_table(DATA / "GSE334497_normalized_counts.csv.gz")
    # First column is an unnamed Ensembl index when read with a header.
    if "ENSMUSG" not in str(raw.columns[0]) and raw.columns[0] in ("Unnamed: 0", ""):
        raw = raw.rename(columns={raw.columns[0]: "ensembl"})
    elif str(raw.iloc[0, 0]).startswith("ENSMUSG") is False and raw.columns[0].startswith("ENSMUSG"):
        raw = raw.reset_index().rename(columns={"index": "ensembl"})
    if "ensembl" not in raw.columns:
        # read_csv without index: the file's first column header is empty and becomes Unnamed: 0
        raw = raw.rename(columns={raw.columns[0]: "ensembl"})
    sample_cols = [c for c in raw.columns if c != "ensembl"]
    raw["ensembl"] = raw["ensembl"].astype(str).str.split(".").str[0]
    raw = raw.drop_duplicates("ensembl", keep="first")
    expr = raw.set_index("ensembl")[sample_cols].apply(pd.to_numeric, errors="coerce")
    symbols = sorted(set(MOUSE_CONTEXT + MOUSE_EFFECTOR + MOUSE_LEUKOCYTE
                          + load_geneset(RES / "hallmark_interferon_alpha_mouse.txt")
                          + load_geneset(RES / "hallmark_interferon_gamma_mouse.txt")))
    return mouse_expr_by_symbol(expr, symbols)


def normcount_matrix(path: Path, animals: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = read_table(path)
    norm_cols = [c for c in raw.columns if c.endswith("RNA_normCounts")]
    if not norm_cols:
        raise RuntimeError(f"no normCounts in {path}")
    # Keep columns whose animal id is in the requested map.
    keep = []
    labels = []
    for col in norm_cols:
        animal = col.split(".")[1]
        if animal in animals:
            keep.append(col)
            labels.append(f"{animals[animal]}_{animal}")
    work = raw[["Feature_name", "biotype", "log2FoldChange", *keep]].copy()
    work = work[work["biotype"] == "protein_coding"]
    renamed = work.rename(columns=dict(zip(keep, labels)))
    expr = collapse_max_mean(renamed, "Feature_name", labels)
    deposited = (
        work.sort_values("log2FoldChange", key=lambda s: s.abs(), ascending=False)
        .drop_duplicates("Feature_name", keep="first")
        .set_index("Feature_name")["log2FoldChange"]
    )
    return expr, deposited


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ensure_downloads()
    ifna_h = load_geneset(RES / "hallmark_interferon_alpha_human.txt")
    ifng_h = load_geneset(RES / "hallmark_interferon_gamma_human.txt")
    ifna_m = load_geneset(RES / "hallmark_interferon_alpha_mouse.txt")
    ifng_m = load_geneset(RES / "hallmark_interferon_gamma_mouse.txt")
    human_sets = {
        "ifn_alpha": ifna_h,
        "ifn_gamma": ifng_h,
        "effector": HUMAN_EFFECTOR,
        "leukocyte": HUMAN_LEUKOCYTE,
    }
    mouse_sets = {
        "ifn_alpha": ifna_m,
        "ifn_gamma": ifng_m,
        "effector": MOUSE_EFFECTOR,
        "leukocyte": MOUSE_LEUKOCYTE,
    }

    summaries = []
    gene_rows = []
    set_rows = []
    corr_rows = []
    deseq_rows = []

    def add(result):
        summary, genes, sets, corrs = result
        summaries.append(summary)
        gene_rows.extend(genes)
        set_rows.extend(sets)
        corr_rows.extend(corrs)

    # GSE245459 SKOV3 TACSTD2 shRNA.
    sk = load_symbol_matrix_from_fpkm(
        DATA / "GSE245459_fpkm.anno.txt.gz",
        ["shNC1", "shNC2", "shNC3", "sh1", "sh2", "sh3",
         "shNCDDP1", "shNCDDP2", "shNCDDP3", "shDDP1", "shDDP2", "shDDP3"],
    )
    add(contrast_record(
        sk, dataset="GSE245459", contrast="shTACSTD2_vs_shNC",
        role="primary_kd", compartment="in_vitro_cell_line", species="human",
        treat=["sh1", "sh2", "sh3"], ctrl=["shNC1", "shNC2", "shNC3"],
        sets=human_sets, context_symbols=HUMAN_CONTEXT,
        target_symbol="TACSTD2", cldn_symbol="CLDN4",
        note="SKOV3 shRNA called sh1-3 versus shNC. Series text says knockout; sample names are shRNA. No cisplatin.",
    ))
    add(contrast_record(
        sk, dataset="GSE245459", contrast="shTACSTD2_DDP_vs_shNC_DDP",
        role="secondary_kd_on_cisplatin", compartment="in_vitro_cell_line", species="human",
        treat=["shDDP1", "shDDP2", "shDDP3"], ctrl=["shNCDDP1", "shNCDDP2", "shNCDDP3"],
        sets=human_sets, context_symbols=HUMAN_CONTEXT,
        target_symbol="TACSTD2", cldn_symbol="CLDN4",
        note="Same shRNA contrast on a cisplatin background. Not independent of the no-cisplatin contrast.",
    ))

    # GSE334497 4T1 Trop2 KO tumors in BALB/c. Library names are GEO sample descriptions.
    print("mapping 4T1 Ensembl ids")
    t4 = load_4t1()
    ko_cols = ["KO162", "KO164", "KO165", "KO172", "RESUB-KO163R"]
    wt_cols = ["RESUB-171R", "RESUB-170R", "RESUB-169R", "RESUB-168R", "control170"]
    add(contrast_record(
        t4, dataset="GSE334497", contrast="Trop2KO_vs_WT",
        role="primary_ko", compartment="bulk_tumor_immunocompetent_mouse", species="mouse",
        treat=ko_cols, ctrl=wt_cols,
        sets=mouse_sets, context_symbols=MOUSE_CONTEXT,
        target_symbol="Tacstd2", cldn_symbol="Cldn4",
        note="4T1 tumors in BALB/c. Column IDs are GEO library names: KO162/164/165/172/RESUB-KO163R versus RESUB-171R/170R/169R/168R/control170.",
    ))

    # GSE289287 T-47D Trop2 KO xenografts. Norm counts live in the DESeq2 table.
    trop_animals = {
        "2808": "WT", "2810": "WT", "2812": "WT",
        "2807": "KO", "2815": "KO", "2817": "KO", "2818": "KO",
    }
    trop_expr, trop_lfc = normcount_matrix(
        DATA / "GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz", trop_animals,
    )
    trop_ctrl = [c for c in trop_expr.columns if c.startswith("WT_")]
    trop_treat = [c for c in trop_expr.columns if c.startswith("KO_")]
    add(contrast_record(
        trop_expr, dataset="GSE289287", contrast="Trop2KO_vs_WT_xenograft",
        role="primary_ko", compartment="bulk_xenograft_NRG_human_genes", species="human",
        treat=trop_treat, ctrl=trop_ctrl,
        sets=human_sets, context_symbols=HUMAN_CONTEXT,
        target_symbol="TACSTD2", cldn_symbol="CLDN4",
        note="T-47D mammary-fat-pad tumors in NRG mice. Recomputed from DESeq2 normalized counts in the deposited table, protein-coding rows only. NRG mice lack adaptive immunity; human gene matrix does not measure mouse infiltrate.",
    ))
    for symbol in ["TACSTD2", "CLDN4", "CLDN7", "CLDN1"]:
        if symbol in trop_lfc.index:
            deseq_rows.append({
                "dataset": "GSE289287",
                "contrast": "Trop2KO_vs_WT_xenograft",
                "symbol": symbol,
                "deposited_deseq2_log2FoldChange": float(trop_lfc.loc[symbol]),
            })

    dsg_animals = {
        "2808": "WT", "2810": "WT", "2812": "WT",
        "4111": "KO", "4113": "KO", "4197": "KO",
        "4121": "KO", "4133": "KO", "4134": "KO", "4135": "KO",
    }
    dsg_expr, dsg_lfc = normcount_matrix(
        DATA / "GSE289287_DESeq2-DSG2KO_tumors_vs_WT.tsv.gz", dsg_animals,
    )
    dsg_ctrl = [c for c in dsg_expr.columns if c.startswith("WT_")]
    dsg_treat = [c for c in dsg_expr.columns if c.startswith("KO_")]
    add(contrast_record(
        dsg_expr, dataset="GSE289287", contrast="DSG2KO_vs_WT_xenograft",
        role="comparator_not_TACSTD2", compartment="bulk_xenograft_NRG_human_genes", species="human",
        treat=dsg_treat, ctrl=dsg_ctrl,
        sets=human_sets, context_symbols=HUMAN_CONTEXT,
        target_symbol="TACSTD2", cldn_symbol="CLDN4",
        note="Same study, DSG2 knockout rather than Trop2. Comparator only. Normalized in a separate DESeq2 model from the Trop2 table, so counts are not pooled with Trop2 KO samples.",
    ))
    for symbol in ["TACSTD2", "CLDN4"]:
        if symbol in dsg_lfc.index:
            deseq_rows.append({
                "dataset": "GSE289287",
                "contrast": "DSG2KO_vs_WT_xenograft",
                "symbol": symbol,
                "deposited_deseq2_log2FoldChange": float(dsg_lfc.loc[symbol]),
            })

    # GSE304294 KYSE30 IMMU132.
    ky = load_symbol_matrix_from_fpkm(
        DATA / "GSE304294_gene_fpkm.txt.gz",
        ["OX1_1", "OX1_2", "OX1_3", "OX2_1", "OX2_2", "OX3_1", "OX3_2", "OX3_3", "OX4_1", "OX4_2", "OX4_3"],
    )
    add(contrast_record(
        ky, dataset="GSE304294", contrast="IMMU132_vs_control",
        role="primary_adc", compartment="in_vitro_cell_line", species="human",
        treat=["OX2_1", "OX2_2"], ctrl=["OX1_1", "OX1_2", "OX1_3"],
        sets=human_sets, context_symbols=HUMAN_CONTEXT,
        target_symbol="TACSTD2", cldn_symbol="CLDN4",
        note="KYSE30, 1 day. Library names OX2 are IMMU132 (n=2) and OX1 are control (n=3). ADC exposure is not a genetic knockout, so TACSTD2 mRNA need not fall.",
    ))
    add(contrast_record(
        ky, dataset="GSE304294", contrast="IACS_vs_control",
        role="same_experiment_not_adc", compartment="in_vitro_cell_line", species="human",
        treat=["OX3_1", "OX3_2", "OX3_3"], ctrl=["OX1_1", "OX1_2", "OX1_3"],
        sets=human_sets, context_symbols=HUMAN_CONTEXT,
        target_symbol="TACSTD2", cldn_symbol="CLDN4",
        note="OXPHOS inhibitor IACS-010759 alone. Not a TROP2 perturbation.",
    ))
    add(contrast_record(
        ky, dataset="GSE304294", contrast="IMMU132_plus_IACS_vs_control",
        role="secondary_adc_combination", compartment="in_vitro_cell_line", species="human",
        treat=["OX4_1", "OX4_2", "OX4_3"], ctrl=["OX1_1", "OX1_2", "OX1_3"],
        sets=human_sets, context_symbols=HUMAN_CONTEXT,
        target_symbol="TACSTD2", cldn_symbol="CLDN4",
        note="IMMU132 plus IACS-010759. Not attributable to the ADC alone.",
    ))

    # GSE312098 CX-1 IMMU132.
    cx = load_symbol_matrix_from_fpkm(
        DATA / "GSE312098_gene_fpkm.txt.gz",
        [f"X_{i}" for i in range(1, 13)],
        encoding="utf-16",
    )
    add(contrast_record(
        cx, dataset="GSE312098", contrast="IMMU132_vs_control",
        role="primary_adc", compartment="in_vitro_cell_line", species="human",
        treat=["X_4", "X_5", "X_6"], ctrl=["X_1", "X_2", "X_3"],
        sets=human_sets, context_symbols=HUMAN_CONTEXT,
        target_symbol="TACSTD2", cldn_symbol="CLDN4",
        note="CX-1 cells, 48 h, 3 ug/ml IMMU132. Library names from GEO: X_1-3 control, X_4-6 IMMU132.",
    ))
    add(contrast_record(
        cx, dataset="GSE312098", contrast="GSK2606414_vs_control",
        role="same_experiment_not_adc", compartment="in_vitro_cell_line", species="human",
        treat=["X_7", "X_8", "X_9"], ctrl=["X_1", "X_2", "X_3"],
        sets=human_sets, context_symbols=HUMAN_CONTEXT,
        target_symbol="TACSTD2", cldn_symbol="CLDN4",
        note="PERK inhibitor GSK2606414 alone.",
    ))
    add(contrast_record(
        cx, dataset="GSE312098", contrast="IMMU132_plus_GSK_vs_control",
        role="secondary_adc_combination", compartment="in_vitro_cell_line", species="human",
        treat=["X_10", "X_11", "X_12"], ctrl=["X_1", "X_2", "X_3"],
        sets=human_sets, context_symbols=HUMAN_CONTEXT,
        target_symbol="TACSTD2", cldn_symbol="CLDN4",
        note="IMMU132 plus GSK2606414. Not attributable to the ADC alone.",
    ))

    # GSE311016 paired CRC PDX.
    pdx_cols = ["C_114", "T_114", "C_36", "T_36", "C_82", "T_82", "C_83", "T_83", "C_196", "T_196"]
    pdx = load_symbol_matrix_from_fpkm(
        DATA / "GSE311016_gene_fpkm.txt.gz", pdx_cols, encoding="utf-16",
    )
    pairs = [("T_114", "C_114"), ("T_36", "C_36"), ("T_82", "C_82"), ("T_83", "C_83"), ("T_196", "C_196")]
    add(contrast_record(
        pdx, dataset="GSE311016", contrast="IMMU132_vs_paired_control",
        role="primary_adc", compartment="bulk_PDX_human_genes", species="human",
        treat=[t for t, _ in pairs], ctrl=[c for _, c in pairs],
        sets=human_sets, context_symbols=HUMAN_CONTEXT,
        target_symbol="TACSTD2", cldn_symbol="CLDN4",
        paired=pairs,
        note="Five CRC PDX models, control versus IMMU132 at day 29. Column prefix C/T matches GEO treatment labels. Matrix is human genes, so mouse stromal immune cells are not in the table. Host strain is not stated in the GEO fields read here.",
    ))

    summary_df = pd.DataFrame(summaries)
    gene_df = pd.DataFrame(gene_rows)
    set_df = pd.DataFrame(set_rows)
    corr_df = pd.DataFrame(corr_rows)
    deseq_df = pd.DataFrame(deseq_rows)

    # Attach recomputed deltas next to deposited DESeq2 log2FC.
    if not deseq_df.empty:
        key = gene_df[["dataset", "contrast", "symbol", "delta_mean_log2p1", "mean_ctrl", "mean_treat", "mw_p"]]
        deseq_df = deseq_df.merge(key, on=["dataset", "contrast", "symbol"], how="left")

    # Cross-contrast sign alignment on prespecified primary contrasts.
    primary = summary_df[summary_df["role"].isin(["primary_kd", "primary_ko", "primary_adc"])].copy()
    align_rows = []
    for _, row in primary.iterrows():
        for score in ("ifn_gamma", "ifn_alpha", "effector", "leukocyte"):
            d_cl = row["cldn_delta_log2p1"]
            d_sc = row[f"{score}_delta"]
            n_used = row[f"{score}_n_genes_used"]
            interpretable = bool(pd.notna(n_used) and n_used >= 5 and row["n_ctrl"] >= 3 and row["n_treat"] >= 3)
            if pd.isna(d_cl) or pd.isna(d_sc) or d_cl == 0 or d_sc == 0:
                relation = "undefined"
            elif np.sign(d_cl) == np.sign(d_sc):
                relation = "same_direction"
            else:
                relation = "opposite_direction"
            align_rows.append({
                "dataset": row["dataset"],
                "contrast": row["contrast"],
                "role": row["role"],
                "compartment": row["compartment"],
                "n_ctrl": row["n_ctrl"],
                "n_treat": row["n_treat"],
                "cldn_delta_log2p1": d_cl,
                "score": score,
                "score_delta": d_sc,
                "score_n_genes_used": n_used,
                "relation": relation,
                "interpretable_program": interpretable,
            })
    align_df = pd.DataFrame(align_rows)

    summary_df.to_csv(OUT / "contrasts.tsv", sep="\t", index=False)
    gene_df.to_csv(OUT / "context_genes.tsv", sep="\t", index=False)
    set_df.to_csv(OUT / "scores.tsv", sep="\t", index=False)
    corr_df.to_csv(OUT / "correlations.tsv", sep="\t", index=False)
    align_df.to_csv(OUT / "shift_alignment.tsv", sep="\t", index=False)
    deseq_df.to_csv(OUT / "deseq2_check.tsv", sep="\t", index=False)

    # Compact printed view for the log.
    cols = [
        "dataset", "contrast", "role", "n_ctrl", "n_treat",
        "target_delta_log2p1", "target_mw_p",
        "cldn_delta_log2p1", "cldn_mw_p", "cldn_mean_ctrl", "cldn_mean_treat",
        "ifn_gamma_delta", "ifn_gamma_mw_p", "ifn_gamma_n_genes_used",
        "ifn_alpha_delta", "effector_delta", "leukocyte_delta",
    ]
    print(summary_df[cols].to_string(index=False))
    print("\nCORRELATIONS")
    print(corr_df.to_string(index=False))
    print("\nALIGNMENT")
    print(align_df.to_string(index=False))
    if not deseq_df.empty:
        print("\nDESEQ2 CHECK")
        print(deseq_df.to_string(index=False))


if __name__ == "__main__":
    main()
