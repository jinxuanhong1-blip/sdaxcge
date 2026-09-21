#!/usr/bin/env python3
"""Concordant-4 malignant ELF3–TACSTD2–CLDN4–CLDN7 coexpression and regulon.

The unit is the locked patient / donor / sample (n=65). Cell counts are not n.
"""

from __future__ import annotations

import argparse
import gzip
import os
import tarfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
LOCKED_PATH = ROOT / "data" / "locked_units.tsv"
DOROTHEA_PATH = ROOT / "resources" / "dorothea_hs_ABC.tsv"
CACHE = ROOT / "results" / "cache"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"

COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
FOCUS = ["ELF3", "GRHL1", "GRHL2", "GRHL3", "TACSTD2", "CLDN4", "CLDN7"]
MODULE4 = ["ELF3", "TACSTD2", "CLDN4", "CLDN7"]
MODULE3 = ["ELF3", "TACSTD2", "CLDN7"]
WATCH = [
    "ELF3", "GRHL1", "GRHL2", "GRHL3", "TACSTD2", "CLDN4", "CLDN7",
    "EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "CDH1", "PTPRC",
]
PAIRS = [
    ("ELF3", "TACSTD2", "primary"),
    ("ELF3", "CLDN4", "primary"),
    ("ELF3", "CLDN7", "primary"),
    ("TACSTD2", "CLDN4", "primary"),
    ("TACSTD2", "CLDN7", "primary"),
    ("CLDN4", "CLDN7", "primary"),
    ("ELF3", "GRHL2", "primary"),
    ("GRHL2", "CLDN4", "primary"),
    ("GRHL2", "TACSTD2", "primary"),
    ("GRHL2", "CLDN7", "primary"),
    ("GRHL1", "ELF3", "secondary"),
    ("GRHL1", "CLDN4", "secondary"),
    ("GRHL3", "ELF3", "secondary"),
    ("GRHL3", "CLDN4", "secondary"),
]
GSE123902_FILES = {
    "LX255B": "GSM3516668_MSK_LX255B_METASTASIS_dense.csv.gz",
    "LX653": "GSM3516662_MSK_LX653_PRIMARY_TUMOUR_dense.csv.gz",
    "LX661": "GSM3516663_MSK_LX661_PRIMARY_TUMOUR_dense.csv.gz",
    "LX666": "GSM3516664_MSK_LX666_METASTASIS_dense.csv.gz",
    "LX675": "GSM3516665_MSK_LX675_PRIMARY_TUMOUR_dense.csv.gz",
    "LX676": "GSM3516667_MSK_LX676_PRIMARY_TUMOUR_dense.csv.gz",
    "LX679": "GSM3516669_MSK_LX679_PRIMARY_TUMOUR_dense.csv.gz",
    "LX680": "GSM3516670_MSK_LX680_PRIMARY_TUMOUR_dense.csv.gz",
    "LX681": "GSM3516671_MSK_LX681_METASTASIS_dense.csv.gz",
    "LX682": "GSM3516672_MSK_LX682_PRIMARY_TUMOUR_dense.csv.gz",
    "LX684": "GSM3516674_MSK_LX684_PRIMARY_TUMOUR_dense.csv.gz",
    "LX699": "GSM3516677_MSK_LX699_METASTASIS_dense.csv.gz",
    "LX701": "GSM3516678_MSK_LX701_METASTASIS_dense.csv.gz",
}
COLORS = {
    "GSE123902": "#1b4f72",
    "GSE131907": "#b85c38",
    "GSE205335": "#1e7f4f",
    "GSE189357": "#6c3483",
}


def say(msg: str) -> None:
    print(msg, flush=True)


def load_locked() -> pd.DataFrame:
    df = pd.read_csv(LOCKED_PATH, sep="\t", dtype={"unit_id": str})
    df["n_cells"] = df["n_cells"].astype(int)
    df["n_malignant"] = df["n_malignant"].astype(int)
    df["n_tnk"] = df["n_tnk"].astype(int)
    return df


def load_dorothea() -> pd.DataFrame:
    df = pd.read_csv(DOROTHEA_PATH, sep="\t")
    df["source"] = df["source"].str.upper()
    df["target"] = df["target"].str.upper()
    df["mor"] = df["mor"].astype(float)
    return df


def spearman_pair(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = int(x.size)
    if n < 4 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan, np.nan, n
    r, p = stats.spearmanr(x, y)
    return float(r), float(p), n


def dl_spearman(rhos, ns):
    rhos = np.asarray(rhos, dtype=float)
    ns = np.asarray(ns, dtype=float)
    ok = np.isfinite(rhos) & np.isfinite(ns) & (ns > 3)
    rhos, ns = rhos[ok], ns[ok]
    if rhos.size == 0:
        return dict(rho=np.nan, p=np.nan, I2=np.nan, ci_lo=np.nan, ci_hi=np.nan, k=0, N=0)
    z = np.arctanh(np.clip(rhos, -0.999999, 0.999999))
    var = 1.0 / (ns - 3.0)
    w = 1.0 / var
    zbar = np.sum(w * z) / np.sum(w)
    q = np.sum(w * (z - zbar) ** 2)
    k = int(rhos.size)
    dfree = k - 1
    cdenom = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (q - dfree) / cdenom) if dfree > 0 and cdenom > 0 else 0.0
    wstar = 1.0 / (var + tau2)
    zre = np.sum(wstar * z) / np.sum(wstar)
    se = float(np.sqrt(1.0 / np.sum(wstar)))
    p = float(2 * stats.norm.sf(abs(zre / se))) if se > 0 else np.nan
    i2 = float(max(0.0, (q - dfree) / q)) if q > 0 else 0.0
    ci = np.tanh(zre + np.array([-1.0, 1.0]) * 1.96 * se)
    return dict(
        rho=float(np.tanh(zre)), p=p, I2=i2, ci_lo=float(ci[0]), ci_hi=float(ci[1]),
        k=k, N=int(ns.sum()),
    )


def within_quartile(x: pd.Series) -> pd.Series:
    r = x.rank(method="average")
    breaks = np.quantile(r.to_numpy(dtype=float), [0, 0.25, 0.5, 0.75, 1.0])
    return pd.cut(r, bins=breaks, include_lowest=True, labels=["Q1", "Q2", "Q3", "Q4"]).astype(str)


def rank_biserial(x4, x1):
    x4 = np.asarray(x4, dtype=float)
    x1 = np.asarray(x1, dtype=float)
    x4, x1 = x4[np.isfinite(x4)], x1[np.isfinite(x1)]
    res = stats.mannwhitneyu(x4, x1, alternative="two-sided", method="asymptotic")
    n4, n1 = int(x4.size), int(x1.size)
    r = 2.0 * float(res.statistic) / (n4 * n1) - 1.0
    return r, float(res.pvalue), n1, n4


def self_check(locked: pd.DataFrame) -> None:
    rows = []
    for ds in COHORTS:
        d = locked[locked.dataset == ds]
        r, p, n = spearman_pair(d.mal_CLDN4_pct, d.frac_tnk)
        rows.append((r, n))
        if ds == "GSE123902" and abs(r + 0.6593406593406593) > 1e-9:
            raise SystemExit(f"self-check GSE123902 rho {r}")
    dl = dl_spearman([r for r, _ in rows], [n for _, n in rows])
    if abs(dl["rho"] + 0.5311678045689989) > 1e-9 or dl["N"] != 65:
        raise SystemExit(f"self-check DL failed: {dl}")
    q = []
    for ds in COHORTS:
        d = locked[locked.dataset == ds]
        qq = within_quartile(d.mal_CLDN4_pct)
        q.append(pd.Series(qq.to_numpy(), index=d.index))
    q = pd.concat(q)
    rb = rank_biserial(locked.loc[q.eq("Q4"), "frac_tnk"], locked.loc[q.eq("Q1"), "frac_tnk"])
    if abs(rb[0] + 0.7236842105263157) > 1e-9 or rb[2:] != (19, 16):
        raise SystemExit(f"self-check quartile failed: {rb}")
    say(f"self-check ok  DL rho={dl['rho']:.6f}  Q4vsQ1 r={rb[0]:.6f}")


def dedupe_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).upper() for c in df.columns]
    return df.loc[:, ~df.columns.duplicated()]


def masks_from_named(get, n: int):
    def g(name):
        v = get(name)
        if v is None:
            return np.zeros(n, dtype=np.float64)
        return np.asarray(v, dtype=np.float64)

    mal = ((g("EPCAM") > 0) | (g("KRT8") > 0) | (g("KRT18") > 0) | (g("KRT19") > 0)) & (g("PTPRC") == 0)
    tnk = (
        (g("CD3D") > 0) | (g("CD3E") > 0) | (g("CD8A") > 0) | (g("NKG7") > 0) | (g("GNLY") > 0) | (g("KLRD1") > 0)
    ) & ~mal
    return mal, tnk


def pair_rows(dataset, unit_id, counts: dict[str, np.ndarray], n_mal: int):
    rows = []
    if n_mal < 30:
        return rows
    for a, b, level in PAIRS:
        if a not in counts or b not in counts:
            continue
        rho, _p, n = spearman_pair(counts[a], counts[b])
        rows.append(dict(dataset=dataset, unit_id=unit_id, gene_a=a, gene_b=b, level=level, rho=rho, n_malignant=n))
    return rows


def copos_row(dataset, unit_id, counts: dict[str, np.ndarray], n_mal: int):
    mats = []
    marg = {}
    for g in MODULE4:
        if g not in counts:
            return None
        pos = np.asarray(counts[g]) > 0
        mats.append(pos)
        marg[g] = float(pos.mean()) if n_mal else np.nan
    all4 = np.logical_and.reduce(mats)
    p_all = float(all4.mean()) if n_mal else np.nan
    prod = float(np.prod([marg[g] for g in MODULE4]))
    ratio = p_all / prod if prod > 0 else np.nan
    row = dict(dataset=dataset, unit_id=str(unit_id), n_malignant=n_mal, p_all4=p_all, ratio=ratio)
    for g in MODULE4:
        row[f"p_{g}"] = marg[g]
    return row


def save_cache(dataset, units, all_mean, mean_cp, focus_pct, focus_mean, cell_rho, copos):
    out = CACHE / dataset
    out.mkdir(parents=True, exist_ok=True)
    units.to_csv(out / "units.tsv", sep="\t", index=False)
    all_mean.to_csv(out / "all_mean_log1p.tsv.gz", sep="\t", compression="gzip")
    mean_cp.to_csv(out / "mean_cp10k.tsv.gz", sep="\t", compression="gzip")
    focus_pct.to_csv(out / "focus_pct.tsv", sep="\t")
    focus_mean.to_csv(out / "focus_mean.tsv", sep="\t")
    cell_rho.to_csv(out / "cell_rho.tsv", sep="\t", index=False)
    copos.to_csv(out / "copos.tsv", sep="\t", index=False)
    say(f"cached {dataset} units={len(units)} genes={all_mean.shape[0]}")


def frames_from_columns(columns: dict[str, dict[str, float]], unit_ids) -> pd.DataFrame:
    df = pd.DataFrame.from_dict(columns, orient="index")
    df = df.reindex(columns=[str(u) for u in unit_ids])
    return df.astype(float).fillna(0.0)


def check_against_locked(dataset, units: pd.DataFrame, focus_pct: pd.DataFrame, focus_mean: pd.DataFrame, locked: pd.DataFrame):
    exp = locked[locked.dataset == dataset].set_index("unit_id")
    got = units.set_index("unit_id")
    missing = [u for u in exp.index if u not in got.index]
    if missing:
        raise SystemExit(f"{dataset} missing units {missing}")
    got = got.loc[exp.index]
    for col in ["n_cells", "n_malignant", "n_tnk"]:
        if not np.array_equal(got[col].to_numpy(), exp[col].to_numpy()):
            bad = got.index[got[col].to_numpy() != exp[col].to_numpy()]
            raise SystemExit(f"{dataset} {col} mismatch {list(bad[:8])}\n{got.loc[bad, col].head()}\n{exp.loc[bad, col].head()}")
    pct = focus_pct.loc["CLDN4", exp.index.astype(str)].to_numpy(dtype=float)
    mean = focus_mean.loc["CLDN4", exp.index.astype(str)].to_numpy(dtype=float)
    if np.nanmax(np.abs(pct - exp.mal_CLDN4_pct.to_numpy())) > 1e-4:
        raise SystemExit(f"{dataset} CLDN4 %pos mismatch max={np.nanmax(np.abs(pct - exp.mal_CLDN4_pct.to_numpy()))}")
    if np.nanmax(np.abs(mean - exp.mal_CLDN4_mean.to_numpy())) > 1e-4:
        raise SystemExit(f"{dataset} CLDN4 mean mismatch max={np.nanmax(np.abs(mean - exp.mal_CLDN4_mean.to_numpy()))}")
    say(f"locked check ok {dataset}")


def summarize_unit_expression(dataset, unit_id, counts, genes_for_cp, lib, n_mal):
    pct, mean, cp = {}, {}, {}
    for g, raw in counts.items():
        raw = np.asarray(raw, dtype=np.float64)
        if g in FOCUS:
            pct[g] = 100.0 * float(np.mean(raw > 0))
            mean[g] = float(np.mean(np.log1p(raw)))
        if g in genes_for_cp:
            scale = np.maximum(np.asarray(lib, dtype=np.float64), 1.0)
            cp[g] = float(np.mean(np.log1p(raw / scale * 1e4)))
    return pct, mean, cp, pair_rows(dataset, unit_id, counts, n_mal), copos_row(dataset, unit_id, counts, n_mal)


def load_gse123902(geo: Path, locked: pd.DataFrame, cp_genes: set[str]):
    dataset = "GSE123902"
    ids = locked.loc[locked.dataset == dataset, "unit_id"].astype(str).tolist()
    tar_path = geo / "gse123902" / "GSE123902_RAW.tar"
    say(f"== {dataset} ==")
    unit_rows, cell_rows, copos_rows = [], [], []
    all_cols, pct_cols, mean_cols, cp_cols = {}, {}, {}, {}
    with tarfile.open(tar_path) as tar:
        for unit in ids:
            member = GSE123902_FILES[unit]
            say(f"  {unit}")
            raw = tar.extractfile(member)
            with gzip.GzipFile(fileobj=raw) as gz:
                df = pd.read_csv(gz, index_col=0)
            df = dedupe_columns(df)
            mal, tnk = masks_from_named(lambda name: df[name].to_numpy() if name in df.columns else None, len(df))
            n_mal = int(mal.sum())
            unit_rows.append(dict(
                dataset=dataset, unit_id=unit, n_cells=int(len(df)),
                n_malignant=n_mal, n_tnk=int(tnk.sum()),
                frac_tnk=float(tnk.mean()),
            ))
            mal_df = df.loc[mal]
            means = np.log1p(mal_df.to_numpy(dtype=np.float64)).mean(axis=0)
            for g, v in zip(mal_df.columns, means):
                all_cols.setdefault(g, {})[unit] = float(v)
            lib = df.to_numpy(dtype=np.float64).sum(axis=1)[mal]
            counts = {g: mal_df[g].to_numpy(dtype=np.float64) for g in set(FOCUS) | cp_genes if g in mal_df.columns}
            pct, mean, cp, prows, crow = summarize_unit_expression(dataset, unit, counts, cp_genes | set(FOCUS), lib, n_mal)
            for g, v in pct.items():
                pct_cols.setdefault(g, {})[unit] = v
            for g, v in mean.items():
                mean_cols.setdefault(g, {})[unit] = v
            for g, v in cp.items():
                cp_cols.setdefault(g, {})[unit] = v
            cell_rows.extend(prows)
            if crow:
                copos_rows.append(crow)
            del df, mal_df
    units = pd.DataFrame(unit_rows)
    return (
        units,
        frames_from_columns(all_cols, ids),
        frames_from_columns(cp_cols, ids),
        frames_from_columns(pct_cols, ids),
        frames_from_columns(mean_cols, ids),
        pd.DataFrame(cell_rows),
        pd.DataFrame(copos_rows),
    )


def _first_symbol_index(symbols: list[str]) -> dict[str, int]:
    seen = {}
    for i, g in enumerate(symbols):
        g = g.upper()
        if g and g not in seen:
            seen[g] = i
    return seen


def load_gse189357(geo: Path, locked: pd.DataFrame, cp_genes: set[str]):
    from scipy.io import mmread

    dataset = "GSE189357"
    ids = locked.loc[locked.dataset == dataset, "unit_id"].astype(str).tolist()
    raw_dir = geo / "gse189357" / "raw"
    if not raw_dir.exists():
        raw_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(geo / "gse189357" / "GSE189357_RAW.tar") as tar:
            tar.extractall(raw_dir)
    say(f"== {dataset} ==")
    unit_rows, cell_rows, copos_rows = [], [], []
    all_cols, pct_cols, mean_cols, cp_cols = {}, {}, {}, {}
    for unit in ids:
        say(f"  {unit}")
        mtx = next(raw_dir.glob(f"*_{unit}_matrix.mtx.gz"))
        feat = next(raw_dir.glob(f"*_{unit}_features.tsv.gz"))
        symbols = []
        with gzip.open(feat, "rt") as handle:
            for line in handle:
                symbols.append(line.rstrip("\n").split("\t")[1].upper())
        seen = _first_symbol_index(symbols)
        with gzip.open(mtx, "rb") as handle:
            mat = mmread(handle).tocsr()
        keep = np.zeros(mat.shape[0], dtype=bool)
        keep[list(seen.values())] = True
        if keep.sum() != mat.shape[0]:
            mat = mat[keep]
        order = [seen[g] for g in seen]
        # seen values are original row indices; after boolean keep, rows are those indices in ascending order
        # rebuild symbol per remaining row from the boolean mask order
        kept_symbols = [s for i, s in enumerate(symbols) if keep[i] and seen.get(s) == i]
        if len(kept_symbols) != mat.shape[0]:
            raise SystemExit(f"{unit} symbol/matrix mismatch {len(kept_symbols)} vs {mat.shape[0]}")
        sym_index = {g: i for i, g in enumerate(kept_symbols)}

        def getter(name, _idx=sym_index, _mat=mat):
            i = _idx.get(name)
            if i is None:
                return None
            return np.asarray(_mat.getrow(i).toarray()).ravel()

        n = mat.shape[1]
        mal, tnk = masks_from_named(getter, n)
        n_mal = int(mal.sum())
        unit_rows.append(dict(
            dataset=dataset, unit_id=unit, n_cells=int(n),
            n_malignant=n_mal, n_tnk=int(tnk.sum()), frac_tnk=float(tnk.mean()),
        ))
        mal_idx = np.flatnonzero(mal)
        piece = mat[:, mal_idx].tocoo()
        acc = np.zeros(mat.shape[0], dtype=np.float64)
        if piece.nnz:
            np.add.at(acc, piece.row, np.log1p(piece.data.astype(np.float64)))
        acc /= max(n_mal, 1)
        for g, i in sym_index.items():
            if acc[i] > 0 or g in FOCUS:
                all_cols.setdefault(g, {})[unit] = float(acc[i])
        lib = np.asarray(mat.sum(axis=0)).ravel()[mal_idx]
        want = [g for g in set(FOCUS) | cp_genes if g in sym_index]
        counts = {g: np.asarray(mat.getrow(sym_index[g]).toarray()).ravel()[mal_idx] for g in want}
        pct, mean, cp, prows, crow = summarize_unit_expression(dataset, unit, counts, cp_genes | set(FOCUS), lib, n_mal)
        for g, v in pct.items():
            pct_cols.setdefault(g, {})[unit] = v
        for g, v in mean.items():
            mean_cols.setdefault(g, {})[unit] = v
        for g, v in cp.items():
            cp_cols.setdefault(g, {})[unit] = v
        cell_rows.extend(prows)
        if crow:
            copos_rows.append(crow)
        del mat
    return (
        pd.DataFrame(unit_rows),
        frames_from_columns(all_cols, ids),
        frames_from_columns(cp_cols, ids),
        frames_from_columns(pct_cols, ids),
        frames_from_columns(mean_cols, ids),
        pd.DataFrame(cell_rows),
        pd.DataFrame(copos_rows),
    )


def load_gse131907(geo: Path, locked: pd.DataFrame, cp_genes: set[str]):
    dataset = "GSE131907"
    ids = locked.loc[locked.dataset == dataset, "unit_id"].astype(str).tolist()
    say(f"== {dataset} ==")
    ann = pd.read_csv(geo / "gse131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t", dtype=str)
    ann = ann[ann["Sample"].isin(ids)]
    n_cells = ann.groupby("Sample").size().reindex(ids).astype(int)
    n_mal_ann = ann["Cell_subtype"].eq("Malignant cells").groupby(ann["Sample"]).sum().reindex(ids).astype(int)
    n_tnk = ann["Cell_type"].isin(["T lymphocytes", "NK cells"]).groupby(ann["Sample"]).sum().reindex(ids).astype(int)
    mat_path = geo / "gse131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    with gzip.open(mat_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    barcodes = header[1:]
    meta = ann.set_index("Index").reindex(barcodes)
    keep = meta["Sample"].isin(ids).to_numpy() & meta["Cell_subtype"].eq("Malignant cells").to_numpy()
    keep = np.where(meta["Sample"].notna().to_numpy(), keep, False)
    ann_mal = set(ann.loc[ann["Cell_subtype"].eq("Malignant cells"), "Index"])
    missing = ann_mal - set(barcodes)
    if missing:
        raise SystemExit(f"GSE131907 malignant cells missing from matrix: {len(missing)}")
    keep_idx = np.flatnonzero(keep)
    unit_of = meta["Sample"].to_numpy()[keep_idx]
    order = np.argsort(unit_of, kind="mergesort")
    unit_of = unit_of[order]
    keep_idx = keep_idx[order]
    uniq, starts, counts = np.unique(unit_of, return_index=True, return_counts=True)
    if list(uniq) != sorted(ids) or not np.array_equal(counts, n_mal_ann.loc[list(uniq)].to_numpy()):
        raise SystemExit(f"GSE131907 malignant counts != annotation {list(zip(uniq, counts))[:5]}")
    starts = starts.astype(np.int64)
    store = set(FOCUS) | cp_genes
    stored = {}
    all_mean = {}
    lib = np.zeros(keep_idx.size, dtype=np.float64)
    seen = set()
    say(f"  streaming malignant cells={keep_idx.size}")
    with gzip.open(mat_path, "rt") as handle:
        handle.readline()
        n_genes = 0
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.upper()
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != len(barcodes):
                raise SystemExit(f"{gene}: {arr.size} != {len(barcodes)}")
            if gene in seen or not gene:
                continue
            seen.add(gene)
            kept = arr[keep_idx]
            lib += kept
            sums = np.add.reduceat(np.log1p(kept.astype(np.float64)), starts)
            if float(sums.sum()) > 0 or gene in FOCUS:
                all_mean[gene] = sums / counts
            if gene in store:
                stored[gene] = kept
            n_genes += 1
            if n_genes % 4000 == 0:
                say(f"  genes={n_genes} stored={len(stored)}")
    say(f"  stream done genes={n_genes}")
    scale = np.maximum(lib, 1.0)
    pct_cols, mean_cols, cp_cols = {}, {}, {}
    cell_rows, copos_rows = [], []
    # regroup stored rows, which follow the sorted unit_of order
    for i, unit in enumerate(uniq):
        sl = slice(int(starts[i]), int(starts[i] + counts[i]))
        counts_u = {g: stored[g][sl] for g in stored}
        n_mal = int(counts[i])
        pct, mean, cp, prows, crow = summarize_unit_expression(
            dataset, unit, counts_u, cp_genes | set(FOCUS), scale[sl], n_mal
        )
        for g, v in pct.items():
            pct_cols.setdefault(g, {})[unit] = v
        for g, v in mean.items():
            mean_cols.setdefault(g, {})[unit] = v
        for g, v in cp.items():
            cp_cols.setdefault(g, {})[unit] = v
        cell_rows.extend(prows)
        if crow:
            copos_rows.append(crow)
    all_df = pd.DataFrame.from_dict({g: pd.Series(v, index=list(uniq)) for g, v in all_mean.items()}, orient="index")
    units = pd.DataFrame({
        "dataset": dataset,
        "unit_id": ids,
        "n_cells": n_cells.loc[ids].to_numpy(),
        "n_malignant": n_mal_ann.loc[ids].to_numpy(),
        "n_tnk": n_tnk.loc[ids].to_numpy(),
    })
    units["frac_tnk"] = units["n_tnk"] / units["n_cells"]
    return (
        units,
        all_df.reindex(columns=ids),
        frames_from_columns(cp_cols, ids),
        frames_from_columns(pct_cols, ids),
        frames_from_columns(mean_cols, ids),
        pd.DataFrame(cell_rows),
        pd.DataFrame(copos_rows),
    )


def ensure_gse205335_extract(geo: Path, cp_genes: set[str]) -> Path:
    out = geo / "gse205335" / "extract"
    need = ["counts.tsv", "all_mean_log1p.tsv", "cell_focus.tsv", "target_cp10k.tsv"]
    if all((out / name).exists() and (out / name).stat().st_size > 0 for name in need):
        return out
    rds = geo / "gse205335" / "GSE205335_Lung_IO_UMI_matrix.rds"
    gz = geo / "gse205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    if not rds.exists():
        say("peeling double-gzip RDS")
        with gzip.open(gz, "rb") as outer:
            magic = outer.read(2)
            outer.seek(0)
            if magic == b"\x1f\x8b":
                with gzip.GzipFile(fileobj=outer) as inner, rds.open("wb") as dest:
                    while True:
                        chunk = inner.read(16 * 1024 * 1024)
                        if not chunk:
                            break
                        dest.write(chunk)
            else:
                with rds.open("wb") as dest:
                    while True:
                        chunk = outer.read(16 * 1024 * 1024)
                        if not chunk:
                            break
                        dest.write(chunk)
    wanted = geo / "gse205335" / "wanted.txt"
    genes = sorted(cp_genes | set(FOCUS) | set(WATCH))
    wanted.write_text("\n".join(genes) + "\n")
    if not (geo / "gse205335" / "map.tsv").exists():
        raise SystemExit("missing GSE205335 map.tsv; re-run the map builder in METHODS")
    env = os.environ.copy()
    env.update({
        "RDS": str(rds),
        "MAP": str(geo / "gse205335" / "map.tsv"),
        "OUT": str(out),
        "WANTED": str(wanted),
    })
    import subprocess
    say("R extract GSE205335")
    subprocess.check_call(["Rscript", str(ROOT / "scripts" / "extract_gse205335.R")], env=env)
    return out


def load_gse205335(geo: Path, locked: pd.DataFrame, cp_genes: set[str]):
    dataset = "GSE205335"
    ids = locked.loc[locked.dataset == dataset, "unit_id"].astype(str).tolist()
    say(f"== {dataset} ==")
    out = ensure_gse205335_extract(geo, cp_genes)
    counts = pd.read_csv(out / "counts.tsv", sep="\t", dtype={"unit_id": str}).set_index("unit_id").loc[ids]
    units = counts.reset_index()
    units.insert(0, "dataset", dataset)
    units["frac_tnk"] = units["n_tnk"] / units["n_cells"]
    all_mean = pd.read_csv(out / "all_mean_log1p.tsv", sep="\t", index_col=0)
    all_mean.columns = all_mean.columns.astype(str)
    all_mean.index = all_mean.index.astype(str).str.upper()
    cp = pd.read_csv(out / "target_cp10k.tsv", sep="\t", index_col=0)
    cp.columns = cp.columns.astype(str)
    cp.index = cp.index.astype(str).str.upper()
    cells = pd.read_csv(out / "cell_focus.tsv", sep="\t", dtype={"unit_id": str})
    cells.columns = [c if c == "unit_id" else str(c).upper() for c in cells.columns]
    pct_cols, mean_cols = {}, {}
    cell_rows, copos_rows = [], []
    for unit, sub in cells.groupby("unit_id", sort=False):
        counts_u = {g: sub[g].to_numpy(dtype=np.float64) for g in FOCUS if g in sub.columns}
        n_mal = int(len(sub))
        for g, raw in counts_u.items():
            pct_cols.setdefault(g, {})[str(unit)] = 100.0 * float(np.mean(raw > 0))
            mean_cols.setdefault(g, {})[str(unit)] = float(np.mean(np.log1p(raw)))
        cell_rows.extend(pair_rows(dataset, str(unit), counts_u, n_mal))
        crow = copos_row(dataset, str(unit), counts_u, n_mal)
        if crow:
            copos_rows.append(crow)
    return (
        units,
        all_mean.reindex(columns=ids),
        cp.reindex(columns=ids),
        frames_from_columns(pct_cols, ids),
        frames_from_columns(mean_cols, ids),
        pd.DataFrame(cell_rows),
        pd.DataFrame(copos_rows),
    )


LOADERS = {
    "GSE123902": load_gse123902,
    "GSE189357": load_gse189357,
    "GSE131907": load_gse131907,
    "GSE205335": load_gse205335,
}


def cache_ready(dataset: str) -> bool:
    out = CACHE / dataset
    return all((out / name).exists() for name in [
        "units.tsv", "all_mean_log1p.tsv.gz", "mean_cp10k.tsv.gz",
        "focus_pct.tsv", "focus_mean.tsv", "cell_rho.tsv", "copos.tsv",
    ])


def run_cohort(dataset, geo, locked, cp_genes, force: bool):
    if cache_ready(dataset) and not force:
        say(f"cache hit {dataset}")
        return
    units, all_mean, mean_cp, focus_pct, focus_mean, cell_rho, copos = LOADERS[dataset](geo, locked, cp_genes)
    for frame in (all_mean, mean_cp, focus_pct, focus_mean):
        frame.columns = frame.columns.astype(str)
    check_against_locked(dataset, units, focus_pct, focus_mean, locked)
    save_cache(dataset, units, all_mean, mean_cp, focus_pct, focus_mean, cell_rho, copos)


def read_cache():
    units, cells, copos = [], [], []
    all_mean, mean_cp, focus_pct, focus_mean = {}, {}, {}, {}
    for ds in COHORTS:
        out = CACHE / ds
        u = pd.read_csv(out / "units.tsv", sep="\t", dtype={"unit_id": str})
        units.append(u)
        cells.append(pd.read_csv(out / "cell_rho.tsv", sep="\t", dtype={"unit_id": str}))
        copos.append(pd.read_csv(out / "copos.tsv", sep="\t", dtype={"unit_id": str}))
        all_mean[ds] = pd.read_csv(out / "all_mean_log1p.tsv.gz", sep="\t", index_col=0)
        mean_cp[ds] = pd.read_csv(out / "mean_cp10k.tsv.gz", sep="\t", index_col=0)
        focus_pct[ds] = pd.read_csv(out / "focus_pct.tsv", sep="\t", index_col=0)
        focus_mean[ds] = pd.read_csv(out / "focus_mean.tsv", sep="\t", index_col=0)
        for frame in (all_mean[ds], mean_cp[ds], focus_pct[ds], focus_mean[ds]):
            frame.index = frame.index.astype(str).str.upper()
            frame.columns = frame.columns.astype(str)
    return (
        pd.concat(units, ignore_index=True),
        pd.concat(cells, ignore_index=True),
        pd.concat(copos, ignore_index=True),
        all_mean, mean_cp, focus_pct, focus_mean,
    )


def z_within(patient: pd.DataFrame, cols: list[str], dataset_col="dataset") -> pd.Series:
    zcols = []
    for col in cols:
        z = pd.Series(np.nan, index=patient.index, dtype=float)
        for ds, idx in patient.groupby(dataset_col).groups.items():
            x = patient.loc[idx, col].astype(float)
            sd = float(x.std(ddof=1))
            if not np.isfinite(sd) or sd == 0:
                continue
            z.loc[idx] = (x - float(x.mean())) / sd
        zcols.append(z)
    return pd.concat(zcols, axis=1).mean(axis=1)


def wmean_score(cp: pd.DataFrame, edges: pd.DataFrame, drop: set[str] | None = None):
    drop = drop or set()
    use = edges[~edges["target"].isin(drop)]
    present = [t for t in use["target"] if t in cp.index]
    if len(present) < 5:
        return pd.Series(np.nan, index=cp.columns), len(present)
    w = use.drop_duplicates("target").set_index("target").loc[present, "mor"].astype(float)
    x = cp.loc[present].astype(float)
    score = x.mul(w, axis=0).sum(axis=0) / float(np.abs(w).sum())
    return score, len(present)


def spearman_against(mean_df: pd.DataFrame, tf: str) -> pd.Series:
    if tf not in mean_df.index:
        return pd.Series(dtype=float)
    x = mean_df.loc[tf].to_numpy(dtype=float)
    mat = mean_df.to_numpy(dtype=float)
    # columns are units; rank within each gene across units
    rx = stats.rankdata(x, method="average")
    rx = rx - rx.mean()
    rmat = stats.rankdata(mat, axis=1, method="average")
    rmat = rmat - rmat.mean(axis=1, keepdims=True)
    denom = np.sqrt((rmat ** 2).sum(axis=1) * np.sum(rx ** 2))
    num = rmat @ rx
    rho = np.full(mat.shape[0], np.nan)
    ok = denom > 0
    rho[ok] = num[ok] / denom[ok]
    return pd.Series(rho, index=mean_df.index)


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4g}"


def fmt_rho(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:.3f}"


def cohort_dl(patient: pd.DataFrame, score: str, outcome: str = "frac_tnk"):
    rows = []
    for ds in COHORTS:
        d = patient[patient.dataset == ds]
        rho, p, n = spearman_pair(d[score], d[outcome])
        rows.append(dict(dataset=ds, n=n, rho=rho, p=p))
    dl = dl_spearman([r["rho"] for r in rows], [r["n"] for r in rows])
    return rows, dl


def summarize(dor: pd.DataFrame):
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    units, cell_rho, copos, all_mean, mean_cp, focus_pct, focus_mean = read_cache()
    patient = units.copy()
    patient["unit_id"] = patient["unit_id"].astype(str)
    for ds in COHORTS:
        pct = focus_pct[ds]
        mean = focus_mean[ds]
        for g in FOCUS:
            if g not in pct.index:
                continue
            mask = patient.dataset.eq(ds)
            patient.loc[mask, f"pct_{g}"] = pct.loc[g, patient.loc[mask, "unit_id"]].to_numpy(dtype=float)
            patient.loc[mask, f"mean_{g}"] = mean.loc[g, patient.loc[mask, "unit_id"]].to_numpy(dtype=float)
        cp = mean_cp[ds]
        elf = dor[dor.source == "ELF3"]
        grh = dor[dor.source == "GRHL2"]
        for name, edges, drop in [
            ("wmean_ELF3", elf, set()),
            ("wmean_GRHL2", grh, set()),
            ("wmean_GRHL2_noCLDN4", grh, {"CLDN4"}),
        ]:
            score, n_t = wmean_score(cp, edges, drop)
            patient.loc[patient.dataset.eq(ds), name] = score.reindex(patient.loc[patient.dataset.eq(ds), "unit_id"]).to_numpy(dtype=float)
            patient.loc[patient.dataset.eq(ds), name + "_ntarget"] = n_t
    patient["module4_pct"] = z_within(patient, [f"pct_{g}" for g in MODULE4])
    patient["module3_pct"] = z_within(patient, [f"pct_{g}" for g in MODULE3])
    patient["module4_mean"] = z_within(patient, [f"mean_{g}" for g in MODULE4])
    patient["module3_mean"] = z_within(patient, [f"mean_{g}" for g in MODULE3])
    patient.to_csv(TABLES / "patient_scores.tsv", sep="\t", index=False)

    # Coexpression across units: mean log1p and %pos.
    co_rows = []
    for a, b, level in PAIRS:
        for scale, prefix in [("mean", "mean_log1p"), ("pct", "pct")]:
            rhos, ns, singles = [], [], []
            for ds in COHORTS:
                d = patient[patient.dataset == ds]
                rho, p, n = spearman_pair(d[f"{scale}_{a}"], d[f"{scale}_{b}"])
                rhos.append(rho)
                ns.append(n)
                singles.append(dict(dataset=ds, rho=rho, p=p, n=n))
            dl = dl_spearman(rhos, ns)
            co_rows.append(dict(
                gene_a=a, gene_b=b, pair_level=level, scale=prefix,
                rho=dl["rho"], p=dl["p"], I2=dl["I2"], ci_lo=dl["ci_lo"], ci_hi=dl["ci_hi"],
                k=dl["k"], N=dl["N"],
                **{f"rho_{s['dataset']}": s["rho"] for s in singles},
            ))
    co_df = pd.DataFrame(co_rows)
    co_df.to_csv(TABLES / "coexpression_pseudobulk.tsv", sep="\t", index=False)

    within_rows = []
    for (a, b, level), sub in cell_rho.groupby(["gene_a", "gene_b", "level"], sort=False):
        for ds in COHORTS + ["ALL"]:
            d = sub if ds == "ALL" else sub[sub.dataset == ds]
            rho = d["rho"].to_numpy(dtype=float)
            rho = rho[np.isfinite(rho)]
            p = np.nan
            if rho.size >= 5 and np.unique(np.round(rho, 6)).size > 1:
                p = float(stats.wilcoxon(rho, alternative="two-sided", zero_method="wilcox").pvalue)
            elif rho.size >= 5 and np.all(rho > 0):
                p = float(stats.wilcoxon(rho, alternative="two-sided", zero_method="wilcox").pvalue)
            within_rows.append(dict(
                gene_a=a, gene_b=b, pair_level=level, dataset=ds,
                n_units=int(rho.size), median_rho=float(np.median(rho)) if rho.size else np.nan,
                q25=float(np.quantile(rho, 0.25)) if rho.size else np.nan,
                q75=float(np.quantile(rho, 0.75)) if rho.size else np.nan,
                frac_pos=float(np.mean(rho > 0)) if rho.size else np.nan,
                wilcoxon_p=p,
            ))
    within_df = pd.DataFrame(within_rows)
    within_df.to_csv(TABLES / "coexpression_within_cell.tsv", sep="\t", index=False)

    copos_sum = []
    for ds in COHORTS + ["ALL"]:
        d = copos if ds == "ALL" else copos[copos.dataset == ds]
        copos_sum.append(dict(
            dataset=ds, n_units=len(d),
            median_p_all4=float(d.p_all4.median()),
            median_ratio=float(d.ratio.median()),
            **{f"median_p_{g}": float(d[f"p_{g}"].median()) for g in MODULE4},
        ))
    copos_sum = pd.DataFrame(copos_sum)
    copos_sum.to_csv(TABLES / "codetection.tsv", sep="\t", index=False)

    # Membership
    mem_rows = []
    for tf in ["ELF3", "GRHL2"]:
        edges = dor[dor.source == tf]
        for gene in ["TACSTD2", "CLDN4", "CLDN7", "CDH1", "GRHL2", "KRT8", "ELF3"]:
            hit = edges[edges.target == gene]
            mem_rows.append(dict(
                tf=tf, gene=gene, in_ABC=bool(len(hit)),
                mor=float(hit.mor.iloc[0]) if len(hit) else np.nan,
                confidence=hit.confidence.iloc[0] if len(hit) else "",
                n_targets_ABC=int(edges.target.nunique()),
            ))
    mem = pd.DataFrame(mem_rows)
    mem.to_csv(TABLES / "regulon_membership.tsv", sep="\t", index=False)

    coh_rows = []
    for tf in ["ELF3", "GRHL2"]:
        targets = sorted(dor.loc[dor.source == tf, "target"].unique())
        for ds in COHORTS:
            rho = spearman_against(all_mean[ds], tf)
            vals = rho.reindex(targets).dropna()
            coh_rows.append(dict(
                tf=tf, dataset=ds, n_targets=int(vals.size),
                median_rho=float(vals.median()) if len(vals) else np.nan,
                frac_pos=float((vals > 0).mean()) if len(vals) else np.nan,
            ))
    coh = pd.DataFrame(coh_rows)
    coh.to_csv(TABLES / "regulon_coherence.tsv", sep="\t", index=False)

    spec_rows = []
    for tf in ["ELF3", "GRHL2"]:
        for ds in COHORTS:
            rho = spearman_against(all_mean[ds], tf)
            finite = rho[np.isfinite(rho)]
            for gene in WATCH:
                if gene not in finite.index:
                    spec_rows.append(dict(tf=tf, dataset=ds, gene=gene, rho=np.nan, percentile=np.nan, n_genes=int(finite.size)))
                    continue
                r = float(finite.loc[gene])
                pct = 100.0 * float(np.mean(finite.to_numpy() <= r))
                spec_rows.append(dict(tf=tf, dataset=ds, gene=gene, rho=r, percentile=pct, n_genes=int(finite.size)))
    spec = pd.DataFrame(spec_rows)
    spec.to_csv(TABLES / "specificity_percentiles.tsv", sep="\t", index=False)

    score_defs = [
        ("pct_ELF3", "ELF3 %pos", "primary"),
        ("pct_TACSTD2", "TACSTD2 %pos", "primary"),
        ("pct_CLDN4", "CLDN4 %pos (pipeline check)", "check"),
        ("pct_CLDN7", "CLDN7 %pos", "primary"),
        ("pct_GRHL2", "GRHL2 %pos", "primary"),
        ("module4_pct", "four-gene %pos module", "primary"),
        ("module3_pct", "ELF3+TACSTD2+CLDN7 %pos (CLDN4 out)", "primary"),
        ("wmean_ELF3", "ELF3 DoRothEA wmean", "primary"),
        ("wmean_GRHL2", "GRHL2 DoRothEA wmean", "primary"),
        ("wmean_GRHL2_noCLDN4", "GRHL2 wmean, CLDN4 out", "primary"),
        ("module4_mean", "four-gene mean-log1p module", "companion"),
        ("module3_mean", "three-gene mean-log1p module", "companion"),
    ]
    tnk_rows = []
    single_rows = []
    q_rows = []
    for key, label, kind in score_defs:
        singles, dl = cohort_dl(patient, key)
        tnk_rows.append(dict(
            score=key, label=label, kind=kind,
            rho=dl["rho"], p=dl["p"], I2=dl["I2"], ci_lo=dl["ci_lo"], ci_hi=dl["ci_hi"],
            k=dl["k"], N=dl["N"],
        ))
        for s in singles:
            single_rows.append(dict(score=key, **s))
        qmap = []
        for ds in COHORTS:
            d = patient[patient.dataset == ds]
            qmap.append(pd.Series(within_quartile(d[key]).to_numpy(), index=d.index))
        q = pd.concat(qmap)
        rb = rank_biserial(patient.loc[q.eq("Q4"), "frac_tnk"], patient.loc[q.eq("Q1"), "frac_tnk"])
        q_rows.append(dict(score=key, r=rb[0], p=rb[1], n_q1=rb[2], n_q4=rb[3]))
    tnk = pd.DataFrame(tnk_rows)
    tnk.to_csv(TABLES / "tnk_association.tsv", sep="\t", index=False)
    pd.DataFrame(single_rows).to_csv(TABLES / "tnk_singles.tsv", sep="\t", index=False)
    qdf = pd.DataFrame(q_rows)
    qdf.to_csv(TABLES / "tnk_q4q1.tsv", sep="\t", index=False)

    write_figures(co_df, within_df, tnk, patient)
    write_finding(patient, co_df, within_df, copos_sum, mem, coh, spec, tnk, qdf)
    say("summarize done")
    return tnk


def forest(df, title, xlabel, path: Path):
    fig_h = 0.42 * len(df) + 1.3
    fig, ax = plt.subplots(figsize=(8.2, fig_h))
    y = np.arange(len(df))[::-1]
    ax.axvline(0, color="#888888", lw=0.8)
    xerr = np.vstack([df["rho"] - df["lo"], df["hi"] - df["rho"]])
    ax.errorbar(df["rho"], y, xerr=xerr, fmt="o", color="#1f4e79", ecolor="#1f4e79", capsize=3, ms=5)
    ax.set_yticks(y)
    ax.set_yticklabels(df["label"])
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.set_xlim(-1.05, 1.05)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_figures(co_df, within_df, tnk, patient):
    sub = co_df[(co_df.pair_level == "primary") & (co_df.scale == "mean_log1p")].copy()
    sub["label"] = sub.gene_a + "–" + sub.gene_b
    sub = sub.rename(columns={"ci_lo": "lo", "ci_hi": "hi"})
    forest(
        sub, "Malignant pseudobulk coexpression",
        "Spearman ρ of mean log1p(UMI), DerSimonian–Laird",
        FIGURES / "fig_coexpression_pseudobulk.png",
    )
    w = within_df[(within_df.pair_level == "primary") & (within_df.dataset == "ALL")].copy()
    w["label"] = w.gene_a + "–" + w.gene_b
    w = w.rename(columns={"median_rho": "rho", "q25": "lo", "q75": "hi"})
    forest(
        w, "Within-malignant-cell coexpression",
        "Median within-unit Spearman ρ (bar = interquartile range)",
        FIGURES / "fig_coexpression_within_cell.png",
    )
    t = tnk[tnk.kind.isin(["primary", "check"])].copy()
    t = t.rename(columns={"ci_lo": "lo", "ci_hi": "hi"})
    forest(
        t, "Patient-level association with T/NK fraction",
        "Spearman ρ versus frac_tnk, DerSimonian–Laird",
        FIGURES / "fig_tnk_forest.png",
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.6), sharey=True)
    for ax, col, title in zip(
        axes,
        ["module4_pct", "module3_pct"],
        ["Four-gene %pos module", "CLDN4 held out"],
    ):
        for ds in COHORTS:
            d = patient[patient.dataset == ds]
            ax.scatter(d[col], d.frac_tnk, s=28, color=COLORS[ds], label=ds, alpha=0.9)
        ax.set_xlabel(title + "\n(within-cohort z mean)")
        ax.axvline(0, color="#cccccc", lw=0.8)
    axes[0].set_ylabel("T/NK fraction")
    axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle("Concordant-4 malignant module versus T/NK", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_module_vs_tnk.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGURES / "fig_module_vs_tnk.pdf", bbox_inches="tight")
    plt.close(fig)


def md_table(headers, rows) -> str:
    line = "| " + " | ".join(headers) + " |"
    sep = "|" + "|".join("---" for _ in headers) + "|"
    body = ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join([line, sep, *body])


def write_finding(patient, co_df, within_df, copos_sum, mem, coh, spec, tnk, qdf):
    n_by = patient.groupby("dataset").size().to_dict()
    check = tnk[tnk.score == "pct_CLDN4"].iloc[0]
    qcheck = qdf[qdf.score == "pct_CLDN4"].iloc[0]

    def co_line(a, b):
        pb = co_df[(co_df.gene_a == a) & (co_df.gene_b == b) & (co_df.scale == "mean_log1p")].iloc[0]
        wc = within_df[(within_df.gene_a == a) & (within_df.gene_b == b) & (within_df.dataset == "ALL")].iloc[0]
        return [
            f"{a}–{b}",
            f"{wc.median_rho:.3f}",
            str(int(wc.n_units)),
            f"{wc.frac_pos:.2f}",
            fmt_p(wc.wilcoxon_p),
            f"{pb.rho:.3f}",
            fmt_p(pb.p),
            f"{pb.I2 * 100:.1f}%",
            f"{pb.ci_lo:.3f} to {pb.ci_hi:.3f}",
        ]

    primary_pairs = [(a, b) for a, b, level in PAIRS if level == "primary"]
    co_md = md_table(
        ["pair", "within-cell median ρ", "units", "fraction ρ>0", "Wilcoxon p", "pseudobulk ρ", "p", "I²", "95% CI"],
        [co_line(a, b) for a, b in primary_pairs],
    )
    sec_pairs = [(a, b) for a, b, level in PAIRS if level == "secondary"]
    sec_md = md_table(
        ["pair", "within-cell median ρ", "units", "pseudobulk ρ", "p"],
        [[
            f"{a}–{b}",
            f"{within_df[(within_df.gene_a==a)&(within_df.gene_b==b)&(within_df.dataset=='ALL')].iloc[0].median_rho:.3f}",
            str(int(within_df[(within_df.gene_a==a)&(within_df.gene_b==b)&(within_df.dataset=='ALL')].iloc[0].n_units)),
            f"{co_df[(co_df.gene_a==a)&(co_df.gene_b==b)&(co_df.scale=='mean_log1p')].iloc[0].rho:.3f}",
            fmt_p(co_df[(co_df.gene_a==a)&(co_df.gene_b==b)&(co_df.scale=='mean_log1p')].iloc[0].p),
        ] for a, b in sec_pairs],
    )
    cop = copos_sum[copos_sum.dataset == "ALL"].iloc[0]
    cop_md = md_table(
        ["cohort", "units", "median P(all four >0)", "median ratio vs independence", "ELF3", "TACSTD2", "CLDN4", "CLDN7"],
        [[
            r.dataset, str(int(r.n_units)), f"{r.median_p_all4:.3f}", f"{r.median_ratio:.2f}",
            f"{r.median_p_ELF3:.2f}", f"{r.median_p_TACSTD2:.2f}", f"{r.median_p_CLDN4:.2f}", f"{r.median_p_CLDN7:.2f}",
        ] for r in copos_sum.itertuples()],
    )
    mem_md = md_table(
        ["TF", "gene", "in DoRothEA ABC", "confidence", "targets in ABC"],
        [[r.tf, r.gene, "yes" if r.in_ABC else "no", r.confidence or "—", str(int(r.n_targets_ABC))]
         for r in mem.itertuples()
         if r.gene in {"TACSTD2", "CLDN4", "CLDN7", "CDH1", "GRHL2", "KRT8"} and r.tf != r.gene],
    )
    coh_md = md_table(
        ["TF", "cohort", "targets with a correlation", "median ρ vs TF", "fraction ρ>0"],
        [[r.tf, r.dataset, str(int(r.n_targets)), f"{r.median_rho:.3f}", f"{r.frac_pos:.2f}"] for r in coh.itertuples()],
    )
    spec_focus = spec[spec.gene.isin(["TACSTD2", "CLDN4", "CLDN7", "GRHL2", "ELF3", "EPCAM", "KRT8", "KRT5", "CDH1", "PTPRC"])]
    spec_md = md_table(
        ["TF", "gene", "median ρ across cohorts", "median percentile"],
        [[tf, gene,
          f"{float(np.nanmedian(spec_focus[(spec_focus.tf==tf)&(spec_focus.gene==gene)].rho)):.3f}",
          f"{float(np.nanmedian(spec_focus[(spec_focus.tf==tf)&(spec_focus.gene==gene)].percentile)):.1f}"]
         for tf in ["ELF3", "GRHL2"] for gene in ["GRHL2", "ELF3", "TACSTD2", "CLDN4", "CLDN7", "CDH1", "EPCAM", "KRT8", "KRT5", "PTPRC"]
         if gene != tf],
    )
    tnk_md = md_table(
        ["score", "kind", "N", "ρ", "p", "I²", "95% CI", "Q4 vs Q1 r (nQ1/nQ4, p)"],
        [[
            r.label, r.kind, str(int(r.N)), f"{r.rho:.3f}", fmt_p(r.p), f"{r.I2*100:.1f}%",
            f"{r.ci_lo:.3f} to {r.ci_hi:.3f}",
            (lambda q: f"{q.r:.3f} ({int(q.n_q1)}/{int(q.n_q4)}, {fmt_p(q.p)})")(qdf[qdf.score == r.score].iloc[0]),
        ] for r in tnk.itertuples()],
    )
    singles = []
    for ds in COHORTS:
        d = patient[patient.dataset == ds]
        rho, p, n_units = spearman_pair(d.pct_CLDN4, d.frac_tnk)
        singles.append(f"| {ds} | {n_units} | {rho:.3f} | {fmt_p(p)} |")
    text = f"""# Concordant-4 ELF3–TACSTD2–CLDN4–CLDN7

ADDITIVE. Malignant cells in the locked concordant-4 only:
GSE123902 + GSE131907 + GSE205335 + GSE189357.
Not GSE148071, GSE127465, GSE207422, or GSE154826.
TACSTD2 is not a gate. The unit is the patient / donor / sample.
Do not quote cell counts as n.

DoRothEA A+B+C is a curated regulon (weighted mean on malignant
log1p CP10k). It is not lung ChIP and not a SCENIC/cisTarget run.
ELF3 → GRHL2 and GRHL2 → CLDN4 are confidence-C edges in that
network. TACSTD2 and CLDN7 are not ELF3 or GRHL2 targets there.

## Pipeline check

Malignant CLDN4 %pos versus T/NK fraction, same 65 units:
ρ = {check.rho:.3f} (p = {fmt_p(check.p)}, I² = {check.I2*100:.1f}%, N = {int(check.N)}).
Stacked within-cohort Q4 versus Q1 rank-biserial r = {qcheck.r:.3f}
({int(qcheck.n_q1)}/{int(qcheck.n_q4)}, p = {fmt_p(qcheck.p)}).
This matches the locked concordant-4 result and is not a new claim.

| cohort | n | ρ | p |
|---|---:|---:|---:|
{chr(10).join(singles)}

## Honest n

- **n_units = {len(patient)}** ({n_by.get('GSE123902',0)} donors + {n_by.get('GSE131907',0)} samples + {n_by.get('GSE205335',0)} patients + {n_by.get('GSE189357',0)} patients).
- Malignant cells in those units: {int(patient.n_malignant.sum())}. Not the test n.
- Within-cell Spearman uses units with at least 30 malignant cells (P4001 is out of that layer only).

## 1. Coexpression in malignant cells

Within-cell ρ is one Spearman per unit, on UMI counts. The Wilcoxon p
treats units as replicates. Pseudobulk ρ is DerSimonian–Laird across
the four cohort Spearmans of malignant mean log1p(UMI).

{co_md}

Secondary GRHL1 / GRHL3 context:

{sec_md}

## 2. Four-gene co-detection

Ratio is P(ELF3, TACSTD2, CLDN4, and CLDN7 all > 0) divided by the
product of the four detection rates. The all-unit median ratio is
{cop.median_ratio:.2f} (median joint detection {cop.median_p_all4:.3f}, {int(cop.n_units)} units).

{cop_md}

## 3. Regulon

{mem_md}

Coherence is the within-cohort Spearman of each ABC target versus its
TF, on malignant mean log1p. A median near zero means the curated
target list does not move together in these malignant cells.

{coh_md}

ELF3 wmean uses {int(patient.wmean_ELF3_ntarget.median())} targets.
GRHL2 wmean uses {int(patient.wmean_GRHL2_ntarget.median())} targets,
{int(patient.wmean_GRHL2_noCLDN4_ntarget.median())} after CLDN4 is removed.

## 4. Where the junction genes sit in the malignant correlation list

Percentile is within cohort, among genes with a finite Spearman against
the TF (100 = the top of the positive tail). The table is the median
of the four cohorts.

{spec_md}

## 5. Patient-level association with T/NK

`frac_tnk = n_tnk / n_cells` on the locked denominator.
The four-gene module is the mean of within-cohort z-scores.
The three-gene module drops CLDN4 so the T/NK number is not only the
pipeline check. Q4 versus Q1 is the stacked within-cohort quartile
contrast (rank-biserial r).

{tnk_md}

Within malignant cells, the four genes are coexpressed: every primary pair
has a positive Spearman in essentially every unit, with median ρ about
0.65–0.77 (63–64 units; P4001 is already out). The same pairs stay
positive on the patient pseudobulk (DL ρ 0.43–0.82). GRHL2 RNA sits
with them, more strongly across patients than inside a single cell
(within-cell median ρ about 0.3). GRHL1 is weaker. GRHL3 is near the
background.

Against T/NK, malignant ELF3 %pos is negative (ρ = −0.459, p = 0.0053),
in the same direction as the locked CLDN4 check. The four-gene %pos
module is also negative (ρ = −0.449, p = 0.032). Dropping CLDN4, the
continuous pool is −0.412 (p = 0.063, CI crosses zero) and the stacked
Q4 versus Q1 contrast is r = −0.480 (p = 0.016). TACSTD2 %pos, GRHL2
%pos, and both DoRothEA wmeans have intervals that include zero.
GSE205335 changes sign for several of those non-CLDN4 scores, which is
why their I² is high. ELF3's own detection tracks low T/NK more clearly
than the average of its DoRothEA targets (52 or 53 symbols, depending on the matrix).

CLDN4 is at the top of the malignant ELF3 correlation list (median
percentile about 99.7), with EPCAM and KRT8. TACSTD2 is in the upper
tail (about the 86th percentile) and is the stronger GRHL2 neighbor
(about the 97th). KRT5 and PTPRC sit in the lower tail of the ELF3 list.

## What this does not say

- A public prior edge is not lung ELF3 or GRHL binding.
- Pseudobulk coexpression is shared malignant-program variation across
  patients. Within-cell ρ is the same-cell measurement.
- GSE205335 mixes histologies. That cohort stays in the locked n.
- No dual-high TACSTD2∩CLDN4 split, no private 8KL matrix, no Visium.

## Reproduce

```bash
bash methods/concordant4_elf3_regulon/scripts/download.sh /tmp/geo_c4
python3 methods/concordant4_elf3_regulon/scripts/analyze.py --geo /tmp/geo_c4
```
"""
    (ROOT / "FINDING.md").write_text(text)
    say("wrote FINDING.md")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--geo", default="/tmp/geo_c4")
    parser.add_argument("--cohort", default="all", choices=["all", "summarize", *COHORTS])
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    locked = load_locked()
    self_check(locked)
    if args.self_check:
        return
    dor = load_dorothea()
    cp_genes = set(dor.loc[dor.source.isin(["ELF3", "GRHL2"]), "target"]) | set(FOCUS)
    geo = Path(args.geo)
    if args.cohort == "summarize":
        summarize(dor)
        return
    todo = COHORTS if args.cohort == "all" else [args.cohort]
    for ds in todo:
        run_cohort(ds, geo, locked, cp_genes, args.force)
    if args.cohort == "all" or args.cohort == "summarize":
        summarize(dor)
    elif all(cache_ready(ds) for ds in COHORTS):
        summarize(dor)


if __name__ == "__main__":
    main()
